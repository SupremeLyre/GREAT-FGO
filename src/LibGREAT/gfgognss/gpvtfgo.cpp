/**
 * @file         gpvtfgo.cpp
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        main code of GNSS RTK besed on factor graph optimization
 * @version      1.0
 * @date         2025-01-01
 *
 * @copyright Copyright (c) 2025, Wuhan University. All rights reserved.
 *
 */

#include "gpvtfgo.h"
#include "gfactor/ginitial_pose_factor.h"
#include "gmodels/gprecisebiasGPP.h"
#include <gutils/gcommon.cpp>

gfgomsf::t_gpvtfgo::t_gpvtfgo(string site, string site_base, t_gsetbase *gset, std::shared_ptr<spdlog::logger> spdlog,
                              t_gallproc *allproc)
    : t_gspp(site, gset, spdlog), t_gpvtflt(site, site_base, gset, spdlog, allproc), t_gfgo(gset), t_gfgo_para(gset)
{
    t_gbiasmodel *precise_bias(new t_gprecisebiasGPP(_allproc, _spdlog, gset));
    _gbias_model = precise_bias;
    for (int i = 0; i <= gwindow_size; i++)
    {
        _Pos[i].setZero();
        _Vel[i].setZero();
        _headers[i] = 0.0;
        _rover_window[i] = nullptr;
        _para_window->delAllParam();
    }
    _win_base_data.resize(gwindow_size + 1);
    cur_sat_prn.clear();
    _all_para_win.delAllParam();
    _DD_msg.clear();
    _vDD_msg.clear();
    shared_ptr<t_gamb_manager> amb_m(new t_gamb_manager(_band_index));
    _amb_manager = amb_m;
    _pos_constrain = false;
    _cntrep = 0;
    int sign = 1;
    if (_sampling > 1)
        _gtime_interval = int(sign * _sampling);
    else
        _gtime_interval = sign * _sampling;
    //_gtime_interval = 20;//?????????????????????????
}

gfgomsf::t_gpvtfgo::~t_gpvtfgo()
{
    delete _last_gnss_info;
    delete _last_gnss_marginalization_info;
    delete _gbias_model;
    for (int i = 0; i <= _rover_count; i++)
    {
        if (_rover_window[i] != nullptr)
        {
            delete _rover_window[i];
        }
    }
}

int gfgomsf::t_gpvtfgo::processBatch(const t_gtime &beg_r, const t_gtime &end_r, bool prtOut)
{
#ifdef BMUTEX
    boost::mutex::scoped_lock lock(_mutex);
#endif
    _gmutex.lock();

    if (_grec == nullptr)
    {
        ostringstream os;
        os << "ERROR: No object found (" << _site << "). Processing terminated!!! " << beg_r.str_ymdhms() << " -> "
           << end_r.str_ymdhms() << endl;
        if (_spdlog)
            SPDLOG_LOGGER_ERROR(_spdlog, string("t_gpvtfgo ") + os.str());
        _gmutex.unlock();
        return -1;
    }

    int sign = 1;
    double subint = 0.1;

    if (!_beg_end)
        sign = -1;

    InitProc(beg_r, end_r, &subint);

    t_gtime now(_beg_time);

    std::cerr << _site << ": Start GNSS Processing: " << now.str_ymdhms() << " " << _end_time.str_ymdhms() << endl;

    bool time_loop = true;
    while (time_loop)
    {
        if (_beg_end && (now < _end_time || now == _end_time))
            time_loop = true;
        else if (_beg_end && now > _end_time)
        {
            time_loop = false;
            break;
        }

        // synchronization
        if (now != _end_time)
        {
            if (!time_sync(now, _sampling, _scale, _spdlog))
            {                                      // now.sod()%(int)_sampling != 0 ){
                now.add_dsec(sign * subint / 100); // add_dsec used for synchronization!

                continue;
            }
            if (_sampling > 1)
                now.reset_dsec();
        }

        bool irc_slip = _slip_detect(now);

        int irc_epo = 0;
        irc_epo = processWindow(now); // sliding window model
        if (irc_epo < 0)
        {
            if (_sampling > 1)
                now.add_secs(int(sign * _sampling)); // =<1Hz data
            else
                now.add_dsec(sign * _sampling); //  >1Hz data

            continue;
        }
        else
            _success = true;

        if (_spdlog)
            SPDLOG_LOGGER_ERROR(_spdlog, string("t_gpvtfgo "), _site + now.str_ymdhms(" processing epoch: "));
        double percent = now.diff(_beg_time) / _end_time.diff(_beg_time) * 100.0;
        std::cerr << "\r" << now.str_ymdhms() << setw(5) << " Q = " << (_amb_state ? 1 : 2) << fixed << setprecision(1)
                  << setw(6) << percent << "%";

        if (_sampling > 1)
            now.add_secs(int(sign * _sampling)); // =<1Hz data
        else
            now.add_dsec(sign * _sampling); //  >1Hz data
    }
    _gmutex.unlock();
    return 1;
}
int gfgomsf::t_gpvtfgo::processWindow(const t_gtime &now, vector<t_gsatdata> *data_rover, vector<t_gsatdata> *data_base)
{
    if (!_get_gdata(now, data_rover, data_base))
        return -1;
    t_gtime runEpoch = _data.begin()->epoch();
    _epoch = runEpoch;
    if (_reset_par > 0)
    {
        if (now.sod() % _reset_par == 0)
        {
            _reset_param();
        }
    }
    // save apriory coordinates
    if (_crd_est != CONSTRPAR::FIX)
        _saveApr(runEpoch, _param, _Qx);
    if (!_crd_xml_valid())
        _sig_init_crd = 100.0;
    // select obs for tb log (rover)
    if (_prepareData() < 0)
        return -1;
    if (_isBase)
    {
        _set_rec_info(_gallobj->obj(_site_base)->crd_arp(_epoch), _vBanc(4), _vBanc_base(4));
    }
    _get_initial_value(runEpoch); // for current epoch
    assert(_data.size() >= _minsat);
    if (_data.size() < _minsat)
    {
        if (_spdlog)
            SPDLOG_LOGGER_ERROR(_spdlog, string("gpvtfgo "), ("Not enough visible satellites!"));
        return -1;
    }
    if (_combine_DD() < 0)
    {
        std::cout << "Epoch: " << runEpoch.sow() << " combine DD wrong!!!" << endl;
        clearWindow();
        return -1;
    }
    _optimization();
    // ambiguity resolution
    if (_last_gnss_info->valid)
    {
        _pre_amb_resolution();
        _amb_resolution();
    }
    _marginalization();
    // if (_rover_count == 1)
    //	_initial_prior = false;

    _slide_window();
    return _amb_state ? 1 : 0;
}

void gfgomsf::t_gpvtfgo::clearWindow()
{
    for (int i = 0; i <= _rover_count; i++)
    {
        _headers[i] = 0.0;
        _Pos[i].setZero();
        if (_rover_window[i] != nullptr)
            delete _rover_window[i];
        _rover_window[i] = nullptr;
        _para_window[i].delAllParam();
    }
    if (_last_gnss_info != nullptr)
        delete _last_gnss_info;
    _last_gnss_info = nullptr;
    if (_last_gnss_marginalization_info != nullptr)
        delete _last_gnss_marginalization_info;
    _last_gnss_marginalization_info = nullptr;
    _win_base_data.clear();
    _win_base_data.resize(gwindow_size + 1);
    cur_sat_prn.clear();
    _DD_msg.clear();
    _vDD_msg.clear();
    _amb_manager->clearState();
    _rover_count = -1;
    _global_sat_id = -1;
    _global_amb_id = -1;
    _initial_prior = true;
    _amb_state = false;
    memset(_para_amb, 0, sizeof(_para_amb));
    memset(_para_CRD, 0, sizeof(_para_CRD));
}

bool gfgomsf::t_gpvtfgo::_get_gdata(const t_gtime &now, vector<t_gsatdata> *data_rover, vector<t_gsatdata> *data_base)
{
    if (/*_data.size()*/ _getData(now, data_rover, false) == 0)
    {
        if (_spdlog)
        {
            SPDLOG_LOGGER_DEBUG(_spdlog, _site + now.str_ymdhms(" no observation found at epoch: "));
        }
        return false;
    }
    // apply dcb
    if (_gallbias)
    {
        for (auto &itdata : _data)
        {
            itdata.apply_bias(_gallbias);
        }
    }

    vector<t_gsatdata>::iterator it = _data.begin();
    string double_freq = "";
    string single_freq = "";

    _sat_freqs.clear();
    while (it != _data.end())
    {
        GOBSBAND b1 = _band_index[it->gsys()][FREQ_1];
        GOBSBAND b2 = _band_index[it->gsys()][FREQ_2];

        auto obsL1 = it->select_phase(b1);
        auto obsL2 = it->select_phase(b2);

        if (obsL1 == GOBS::X && obsL2 != GOBS::X || obsL1 != GOBS::X && obsL2 == GOBS::X)
        {
            single_freq += "  " + it->sat();
            _sat_freqs[it->sat()] = "1";
        }

        if (obsL1 != GOBS::X && obsL2 != GOBS::X)
        {
            double_freq += "  " + it->sat();
            _sat_freqs[it->sat()] = "2";
        }

        ++it;
    }

    if (_isBase)
    {

        if (/*_data_base.size()*/ _getData(now, data_base, true) == 0)
        {
            if (_spdlog)
            {
                SPDLOG_LOGGER_DEBUG(_spdlog, _site_base + now.str_ymdhms(" no observation found at epoch: "));
            }
            return false;
        }
        // apply dcb
        if (_gallbias)
        {
            for (auto &itdata_base : _data_base)
            {
                itdata_base.apply_bias(_gallbias);
            }
        }
    }

    if (_gallobj != nullptr)
    {
        auto it_data = _data.begin();
        while (it_data != _data.end())
        {
            string sat_id = it_data->sat();
            shared_ptr<t_gobj> sat_obj = _gallobj->obj(sat_id);

            if (sat_obj == nullptr)
            {
                if (_spdlog)
                {
                    SPDLOG_LOGGER_DEBUG(_spdlog, "remove satellite " + sat_id + " due to missing object");
                }
                it_data = _data.erase(it_data);
            }
            else
            {
                shared_ptr<t_gpcv> sat_pcv = sat_obj->pcv(now);
                if (sat_pcv == nullptr)
                {
                    if (_spdlog)
                    {
                        SPDLOG_LOGGER_DEBUG(_spdlog, "remove satellite " + sat_id + " due to missing PCV data");
                    }
                    it_data = _data.erase(it_data);
                }
                else
                {
                    ++it_data;
                }
            }
        }

        if (_isBase && data_base != nullptr)
        {
            auto it_base = data_base->begin();
            while (it_base != data_base->end())
            {
                string sat_id = it_base->sat();
                shared_ptr<t_gobj> sat_obj = _gallobj->obj(sat_id);

                if (sat_obj == nullptr)
                {
                    it_base = data_base->erase(it_base);
                }
                else
                {
                    shared_ptr<t_gpcv> sat_pcv = sat_obj->pcv(now);
                    if (sat_pcv == nullptr)
                    {
                        it_base = data_base->erase(it_base);
                    }
                    else
                    {
                        ++it_base;
                    }
                }
            }
        }
    }

    return true;
}

void gfgomsf::t_gpvtfgo::_get_initial_value(const t_gtime &runEpoch)
{

    // Predict ambiguity
    // Add or remove ambiguity parameter and appropriate rows/columns covar. matrix
    if (_phase)
    {
        _syncAmb();
    }
    _predictCrd();
    _predictAmb();
    _set_initial_value(runEpoch);
    _initialized = true;
}

void gfgomsf::t_gpvtfgo::_set_initial_value(const t_gtime &runEpoch)
{
    _rover_count++;
    if (_rover_count >= 1)
    {
        // detect time gap
        if (runEpoch.sow() - _headers[_rover_count - 1] > _gtime_interval)
        {
            cout << "[" << _headers[_rover_count - 1] << "]--" << "[" << runEpoch.sow() << "]" << ": time gap !!!]"
                 << endl;
            clearWindow();
            _rover_count++;
        }
    }
    t_gtriple xyz;
    _param.getCrdParam(_site, xyz); //!!!need to check the solution!!!
    _Pos[_rover_count] = Eigen::Vector3d(xyz.crd(0), xyz.crd(1), xyz.crd(2));
    _headers[_rover_count] = runEpoch.sow();
    _rover_window[_rover_count] = new t_grover_msg(runEpoch, _Pos[_rover_count]);
    _win_base_data[_rover_count] = _data_base;
    t_gallpar params_add;
    _gtemp_params(_param, params_add);
    _para_window[_rover_count] = params_add;
    vector<string> cur_sat_list;
    for (auto it : _data)
    {
        cur_sat_list.push_back(it.sat()); // Current list of sats
    }
    map<string, int>::iterator iter_cur_prn = cur_sat_prn.begin(); // Last list of sats
    for (; iter_cur_prn != cur_sat_prn.end();)
    {
        string sat = iter_cur_prn->first;
        auto it_find = find_if(cur_sat_list.begin(), cur_sat_list.end(),
                               [sat](const string &sat_name) { return sat_name == sat; });
        if (it_find == cur_sat_list.end()) // If a satellite present in the previous epoch is absent in the current
                                           // epoch, it indicates a loss of lock.
        {
            iter_cur_prn = cur_sat_prn.erase(iter_cur_prn);
            cout << "sat : " << sat << " lost tracking!!!" << endl;
        }
        else
            iter_cur_prn++;
    }

    for (auto it : _data)
    {
        string sat_name = it.sat();
        shared_ptr<t_gobj> sat_obj2 = _gallobj->obj(sat_name);
        shared_ptr<t_gpcv> sat_pcv2;
        sat_pcv2 = sat_obj2->pcv(runEpoch);
        if (sat_pcv2 == nullptr)
        {
            cur_sat_prn.erase(sat_name);
            continue;
        }
        auto it_find = cur_sat_prn.find(sat_name);

        if (!_amb_manager->is_sat_tracking(sat_name))
        {

            // new sat
            _global_sat_id++;
            cur_sat_prn[sat_name] = _global_sat_id;
            _amb_manager->addNewSat(runEpoch, _rover_count, _global_sat_id, _global_amb_id, it, _param);
        }
        else
        {
            bool is_slip = false;
            vector<GOBSBAND> band = dynamic_cast<t_gsetgnss *>(_set)->band(it.gsys());
            GSYS gs = it.gsys();
            int nf = 5;
            if (band.size())
                nf = band.size();
            for (FREQ_SEQ f = FREQ_1; f <= nf; f = (FREQ_SEQ)(f + 1))
            {
                GOBSBAND b;
                if (f > _frequency)
                    continue; // modified by lvhb in 20201211
                if (band.size() >= f)
                { // automatic dual band selection -> for Anubis purpose
                    b = band[f - 1];
                }
                else
                    b = t_gsys::band_priority(gs, f);

                t_gobs gobs(it.select_phase(b));
                if (true && it.getlli(gobs.gobs()) >= 1)
                {
                    // slip occured
                    cur_sat_prn.erase(sat_name);
                    _global_sat_id++;
                    cur_sat_prn[sat_name] = _global_sat_id;
                    _amb_manager->addNewSat(runEpoch, _rover_count, _global_sat_id, _global_amb_id, it, _param);
                    is_slip = true;
                }
            }
            if (!is_slip)
            {
                //_amb_manager->addRover(sat_name, _rover_count);
                _amb_manager->addRover(_epoch.sow(), sat_name, _rover_count);
            }
        }
    }
    _amb_manager->get_last_epoch_sats(_epoch.sow());
    assert(_amb_manager->cur_sats.size() == _data.size());
}

bool gfgomsf::t_gpvtfgo::_gtemp_params(t_gallpar &params, t_gallpar &params_temp)
{
    params_temp = params;
    t_gpar par_x_base;
    par_x_base.site = _site_base;
    par_x_base.parType = par_type::CRD_X;
    par_x_base.value(_gcrd_base[0]);
    t_gpar par_y_base;
    par_y_base.site = _site_base;
    par_y_base.parType = par_type::CRD_Y;
    par_y_base.value(_gcrd_base[1]);
    t_gpar par_z_base;
    par_z_base.site = _site_base;
    par_z_base.parType = par_type::CRD_Z;
    par_z_base.value(_gcrd_base[2]);
    t_gpar par_clk_rover;
    par_clk_rover.site = _site;
    par_clk_rover.parType = par_type::CLK;
    par_clk_rover.value(_gclk_rover);
    t_gpar par_clk_base;
    par_clk_base.site = _site_base;
    par_clk_base.parType = par_type::CLK;
    par_clk_base.value(_gclk_base);
    params_temp.addParam(par_x_base);
    params_temp.addParam(par_y_base);
    params_temp.addParam(par_z_base);
    params_temp.addParam(par_clk_rover);
    params_temp.addParam(par_clk_base);
    params_temp.reIndex();
    return true;
}

void gfgomsf::t_gpvtfgo::_set_rec_info(const t_gtriple &xyz_base, double clk_rover, double clk_base)
{
    _gcrd_base = xyz_base;
    _gclk_rover = clk_rover;
    _gclk_base = clk_base;
}

int gfgomsf::t_gpvtfgo::_combine_DD()
{
    int flag = -1;
    _select_ref_sat();
    auto it_dd = _DD_msg.begin();
    while (it_dd != _DD_msg.end())
    {
        if (!_get_DD_data(*it_dd, _data_base))
        {
            cout << "[" << it_dd->rover_ref_sat.sat() << "-" << it_dd->rover_nonref_sat.sat() << it_dd->obs_type << ","
                 << it_dd->freq << "]" << "DD msg wrong!" << endl;
            _DD_msg.erase(it_dd);
            continue;
        }

        it_dd++;
    }
    if (!_DD_msg.empty())
    {
        _vDD_msg.push_back(_DD_msg);
        flag = 1;
    }
    return flag;
}

void gfgomsf::t_gpvtfgo::_prior_factor(ceres::Problem &problem)
{
    // prior
    if (_last_gnss_marginalization_info && _last_gnss_marginalization_info->valid)
    {
        MarginalizationGNSSFactor *marginalization_gnss_factor =
            new MarginalizationGNSSFactor(_last_gnss_marginalization_info);
        problem.AddResidualBlock(marginalization_gnss_factor, NULL, _last_gnss_marginalization_para_blocks);
    }
}

void gfgomsf::t_gpvtfgo::_double_to_vector()
{
    all_parameter_block.clear();
    for (int i = 0; i <= _rover_count; i++)
    {
        _Pos[i] = Eigen::Vector3d(_para_CRD[i][0], _para_CRD[i][1], _para_CRD[i][2]);
    }
    for (int i = 0; i < _amb_manager->ambiguity_ids.size(); i++)
    {
        int amb_id = _amb_manager->ambiguity_ids[i];
        _amb_manager->updateAmb(amb_id, _para_amb[amb_id][0]);

        // cout << amb_id <<": "<< _para_amb[amb_id][0] << "   " ;
    }
}

void gfgomsf::t_gpvtfgo::_vector_to_double()
{
    for (int i = 0; i <= _rover_count; i++)
    {
        _para_CRD[i][0] = _Pos[i].x();
        _para_CRD[i][1] = _Pos[i].y();
        _para_CRD[i][2] = _Pos[i].z();
    }
    for (int i = 0; i < _amb_manager->ambiguity_ids.size(); i++)
    {
        int amb_id = _amb_manager->ambiguity_ids[i];
        _para_amb[amb_id][0] = _amb_manager->getAmb(amb_id);
    }
    _amb_manager->generateAmbSearchIndex();
}

void gfgomsf::t_gpvtfgo::_select_ref_sat()
{
    _DD_msg.clear();
    //_obs_index.clear();
    _sat_ref.clear();
    set<string> sysall = dynamic_cast<t_gsetgen *>(_set)->sys();
    bool isSetRefSat = dynamic_cast<t_gsetamb *>(_set)->isSetRefSat();
    bool isPhaseProcess = true;
    for (int obslevel = _obs_level; obslevel <= 3; obslevel++)
    {
        for (auto sys_iter = sysall.begin(); sys_iter != sysall.end(); sys_iter++)
        {
            enum GSYS sys = t_gsys::str2gsys(*sys_iter);
            vector<GOBSBAND> band = dynamic_cast<t_gsetgnss *>(_set)->band(sys);
            int nf = 5;
            if (band.size())
                nf = band.size();
            if (_observ == OBSCOMBIN::IONO_FREE)
                nf = 1;
            FREQ_SEQ f;
            string sat_ref;
            for (FREQ_SEQ freq = FREQ_1; freq <= 2 * nf; freq = (FREQ_SEQ)(freq + 1))
            {
                if (freq <= nf)
                { // phase equations
                    isPhaseProcess = true;
                    f = freq;
                }
                else
                { // code equations
                    isPhaseProcess = false;
                    f = (FREQ_SEQ)(freq - nf);
                }
                if (f > _frequency)
                    continue;
                string sat;
                t_gsatdata obs_sat_ref;
                enum GSYS gs;
                if (!isSetRefSat || (_observ == OBSCOMBIN::RAW_MIX && !isPhaseProcess))
                    sat_ref.clear();
                // sat_ref.empty();
                if (sat_ref.empty())
                {
                    for (auto it = _data.begin(); it != _data.end(); ++it)
                    {
                        std::string sat_id = it->sat(); // Obtain sat ID
                        double elevation = it->ele();   // Obtain the elevation angle
                    }
                    for (auto it = _data.begin(); it != _data.end(); it++)
                    {
                        string sat = it->sat();

                        gs = it->gsys();
                        if (gs == QZS)
                            gs = GPS;
                        if (gs != sys)
                            continue;
                        if ((gs == BDS) && t_gsys::bds_geo(sat))
                            continue;
                        if (/*(_observ == RAW_ALL||_observ == IONO_FREE||_observ == RAW_DOUBLE)
                        //	 &&*/
                            !_reset_amb && !_reset_par && it->islip())
                            continue;

                        GOBSBAND b = _band_index[gs][f];
                        if (isPhaseProcess)
                        {
                            if (!it->band_avail(true).count(b))
                            {
                                continue;
                            }
                        }
                        else
                        {
                            if (!it->band_avail(true).count(b) || !it->band_avail(false).count(b))
                            {
                                continue;
                            }
                        }
                        int base_flag = 0;
                        for (auto it_base = _data_base.begin(); it_base != _data_base.end(); it_base++)
                        {
                            if (it_base->sat() != sat)
                            {
                                continue;
                            }

                            if (isPhaseProcess)
                            {
                                if (it_base->band_avail(true).count(b))
                                {
                                    base_flag = 1;
                                }
                            }
                            else
                            {
                                if (it_base->band_avail(true).count(b) && it_base->band_avail(false).count(b))
                                {
                                    base_flag = 1;
                                }
                            }
                        }
                        if (!base_flag)
                        {
                            continue;
                        }

                        if (sat_ref.empty())
                        {
                            sat_ref = sat;
                            obs_sat_ref = *it;
                            continue;
                        }

                        double e = it->ele_deg();
                        double e2 = obs_sat_ref.ele_deg();
                        if (e >= e2)
                        {
                            sat_ref = sat;
                            sat_ref = it->sat();
                            obs_sat_ref = *it;
                        }

                    } // end select sat_ref
                }
                if (sat_ref.empty())
                    continue;
                if (_ipSatRep[sys] != "" && sat_ref != _ipSatRep[sys])
                {
                    cerr << "refsat bug" << endl;
                    continue;
                }
                if (freq == FREQ_1)
                    _sat_ref.insert(sat_ref);
                t_gsatdata ref_sat_data;
                auto it = find_if(_data.begin(), _data.end(), [sat_ref](t_gsatdata it) { return it.sat() == sat_ref; });
                assert(it != _data.end());
                ref_sat_data = *it;
                for (auto it = _data.begin(); it != _data.end(); it++)
                {
                    sat = it->sat();
                    if (sat == sat_ref)
                    {
                        continue;
                    }
                    gs = it->gsys();
                    if (gs == QZS)
                        gs = GPS;
                    if (gs != sys)
                        continue;
                    GOBSTYPE obstype = TYPE_C;
                    if (isPhaseProcess)
                        obstype = GOBSTYPE::TYPE_L;
                    //_obs_index.push_back(make_pair(it->sat(), make_pair(f, obstype)));
                    // delete unrecorded observations, added by hyChang
                    GOBSBAND b = _band_index[gs][f];
                    if (isPhaseProcess)
                    {
                        if (!it->band_avail(true).count(b))
                        {
                            continue;
                        }
                    }
                    else
                    {
                        if (!it->band_avail(true).count(b) || !it->band_avail(false).count(b))
                        {
                            continue;
                        }
                    }
                    int base_flag = 0;
                    for (auto it_base = _data_base.begin(); it_base != _data_base.end(); it_base++)
                    {
                        if (it_base->sat() != sat)
                        {
                            continue;
                        }

                        if (isPhaseProcess)
                        {
                            if (it_base->band_avail(true).count(b))
                            {
                                base_flag = 1;
                            }
                        }
                        else
                        {
                            if (it_base->band_avail(true).count(b) && it_base->band_avail(false).count(b))
                            {
                                base_flag = 1;
                            }
                        }
                    }
                    if (!base_flag)
                    {
                        continue;
                    }

                    DDEquMsg dd_msg(ref_sat_data, *it, obstype, f);
                    // if (gins_window_size>0)
                    //{
                    //	dd_msg.ref_sat_global_id = cur_sat_prn[ref_sat_data.sat()];
                    //	dd_msg.nonref_sat_global_id = cur_sat_prn[it->sat()];
                    // }
                    // else
                    //{
                    dd_msg.ref_sat_global_id = _amb_manager->get_sat_id(ref_sat_data.sat());
                    // dd_msg.ref_sat_global_id = cur_sat_prn[ref_sat_data.sat()];
                    dd_msg.nonref_sat_global_id = _amb_manager->get_sat_id(it->sat());
                    //}

                    assert(dd_msg.ref_sat_global_id != -1);
                    assert(dd_msg.nonref_sat_global_id != -1);

                    _DD_msg.push_back(dd_msg);

                } // end sat
            } // end f
        } // end sys
    }
}

bool gfgomsf::t_gpvtfgo::_get_DD_data(DDEquMsg &dd_msg, vector<t_gsatdata> base_sat_data)
{
    t_gsatdata rover_ref_sat = dd_msg.rover_ref_sat;
    t_gsatdata rover_nonref_sat = dd_msg.rover_nonref_sat;
    FREQ_SEQ freq = dd_msg.freq;
    GOBSTYPE obstype = dd_msg.obs_type;
    t_gsatdata base_ref_sat, base_nonref_sat;
    if (base_sat_data.empty())
        return false;
    int ref_i = 0;
    int nonref_i = 0;

    for (auto it : base_sat_data)
    {
        if (it.sat() == rover_ref_sat.sat())
        {
            ref_i = 1;
            base_ref_sat = it;
        }
        if (it.sat() == rover_nonref_sat.sat())
        {
            nonref_i = 1;
            base_nonref_sat = it;
        }
    }
    if (!ref_i || !nonref_i)
    {
        cout << "no common view between base and rover!!!" << endl;
        return false;
    }
    map<FREQ_SEQ, GOBSBAND> crt_bands = _band_index[rover_ref_sat.gsys()];
    if (crt_bands.empty())
        return false;
    if (freq > _frequency)
        return false;
    dd_msg.band = crt_bands[freq];
    dd_msg.base_ref_sat = base_ref_sat;
    dd_msg.base_nonref_sat = base_nonref_sat;
    dd_msg.time = rover_ref_sat.epoch();
    dd_msg.base_site = _site_base;
    dd_msg.rover_site = _site;
    return true;
}

bool gfgomsf::t_gpvtfgo::_remove_outlier_sat(const pair<string, int> &outlier)
{
    if (outlier.first != " ")
    {

        if (_amb_manager->cur_sats.size() > _minsat)
        {
            pair<string, int> sat_id = outlier;
            auto it_DD = _vDD_msg[_rover_count].begin();
            int sat_global_id = outlier.second;
            /*
            for (int i = _rover_count; i <= _rover_count; i++)
            {
                for (vector<DDEquMsg>::iterator it_DD = _vDD_msg[i].begin(); it_DD != _vDD_msg[i].end(); )
                {
                    if (it_DD->nonref_sat_global_id == sat_id.second)
                    {
                        it_DD = _vDD_msg[i].erase(it_DD);
                    }
                    else
                    {
                        ++it_DD;
                    }
                }
            }
            */

            while (it_DD != _vDD_msg[_rover_count].end())
            {
                if (it_DD->ref_sat_global_id == sat_global_id || it_DD->nonref_sat_global_id == sat_global_id)
                {
                    _vDD_msg[_rover_count].erase(it_DD);
                    continue;
                }

                it_DD++;
            }

            _amb_manager->removeSat(sat_global_id, _rover_count);
            _amb_manager->get_last_epoch_sats(_epoch.sow()); // update cur sat map

            // update cur sat map
            // auto it = cur_sat_prn.find(sat_id.first);
            // if (it != cur_sat_prn.end())
            //	cur_sat_prn.erase(sat_id.first);
            _last_gnss_info->valid = false;
            auto it = cur_sat_prn.find(sat_id.first);
            if (it != cur_sat_prn.end())
                cur_sat_prn.erase(sat_id.first);
            return true;
        }
        else
            return false;
    }

    return true;
}

bool gfgomsf::t_gpvtfgo::_pre_amb_resolution()
{
    t_gallpar construct_para = _all_para_win;
    int nobs_total, npar_number;
    Matrix A_fgo;
    SymmetricMatrix P_fgo;
    ColumnVector l_fgo, dx_fgo;
    SymmetricMatrix Qx0_fgo, Qx_fgo;
    double vtpv_fgo;
    nobs_total = _last_gnss_info->linearized_jacobians.rows();
    npar_number = _last_gnss_info->linearized_jacobians.cols();
    assert(npar_number == construct_para.parNumber());
    // cout << "nobs_total：" << nobs_total << endl;
    // cout << "npar_number：" << npar_number << endl;
    A_fgo.ReSize(nobs_total, npar_number);
    A_fgo = 0.0;
    P_fgo.ReSize(nobs_total);
    P_fgo = 0.0;
    l_fgo.ReSize(nobs_total);
    l_fgo = 0.0;
    dx_fgo.ReSize(npar_number);
    dx_fgo = 0.0;
    Qx0_fgo.ReSize(npar_number);
    Qx0_fgo = 0.0;
    Qx_fgo.ReSize(npar_number);
    Qx_fgo = 0.0;
    _sig_unit = _last_gnss_info->sig_unit;
    vtpv_fgo = _last_gnss_info->vtpv;
    // for Qx
    for (int i = 0; i < npar_number; i++)
    {
        for (int j = 0; j < npar_number; j++)
        {
            // Qx_fgo(i + 1, j + 1) = _last_gnss_info->Qx(i, j);
            Qx0_fgo(i + 1, j + 1) = _last_gnss_info->Qx(i, j);
        }
    }
    Qx_fgo = Qx0_fgo;
    // for A
    for (int i = 0; i < nobs_total; i++)
    {
        for (int j = 0; j < npar_number; j++)
        {
            A_fgo(i + 1, j + 1) = _last_gnss_info->linearized_jacobians(i, j);
        }
    }
    // for P
    for (int i = 0; i < nobs_total; i++)
    {
        for (int j = 0; j < nobs_total; j++)
        {
            P_fgo(i + 1, j + 1) = _last_gnss_info->weight(i, j);
        }
    }
    // for l
    for (int i = 0; i < nobs_total; i++)
    {
        l_fgo(i + 1) = _last_gnss_info->linearized_residuals(i);
    }
    _filter->add_data(construct_para, dx_fgo, Qx_fgo, _sig_unit, Qx0_fgo);
    _filter->add_data(A_fgo, P_fgo, l_fgo);
    _filter->add_data(vtpv_fgo, nobs_total, npar_number);

    return true;
    ;
}

int gfgomsf::t_gpvtfgo::_optimization()
{
    _removed_sats.clear();
    t_tictoc fgo_gnss;
    int count = 0;
    bool iter_flag = false;
    pair<string, int> outlier = make_pair(" ", -1);
    do
    {
        count++;
        if (!_remove_outlier_sat(outlier))
            break;
        _vector_to_double();
        ceres::Problem problem;
        ceres::LossFunction *loss_function;
        // loss_function = new ceres::HuberLoss(_loss_func_value);
        loss_function = new ceres::CauchyLoss(2);
        for (int i = 0; i < _rover_count + 1; i++)
        {
            problem.AddParameterBlock(_para_CRD[i], 3);
        }

        _prior_factor(problem);
        //_gnss_obs_index.clear();
        assert(_vDD_msg.size() == _rover_count + 1);
        for (int i = 0; i <= _rover_count; i++)
        {
            t_gallpar params_temp(_para_window[i]);
            vector<DDEquMsg> DD_tmp = _vDD_msg[i];
            vector<t_gsatdata> b_sat_data = _win_base_data[i];
            int dd_equ_count = 0;
            for (auto &dd_iter : DD_tmp)
            {
                if (!_get_DD_data(dd_iter, b_sat_data))
                    continue;
                pair<string, string> base_rover_site = make_pair(dd_iter.base_site, dd_iter.rover_site);
                pair<FREQ_SEQ, GOBSBAND> freq_band = make_pair(dd_iter.freq, dd_iter.band);
                vector<pair<t_gsatdata, t_gsatdata>> DD_sat_data;
                DD_sat_data.push_back(make_pair(dd_iter.base_ref_sat, dd_iter.rover_ref_sat));       // reference sat
                DD_sat_data.push_back(make_pair(dd_iter.base_nonref_sat, dd_iter.rover_nonref_sat)); // unreference sat
                GOBSTYPE obstype = dd_iter.obs_type;
                //_gnss_obs_index.push_back(make_pair(make_pair(dd_iter.rover_nonref_sat.sat(),
                //dd_iter.nonref_sat_global_id), make_pair(dd_iter.freq, obstype)));
                if (obstype == GOBSTYPE::TYPE_C)
                {
                    PseudorangeDDFactor *pf =
                        new PseudorangeDDFactor(dd_iter.time, base_rover_site, params_temp, DD_sat_data, _gbias_model,
                                                freq_band); // DD_Pseudorange
                    problem.AddResidualBlock(pf, NULL, _para_CRD[i]);
                }
                if (obstype == GOBSTYPE::TYPE_L)
                {
                    int id1, id2;
                    id1 = _amb_manager->getAmbSearchIndex(make_pair(dd_iter.ref_sat_global_id, dd_iter.freq));
                    id2 = _amb_manager->getAmbSearchIndex(make_pair(dd_iter.nonref_sat_global_id, dd_iter.freq));

                    if (id1 == -1 || id2 == -1)
                        continue;

                    CarrierphaseDDFactor *lf =
                        new CarrierphaseDDFactor(dd_iter.time, base_rover_site, params_temp, DD_sat_data, _gbias_model,
                                                 freq_band); // DD_Carrierphase
                    problem.AddResidualBlock(lf, NULL, _para_CRD[i], _para_amb[id1], _para_amb[id2]);
                }
                dd_equ_count++;
            }
            assert(dd_equ_count == DD_tmp.size());
        }
        // if (_initial_prior)

        if (1)
        {
            for (int i = 0; i < _amb_manager->ambiguity_ids.size(); i++)
            {
                int amb_id = _amb_manager->ambiguity_ids[i];

                double initial_amb = _amb_manager->getInitialAmb(amb_id);
                InitialGnssAMB *amb_prior = new InitialGnssAMB(initial_amb);
                problem.AddResidualBlock(amb_prior, NULL, _para_amb[amb_id]);
            }
        }
        for (int i = 0; i < _rover_count + 1; i++)
        {

            InitialPosFactor *initial_pos =
                new InitialPosFactor(Eigen::Vector3d(_para_CRD[i])); // initial position factor
            // initial_pos->sqrt_info = 1e-4 * Eigen::Matrix<double, 3, 3>::Identity();
            problem.AddResidualBlock(initial_pos, NULL, _para_CRD[i]);
        }
        // ceres solver
        ceres::Solver::Options options;
        options.linear_solver_type = ceres::DENSE_QR;
        options.max_num_iterations = 10;
        options.trust_region_strategy_type = ceres::DOGLEG;
        ceres::Solver::Summary summary;
        ceres::Solve(options, &problem, &summary);

        _posteriori_test(problem);
        double idx = _gobs_outlier_detection(outlier);
        if (idx >= 0)
            iter_flag = true;
        else
            iter_flag = false;

    } while (iter_flag);

    std::cout << "Epoch: " << _headers[_rover_count] << " time cost (ms): " << fgo_gnss.toc()
              << " iter times: " << count << endl;
    if (_last_gnss_info->valid)
    {
        std::cout << "cur epoch: " << "sigma: " << _last_gnss_info->sig_unit << "  vtpv: " << _last_gnss_info->vtpv
                  << std::endl;
    }
    else
    {
        cout << "soving failed in this epoch!!!" << endl;
    }
    _double_to_vector();
    return 1;
}

void gfgomsf::t_gpvtfgo::_slide_window()
{
    if (_rover_count == gwindow_size)
    {
        for (int i = 0; i < gwindow_size; i++)
        {
            _headers[i] = _headers[i + 1];
            _Pos[i].swap(_Pos[i + 1]);
            std::swap(_rover_window[i], _rover_window[i + 1]);
            _win_base_data[i].swap(_win_base_data[i + 1]);
            _para_window[i] = _para_window[i + 1];
        }
        _vDD_msg.erase(_vDD_msg.begin());
        delete _rover_window[gwindow_size];
        _rover_window[gwindow_size] = nullptr;
        _amb_manager->slidingWindow();
        _rover_count--;
    }

    _amb_manager->get_last_epoch_sats(_epoch.sow());
}

void gfgomsf::t_gpvtfgo::_marginalization()
{

    // The residual factors added during marginalization include: the marginalization factor retained from the previous
    // iteration, the absolute observation factor of the current dual-difference pseudorange and phase, and the initial
    // ambiguity constraint factor.
    if (_rover_count == gwindow_size)
    {
        _vector_to_double();
        ceres::LossFunction *loss_function;
        loss_function = new ceres::CauchyLoss(2);
        // loss_function = new  ceres::CauchyLoss(_loss_func_value);
        GNSSInfo *gnss_marginalization_info = new GNSSInfo();
        if (_last_gnss_marginalization_info && _last_gnss_marginalization_info->valid)
        {
            vector<int> drop_set;
            vector<int> amb_margin = _amb_manager->getMarginAmb();
            // preevious para_blocks
            for (int i = 0; i < static_cast<int>(_last_gnss_marginalization_para_blocks.size()); i++)
            {
                for (int j = 0; j < amb_margin.size(); j++)
                {
                    if (_last_gnss_marginalization_para_blocks[i] == _para_amb[amb_margin[j]])
                        drop_set.push_back(i);
                }
            }
            // construct new marginlization_factor
            MarginalizationGNSSFactor *marginalization_gnss_factor =
                new MarginalizationGNSSFactor(_last_gnss_marginalization_info);
            GNSSResidualBlockInfo *residual_block_info = new GNSSResidualBlockInfo(
                marginalization_gnss_factor, NULL, _last_gnss_marginalization_para_blocks, drop_set);
            gnss_marginalization_info->addResidualBlockInfo(residual_block_info);
        }
        /*
        for (int i = 0; i <= _rover_count; i++)
        {
            if (i == 0)
            {
                t_gallpar params_temp = _para_window[i];
                vector<DDEquMsg> dd_msg = _vDD_msg[i];
                vector<t_gsatdata> b_sat_data = _win_base_data[i];
                for (auto &dd_iter : dd_msg)
                {
                    if (!_get_DD_data(dd_iter, b_sat_data)) continue;
                    pair<string, string> base_rover_site = make_pair(dd_iter.base_site, dd_iter.rover_site);
                    pair<FREQ_SEQ, GOBSBAND> freq_band = make_pair(dd_iter.freq, dd_iter.band);
                    vector<pair<t_gsatdata, t_gsatdata>> DD_sat_data;
                    DD_sat_data.push_back(make_pair(dd_iter.base_ref_sat, dd_iter.rover_ref_sat));
                    DD_sat_data.push_back(make_pair(dd_iter.base_nonref_sat, dd_iter.rover_nonref_sat));
                    GOBSTYPE  obstype = dd_iter.obs_type;
                    if (obstype == GOBSTYPE::TYPE_C)
                    {
                        PseudorangeDDFactor *pf = new PseudorangeDDFactor(dd_iter.time, base_rover_site, params_temp,
        DD_sat_data, _gbias_model, freq_band); GNSSResidualBlockInfo *residual_block_info = new
        GNSSResidualBlockInfo(pf, NULL, vector<double *> {_para_CRD[i]}, vector<int>{0});
                        gnss_marginalization_info->addResidualBlockInfo(residual_block_info);
                    }
                    if (obstype == GOBSTYPE::TYPE_L)
                    {
                        int id1, id2;
                        id1 = _amb_manager->getAmbSearchIndex(make_pair(dd_iter.ref_sat_global_id, dd_iter.freq));
                        id2 = _amb_manager->getAmbSearchIndex(make_pair(dd_iter.nonref_sat_global_id, dd_iter.freq));

                        if (id1 == -1 || id2 == -1)
                            continue;

                        vector<int> drop_set{ 0 };
                        if (_amb_manager->getAmbStartRoverID(id1) == 0 && _amb_manager->getAmbStartRoverID(id2) == 0)
                        {
                            if (_amb_manager->getAmbEndRoverID(id1) == 0)
                                drop_set.push_back(1);
                            if (_amb_manager->getAmbEndRoverID(id2) == 0)
                                drop_set.push_back(2);

                            CarrierphaseDDFactor *lf = new CarrierphaseDDFactor(dd_iter.time, base_rover_site,
        params_temp, DD_sat_data, _gbias_model, freq_band); GNSSResidualBlockInfo *residual_block_info = new
        GNSSResidualBlockInfo(lf, NULL, vector<double *> {_para_CRD[i], _para_amb[id1], _para_amb[id2]}, drop_set);
                            gnss_marginalization_info->addResidualBlockInfo(residual_block_info);
                        }

                    }
                }
            }
        }	*/

        for (int i = 0; i <= _rover_count; i++)
        {
            if (i == 0)
            {
                t_gallpar params_temp = _para_window[i];
                for (auto &dd_iter : _vDD_msg[i])
                {
                    pair<string, string> base_rover_site = make_pair(dd_iter.base_site, dd_iter.rover_site);
                    pair<FREQ_SEQ, GOBSBAND> freq_band = make_pair(dd_iter.freq, dd_iter.band);
                    vector<pair<t_gsatdata, t_gsatdata>> DD_sat_data;
                    DD_sat_data.push_back(make_pair(dd_iter.base_ref_sat, dd_iter.rover_ref_sat));
                    DD_sat_data.push_back(make_pair(dd_iter.base_nonref_sat, dd_iter.rover_nonref_sat));
                    GOBSTYPE obstype = dd_iter.obs_type;
                    if (obstype == GOBSTYPE::TYPE_C)
                    {
                        PseudorangeDDFactor *pf = new PseudorangeDDFactor(dd_iter.time, base_rover_site, params_temp,
                                                                          DD_sat_data, _gbias_model, freq_band);
                        GNSSResidualBlockInfo *residual_block_info = new GNSSResidualBlockInfo(
                            pf, loss_function, vector<double *>{_para_CRD[i]}, vector<int>{0});
                        gnss_marginalization_info->addResidualBlockInfo(residual_block_info);
                    }
                    if (obstype == GOBSTYPE::TYPE_L)
                    {
                        int id1, id2;
                        id1 = _amb_manager->getAmbSearchIndex(make_pair(dd_iter.ref_sat_global_id, dd_iter.freq));
                        id2 = _amb_manager->getAmbSearchIndex(make_pair(dd_iter.nonref_sat_global_id, dd_iter.freq));

                        if (id1 == -1 || id2 == -1)
                            continue;

                        vector<int> drop_set{0};
                        if (_amb_manager->getAmbStartRoverID(id1) == 0 && _amb_manager->getAmbStartRoverID(id2) == 0)
                        {
                            if (_amb_manager->getAmbEndRoverID(id1) == 0)
                                drop_set.push_back(1);
                            if (_amb_manager->getAmbEndRoverID(id2) == 0)
                                drop_set.push_back(2);

                            CarrierphaseDDFactor *lf = new CarrierphaseDDFactor(
                                dd_iter.time, base_rover_site, params_temp, DD_sat_data, _gbias_model, freq_band);
                            GNSSResidualBlockInfo *residual_block_info = new GNSSResidualBlockInfo(
                                lf, NULL, vector<double *>{_para_CRD[i], _para_amb[id1], _para_amb[id2]}, drop_set);
                            gnss_marginalization_info->addResidualBlockInfo(residual_block_info);
                        }
                    }
                }
            }
        }

        gnss_marginalization_info->preMarginalize();
        gnss_marginalization_info->marginalize();

        std::unordered_map<long, double *> addr_shift;

        vector<int> cur_amb = _amb_manager->getCurWinAmb();
        for (int i = 0; i < cur_amb.size(); i++)
        {
            addr_shift[reinterpret_cast<long>(_para_amb[cur_amb[i]])] = _para_amb[cur_amb[i]];
        }

        vector<double *> parameter_blocks = gnss_marginalization_info->getParameterBlocks(addr_shift);
        if (_last_gnss_marginalization_info)
            delete _last_gnss_marginalization_info;
        // save marginalization info
        _last_gnss_marginalization_info = gnss_marginalization_info;
        _last_gnss_marginalization_para_blocks = parameter_blocks;
    }
}

int gfgomsf::t_gpvtfgo::_gobs_outlier_detection(pair<string, int> &outlier)
{
    pair<string, int> sat_id;
    int idx = -1;
    if (_last_gnss_info->valid)
    {
        Eigen::VectorXd v = _last_gnss_info->v_norm;
        // std::cout << v << endl;
        int nobs = v.rows();
        double max = 0.0;

        for (int i = 0; i < nobs; i++)
        {
            if (fabs(v(i)) > max && fabs(v(i)) > _max_res_norm)
            {
                max = fabs(v(i));
                idx = i;
            }
        }
        if (idx >= 0)
        {
            sat_id = _gnss_obs_index[idx].first;
            outlier = sat_id;
            int id = sat_id.second;
            auto it_find = find_if(_removed_sats.begin(), _removed_sats.end(), [id](pair<string, int> &sat_id) {
                return sat_id.second == id;
            });
            if (it_find == _removed_sats.end())
                _removed_sats.push_back(sat_id);

            string obsType = gobstype2str(_gnss_obs_index[idx].second.second);
            ostringstream os;
            /*os << _site << " outlier (" << obsType << _obs_index[idx].second.first << ") " << sat
                << " v: " << fixed << setw(16) << right << setprecision(3) << max;
            if (_log)
                _log->comment(1, "gpvtfgo", _epoch.str_ymdhms(" epoch ") + os.str());*/

            std::cout << _site << " outlier (" << obsType << _gnss_obs_index[idx].second.first << ") " << sat_id.first
                      << " v: " << fixed << setw(16) << right << setprecision(3) << max << endl;

            if (_spdlog)
                SPDLOG_LOGGER_ERROR(_spdlog, string("gpvtfgo "), _epoch.str_ymdhms(" epoch ") + os.str());
        }

        if (_vDD_msg[_rover_count].size() - _removed_sats.size() <= 2)
            idx = -1;
    }

    return idx;
}

bool gfgomsf::t_gpvtfgo::_check_outlier(const string &sat)
{
    auto it_find = find_if(outlier_sats.begin(), outlier_sats.end(), [sat](const string &sat_name) {
        return sat_name == sat;
    });
    if (it_find != outlier_sats.end())
        return true;
    else
        return false;
}

void gfgomsf::t_gpvtfgo::_posteriori_test(ceres::Problem &problem)
{
    cout << "begin to gnss posteriori_test" << endl;
    // for current gnss epoch
    _all_para_win.delAllParam();
    ceres::LossFunction *loss_function;
    // loss_function = new ceres::HuberLoss(_loss_func_value);
    loss_function = new ceres::CauchyLoss(_loss_func_value);
    GNSSInfo *gnss_info = new GNSSInfo();
    if (_vDD_msg[_rover_count].size() > 2)
    {
        // constrcut window para_index
        vector<vector<int>> crd_para_col_index;
        vector<vector<int>> amb_para_col_index;
        vector<double *> _parameter_blocks;
        int total_para_size = 0;
        vector<int> crdi;
        vector<par_type> partype{par_type::CRD_X, par_type::CRD_Y, par_type::CRD_Z};
        _parameter_blocks.push_back(_para_CRD[_rover_count]);
        for (int j = 0; j < 3; j++)
        {
            t_gpar par_crd;
            par_crd.site = _site;
            par_crd.parType = partype[j];
            par_crd.value(_para_CRD[_rover_count][j]);
            par_crd.beg = _rover_window[_rover_count]->cur_time;
            par_crd.end = _rover_window[_rover_count]->cur_time;
            par_crd.index = j + 1;
            _all_para_win.addParam(par_crd);
            int id = j;
            crdi.push_back(id);
            total_para_size = total_para_size + 1;
        }
        crd_para_col_index.push_back(crdi);

        int amb_index_start = 3;
        map<int, int> amb_col_id;
        int amb_size = -1;
        _gnss_obs_index.clear();
        t_gallpar params_temp = _para_window[_rover_count];
        vector<DDEquMsg> dd_msg = _vDD_msg[_rover_count];
        vector<t_gsatdata> b_sat_data = _win_base_data[_rover_count];
        for (auto &dd_iter : dd_msg)
        {
            if (!_get_DD_data(dd_iter, b_sat_data))
                continue;
            pair<string, string> base_rover_site = make_pair(dd_iter.base_site, dd_iter.rover_site);
            pair<FREQ_SEQ, GOBSBAND> freq_band = make_pair(dd_iter.freq, dd_iter.band);
            vector<pair<t_gsatdata, t_gsatdata>> DD_sat_data;
            DD_sat_data.push_back(make_pair(dd_iter.base_ref_sat, dd_iter.rover_ref_sat));
            DD_sat_data.push_back(make_pair(dd_iter.base_nonref_sat, dd_iter.rover_nonref_sat));
            GOBSTYPE obstype = dd_iter.obs_type;
            _gnss_obs_index.push_back(make_pair(make_pair(dd_iter.rover_nonref_sat.sat(), dd_iter.nonref_sat_global_id),
                                                make_pair(dd_iter.freq, obstype)));
            if (obstype == GOBSTYPE::TYPE_C)
            {
                PseudorangeDDFactor *pf = new PseudorangeDDFactor(dd_iter.time, base_rover_site, params_temp,
                                                                  DD_sat_data, _gbias_model, freq_band);
                GNSSResidualBlockInfo *residual_block_info =
                    new GNSSResidualBlockInfo(pf, NULL, vector<double *>{_para_CRD[_rover_count]});
                map<long, vector<int>> para_index;
                para_index[reinterpret_cast<long>(_para_CRD[_rover_count])] = crd_para_col_index[0];
                gnss_info->addResidualBlockInfo(residual_block_info, para_index);
            }
            if (obstype == GOBSTYPE::TYPE_L)
            {
                vector<int> amb_id12(2);
                amb_id12[0] = _amb_manager->getAmbSearchIndex(make_pair(dd_iter.ref_sat_global_id, dd_iter.freq));
                amb_id12[1] = _amb_manager->getAmbSearchIndex(make_pair(dd_iter.nonref_sat_global_id, dd_iter.freq));

                if (amb_id12[0] == -1 || amb_id12[1] == -1)
                    continue;

                for (int i = 0; i < amb_id12.size(); i++)
                {
                    int amb_id = amb_id12[i];
                    if (amb_col_id.find(amb_id) == amb_col_id.end())
                    {
                        amb_size++;
                        int amb_para_id = amb_index_start + amb_size;
                        amb_col_id[amb_id] = amb_para_id;
                        amb_para_col_index.push_back(vector<int>{amb_para_id});
                        _amb_manager->addGpara(_all_para_win, amb_id, _para_amb[amb_id][0]);
                        _parameter_blocks.push_back(_para_amb[amb_id]);
                    }
                }
                map<long, vector<int>> para_index;
                // cout << "para_col_index: " << para_col_index.size() << " id1: " << _rover_count + 1+ id1 << " id2: "
                // << _rover_count + 1+id2 << endl;
                double *addr1 = _para_amb[amb_id12[0]];
                double *addr2 = _para_amb[amb_id12[1]];
                para_index[reinterpret_cast<long>(_para_CRD[_rover_count])] = crd_para_col_index[0];
                para_index[reinterpret_cast<long>(addr1)] = vector<int>{amb_col_id[amb_id12[0]]};
                para_index[reinterpret_cast<long>(addr2)] = vector<int>{amb_col_id[amb_id12[1]]};
                CarrierphaseDDFactor *lf = new CarrierphaseDDFactor(dd_iter.time, base_rover_site, params_temp,
                                                                    DD_sat_data, _gbias_model, freq_band);
                GNSSResidualBlockInfo *residual_block_info = new GNSSResidualBlockInfo(
                    lf, NULL,
                    vector<double *>{_para_CRD[_rover_count], _para_amb[amb_id12[0]], _para_amb[amb_id12[1]]});
                gnss_info->addResidualBlockInfo(residual_block_info, para_index);
            }
        }
        total_para_size = total_para_size + amb_size + 1;
        // construct variances by ceres solver
        ceres::Covariance::Options options_co;
        // options_co.algorithm_type = ceres::SPARSE_QR;
        // options_co.algorithm_type = ceres::DENSE_SVD;
        // options_co.sparse_linear_algebra_library_type = ceres::SparseLinearAlgebraLibraryType::SUITE_SPARSE;
        // options_co.apply_loss_function = false; //optional, true or false is depended on the reliability of
        // covariance
        ceres::Covariance covariance(options_co);
        std::vector<const double *> covariance_blocks; // all parameter_blocks
        for (int i = 0; i < _parameter_blocks.size(); i++)
        {
            covariance_blocks.push_back(_parameter_blocks[i]);
        }
        Eigen::MatrixXd Qx = Eigen::MatrixXd::Zero(total_para_size, total_para_size);
        try
        {
            covariance.Compute(covariance_blocks, &problem);
            covariance.GetCovarianceMatrix(covariance_blocks, Qx.data());
        }
        catch (...)
        {
            cout << "Covariance Compute Failed" << endl;
        }
        gnss_info->constructEqu_fromCeres(Qx);
    }
    if (_last_gnss_info)
        delete _last_gnss_info;
    _last_gnss_info = gnss_info;
}

gfgomsf::t_gpvtfgo::DDEquMsg::DDEquMsg(const t_gsatdata &_ref_sat, const t_gsatdata &_nonref_sat,
                                       const GOBSTYPE &_obs_type, const FREQ_SEQ &_freq)
    : rover_ref_sat(_ref_sat), rover_nonref_sat(_nonref_sat), obs_type(_obs_type), freq(_freq), ref_sat(_ref_sat.sat()),
      nonref_sat(_nonref_sat.sat())
{
    this->time = _nonref_sat.epoch();
}

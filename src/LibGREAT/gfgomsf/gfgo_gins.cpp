/**
 * @file         gfgo_gins.cpp
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        Class for factor graph-based GNSS/SINS integrated solution
 * @version      1.0
 * @date         2025-11-04
 *
 * @copyright Copyright (c) 2025, Wuhan University. All rights reserved.
 *
 */

#include "gfgo_gins.h"
#include "gfactor/carrierphase_DD_integration_factor.h"
#include "gfactor/ginitial_bias_factor.h"
#include "gfactor/ginitial_pose_factor.h"
#include "gfactor/pseudorange_DD_integration_factor.h"
#include "gfgo/gutility.h"
#include "gset/gsetsensors.h"
using namespace gfgomsf;

gfgomsf::t_gfgo_gins::t_gfgo_gins(string site, string site_base, t_gsetbase *gset,
                                  std::shared_ptr<spdlog::logger> spdlog, t_gallproc *allproc)
    : t_gspp(site, gset, spdlog), t_gpvtfgo(site, site_base, gset, spdlog, allproc),
      t_gpvtflt(site, site_base, gset, spdlog, allproc), t_gsinskf(gset, spdlog, site), t_gfgo_para(gset), t_gfgo(gset),
      t_gpara(gset)
{
    _solver_flag = INITIAL;
    _c_gnss_factor = false;
    _opt_valid = 0;
    _imudata = nullptr;
    _fgoflag = 1;
    _gtime_interval = 10.0;
    for (int i = 0; i <= WINDOW_SIZE + 1; i++)
        _pre_integrations[i] = nullptr;
    _cntrep = 0;
    _publish_result.Initialize();
}

int gfgomsf::t_gfgo_gins::ProcessBatch(const t_gtime &beg, const t_gtime &end, bool beg_end)
{
#ifdef BMUTEX
    boost::mutex::scoped_lock lock(_mutex);
#endif
    _gmutex.lock();
    if (_grec == nullptr)
    {
        ostringstream os;
        os << "ERROR: No object found (" << _site << "). Processing terminated!!! " << beg.str_ymdhms() << " -> "
           << end.str_ymdhms() << endl;
        if (_spdlog)
            SPDLOG_LOGGER_ERROR(_spdlog, string("gintegration:  ") + os.str());
        _gmutex.unlock();
        return -1;
    }
    int sign = 1;
    double subint = 0.1;
    bool prtOut = true;
    _beg_end = beg_end;
    t_gpvtflt::InitProc(beg, end);
    if (!this->_init())
        return -1;
    // need execute frequently in real-time mode
    _gobs->setepoches(_site);
    if (_beg_end && _imu_beg < _gnss_beg)
    {
        _imu_crt = _imudata->erase_bef(_gnss_beg, _beg_end);
    }
    else
    {
        _gmutex.unlock();
        if (_beg_end)
            t_gpvtflt::processBatch(_gnss_crt, _imu_crt, false);
        if (_beg_end && _gnss_crt < _imu_crt)
        {
            int nEpo = round(_imu_crt.diff(_gnss_crt) / (_fgoflag * _sampling) + 0.5);
            if (_sampling > 1)
            {
                _gnss_crt.add_secs(int(_fgoflag * _sampling * nEpo)); //  < 1Hz data
            }
            else
            {
                _gnss_crt.add_dsec(_fgoflag * _sampling * nEpo); //  >=1Hz data
            }
        }
    }

    std::cerr << "\n"
              << _site << ": Start GNSS/SINS Integration Processing:" << _imu_beg.str_ymdhms() << " "
              << _imu_end.str_ymdhms() << endl;
    std::cout << "EPOCH " << _imu_crt.str_ymdhms() << endl;
    t_gposdata::data_pos posdata;
    bool time_loop = true;
    while (time_loop)
    {
        if (_beg_end && (_imu_crt > _imu_end || _gnss_crt > _gnss_end))
        {
            cerr << " Processing Finished!" << endl;
            if (_spdlog)
                SPDLOG_LOGGER_INFO(_spdlog, string("gintegration:  ") + "Processing Finished!");
            time_loop = false;
            break;
        }
        if (!_imudata->load(_wm, _vm, _shm.t, _shm.ts, _shm.nSamples, _beg_end))
            break;
        _imu_crt = t_gtime(_gnss_crt.gwk(), _shm.t);
        if (!_aligned)
        {
            if (t_gfgo_gins::_time_valid(_gnss_crt, _imu_crt))
                Flag = _msf_flag = t_gfgo_gins::_getPOS(posdata);
            _aligned = cascaded_align(posdata.pos, posdata.vn);
        }
        else
        {
            sins.Update(_wm, _vm, _shm);
            _gravity_update();
            _input_imu_data(_shm.t, _wm, _vm); // input imu data
            int imeas = _getMeas(), irc = 0;
            set<MEAS_TYPE>::const_iterator it = _Meas_Type.begin();
            if (_Meas_Type.size() > 0)
                std::cout << endl << "imu_time:" << setprecision(10) << _shm.t << endl;
            _opt_valid = 0;

            while (it != _Meas_Type.end())
            {
                switch (*it)
                {
                case GNSS_MEAS: {
                    _cur_node_time = _gnss_crt.sow();
                    _gnss_processing(); // process GNSS data
                    if (_c_gnss_factor)
                    {
                        _set_frame_pose(_rover_count);
                        _imu_pre_integration(_last_pre_integration_time, _cur_node_time); // IMU pre_integration
                        _last_pre_integration_time = _cur_node_time;
                    }
                    if (_msf_type == MSF_TYPE::GINS_TC_MODE)
                        _gins_processing(); // GNSS and IMU Fusion Optimization
                }
                break;
                default:
                    break;
                }
                ++it;
            }

            while (abs(_imu_crt.diff(_gnss_crt)) > _sampling)
            {
                if (_sampling > 1)
                    _gnss_crt.add_secs(_fgoflag * int(_sampling)); // =<1Hz data
                else
                    _gnss_crt.add_dsec(_fgoflag * _sampling); //  >1Hz data
            }
            if (_opt_flag)
            {
                _gins_feedback();
                cerr << "\r" << "Current Epoch: " << (_imu_crt).str_ymdhms();
            }
            if (fabs(sins.t - int(sins.t)) < _shm.delay)
                this->write();
            if (!_initial_gnss_pos)
            {
                _imu_state.orientation = Eigen::Quaterniond(t_gbase::q2mat(sins.qnb));
                _imu_state.position = Geod2Cart(sins.pos, false);
                _publish_result.UpdateNewState(_imu_state);
            }
            _opt_flag = false;
        }
    }
    _gmutex.unlock();
    return 0;
}
void gfgomsf::t_gfgo_gins::_input_imu_data(const double &t, const vector<Eigen::Vector3d> &wm,
                                           const vector<Eigen::Vector3d> &vm)
{
    assert(wm.size() == vm.size());
    /* only one sample can be processed successfully */
    assert(wm.size() == 1);
    for (int i = 0; i < wm.size(); i++)
    {
        _accBuf.push(make_pair(t, vm[i] / _imu_ts));
        _gyrBuf.push(make_pair(t, wm[i] / _imu_ts));
    }
}

int gfgomsf::t_gfgo_gins::_init()
{

    if (!_ins_init())
    {

        if (_spdlog)
            SPDLOG_LOGGER_ERROR(_spdlog, string("gintegration:  ") + "initialize integration filter failed!!! ");
        return -1;
    }
    _imu_crt = _gnss_crt = t_gtime::current_time(t_gtime::GPS);
    _gnss_beg = _beg_time;
    _gnss_end = _end_time;

    if (_imudata)
    {
        double start = _imudata->beg_obs(_beg_end);
        double end = _imudata->end_obs(_beg_end);
        double start_set = dynamic_cast<t_gsetins *>(_set)->start();
        double end_set = dynamic_cast<t_gsetins *>(_set)->end();
        start = (start < start_set) ? start_set : start;
        end = (end < end_set) ? end : end_set;
        if (_beg_end)
        {
            _gnss_crt = _gobs->beg_obs(_site);
            _imu_beg = t_gtime(_gnss_beg.gwk(), start);
            _imu_end = t_gtime(_gnss_end.gwk(), end);
        }
        else
        {
            _gnss_crt = _gobs->end_obs(_site);
            _imu_beg = t_gtime(_gnss_beg.gwk(), end);
            _imu_end = t_gtime(_gnss_end.gwk(), start);
        }
        _imu_crt = _imu_beg;
    }
    else
    {
        /*if (_log) _log->comment(1, "t_gfgo_gins", "IMU observation is not existing!!! ");*/

        if (_spdlog)
            SPDLOG_LOGGER_ERROR(_spdlog, string("t_gfgo_gins "), "IMU observation is not existing!!! ");
        return -1;
    }

    return 1;
}

void gfgomsf::t_gfgo_gins::_set_frame_pose(const int &framecount)
{
    if (_initial_gnss_pos)
    {
        _initial_gnss_pos = false;
        _initial_pos = sins.pos_ecef;
    }

    Eigen::Quaterniond e_q = Eigen::Quaterniond(sins.qeb.q0, sins.qeb.q1, sins.qeb.q2, sins.qeb.q3);
    e_q.normalized();
    assert(framecount <= WINDOW_SIZE);
    t_gfgo_gins::_headers[framecount] = _cur_node_time;
    _Rs[framecount] = e_q.toRotationMatrix();
    _Ps[framecount] = sins.pos_ecef;
    _Vs[framecount] = sins.ve;
    _Bas[framecount] = sins.db;
    _Bgs[framecount] = sins.eb;
}

bool gfgomsf::t_gfgo_gins::_time_valid(t_gtime gt, t_gtime inst)
{
    //_mutex.lock();
    _amb_state = false;
    bool res_valid = false;
    double crt = inst.sow() + inst.dsec();
    _data.erase(_data.begin(), _data.end());
    t_gtime obsEpo = _gobs->load(_site, crt);

    if ((abs(inst.diff(obsEpo)) < _shm.delay && inst >= obsEpo) || abs(inst.diff(obsEpo)) < 1e-6)
    {

        // detect cycle slip by lvhb
        _gpre->ProcessBatch(_site, obsEpo, obsEpo, _sampling, false);
        if (_isBase)
        {
            _gpre->ProcessBatch(_site_base, obsEpo, obsEpo, _sampling, false);
        }

        /*if (!_isClient)*/ _data = _gobs->obs(_site, obsEpo);

        if (_data.size() > 0)
        {
            res_valid = true;
            // apply dcb
            if (_gallbias)
            {
                for (auto &itdata : _data)
                {
                    itdata.apply_bias(_gallbias);
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
                        shared_ptr<t_gpcv> sat_pcv = sat_obj->pcv(obsEpo);
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
            }
        }
        else
        {

            if (_spdlog)
                SPDLOG_LOGGER_ERROR(_spdlog, string("gintegration "),
                                    _site + gt.str_ymdhms(" no observation found at epoch: "));

            res_valid = false;
        }

        if (_isBase)
        {
            _data_base.erase(_data_base.begin(), _data_base.end());
            _data_base = _gobs->obs(_site_base, obsEpo);

            if (_data_base.size() > 0)
            {
                res_valid = true;
                // apply dcb
                if (_gallbias)
                {
                    for (auto &itdata : _data_base)
                    {
                        itdata.apply_bias(_gallbias);
                    }
                }

                if (_gallobj != nullptr)
                {
                    auto it_base = _data_base.begin();
                    while (it_base != _data_base.end())
                    {
                        string sat_id = it_base->sat();
                        shared_ptr<t_gobj> sat_obj = _gallobj->obj(sat_id);

                        if (sat_obj == nullptr)
                        {
                            it_base = _data_base.erase(it_base);
                        }
                        else
                        {
                            shared_ptr<t_gpcv> sat_pcv = sat_obj->pcv(obsEpo);
                            if (sat_pcv == nullptr)
                            {
                                if (_spdlog)
                                {
                                    SPDLOG_LOGGER_DEBUG(_spdlog,
                                                        "remove base satellite " + sat_id + " due to missing PCV data");
                                }
                                it_base = _data_base.erase(it_base);
                            }
                            else
                            {
                                ++it_base;
                            }
                        }
                    }
                }
            }
            else
            {

                if (_spdlog)
                    SPDLOG_LOGGER_ERROR(_spdlog, string("gintegration "),
                                        _site_base + gt.str_ymdhms(" no base observation found at epoch: "));
                res_valid = false;
            }
        }
    }

    if (_data.size() == 0)
    {
        res_valid = false;
    }

    if (res_valid)
        _gnss_crt = obsEpo;
    //_mutex.unlock();
    return res_valid;
}

int gfgomsf::t_gfgo_gins::_getMeas()
{
#ifdef BMUTEX
    boost::mutex::scoped_lock lock(_mutex);
#endif

    _Meas_Type.clear();
    _msf_flag = NO_MEAS;
    double t = _imu_crt.sow() + _imu_crt.dsec();

    MEAS_TYPE meas_type = meas_state();
    if (double_eq(fabs(t - int(t)), 0.0))
    {
        if (meas_type == ZUPT_MEAS)
        {
            Eigen::Quaterniond e_q = Eigen::Quaterniond(sins.qeb.q0, sins.qeb.q1, sins.qeb.q2, sins.qeb.q3);
            _map_motion.insert(make_pair(_imu_crt, MOTION_TYPE::m_static));
            _map_pose.insert(make_pair(_imu_crt, make_pair(sins.pos_ecef, e_q)));
        }
        else if (meas_type == NHC_MEAS)
        {
            _map_motion.insert(make_pair(_imu_crt, MOTION_TYPE::m_straight));
        }
        else
        {
            _map_motion.insert(make_pair(_imu_crt, MOTION_TYPE::m_default));
        }
        if (meas_type != ZUPT_MEAS)
            _map_yaw.clear();
    }

    if (_time_valid(_gnss_crt, _imu_crt))
        _Meas_Type.insert(MEAS_TYPE::GNSS_MEAS);

    return _Meas_Type.size();
}

MEAS_TYPE gfgomsf::t_gfgo_gins::_getPOS(t_gposdata::data_pos &pos)
{
    //_mutex.lock();
    MEAS_TYPE res_type;
    double crt = _imu_crt.sow() + _imu_crt.dsec();
    t_gtime runEpoch = _gobs->load(_site, crt);
    int irc = ProcessOneEpoch(runEpoch);
    if (irc < 0)
    {
        return MEAS_TYPE::NO_MEAS;
    }
    _get_result(runEpoch, pos);

    MeasVel = pos.vn;
    MeasPos = pos.pos;
    tmeas = pos.t;
    _Cov_MeasVn = pos.Rvn;
    _Cov_MeasPos = pos.Rpos;

    if (!_aligned)
    {
        sins.pos = Cart2Geod(MeasPos, false);
        sins.vn = t_gbase::Cen(sins.pos).transpose() * MeasVel;
    }

    res_type = MEAS_TYPE::POS_MEAS;
    if (!double_eq(MeasVel.norm(), 0.0))
        res_type = POS_VEL_MEAS;

    if (pos.PDOP > _shm.max_pdop)
        res_type = NO_MEAS;
    if (pos.nSat < _shm.min_sat)
        res_type = NO_MEAS;

    return res_type;
}

MEAS_TYPE gfgomsf::t_gfgo_gins::_get_gnss_measurements(const t_gtime &runEpoch)
{

    _epoch = runEpoch;

    t_gtriple crdapr = _grec->crd_arp(_epoch);
    if (double_eq(crdapr[0], 0.0) && double_eq(crdapr[1], 0.0) && double_eq(crdapr[2], 0.0))
    {
        _valid_crd_xml = false;
    }
    else
    {
        _valid_crd_xml = true;
    }
    if (!_valid_crd_xml)
        _sig_init_crd = 100.0;

    Eigen::Vector3d XYZ_INS = sins.pos_ecef + sins.Ceb * lever;
    t_gtriple XYZ(XYZ_INS(0), XYZ_INS(1), XYZ_INS(2));
    _external_pos(XYZ, t_gtriple());

    if (_prepareData() < 0)
    {

        if (_spdlog)
            SPDLOG_LOGGER_ERROR(_spdlog, string("t_gfgo_gins "), ("_prepareData Failed!"));
        cur_sat_prn.clear();
        _initial_prior = true;
        return NO_MEAS;
    }
    if (_isBase)
        _set_rec_info(_gallobj->obj(_site_base)->crd_arp(_epoch), _vBanc(4), _vBanc_base(4));

    _get_initial_value(runEpoch); // for current epoch

    if (_data.size() < _minsat)
    {
        _rover_count--;

        if (_spdlog)
            SPDLOG_LOGGER_ERROR(_spdlog, string("t_gfgo_gins "), ("Not enough visible satellites!"));
        return NO_MEAS;
    }
    if (_combine_DD() < 0)
    {
        _rover_count--;

        if (_spdlog)
            SPDLOG_LOGGER_ERROR(_spdlog, string("t_gfgo_gins "), ("Combining the Double-Difference Pairs Failed!"));
        cur_sat_prn.clear();
        _initial_prior = true;
        return NO_MEAS;
    }

    return GNSS_MEAS;
}

bool gfgomsf::t_gfgo_gins::_get_imu_interval(double t0, double t1, vector<pair<double, Eigen::Vector3d>> &accVector,
                                             vector<pair<double, Eigen::Vector3d>> &gyrVector)
{

    if (_accBuf.empty())
    {
        printf("not receive imu\n");
        return false;
    }

    // if (fabs(t1 - _accBuf.back().first) <= _shm.delay)

    while (_accBuf.front().first <= t0)
    {
        _accBuf.pop();
        _gyrBuf.pop();
    }
    while (_accBuf.front().first < t1)
    {
        accVector.push_back(_accBuf.front());
        _accBuf.pop();
        gyrVector.push_back(_gyrBuf.front());
        _gyrBuf.pop();
        if (_accBuf.empty() || _gyrBuf.empty())
            break;
    }
    if (!_accBuf.empty() && !_gyrBuf.empty())
    {
        accVector.push_back(_accBuf.front());
        gyrVector.push_back(_gyrBuf.front());
    }

    assert(accVector.size() == gyrVector.size());
    return true;
}

bool gfgomsf::t_gfgo_gins::_imu_pre_integration(double t1, double t2)
{
    // NOTE: t2 > t1
    vector<pair<double, Eigen::Vector3d>> accVector, gyrVector;
    if (!_get_imu_interval(t1, t2, accVector, gyrVector))
        return false;

    for (size_t i = 0; i < accVector.size(); i++)
    {
        double dt;
        if (i == 0)
            dt = accVector[i].first - t1;
        else if (i == accVector.size() - 1)
            dt = t2 - accVector[i - 1].first;
        else
            dt = accVector[i].first - accVector[i - 1].first;

        double imu_ti = accVector[i].first;

        Eigen::Vector3d linear_acceleration = accVector[i].second;
        Eigen::Vector3d angular_velocity = gyrVector[i].second;

        if (!_first_imu)
        {
            _first_imu = true;
            _acc_0 = linear_acceleration;
            _gyr_0 = angular_velocity;
        }
        if (!_pre_integrations[_rover_count])
        {
            _pre_integrations[_rover_count] =
                new IntegrationBase{_acc_0, _gyr_0, _Bas[_rover_count], _Bgs[_rover_count]};
            _pre_integrations[_rover_count]->init_ins(_acc_n, _acc_w, _gyr_n, _gyr_w, _gravity);
        }
        if (_rover_count != 0)
        {
            _pre_integrations[_rover_count]->push_back(dt, linear_acceleration, angular_velocity);
            _tmp_pre_integration->push_back(dt, linear_acceleration, angular_velocity);
            _dt_buf[_rover_count].push_back(dt);
            _linear_acceleration_buf[_rover_count].push_back(linear_acceleration);
            _angular_velocity_buf[_rover_count].push_back(angular_velocity);
        }
        _acc_0 = linear_acceleration;
        _gyr_0 = angular_velocity;
    }

    return true;
}

int gfgomsf::t_gfgo_gins::_gnss_processing()
{
    std::cout << "gnss_now:" << setprecision(10) << (_gnss_crt.sow() + _gnss_crt.dsec()) << endl;

    t_gposdata::data_pos posdata;

    if (_msf_type == MSF_TYPE::GINS_TC_MODE)
    {
        _msf_flag = _get_gnss_measurements(_gnss_crt);
        if (_msf_flag != NO_MEAS)
            _c_gnss_factor = true;
        else
            _c_gnss_factor = false;
    }
    if (_msf_flag == NO_MEAS)
        return -1;
    return 1;
}

int gfgomsf::t_gfgo_gins::_gins_processing()
{
    _opt_valid = 0;
    if (_c_gnss_factor)
    {

        t_gnss_node gnssNode;
        gnssNode.pre_integration = _tmp_pre_integration;
        _all_gnss_node.insert(make_pair(_cur_node_time, gnssNode));
        _tmp_pre_integration = new IntegrationBase{_acc_0, _gyr_0, _Bas[_rover_count], _Bgs[_rover_count]};
        _tmp_pre_integration->init_ins(_acc_n, _acc_w, _gyr_n, _gyr_w, _gravity);
        if (1 && _solver_flag == INITIAL)
        {
            if (_rover_count == gins_window_size)
            {
                map<double, t_gnss_node>::iterator node_it;
                int i = 0;
                for (node_it = _all_gnss_node.begin(); node_it != _all_gnss_node.end(); ++node_it)
                {
                    node_it->second.R = _Rs[i];
                    node_it->second.T = _Ps[i];
                    i++;
                }

                for (int i = 0; i < _rover_count; i++)
                {
                    _pre_integrations[i]->repropagate(Eigen::Vector3d::Zero(), _Bgs[i]);
                }
                _gins_optimization();
                _gins_marginalization();
                _slide_gins();
                _solver_flag = NON_LINEAR;
            }
        }
        else
        {
            _gins_optimization();
            if (_c_gnss_factor)
                _gnss_amb_resolution();
            _gins_marginalization();
            _slide_gins();
        }
        _c_gnss_factor = false;
    }

    return 0;
}

void gfgomsf::t_gfgo_gins::_gins_feedback()
{
    if (_solver_flag == NON_LINEAR)
    {
        t_gposdata::data_pos pos;
        _get_result(_imu_crt, pos);

        cout << _imu_crt.str_hms() << _gnss_crt.str_hms() << _epoch.str_hms() << endl;
        cout << pos.pos.transpose() << endl;
        cout << _Ps[_rover_count].transpose() << endl;

        MEAS_TYPE meas_type = meas_state();

        Eigen::Vector3d mean_ba = _Bas[_rover_count];
        Eigen::Vector3d mean_bg = _Bgs[_rover_count];
        Eigen::Quaterniond e_q = Eigen::Quaterniond(_Rs[_rover_count]);
        e_q.normalized();
        sins.qeb = t_gquat(e_q.w(), e_q.x(), e_q.y(), e_q.z());
        sins.Ceb = t_gbase::q2mat(sins.qeb);
        sins.ve = _Vs[_rover_count];
        sins.pos_ecef = _Ps[_rover_count];
        Eigen::Vector3d robustpos;
        if (_getRobustFixedPosition(robustpos))
            sins.pos_ecef = robustpos;
        // sins.pos_ecef = pos.pos;
        sins.eb = mean_bg;
        sins.db = mean_ba;
        sins.qnb = t_gbase::m2qua(sins.eth.Cne) * sins.qeb;
        sins.Cnb = t_gbase::q2mat(sins.qnb);
        sins.vn = sins.eth.Cne * sins.ve;
        sins.pos = Cart2Geod(sins.pos_ecef, false);
        sins.att = t_gbase::q2att(sins.qnb);
        sins.xyz_out = pos.pos;
    }
}

void gfgomsf::t_gfgo_gins::_gravity_update()
{
    Eigen::Vector3d g_n = -sins.eth.gcc;
    Eigen::Vector3d g_e = sins.eth.Cen * g_n;
    // cout << " gravity_n: " << g_n << endl;
    // cout << " gravity_e: " << g_e << endl;
    if (_nav_frame == NAV_REFERENCE_FRAME::E_F)
        _gravity = g_e;
    else
        _gravity = g_n;
}
void gfgomsf::t_gfgo_gins::get_imu(great::t_gimudata *_imu)
{
    _imudata = _imu;
}

void gfgomsf::t_gfgo_gins::_gins_double_to_vector()
{
    for (int i = 0; i <= _rover_count; i++)
    {
        _Rs[i] = Eigen::Quaterniond(_para_pose[i][6], _para_pose[i][3], _para_pose[i][4], _para_pose[i][5])
                     .normalized()
                     .toRotationMatrix();
        _Ps[i] = Eigen::Vector3d(_para_pose[i][0], _para_pose[i][1], _para_pose[i][2]);
        if (_imu_enable)
        {
            _Vs[i] = Eigen::Vector3d(_para_speed_bias[i][0], _para_speed_bias[i][1], _para_speed_bias[i][2]);

            _Bas[i] = Eigen::Vector3d(_para_speed_bias[i][3], _para_speed_bias[i][4], _para_speed_bias[i][5]);

            _Bgs[i] = Eigen::Vector3d(_para_speed_bias[i][6], _para_speed_bias[i][7], _para_speed_bias[i][8]);
        }
    }
    //}
    for (int i = 0; i < _amb_manager->ambiguity_ids.size(); i++)
    {
        int amb_id = _amb_manager->ambiguity_ids[i];
        _amb_manager->updateAmb(amb_id, _para_amb[amb_id][0]);
    }
}

void gfgomsf::t_gfgo_gins::_gins_vector_to_double()
{
    for (int i = 0; i <= _rover_count; i++)
    {
        _para_pose[i][0] = _Ps[i].x();
        _para_pose[i][1] = _Ps[i].y();
        _para_pose[i][2] = _Ps[i].z();
        Eigen::Quaterniond q{_Rs[i]};
        _para_pose[i][3] = q.x();
        _para_pose[i][4] = q.y();
        _para_pose[i][5] = q.z();
        _para_pose[i][6] = q.w();

        if (_imu_enable)
        {
            _para_speed_bias[i][0] = _Vs[i].x();
            _para_speed_bias[i][1] = _Vs[i].y();
            _para_speed_bias[i][2] = _Vs[i].z();

            _para_speed_bias[i][3] = _Bas[i].x();
            _para_speed_bias[i][4] = _Bas[i].y();
            _para_speed_bias[i][5] = _Bas[i].z();

            _para_speed_bias[i][6] = _Bgs[i].x();
            _para_speed_bias[i][7] = _Bgs[i].y();
            _para_speed_bias[i][8] = _Bgs[i].z();
        }
    }

    for (int i = 0; i < _amb_manager->ambiguity_ids.size(); i++)
    {
        int amb_id = _amb_manager->ambiguity_ids[i];
        _para_amb[amb_id][0] = _amb_manager->getAmb(amb_id);
    }
    _amb_manager->generateAmbSearchIndex();
}

void gfgomsf::t_gfgo_gins::_gins_optimization()
{
    /*if (_cur_node_time == 200340.0)
        cerr << endl;*/

    if (_rover_count < 1)
    {
        _c_gnss_factor = false;
        return;
    }
    t_tictoc fgo_gins;
    _removed_sats.clear();
    int count = -1;
    bool iter_flag = false;
    pair<string, int> outlier = make_pair(" ", -1);
    do
    {
        count++;
        _remove_outlier_sat(outlier);
        _gins_vector_to_double();
        // if (_vDD_msg[_rover_count].size() <= 1)
        if (_vDD_msg[_rover_count].size() < _minsat * 4)
        {
            _last_gnss_info->valid = false;
            _opt_valid = 0;
            break;
        }
        ceres::Problem problem;
        ceres::LossFunction *loss_function;
        loss_function = new ceres::HuberLoss(_loss_func_value);

        for (int i = 0; i <= _rover_count; i++)
        {
            ceres::LocalParameterization *local_parameterization = new PoseLocalParameterization();
            problem.AddParameterBlock(_para_pose[i], SIZE_POSE, local_parameterization);
            if (_imu_enable)
                problem.AddParameterBlock(_para_speed_bias[i], SIZE_SPEEDBIAS);

            Eigen::Vector3d pos(_para_pose[i][0], _para_pose[i][1], _para_pose[i][2]);
            Eigen::Quaterniond quat(_para_pose[i][6], _para_pose[i][3], _para_pose[i][4], _para_pose[i][5]);
            InitialPoseFactor *initial_pose = new InitialPoseFactor(pos, quat);
            initial_pose->sqrt_info = 1e-7 * Eigen::Matrix<double, 6, 6>::Identity();
            problem.AddResidualBlock(initial_pose, NULL, _para_pose[i]);

            InitialVelBiasFactor *initial_bias = new InitialVelBiasFactor(
                Eigen::Vector3d(_para_speed_bias[i][0], _para_speed_bias[i][1], _para_speed_bias[i][2]),
                Eigen::Vector3d(_para_speed_bias[i][3], _para_speed_bias[i][4], _para_speed_bias[i][5]),
                Eigen::Vector3d(_para_speed_bias[i][6], _para_speed_bias[i][7], _para_speed_bias[i][8]));
            initial_bias->sqrt_info = 1e-7 * Eigen::Matrix<double, 9, 9>::Identity();
            problem.AddResidualBlock(initial_bias, NULL, _para_speed_bias[i]);
        }

        /*for (int i = 0; i < _amb_manager->ambiguity_ids.size(); i++)
        {
            int amb_id = _amb_manager->ambiguity_ids[i];
            problem.AddParameterBlock(_para_amb[amb_id], SIZE_AMB);
        }*/

        // prior

        if (_last_marginalization_info && _last_marginalization_info->valid)
        {
            // Constructing Marginalization Factors
            MarginalizationFactor *marginalization_factor = new MarginalizationFactor(_last_marginalization_info);
            problem.AddResidualBlock(marginalization_factor, NULL, _last_marginalization_parameter_blocks);
        }

        // IMU

        if (_rover_count > 0)
        {
            for (int i = 0; i < _rover_count; i++)
            {
                int j = i + 1;
                if (_pre_integrations[j]->sum_dt > _gtime_interval)
                    continue;
                auto preint = _pre_integrations[j];
                IMUFactor *imu_factor = new IMUFactor(_pre_integrations[j]);
                problem.AddResidualBlock(imu_factor, NULL, _para_pose[i], _para_speed_bias[i], _para_pose[j],
                                         _para_speed_bias[j]);
            }
        }

        // gnss
        _obs_index.clear();
        assert(_vDD_msg.size() == _rover_count + 1);
        for (int i = 0; i <= _rover_count; i++)
        {
            vector<DDEquMsg> DD_tmp = _vDD_msg[i];
            t_gtime crt = DD_tmp.begin()->time;
            map<t_gtime, vector<t_gsatdata>>::const_iterator base_iter = _map_basedata.find(crt);
            map<t_gtime, t_gallpar>::const_iterator param_iter = _map_param.find(crt);
            if (base_iter == _map_basedata.end() || param_iter == _map_param.end())
                continue;
            int dd_equ_count = 0;
            for (auto &dd_iter : DD_tmp)
            {
                /*if (_cur_node_time == 200340.0 && i == 2 && dd_iter.ref_sat == "G10" && dd_iter.nonref_sat == "G25")
                    cerr << endl;*/

                if (!_get_DD_data(dd_iter, base_iter->second))
                    continue;
                pair<string, string> base_rover_site = make_pair(dd_iter.base_site, dd_iter.rover_site);
                pair<FREQ_SEQ, GOBSBAND> freq_band = make_pair(dd_iter.freq, dd_iter.band);
                vector<pair<t_gsatdata, t_gsatdata>> DD_sat_data;
                DD_sat_data.push_back(make_pair(dd_iter.base_ref_sat, dd_iter.rover_ref_sat));
                DD_sat_data.push_back(make_pair(dd_iter.base_nonref_sat, dd_iter.rover_nonref_sat));

                // cout << "base_ref_sat:" << dd_iter.base_ref_sat.epoch().sow() <<
                //	"  rover_ref_sat:" << dd_iter.rover_ref_sat.epoch().sow() << endl;
                // cout << "base_nonref_sat:" << dd_iter.base_nonref_sat.epoch().sow() <<
                //	"  rover_nonref_sat:" << dd_iter.rover_nonref_sat.epoch().sow() << endl;

                GOBSTYPE obstype = dd_iter.obs_type;
                _obs_index.push_back(make_pair(dd_iter.rover_nonref_sat.sat(), make_pair(dd_iter.freq, obstype)));
                if (obstype == GOBSTYPE::TYPE_C)
                {
                    PseudorangeDDINGFactor *pinsf = new PseudorangeDDINGFactor(
                        dd_iter.time, base_rover_site, param_iter->second, DD_sat_data, _gbias_model, freq_band, lever);
                    problem.AddResidualBlock(pinsf, NULL, _para_pose[i]);

                    /*double **para = new double *[3];
                    para[0] = _para_pose[i];
                    pinsf->check(para);*/
                }
                if (obstype == GOBSTYPE::TYPE_L)
                {
                    int id1, id2;
                    id1 = _amb_manager->getAmbSearchIndex(make_pair(dd_iter.ref_sat_global_id, dd_iter.freq));
                    id2 = _amb_manager->getAmbSearchIndex(make_pair(dd_iter.nonref_sat_global_id, dd_iter.freq));
                    if (id1 == -1 || id2 == -1)
                        continue;
                    CarrierphaseDDINGFactor *linsf = new CarrierphaseDDINGFactor(
                        dd_iter.time, base_rover_site, param_iter->second, DD_sat_data, _gbias_model, freq_band, lever);
                    problem.AddResidualBlock(linsf, NULL, _para_pose[i], _para_amb[id1], _para_amb[id2]);

                    /*double **para = new double *[3];
                    para[0] = _para_pose[i];
                    para[1] = _para_amb[id1];
                    para[2] = _para_amb[id2];
                    linsf->check(para);*/
                }
                dd_equ_count++;
            }
            assert(dd_equ_count == DD_tmp.size());
        }

        // static constraints
        // for (int i = 0; i <= _rover_count; i++)
        //{
        //	vector<DDEquMsg> DD_tmp = _vDD_msg[i];
        //	t_gtime crt = DD_tmp.begin()->time;
        //	//if (crt.sow() + crt.dsec() > 449600 && crt.sow() + crt.dsec() < 449640)
        //	//	continue;
        //	map<t_gtime, MOTION_TYPE>::const_iterator it = _map_motion.find(crt);
        //	if (it->second == MOTION_TYPE::m_static)
        //	{
        //		for (it; it != _map_motion.begin(); --it)
        //		{
        //			if (it->second != MOTION_TYPE::m_static)
        //			{
        //				++it; break;
        //			}
        //		}
        //		map<t_gtime, pair<Eigen::Vector3d, Eigen::Quaterniond>>::const_iterator it_pose =
        //_map_pose.find(it->first); 		Eigen::Vector3d pos = it_pose->second.first; 		Eigen::Quaterniond quat =
        //it_pose->second.second; 		InitialPoseFactor* pose_factor = new InitialPoseFactor(pos, quat);
        //		pose_factor->sqrt_info = 1e4 * Eigen::Matrix<double, 6, 6>::Identity();
        //		//ceres::CostFunction* pose_factor = PoseError::Create(pos.x(), pos.y(), pos.z(),
        //		//	quat.w(), quat.x(), quat.y(), quat.z(), 1e-6, 1e-9);
        //		problem.AddResidualBlock(pose_factor, NULL, _para_pose[i]);
        //	}
        //}

        // if (_initial_prior)

        if (1)
        {
            for (int i = 0; i < _amb_manager->ambiguity_ids.size(); i++)
            {
                int amb_id = _amb_manager->ambiguity_ids[i];
                double initial_amb = _para_amb[amb_id][0];
                InitialGnssAMB *amb_prior = new InitialGnssAMB(initial_amb);
                problem.AddResidualBlock(amb_prior, NULL, _para_amb[amb_id]);
            }
        }

        ///*if (_cur_node_time == 200340.0)
        //{*/
        //	Eigen::MatrixXd jacobian = GetFullJacobian(problem);
        //	std::cout << jacobian << endl;
        ////}
        // AnalyzeFullJacobian(problem);

        // optimization
        // ceres::Solver::Options options;
        // options.minimizer_type = ceres::LINE_SEARCH;
        ////options.trust_region_strategy_type = ceres::DOGLEG;
        // options.max_num_iterations = 1;

        ceres::Solver::Options options;
        options.linear_solver_type = ceres::DENSE_SCHUR;
        options.trust_region_strategy_type = ceres::DOGLEG;
        // options.trust_region_strategy_type = ceres::LEVENBERG_MARQUARDT;
        // options.use_nonmonotonic_steps = false;
        // options.min_trust_region_radius = options.max_trust_region_radius = 1e6;
        options.max_num_iterations = 5;
        ceres::Solver::Summary summary;
        ceres::Solve(options, &problem, &summary);
        std::cout << summary.BriefReport() << endl;
        _opt_valid = 1;
        //_gins_posteriori_test();
        _gins_posteriori_test(problem);
        if (_gobs_outlier_detection(outlier) >= 0)
            iter_flag = true;
        else
            iter_flag = false;
        //_detect_outlier()
    } while (iter_flag);

    std::cout << "Epoch: " << _cur_node_time << " time cost (ms): " << fgo_gins.toc() << " iter times: " << count
              << endl;
    if (_last_gnss_info->valid)
    {
        std::cout << "sigma: " << _last_gnss_info->sig_unit << "  vtpv: " << _last_gnss_info->vtpv << std::endl;
        _gins_double_to_vector();
        _opt_flag = true;

        if (_spdlog)
            _spdlog->info("Epoch: {} GNSS/SINS Integrated Navigation Processing Succeed at Epoch: {}", _cur_node_time,
                          _imu_crt.str_ymdhms(""));
    }
    else
    {

        if (_spdlog)
            _spdlog->warn("Epoch: {} GNSS/SINS Integrated Navigation Processing Failed at Epoch: {}", _cur_node_time,
                          _imu_crt.str_ymdhms(""));
    }
}

void gfgomsf::t_gfgo_gins::_gins_marginalization()
{
    // marginalization
    if (_rover_count == gins_window_size)
    {
        ceres::LossFunction *loss_function;
        loss_function = new ceres::HuberLoss(_loss_func_value);
        MarginalizationInfo *marginalization_info = new MarginalizationInfo();
        if (_vDD_msg[_rover_count].size() > 1)
        {
            _gins_vector_to_double();
            // marginalization prior
            if (_last_marginalization_info && _last_marginalization_info->valid)
            {
                vector<int> drop_set;
                vector<int> amb_margin = _amb_manager->getMarginAmb();
                for (int i = 0; i < static_cast<int>(_last_marginalization_parameter_blocks.size()); i++)
                {
                    if (_last_marginalization_parameter_blocks[i] == _para_pose[0] ||
                        _last_marginalization_parameter_blocks[i] == _para_speed_bias[0])
                        drop_set.push_back(i);

                    for (int j = 0; j < amb_margin.size(); j++)
                    {
                        if (_last_marginalization_parameter_blocks[i] == _para_amb[amb_margin[j]])
                            drop_set.push_back(i);
                    }
                }
                // construct new marginlization_factor
                MarginalizationFactor *marginalization_factor = new MarginalizationFactor(_last_marginalization_info);
                ResidualBlockInfo *residual_block_info = new ResidualBlockInfo(
                    marginalization_factor, NULL, _last_marginalization_parameter_blocks, drop_set);
                marginalization_info->addResidualBlockInfo(residual_block_info);
            }
            // marginalization INS
            if (_imu_enable)
            {
                if (_pre_integrations[1]->sum_dt <= _gtime_interval)
                {
                    IMUFactor *imu_factor = new IMUFactor(_pre_integrations[1]);
                    ResidualBlockInfo *residual_block_info = new ResidualBlockInfo(
                        imu_factor, NULL,
                        vector<double *>{_para_pose[0], _para_speed_bias[0], _para_pose[1], _para_speed_bias[1]},
                        vector<int>{0, 1});
                    marginalization_info->addResidualBlockInfo(residual_block_info);
                }
            }

            ////marginalization gnss
            for (int i = 0; i <= _rover_count; i++)
            {
                if (i == 0)
                {
                    vector<DDEquMsg> dd_msg = _vDD_msg[i];
                    t_gtime crt = dd_msg.begin()->time;
                    map<t_gtime, vector<t_gsatdata>>::const_iterator base_iter = _map_basedata.find(crt);
                    map<t_gtime, t_gallpar>::const_iterator param_iter = _map_param.find(crt);
                    if (base_iter == _map_basedata.end() || param_iter == _map_param.end())
                        continue;
                    for (auto &dd_iter : dd_msg)
                    {
                        if (!_get_DD_data(dd_iter, base_iter->second))
                            continue;
                        pair<string, string> base_rover_site = make_pair(dd_iter.base_site, dd_iter.rover_site);
                        pair<FREQ_SEQ, GOBSBAND> freq_band = make_pair(dd_iter.freq, dd_iter.band);
                        vector<pair<t_gsatdata, t_gsatdata>> DD_sat_data;
                        DD_sat_data.push_back(make_pair(dd_iter.base_ref_sat, dd_iter.rover_ref_sat));
                        DD_sat_data.push_back(make_pair(dd_iter.base_nonref_sat, dd_iter.rover_nonref_sat));

                        GOBSTYPE obstype = dd_iter.obs_type;
                        if (obstype == GOBSTYPE::TYPE_C)
                        {
                            PseudorangeDDINGFactor *pinsf =
                                new PseudorangeDDINGFactor(dd_iter.time, base_rover_site, param_iter->second,
                                                           DD_sat_data, _gbias_model, freq_band, lever);
                            ResidualBlockInfo *residual_block_info =
                                new ResidualBlockInfo(pinsf, NULL, vector<double *>{_para_pose[0]}, vector<int>{0});
                            marginalization_info->addResidualBlockInfo(residual_block_info);
                        }
                        if (obstype == GOBSTYPE::TYPE_L)
                        {
                            int id1, id2;
                            id1 = _amb_manager->getAmbSearchIndex(make_pair(dd_iter.ref_sat_global_id, dd_iter.freq));
                            id2 =
                                _amb_manager->getAmbSearchIndex(make_pair(dd_iter.nonref_sat_global_id, dd_iter.freq));
                            if (id1 == -1 || id2 == -1)
                                continue;
                            vector<int> drop_set{0};
                            if (_amb_manager->getAmbStartRoverID(id1) == 0 &&
                                _amb_manager->getAmbStartRoverID(id2) == 0)
                            {
                                if (_amb_manager->getAmbEndRoverID(id1) == 0)
                                    drop_set.push_back(1);
                                if (_amb_manager->getAmbEndRoverID(id2) == 0)
                                    drop_set.push_back(2);

                                CarrierphaseDDINGFactor *linsf =
                                    new CarrierphaseDDINGFactor(dd_iter.time, base_rover_site, param_iter->second,
                                                                DD_sat_data, _gbias_model, freq_band, lever);
                                ResidualBlockInfo *residual_block_info = new ResidualBlockInfo(
                                    linsf, NULL, vector<double *>{_para_pose[0], _para_amb[id1], _para_amb[id2]},
                                    drop_set);
                                marginalization_info->addResidualBlockInfo(residual_block_info);
                            }
                        }
                    }
                }
            }

            // static constraints
            for (int i = 0; i <= _rover_count; i++)
            {
                if (i == 0)
                {
                    vector<DDEquMsg> DD_tmp = _vDD_msg[i];
                    t_gtime crt = DD_tmp.begin()->time;
                    map<t_gtime, MOTION_TYPE>::const_iterator it = _map_motion.find(crt);
                    if (it->second == MOTION_TYPE::m_static)
                    {
                        for (it; it != _map_motion.begin(); --it)
                        {
                            if (it->second != MOTION_TYPE::m_static)
                            {
                                ++it;
                                break;
                            }
                        }
                        map<t_gtime, pair<Eigen::Vector3d, Eigen::Quaterniond>>::const_iterator it_pose =
                            _map_pose.find(it->first);
                        Eigen::Vector3d pos = it_pose->second.first;
                        Eigen::Quaterniond quat = it_pose->second.second;
                        InitialPoseFactor *pose_factor = new InitialPoseFactor(pos, quat);
                        pose_factor->sqrt_info = 1e4 * Eigen::Matrix<double, 6, 6>::Identity();
                        // ceres::CostFunction* pose_factor = PoseError::Create(pos.x(), pos.y(), pos.z(),
                        //	quat.w(), quat.x(), quat.y(), quat.z(), 1e-6, 1e-9);
                        ResidualBlockInfo *residual_block_info =
                            new ResidualBlockInfo(pose_factor, NULL, vector<double *>{_para_pose[0]}, vector<int>{0});
                        marginalization_info->addResidualBlockInfo(residual_block_info);
                    }
                }
            }

            marginalization_info->preMarginalize();
            // printf("pre marginalization %f ms\n", t_pre_margin.toc());
            // TicToc t_margin;
            marginalization_info->marginalize();
            // printf("marginalization %f ms\n", t_margin.toc());
            std::map<long, double *> addr_shift;
            for (int i = 1; i <= _rover_count; i++)
            {
                addr_shift[reinterpret_cast<long>(_para_pose[i])] = _para_pose[i - 1];
                if (_imu_enable)
                    addr_shift[reinterpret_cast<long>(_para_speed_bias[i])] = _para_speed_bias[i - 1];
            }
            // gnss
            vector<int> cur_amb = _amb_manager->getCurWinAmb();
            for (int i = 0; i < cur_amb.size(); i++)
            {
                addr_shift[reinterpret_cast<long>(_para_amb[cur_amb[i]])] = _para_amb[cur_amb[i]];
            }

            vector<double *> parameter_blocks =
                marginalization_info->getParameterBlocks(addr_shift); /// reserved parameters
            _last_marginalization_parameter_blocks = parameter_blocks;
        }
        if (_last_marginalization_info)
            delete _last_marginalization_info;
        _last_marginalization_info = marginalization_info;

        if (_last_marginalization_info->factors.size() == 0)
            _last_marginalization_info->valid = false;
        _initial_prior = false;
    }
}

void gfgomsf::t_gfgo_gins::_slide_gins()
{
    double t_0 = t_gfgo_gins::_headers[0];
    if (_rover_count == gins_window_size)
    {
        for (int i = 0; i < _rover_count; i++)
        {
            t_gfgo_gins::_headers[i] = t_gfgo_gins::_headers[i + 1];
            _Rs[i].swap(_Rs[i + 1]);
            _Ps[i].swap(_Ps[i + 1]);
            if (_imu_enable)
            {
                std::swap(_pre_integrations[i], _pre_integrations[i + 1]);

                _dt_buf[i].swap(_dt_buf[i + 1]);
                _linear_acceleration_buf[i].swap(_linear_acceleration_buf[i + 1]);
                _angular_velocity_buf[i].swap(_angular_velocity_buf[i + 1]);

                _Vs[i].swap(_Vs[i + 1]);
                _Bas[i].swap(_Bas[i + 1]);
                _Bgs[i].swap(_Bgs[i + 1]);
            }
        }

        if (_imu_enable)
        {
            delete _pre_integrations[_rover_count];
            _pre_integrations[_rover_count] =
                new IntegrationBase{_acc_0, _gyr_0, _Bas[_rover_count], _Bgs[_rover_count]};
            _pre_integrations[_rover_count]->init_ins(_acc_n, _acc_w, _gyr_n, _gyr_w, _gravity);
            /*delete _pre_integrations[WINDOW_SIZE];
            _pre_integrations[WINDOW_SIZE] = new IntegrationBase{ _acc_0, _gyr_0, _Bas[WINDOW_SIZE], _Bgs[WINDOW_SIZE]
            };*/
            _dt_buf[_rover_count].clear();
            _linear_acceleration_buf[_rover_count].clear();
            _angular_velocity_buf[_rover_count].clear();
        }

        if (true || _solver_flag == INITIAL)
        {
            // pre_integration = nullptr;
            map<double, t_gnss_node>::iterator it_0 = _all_gnss_node.find(t_0);

            delete it_0->second.pre_integration;
            _all_gnss_node.erase(_all_gnss_node.begin(), it_0);
        }

        while (!_map_motion.empty() && _map_motion.rbegin()->first - _map_motion.begin()->first > 60 * 10)
        {
            _map_motion.erase(std::begin(_map_motion));
        }
        // slide gnss
        _vDD_msg.erase(_vDD_msg.begin());
        _amb_manager->slidingWindow();
        _rover_count--;
    }
}

void gfgomsf::t_gfgo_gins::_gins_posteriori_test(ceres::Problem &problem)
{
    _all_para_win.delAllParam();
    ceres::LossFunction *loss_function;
    loss_function = new ceres::HuberLoss(_loss_func_value);
    // loss_function = new ceres::CauchyLoss(_loss_func_value);
    GNSSInfo *gnss_info = new GNSSInfo();
    if (_vDD_msg[_rover_count].size() >= 1)
    {
        // constrcut window para_index
        vector<vector<int>> pose_para_col_index;
        vector<vector<int>> amb_para_col_index;
        int total_para_size = 0;
        vector<double *> _parameter_blocks;
        vector<par_type> crd_partype{par_type::CRD_X, par_type::CRD_Y, par_type::CRD_Z};
        vector<par_type> att_partype{par_type::ATT_X, par_type::ATT_Y, par_type::ATT_Z};
        vector<int> posei;
        t_gquat qi = t_gquat(_para_pose[_rover_count][6], _para_pose[_rover_count][3], _para_pose[_rover_count][4],
                             _para_pose[_rover_count][5]);
        // qi.normlize(qi);
        Eigen::Vector3d rv = t_gbase::q2rv(qi);
        _parameter_blocks.push_back(_para_pose[_rover_count]);
        for (int j = 0; j < 6; j++)
        {
            if (j < 3)
            {
                t_gpar par_crd;
                par_crd.site = _site;
                par_crd.parType = crd_partype[j];
                par_crd.value(_para_pose[_rover_count][j]);
                par_crd.beg = _gnss_crt;
                par_crd.end = _gnss_crt;
                // par_crd.index = i * 6 + j + 1;
                par_crd.index = j + 1;
                _all_para_win.addParam(par_crd);
            }
            if (j >= 3)
            {
                t_gpar par_att;
                par_att.site = _site;
                par_att.parType = att_partype[j - 3];
                par_att.value(rv(j - 3));
                par_att.beg = _gnss_crt;
                par_att.end = _gnss_crt;
                // par_att.index = i * 6 + j + 1;
                par_att.index = j + 1;
                _all_para_win.addParam(par_att);
            }
            int id = j;
            posei.push_back(id);
            total_para_size = total_para_size + 1;
        }
        pose_para_col_index.push_back(posei);

        int amb_index_start = 6;

        map<int, int> amb_col_id;
        int amb_size = -1;
        _gnss_obs_index.clear();
        int unused_DD_nums = 0;
        for (int i = _rover_count; i <= _rover_count; i++)
        {
            vector<DDEquMsg> dd_msg = _vDD_msg[i];
            t_gtime crt = dd_msg.begin()->time;
            map<t_gtime, vector<t_gsatdata>>::const_iterator base_iter = _map_basedata.find(crt);
            map<t_gtime, t_gallpar>::const_iterator param_iter = _map_param.find(crt);
            if (base_iter == _map_basedata.end() || param_iter == _map_param.end())
                continue;
            for (auto &dd_iter : dd_msg)
            {
                if (!_get_DD_data(dd_iter, base_iter->second))
                {
                    unused_DD_nums++;
                    continue;
                }
                pair<string, string> base_rover_site = make_pair(dd_iter.base_site, dd_iter.rover_site);
                pair<FREQ_SEQ, GOBSBAND> freq_band = make_pair(dd_iter.freq, dd_iter.band);
                vector<pair<t_gsatdata, t_gsatdata>> DD_sat_data;
                DD_sat_data.push_back(make_pair(dd_iter.base_ref_sat, dd_iter.rover_ref_sat));
                DD_sat_data.push_back(make_pair(dd_iter.base_nonref_sat, dd_iter.rover_nonref_sat));

                GOBSTYPE obstype = dd_iter.obs_type;
                _gnss_obs_index.push_back(
                    make_pair(make_pair(dd_iter.rover_nonref_sat.sat(), dd_iter.nonref_sat_global_id),
                              make_pair(dd_iter.freq, obstype)));
                if (obstype == GOBSTYPE::TYPE_C)
                {
                    PseudorangeDDINGFactor *pinsf = new PseudorangeDDINGFactor(
                        dd_iter.time, base_rover_site, param_iter->second, DD_sat_data, _gbias_model, freq_band, lever);
                    GNSSResidualBlockInfo *residual_block_info =
                        new GNSSResidualBlockInfo(pinsf, NULL, vector<double *>{_para_pose[i]});
                    map<long, vector<int>> para_index;
                    para_index[reinterpret_cast<long>(_para_pose[i])] = pose_para_col_index[0];
                    gnss_info->addResidualBlockInfo(residual_block_info, para_index);
                }
                if (obstype == GOBSTYPE::TYPE_L)
                {
                    int index = dd_iter.ref_sat_global_id;
                    int nonindex = dd_iter.nonref_sat_global_id;
                    vector<int> amb_id12(2);
                    amb_id12[0] = _amb_manager->getAmbSearchIndex(make_pair(dd_iter.ref_sat_global_id, dd_iter.freq));
                    amb_id12[1] =
                        _amb_manager->getAmbSearchIndex(make_pair(dd_iter.nonref_sat_global_id, dd_iter.freq));

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
                    // cout << "para_col_index: " << para_col_index.size() << " id1: " << _rover_count + 1+ id1 << "
                    // id2: " << _rover_count + 1+id2 << endl;
                    double *addr1 = _para_amb[amb_id12[0]];
                    double *addr2 = _para_amb[amb_id12[1]];
                    para_index[reinterpret_cast<long>(_para_pose[i])] = pose_para_col_index[0];
                    para_index[reinterpret_cast<long>(addr1)] = vector<int>{amb_col_id[amb_id12[0]]};
                    para_index[reinterpret_cast<long>(addr2)] = vector<int>{amb_col_id[amb_id12[1]]};
                    CarrierphaseDDINGFactor *linsf = new CarrierphaseDDINGFactor(
                        dd_iter.time, base_rover_site, param_iter->second, DD_sat_data, _gbias_model, freq_band, lever);
                    GNSSResidualBlockInfo *residual_block_info = new GNSSResidualBlockInfo(
                        linsf, NULL, vector<double *>{_para_pose[i], _para_amb[amb_id12[0]], _para_amb[amb_id12[1]]});
                    gnss_info->addResidualBlockInfo(residual_block_info, para_index);
                }
            }
        }

        total_para_size = total_para_size + amb_size + 1;
        // construct variances by ceres solver
        ceres::Covariance::Options options_co;
        options_co.algorithm_type = ceres::DENSE_SVD;
        options_co.min_reciprocal_condition_number = 1e-50;
        // options_co.sparse_linear_algebra_library_type = ceres::SparseLinearAlgebraLibraryType::SUITE_SPARSE;
        // options_co.apply_loss_function = false; //optional, true or false is depended on the reliability of
        // covariance
        ceres::Covariance covariance(options_co);
        std::vector<const double *> covariance_blocks; // all parameter_blocks
        for (int i = 0; i < _parameter_blocks.size(); i++)
        {
            covariance_blocks.push_back(_parameter_blocks[i]);
        }
        Eigen::MatrixXd Qx = Eigen::MatrixXd::Zero(total_para_size + 1, total_para_size + 1);
        try
        {
            covariance.Compute(covariance_blocks, &problem);
            covariance.GetCovarianceMatrix(covariance_blocks, Qx.data());
            // std::cout.precision(15);
            t_gbase::delrowcol(Qx, 6);
            // cout << "Qx: " << endl;
            // std::cout << Qx << std::endl;

            // if (_vDD_msg.begin()->begin()->time.sow() > 447786)
            //{
            //	covariance.Compute(covariance_blocks, &problem);
            //	covariance.GetCovarianceMatrix(covariance_blocks, Qx.data());
            //	//std::cout.precision(15);
            //	t_gbase::delrowcol(Qx, 6);
            //	//cout << "Qx: " << endl;
            //	//std::cout << Qx << std::endl;
            // }
        }
        catch (...)
        {
            cout << "Covariance Compute Failed" << endl;
        }

        if (_vDD_msg[_rover_count].size() - unused_DD_nums >= 1)
        {
            gnss_info->constructEqu_fromCeres(Qx);
        }
        else
        {
            cout << "DD_msgs: " << _vDD_msg[_rover_count].size() << "   unused_DD_nums: " << unused_DD_nums << endl;
        }
    }
    if (_last_gnss_info)
        delete _last_gnss_info;
    _last_gnss_info = gnss_info;
}

bool gfgomsf::t_gfgo_gins::_getRobustFixedPosition(Eigen::Vector3d &pos)
{
    set<string> ambs = _param.amb_prns();
    int nsat = ambs.size();
    if (_amb_state && nsat > 10 && _ambfix->get_ratio() > 5)
    {
        t_gtriple xyz;
        _param_fixed.getCrdParam(_site, xyz);
        pos = Eigen::Vector3d(xyz.crd_array());
        return true;
    }
    else
    {
        pos = Eigen::Vector3d();
        return false;
    }
}

void gfgomsf::t_gfgo_gins::_gnss_amb_resolution()
{
    if (_last_gnss_info->valid)
    {
        _pre_amb_resolution();
        _amb_resolution();
    }
}
void gfgomsf::t_gfgo_gins::_get_initial_value(const t_gtime &runEpoch)
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
void gfgomsf::t_gfgo_gins::_set_initial_value(const t_gtime &runEpoch)
{
    _rover_count++;
    if (_rover_count >= 1)
    {
        if (_msf_type == MSF_TYPE::GINS_TC_MODE)
        {
            double time_gap = runEpoch.sow() - _last_pre_integration_time;
            // detect time gap
            if (time_gap > _gtime_interval)
            {
                cout << "[" << _last_pre_integration_time << "]--" << "[" << runEpoch.sow() << "]"
                     << ":gnss time gap !!!]" << endl;
                clearGNSSmsg();
                _rover_count++;
                _last_pre_integration_time = -1;
                while (!_accBuf.empty())
                {
                    _accBuf.pop();
                    _gyrBuf.pop();
                }
                if (_last_gnss_marginalization_info)
                    _last_marginalization_info->valid = false;
            }
        }
    }
    t_gallpar params_add;
    _gtemp_params(_param, params_add);
    _map_basedata.insert(make_pair(runEpoch, _data_base));
    _map_param.insert(make_pair(runEpoch, params_add));
    //_rover_window[_rover_count]->setSatData(_data);
    vector<string> cur_sat_list;
    for (auto it : _data)
        cur_sat_list.push_back(it.sat());

    if (_rover_count >= 1 && cur_sat_prn.empty())
        cout << "Refill current satellite list !" << endl;

    map<string, int>::iterator iter_cur_prn = cur_sat_prn.begin();
    for (; iter_cur_prn != cur_sat_prn.end();)
    {
        string sat = iter_cur_prn->first;
        auto it_find = find_if(cur_sat_list.begin(), cur_sat_list.end(),
                               [sat](const string &sat_name) { return sat_name == sat; });
        if (it_find == cur_sat_list.end())
        {
            iter_cur_prn = cur_sat_prn.erase(iter_cur_prn);
            cout << "the sat : " << sat << " lost tracking !" << endl;
        }
        else
            iter_cur_prn++;
    }

    for (auto it : _data)
    {
        string sat_name = it.sat();
        auto it_find = cur_sat_prn.find(sat_name);
        if (it_find == cur_sat_prn.end())
        {
            _global_sat_id++;
            cur_sat_prn[sat_name] = _global_sat_id;
            _amb_manager->addNewSat(runEpoch, _rover_count, _global_sat_id, _global_amb_id, it, _param);
        }
        else
        {
            if (it.islip())
            {
                cur_sat_prn.erase(sat_name);
                _global_sat_id++;
                cur_sat_prn[sat_name] = _global_sat_id;
                _amb_manager->addNewSat(runEpoch, _rover_count, _global_sat_id, _global_amb_id, it, _param);
            }
            else
                _amb_manager->addRover(_epoch.sow(), sat_name, _rover_count);
        }
    }
    _amb_manager->get_last_epoch_sats(_epoch.sow());
    assert(cur_sat_prn.size() == _data.size());
}

void gfgomsf::t_gfgo_gins::clearGNSSmsg()
{
    if (_last_gnss_info != nullptr)
        delete _last_gnss_info;
    _last_gnss_info = nullptr;
    cur_sat_prn.clear();
    _DD_msg.clear();
    _vDD_msg.clear();
    _amb_manager->clearState();
    _rover_count = -1;
    _global_sat_id = -1;
    _global_amb_id = -1;
    _amb_state = false;
    memset(_para_amb, 0, sizeof(_para_amb));
    _initial_prior = true;
}

void gfgomsf::t_gfgo_gins::write(string info)
{
    set<string> ambs = _param.amb_prns();
    int nsat = ambs.size();
    if (_Meas_Type.find(MEAS_TYPE::GNSS_MEAS) == _Meas_Type.end())
        nsat = 0;
    // get amb status
    string amb = "Float";
    if (_amb_state)
        amb = "Fixed";
    string meas = "INS";
    if (_Meas_Type.size())
        meas = meas2str(*_Meas_Type.begin());
    ostringstream os;
    os << fixed << _opt_valid;
    int q = prt_ins_kml();
    sins.prt_sins(os);
    os << fixed << setprecision(0) << " " << setw(10) << meas                                   // meas
       << " " << setw(5) << nsat                                                                // nsat
       << fixed << setprecision(2) << " " << setw(7) << (_dop.pdop() < 100 ? _dop.pdop() : 0.0) // pdop
       << fixed << setprecision(2) << " " << setw(8) << amb << setw(8) << (_amb_state ? _ambfix->get_ratio() : 0.0);
    os << endl;
    _fins->write(os.str().c_str(), os.str().size());

    double Xrms = 0.0, Yrms = 0.0, Zrms = 0.0, Vxrms = 0.0, Vyrms = 0.0, Vzrms = 0.0;
    ostringstream osflt("");
    osflt << os.str().substr(0, 104);
    osflt << fixed << setprecision(4) << " " << setw(9) << Xrms // [m]
          << " " << setw(9) << Yrms                             // [m]
          << " " << setw(9) << Zrms                             // [m]
          << " " << setw(9) << Vxrms                            // [m/s]
          << " " << setw(9) << Vyrms                            // [m/s]
          << " " << setw(9) << Vzrms                            // [m/s]
          << fixed << setprecision(0) << " " << setw(5) << nsat // nsat
          << fixed
          << setprecision(1)
          //<< " " << setw(3) << gdop            // gdop
          << fixed << setprecision(2) << " " << setw(5) << _dop.pdop() // pdop
          << fixed << setprecision(2) << " " << setw(8) << 0           // m0
          << " " << setw(8) << amb << setprecision(2) << " " << fixed << setw(10)
          << (_amb_state ? _ambfix->get_ratio() : 0.0);
    osflt << endl;
    if ((!_opt_flag) && _flt)
        _flt->write(osflt.str().c_str(), osflt.str().size());

    os.str("");
    osflt.str("");
}

int gfgomsf::t_gfgo_gins::prt_ins_kml()
{
    Eigen::Vector3d Geo_pos = Eigen::Vector3d(sins.pos(0) / t_gglv::deg, sins.pos(1) / t_gglv::deg, sins.pos(2));
    bool ins = dynamic_cast<t_gsetinp *>(_set)->input_size("imu") > 0 ? true : false;

    double crt = _ins_crt.sow() + _ins_crt.dsec();
    Eigen::Vector3d Qpos = Pk.block(6, 6, 3, 3).diagonal(), Qvel = Pk.block(3, 3, 3, 3).diagonal();
    if (_amb_state && _ign_type == IGN_TYPE::TCI)
    {
        Qpos = _global_variance.block(6, 6, 3, 3).diagonal();
        Qvel = _global_variance.block(3, 3, 3, 3).diagonal();
    }
    Eigen::Vector3d position = sins.pos_ecef, velocity = sins.ve;
    t_gposdata::data_pos posdata =
        t_gposdata::data_pos{crt, position, velocity, Qpos, Qvel, 1.3, int(_data.size()), _amb_state};

    // write kml
    if (_kml && ins)
    {
        ostringstream out;
        out << fixed << setprecision(11) << " " << setw(0) << Geo_pos[1] << ',' << Geo_pos[0];
        string val = out.str();

        xml_node root = _doc;
        xml_node node = this->_default_node(root, _root.c_str());
        xml_node document = node.child("Document");
        xml_node last_child = document.last_child();
        xml_node placemark = document.insert_child_after("Placemark", last_child);
        string q = "#P" + _quality_grade(posdata);
        this->_default_node(placemark, "styleUrl", q.c_str());
        this->_default_node(placemark, "time", int2str(_ins_crt.sow()).c_str());
        xml_node point = this->_default_node(placemark, "Point");
        this->_default_node(point, "coordinates", val.c_str()); // for point
        xml_node description = placemark.append_child("description");
        description.append_child(pugi::node_cdata).set_value(_gen_kml_description(_ins_crt, posdata).c_str());
        xml_node TimeStamp = placemark.append_child("TimeStamp");
        string time = trim(_ins_crt.str_ymd()) + "T" + trim(_ins_crt.str_hms()) + "Z";
        this->_default_node(TimeStamp, "when", time.c_str());

        xml_node Placemark = document.child("Placemark");
        xml_node LineString = Placemark.child("LineString");
        this->_default_node(LineString, "coordinates", val.c_str(), false); // for line
    }

    return 1;
}

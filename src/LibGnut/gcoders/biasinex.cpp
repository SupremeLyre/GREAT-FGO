/* ----------------------------------------------------------------------
 * G-Nut - GNSS software development library
 *
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
  This file is part of the G-Nut C++ library.
-*/

#include <cstring>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>

#include "gall/gallbias.h"
#include "gcoders/biasinex.h"
#include "gmodels/gbias.h"
#include "gutils/gtypeconv.h"

using namespace std;

namespace gnut
{

t_biasinex::t_biasinex(t_gsetbase *s, string version, int sz, string id) : t_sinex(s, version, sz, id)
{

    _allbias = 0;
}

int t_biasinex::_decode_comm()
{
    std::string::size_type idx;
    if ((idx = _line.find("BIAS")) != string::npos)
        _mapidx["TYP"] = make_pair(idx, 4);
    if ((idx = _line.find("PRN")) != string::npos)
        _mapidx["SAT"] = make_pair(idx, 3);
    if ((idx = _line.find("OBS1")) != string::npos)
        _mapidx["OBS1"] = make_pair(idx, 4);
    if ((idx = _line.find("OBS2")) != string::npos)
        _mapidx["OBS2"] = make_pair(idx, 4);

    if ((idx = _line.find("BIAS_START____")) != string::npos)
        _mapidx["BEG"] = make_pair(idx, 14);
    if ((idx = _line.find("BIAS_END______")) != string::npos)
        _mapidx["END"] = make_pair(idx, 14);

    if ((idx = _line.find("__ESTIMATED_VALUE____")) != string::npos)
        _mapidx["EST"] = make_pair(idx, 21);
    if ((idx = _line.find("_STD_DEV___")) != string::npos)
        _mapidx["STD"] = make_pair(idx, 11);
    if ((idx = _line.find("UNIT")) != string::npos)
        _mapidx["UNIT"] = make_pair(idx, 4);

    return 1;
}

int t_biasinex::_decode_block()
{

    t_sinex::_decode_block();

    // -------- "BIAS/DESCRIPTION" --------
    if (_block.find("BIAS/DESCRIPTION") != string::npos)
    {

        if (_line.find(" OBSERVATION_SAMPLING ") != string::npos)
        {
            if (_spdlog)
                SPDLOG_LOGGER_DEBUG(_spdlog, "Read BIAS/SMP: " + cut_crlf(_line.substr(31)));
        }
        else if (_line.find(" PARAMETER_SPACING ") != string::npos)
        {
            if (_spdlog)
                SPDLOG_LOGGER_DEBUG(_spdlog, "Read BIAS/SPC: " + cut_crlf(_line.substr(31)));
        }
        else if (_line.find(" DETERMINATION_METHOD ") != string::npos)
        {
            if (_spdlog)
                SPDLOG_LOGGER_DEBUG(_spdlog, "Read BIAS/MTD: " + cut_crlf(_line.substr(31)));
        }
        else if (_line.find(" BIAS_MODE ") != string::npos)
        {
            if (_spdlog)
                SPDLOG_LOGGER_DEBUG(_spdlog, "Read BIAS/MOD: " + cut_crlf(_line.substr(31)));
            if (_line.find("ABSOLUTE") != string::npos)
            {
                if (_allbias)
                    _allbias->set_osb(true);
                if (_spdlog)
                    SPDLOG_LOGGER_INFO(_spdlog, "BIAS_MODE: ABSOLUTE");
            }
        }
        else if (_line.find(" TIME_SYSTEM ") != string::npos)
        {
            if (_spdlog)
                SPDLOG_LOGGER_DEBUG(_spdlog, "Read BIAS/TSY: " + cut_crlf(_line.substr(31)));
        }
    }
    else if (_block.find("BIAS/SOLUTION") != string::npos)
    {

        string typ = "";
        string prn = "";
        string unit = "ns";
        GOBS gobs1, gobs2;
        gobs1 = gobs2 = X;
        t_gtime beg = FIRST_TIME;
        t_gtime end = LAST_TIME;
        double dcb = 0.0;

        for (auto it = _mapidx.begin(); it != _mapidx.end(); it++)
        {
            size_t pos = it->second.first;
            size_t len = it->second.second;
            if (it->first == "TYP")
            {
                typ = _line.substr(pos, len);
                typ.erase(typ.find_last_not_of(" \n\r\t*") + 1);
                typ.erase(0, typ.find_first_not_of(" \n\r\t*"));
            }
            if (it->first == "SAT")
                prn = _line.substr(pos, len);
            if (it->first == "OBS1")
            {
                string obs_str = _line.substr(pos, len);
                obs_str.erase(obs_str.find_last_not_of(" \n\r\t") + 1);
                obs_str.erase(0, obs_str.find_first_not_of(" \n\r\t"));
                if (!obs_str.empty())
                    gobs1 = str2gobs(obs_str);
            }
            if (it->first == "OBS2")
            {
                string obs_str = _line.substr(pos, len);
                obs_str.erase(obs_str.find_last_not_of(" \n\r\t") + 1);
                obs_str.erase(0, obs_str.find_first_not_of(" \n\r\t"));
                if (!obs_str.empty())
                    gobs2 = str2gobs(obs_str);
            }
            if (it->first == "BEG")
                beg.from_str("%Y:%j:%s", _line.substr(pos, len));
            if (it->first == "END" && _line.substr(pos, len) != "0000:000:00000")
                end.from_str("%Y:%j:%s", _line.substr(pos, len));
            if (it->first == "UNIT")
            {
                unit = _line.substr(pos, len);
                unit.erase(unit.find_last_not_of(" \n\r\t") + 1);
                unit.erase(0, unit.find_first_not_of(" \n\r\t"));
            }
            if (it->first == "EST")
                dcb = str2dbl(_line.substr(pos, len));
        }
        shared_ptr<t_gbias> p_bias;

        if (_allbias)
        {
            if (typ == "OSB")
            {
                _allbias->set_osb(true);
            }

            p_bias = make_shared<t_gbias>(_spdlog);
            p_bias->set(beg, end, dcb, gobs1, gobs2);

            _allbias->add(_ac, beg, prn, p_bias);
        }
    }

    return 1;
}

void t_biasinex::_add_data(const string &id, t_gdata *pt_data)
{
    if (_spdlog)
    {
        SPDLOG_LOGGER_INFO(_spdlog, "Add data to biasinex: {}", id);
    }
    if (pt_data == 0)
        return;

    // ALL OBJECTS
    if (pt_data->id_type() == t_gdata::ALLBIAS)
    {
        if (!_allbias)
        {
            _allbias = dynamic_cast<t_gallbias *>(pt_data);
        }
    }

    return;
}

} // namespace gnut

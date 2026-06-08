/* ----------------------------------------------------------------------
 * G-Nut - GNSS software development library
 *
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
  This file is part of the G-Nut C++ library.
-*/

#include <cstring>
#include <algorithm>
#include <cctype>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>

#include "gall/gallbias.h"
#include "gcoders/biasinex.h"
#include "gmodels/gbias.h"
#include "gutils/gconst.h"
#include "gutils/gsys.h"
#include "gutils/gtypeconv.h"

using namespace std;

namespace gnut
{
namespace
{
string _trim(const string &value)
{
    const string blank = " \n\r\t";
    const size_t first = value.find_first_not_of(blank);
    if (first == string::npos)
        return "";
    const size_t last = value.find_last_not_of(blank);
    return value.substr(first, last - first + 1);
}

string _upper(string value)
{
    transform(value.begin(), value.end(), value.begin(), [](unsigned char c) { return static_cast<char>(toupper(c)); });
    return value;
}

double _frequency_from_prn_obs(const string &prn, const GOBS &gobs)
{
    const string obs = gobs2str(gobs);
    if (prn.empty() || obs.size() < 2)
        return 0.0;

    const GOBSBAND band = str2gobsband(obs);
    switch (prn[0])
    {
    case 'G':
        if (band == BAND_1)
            return G01_F;
        if (band == BAND_2)
            return G02_F;
        if (band == BAND_5)
            return G05_F;
        break;
    case 'E':
        if (band == BAND_1)
            return E01_F;
        if (band == BAND_5)
            return E05_F;
        if (band == BAND_6)
            return E06_F;
        if (band == BAND_7)
            return E07_F;
        if (band == BAND_8)
            return E08_F;
        break;
    case 'C':
        if (band == BAND_1)
            return C01_F;
        if (band == BAND_2)
            return C02_F;
        if (band == BAND_5)
            return C05_F;
        if (band == BAND_6)
            return C06_F;
        if (band == BAND_7)
            return C07_F;
        if (band == BAND_8)
            return C08_F;
        if (band == BAND_9)
            return C09_F;
        break;
    case 'J':
        if (band == BAND_1)
            return J01_F;
        if (band == BAND_2)
            return J02_F;
        if (band == BAND_5)
            return J05_F;
        if (band == BAND_6)
            return J06_F;
        break;
    case 'S':
        if (band == BAND_1)
            return S01_F;
        if (band == BAND_5)
            return S05_F;
        break;
    case 'I':
        if (band == BAND_5)
            return I05_F;
        break;
    default:
        break;
    }

    return 0.0;
}

bool _bias_value_to_meters(const string &unit, const string &prn, const GOBS &gobs, const double &value, double &meters)
{
    const string normalized_unit = _upper(_trim(unit));

    if (normalized_unit.empty() || normalized_unit == "NS" || normalized_unit == "NSEC")
    {
        meters = value * CLIGHT * 1e-9;
        return true;
    }
    if (normalized_unit == "PS" || normalized_unit == "PSEC")
    {
        meters = value * CLIGHT * 1e-12;
        return true;
    }
    if (normalized_unit == "US" || normalized_unit == "USEC")
    {
        meters = value * CLIGHT * 1e-6;
        return true;
    }
    if (normalized_unit == "MS" || normalized_unit == "MSEC")
    {
        meters = value * CLIGHT * 1e-3;
        return true;
    }
    if (normalized_unit == "S" || normalized_unit == "SEC")
    {
        meters = value * CLIGHT;
        return true;
    }
    if (normalized_unit == "M" || normalized_unit == "MET" || normalized_unit == "METER" || normalized_unit == "METERS")
    {
        meters = value;
        return true;
    }
    if (normalized_unit == "CYC" || normalized_unit == "CYCL" || normalized_unit == "CYCLE" ||
        normalized_unit == "CYCLES")
    {
        const double frequency = _frequency_from_prn_obs(prn, gobs);
        if (frequency <= 0.0)
            return false;
        meters = value * CLIGHT / frequency;
        return true;
    }

    return false;
}

} // namespace

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
            double bias_meter = 0.0;
            if (!_bias_value_to_meters(unit, prn, gobs1, dcb, bias_meter))
            {
                if (_spdlog)
                {
                    SPDLOG_LOGGER_WARN(_spdlog, "Skip BIA bias with unsupported UNIT: ac={}, prn={}, obs={}, unit={}",
                                       _ac.c_str(), prn.c_str(), gobs2str(gobs1).c_str(), unit.c_str());
                }
                return 1;
            }

            // t_gbias is the internal product container and stores meters.
            // Observation application decides the observation domain later:
            // code subtracts meters directly; phase subtracts meters/lambda
            // because carrier phase observations are stored in cycles.
            p_bias->set(beg, end, bias_meter, gobs1, gobs2);

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

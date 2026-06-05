/* ----------------------------------------------------------------------
 * G-Nut - GNSS software development library
 *
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
  This file is part of the G-Nut C++ library.
-*/

#include <cstring>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <memory>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gall/gallnav.h"
#include "gall/gallobj.h"
#include "gcoders/rinexn.h"
#include "gdata/geph.h"
#include "gdata/gnav.h"
#include "gdata/grxnhdr.h"
#include "gdata/gtrn.h"
#include "gutils/gfileconv.h"
#include "gutils/gsys.h"
#include "gutils/gtime.h"
#include "gutils/gtypeconv.h"

#include "gdata/gnavbds.h"
#include "gdata/gnavgal.h"
#include "gdata/gnavglo.h"
#include "gdata/gnavgps.h"
#include "gdata/gnavirn.h"
#include "gdata/gnavqzs.h"
#include "gdata/gnavsbs.h"

using namespace std;

namespace gnut
{

namespace
{

int _rinex4_eph_rec_count(char sys, const string &msg_type)
{
    switch (sys)
    {
    case 'G':
        if (msg_type == "LNAV")
            return MAX_RINEXN_REC_GPS;
        if (msg_type == "CNAV")
            return MAX_RINEXN_REC_GPS_CNAV;
        if (msg_type == "CNV2")
            return MAX_RINEXN_REC_GPS_CNV2;
        break;
    case 'E':
        if (msg_type == "INAV" || msg_type == "FNAV")
            return MAX_RINEXN_REC_GAL;
        break;
    case 'R':
        if (msg_type == "FDMA")
            return MAX_RINEXN_REC_GLO;
        if (msg_type == "L1OC")
            return MAX_RINEXN_REC_GLO_L1OC;
        if (msg_type == "L3OC")
            return MAX_RINEXN_REC_GLO_L3OC;
        break;
    case 'C':
        if (msg_type == "D1" || msg_type == "D2")
            return MAX_RINEXN_REC_BDS;
        if (msg_type == "CNV1")
            return MAX_RINEXN_REC_BDS_CNV1;
        if (msg_type == "CNV2")
            return MAX_RINEXN_REC_BDS_CNV2;
        if (msg_type == "CNV3")
            return MAX_RINEXN_REC_BDS_CNV3;
        break;
    case 'J':
        if (msg_type == "LNAV")
            return MAX_RINEXN_REC_QZS;
        if (msg_type == "CNAV")
            return MAX_RINEXN_REC_QZS_CNAV;
        if (msg_type == "CNV2")
            return MAX_RINEXN_REC_QZS_CNV2;
        break;
    case 'S':
        if (msg_type == "SBAS")
            return MAX_RINEXN_REC_SBS;
        break;
    case 'I':
        if (msg_type == "LNAV")
            return MAX_RINEXN_REC_IRN;
        if (msg_type == "L1NV")
            return MAX_RINEXN_REC_IRN_L1NV;
        break;
    default:
        break;
    }
    return -1;
}

int _rinex4_non_eph_line_count(const string &nav_type, char sys, const string &msg_type, const string &msg_subtype)
{
    if (nav_type == "STO")
        return 2;
    if (nav_type == "EOP")
        return 3;
    if (nav_type != "ION")
        return -1;

    if (sys == 'I' && msg_type == "L1NV")
    {
        if (msg_subtype == "NEQN")
            return 7;
        if (msg_subtype == "KLOB")
            return 4;
    }
    if (sys == 'E' && msg_type == "IFNV")
        return 2;
    if (sys == 'R')
        return 1;

    return 3;
}

void _set_rinex4_epoch_tsys(char sys, t_gtime &epoch)
{
    switch (sys)
    {
    case 'R':
        epoch.tsys(t_gtime::UTC);
        break;
    case 'E':
        epoch.tsys(t_gtime::GAL);
        break;
    case 'C':
        epoch.tsys(t_gtime::BDS);
        break;
    default:
        epoch.tsys(t_gtime::GPS);
        break;
    }
}

int _rinex4_week(char sys, const t_gtime &epoch)
{
    return sys == 'C' ? epoch.bwk() : epoch.gwk();
}

bool _rinex4_numeric_field(const string &line, size_t pos, int len, double &value)
{
    if (line.size() <= pos)
        return false;

    string field = trim(line.substr(pos, len));
    if (field.empty())
        return false;

    char *end = nullptr;
    value = strtod(field.c_str(), &end);
    return end != field.c_str();
}

void _append_rinex4_values(const string &line, int start, int fields, vector<double> &data)
{
    for (int i = 0; i < fields && data.size() < MAX_RINEX4_SEI; ++i)
    {
        double value = 0.0;
        if (_rinex4_numeric_field(line, start + i * 19, 19, value))
            data.push_back(value);
    }
}

bool _rinex4_klobuchar_iono(char sys, const string &msg_type, const string &msg_subtype,
                            const vector<double> &data, IONO_CORR &alpha, IONO_CORR &beta,
                            int &alpha_idx, int &beta_idx)
{
    alpha_idx = 0;
    beta_idx = 4;

    if (sys == 'G' && (msg_type == "LNAV" || msg_type == "CNVX"))
    {
        alpha = IO_GPSA;
        beta = IO_GPSB;
    }
    else if (sys == 'C' && msg_type == "D1D2")
    {
        alpha = IO_BDSA;
        beta = IO_BDSB;
    }
    else if (sys == 'J' && (msg_type == "LNAV" || msg_type == "CNVX"))
    {
        alpha = IO_QZSA;
        beta = IO_QZSB;
    }
    else if (sys == 'I' && msg_type == "LNAV")
    {
        alpha = IO_IRNA;
        beta = IO_IRNB;
    }
    else if (sys == 'I' && msg_type == "L1NV" && msg_subtype == "KLOB")
    {
        alpha = IO_IRNA;
        beta = IO_IRNB;
        alpha_idx = 1;
        beta_idx = 5;
    }
    else
    {
        return false;
    }

    return static_cast<int>(data.size()) >= beta_idx + 4;
}

bool _rinex4_default_nav_message(GSYS gs, const string &msg_type)
{
    switch (gs)
    {
    case GPS:
    case QZS:
    case IRN:
        return msg_type == "LNAV";
    case GLO:
        return msg_type == "FDMA";
    case BDS:
        return msg_type == "D1" || msg_type == "D2";
    case SBS:
        return msg_type == "SBAS";
    default:
        return false;
    }
}

bool _rinex4_msg_selected(GSYS gs, const string &msg_type, const set<string> &nav)
{
    if (nav.empty() || nav.find(msg_type) != nav.end())
        return true;

    bool bds_cnav = gs == BDS && (msg_type == "CNV1" || msg_type == "CNV2" || msg_type == "CNV3");
    if (bds_cnav && (nav.find("CNAV") != nav.end() || nav.find("CNV") != nav.end()))
        return true;

    bool legacy = _rinex4_default_nav_message(gs, msg_type);
    if (legacy && nav.find("NAV") != nav.end())
        return true;

    bool bds_d = gs == BDS && (msg_type == "D1" || msg_type == "D2");
    if (bds_d && nav.find("D1D2") != nav.end())
        return true;

    return false;
}

GNAVTYPE _rinex4_msg_gnavtype(const string &msg_type)
{
    if (msg_type == "CNAV")
        return CNAV;
    if (msg_type == "D1")
        return GNAV_D1;
    if (msg_type == "D2")
        return GNAV_D2;
    if (msg_type == "CNV1")
        return GNAV_CNV1;
    if (msg_type == "CNV2")
        return GNAV_CNV2;
    if (msg_type == "CNV3")
        return GNAV_CNV3;
    if (msg_type == "INAV")
        return INAV;
    if (msg_type == "FNAV")
        return FNAV;
    return NAV;
}

void _clear_navdata(t_gnavdata &data)
{
    for (int i = 0; i < MAX_RINEXN_REC; ++i)
        data[i] = 0.0;
}

void _copy_navdata(const t_gnavdata &src, t_gnavdata &dst)
{
    for (int i = 0; i < MAX_RINEXN_REC; ++i)
        dst[i] = src[i];
}

void _copy_adot_kepler_to_legacy(const t_gnavdata &raw, t_gnavdata &data)
{
    data[0] = raw[0];
    data[1] = raw[1];
    data[2] = raw[2];
    data[4] = raw[4];
    data[5] = raw[5];
    data[6] = raw[6];
    data[7] = raw[7];
    data[8] = raw[8];
    data[9] = raw[9];
    data[10] = raw[10];
    data[11] = raw[11];
    data[12] = raw[12];
    data[13] = raw[13];
    data[14] = raw[14];
    data[15] = raw[15];
    data[16] = raw[16];
    data[17] = raw[17];
    data[18] = raw[18];
    data[19] = raw[19];
}

void _normalize_gps_qzs_cnav(const t_gnavdata &raw, t_gnavdata &data, const t_gtime &epoch, int ttm_idx, int week_idx)
{
    _copy_adot_kepler_to_legacy(raw, data);
    data[21] = raw[week_idx] > 0.0 ? raw[week_idx] : epoch.gwk();
    data[23] = raw[23]; // URAI_ED, retained as the closest accuracy indicator.
    data[24] = raw[24];
    data[25] = raw[25];
    data[27] = raw[ttm_idx];
}

void _normalize_bds_cnv12(const t_gnavdata &raw, t_gnavdata &data, const t_gtime &epoch)
{
    _copy_adot_kepler_to_legacy(raw, data);
    data[3] = raw[38];      // IODE
    data[21] = epoch.bwk(); // CNV1/CNV2 carry full Toc; derive BDT week from it.
    data[23] = raw[23];     // SISAI_oe, retained as the closest accuracy indicator.
    data[24] = raw[32];
    data[25] = raw[29];
    data[26] = raw[30];
    data[27] = raw[35];
    data[28] = raw[34]; // IODC
    data[29] = raw[3];  // A DOT
    data[30] = raw[20]; // Delta n0 dot
}

void _normalize_bds_cnv3(const t_gnavdata &raw, t_gnavdata &data, const t_gtime &epoch)
{
    _copy_adot_kepler_to_legacy(raw, data);
    data[21] = epoch.bwk();
    data[23] = raw[23]; // SISAI_oe
    data[24] = raw[28];
    data[25] = raw[30];
    data[27] = raw[31];
    data[29] = raw[3];  // A DOT
    data[30] = raw[20]; // Delta n0 dot
}

void _normalize_glo_cdma(const t_gnavdata &raw, t_gnavdata &data, const t_gtime &epoch)
{
    data[0] = raw[0];
    data[1] = raw[1];
    data[2] = raw[34] > 0.0 ? raw[34] : epoch.sow();

    data[3] = raw[3];
    data[4] = raw[4];
    data[5] = raw[5];
    data[6] = (raw[6] != 0.0 || raw[10] != 0.0) ? 1.0 : 0.0;

    data[7] = raw[7];
    data[8] = raw[8];
    data[9] = raw[9];
    data[10] = 0.0; // CDMA messages do not carry the FDMA frequency channel.

    data[11] = raw[11];
    data[12] = raw[12];
    data[13] = raw[13];
    data[14] = raw[17]; // AODE (EE), closest analogue to age of operation.
}

void _normalize_irn_l1nv(const t_gnavdata &raw, t_gnavdata &data, const t_gtime &epoch)
{
    data[0] = raw[0];
    data[1] = raw[1];
    data[2] = raw[2];
    data[3] = raw[11]; // IODEC
    data[4] = raw[4];
    data[5] = raw[5];
    data[6] = raw[6];
    data[7] = raw[7];
    data[8] = raw[8];
    data[9] = raw[9];
    data[10] = raw[10];
    data[11] = epoch.sow(); // L1NV has no separate Toe field; use Toc per full epoch.
    data[12] = raw[12];
    data[13] = raw[13];
    data[14] = raw[14];
    data[15] = raw[15];
    data[16] = raw[16];
    data[17] = raw[17];
    data[18] = raw[18];
    data[19] = raw[19];
    data[21] = epoch.gwk();
    data[23] = raw[23];
    data[24] = raw[24];
    data[25] = raw[25] != 0.0 ? raw[25] : raw[26];
    data[26] = raw[11];
    data[27] = raw[31];
}

bool _normalize_rinex4_eph(char sys, const string &msg_type, const t_gtime &epoch, t_gnavdata &data)
{
    bool map_adot_cnav = (sys == 'G' || sys == 'J') && (msg_type == "CNAV" || msg_type == "CNV2");
    bool map_bds_cnv12 = sys == 'C' && (msg_type == "CNV1" || msg_type == "CNV2");
    bool map_bds_cnv3 = sys == 'C' && msg_type == "CNV3";
    bool map_glo_cdma = sys == 'R' && (msg_type == "L1OC" || msg_type == "L3OC");
    bool map_irn_l1nv = sys == 'I' && msg_type == "L1NV";

    if (!map_adot_cnav && !map_bds_cnv12 && !map_bds_cnv3 && !map_glo_cdma && !map_irn_l1nv)
        return true;

    t_gnavdata raw;
    _copy_navdata(data, raw);
    _clear_navdata(data);

    if (map_adot_cnav)
        _normalize_gps_qzs_cnav(raw, data, epoch, msg_type == "CNAV" ? 31 : 33, msg_type == "CNAV" ? 32 : 34);
    else if (map_bds_cnv12)
        _normalize_bds_cnv12(raw, data, epoch);
    else if (map_bds_cnv3)
        _normalize_bds_cnv3(raw, data, epoch);
    else if (map_glo_cdma)
        _normalize_glo_cdma(raw, data, epoch);
    else if (map_irn_l1nv)
        _normalize_irn_l1nv(raw, data, epoch);

    return true;
}

} // namespace

t_rinexn::t_rinexn(t_gsetbase *s, string version, int sz) : t_gcoder(s, version, sz)
{
    _gnsssys = 'G';
    _rxnhdr_all = true;
    _check_dt = FIRST_TIME;

    _gset(s); // HAVE TO BE EXPLICITLY CALLED HERE (AT THE END OF CONSTRUCTOR)
}

int t_rinexn::decode_head(char *buff, int sz, vector<string> &errmsg)
{

    _mutex.lock();

    if (t_gcoder::_add2buffer(buff, sz) == 0)
    {
        _mutex.unlock();
        return 0;
    };

    string line;
    int consume = 0;
    int tmpsize = 0;
    while ((tmpsize = t_gcoder::_getline(line)) >= 0)
    {

        consume += tmpsize;

        // -------- "RINEX VERSION" --------
        if (line.find("RINEX VERSION", 60) != string::npos)
        { // first line
            switch (line[20])
            {
            case 'N':
                _gnsssys = 'G';
                break; // Navigation data - according to Rinex specification
            case 'G':
                _gnsssys = 'R';
                break; // GLONASS NAVIGATION - occures sometimes in brdc
            default: {
                string lg("warning: not rinex navigation data file");

                if (_spdlog)
                    SPDLOG_LOGGER_INFO(_spdlog, lg);
                _mutex.unlock();
                return -1;
            }
            }

            switch (line[40])
            {
            case 'G':
                _gnsssys = 'G';
                break; // GPS
            case 'R':
                _gnsssys = 'R';
                break; // GLONASS
            case 'E':
                _gnsssys = 'E';
                break; // GALILEO
            case 'J':
                _gnsssys = 'J';
                break; // QZSS
            case 'S':
                _gnsssys = 'S';
                break; // SBAS
            case 'I':
                _gnsssys = 'I';
                break; // IRNSS
            case 'M':
                _gnsssys = 'M';
                break; // MIXED
            case ' ': {
                if (line[20] == 'N')
                    _gnsssys = 'G';
                if (line[20] == 'G')
                    _gnsssys = 'R';

                if (_spdlog)
                    SPDLOG_LOGGER_WARN(_spdlog,
                                       "warning - RINEXN system not defined, used " + t_gsys::char2str(_gnsssys));
                break;
            }
            default: {
                string lg("warning: not supported satellite system " + line.substr(40, 1));

                if (_spdlog)
                    SPDLOG_LOGGER_INFO(_spdlog, lg);
            }
            }

            _version = trim(line.substr(0, 9));

            _rxnhdr.path(_fname);
            if (_spdlog && substitute(_version, " ", "") > 0)
                if (_spdlog)
                    SPDLOG_LOGGER_INFO(_spdlog, "reading VER: " + _version + " SYS: " + string(1, _gnsssys));
        }
        else if (line.find("PGM / RUN BY / DATE", 60) != string::npos)
        {
            _rxnhdr.program(trim(line.substr(0, 20)));
            _rxnhdr.runby(trim(line.substr(20, 20)));
            t_gtime gtime(t_gtime::UTC);
            if (line.substr(56, 3) != "UTC")
                gtime.tsys(t_gtime::LOC);

            if (gtime.from_str("%Y%m%d %H%M%S", line.substr(40, 15)) == 0)
                ;
            else if (gtime.from_str("%Y-%m-%d %H-%M-%S", line.substr(40, 20)) == 0)
                ;
            else
            {
                gtime = FIRST_TIME;
            }
            _rxnhdr.gtime(gtime);

            if (_spdlog)
                SPDLOG_LOGGER_DEBUG(_spdlog, "PGM / RUN BY / DATE: " + _rxnhdr.program() + " " + _rxnhdr.runby() + " " +
                                                 _rxnhdr.gtime().str_ymdhms());

            // -------- "IONOSPHERIC CORR" --------
        }
        else if (line.find("IONOSPHERIC CORR", 60) != string::npos)
        {

            IONO_CORR IO = str2iono_corr(line.substr(0, 4));
            t_iono_corr io;

            io.x0 = strSci2dbl(line.substr(5 + 0, 12));
            io.x1 = strSci2dbl(line.substr(5 + 12, 12));
            io.x2 = strSci2dbl(line.substr(5 + 24, 12));
            io.x3 = strSci2dbl(line.substr(5 + 36, 12));

            _rxnhdr.iono_corr(IO, io);

            map<string, t_gdata *>::iterator it = _data.begin();
            while (it != _data.end())
            {
                if (it->second->id_type() == t_gdata::ALLNAV || it->second->id_group() == t_gdata::GRP_EPHEM)
                {

                    ((t_gallnav *)it->second)->add_iono_corr(IO, io);
                }
                it++;
            }

            if (_spdlog)
                SPDLOG_LOGGER_DEBUG(_spdlog, "IONOSPHERIC CORR " + iono_corr2str(IO) + dbl2str(io.x0) + dbl2str(io.x1) +
                                                 dbl2str(io.x2) + dbl2str(io.x3) + " " + base_name(_fname));

            // -------- "ION ALPHA" --------
        }
        else if (line.find("ION ALPHA", 60) != string::npos)
        {

            IONO_CORR IO = IO_GPSA;
            t_iono_corr io;

            io.x0 = strSci2dbl(line.substr(2 + 0, 12));
            io.x1 = strSci2dbl(line.substr(2 + 12, 12));
            io.x2 = strSci2dbl(line.substr(2 + 24, 12));
            io.x3 = strSci2dbl(line.substr(2 + 36, 12));

            _rxnhdr.iono_corr(IO, io);
            map<string, t_gdata *>::iterator it = _data.begin();
            while (it != _data.end())
            {
                if (it->second->id_type() == t_gdata::ALLNAV || it->second->id_group() == t_gdata::GRP_EPHEM)
                {

                    ((t_gallnav *)it->second)->add_iono_corr(IO, io);
                }
                it++;
            }

            if (_spdlog)
                SPDLOG_LOGGER_DEBUG(_spdlog, "ION ALPHA " + iono_corr2str(IO) + dbl2str(io.x0) + dbl2str(io.x1) +
                                                 dbl2str(io.x2) + dbl2str(io.x3) + " " + base_name(_fname));

            // -------- "ION BETA" --------
        }
        else if (line.find("ION BETA", 60) != string::npos)
        {

            IONO_CORR IO = IO_GPSB;
            t_iono_corr io;

            io.x0 = strSci2dbl(line.substr(2 + 0, 12));
            io.x1 = strSci2dbl(line.substr(2 + 12, 12));
            io.x2 = strSci2dbl(line.substr(2 + 24, 12));
            io.x3 = strSci2dbl(line.substr(2 + 36, 12));

            _rxnhdr.iono_corr(IO, io);
            map<string, t_gdata *>::iterator it = _data.begin();
            while (it != _data.end())
            {
                if (it->second->id_type() == t_gdata::ALLNAV || it->second->id_group() == t_gdata::GRP_EPHEM)
                {

                    ((t_gallnav *)it->second)->add_iono_corr(IO, io);
                }
                it++;
            }

            if (_spdlog)
                SPDLOG_LOGGER_DEBUG(_spdlog, "ION BETA " + iono_corr2str(IO) + dbl2str(io.x0) + dbl2str(io.x1) +
                                                 dbl2str(io.x2) + dbl2str(io.x3) + " " + base_name(_fname));

            // -------- "TIME SYSTEM CORR" --------
        }
        else if (line.find("TIME SYSTEM CORR", 60) != string::npos)
        {

            TSYS_CORR TS = str2tsys_corr(line.substr(0, 4));
            t_tsys_corr ts;

            if (line.substr(5, 2) != "  ")
            { // ELIMINATE INCORRECT HEADERS FROM SOME RECEIVERS!

                ts.a0 = strSci2dbl(line.substr(5, 17));
                ts.a1 = strSci2dbl(line.substr(22, 16));
                ts.T = str2int(line.substr(38, 7));
                ts.W = str2int(line.substr(45, 5));

                _rxnhdr.tsys_corr(TS, ts);

                if (_spdlog)
                    SPDLOG_LOGGER_DEBUG(_spdlog, "TIME SYSTEM CORR " + tsys_corr2str(TS) + dbl2str(ts.a0 * 1e9, 6) +
                                                     dbl2str(ts.a1 * 1e9, 6) + " " + int2str(ts.T) + " " +
                                                     int2str(ts.W) + " " + base_name(_fname));
            }

            // -------- "DELTA-UTC" --------
        }
        else if (line.find("DELTA-UTC: A0,A1,T,W", 60) != string::npos)
        {

            TSYS_CORR TS = TS_GPUT;
            t_tsys_corr ts;

            ts.a0 = strSci2dbl(line.substr(3, 19));
            ts.a1 = strSci2dbl(line.substr(22, 19));
            ts.T = str2int(line.substr(41, 9));
            ts.W = str2int(line.substr(50, 9));

            _rxnhdr.tsys_corr(TS, ts);

            if (_spdlog)
                SPDLOG_LOGGER_DEBUG(_spdlog, "DELTA-UTC: A0,A1,T,W" + tsys_corr2str(TS) + dbl2str(ts.a0 * 1e9, 6) +
                                                 dbl2str(ts.a1 * 1e9, 6) + " " + int2str(ts.T) + " " + int2str(ts.W) +
                                                 " " + base_name(_fname));

            // -------- "CORR TO SYSTEM" --------
        }
        else if (line.find("CORR TO SYSTEM", 60) != string::npos)
        { // <= RINEX v2.10

            TSYS_CORR TS = TS_GLUT;
            t_tsys_corr ts;

            ts.a0 = strSci2dbl(line.substr(21, 19));
            ts.a1 = 0.0;
            ts.T = 0;
            ts.W = 0;

            _rxnhdr.tsys_corr(TS, ts);

            if (_spdlog)
                SPDLOG_LOGGER_DEBUG(_spdlog, "TIME SYSTEM CORR " + tsys_corr2str(TS) + dbl2str(ts.a0 * 1e9, 6) +
                                                 dbl2str(ts.a1 * 1e9, 6) + " " + int2str(ts.T) + " " + int2str(ts.W) +
                                                 " " + base_name(_fname));

            if (_spdlog)
                SPDLOG_LOGGER_DEBUG(_spdlog, "reading CORR TO SYSTEM");

            // -------- "D-UTC" --------
        }
        else if (line.find("D-UTC A0,A1,T,W,S,U", 60) != string::npos)
        { // == RINEX v2.11

            TSYS_CORR TS = TS_SBUT;
            t_tsys_corr ts;

            ts.a0 = strSci2dbl(line.substr(0, 19));
            ts.a1 = strSci2dbl(line.substr(20, 19));
            ts.T = str2int(line.substr(40, 7));
            ts.W = str2int(line.substr(48, 5));

            _rxnhdr.tsys_corr(TS, ts);

            if (_spdlog)
                SPDLOG_LOGGER_DEBUG(_spdlog, "D-UTC: A0,A1,T,W,S,U" + tsys_corr2str(TS) + dbl2str(ts.a0 * 1e9, 6) +
                                                 dbl2str(ts.a1 * 1e9, 6) + " " + int2str(ts.T) + " " + int2str(ts.W) +
                                                 " " + base_name(_fname));

            // -------- "LEAP SECONDS" --------
        }
        else if (line.find("LEAP SECONDS", 60) != string::npos)
        {
            _rxnhdr.leapsec(str2int(line.substr(0, 6)));

            if (_spdlog)
                SPDLOG_LOGGER_DEBUG(_spdlog, "reading LEAP SECONDS");

            // -------- "COMMENT" --------
        }
        else if (line.find("COMMENT", 60) != string::npos || line.find("DOI", 60) != string::npos ||
                 line.find("MERGED FILE", 60) != string::npos)
        {
            _comment.push_back(line.substr(0, 60));

            if (_spdlog)
                SPDLOG_LOGGER_DEBUG(_spdlog, "reading COMMENT");

            // -------- "END OF HEADER" --------
        }
        else if (line.find("END OF HEADER", 60) != string::npos)
        {
            _fill_head();

            if (_spdlog)
                SPDLOG_LOGGER_DEBUG(_spdlog, "reading END OF HEADER ");
            t_gcoder::_consume(tmpsize);
            _mutex.unlock();
            return -1;
        }
        else
        {
            // t_gcoder::_consume(tmpsize);
            _mutex.unlock();
            return -1;
        }
        t_gcoder::_consume(tmpsize);
    }

    _mutex.unlock();
    return consume;
}

int t_rinexn::decode_data(char *buff, int sz, int &cnt, vector<string> &errmsg)
{

    _mutex.lock();

    if (t_gcoder::_add2buffer(buff, sz) == 0)
    {
        _mutex.unlock();
        return 0;
    };

    t_gtime epoch;
    int b, e, s, l, i;
    int maxrec = MAX_RINEXN_REC; // implicite
    string timstr;
    t_gnavdata data = {0};

    // RINEX v2.xx
    if (_version[0] == '2')
    {
        b = 0;
        e = 22;
        l = 19;
        s = 3; // timstr = "%2s %2d %02d %02d %02d %02d %02d";
               // RINEX v3.xx
    }
    else if (_version[0] == '3')
    {
        b = 0;
        e = 23;
        l = 19;
        s = 4; // timstr = "%3s %4d %02d %02d %02d %02d %02d";
               // RINEX ???
    }
    else
    {
        b = 0;
        e = 23;
        l = 19;
        s = 4; // timstr = "%3s %4d %02d %02d %02d %02d %02d";
    }

    string line;
    int consume = 0;
    int tmpsize = 0;
    int recsize = 0;
    string msg_type_str = "";
    string msg_subtype_str = "";

    while ((tmpsize = t_gcoder::_getline(line, 0)) >= 0)
    {

        consume += tmpsize;
        recsize += tmpsize;

        string prn;
        int yr, mn, dd, hr, mi;
        int svn = 0;
        int irc = 0;
        float sec = 0.0;
        char tmpbuff[83]; // RINEX 2 (80) + RINEX 3 (81)
        int min_sz = 22;  // minimum size to be decoded (timestamp)
        maxrec = MAX_RINEXN_REC;
        _clear_navdata(data);

        switch (_version[0])
        {
        case '2':
            min_sz = 22;
            break; // RINEX 2
        case '3':
            min_sz = 23;
            break; // RINEX 3
        case '4':
            min_sz = 14;
            break; // RINEX 4
        default:
            break;
        }

        if (line.size() > 82 || _decode_buffer.size() <= min_sz)
        {
            t_gcoder::_consume(tmpsize);
            recsize = consume = 0;
            break; // read buffer
        }

        int copy_sz = line.size() < static_cast<size_t>(min_sz) ? static_cast<int>(line.size()) : min_sz;
        strncpy(tmpbuff, line.c_str(), copy_sz);
        tmpbuff[copy_sz] = '\0';

        if (_version[0] == '2')
        {
            irc = sscanf(tmpbuff, "%2d%*[ ] %2d%*[ ] %2d%*[ ] %2d%*[ ] %2d%*[ ] %2d%*[ ] %5f", &svn, &yr, &mn, &dd, &hr,
                         &mi, &sec);

            if (irc < 7)
            {
                t_gcoder::_consume(tmpsize);
                recsize = consume = 0;
                continue;
            }
            prn = t_gsys::eval_sat(svn, t_gsys::char2gsys(_gnsssys));
        }
        else if (_version[0] == '3')
        {
            int isec = 0;
            char sat[3 + 1];
            sat[3] = '\0'; // don't use '%2i' in scan, but '%2d' instead !
            irc = sscanf(tmpbuff, "%c%2d%*[ ] %4d%*[ ] %2d%*[ ] %2d%*[ ]%2d%*[ ]%2d%*[ ]%2d", &sat[0], &svn, &yr, &mn,
                         &dd, &hr, &mi, &isec);

            if (irc < 8)
            {
                t_gcoder::_consume(tmpsize);
                recsize = consume = 0;
                continue;
            }
            prn = t_gsys::eval_sat(svn, t_gsys::char2gsys(sat[0]));

            sec = isec;
        }
        else if (_version[0] == '4')
        {
            if (line.empty() || line[0] != '>')
            {
                t_gcoder::_consume(tmpsize);
                recsize = consume = 0;
                continue;
            }

            char marker = '\0';
            string s_nav_type;
            string sat_id;
            istringstream hdr(line);
            hdr >> marker >> s_nav_type >> sat_id >> msg_type_str;
            if (marker != '>' || s_nav_type.empty() || sat_id.empty() || msg_type_str.empty())
            {
                t_gcoder::_consume(tmpsize);
                recsize = consume = 0;
                continue;
            }

            msg_subtype_str = "";
            hdr >> msg_subtype_str;
            if (s_nav_type == "EPH")
            {
                tmpsize = t_gcoder::_getline(line, recsize);
                if (tmpsize < 0)
                    break;

                consume += tmpsize;
                recsize += tmpsize;

                int isec = 0;
                char sat[3 + 1];
                sat[3] = '\0';

                // For EPH, the second line format matches RINEX 3 exactly
                irc = sscanf(line.c_str(), "%c%2d%*[ ] %4d%*[ ] %2d%*[ ]%2d%*[ ]%2d%*[ ]%2d%*[ ]%2d", &sat[0], &svn,
                             &yr, &mn, &dd, &hr, &mi, &isec);

                if (irc < 8)
                {
                    t_gcoder::_consume(recsize);
                    recsize = consume = 0;
                    continue;
                }
                prn = t_gsys::eval_sat(svn, t_gsys::char2gsys(sat[0]));
                sec = isec;
            }
            else
            {
                int non_eph = _decode_v4_non_eph(s_nav_type, sat_id, msg_type_str, msg_subtype_str, recsize, consume);
                if (non_eph < 0)
                {
                    break;
                }
                cnt += non_eph;
                continue;
            }
        }

        shared_ptr<t_gnav> geph = make_shared<t_gnav>(_spdlog);

        if (prn[0] == 'G')
        {
            maxrec = (_version[0] == '4') ? _rinex4_eph_rec_count(prn[0], msg_type_str) : MAX_RINEXN_REC_GPS;
            geph = make_shared<t_gnavgps>(_spdlog);
            geph->spdlog(_spdlog);
            epoch.tsys(t_gtime::GPS);
        }
        else if (prn[0] == 'R')
        {
            maxrec = (_version[0] == '4') ? _rinex4_eph_rec_count(prn[0], msg_type_str) : MAX_RINEXN_REC_GLO;
            geph = make_shared<t_gnavglo>(_spdlog);
            geph->spdlog(_spdlog);
            epoch.tsys(t_gtime::UTC);
        }
        else if (prn[0] == 'E')
        {
            maxrec = (_version[0] == '4') ? _rinex4_eph_rec_count(prn[0], msg_type_str) : MAX_RINEXN_REC_GAL;
            geph = make_shared<t_gnavgal>(_spdlog);
            geph->spdlog(_spdlog);
            epoch.tsys(t_gtime::GAL);
        }
        else if (prn[0] == 'J')
        {
            maxrec = (_version[0] == '4') ? _rinex4_eph_rec_count(prn[0], msg_type_str) : MAX_RINEXN_REC_QZS;
            geph = make_shared<t_gnavqzs>(_spdlog);
            geph->spdlog(_spdlog);
            epoch.tsys(t_gtime::GPS); // QZSS is equal to GPS time
        }
        else if (prn[0] == 'S')
        {
            maxrec = (_version[0] == '4') ? _rinex4_eph_rec_count(prn[0], msg_type_str) : MAX_RINEXN_REC_SBS;
            geph = make_shared<t_gnavsbs>(_spdlog);
            geph->spdlog(_spdlog);
            epoch.tsys(t_gtime::GPS); // SBAS is equal to GPS time
        }
        else if (prn[0] == 'C')
        {
            maxrec = (_version[0] == '4') ? _rinex4_eph_rec_count(prn[0], msg_type_str) : MAX_RINEXN_REC_BDS;
            geph = make_shared<t_gnavbds>(_spdlog);
            geph->spdlog(_spdlog);
            epoch.tsys(t_gtime::BDS);
        }
        else if (prn[0] == 'I')
        {
            maxrec = (_version[0] == '4') ? _rinex4_eph_rec_count(prn[0], msg_type_str) : MAX_RINEXN_REC_IRN;
            geph = make_shared<t_gnavirn>(_spdlog);
            geph->spdlog(_spdlog);
            epoch.tsys(t_gtime::GPS); // ??
        }
        else
        {
            string lg("Warning: not supported satellite satellite system: " + prn);
            if (_spdlog)
                SPDLOG_LOGGER_WARN(_spdlog, lg);
        }

        if (_version[0] == '4' && maxrec < 0)
        {
            if (_spdlog)
                SPDLOG_LOGGER_WARN(_spdlog, "skip unsupported RINEX 4 EPH message type: " + prn + " " + msg_type_str);
            t_gcoder::_consume(recsize);
            recsize = consume = 0;
            continue;
        }

        epoch.from_ymdhms(yr, mn, dd, hr, mi, sec);

        if (fabs(epoch - _check_dt) > 7 * 86400 && _check_dt != FIRST_TIME)
        {
            string lg(prn + " strange epoch [" + epoch.str_ymdhms() + "] or corrupted file [" + base_name(_fname) +
                      "]");
            mesg(GWARNING, lg);
            t_gcoder::_consume(tmpsize);
            recsize = consume = 0;
            continue;
        }

        if (tmpsize < 57 + s)
            break;

        data[0] = strSci2dbl(line.substr(19 + s, l));
        data[1] = strSci2dbl(line.substr(38 + s, l));
        data[2] = strSci2dbl(line.substr(57 + s, l));

        i = 2;
        while (i < MAX_RINEXN_REC)
        {

            // incomplete record
            if ((tmpsize = t_gcoder::_getline(line, recsize)) < 0)
            {
                break;
            }

            consume += tmpsize;
            recsize += tmpsize;
            if (++i < maxrec)
            {
                if (tmpsize > s)
                    data[i] = strSci2dbl(line.substr(s, l));
            }
            if (++i < maxrec)
            {
                if (tmpsize > 19 + s)
                    data[i] = strSci2dbl(line.substr(19 + s, l));
            }
            if (++i < maxrec)
            {
                if (tmpsize > 38 + s)
                    data[i] = strSci2dbl(line.substr(38 + s, l));
            }
            if (++i < maxrec)
            {
                if (tmpsize > 57 + s)
                    data[i] = strSci2dbl(line.substr(57 + s, l));
            }

            // is record complete and filter-out GNSS systems
            if (geph && i + 1 >= maxrec)
            {
                t_gcoder::_consume(recsize);
                recsize = 0;

                // filter GNSS and SAT
                if (!_filter_gnss(prn))
                {

                    if (_spdlog)
                        SPDLOG_LOGGER_DEBUG(_spdlog, "skip " + prn);
                    break;
                }

                if (epoch < _beg - MAX_NAV_TIMEDIFF || epoch > _end + MAX_NAV_TIMEDIFF)
                {

                    if (_spdlog)
                        SPDLOG_LOGGER_DEBUG(_spdlog, "skip " + prn + " " + epoch.str_ymdhms());
                    break;
                }
                if (_version[0] == '4' && !_normalize_rinex4_eph(prn[0], msg_type_str, epoch, data))
                    break;

                if (_version[0] == '4' && prn[0] == 'C')
                    dynamic_pointer_cast<t_gnavbds>(geph)->gnavtype(_rinex4_msg_gnavtype(msg_type_str));

                geph->data2nav(prn, epoch, data);
                geph->gio(_gio_ptr.lock());

                // reset check_dt (already filtered sat)
                if (geph->healthy())
                {
                    _check_dt = geph->epoch();
                };

                if (!_filter_gnav(geph, prn, msg_type_str))
                {

                    if (_spdlog)
                        SPDLOG_LOGGER_DEBUG(_spdlog,
                                            "skip " + prn + " nav type [" + gnavtype2str(geph->gnavtype()) + "]");
                    break;
                }

                // collect 1-line messages
                if (_spdlog)
                    SPDLOG_LOGGER_DEBUG(_spdlog, geph->linefmt() + " " + base_name(_fname));

                // fill data
                map<string, t_gdata *>::iterator it = _data.begin();
                while (it != _data.end())
                {

                    if (it->second->id_type() == t_gdata::ALLNAV || it->second->id_group() == t_gdata::GRP_EPHEM)
                    {
                        ((t_gallnav *)it->second)->add(geph);
                    }
                    else if (it->second->id_type() == t_gdata::ALLOBJ)
                    {
                        t_gallobj *all_obj = (t_gallobj *)it->second;
                        shared_ptr<t_gobj> one_obj = all_obj->obj(geph->sat());
                        if (one_obj != 0)
                            one_obj->spdlog(_spdlog);
                        if (geph->id_type() == t_gdata::EPHGLO)
                        {
                            int ch = dynamic_pointer_cast<t_gnavglo>(geph)->channel();
                            if (one_obj != 0)
                                one_obj->channel(ch);
                        }
                    }
                    it++;
                }
                cnt++;
                break;
            }
        }
        if (recsize != 0)
            break; // break the initialization loop if read not finished correctly
    }
    _mutex.unlock();
    return consume;
}

int t_rinexn::_fill_head()
{

    int cnt = 0;

    for (auto itDAT = _data.begin(); itDAT != _data.end(); ++itDAT)
    {
        if (itDAT->second->id_type() == t_gdata::ALLOBJ)
        {
            t_gallobj *all_obj = (t_gallobj *)itDAT->second;

            shared_ptr<t_gtrn> obj = dynamic_pointer_cast<t_gtrn>(all_obj->obj(ID_ALLSAT));

            if (obj == 0)
            {
                obj = make_shared<t_gtrn>(_spdlog);
                obj->id(ID_ALLSAT);
                obj->name(ID_ALLSAT);
                obj->spdlog(_spdlog);
                obj->header(_rxnhdr, _fname);
                all_obj->add(obj);

                if (_spdlog)
                    SPDLOG_LOGGER_DEBUG(_spdlog, "Object created: " + obj->id());
            }
            else
            {

                if (_spdlog)
                    SPDLOG_LOGGER_DEBUG(_spdlog, "Object completed: " + obj->id());
                obj->header(_rxnhdr, _fname);
            }
            ++cnt;
        }
    }

    return cnt;
}

bool t_rinexn::_filter_gnav(shared_ptr<t_gnav> geph, const string &prn, const string &msg_type)
{

    bool ret = true;

    GSYS gs = t_gsys::char2gsys(prn[0]);
    if (_nav[gs].empty())
        return true;

    if (_version[0] == '4' && !msg_type.empty())
        return _rinex4_msg_selected(gs, msg_type, _nav[gs]);

    // Galileo has source-specific INAV/FNAV handling in existing navigation classes.
    if (gs == GAL)
    {

        GNAVTYPE gnav;
        if (_nav[gs].find(gnavtype2str(INAV)) != _nav[gs].end())
        {
            gnav = dynamic_pointer_cast<t_gnavgal>(geph)->gnavtype(false);
        } // any-source fit in INAV request
        else
        {
            gnav = dynamic_pointer_cast<t_gnavgal>(geph)->gnavtype(true);
        } // exact fit for NAV source

        if (_nav[gs].find(gnavtype2str(gnav)) != _nav[gs].end())
        {
            ret = true;
        }
        else
        {
            ret = false;
        }
    }
    else if (_nav[gs].find(gnavtype2str(geph->gnavtype(true))) == _nav[gs].end() &&
             _nav[gs].find(gnavtype2str(geph->gnavtype(false))) == _nav[gs].end())
    {
        ret = false;
    }

    return ret; // ALL OTHER SYSTEMS
}

int t_rinexn::_decode_v4_non_eph(string nav_type, string sat_id, string msg_type, string msg_subtype, int &recsize, int &consume)
{
    string line;
    int tmpsize;

    char sys = sat_id.empty() ? '\0' : sat_id[0];
    int expected_lines = _rinex4_non_eph_line_count(nav_type, sys, msg_type, msg_subtype);
    if (expected_lines < 0)
    {
        if (_spdlog)
            SPDLOG_LOGGER_WARN(_spdlog, "skip unsupported RINEX 4 " + nav_type + " message type: " + sat_id + " " + msg_type + " " + msg_subtype);

        for (int i = 0; i < MAX_RINEX4_SEI; ++i)
        {
            tmpsize = t_gcoder::_getline(line, recsize);
            if (tmpsize < 0 || (!line.empty() && line[0] == '>'))
                break;

            consume += tmpsize;
            recsize += tmpsize;
        }
        t_gcoder::_consume(recsize);
        recsize = 0;
        return 0;
    }

    t_rinex4_navmsg msg;
    msg.nav_type = nav_type;
    msg.sat = sat_id;
    msg.msg_type = msg_type;
    msg.msg_subtype = msg_subtype;

    for (int line_idx = 0; line_idx < expected_lines; ++line_idx)
    {
        tmpsize = t_gcoder::_getline(line, recsize);
        if (tmpsize < 0)
            return -1;

        consume += tmpsize;
        recsize += tmpsize;

        if (line_idx == 0)
        {
            int yr = 0, mn = 0, dd = 0, hr = 0, mi = 0, isec = 0;
            int irc = sscanf(line.c_str(), "%4d %02d %02d %02d %02d %02d", &yr, &mn, &dd, &hr, &mi, &isec);
            if (irc < 6)
            {
                t_gcoder::_consume(recsize);
                recsize = 0;
                return 0;
            }

            _set_rinex4_epoch_tsys(sys, msg.epoch);
            msg.epoch.from_ymdhms(yr, mn, dd, hr, mi, isec);

            if (nav_type == "STO")
            {
                string tail = line.size() > 23 ? trim(line.substr(23)) : "";
                istringstream sto_info(tail);
                sto_info >> msg.corr_type >> msg.time_ref;
            }
            else
            {
                _append_rinex4_values(line, 23, 3, msg.data);
            }
        }
        else
        {
            _append_rinex4_values(line, 4, 4, msg.data);
        }
    }

    _rxnhdr.rinex4_navmsg(msg);

    if (nav_type == "STO" && msg.data.size() >= 3)
    {
        TSYS_CORR TS = str2tsys_corr(msg.corr_type);
        if (TS != TS_NONE)
        {
            t_tsys_corr ts;
            ts.a0 = msg.data[1];
            ts.a1 = msg.data[2];
            ts.T = static_cast<int>(msg.data[0]);
            ts.W = _rinex4_week(sys, msg.epoch);
            _rxnhdr.tsys_corr(TS, ts);
        }
    }
    else if (nav_type == "ION")
    {
        auto store_iono = [this](IONO_CORR IO, const t_iono_corr &io) {
            _rxnhdr.iono_corr(IO, io);
            map<string, t_gdata *>::iterator it = _data.begin();
            while (it != _data.end())
            {
                if (it->second->id_type() == t_gdata::ALLNAV || it->second->id_group() == t_gdata::GRP_EPHEM)
                    ((t_gallnav *)it->second)->add_iono_corr(IO, io);
                it++;
            }
        };

        IONO_CORR alpha = IO_NONE;
        IONO_CORR beta = IO_NONE;
        int alpha_idx = 0;
        int beta_idx = 0;
        if (_rinex4_klobuchar_iono(sys, msg_type, msg_subtype, msg.data, alpha, beta, alpha_idx, beta_idx))
        {
            t_iono_corr io_alpha;
            io_alpha.x0 = msg.data[alpha_idx + 0];
            io_alpha.x1 = msg.data[alpha_idx + 1];
            io_alpha.x2 = msg.data[alpha_idx + 2];
            io_alpha.x3 = msg.data[alpha_idx + 3];
            io_alpha.T = msg.epoch.sow();
            io_alpha.sat = sat_id.size() > 1 ? str2int(sat_id.substr(1)) : 0;
            store_iono(alpha, io_alpha);

            t_iono_corr io_beta;
            io_beta.x0 = msg.data[beta_idx + 0];
            io_beta.x1 = msg.data[beta_idx + 1];
            io_beta.x2 = msg.data[beta_idx + 2];
            io_beta.x3 = msg.data[beta_idx + 3];
            io_beta.T = msg.epoch.sow();
            io_beta.sat = io_alpha.sat;
            store_iono(beta, io_beta);
        }
        else if (sys == 'E' && msg_type == "IFNV" && msg.data.size() >= 3)
        {
            t_iono_corr io;
            io.x0 = msg.data[0];
            io.x1 = msg.data[1];
            io.x2 = msg.data[2];
            io.x3 = msg.data.size() > 3 ? msg.data[3] : 0.0;
            io.T = msg.epoch.sow();
            io.sat = sat_id.size() > 1 ? str2int(sat_id.substr(1)) : 0;
            store_iono(IO_GAL, io);
        }
    }

    if (_spdlog)
        SPDLOG_LOGGER_DEBUG(_spdlog, "RINEX 4 " + nav_type + " " + sat_id + " " + msg_type + " " + msg_subtype + " " + base_name(_fname));

    _fill_head();

    t_gcoder::_consume(recsize);
    recsize = 0;
    return 1;
}

} // namespace gnut

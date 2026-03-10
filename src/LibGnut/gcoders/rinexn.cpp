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

        strncpy(tmpbuff, line.c_str(), min_sz);
        tmpbuff[line.size()] = '\0';

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
            if (line[0] != '>')
            {
                t_gcoder::_consume(tmpsize);
                recsize = consume = 0;
                continue;
            }

            char nav_type[4];
            char sat_id[4];
            char msg_type[5];
            irc = sscanf(tmpbuff, "> %3s %3s %4s", nav_type, sat_id, msg_type);
            if (irc < 3)
            {
                t_gcoder::_consume(tmpsize);
                recsize = consume = 0;
                continue;
            }

            msg_type_str = string(msg_type);
            string s_nav_type(nav_type);
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
                if (_decode_v4_non_eph(s_nav_type, sat_id, msg_type, recsize, consume) < 0)
                {
                    break;
                }
                continue;
            }
        }

        shared_ptr<t_gnav> geph = make_shared<t_gnav>(_spdlog);

        if (prn[0] == 'G')
        {
            if (msg_type_str == "CNAV") {
                maxrec = MAX_RINEXN_REC_GPS_CNAV;
            } else if (msg_type_str == "CNV2") {
                maxrec = MAX_RINEXN_REC_GPS_CNV2;
            } else {
                maxrec = MAX_RINEXN_REC_GPS;
            }
            geph = make_shared<t_gnavgps>(_spdlog);
            geph->spdlog(_spdlog);
            epoch.tsys(t_gtime::GPS);
        }
        else if (prn[0] == 'R')
        {
            if (msg_type_str == "L1OC") {
                maxrec = MAX_RINEXN_REC_GLO_L1OC;
            } else if (msg_type_str == "L3OC") {
                maxrec = MAX_RINEXN_REC_GLO_L3OC;
            } else {
                maxrec = MAX_RINEXN_REC_GLO;
            }
            geph = make_shared<t_gnavglo>(_spdlog);
            geph->spdlog(_spdlog);
            epoch.tsys(t_gtime::UTC);
        }
        else if (prn[0] == 'E')
        {
            maxrec = MAX_RINEXN_REC_GAL;
            geph = make_shared<t_gnavgal>(_spdlog);
            geph->spdlog(_spdlog);
            epoch.tsys(t_gtime::GAL);
        }
        else if (prn[0] == 'J')
        {
            if (msg_type_str == "CNAV") {
                maxrec = MAX_RINEXN_REC_QZS_CNAV;
            } else if (msg_type_str == "CNV2") {
                maxrec = MAX_RINEXN_REC_QZS_CNV2;
            } else {
                maxrec = MAX_RINEXN_REC_QZS;
            }
            geph = make_shared<t_gnavqzs>(_spdlog);
            geph->spdlog(_spdlog);
            epoch.tsys(t_gtime::GPS); // QZSS is equal to GPS time
        }
        else if (prn[0] == 'S')
        {
            maxrec = MAX_RINEXN_REC_SBS;
            geph = make_shared<t_gnavsbs>(_spdlog);
            geph->spdlog(_spdlog);
            epoch.tsys(t_gtime::GPS); // SBAS is equal to GPS time
        }
        else if (prn[0] == 'C')
        {
            if (msg_type_str == "CNV1") {
                maxrec = MAX_RINEXN_REC_BDS_CNV1;
            } else if (msg_type_str == "CNV2") {
                maxrec = MAX_RINEXN_REC_BDS_CNV2;
            } else if (msg_type_str == "CNV3") {
                maxrec = MAX_RINEXN_REC_BDS_CNV3;
            } else {
                maxrec = MAX_RINEXN_REC_BDS;
            }
            geph = make_shared<t_gnavbds>(_spdlog);
            geph->spdlog(_spdlog);
            epoch.tsys(t_gtime::BDS);
        }
        else if (prn[0] == 'I')
        {
            if(msg_type_str == "L1NV") {
                maxrec = MAX_RINEXN_REC_IRN_L1NV;
            } else {
                maxrec = MAX_RINEXN_REC_IRN;
            }
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
                // only process default EPH for RINEX v4, but allow all for other versions (for backward compatibility)
                if (_version[0] == '4')
                {
                    bool allowed = false;
                    if (prn[0] == 'G' && msg_type_str == "LNAV") allowed = true;
                    else if (prn[0] == 'E' && (msg_type_str == "INAV" || msg_type_str == "FNAV")) allowed = true;
                    else if (prn[0] == 'R' && msg_type_str == "FDMA") allowed = true;
                    else if (prn[0] == 'C' && (msg_type_str == "D1" || msg_type_str == "D2")) allowed = true;
                    else if (prn[0] == 'J' && msg_type_str == "LNAV") allowed = true;
                    else if (prn[0] == 'I' && msg_type_str == "LNAV") allowed = true;
                    else if (prn[0] == 'S') allowed = true;

                    if (!allowed)
                    {
                        break;
                    }
                }

                geph->data2nav(prn, epoch, data);
                geph->gio(_gio_ptr.lock());

                // reset check_dt (already filtered sat)
                if (geph->healthy())
                {
                    _check_dt = geph->epoch();
                };

                if (!_filter_gnav(geph, prn))
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

bool t_rinexn::_filter_gnav(shared_ptr<t_gnav> geph, const string &prn)
{

    bool ret = true;

    GSYS gs = t_gsys::char2gsys(prn[0]);

    // currently only GALILEO can be filtered!
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

        if (_nav[gs].size() == 0 || _nav[gs].find(gnavtype2str(gnav)) != _nav[gs].end())
        {
            ret = true;
        }
        else
        {
            ret = false;
        }
    }

    return ret; // ALL OTHER SYSTEMS
}

int t_rinexn::_decode_v4_non_eph(string nav_type, string sat_id, string msg_type, int &recsize, int &consume)
{
    string line;
    int tmpsize;

    // Read the time line (second line)
    tmpsize = t_gcoder::_getline(line, recsize);
    if (tmpsize < 0)
        return -1;

    consume += tmpsize;
    recsize += tmpsize;

    int yr, mn, dd, hr, mi;
    int isec = 0;
    // other message format "    2025 03 18 01 20 00 BDGA"
    int irc = sscanf(line.c_str(), "%4d %02d %02d %02d %02d %02d", &yr, &mn, &dd, &hr, &mi, &isec);

    t_gtime epoch;
    if (irc >= 6)
    {
        epoch.tsys(t_gtime::GPS); // Assume GPS time for now
        epoch.from_ymdhms(yr, mn, dd, hr, mi, isec);
    }

    // Read up to 6 more lines (7 lines total counting the time line, max 27 data)
    int s = 4;
    int l = 19;
    double data[27] = {0};
    int data_idx = 0;
    int lines_read = 0;
    if (tmpsize < 57 + s)
    {
    }
    else
    {
        data[data_idx++] = strSci2dbl(line.substr(19 + s, l));
        data[data_idx++] = strSci2dbl(line.substr(38 + s, l));
        data[data_idx++] = strSci2dbl(line.substr(57 + s, l));
    }

    while (data_idx < 27 && lines_read < 6)
    {
        tmpsize = t_gcoder::_getline(line, recsize);
        if (tmpsize < 0)
            return -1;

        if (line.size() > 0 && line[0] == '>')
        {
            // Reached next block safely
            break;
        }

        consume += tmpsize;
        recsize += tmpsize;
        lines_read++;

        if (tmpsize > s)
            data[data_idx++] = strSci2dbl(line.substr(s, l));
        if (tmpsize > 19 + s && data_idx < 39)
            data[data_idx++] = strSci2dbl(line.substr(19 + s, l));
        if (tmpsize > 38 + s && data_idx < 39)
            data[data_idx++] = strSci2dbl(line.substr(38 + s, l));
        if (tmpsize > 57 + s && data_idx < 39)
            data[data_idx++] = strSci2dbl(line.substr(57 + s, l));
    }

    t_gcoder::_consume(recsize);
    recsize = 0;
    return 0;
}

} // namespace gnut

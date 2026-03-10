#include "odofile.h"
#include "gdata/gododata.h"
#include "gset/gsetign.h"
using namespace Eigen;

great::t_odofile::t_odofile(t_gsetbase *s, string version, int sz) : t_gcoder(s, version, sz)
{
    _info = dynamic_cast<t_gsetign *>(s)->wheelraduis() * 6.2832;
    _format = dynamic_cast<t_gsetign *>(s)->odo();
    _ts = dynamic_cast<t_gsetign *>(s)->TS_odo();
    _sign = dynamic_cast<t_gsetign *>(s)->odo_inst() == "Right" ? 1 : -1;
}

int great::t_odofile::decode_head(char *buff, int sz, vector<string> &errmsg)
{
#ifdef BMUTEX
    boost::mutex::scoped_lock lock(_mutex);
#endif
    _mutex.lock();

    // no header expected, but fill the buffer
    t_gcoder::_add2buffer(buff, sz);
    _mutex.unlock();
    return -1;
}

int great::t_odofile::decode_data(char *buff, int sz, int &cnt, vector<string> &errmsg)
{

    _mutex.lock();

    if (t_gcoder::_add2buffer(buff, sz) == 0)
    {
        _mutex.unlock();
        return 0;
    };

    int tmpsize = 0;

    string line;
    while ((tmpsize = t_gcoder::_getline(line, 0)) >= 0)
    {
        if (line[0] == '#' || line[0] == 'T')
        {
            t_gcoder::_consume(tmpsize);
            continue;
        }

        for (int i = 0; i < line.size(); i++)
        {
            if (line[i] == ',' || line[i] == '*' || line[i] == ';')
                line[i] = ' ';
        }

        stringstream ss(line);
        string tmp1;
        double tmp2;
        double _tt, val;
        int flag;

        if (_format == "Raw")
        {
            ss >> tmp1 >> tmp2 >> _tt;
            for (int i = 4; i <= 10; i++)
                ss >> tmp2;
            ss >> val >> tmp2 >> flag;
            // fprintf(_f_out, "%15.5lf %15.5lf %15.5lf %d\n", _tt, val, tmp2, flag);
            val = val * _info * (flag == 0 ? _sign : -_sign) / _ts;
        }
        else if (_format == "Pulse")
        {
            ss >> _tt >> val >> tmp2 >> flag;
            val = val * _info * (flag == 0 ? _sign : -_sign) / _ts;
            // fprintf(_f_out, "%15.5lf%15.6lf\n", _tt, val);
        }
        else if (_format == "Velocity")
        {
            ss >> _tt >> val;
        }

        map<string, t_gdata *>::iterator it = _data.begin();
        while (it != _data.end())
        {
            if (it->second->id_type() == t_gdata::ODODATA)
            {
                ((t_gododata *)it->second)->add_odo(_tt, -val);
            }
            ++it;
        }
        if (ss.fail())
        {
            if (_spdlog)
                SPDLOG_LOGGER_DEBUG(_spdlog, "odofile", "warning: incorrect ODO data record: " + ss.str());
            t_gcoder::_consume(tmpsize);
            _mutex.unlock();
            return -1;
        }
        t_gcoder::_consume(tmpsize);
        cnt++;
    }

    _mutex.unlock();
    return 0;
}
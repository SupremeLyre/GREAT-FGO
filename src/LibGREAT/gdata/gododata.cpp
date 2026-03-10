#include "gododata.h"

great::t_gododata::t_gododata()
{

    id_type(ID_TYPE::ODODATA);
    _last_v = 0.0;
    _last_dv = 0.0;
}
great::t_gododata::t_gododata(t_spdlog spdlog) : t_gdata(spdlog)
{
    id_type(ID_TYPE::ODODATA);
    _last_v = 0.0;
    _last_dv = 0.0;
}

great::t_gododata::~t_gododata()
{
}

int great::t_gododata::add_odo(const double &t, const double &val)
{
    // int it = int(t);
    // for (auto i = 1; i <= 10; i++) {
    // if (double_eq(it + i / 10.0 - t, 0.0)) {
    // if (double_eq(t - it, 0.5) || double_eq(t - it, 0)) {
    _mapodo.insert(make_pair(t, val));
    return 1;
    //    }
    //}
    return -1;
}

bool great::t_gododata::load(const t_gtime &crt, double &val, double dt)
{
    double t = crt.sow() + crt.dsec();

    map<double, double>::iterator it = _mapodo.lower_bound(t);

    if (it != _mapodo.end() && fabs(t - it->first) <= dt)
    {
        _last_dv = it->second - _last_v;
        _last_v = val = it->second;
        return true;
    }
    else
        val = -1e100;
    return false;
}

bool great::t_gododata::avaliable(const t_gtime &crt, double dt)
{
    // std::cout << "[available] this=" << this << " mapodo.addr=" << &_mapodo << "\n";
    double t = crt.sow() + crt.dsec();
    // map<double, double>::iterator it = _mapodo.lower_bound(t);
    // if (it != _mapodo.end() && fabs(t - it->first) <= dt)
    //     return true;
    // else
    //     return false;

    map<double, double>::iterator it = _mapodo.begin();
    while (it != _mapodo.end() && it->first - t <= dt)
    {
        if (t - it->first <= dt)
            return true;
        it = _mapodo.erase(it);
        // it++;
    }
    return false;
}

double great::t_gododata::beg_obs()
{
    if (_mapodo.size() == 0)
        return 0;
    return _mapodo.begin()->first;
}

double great::t_gododata::end_obs()
{
    if (_mapodo.size() == 0)
        return 0;
    return _mapodo.rbegin()->first;
}
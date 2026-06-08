/* ----------------------------------------------------------------------
 * G-Nut - GNSS software development library
 *
  (c) 2018 G-Nut Software s.r.o. (software@gnutsoftware.com)
  This file is part of the G-Nut C++ library.
-*/

#include "gall/gallbias.h"
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <stdlib.h>

using namespace std;

namespace gnut
{
namespace
{
GOBS _make_gobs(const char type, const char band, const char attr)
{
    string obs;
    obs.push_back(type);
    obs.push_back(band);
    obs.push_back(attr);
    return str2gobs(obs);
}

bool _valid_gobs(const GOBS &gobs)
{
    return gobs != X;
}

} // namespace

t_gallbias::t_gallbias() : t_gdata(), _isOverWrite(false)
{
    id_type(t_gdata::ALLBIAS);
    id_group(t_gdata::GRP_MODEL);

    _acOrder["COD_A"] = 1;
    _acOrder["CAS_A"] = 2;
    _acOrder["WHU_A"] = 3;
    _acOrder["DLR_A"] = 4;
    _acOrder["CAS_R"] = 5;
    _acOrder["COD_R"] = 6;
    _acOrder["WHU_R"] = 7;
    _acOrder["DLR_R"] = 8;
    _acOrder["CNT_A"] = 9;
    _acOrder["RTB_A"] = 10;
}

t_gallbias::t_gallbias(t_spdlog spdlog) : t_gdata(spdlog), _isOverWrite(false)
{
    id_type(t_gdata::ALLBIAS);
    id_group(t_gdata::GRP_MODEL);

    _acOrder["COD_A"] = 1;
    _acOrder["CAS_A"] = 2;
    _acOrder["WHU_A"] = 3;
    _acOrder["DLR_A"] = 4;
    _acOrder["CAS_R"] = 5;
    _acOrder["COD_R"] = 6;
    _acOrder["WHU_R"] = 7;
    _acOrder["DLR_R"] = 8;
    _acOrder["CNT_A"] = 9;
    _acOrder["RTB_A"] = 10;
}

t_gallbias::~t_gallbias()
{
    _mapBias.clear();
}

void t_gallbias::add(const string &ac, const t_gtime &epo, const string &obj, t_spt_bias pt_cb)
{

    _gmutex.lock();

    if (pt_cb == nullptr)
    {
        _gmutex.unlock();
        return;
    }

    if (pt_cb->ref() == X)
    {
        _mapBias[ac][epo][obj][pt_cb->gobs()] = pt_cb;
    }
    else
    {
        if (_mapBias[ac][epo][obj].size() == 0)
        {
            shared_ptr<t_gbias> pt_ref = make_shared<t_gbias>(_spdlog); // create first reference bias
            pt_ref->set(pt_cb->beg(), pt_cb->end(), 0.0, pt_cb->ref(),
                        pt_cb->ref());                       // reference bias is set up to zero
            _mapBias[ac][epo][obj][pt_ref->gobs()] = pt_ref; // store new bias (reference)
            _mapBias[ac][epo][obj][pt_cb->gobs()] = pt_cb;   // store new bias
        }
        else
        {
            t_spt_bias pt_obs1 = _find(ac, epo, obj, pt_cb->gobs());
            t_spt_bias pt_obs2 = _find(ac, epo, obj, pt_cb->ref());
            if (pt_obs1 != nullptr && pt_obs2 == nullptr)
            { // connection with first signal
                _connect_first(pt_obs1, pt_cb);
                _mapBias[ac][epo][obj][pt_cb->gobs()] = pt_cb; // store modified bias
            }
            else if (pt_obs1 == nullptr && pt_obs2 != nullptr)
            { // connection with second signal
                _connect_second(pt_obs2, pt_cb);
                _mapBias[ac][epo][obj][pt_cb->gobs()] = pt_cb; // store modified bias
            }
            else if (pt_obs1 != nullptr && pt_obs2 != nullptr)
            {
                // connectin two groups with different reference signal
                // connection with first signal
                _connect_first(pt_obs1, pt_cb);
                // all biases connected with second signal need to be consolidated
                _consolidate(ac, obj, pt_cb, pt_obs2);
            }
            else
            {
                // glfeng add for GAL to store Q & X DCB bias
                shared_ptr<t_gbias> pt_ref = make_shared<t_gbias>(_spdlog); // create first reference bias
                pt_ref->set(pt_cb->beg(), pt_cb->end(), 0.0, pt_cb->ref(),
                            pt_cb->ref());                       // reference bias is set up to zero
                _mapBias[ac][epo][obj][pt_ref->gobs()] = pt_ref; // store new bias (reference)
                _mapBias[ac][epo][obj][pt_cb->gobs()] = pt_cb;   // store new bias
            }
        }
    }

    _gmutex.unlock();
    return;
}

double t_gallbias::get(const string &prd, const t_gtime &epo, const string &prn, const GOBS &gobs, const bool &meter)
{
    _gmutex.lock();

    double bias = 999.0;
    auto itAC = _mapBias.find(prd);
    if (itAC != _mapBias.end())
    {
        auto itEPO = itAC->second.upper_bound(epo);
        if (itEPO != itAC->second.begin() && itEPO != itAC->second.end())
            itEPO--; // between epochs
        if (itEPO == itAC->second.end() && itAC->second.size() != 0)
            itEPO--; // no epochs

        if (itEPO != itAC->second.end())
        {
            auto itSAT = itEPO->second.find(prn);
            if (itSAT != itEPO->second.end() && itSAT->second.find(gobs) != itSAT->second.end())
            {
                t_spt_bias pobs1 = itSAT->second.find(gobs)->second;
                bias = pobs1->bias();
            }
        }
    }

    _gmutex.unlock();
    return bias;
}

bool t_gallbias::get_osb(const t_gtime &epo, const string &obj, const GOBS &gobs, double &bias)
{
    _gmutex.lock();

    const string ac = _select_ac_unlocked();
    const bool found = _get_osb_unlocked(ac, epo, obj, gobs, bias);

    _gmutex.unlock();
    return found;
}

bool t_gallbias::has_osb(const t_gtime &epo, const string &obj, const GOBS &gobs)
{
    double bias = 0.0;
    return get_osb(epo, obj, gobs, bias);
}

vector<string> t_gallbias::get_ac()
{
    vector<string> ac_list;
    for (const auto &item : _mapBias)
    {
        ac_list.push_back(item.first);
    }
    return ac_list;
}

string t_gallbias::get_ac_priority()
{
    string used_ac;
    int loc = 999;
    const map<string, int> ac_order{{"COD_A", 1}, {"CAS_A", 2}, {"WHU_A", 3}, {"DLR_A", 4},  {"CAS_R", 5}, {"COD_R", 6},
                                    {"WHU_R", 7}, {"DLR_R", 8}, {"CNT_A", 9}, {"RTB_A", 10}, {"SGG_A", 11}};
    for (auto item : _mapBias)
    {
        string ac = (item.first == "WHU_A_PHASE") ? "WHU_A" : item.first;
        auto it = ac_order.find(ac);
        int current_loc = (it != ac_order.end()) ? it->second : 999;

        if (current_loc < loc)
        {
            used_ac = ac;
            loc = current_loc;
        }
        else if (used_ac.empty())
        {
            used_ac = ac; // 如果全都没命中预设优先级，至少保证返回一个
        }
    }
    return used_ac;
}

string t_gallbias::get_used_ac()
{
    if (_acUsed.empty())
        _acUsed = get_ac_priority();
    return _acUsed;
}

string t_gallbias::_select_ac_unlocked(const string &tmp)
{
    string ac(tmp);
    if (ac == "" && _isOrdered == true)
        ac = _acPri;

    if (ac == "" && _isOrdered == false)
    {
        int loc = 999;
        for (const auto &item : _mapBias)
        {
            auto it = _acOrder.find(item.first);
            int order = (it != _acOrder.end()) ? it->second : 999;

            if (order < loc)
            {
                ac = item.first;
                loc = order;
            }
            else if (ac.empty())
            {
                ac = item.first;
            }
        }

        _isOrdered = true;
        _acPri = ac;
    }

    return ac;
}

bool t_gallbias::_get_osb_unlocked(const string &ac, const t_gtime &epo, const string &obj, const GOBS &gobs, double &bias)
{
    if (ac.empty() || obj.empty() || gobs == X)
        return false;

    t_gobs normalized(gobs);
    normalized.gobs2to3(t_gsys::char2gsys(obj[0]));
    GOBS obs = normalized.gobs();

    t_spt_bias exact = _find_exact(ac, epo, obj, obs);
    if (exact != nullptr)
    {
        bias = exact->bias();
        return true;
    }

    const string obs_str = gobs2str(obs);
    if (obs_str.size() < 3)
        return false;

    const char sys = obj[0];
    const char type = obs_str[0];
    const char band = obs_str[1];
    const char attr = obs_str[2];

    auto exact_bias = [&](const GOBS &candidate, double &value) -> bool {
        if (!_valid_gobs(candidate) || candidate == obs)
            return false;
        t_spt_bias pt = _find_exact(ac, epo, obj, candidate);
        if (pt == nullptr)
            return false;
        value = pt->bias();
        return true;
    };

    auto average_bias = [&](const GOBS &candidate1, const GOBS &candidate2, double &value) -> bool {
        double bias1 = 0.0;
        double bias2 = 0.0;
        if (!exact_bias(candidate1, bias1) || !exact_bias(candidate2, bias2))
            return false;
        value = 0.5 * (bias1 + bias2);
        return true;
    };

    if (type == 'C')
    {
        // PRIDE read_bias completes vacant GPS/QZSS composite code OSBs by
        // averaging their two component tracking attributes. We keep the
        // product unchanged and synthesize only at lookup time so missing
        // raw product entries are still detectable by callers.
        if ((sys == 'G' || sys == 'J') && attr == 'X')
        {
            if (band == '1' || band == '2')
                return average_bias(_make_gobs('C', band, 'S'), _make_gobs('C', band, 'L'), bias);
            if (band == '5')
                return average_bias(_make_gobs('C', band, 'I'), _make_gobs('C', band, 'Q'), bias);
        }

        // PRIDE treats Galileo E6 C and Q code biases as equivalent to X;
        // if X is absent but C is present, X can also use C. This covers
        // files that publish only the pilot/data/common tracking attribute.
        if (sys == 'E')
        {
            if (band == '6' && attr == 'X' && exact_bias(C6C, bias))
                return true;
            if ((attr == 'C' || attr == 'Q') && exact_bias(_make_gobs('C', band, 'X'), bias))
                return true;
        }

        // Keep the historical GREAT fallback for C6X->C6C. It is applied
        // after the PRIDE rules and only for OSB single-signal lookup.
        if (obs == C6X && exact_bias(C6C, bias))
            return true;
    }
    else if (type == 'L')
    {
        // PRIDE assumes phase OSBs on the same frequency are identical
        // across signal attributes. This is necessary for AR because phase
        // products are often sparser than the RINEX observation attributes.
        const string attrs = "ABCDILMNPQSWXYZ";
        for (char fallback_attr : attrs)
        {
            if (fallback_attr == attr)
                continue;
            if (exact_bias(_make_gobs('L', band, fallback_attr), bias))
                return true;
        }
    }

    return false;
}

double t_gallbias::get(const t_gtime &epo, const string &obj, const GOBS &gobs1, const GOBS &gobs2, const string &tmp)
{

    _gmutex.lock();

    double dcb = 0.0;
    string ac(tmp);
    if (ac == "" && _isOrdered == true)
        ac = _acPri;

    if (ac == "" && _isOrdered == false)
    {
        int loc = 999;
        for (const auto &item : _mapBias)
        {
            auto it = _acOrder.find(item.first);
            int order = (it != _acOrder.end()) ? it->second : 999;

            if (order < loc)
            {
                ac = item.first;
                loc = order;
            }
            else if (ac.empty())
            {
                ac = item.first;
            }
        }

        _isOrdered = true;
        _acPri = ac;
    }

    if (gobs2 == gobs1)
    {
        GOBS gobs1_convert = gobs1;
        this->_convert_obstype(ac, obj, gobs1_convert);
        t_spt_bias pobs1 = _find(ac, epo, obj, gobs1_convert);
        if (pobs1 != nullptr)
        {
            dcb = pobs1->bias();
        }
    }
    else
    {
        GOBS gobs1_convert = gobs1;
        GOBS gobs2_convert = gobs2;

        // align obstype to ac(code or cas)
        this->_convert_obstype(ac, obj, gobs1_convert);
        this->_convert_obstype(ac, obj, gobs2_convert);

        t_spt_bias pobs1 = _find(ac, epo, obj, gobs1_convert);
        t_spt_bias pobs2 = _find(ac, epo, obj, gobs2_convert);

        if (pobs1 != nullptr && pobs2 != nullptr && pobs1->ref() == pobs2->ref())
        {
            dcb = pobs1->bias() - pobs2->bias();
        }
        if (pobs1 != nullptr && pobs2 != nullptr && pobs2->gobs() == C6C)
        {
            dcb = pobs1->bias() - pobs2->bias();
        }
        
    }

    _gmutex.unlock();
    return dcb;
}

t_spt_bias t_gallbias::_find(const string &ac, const t_gtime &epo, const string &obj, const GOBS &gobs)
{
    t_spt_bias pt_bias = nullptr;
    GOBS obs = gobs;
    auto itAC = _mapBias.find(ac);
    if (itAC != _mapBias.end())
    {
        auto itEPO = itAC->second.upper_bound(epo);
        if (itEPO != itAC->second.begin() && itEPO != itAC->second.end())
            itEPO--; // between epochs
        if (itEPO == itAC->second.end() && itAC->second.size() != 0)
            itEPO--; // no epochs

        if (itEPO != itAC->second.end())
        {
            auto itOBJ = itEPO->second.find(obj);
            if (itOBJ != itEPO->second.end())
            {

                auto itGOBS = itOBJ->second.find(obs);
                if (itGOBS == itOBJ->second.end() && obs == C6X)
                {
                    obs = C6C;
                    itGOBS = itOBJ->second.find(obs);
                }
                if (itGOBS != itOBJ->second.end())
                {
                    if (itGOBS->second->valid(epo))
                    {
                        pt_bias = itGOBS->second;
                    }
                }
            }
        }
    }

    return pt_bias;
}

t_spt_bias t_gallbias::_find_exact(const string &ac, const t_gtime &epo, const string &obj, const GOBS &gobs)
{
    t_spt_bias pt_bias = nullptr;

    auto itAC = _mapBias.find(ac);
    if (itAC != _mapBias.end())
    {
        auto itEPO = itAC->second.upper_bound(epo);
        if (itEPO != itAC->second.begin() && itEPO != itAC->second.end())
            itEPO--; // between epochs
        if (itEPO == itAC->second.end() && itAC->second.size() != 0)
            itEPO--; // no epochs

        if (itEPO != itAC->second.end())
        {
            auto itOBJ = itEPO->second.find(obj);
            if (itOBJ != itEPO->second.end())
            {
                auto itGOBS = itOBJ->second.find(gobs);
                if (itGOBS != itOBJ->second.end() && itGOBS->second->valid(epo))
                {
                    pt_bias = itGOBS->second;
                }
            }
        }
    }

    return pt_bias;
}

vector<t_spt_bias> t_gallbias::_find_ref(const string &ac, const t_gtime &epo, const string &obj, const GOBS &ref)
{
    vector<t_spt_bias> vec_bias;

    auto itAC = _mapBias.find(ac);
    if (itAC != _mapBias.end())
    {
        auto itEPO = itAC->second.find(epo);
        if (itEPO != itAC->second.end())
        {
            auto itOBJ = itEPO->second.find(obj);
            if (itOBJ != itEPO->second.end())
            {
                for (auto itGOBS = itOBJ->second.begin(); itGOBS != itOBJ->second.end(); itGOBS++)
                {
                    if (itGOBS->second->ref() == ref)
                        vec_bias.push_back(itGOBS->second);
                }
            }
        }
    }

    return vec_bias;
}

void t_gallbias::_convert_obstype(const string &ac, const string &obj, GOBS &obstype)
{
    if (ac == "COD_R")
    {
        if (obj[0] == 'G' || obj[0] == 'R' || obj[0] == '2' || obj[0] == '3' || obj[0] == '4' || obj[0] == '5')
        {
            switch (obstype)
            {
            case C1C:
                obstype = C1;
                break;
            case C1P:
            case C1Y:
            case C1W:
                obstype = P1;
                break;
            case C2C:
                obstype = C2;
                break;
            case C2P:
            case C2Y:
            case C2W:
                obstype = P2;
                break;
            default:
                break;
            }
        }
    }
    else if (ac == "CAS_R")
    {
        if (obj[0] == 'G' || obj[0] == '2' || obj[0] == '3' || obj[0] == '4' || obj[0] == '5')
        {
            switch (obstype)
            {
            case P1:
                obstype = C1W;
                break;
            case P2:
                obstype = C2W;
                break;
            case C1:
                obstype = C1C;
                break;
            case C2:
                obstype = C2C;
                break;
            default:
                break;
            }
        }
        if (obj[0] == 'R')
        {
            switch (obstype)
            {
            case P1:
                obstype = C1P;
                break;
            case P2:
                obstype = C2P;
                break;
            case C1:
                obstype = C1C;
                break;
            case C2:
                obstype = C2C;
                break;
            default:
                break;
            }
        }
    }
}

void t_gallbias::_connect_first(const t_spt_bias &pt_cb1, const t_spt_bias &pt_cb2)
{
    double newval = pt_cb1->bias() - pt_cb2->bias();
    pt_cb2->set(newval, pt_cb2->ref(), pt_cb1->ref());
}

void t_gallbias::_connect_second(const t_spt_bias &pt_cb1, const t_spt_bias &pt_cb2)
{
    double newval = pt_cb1->bias() + pt_cb2->bias();
    pt_cb2->set(newval, pt_cb2->gobs(), pt_cb1->ref());
}

void t_gallbias::_consolidate(const string &ac, const string &obj, const t_spt_bias &pt_cb1, const t_spt_bias &pt_cb2)
{
    double diff = pt_cb2->val() - pt_cb1->val();
    t_gtime epo = pt_cb1->beg();
    vector<t_spt_bias> vec = _find_ref(ac, epo, obj, pt_cb2->ref());

    for (auto itSPT = vec.begin(); itSPT != vec.end(); itSPT++)
    {
        double newval = (*itSPT)->bias() - diff;
        GOBS gobs = (*itSPT)->gobs();
        GOBS newref = pt_cb1->ref();
        (*itSPT)->set(newval, gobs, newref);
    }
}

} // namespace gnut

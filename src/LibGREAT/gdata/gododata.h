/**
 * @file         gimudata.h
 * @author       GREAT-WHU (https://github.com/GREAT-WHU)
 * @brief        ODO data structure for storing odo data
 * @version      1.0
 * @date         2025-06-01
 *
 * @copyright Copyright (c) 2025, Wuhan University. All rights reserved.
 *
 */

#ifndef GODODATA_H
#define GODODATA_H

#include "gdata/gdata.h"
#include "gexport/ExportLibGREAT.h"
#include "gins/gutility.h"
using namespace gnut;

namespace great
{
/**
 * @brief t_gododata Class for storing the odometry data
 *
 */
class LibGREAT_LIBRARY_EXPORT t_gododata : public t_gdata
{
  public:
    /**
     * @brief default constructor.
     */
    explicit t_gododata();

    /**
     * @brief default destructor.
     */
    explicit t_gododata(t_spdlog spdlog);
    ~t_gododata();

    /**
     * @brief add one data into maps
     * @note data includes time and velocity
     * @param[in] t            time
     * @param[in] val        forward velocity
     * @return int            value reprents success or failure of add
     */
    int add_odo(const double &t, const double &val);

    /**
     * @brief load data into pos and v and t accroding to direction
     * @note direction stands for loading direction
     * @param[in] t   current epoch(int sec)
     * @param[in] dt  fractional second
     * @param[out] val   data value
     * @return bool
     */
    bool load(const t_gtime &t, double &val, double dt);

    /**
    * @brief judge whether IMU data is available
    * @param[in] t   current epoch(int sec)
    * @param[in] dt  fractional second
    * @return
        @retval true   available
        @retval false  unavailable
    */
    bool avaliable(const t_gtime &t, double dt);

    /**
     * @brief get the begin time of imu data
     * @param[in] '_beg_end' is the direction
     * @return double  first epoch
     */
    virtual double beg_obs();

    /**
     * @brief get the end time of imu data
     * @param[in] '_beg_end' is the direction
     * @return double  end epoch
     */
    virtual double end_obs();

  private:
    map<double, double> _mapodo;
    // to judge whether the tire is sliding or the road is bumpy
    double _last_v;
    double _last_dv;
};

} // namespace great

#endif
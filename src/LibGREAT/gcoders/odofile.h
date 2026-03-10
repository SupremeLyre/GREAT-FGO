#ifndef ODOFILE_H
#define ODOFILE_H

/**
* @file odofile.h
* @brief read and decode odometry file

* @verbatim
    History
    -1.0    ShenZhiheng          2021-05-02    creat the file.
    -1.1    WuZongzhou          2021-07-26    Modify Doxygen Style Code Remarks
  @endverbatim

* @author ShenZhiheng
* @version 1.0
* @date 2021-05-02
* @license
*/

#include "gcoders/gcoder.h"
#include "gins/gutility.h"
#include "gutils/gstring.h"
#include <Eigen/Eigen>
using namespace std;
using namespace gnut;
using namespace great;

namespace great
{

/**
 * @class t_odofile
 * @brief t_odofile Class for decoding and encoding the ODO file
 *
 */
class LibGREAT_LIBRARY_EXPORT t_odofile : public t_gcoder
{
  public:
    /**
     * @brief constructor.
     *
     * @param[in]  s        setbase control
     * @param[in]  version  version of the gcoder
     * @param[in]  sz       size of the buffer
     */
    explicit t_odofile(t_gsetbase *s, string version = "", int sz = DEFAULT_BUFFER_SIZE);

    /**
     * @brief destructor.
     */
    ~t_odofile()
    {
    }

    /**
    * @brief decode the header of the ODO data file.
    * odo file doesn't include head block
    *
    * @param[in]  buff        buffer of the data
    * @param[in]  sz          buffer size of the data
    * @param[in]  errmsg      error message of the data decoding
    * @return
        @retval >=0 consume size of header decoding
        @retval <0  finish reading
    */
    virtual int decode_head(char *buff, int sz, vector<string> &errmsg);

    /**
    * @brief decode the data body of the ODO data file.
    *
    *
    * @param[in]  buff        buffer of the data
    * @param[in]  sz          buffer size of the data
    * @param[in]  cnt         size of lines successfully decoded
    * @param[in]  errmsg      error message of the data decoding
    * @return
        @retval >=0 consume size of header decoding
        @retval <0  finish reading
    */
    virtual int decode_data(char *buff, int sz, int &cnt, vector<string> &errmsg);

  private:
    double _tt;
    double _info; /// info=displacement/pulse_count
    double _ts;
    string _format;
    int _sign;
};
} // namespace great

#endif
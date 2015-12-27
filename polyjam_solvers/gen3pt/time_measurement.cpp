#include "time_measurement.hpp"

#ifdef WIN32

#include <sys\timeb.h>

void
gettimeofday( timeval * timeofday, int dummy)
{
  struct timeb time;
  ftime(&time);
  timeofday->tv_sec = (int) time.time;
  timeofday->tv_usec = 1000 * (int) time.millitm;
}

#endif

timeval
timeval_minus( const struct timeval &t1, const struct timeval &t2 )
{
  timeval ret;
  ret.tv_sec = t1.tv_sec - t2.tv_sec;
  if( t1.tv_usec < t2.tv_usec )
  {
    ret.tv_sec--;
    ret.tv_usec = t1.tv_usec - t2.tv_usec + 1000000;
  }
  else
    ret.tv_usec = t1.tv_usec - t2.tv_usec;

  return ret;
}


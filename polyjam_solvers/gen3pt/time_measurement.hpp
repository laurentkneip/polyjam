#ifndef OPENGV_TIME_MEASUREMENT_HPP_
#define OPENGV_TIME_MEASUREMENT_HPP_

#include <stdlib.h>
#include <stdio.h>

#ifdef WIN32
  struct timeval {
    int tv_sec;
    int tv_usec;
  };

  void gettimeofday( timeval * timeofday, int dummy);
#else
  #include <sys/time.h>
#endif

#define TIMETODOUBLE(x) ( x.tv_sec + x.tv_usec * 1e-6 )

timeval
timeval_minus( const struct timeval &t1, const struct timeval &t2 );

#endif /* OPENGV_TIME_MEASUREMENT_HPP_ */

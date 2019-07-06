#ifndef __BASEVARC_STATS_H__
#define __BASEVARC_STATS_H__

#include "BamProcess.h"
#include "algorithm.h"
#include "htslib/kfunc.h"

static const char STRAND[2] = {'-', '+'};
static const char BASE2CHAR[4] = {'A', 'C', 'G', 'T'};
static const char BASE_INT8_TABLE[128] = {
    -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1,
    -1,  0, -1,  1, -1, -1, -1,  2,
    -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1,  3, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1,
};

struct Stat
{
    double phred_qual;
    double phred_mapq;
    double phred_rpr;
    double fs;
    double sor;
    int ref_fwd;
    int ref_rev;
    int alt_fwd;
    int alt_rev;
};


inline double chisf(double x, double k)
{
    // sf = 1 - cdf
    return kf_gammaq(k/2.0, x/2.0);
}

inline double normsf(double x)
{
    // mean = 0, sd = 1
    return 1 - (1 + erf(x/sqrt(2.0)))/2;
}


#endif

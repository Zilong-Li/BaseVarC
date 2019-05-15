#ifndef __BASEVARC_STATS_H__
#define __BASEVARC_STATS_H__

#include "BamProcess.h"
#include "ranksum.h"
#include "em.h"
#include "htslib/kfunc.h"
#include "htslib/bgzf.h"

typedef std::vector<double> ProbV;
typedef std::vector<ProbV> FreqV;
typedef std::vector<int8_t> BaseV;
typedef std::vector<BaseV> CombV;
typedef std::unordered_map<int32_t, int32_t> DepM;
typedef std::string String;

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


inline double chisf(const double& x, const double& k)
{
    // sf = 1 - cdf
    return 1 - kf_gammap(k/2, x/2);
}

inline double normsf(const double& x)
{
    // mean = 0, sd = 1
    return 1 - (1 + erf(x/sqrt(2)))/2;
}


#endif

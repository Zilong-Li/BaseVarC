#ifndef __BASEVARC_BASE_TYPE_H__
#define __BASEVARC_BASE_TYPE_H__

#include <algorithm>
#include <complex>
#include <unordered_map>
#include <vector>
#include "em.h"
#include "ranksum.h"
#include "htslib/kfunc.h"

#define LRT_THRESHOLD 24    // chi-pvalue of 10^-6
#define QUAL_THRESHOLD 60   // -10 * lg(10^-6)
#define MLN10TO10 -0.23025850929940458    // -log(10)/10
#define MINAF 0.001                       // base freqence threshold
#define NTYPE 4

const static uint8_t BASE[4] = {0, 1, 2, 3};

typedef std::vector<float> ProbV;
typedef std::vector<ProbV> FreqV;
typedef std::vector<uint8_t> BaseV;
typedef std::vector<BaseV> CombV;
typedef std::unordered_map<int, int> DepM;

void combs(const BaseV& bases, CombV& comb_v, int32_t k)
{
    // k <= n
    int32_t n = bases.size();
    std::string bitmask(k, 1); // K leading 1's
    bitmask.resize(n, 0); // N-K trailing 0's

    BaseV bv;
    comb_v.clear();
    bv.reserve(k);
    do {
        for (int32_t i = 0; i < n; ++i) // [0..N-1] integers
        {
            if (bitmask[i]) bv.push_back(bases[i]);
        }
        comb_v.push_back(bv);
        bv.clear();
    } while (std::prev_permutation(bitmask.begin(), bitmask.end()));
}

inline double chisf(double x, double k)
{
    // sf = 1 - cdf
    return 1 - kf_gammap(k/2, x/2);
}

inline double normsf(double x)
{
    // mean = 0, sd = 1
    return 1 - (1 + erf(x/sqrt(2)))/2;
}

class BaseType
{
 public:
    // REMINDME: consider using AlleleInfoVector directly
    BaseType(BaseV base, BaseV qual, int ref, double minaf);
    ~BaseType() {
        delete []ind_allele_likelihood;
    }
    void LRT();
    double var_qual;
    BaseV alt_bases;
    ProbV af_lrt;
    // std::unordered_map<uint8_t, float> af_lrt;

 private:
    BaseV bases;
    BaseV quals;
    const int ref_base;
    const double min_af;
    double depth_total = 0;
    int32_t nind;
    DepM depth{ {0, 0},{1, 0},{2, 0},{3, 0}};
    double *ind_allele_likelihood;
    double init_allele_freq[NTYPE];

    void SetAlleleFreq(const BaseV& bases);
    void UpdateF(const BaseV& bases, CombV& bc, ProbV& lr, FreqV& bp, int32_t k);
};


#endif
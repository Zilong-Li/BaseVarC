#ifndef __BASEVARC_BASE_TYPE_H__
#define __BASEVARC_BASE_TYPE_H__

#include <algorithm>
#include <complex>
#include <sstream>
#include "Stats.h"

#define LRT_THRESHOLD 24    // chi-pvalue of 10^-6
#define MLN10TO10 -0.23025850929940458    // -log(10)/10
#define MINAF 0.001                       // base freqence threshold
#define QUAL_THRESHOLD 60   // -10 * lg(10^-6)
#define NTYPE 4

const static int8_t BASE[4] = {0, 1, 2, 3};

class BaseType
{
 public:
    // REMINDME: consider using AlleleInfoVector directly
    BaseType(BaseV base, BaseV qual, int ref, double minaf);
    ~BaseType() {
        delete []ind_allele_likelihood;
        delete []init_allele_freq;
    }
    void LRT();
    void stats(int8_t ref_base, const BaseV& alt_bases, const AlleleInfoVector& aiv, Stat& s);
    void writeVcf(const String& chr, int32_t pos, int8_t ref_base, const BaseType& bt, const AlleleInfoVector& aiv, const DepM& idx, int32_t N);
    double var_qual;
    BaseV alt_bases;
    ProbV af_lrt;
    DepM depth{ {0, 0},{1, 0},{2, 0},{3, 0} };
    double depth_total = 0;

 private:
    BaseV bases;
    BaseV quals;
    const int ref_base;
    const double min_af;
    int32_t nind;
    double *ind_allele_likelihood;
    double *init_allele_freq;

    void combs(const BaseV& bases, CombV& comb_v, int32_t k);
    void SetAlleleFreq(const BaseV& bases);
    void UpdateF(const BaseV& bases, CombV& bc, ProbV& lr, FreqV& bp, int32_t k);
};


#endif
#ifndef __BASEVARC_BASE_TYPE_H__
#define __BASEVARC_BASE_TYPE_H__

#include <algorithm>
#include "Algorithm.h"
#include "BamProcess.h"
#include "robin_hood.h"

#define LRT_THRESHOLD 24.0    // chi-pvalue of 10^-6
#define MLN10TO10 -0.23025850929940458    // -log(10)/10
#define MINAF 0.001                       // base freqence threshold
#define QUAL_THRESHOLD 60   // -10 * lg(10^-6)
#define NTYPE 4

typedef std::string String;
typedef std::vector<double> ProbV;
typedef std::vector<ProbV> FreqV;
typedef std::vector<int8_t> BaseV;
typedef std::vector<BaseV> CombV;
typedef std::map<String, String> InfoM;
typedef robin_hood::unordered_map<int32_t, int32_t> DepM;
static const int BASE[4] = {0, 1, 2, 3};
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
    int ref_fwd = 0;
    int ref_rev = 0;
    int alt_fwd = 0;
    int alt_rev = 0;
};

void combs_(const BaseV& bases, CombV& comb_v, int32_t k);

class BaseType
{
    friend String WriteVcf(const BaseType& bt, const String& chr, int32_t pos, int8_t ref_base, const AlleleInfoVector& aiv, const DepM& idx, InfoM& info, int32_t N);

 public:
    BaseType(BaseV base, BaseV qual, int8_t ref, double minaf);
    ~BaseType() {}

    void SetBase (const BaseV& v) { this->base_comb = v; }
    bool LRT();

    double var_qual;
    double depth_total;
    BaseV alt_bases;
    DepM depth;
    robin_hood::unordered_map<int8_t, double> af_lrt;

 private:
    BaseV bases;
    BaseV quals;
    BaseV base_comb{0, 1, 2, 3};
    const int8_t ref_base;
    const double min_af;
    const int32_t nind;
    ProbV ind_allele_likelihood;
    ProbV init_allele_freq;

    void SetAlleleFreq(const BaseV& bases);

    void UpdateF(const BaseV& bases, CombV& bc, ProbV& lr, FreqV& bp, int32_t k);
};


#endif

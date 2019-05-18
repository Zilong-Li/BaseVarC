#ifndef __BASEVARC_BASE_TYPE_H__
#define __BASEVARC_BASE_TYPE_H__

#include <algorithm>
#include <complex>
#include <sstream>
#include "Stats.h"
#include "htslib/bgzf.h"

#define LRT_THRESHOLD 24    // chi-pvalue of 10^-6
#define MLN10TO10 -0.23025850929940458    // -log(10)/10
#define MINAF 0.001                       // base freqence threshold
#define QUAL_THRESHOLD 60   // -10 * lg(10^-6)
#define NTYPE 4

typedef std::string String;
typedef std::vector<double> ProbV;
typedef std::vector<ProbV> FreqV;
typedef std::vector<int8_t> BaseV;
typedef std::vector<BaseV> CombV;
typedef std::unordered_map<int32_t, int32_t> DepM;
typedef std::unordered_map<int32_t, String> SamM;

static const int8_t BASE[4] = {0, 1, 2, 3};

static const char* CVG_HEADER =
"##fileformat=CVGv1.0\n"
"##Group information is the depth of A:C:G:T\n"
"#CHROM\tPOS\tREF\tDepth\tA\tC\tG\tT\tFS\tSOR\tStrand_Coverage(REF_FWD,REF_REV,ALT_FWD,ALT_REV)\n";

static const char* VCF_HEADER =
"##fileformat=VCFv4.2\n"
"##FILTER=<ID=LowQual,Description=\"Low quality (QUAL < 60)\">\n"
"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
"##FORMAT=<ID=AB,Number=1,Type=String,Description=\"Allele Base\">\n"
"##FORMAT=<ID=SO,Number=1,Type=String,Description=\"Strand orientation of the mapping base. Marked as + or -\">\n"
"##FORMAT=<ID=BP,Number=1,Type=String,Description=\"Base Probability which calculate by base quality\">\n"
"##INFO=<ID=CM_AF,Number=A,Type=Float,Description=\"An ordered, comma delimited list of allele frequencies base on LRT algorithm\">\n"
"##INFO=<ID=CM_CAF,Number=A,Type=Float,Description=\"An ordered, comma delimited list of allele frequencies just base on read count\">\n"
"##INFO=<ID=CM_AC,Number=A,Type=Integer,Description=\"An ordered, comma delimited allele depth in CMDB\">\n"
"##INFO=<ID=CM_DP,Number=A,Type=Integer,Description=\"Total Depth\">\n"
"##INFO=<ID=SB_REF,Number=A,Type=Integer,Description=\"Read number support REF: Forward,Reverse\">\n"
"##INFO=<ID=SB_ALT,Number=A,Type=Integer,Description=\"Read number support ALT: Forward,Reverse\">\n"
"##INFO=<ID=FS,Number=1,Type=Float,Description=\"Phred-scaled p-value using Fisher's exact test to detect strand bias\">\n"
"##INFO=<ID=BaseQRankSum,Number=1,Type=Float,Description=\"Phred-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities\">\n"
"##INFO=<ID=SOR,Number=1,Type=Float,Description=\"Symmetric Odds Ratio of 2x2 contingency table to detect strand bias\">\n"
"##INFO=<ID=MQRankSum,Number=1,Type=Float,Description=\"Phred-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities\">\n"
"##INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description=\"Phred-score from Wilcoxon rank sum test of Alt vs. Ref read position bias\">\n"
"##INFO=<ID=QD,Number=1,Type=Float,Description=\"Variant Confidence Quality by Depth\">\n";

class BaseType
{
 public:
    // REMINDME: consider using AlleleInfoVector directly
    BaseType(BaseV base, BaseV qual, int8_t ref, double minaf);
    ~BaseType() {
        delete []ind_allele_likelihood;
        delete []init_allele_freq;
    }
    bool LRT();
    void WriteVcf(BGZF* fpv, const BaseType& bt, const String& chr, int32_t pos, int8_t ref_base, const AlleleInfoVector& aiv, const DepM& idx, int32_t N);

    double var_qual;
    BaseV alt_bases;
    DepM depth;
    double depth_total;
    std::unordered_map<int8_t, double> af_lrt;

 private:
    BaseV bases;
    BaseV quals;
    const int8_t ref_base;
    const double min_af;
    const int32_t nind;
    double *ind_allele_likelihood;
    double *init_allele_freq;

    void stats(int8_t ref_base, const BaseV& alt_bases, const AlleleInfoVector& aiv, Stat& s);
    void combs(const BaseV& bases, CombV& comb_v, int32_t k);
    void SetAlleleFreq(const BaseV& bases);
    void UpdateF(const BaseV& bases, CombV& bc, ProbV& lr, FreqV& bp, int32_t k);
};


#endif
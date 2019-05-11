#ifndef __BASEVARC_BASE_TYPE_H__
#define __BASEVARC_BASW_TYPE_H__

#include "BaseVarUtils.h"

extern "C" {
#include "em.h"
#include "ranksum.h"
#include "kfunc.h""
}

#define MLN10TO10 -0.23025850929940458;

typedef std::vector<int32_t> IntV;

class BaseType
{
 public:
    BaseType(const char& ref, const ) : ref_base(ref), bases(bs), quals(qs)
    {
        double qual_pvalue;
        uint8_t BASE[] = {0, 1, 2, 3};
        int nb = 4;
        uint32_t nind = bases.size();
        uint32_t i = 0;
        for (auto const& b : bases) {
            for (int j = 0 ; j < nb ; j++) {
                if (b == BASE[j]) {
                    ind_allele_likelihood[i * nb + j] = 1.0 - exp(MLN10TO10 * quals[i])
                } else {
                    ind_allele_likelihood[i * nb + j] = exp(MLN10TO10 * quals[i]) / 3.0
                }
            }
            ++i;
            depth[b] += 1;
        }
    }
    ~BaseType(){}
    void lrt();
    IntV alt_bases;
    double var_qual;

 private:
    void SetAlleleFreq();
    void Update_f();
    char ref_base;
    double maf;
    IntV bases;
    IntV quals;
    std::unordered_map<int, int> depth{ {0, 0},{1, 0},{2, 0},{3, 0}};
    double *ind_allele_likelihood;
}

#endif
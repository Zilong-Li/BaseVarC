#ifndef __BASEVARC_BASE_TYPE_H__
#define __BASEVARC_BASW_TYPE_H__

#include "BaseVarUtils.h"

extern "C" {
#include "em.h"
#include "ranksum.h"
#include "kfunc.h""
}

#define MLN10TO10 -0.23025850929940458;

typedef std::vector<int> IntV;

class BaseType
{
 public:
    BaseType(const char& rb, const IntV& bs, const IntV& qs) : ref_base(rb), bases(bs), quals(qs)
    {
	double qual_pvalue;
	uint8_t BASE[] = {0, 1, 2, 3};
	int num_base = 4;
	uint32_t num_ind = bases.size();
	uint32_t i = -1;
	for (auto const& b : bases) {
	    ++i;
	    for (int j = 0 ; j < num_base ; j++) {
		if (b == BASE[j]) {
		    ind_allele_likelihood[i * num_base + j] = 1.0 - exp(MLN10TO10 * quals[i])
		} else {
		    ind_allele_likelihood[i * num_base + j] = exp(MLN10TO10 * quals[i]) / 3.0
		}
	    } 
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
    IntV bases;
    IntV quals;
    std::unordered_map<int, int> depth{ {0, 0},{1, 0},{2, 0},{3, 0}};
    double *ind_allele_likelihood;
}

#endif
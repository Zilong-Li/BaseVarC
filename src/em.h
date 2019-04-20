#ifndef BASEVAR_EM_H
#define BASEVAR_EM_H

#ifdef __cplusplus
extern "C" {
#endif

void EM(double *init_allele_freq, double *ind_allele_likelihood, double *marginal_likelihood, double *expect_allele_prob, int nsample, int ntype, int iter_num, double epsilon);

#ifdef __cplusplus
}
#endif

#endif

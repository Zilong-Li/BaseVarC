#ifndef BASEVAR_ALGORITHM_H
#define BASEVAR_ALGORITHM_H

#ifdef __cplusplus
extern "C" {
#endif

double RankSumTest(double *x, int n1, double *y, int n2);

void EM(double *init_allele_freq, double *ind_allele_likelihood, double *marginal_likelihood, double *expect_allele_prob, int nsample, int ntype, int iter_num, double epsilon);

#ifdef __cplusplus
}
#endif

#endif

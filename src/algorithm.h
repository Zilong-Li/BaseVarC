#ifndef BASEVARC_ALGORITHM_H
#define BASEVARC_ALGORITHM_H

#include <cmath>
#include "htslib/kfunc.h"

double chisf(double x, double k);

double normsf(double x);

double bt_fisher_exact(int n11, int n12, int n21, int n22);


double RankSumTest(std::vector<double>& x, std::vector<double>& y);

void EM(std::vector<double>& init_allele_freq, const std::vector<double>& ind_allele_likelihood, std::vector<double>& marginal_likelihood, std::vector<double>& expect_allele_prob, int32_t nsample, int ntype, int iter_num, double epsilon);

#endif

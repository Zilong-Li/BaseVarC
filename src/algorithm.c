/*
*   lizilong@bgi.com 201903
*/
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "algorithm.h"

static int compare_floats(const void* a, const void* b)
{
    double arg1 = *(const double*)a;
    double arg2 = *(const double*)b;
    if (arg1 < arg2) return -1;
    if (arg1 > arg2) return 1;
    return 0;
}

static double rankR1(double *x, int n1, double *y, int n2)
{
    // merge and store the index of sample1
    int ia = n1 - 1;
    int ib = n2 - 1;
    int i, j;
    int *index = malloc(n1 * sizeof(int));
    for (i = n1 + n2 - 1; i>=0; i--) {
        if (ia >= 0 && ib < 0) {
            while(ia >= 0){
                index[ia] = ia;
                ia--;
            }
            break;
        }
        if (ia < 0 && ib >= 0) {
            x[i] = y[ib--];
        }
        if (ia >= 0 && ib >= 0) {
            if (x[ia] > y[ib]) {
                index[ia] = i;
                x[i] = x[ia--];
            }else{
                x[i] = y[ib--];
            }
        }
    }
    // average method
    int k = 0, n = 0;
    for (i = 0; i < n1 + n2; i++){
        if (x[i] == x[i+1]) {
            k += i + 1;
            n++;
        }else{
            if (k > 0) {
                k += i + 1;
                double avg = (double)k / (n + 1);
                for (j = i; j >= i - n; j--) {
                    x[j] = avg;
                }
                k = 0;
                n = 0;
            }else{
                x[i] = i + 1;
            }
        }
    }

    //return R1 - the sum of the ranks in Sample1
    double r1 = 0;
    for (i = 0; i < n1; i++) {
        r1 += x[index[i]];
    }

    free(index);
    return r1;
}

double RankSumTest(double *x, int n1, double *y, int n2)
{
    int n = n1 + n2;
    double *xx = malloc(n * sizeof(double));
    qsort(x, n1, sizeof(double), compare_floats);
    qsort(y, n2, sizeof(double), compare_floats);
    memcpy(xx, x, n1 * sizeof(double));
    double r1 = rankR1(xx, n1, y, n2);
    double expected = (double)(n1 * (n1 + n2 + 1)) / 2.0;
    double z = (r1 - expected) / sqrt((double)(n1*n2*(n1+n2+1))/12.0);

    free(xx);
    return z;
}

static void singleEM(double *allele_freq, double *ind_allele_likelihood, double *marginal_likelihood, double *expect_allele_prob, int nsample, int ntype)
{
    int i, j;
    double *likelihood = calloc(ntype, sizeof(double));
    double *ind_allele_prob = calloc(ntype * nsample, sizeof(double)); // covert ind_allele_prob to col major 
    // step E 
    for(i=0; i<nsample; ++i){
        for(j=0; j<ntype; ++j){
            likelihood[j] = allele_freq[j] * ind_allele_likelihood[i * ntype + j];
            marginal_likelihood[i] += likelihood[j];
        }
        for(j=0; j<ntype; ++j){
            /* col major may be fast for step M
             * need to deal with marginal_likelihood[i] is close to zero */
            ind_allele_prob[j * nsample + i] = likelihood[j] / marginal_likelihood[i]; 
        }
    }
    free(likelihood);

    // step M
    for(j=0; j<ntype; ++j){
        for(i=0; i<nsample; ++i){
            expect_allele_prob[j] += ind_allele_prob[j * nsample + i]; 
        }
        expect_allele_prob[j] = expect_allele_prob[j] / nsample;
    }
    free(ind_allele_prob);
}

static void update_allele_freq(double *allele_freq, double *expect_allele_prob, int ntype)
{
    for(int j=0; j<ntype; ++j){
        allele_freq[j] = expect_allele_prob[j];
        expect_allele_prob[j] = 0.0;
    }
}

static double delta_bylog(double *bf, double *af, int n)
{
    double delta = 0.0;
    for(int i=0; i<n; ++i){
        // need to deal with log(0) == inf;
        delta += fabs(log(af[i]) - log(bf[i]));
        bf[i] = af[i];
        af[i] = 0.0;
    }
    return delta;
}

void EM(double *init_allele_freq, double *ind_allele_likelihood, double *marginal_likelihood, double *expect_allele_prob, int nsample, int ntype, int iter_num, double epsilon)
{
    double *af_marginal_likelihood = calloc(nsample, sizeof(double));
    double *allele_freq = malloc(ntype * sizeof(double));
    double delta;
    int i, j;

    /*
     copy allele_freq in case that init_allele_freq be modified; 
    */
    for(j = 0; j < ntype; ++j){
        allele_freq[j] = init_allele_freq[j];
    }
    singleEM(allele_freq, ind_allele_likelihood, marginal_likelihood, expect_allele_prob, nsample, ntype);

    for(i=0; i<iter_num; ++i){
        update_allele_freq(allele_freq, expect_allele_prob, ntype);
        singleEM(allele_freq, ind_allele_likelihood, af_marginal_likelihood, expect_allele_prob, nsample, ntype);
        delta = delta_bylog(marginal_likelihood, af_marginal_likelihood, nsample);
        if(delta < epsilon){
            break;
        }
    }

    free(af_marginal_likelihood);
    free(allele_freq);
}


/*
**      Author:     Tapas Kanungo, kanungo@cfar.umd.edu
**      Date:       15 December 1997
**      File:       backward.c
**      Purpose:    Backward algorithm for computing the probabilty
**                  of observing a sequence given a HMM model parameter.
**      Organization: University of Maryland
**
**      $Id: backward.c,v 1.3 1998/02/23 07:56:05 kanungo Exp kanungo $
*/

/*  Solution to the evaluation problem, ie:
    Given the sequence of observations O=O1,O2,...,OT and the hidden Markov model
    lambda = (A, B, pi), how to can be calculated efficiently P(O|lambda),
    the probability of the sequence of observations, given the lambda model. */

#include <stdio.h>
#include "hmm.h"
static char rcsid[] = "$Id: backward.c,v 1.3 1998/02/23 07:56:05 kanungo Exp kanungo $";

/*  Where beta is the auxiliary variable backward, defined as the probability of
    partial sequence of observations from t+1 to the end, given the state Si at time t,
    and the model lambda=(A, B, pi) of the hidden Markov model.
    beta_t(i)=P(Ot+1 Ot+2...OT, q_t=S_i|lambda) */

void Backward(HMM *phmm, int T, int *O, double **beta, double *pprob)
{
    /*  Indexes of the state */
    int i, j;
    /*  Indexes of the time */
    int t;
    /*  Partial sum */
    double sum;

    /*  1. Initialization */
    /*  beta_T(i) = 1, with 1 <= i <= N */
    for (i = 1; i <= phmm->N; i++)
        beta[T][i] = 1.0;

    /*  2. Induction */
    /*  for t=T-1,T-2...1, with 1 <= i <= N
        beta_t+1(j) = {sum a_ij*b_j(Ot+1)}*beta_t+1(j) */
    for (t = T - 1; t >= 1; t--)
    {
        for (i = 1; i <= phmm->N; i++)
        {
            sum = 0.0;
            for (j = 1; j <= phmm->N; j++)
                sum += phmm->A[i][j] * (phmm->B[j][O[t+1]])*beta[t+1][j];
            beta[t][i] = sum;
        }
    }

    /* 3. Termination
         This phase in the Backward there is not. */
    /*
    // It may be like:
    *pprob = 0.0;
    for (i = 1; i <= phmm->N; i++)
        *pprob += beta[1][i];
    */
}

void BackwardWithScale(HMM *phmm, int T, int *O, double **beta,	double *scale, double *pprob)
{
    int i, j;
    int t;
	double sum;

    /* 1. Initialization */
    for (i = 1; i <= phmm->N; i++)
        beta[T][i] = 1.0/scale[T];

    /* 2. Induction */
    for (t = T - 1; t >= 1; t--)
    {
        for (i = 1; i <= phmm->N; i++)
        {
			sum = 0.0;
            for (j = 1; j <= phmm->N; j++)
                sum += phmm->A[i][j] * (phmm->B[j][O[t+1]])*beta[t+1][j];
            beta[t][i] = sum/scale[t];

        }
    }

    /* 3. Termination
        This phase in the Backward there is not.  */
}

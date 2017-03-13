/*
**      Author: Tapas Kanungo, kanungo@cfar.umd.edu
**      Date:   15 December 1997
**      File:   viterbi.c
**      Purpose: Viterbi algorithm for computing the maximum likelihood
**		state sequence and probablity of observing a sequence
**		given the model.
**      Organization: University of Maryland
**
**      $Id: viterbi.c,v 1.1 1999/05/06 05:25:37 kanungo Exp kanungo $
*/

/*
        Change
        Author:     Pietro Russo
        Date:       12/03/2011
        File:       viterbi.c
        Purpose:    Added the ViterbiLog_C function that adapts the Viterbi algorithm
                    to the speech recognition problem.
*/

/*  Solution to the problem of decoding:
    given the sequence of observations O = O1O2 ... OT and the Hidden Markovian Model lambda = (A, B, pi),
    how to choose a corresponding state sequence Q = q1 q2 ... qT that is optimal,
    to justify the better observations? */
#include <math.h>
#include "hmm.h"
#include "nrutil.h"
static char rcsid[] = "$Id: viterbi.c,v 1.1 1999/05/06 05:25:37 kanungo Exp kanungo $";

#define VITHUGE  100000000000.0
#define VITHUGE2  100000000.0
//#define VITHUGE  -100000000000.0
/*  Want to find Q (sequences of states) that maximizes P(O,Q|lambda).
    It defines the quantity:
    delta_t(i) = max P(q1 q2 ... qt = i, O2 O1...Ot|lambda)
    which is the highest probability of a single path, at time t, which justifies
    the first t observations and ends in the state Si.
    By induction we have:
    delta_t+1(j) = {max delta_t(s)*a_ij}+b_j(Ot+1)
    It keeps track of the argument that maximizes the above formula for each t e j
    in psi_t (j)
*/
/*  Function of the Viterbi algorithm */
void Viterbi(HMM *phmm, int T, int *O, double **delta, int **psi, int *q, double *pprob)
{
    /*  Indes of the states */
	int i, j;
	/*  Index of the time */
	int  t;

	int	maxvalind;
	double	maxval, val;

	/*  1. Initialisation */
	/*  delta_1(i) = pi_i*b_i(O1), with 1 <= i <= N
        psi_1(i)=0 */
	for (i = 1; i <= phmm->N; i++)
	{
		delta[1][i] = phmm->pi[i] * (phmm->B[i][O[1]]);
		psi[1][i] = 0;
	}

	/*  2. Recursion */
	/*  delta_t(j) = max{delta_t-1(i)*a_ij}*b_j(Ot), with 2 <= t <= T e 1 <= j <= N
        psi_t(j) = arg max{delta_t-1(i)*a_ij}, with 2 <= t <= T e 1 <= j <= N */
	for (t = 2; t <= T; t++)
	{
		for (j = 1; j <= phmm->N; j++)
		{
			maxval = 0.0;
			maxvalind = 1;
			for (i = 1; i <= phmm->N; i++)
			{
				val = delta[t-1][i]*(phmm->A[i][j]);
				//printf("Val: %g\t", val);
				if (val > maxval)
				{
				    maxval = val;
					maxvalind = i;
				}
			}
			delta[t][j] = maxval*(phmm->B[j][O[t]]);
			psi[t][j] = maxvalind;
		}
	}

	/*  3. Termination */
	/*  P* = max delta_t(i)
        q*_T = arg max delta_T(i) */
	*pprob = 0.0;
	q[T] = 1;
	for (i = 1; i <= phmm->N; i++)
	{
	    if (delta[T][i] > *pprob)
	    {
			*pprob = delta[T][i];
			q[T] = i;
		}
	}

    /*  4. Backtracking on the path (sequence of states) */
	/*  q*_t = psi_t+1 (q*_t+1), for t=T-1, T-2, ..., 1*/
	for (t = T - 1; t >= 1; t--)
		q[t] = psi[t+1][q[t+1]];
}

/*  Function with the Viterbi algorithm, with operations with logarithms. */
void ViterbiLog(HMM *phmm, int T, int *O, double **delta, int **psi, int *q, double *pprob)
{
    int i, j;
    int t;

    int     maxvalind;
    double  maxval, val;

	double  **biot;

	/*  0. Preprocessing */

	for (i = 1; i <= phmm->N; i++)
		phmm->pi[i] = log(phmm->pi[i]);
	for (i = 1; i <= phmm->N; i++)
		for (j = 1; j <= phmm->N; j++)
		{
			phmm->A[i][j] = log(phmm->A[i][j]);
		}

	biot = dmatrix(1, phmm->N, 1, T);
	for (i = 1; i <= phmm->N; i++)
		for (t = 1; t <= T; t++)
		{
			biot[i][t] = log(phmm->B[i][O[t]]);
		}

    /*  1. Initialisation  */

    for (i = 1; i <= phmm->N; i++)
    {
        delta[1][i] = phmm->pi[i] + biot[i][1];
        psi[1][i] = 0;
    }

    /*  2. Recursion */
    for (t = 2; t <= T; t++)
    {
        for (j = 1; j <= phmm->N; j++)
        {
            maxval = -VITHUGE;
            maxvalind = 1;
            for (i = 1; i <= phmm->N; i++)
            {
                val = delta[t-1][i] + (phmm->A[i][j]);
                if (val > maxval)
                {
                    maxval = val;
                    maxvalind = i;
                }
            }
            delta[t][j] = maxval + biot[j][t];
            psi[t][j] = maxvalind;
        }
    }

    /*  3. Termination */
    *pprob = -VITHUGE;
    q[T] = 1;
    for (i = 1; i <= phmm->N; i++)
    {
        if (delta[T][i] > *pprob)
        {
            *pprob = delta[T][i];
            q[T] = i;
        }
    }

    /* 4. Backtracking on the path (sequence of states) */
    for (t = T - 1; t >= 1; t--)
		q[t] = psi[t+1][q[t+1]];
}

/*  Function with the Viterbi algorithm, adapted to the problem of recognition. */
void ViterbiLog_C(HMM *phmm, int T, int *O, double **delta, int **psi, int *q, double *pprob)
{
    int i, j;
    int t;

    int     maxvalind;
    double  maxval, val;

    /*  1. Initialisation  */

   	for (i = 1; i <= phmm->N; i++)
	{
	    if(phmm->pi[i]*phmm->B[i][O[1]] != 0)
            delta[1][i] = log (phmm->pi[i] * (phmm->B[i][O[1]]));
        else
            delta[1][i] = -VITHUGE2;
        psi[1][i] = 0;
	}

    /*  2. Recursion */
    for (t = 2; t <= T; t++)
    {
        for (j = 1; j <= phmm->N; j++)
        {
            maxval = -VITHUGE;
            maxvalind = 1;
            for (i = 1; i <= phmm->N; i++)
            {
                if(phmm->A[i][j]*(phmm->B[j][O[t]]) != 0)
                    val = delta[t-1][i] + log(phmm->A[i][j]*(phmm->B[j][O[t]]));
                else
                    val = delta[t-1][i];
                if (val > maxval)
                {
                    maxval = val;
                    maxvalind = i;
                }
            }
            delta[t][j] = maxval;
            psi[t][j] = maxvalind;
        }
    }

    /*  3. Termination */
    *pprob = -VITHUGE;
    q[T] = 1;
    for (i = 1; i <= phmm->N; i++)
    {
        if (delta[T][i] > *pprob)
        {
            *pprob = delta[T][i];
            q[T] = i;
        }
    }

    /* 4. Backtracking on the path (sequence of states) */
    for (t = T - 1; t >= 1; t--)
		q[t] = psi[t+1][q[t+1]];
}

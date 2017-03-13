/*
**      Author: Tapas Kanungo, kanungo@cfar.umd.edu
**      Date:   15 December 1997
**      File:   baumwelch.c
**      Purpose: Baum-Welch algorithm for estimating the parameters
**              of a HMM model, given an observation sequence.
**      Organization: University of Maryland
**
**	Update:
**	Author: Tapas Kanungo
**	Date:	19 April 1999
**	Purpose: Changed the convergence criterion from ratio
**		to absolute value.
**
**      $Id: baumwelch.c,v 1.6 1999/04/24 15:58:43 kanungo Exp kanungo $
*/

/*
        Change
        Author:     Pietro Russo
        Date:       12/03/2011
        File:       baum.h
        Purpose:    Added the BaumWelch_C function that fits the Baum-Welch algorithm
                    for the speech recognition problem.
*/


/*  Solution to the problem of learning:
    how to adjust the parameters of the hidden Markov model
    lambda = (A, B, pi), so as to maximize P(O|lambda)
*/
#include <stdio.h>
#include "nrutil.h"
#include "hmm.h"
#include <math.h>

static char rcsid[] = "$Id: baumwelch.c,v 1.6 1999/04/24 15:58:43 kanungo Exp kanungo $";

#define DELTA 0.001
#define EPSILON 0.0001

/*  Function that invokes the Baum-Welch algorithm */
void BaumWelch(HMM *phmm, int T, int *O, double **alpha, double **beta,	double **gamma, int *pniter, double *plogprobinit, double *plogprobfinal)
{
	int	i, j, k;
	int	t, l = 0;

	double	logprobf, logprobb,  threshold;
	double	numeratorA, denominatorA;
	double	numeratorB, denominatorB;

	double ***xi, *scale;
	double delta, deltaprev, logprobprev;

	deltaprev = 10e-70;

	xi = AllocXi(T, phmm->N);
	scale = dvector(1, T);

	ForwardWithScale(phmm, T, O, alpha, scale, &logprobf);
	*plogprobinit = logprobf; /* log P(O |intial model) */
	BackwardWithScale(phmm, T, O, beta, scale, &logprobb);
	ComputeGamma(phmm, T, alpha, beta, gamma);
	ComputeXi(phmm, T, O, alpha, beta, xi);
	logprobprev = logprobf;

	do
	{
        /*  Re-estimates the frequency of the state i over the time t=1 */
        /*  The robability of being in state S_i at time t */
		for (i = 1; i <= phmm->N; i++)
			phmm->pi[i] = .001 + .999*gamma[1][i];

        /*  Re-estimating the transition probability matrix and the different
            observation symbols per state */

		for (i = 1; i <= phmm->N; i++)
		{   /*  a^_ij = sum xi_t(i,j) / sum gamma_t(i) */
			denominatorA = 0.0;
			for (t = 1; t <= T - 1; t++)
				denominatorA += gamma[t][i];

			for (j = 1; j <= phmm->N; j++)
			{
				numeratorA = 0.0;
				for (t = 1; t <= T - 1; t++)
					numeratorA += xi[t][i][j];
				phmm->A[i][j] = .001 + .999*numeratorA/denominatorA;
			}

            /*  b^_j(k) = sum gamma_t(j) / sum gamma_t(j) */
			denominatorB = denominatorA + gamma[T][i];
			for (k = 1; k <= phmm->M; k++)
			{
				numeratorB = 0.0;
				for (t = 1; t <= T; t++)
				{
					if (O[t] == k)
						numeratorB += gamma[t][i];
				}
				phmm->B[i][k] = .001 + .999*numeratorB/denominatorB;
			}
		}

		ForwardWithScale(phmm, T, O, alpha, scale, &logprobf);
		BackwardWithScale(phmm, T, O, beta, scale, &logprobb);
		ComputeGamma(phmm, T, alpha, beta, gamma);
		ComputeXi(phmm, T, O, alpha, beta, xi);

        /*  Calculates the difference between the logarithmic probability
            of two iterations, the current and the previous */
		delta = logprobf - logprobprev;
		logprobprev = logprobf;
        /*  Increments the counter that tracks the number of iterations */
		l++;
	}while (delta > DELTA);
	/*  Exits from the do while if the logarithmic probability does not change much */

	*pniter = l;
	/*  log P(O|estimated models) */
	*plogprobfinal = logprobf;
	FreeXi(xi, T, phmm->N);
	free_dvector(scale, 1, T);
}

/*  Function that calculates the gamma as:
    gamma_t(i) = { alpha_t(i) * beta_t(i) } / P(O|lambda) */
void ComputeGamma(HMM *phmm, int T, double **alpha, double **beta, double **gamma)
{
	int 	i, j;
	int	t;
	double	denominator;

	for (t = 1; t <= T; t++)
	{
		denominator = 0.0;
		for (j = 1; j <= phmm->N; j++)
		{
			gamma[t][j] = alpha[t][j]*beta[t][j];
			denominator += gamma[t][j];
		}

		for (i = 1; i <= phmm->N; i++)
			gamma[t][i] = gamma[t][i]/denominator;
	}
}

/*  Function that calculates the xi as:
    xi_t(i,j) = { alpha_t(i) * beta_t+1(j) * a_ij * b_j } / P(O|lambda) */
void ComputeXi(HMM* phmm, int T, int *O, double **alpha, double **beta,	double ***xi)
{
	int i, j;
	int t;
	double sum;

	for (t = 1; t <= T - 1; t++)
	{
		sum = 0.0;
		for (i = 1; i <= phmm->N; i++)
			for (j = 1; j <= phmm->N; j++)
			{
				xi[t][i][j] = alpha[t][i]*beta[t+1][j]*(phmm->A[i][j])*(phmm->B[j][O[t+1]]);
				sum += xi[t][i][j];
			}

		for (i = 1; i <= phmm->N; i++)
			for (j = 1; j <= phmm->N; j++)
				xi[t][i][j] /= sum;
	}
}

/*  Allocates the array of 3 dimensions xi */
double *** AllocXi(int T, int N)
{
	int t;
	double ***xi;

	xi = (double ***) malloc(T*sizeof(double **));

	xi --;

	for (t = 1; t <= T; t++)
		xi[t] = dmatrix(1, N, 1, N);
	return xi;
}

/*  Free the allocated memory space */
void FreeXi(double *** xi, int T, int N)
{
	int t;

	for (t = 1; t <= T; t++)
		free_dmatrix(xi[t], 1, N, 1, N);

	xi ++;
	free(xi);
}

/*  Added the function to re-estimate the values of HMM, according to the problem of speech recognition */
void BaumWelch_C(HMM *phmm, int T, int *O, double **alpha, double **beta,	double **gamma, int *pniter, double *plogprobinit, double *plogprobfinal)
{
	int	i, j, k;
	int	t, l = 0;

	double	logprobf, logprobb,  threshold;
	double	numeratorA, denominatorA;
	double	numeratorB, denominatorB;

	double ***xi, *scale;
	double delta, deltaprev, logprobprev;

    int temp=0;
    double somma=0;
    double val;

	deltaprev = 10e-70; // non viene utilizzata?!?!?

	xi = AllocXi(T, phmm->N);
	scale = dvector(1, T);

	Forward(phmm, T, O, alpha, &logprobf);
	*plogprobinit =logprobf; /* log P(O |intial model) */
	Backward(phmm, T, O, beta, &logprobb);
	ComputeGamma(phmm, T, alpha, beta, gamma);
	ComputeXi(phmm, T, O, alpha, beta, xi);
	logprobprev = logprobf;

	do
	{
        /*  Re-estimating the frequency of the state i over the time t=1 */
        /*  The probability of being in state S_i at time t */
		for (i = 1; i <= phmm->N; i++)
			//phmm->pi[i] = .001 + .999*gamma[1][i];
			phmm->pi[i] = gamma[1][i];

        /*  Re-estimating the transition probability matrix and the different observation symbols per state */

		for (i = 1; i <= phmm->N; i++)
		{   /*  a^_ij = sum xi_t(i,j) / sum gamma_t(i) */
			denominatorA = 0.0;
			for (t = 1; t <= T - 1; t++)
				denominatorA += gamma[t][i];

			for (j = 1; j <= phmm->N; j++)
			{
				numeratorA = 0.0;
				for (t = 1; t <= T - 1; t++)
					numeratorA += xi[t][i][j];
				phmm->A[i][j] = numeratorA/denominatorA;
			}

            /*  b^_j(k) = sum gamma_t(j) / sum gamma_t(j) */
			denominatorB = denominatorA + gamma[T][i];
			for (k = 1; k <= phmm->M; k++)
			{
				numeratorB = 0.0;
				for (t = 1; t <= T; t++)
				{
					if (O[t] == k)
						numeratorB += gamma[t][i];
				}
				phmm->B[i][k] = numeratorB/denominatorB;
				if(phmm->B[i][k] < EPSILON)
				{
				    phmm->B[i][k] = EPSILON;
				    ++temp;
				}
				else
                    somma+=phmm->B[i][k];
			}

            if(temp != phmm->M && temp !=0)
			{
			    val = 1-(EPSILON*temp);
			    for (k = 1; k <= phmm->M; k++)
                {
                    if(phmm->B[i][k] != EPSILON)
                    {
                        phmm->B[i][k] = phmm->B[i][k] * val / somma ;
                    }
                }
			}

			if(temp == phmm->M)
			{
			    somma=temp*EPSILON;
			    for (k = 1; k <= phmm->M; k++)
                    phmm->B[i][k] = phmm->B[i][k] / somma;
			}
		}

        Forward(phmm, T, O, alpha, &logprobf);
		Backward(phmm, T, O, beta, &logprobb);
		ComputeGamma(phmm, T, alpha, beta, gamma);
		ComputeXi(phmm, T, O, alpha, beta, xi);

        /*  Calculates the difference between the logarithmic probability of two iterations,
            the current and the previous */
        delta = log(logprobf) - log(logprobprev);
		logprobprev = logprobf;
        /*  Increments the counter that tracks the number of iterations */

		l++;

	}while (delta > DELTA);
	/*  Exits the do while if the logarithmic probability does not change much */

	*pniter = l;
	/*  log P(O|estimated models) */
	*plogprobfinal = logprobf;
	FreeXi(xi, T, phmm->N);
	free_dvector(scale, 1, T);
}

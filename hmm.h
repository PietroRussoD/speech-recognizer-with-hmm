/*
**      Author:     Tapas Kanungo, kanungo@cfar.umd.edu
**      Date:       15 December 1997
**      File:       hmm.h
**      Purpose:    datastructures used for HMM.
**      Organization: University of Maryland
**
**	    Update:
**	    Author:     Tapas Kanungo
**	    Purpose:    include <math.h>. Not including this was
**		            creating a problem with forward.c
**      $Id: hmm.h,v 1.9 1999/05/02 18:38:11 kanungo Exp kanungo $
*/

/*
        Change
        Author:     Pietro Russo
        Date:       11/03/2011
        File:       hmm.h
        Purpose:    Added functions that allows to save the HMM on a binary file.
                    Added functions that initializes the matrices of the model
                    according to a predefined sequence of states.
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "lista.h"

/*  Structure of the Hidden Markovian Model HMM */
typedef struct {
    /* N number of states;  Q={1,2,...,N} */
	int N;
	/* M number of symbols of the observation; V={1,2,...,M}*/
	int M;
	/*  A[1..N][1..N]. a[i][j] is the distribution of the transition probability of change
        from state i at time t to state j at time t + 1 */
	double	**A;
	/*  B[1..N][1..M]. b[j][k] is the probability distribution of observation k symbols in state j */
	double	**B;
	/*  pi[1..N] pi[i] is the initial distribution of states */
	double	*pi;
} HMM;

void ReadHMM(FILE *fp, HMM *phmm);
void PrintHMM(FILE *fp, HMM *phmm);
void InitHMM(HMM *phmm, int N, int M, int seed);
void CopyHMM(HMM *phmm1, HMM *phmm2);
void FreeHMM(HMM *phmm);


void InitHMM_SC1(HMM *phmm, int N, int M, int seed);
void InitHMM_SC2(HMM *phmm, int N, int M, int seed);

void Salva_HMM(char nome_file[], HMM *phmm);
void Carica_HMM(char nome_file[], HMM *phmm, int N, int M);

void ReadSequence(FILE *fp, int *pT, int **pO);


void PrintSequence(FILE *fp, int T, int *O);
void GenSequenceArray(HMM *phmm, int seed, int T, int *O, int *q);
int GenInitalState(HMM *phmm);
int GenNextState(HMM *phmm, int q_t);
int GenSymbol(HMM *phmm, int q_t);

void Forward(HMM *phmm, int T, int *O, double **alpha, double *pprob);
void ForwardWithScale(HMM *phmm, int T, int *O, double **alpha,
        double *scale, double *pprob);
void Backward(HMM *phmm, int T, int *O, double **beta, double *pprob);
void BackwardWithScale(HMM *phmm, int T, int *O, double **beta,
        double *scale, double *pprob);
void BaumWelch(HMM *phmm, int T, int *O, double **alpha, double **beta,
        double **gamma, int *niter,
	double *plogprobinit, double *plogprobfinal);

void BaumWelch_C(HMM *phmm, int T, int *O, double **alpha, double **beta,
        double **gamma, int *niter,
	double *plogprobinit, double *plogprobfinal);

double *** AllocXi(int T, int N);
void FreeXi(double *** xi, int T, int N);
void ComputeGamma(HMM *phmm, int T, double **alpha, double **beta,
        double **gamma);
void ComputeXi(HMM* phmm, int T, int *O, double **alpha, double **beta,
        double ***xi);
void Viterbi(HMM *phmm, int T, int *O, double **delta, int **psi,
        int *q, double *pprob);
void ViterbiLog(HMM *phmm, int T, int *O, double **delta, int **psi,
        int *q, double *pprob);

void ViterbiLog_C(HMM *phmm, int T, int *O, double **delta, int **psi,
        int *q, double *pprob);

/*  Functions belonging to the random number generator */
unsigned int hmmgetseed(void);
void hmmsetseed(unsigned int seed);
double hmmgetrand(void);

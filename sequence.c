/*
**      Author: Tapas Kanungo, kanungo@cfar.umd.edu
**      Date:   22 February 1998
**      File:   sequence.c
**      Purpose: Routines for generating, reading and writing sequence of
**		observation symbols.
**      Organization: University of Maryland
**
**	Update:
**	Author: Tapas Kanungo
**	Purpose: To make calls to generic random number generators
**		and to change the seed everytime the software is executed.
**
**      $Id: sequence.c,v 1.2 1998/02/23 06:19:41 kanungo Exp kanungo $
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "nrutil.h"
#include "hmm.h"

static char rcsid[] = "$Id: sequence.c,v 1.2 1998/02/23 06:19:41 kanungo Exp kanungo $";

/*  Function that generates a sequence of random observations. */
void GenSequenceArray(HMM *phmm, int seed, int T, int *O, int *q)
{
    /*  T is the number of the sequence of observations O.
        Each observation is denoted by Ot */
    /*  q is the vector containing the initial values of the states */
    int     t = 1;

	hmmsetseed(seed);

        q[1] = GenInitalState(phmm);
        O[1] = GenSymbol(phmm, q[1]);

        for (t = 2; t <= T; t++) {
                q[t] = GenNextState(phmm, q[t-1]);
                O[t] = GenSymbol(phmm, q[t]);
        }
}

/*  Select a random initial state */
int GenInitalState(HMM *phmm)
{
    double val, accum;
    int i, q_t;
    /*  Attache to val a random value between [0,1) */
    val = hmmgetrand();
    accum = 0.0;
    /*  Set to q_t the number of states */
    q_t = phmm->N;
    /*  For each state */
    for (i = 1; i <= phmm->N; i++)
    {
        /*  If the randomly generated value is less than the sum of the i-th initial value
            and the sum of the values of the initial up to the i-th */
        if (val < phmm->pi[i] + accum)
        {
            q_t = i;
            /*  Exit from the for cycle */
            break;
        }
        else
        {
            accum += phmm->pi[i];
        }
    }
    return q_t;
}
/*  Select the next state */
int GenNextState(HMM *phmm, int q_t)
{
    double val, accum;
    int j, q_next;

    val = hmmgetrand();
    accum = 0.0;
    q_next = phmm->N;
    /*  For each column j of the probability distribution A of transition to change from q_t state to j state */
    for (j = 1; j <= phmm->N; j++)
    {
        if ( val < phmm->A[q_t][j] + accum )
        {
            q_next = j;
            break;
        }
        else
            accum += phmm->A[q_t][j];
    }
    return q_next;
}
/*  Select a symbol of observation */
int GenSymbol(HMM *phmm, int q_t)
{
    double val, accum;
    int j, o_t;

    val = hmmgetrand();
    accum = 0.0;
    o_t = phmm->M;
    /*  For each column j of the probability distribution B of observation symbols of the state in q_t */
    for (j = 1; j <= phmm->M; j++)
    {
        /*  If val is less than the sum of the value and the sum of the values examined until now */
        if ( val < phmm->B[q_t][j] + accum )
        {
            o_t = j;
            break;
        }
        else
            accum += phmm->B[q_t][j];
    }
    return o_t;
}
/*  Given a file reads the sequence of observations */
void ReadSequence(FILE *fp, int *pT, int **pO)
{
        int *O;
        int i;

        fscanf(fp, "T= %d\n", pT);
        O = ivector(1,*pT);
        for (i=1; i <= *pT; i++)
                fscanf(fp,"%d", &O[i]);
        *pO = O;
}
/*  Save in a file the sequence of observations */
void PrintSequence(FILE *fp, int T, int *O)
{
        int i;

        fprintf(fp, "T= %d\n", T);
        for (i=1; i <= T; i++)
                fprintf(fp,"%d ", O[i]);
	printf("\n");

}

/*
**      Author:     Tapas Kanungo, kanungo@cfar.umd.edu
**      Date:       15 December 1997
**      File:       hmmutils.c
**      Purpose:    utilities for reading, writing HMM stuff.
**      Organization: University of Maryland
**
**      $Id: hmmutils.c,v 1.4 1998/02/23 07:51:26 kanungo Exp kanungo $
*/

/*
        Change
        Author:     Pietro Russo
        Date:       11/03/2011
        File:       hmmutils.c
        Purpose:    Addition of functions that saves the HMM on a binary file.
                    Addition of the functions that initializes the pattern matrices
                    according to a predefined sequence of states.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "nrutil.h"
#include "hmm.h"
static char rcsid[] = "$Id: hmmutils.c,v 1.4 1998/02/23 07:51:26 kanungo Exp kanungo $";

/*  Opens the file, it reads the data inside it and inserts it in the pointer to HMM */
void ReadHMM(FILE *fp, HMM *phmm)
{
	int i, j, k;
    /*  Reads the number of observation symbols */
	fscanf(fp, "M= %d\n", &(phmm->M));
    /*  Reads the number of states */
	fscanf(fp, "N= %d\n", &(phmm->N));
    /*  Reads the values of the distribution of transition probabilities */
	fscanf(fp, "A:\n");
	phmm->A = (double **) dmatrix(1, phmm->N, 1, phmm->N);
	for (i = 1; i <= phmm->N; i++) {
		for (j = 1; j <= phmm->N; j++) {
			fscanf(fp, "%lf", &(phmm->A[i][j]));
		}
		fscanf(fp,"\n");
	}
    /*  Reads the values of the distribution probability of observing k symbols in state j */
	fscanf(fp, "B:\n");
	phmm->B = (double **) dmatrix(1, phmm->N, 1, phmm->M);
	for (j = 1; j <= phmm->N; j++) {
		for (k = 1; k <= phmm->M; k++) {
			fscanf(fp, "%lf", &(phmm->B[j][k]));
		}
		fscanf(fp,"\n");
	}
    /*  Reads the values of the initial distribution of the states */
	fscanf(fp, "pi:\n");
	phmm->pi = (double *) dvector(1, phmm->N);
	for (i = 1; i <= phmm->N; i++)
		fscanf(fp, "%lf", &(phmm->pi[i]));

}
/*  Free the allocated dynamically memory */
void FreeHMM(HMM *phmm)
{
	free_dmatrix(phmm->A, 1, phmm->N, 1, phmm->N);
	free_dmatrix(phmm->B, 1, phmm->N, 1, phmm->M);
	free_dvector(phmm->pi, 1, phmm->N);
}

/*  InitHMM() This function initializes the matrices A, B and the vector pi with random values.
    If you do the BaumWelch algorithm might behave oddly. */
void InitHMM(HMM *phmm, int N, int M, int seed)
{
	int i, j, k;
	double sum;

	/*  Initialize the random number generator */
	hmmsetseed(seed);
    /*  Set the number of observation symbols */
    phmm->M = M;
    /*  Set the number of states */
    phmm->N = N;
    /*  Allocate the matrix A, that is the distribution of the transition probability
        from state i at time t to state j at time t+1 */
    phmm->A = (double **) dmatrix(1, phmm->N, 1, phmm->N);
    /*  Fill the array with random numbers and then normalizes them,
        in such a way that the sum of each row is 1. */
    /*  For each state i */
    for (i = 1; i <= phmm->N; i++)
    {
		sum = 0.0;
        for (j = 1; j <= phmm->N; j++)
        {
            phmm->A[i][j] = hmmgetrand();
			sum += phmm->A[i][j];
		}
        for (j = 1; j <= phmm->N; j++)
            phmm->A[i][j] /= sum;
	}
    /*  Allocate the matrix B, which is the probability distribution of observation k symbols in state j */
    phmm->B = (double **) dmatrix(1, phmm->N, 1, phmm->M);
    /*  For each state j */
    for (j = 1; j <= phmm->N; j++)
    {
        sum = 0.0;
        for (k = 1; k <= phmm->M; k++)
        {
            phmm->B[j][k] = hmmgetrand();
			sum += phmm->B[j][k];
		}
        for (k = 1; k <= phmm->M; k++)
			phmm->B[j][k] /= sum;
	}
    /*  Allocate the vector, which is the initial distribution of the states */
    phmm->pi = (double *) dvector(1, phmm->N);
    sum = 0.0;
    for (i = 1; i <= phmm->N; i++)
    {
        phmm->pi[i] = hmmgetrand();
		sum += phmm->pi[i];
	}
    for (i = 1; i <= phmm->N; i++)
		phmm->pi[i] /= sum;
}

/*  This function copies a hidden Markov model to another */
void CopyHMM(HMM *phmm1, HMM *phmm2)
{
    int i, j, k;

    phmm2->M = phmm1->M;
    phmm2->N = phmm1->N;

    phmm2->A = (double **) dmatrix(1, phmm2->N, 1, phmm2->N);
    for (i = 1; i <= phmm2->N; i++)
        for (j = 1; j <= phmm2->N; j++)
            phmm2->A[i][j] = phmm1->A[i][j];

    phmm2->B = (double **) dmatrix(1, phmm2->N, 1, phmm2->M);
        for (j = 1; j <= phmm2->N; j++)
            for (k = 1; k <= phmm2->M; k++)
                phmm2->B[j][k] = phmm1->B[j][k];

    phmm2->pi = (double *) dvector(1, phmm2->N);
        for (i = 1; i <= phmm2->N; i++)
            phmm2->pi[i] = phmm1->pi[i];
}

/*  This function writes a file with the hidden Markov model values passed into the function */
void PrintHMM(FILE *fp, HMM *phmm)
{
    int i, j, k;

	fprintf(fp, "M= %d\n", phmm->M);
	fprintf(fp, "N= %d\n", phmm->N);

	fprintf(fp, "A:\n");
    for (i = 1; i <= phmm->N; i++)
    {
        for (j = 1; j <= phmm->N; j++)
        {
            fprintf(fp, "%f ", phmm->A[i][j] );
		}
		fprintf(fp, "\n");
	}

	fprintf(fp, "B:\n");
    for (j = 1; j <= phmm->N; j++)
    {
        for (k = 1; k <= phmm->M; k++)
        {
            fprintf(fp, "%f ", phmm->B[j][k]);
		}
		fprintf(fp, "\n");
	}

	fprintf(fp, "pi:\n");
    for (i = 1; i <= phmm->N; i++)
    {
		fprintf(fp, "%f ", phmm->pi[i]);
	}
	fprintf(fp, "\n\n");
}

/*  Function that saves the hidden Markov model of a binary file */
void Salva_HMM(char nome_file[], HMM *phmm)
{
    FILE *fp;
    int i;
    fp=fopen(nome_file,"wb");
    if(fp==NULL)
         errore("Errore di apertura del file.");
    /*  Memorize the first matrix A of size NxN
        that is the distribution of the transition probability
        from state i at time t to state j at time t+1 */
    size_t dimensione = sizeof(double);
    size_t elementi = phmm->N;
    size_t ritorno;

    for(i=1; i<phmm->N; i++)
    {
        ritorno = fwrite(phmm->A[i], dimensione, elementi, fp);
        if (ritorno != elementi)
            errore("Errore scrittura su file.");
    }


    /*  Memorize the matrix B of size NxM
        that is the probability distribution of observation k symbols in state j */
    elementi = phmm->M;

    for(i=1; i<phmm->N; i++)
    {
        ritorno = fwrite(phmm->B[i], dimensione, elementi, fp);
        if (ritorno != elementi)
            errore("Errore scrittura su file.");
    }

    /*  Finally memorize the vector pi of size N which is the initial distribution of the states */
    elementi = phmm->N;

    ritorno = fwrite(phmm->pi, dimensione, elementi, fp);
        if (ritorno != elementi)
            errore("Errore scrittura su file.");

    fclose(fp);
}

/*  Function that loads a binary file on hidden Markov model */
void Carica_HMM(char nome_file[], HMM *phmm, int N, int M)
{
    FILE *fp;
    int i;

    phmm->N = N;
    phmm->M = M;

    phmm->A = (double **) dmatrix(1, phmm->N, 1, phmm->N);
    phmm->B = (double **) dmatrix(1, phmm->N, 1, phmm->M);
    phmm->pi = (double *) dvector(1, phmm->N);

    fp=fopen(nome_file,"rb");
    if(fp==NULL)
         errore("Errore di apertura del file.");
    /*  Memorize the first matrix A of size NxN
        that is the distribution of the transition probability
        from state i at time t to state j at time t+1 */
    size_t dimensione = sizeof(double);
    size_t elementi = phmm->N;
    size_t ritorno;

    for(i=1; i<phmm->N; i++)
    {
        ritorno = fread(phmm->A[i], dimensione, elementi, fp);
        if (ritorno != elementi)
            errore("Errore scrittura su file.");
    }

    /*  Memorize the matrix B of size NxM
        that is the probability distribution of observation k symbols in state j */
    elementi = phmm->M;

    for(i=1; i<phmm->N; i++)
    {
        ritorno = fread(phmm->B[i], dimensione, elementi, fp);
        if (ritorno != elementi)
            errore("Errore scrittura su file.");
    }

    /*  Finally memorize the vector pi of size N which is the initial distribution of the states */
    elementi = phmm->N;

    ritorno = fread(phmm->pi, dimensione, elementi, fp);
        if (ritorno != elementi)
            errore("Errore scrittura su file.");

    fclose(fp);
}

/*  Function that initializes the hidden Markov model based on the first type of serial constraint */
void InitHMM_SC1(HMM *phmm, int N, int M, int seed)
{
	int i, j, k;
	double sum;

	/*  Initialize the random number generator */
	hmmsetseed(seed);
    /*  Set the number of observation symbols */
    phmm->M = M;
    /*  Set the number of states */
    phmm->N = N;
    /*  Allocate the matrix A, that is the distribution of the transition probability
        from state i at time t to state j at time t+1 */
    phmm->A = (double **) dmatrix(1, phmm->N, 1, phmm->N);
    /*  Fill the array with random numbers and then normalizes them,
        in such a way that the sum of each is 1. */
    /*  For each state i */
    for (i = 1; i <= phmm->N; i++)
    {
		sum = 0.0;
        for (j = 1; j <= phmm->N; j++)
        {
            if(j<i || j>=i+3)
                phmm->A[i][j] = 0.0;
            else
                phmm->A[i][j] = hmmgetrand();
			sum += phmm->A[i][j];
		}
        for (j = 1; j <= phmm->N; j++)
            phmm->A[i][j] /= sum;
	}
    /*  Allocate the matrix B, which is the probability distribution of observation k symbols in state j */
    phmm->B = (double **) dmatrix(1, phmm->N, 1, phmm->M);
    /*  Per ogni stato j */
    for (j = 1; j <= phmm->N; j++)
    {
        sum = 0.0;
        for (k = 1; k <= phmm->M; k++)
        {
            phmm->B[j][k] = hmmgetrand();
			sum += phmm->B[j][k];
		}
        for (k = 1; k <= phmm->M; k++)
			phmm->B[j][k] /= sum;
	}
    /*  Allocate the vector pi, which is the initial distribution of the states */
    phmm->pi = (double *) dvector(1, phmm->N);
    sum = 0.0;
    phmm->pi[1] = 1;
    for (i = 2; i <= phmm->N; i++)
    {
        phmm->pi[i] = 0;
	}
}

/*  Function that initializes the hidden Markov model based on the second type of serial constraint */
void InitHMM_SC2(HMM *phmm, int N, int M, int seed)
{
	int i, j, k;
	double sum;

	/*  Initialize the random number generator */
	hmmsetseed(seed);
    /*  Set the number of observation symbols */
    phmm->M = M;
    /*  Set the number of states */
    phmm->N = N;
    /*  Allocate the matrix A, that is the distribution of the transition probability
        from state i at time t to state j at time t+1 */
    phmm->A = (double **) dmatrix(1, phmm->N, 1, phmm->N);
    /*  Fill the array with random numbers and then normalizes them,
        in such a way that the sum of each is 1. */
    /*  For each state i */
    for (i = 1; i <= phmm->N; i++)
    {
		sum = 0.0;
		/*  Ho un certo valore se vado nello stato j */
        for (j = 1; j <= phmm->N; j++)
        {
            if(j<i || j>=i+2)
                phmm->A[i][j] = 0.0;
            else
                phmm->A[i][j] = hmmgetrand();
			sum += phmm->A[i][j];
		}
        for (j = 1; j <= phmm->N; j++)
            phmm->A[i][j] /= sum;
	}
    /*  Allocate the matrix B, which is the probability distribution of observation k symbols in state j */
    phmm->B = (double **) dmatrix(1, phmm->N, 1, phmm->M);
    /*  Per ogni stato j */
    for (j = 1; j <= phmm->N; j++)
    {
        sum = 0.0;
        for (k = 1; k <= phmm->M; k++)
        {
            phmm->B[j][k] = hmmgetrand();
			sum += phmm->B[j][k];
		}
        for (k = 1; k <= phmm->M; k++)
			phmm->B[j][k] /= sum;
	}
    /*  Allocate the vector pi, which is the initial distribution of the states */
    phmm->pi = (double *) dvector(1, phmm->N);
    sum = 0.0;
    phmm->pi[1] = 1;
    for (i = 2; i <= phmm->N; i++)
    {
        phmm->pi[i] = 0;
	}
}

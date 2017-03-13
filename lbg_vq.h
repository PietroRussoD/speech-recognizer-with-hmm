/*
        Author:     Pietro Russo
        Date:       18/02/2011
        File:       lbg_vq.h
        Purpose:    Functions that are used to determine the codebook and the codevectors.
*/
#include "matr_vett.h"

double **lbg(double **training, int M, int K, int N);
int *vq(double **codebook, double **dati, int M, int K, int N);

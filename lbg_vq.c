/*
        Change
        Author:     Pietro Russo
        Date:       18/02/2011
        File:       lbg_vq.c
        Purpose:    Functions that are used to determine the codebook and the codevectors.
*/

#include "lbg_vq.h"

/*  M is the number of source vectors of the training sequence
    T = { x_1, x_2, ... x_M }*/
/*  K is the size of the source vectors of the training sequence
    x_m = (x_m,1, x_m,2, ..., x_m,K)    with m = 1, 2, ..., M */
/*  N is the desired number of codevectors, which together form
    the codebook C = {c_1, c_2, ..., c_N} where every codevector has the dimension K
    c_n = (c_n,1, c_n,2, ..., c_n,K)    with n = 1, 2, ..., N */
/*  The relationship M >= 1000*N should normally be valid */

/*  Function that creates the codebook starting from a training set */
double **lbg(double **training, int M, int K, int N)
{
    double epsilon= 0.001;;

    int k,m,n;
    int n_codevectors;
    double **codebook = matrice_d(N, K);
    double **codebook_star = matrice_d(N, K);
    int quanti[N], indice;
    int flag_epsilon=0;
    double distanza_min, distanza_temp;
    double distorsione, distorsione_star, rapporto_distorsione;

    /*  INITIALIZATION */
    /*  The first codevectors calculation */
    n_codevectors=1;

    for (k=1; k<K; k++)
    {
        codebook_star[1][k]=0.0;

        for (m=1; m<M; m++)
            codebook_star[1][k] += training[m][k];

        codebook_star[1][k] /= (float) M;
    }

    /*  Until is not reached the desired number of codevectors */
    while (n_codevectors < N-1)
    {
        /* DIVISION */
        for(n=1; n<=n_codevectors; n++)
            for(k=1; k<K; k++)
            {
                codebook[n][k] = (1+epsilon)*codebook_star[n][k];
                codebook[n+n_codevectors][k] = (1-epsilon)*codebook_star[n][k];
            }

        n_codevectors *= 2;

        /*  Set a high value to distorsione_star so the cycle is repeated more than once */
        distorsione_star = 10.0E15;

        flag_epsilon = 0;
        /*  Until rapporto_distorsione is greater than epsilon */
        while ( flag_epsilon == 0)
        {
            /*  ITERATION */
            distorsione=0.0;

            /*  Preparation of vectors */
            for(n=1; n<=n_codevectors; n++)
            {
                for(k=1; k<K; k++)
                    codebook_star[n][k]=0.0;
                quanti[n]=0;
            }

            for(m=1; m<M; m++)
            {
                distanza_min=0.0;

                /*  Nearest Neighbor Condition
                    Associate the closest codevector */

                /*  Because there is a comparison later, it imposes the minimum distance values for n=1 */
                for(k=1; k<K; k++)
                {
                    /*  Take in consideration codebook[n][k] with n = 1 */
                    distanza_min += pow( (training[m][k] - codebook[1][k]), 2 );
                }
                /*  Store the index n=1 */
                indice=1;

                /*  Start from n=2, because n=1 it has already processed before */
                for(n=2; n<=n_codevectors; n++)
                {
                    distanza_temp=0.0;
                    for(k=1; k<K; k++)
                        distanza_temp += pow( (training[m][k] - codebook[n][k]), 2 );

                    if(distanza_temp < distanza_min)
                    {
                        distanza_min = distanza_temp;
                        indice = n;
                    }
                }
                /*  Increase of one the n-th position of the codevector that satisfied the previous criteria */
                ++quanti[indice];

                /*  Start to calculate the codebook that will be updated */
                for(k=1; k<K; k++)
                    codebook_star[indice][k] += training[m][k];

                distorsione = distorsione+distanza_min;
            }

            /*  Centroid Condition
                The codevector c_n should be the average of all vectors of the training that are in the codified region S_n */
            for(n=1; n<=n_codevectors; n++)
            {
                if(quanti[n] > 0)
                {
                    for(k=1; k<K; k++)
                    {
                        codebook_star[n][k] /= (float) quanti[n];
                        codebook[n][k] = codebook_star[n][k];
                    }
                }
            }

            distorsione /= (float) M*K;

            rapporto_distorsione = (distorsione_star-distorsione)/distorsione_star;

            flag_epsilon = 1;

            if(rapporto_distorsione > 0.001)
            {
                flag_epsilon = 0;
                distorsione_star = distorsione;
            }
        }
    }

    return codebook;
}

/*  Function that quantises the data passed on the function of codebook */
int *vq(double **codebook, double **dati, int n_fin, int K, int N)
{
    double distanza_min, distanza_temp;
    int m, k, n, indice;
    int *quali = vettore_i(n_fin+1);

    for(m=0; m<n_fin; m++)
    {
        /*  Nearest Neighbor Condition
            Associate the closest codevector */
        distanza_min=0.0;

        for(k=1; k<K; k++)
        {
            distanza_min += pow( (dati[m][k-1]-codebook[1][k]), 2);
        }
        indice=1;

        for(n=2; n<N; n++)
        {
            distanza_temp=0.0;
            for(k=1; k<K; k++)
                distanza_temp += pow( (dati[m][k-1]-codebook[n][k]), 2);

            if(distanza_temp < distanza_min)
            {
                distanza_min = distanza_temp;
                indice = n;
            }
        }
        quali[m+1] = indice;
    }
    return quali;
}


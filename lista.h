/*
        Author:     Pietro Russo
        Date:       03/03/2011
        File:       lista.h
        Purpose:    Functions that are used to implement a list containing
                    voice signals used in the training set.
                    The signals have associated their HMM that will serve
                    in the recognition process to determine the word.
*/
#include "matr_vett.h"

/*  The windows is defined in milliseconds */
#define DIM_FINESTRA 25
#define DIM_PASSO 5
/*  The frequency is defined in Hertz */
#define FREQUENZA_MAX 5000
/*  Number of considered cepstral coefficients */
#define N_COEFF_CEP 13

/*  Data structure for the list of words */
typedef struct
{
    char nome_file_audio[32];
    int n_finestre;
    int frequenza_campionamento;
    int campioni_segnale;
    char nome_file_car[32];
    char nome_file_seq[32];
    char nome_file_hmm[32];
    struct lista_parole *prossimo;
}lista_parole;

lista_parole *aggiungi(char nome_file[32], lista_parole *lista);
void mostra(lista_parole *lista);

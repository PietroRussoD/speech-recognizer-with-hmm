/*
        Author:     Pietro Russo
        Data:       15/02/2011
        File:       analisi_caratteristiche.h
        Purpose:    Functions that determines the values uses to proceed
                    with the speech recognition.
*/

#include "matr_vett.h"
#include <math.h>
#include "sndfile.h"
#include "fftw3.h"
#include <complex.h>

/*
    Definition of the structure of the audio files to be analysed.
    typedef allows to create a synonym for the type of structure.
*/
typedef struct
{
  SNDFILE *audio;
  SF_INFO info;
  char *nome;
}dato_audio;

void info_dato(dato_audio *dato);
double **preleva_dati_finestra(double x[], int campioni_segnale, int freq_camp, int dim_finestra, int dim_passo, int *camp_fin, int *camp_passo, int *n_fin);
double *pre_enfasi(double *x, int campioni_segnale, double alpha);
double **fft(double **x, int n_finestre, int campioni_finestra, int frequenza_campionamento, int freq_max, int *n_freq, double *freq_passo);
double freq_mel(double freq);
double mel_freq(double mel);
double **filtri_mel(int n_filtri, int n_frequenze, int freq_max, double frequenza_passo);
double **mfcc(double **absX, int n_finestre, int n_frequenze, int campioni_finestra, int n_filtri, double **banco_filtri, int n_coeff_cepstrali, double *energia_log);
double *en_log(double **x, int n_finestre, int campioni_finestra);
double **regressione_lineare(double **x, int n_finestre, int n_coeff);
double **estrai_caratteristiche(int n_finestre, int n_coeff, double **mfcc, double **delta, double **doppio_delta);
double **analisi(double *segnale, int frequenza_campionamento, int campioni_segnale, int freq_max, int dim_finestra, int dim_passo, int n_coeff_cepstrali, int *n_fin);
double *apertura_wav(char nome_file[], int *freq_campionamento, int *camp_segnale);

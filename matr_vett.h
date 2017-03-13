/*
        Author:     Pietro Russo
        Date:       15/02/2011
        File:       matr_vett.h
        Purpose:    Functions to allocate memory arrays and matrices.
                    Useful functions for arrays.
*/

#include <stdio.h>

void errore(char testo[]);
double *vettore_d(int dimensione);
int *vettore_i(int dimensione);
double **matrice_d(int righe, int colonne);
int **matrice_i(int righe, int colonne);
double **trasposta_d(double **x, int righe, int colonne);
double **prodotto_matrici_d(double **A, int righe_A, int colonne_A, double **B, int colonne_B);
void salva_matrice_d(char nome_file[], double **dati, int righe, int colonne);
void carica_matrice_d(char nome_file[], double **dati, int righe, int colonne);
void visualizza_matrice_d(double **dati, int righe, int colonne);
void salva_vettore_i(char nome_file[], int *dati, int elementi);
void carica_vettore_i(char nome_file[], int *dati, int elementi);

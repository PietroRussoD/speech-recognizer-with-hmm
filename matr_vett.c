/*
        Author:     Pietro Russo
        Date:       15/02/2011
        File:       matr_vett.c
        Purpose:    Functions to allocate memory arrays and matrices.
                    Useful functions for arrays.
*/
#include "matr_vett.h"

/*   Function that prints the text string on stderr stream */
void errore(char testo[])
{
    fprintf(stderr,"%s\n",testo);
	fprintf(stderr,"Uscita dal programma...\n");
	exit(1);
}

/*  Function that allocates space in memory to store a double type vector */
double *vettore_d(int dimensione)
{
	double *v;
	int i;

    v = (double *)malloc(sizeof(double) * dimensione);
	if (v==NULL)
        errore("Allocazione del vettore double tramite funzione fallita");

	for(i=0; i<dimensione; i++)
        v[i]=0.0;

    return v;
}

/*  Function that allocates space in memory to store a int type vector */
int *vettore_i(int dimensione)
{
	int *v, i;

    v = (int *)malloc(sizeof(int) * dimensione);
	if (v==NULL)
        errore("Allocazione del vettore int tramite funzione fallita");

    for(i=0; i<dimensione; i++)
        v[i]=0;

	return v;
}

/*  Function that allocates space in memory to store a double type matrix */
double **matrice_d(int righe, int colonne)
{
	int i, j;
	double **m;

	m = (double **) malloc(sizeof(double*) * righe);
	if (m==NULL)
        errore("Allocazione delle righe della matrice double tramite funzione fallita");

	for(i=0; i<righe; i++)
	{
		m[i] = (double *) malloc(sizeof(double) * colonne);
		if (m[i]==NULL)
            errore("Allocazione delle colonne della matrice double tramite funzione fallita");

	}

    for(i=0; i<righe; i++)
        for(j=0; j<colonne; j++)
            m[i][j]=0.0;

	return m;
}

/*  Function that allocates space in memory to store a int type matrix */
int **matrice_i(int righe, int colonne)
{
	int i, j;
	int **m;

	m = (int **) malloc(sizeof(int*) * righe);
	if (m==NULL)
        errore("Allocazione delle righe della matrice int tramite funzione fallita");

	for(i=0; i<righe; i++)
	{
		m[i] = (int *) malloc(sizeof(int) * colonne);
		if (m[i]==NULL)
            errore("Allocazione delle colonne della matrice int tramite funzione fallita");

	}

    for(i=0; i<righe; i++)
        for(j=0; j<colonne; j++)
            m[i][j]=0.0;

	return m;
}

/*  Function that calculates the transpose of a matrix of type double */
double **trasposta_d(double **x, int righe, int colonne)
{
    double **y = matrice_d(colonne, righe);
    int n,m;

    for(n=0; n<colonne; n++)
        for(m=0; m<righe; m++)
            y[n][m]=x[m][n];

    return y;
}

/*  Function that calculates the product of matrices of type double */
double **prodotto_matrici_d(double **A, int righe_A, int colonne_A, double **B, int colonne_B)
{
    /*  If A is MxN and B id NxP then C is MxP */
    double temp;
    int n,m,i;
    double **C = matrice_d(righe_A, colonne_B);

    for(n=0; n<righe_A; n++)
        for(m=0; m<colonne_B; m++)
            {
                temp = 0.0;
                for(i=0; i<colonne_A; i++)
                    temp += A[n][i]*B[i][m];
                C[n][m]=temp;
            }

    return C;
}

/*  Function that saves on a binary file a double type matrix */
void salva_matrice_d(char nome_file[], double **dati, int righe, int colonne)
{
    FILE *fp;
    int i;
    fp=fopen(nome_file,"wb") ;
    if(fp==NULL)
         errore("Errore di apertura del file.");

    size_t dimensione = sizeof(double);
    size_t elementi = colonne ;
    size_t ritorno;

    for(i=0; i<righe; i++)
    {
        ritorno = fwrite(dati[i], dimensione, elementi, fp);
        if (ritorno != elementi)
            errore("Errore scrittura su file.");
    }

    fclose(fp);
}

/*  Function that load from a binary file a double type matrix */
void carica_matrice_d(char nome_file[], double **dati, int righe, int colonne)
{
    FILE *fp;
    int i;
    fp=fopen(nome_file,"rb") ;

    if(fp==NULL)
         errore("Errore di apertura del file.");

    size_t dimensione = sizeof(double);
    size_t elementi = colonne ;
    size_t ritorno;

    for(i=0; i<righe; i++)
    {
        ritorno = fread(dati[i], dimensione, elementi, fp);

        if (ritorno != colonne)
            errore("Errore scrittura su file.");
    }

    fclose(fp);
}

/*  Function that prints on the screen a double type matrix */
void visualizza_matrice_d(double **dati, int righe, int colonne)
{
    int n,m;
    for(n=0; n<righe; n++)
        for(m=0; m<colonne; m++)
            printf("%d,%d\t%f\n", n, m, dati[n][m]);
}

/*  Function that load from a binary file a double type vector */
void salva_vettore_i(char nome_file[], int *dati, int elementi)
{
    FILE *fp;
    int i;
    fp=fopen(nome_file,"wb") ;
    if(fp==NULL)
         errore("Errore di apertura del file.");

    size_t dimensione = sizeof(int);
    size_t ritorno;

    ritorno = fwrite(dati, dimensione, elementi, fp);
        if (ritorno != elementi)
            errore("Errore scrittura su file.");

    fclose(fp);
}

/*  Function that load from a binary file a double type vector */
void carica_vettore_i(char nome_file[], int *dati, int elementi)
{
    FILE *fp;
    int i;
    fp=fopen(nome_file,"rb") ;
    if(fp==NULL)
         errore("Errore di apertura del file.");

    size_t dimensione = sizeof(int);
    size_t ritorno;

    ritorno = fread(dati, dimensione, elementi, fp);
        if (ritorno != elementi)
            errore("Errore scrittura su file.");

    fclose(fp);
}

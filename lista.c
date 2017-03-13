/*
        Author:     Pietro Russo
        Date:       03/03/2011
        File:       lista.c
        Purpose:    Functions that are used to implement a list containing
                    voice signals used in the training set.
                    The signals have associated their HMM that will serve
                    in the recognition process to determine the word.
*/

#include "lista.h"

/*  Function that adds the wav file name, which determines signal characteristics. */
lista_parole *aggiungi(char nome_file[32], lista_parole *lista)
{
    lista_parole *nuovo = NULL;
    lista_parole *temporaneo = NULL;
    lista_parole *precedente = NULL;

    double *segnale=NULL, **caratteristiche=NULL;

    nuovo = (lista_parole *) malloc(sizeof(lista_parole));
    if(nuovo == NULL)
    {
        printf("\nImpossibile allocare la memoria.\n");
        exit(1);
    }

    /*  Copy the file name */
    strcpy(nuovo->nome_file_audio , nome_file);
    strcpy(nuovo->nome_file_car , "");
    strcpy(nuovo->nome_file_seq , "");
    strcpy(nuovo->nome_file_hmm , "");
    /*  Set to NULL the next element pointer */
    nuovo->prossimo = NULL;

    segnale=apertura_wav(nuovo->nome_file_audio, &nuovo->frequenza_campionamento, &nuovo->campioni_segnale);
    caratteristiche=analisi(segnale, nuovo->frequenza_campionamento, nuovo->campioni_segnale, FREQUENZA_MAX, DIM_FINESTRA, DIM_PASSO, N_COEFF_CEP, &nuovo->n_finestre);

    char *nome=strdup(nuovo->nome_file_audio);
    /*  Search for the last occurrence of the character '.' */
    char *pos=strrchr(nome, '.');

    int i=pos-nome;
    /*  After the '.', it changes the extension */
    nome[i+1]='c';
    nome[i+2]='a';
    nome[i+3]='r';

    strcpy(nuovo->nome_file_car, nome);

    salva_matrice_d(nuovo->nome_file_car, caratteristiche, nuovo->n_finestre, 3*N_COEFF_CEP);

    /*  If it is the first element, it adds on the top */
    if(lista == NULL)
    {
        lista = nuovo;
        nuovo->prossimo = NULL;
    }
    /*  If it is not the first */
    else
    {
        /*  Check if it goes before the first element
            strcmp (str1, SRT2) returns a value:
                < 0 if str1<str2
                = 0 if str1=str2
                > 0 if str1>str2 */
        if(strcmp(nuovo->nome_file_audio, lista->nome_file_audio) < 0)
        {
            nuovo->prossimo = lista;
            lista = nuovo;
        }
        /*  If it goes after the first */
        else
        {
            temporaneo = lista->prossimo;
            precedente = lista;
            /*  If it is the second */
            if(temporaneo == NULL)
            {
                precedente->prossimo = nuovo;
            }

            else
            {
                /*  If it goes before the end */
                while(temporaneo->prossimo != NULL)
                {
                    if(strcmp(nuovo->nome_file_audio, lista->nome_file_audio) < 0)
                    {
                        nuovo->prossimo = temporaneo;
                        if(nuovo->prossimo != precedente->prossimo)
                        {
                            printf("Errore");
                            exit(1);
                        }
                        precedente->prossimo = nuovo;
                        /*  Added item, and then exits the while */
                        break;
                    }
                    else
                    {
                        /*  Scan the list */
                        temporaneo = temporaneo->prossimo;
                        precedente = precedente->prossimo;
                    }
                }
                /*  If it is near the end */
                if(temporaneo->prossimo == NULL)
                {
                    /*  It is the penultimate */
                    if(strcmp(nuovo->nome_file_audio, lista->nome_file_audio) < 0)
                    {
                        nuovo->prossimo = temporaneo;
                        precedente->prossimo = nuovo;
                    }
                    /*  It is the last */
                    else
                    {
                        temporaneo->prossimo = nuovo;
                        nuovo->prossimo = NULL;
                    }
                }
            }
        }
    }
    return lista;
}

/*  Function that prints a list on the screen */
void mostra(lista_parole *lista)
{
    lista_parole *corrente;
    int i=1;

    corrente = lista;
    /*  Start from the top and scan the list */
    while(corrente != NULL)
    {
        printf("Parola numero:%3i\nNome File: %s\n", i++, corrente->nome_file_audio);
        printf("N.finestre: %d\nFreq. camp.: %d\n", corrente->n_finestre, corrente->frequenza_campionamento);
        printf("N.campioni: %d\n", corrente->campioni_segnale);
        printf("Nome Carat: %s\n", corrente->nome_file_car);
        printf("Nome Seq: %s\n", corrente->nome_file_seq);
        printf("Nome HMM: %s\n", corrente->nome_file_hmm);
        printf("\n\n");
        corrente = corrente->prossimo;
    }
}

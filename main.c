#include "analisi_caratteristiche.h"
#include "hmm.h"
#include "nrutil.h"

/*  Number of codevectors, allowed only power of 2 */
//#define NVQ 32
#define NVQ 64
//#define NVQ 128

/*
    If SET is 0 it compares with words not used in the training
    If SET is 1 it compares with words used in the training
    If SET is 2 it compares with words not used in the training contaminated with Gaussian noise
    If SET is 3 it compares with words used in the training contaminated with Gaussian noise */
#define SET 0

/*  If SC is 1 it uses the model with serial constraint with single and double transitions
    If SC is 2 it uses the model with serial constraint with single transitions */
#define SC 1

int main()
{
    int n,m,k, n_voci=0, i, j;
    lista_parole *lista=NULL;

    /*  Fill the list with the audio file names that will represent words of the training. */
    lista=aggiungi("00_01_giuseppe.s.wav", lista);
    lista=aggiungi("00_02_giuseppe.s.wav", lista);
    lista=aggiungi("00_03_giuseppe.s.wav", lista);
    lista=aggiungi("00_01_pietro.r.wav", lista);
    lista=aggiungi("00_02_pietro.r.wav", lista);
    lista=aggiungi("00_03_pietro.r.wav", lista);
    lista=aggiungi("01_01_giuseppe.s.wav", lista);
    lista=aggiungi("01_02_giuseppe.s.wav", lista);
    lista=aggiungi("01_03_giuseppe.s.wav", lista);
    lista=aggiungi("01_01_pietro.r.wav", lista);
    lista=aggiungi("01_02_pietro.r.wav", lista);
    lista=aggiungi("01_03_pietro.r.wav", lista);
    lista=aggiungi("02_01_giuseppe.s.wav", lista);
    lista=aggiungi("02_02_giuseppe.s.wav", lista);
    lista=aggiungi("02_03_giuseppe.s.wav", lista);
    lista=aggiungi("02_01_pietro.r.wav", lista);
    lista=aggiungi("02_02_pietro.r.wav", lista);
    lista=aggiungi("02_03_pietro.r.wav", lista);
    lista=aggiungi("03_01_giuseppe.s.wav", lista);
    lista=aggiungi("03_02_giuseppe.s.wav", lista);
    lista=aggiungi("03_03_giuseppe.s.wav", lista);
    lista=aggiungi("03_01_pietro.r.wav", lista);
    lista=aggiungi("03_02_pietro.r.wav", lista);
    lista=aggiungi("03_03_pietro.r.wav", lista);
    lista=aggiungi("04_01_giuseppe.s.wav", lista);
    lista=aggiungi("04_02_giuseppe.s.wav", lista);
    lista=aggiungi("04_03_giuseppe.s.wav", lista);
    lista=aggiungi("04_01_pietro.r.wav", lista);
    lista=aggiungi("04_02_pietro.r.wav", lista);
    lista=aggiungi("04_03_pietro.r.wav", lista);
    lista=aggiungi("05_01_giuseppe.s.wav", lista);
    lista=aggiungi("05_02_giuseppe.s.wav", lista);
    lista=aggiungi("05_03_giuseppe.s.wav", lista);
    lista=aggiungi("05_01_pietro.r.wav", lista);
    lista=aggiungi("05_02_pietro.r.wav", lista);
    lista=aggiungi("05_03_pietro.r.wav", lista);
    lista=aggiungi("06_01_giuseppe.s.wav", lista);
    lista=aggiungi("06_02_giuseppe.s.wav", lista);
    lista=aggiungi("06_03_giuseppe.s.wav", lista);
    lista=aggiungi("06_01_pietro.r.wav", lista);
    lista=aggiungi("06_02_pietro.r.wav", lista);
    lista=aggiungi("06_03_pietro.r.wav", lista);
    lista=aggiungi("07_01_giuseppe.s.wav", lista);
    lista=aggiungi("07_02_giuseppe.s.wav", lista);
    lista=aggiungi("07_03_giuseppe.s.wav", lista);
    lista=aggiungi("07_01_pietro.r.wav", lista);
    lista=aggiungi("07_02_pietro.r.wav", lista);
    lista=aggiungi("07_03_pietro.r.wav", lista);
    lista=aggiungi("08_01_giuseppe.s.wav", lista);
    lista=aggiungi("08_02_giuseppe.s.wav", lista);
    lista=aggiungi("08_03_giuseppe.s.wav", lista);
    lista=aggiungi("08_01_pietro.r.wav", lista);
    lista=aggiungi("08_02_pietro.r.wav", lista);
    lista=aggiungi("08_03_pietro.r.wav", lista);
    lista=aggiungi("09_01_giuseppe.s.wav", lista);
    lista=aggiungi("09_02_giuseppe.s.wav", lista);
    lista=aggiungi("09_03_giuseppe.s.wav", lista);
    lista=aggiungi("09_01_pietro.r.wav", lista);
    lista=aggiungi("09_02_pietro.r.wav", lista);
    lista=aggiungi("09_03_pietro.r.wav", lista);

    /*  M_vq represents the number of source vectors in the training set */
    int M_vq=1;

    lista_parole *corrente;
    corrente = lista;
    /*  Starts from the beginning and through the list */
    /*  Calculating the total number of windows on the list and the number of words included in the list */
    while(corrente != NULL)
    {
        M_vq += corrente->n_finestre;
        corrente = corrente->prossimo;
        n_voci++;
    }

    char **set=NULL;
    int lettere = 25;

    set = (char **) malloc(sizeof(char*) * n_voci+1);

	if (set==NULL)
        errore("Allocazione fallita");

	for(i=0; i<n_voci+1; i++)
	{
		set[i] = (char *) malloc(sizeof(char) * lettere);
		if (set[i]==NULL)
            errore("Allocazione fallita");
	}


    /*  K is the dimension of the vectors of the training set */
    int K = N_COEFF_CEP*3 + 1;

    /*  Declare the matrix that contains the set of the training set vectors */
    double **training = matrice_d(M_vq, K);

    corrente = lista;
    double **caratteristiche=NULL;
    m=1;
    /*  Fill the training matrix with the characteristics of the audio signals of the training set
        and keep track of the names in another variable */
    i=0;
    while(corrente != NULL)
    {
        caratteristiche=matrice_d(corrente->n_finestre, K-1);
        carica_matrice_d(corrente->nome_file_car, caratteristiche, corrente->n_finestre, K-1);

        for(n=0; n<corrente->n_finestre; n++, m++)
            for(k=0; k< K-1 ; k++)
            {
                training[m][k+1] = caratteristiche[n][k];
            }

        strcpy(set[i],corrente->nome_file_audio);
        i++;
        corrente = corrente->prossimo;
    }

    set[n_voci+1]="nessuno";

    /*  Declare the matrix of codebook */
    double **codebook=NULL;
    /*  In the vector quantization:
            - M is the total number of vectors of the training set
            - N is the number of codevectors
            - K is the vector dimension of training set */
    int N_vq = NVQ + 1;

    /*  Determine the codebook */
    codebook=lbg(training, M_vq, K, N_vq);

    FILE *codebook_file=NULL;
    codebook_file=fopen("codebook.txt","w");
    for(n=1; n<N_vq; n++)
    {
        fprintf(codebook_file,"%d\n", n);
        for(k=1; k<K; k++)
            fprintf(codebook_file,"%f\t", codebook[n][k]);
        fprintf(codebook_file,"\n");
    }

    corrente = lista;
    char *nome=NULL;
    int *sequenza=NULL;
    char *pos=NULL;

    FILE *tutte_caratteristiche=NULL;
    tutte_caratteristiche = fopen("tutte_caratteristiche.txt","w");

    FILE *tutte_sequenze=NULL;
    tutte_sequenze = fopen("tutte_sequenze.txt","w");

    /*  After the codebook, it determines the vector quantization for every vector of the training set */
    while(corrente != NULL)
    {
        caratteristiche=matrice_d(corrente->n_finestre, K-1);
        carica_matrice_d(corrente->nome_file_car, caratteristiche, corrente->n_finestre, K-1);

        /*  The quantized vector is the vector of sequences that will be used in the Hidden Markov Models */
        sequenza = vq(codebook, caratteristiche, corrente->n_finestre, K, N_vq);
        nome=strdup(corrente->nome_file_audio);
        /*  Find the last character '.' */
        pos=strrchr(nome, '.' );
        i=pos-nome;
        /*  After '.' changes the extension */
        nome[i+1]='s';
        nome[i+2]='e';
        nome[i+3]='q';

        fprintf(tutte_sequenze,"%s\n", corrente->nome_file_audio);
        for(n=1; n<=corrente->n_finestre; n++)
        {
            fprintf(tutte_sequenze, "%d\t", sequenza[n]);
        }
        fprintf(tutte_sequenze, "\n\n");


        strcpy(corrente->nome_file_seq, nome);
        salva_vettore_i(corrente->nome_file_seq, sequenza, corrente->n_finestre +1);

        corrente = corrente->prossimo;
    }

    /*  Set the seed of the random. */
    unsigned int seed=333;

    corrente = lista;

    HMM hmm;
    /*  In the Hidden Markov Models:
            - M is the number of the observed symbols
            - N is the number of states */
    int M=N_vq-1;
    int N=5;

    double	logprobinit, logprobfinal;
    double 	**alpha;
	double	**beta;
	double	**gamma;

    int *O=NULL;
    int T, niter, T_max=0;

    while(corrente != NULL)
    {
        nome=strdup(corrente->nome_file_audio);
        pos=strrchr(nome, '.' );
        i=pos-nome;
        nome[i+1]='h';
        nome[i+2]='m';
        nome[i+3]='m';

        strcpy(corrente->nome_file_hmm, nome);

        O=vettore_i(corrente->n_finestre +1);

        carica_vettore_i(corrente->nome_file_seq, O, corrente->n_finestre +1);

        T=corrente->n_finestre;

        if(T_max<T)
            T_max=T;

        #if SC==1
            InitHMM_SC1(&hmm, N, M, seed);
        #elif SC==2
            InitHMM_SC2(&hmm, N, M, seed);
        #endif

        alpha = dmatrix(1, T, 1, hmm.N);
        beta = dmatrix(1, T, 1, hmm.N);
        gamma = dmatrix(1, T, 1, hmm.N);

        BaumWelch_C(&hmm, T, O, alpha, beta, gamma, &niter, &logprobinit, &logprobfinal);

        Salva_HMM(corrente->nome_file_hmm, &hmm);

        corrente = corrente->prossimo;
    }

    double 	proba[n_voci];
    double somma_proba=0;

    int	*q;
    double **delta;
    int	**psi;

    int frequenza_campionamento, campioni_segnale, n_finestre;

    double *segnale=NULL;
    char **analizza=NULL;

    #if (SET==0 || SET ==2)
        int parole = 40;
    #elif (SET==1 || SET ==3)
        int parole = 60;
    #endif

    analizza = (char **) malloc(sizeof(char*) * parole);

	if (analizza==NULL)
        errore("Allocazione fallita");

	for(i=0; i<parole; i++)
	{
		analizza[i] = (char *) malloc(sizeof(char) * lettere);
		if (analizza[i]==NULL)
            errore("Allocazione fallita");
	}

    #if SET==0
        analizza[0]="00_04_giuseppe.s.wav";
        analizza[1]="00_05_giuseppe.s.wav";
        analizza[2]="00_04_pietro.r.wav";
        analizza[3]="00_05_pietro.r.wav";
        analizza[4]="01_04_giuseppe.s.wav";
        analizza[5]="01_05_giuseppe.s.wav";
        analizza[6]="01_04_pietro.r.wav";
        analizza[7]="01_05_pietro.r.wav";
        analizza[8]="02_04_giuseppe.s.wav";
        analizza[9]="02_05_giuseppe.s.wav";
        analizza[10]="02_04_pietro.r.wav";
        analizza[11]="02_05_pietro.r.wav";
        analizza[12]="03_04_giuseppe.s.wav";
        analizza[13]="03_05_giuseppe.s.wav";
        analizza[14]="03_04_pietro.r.wav";
        analizza[15]="03_05_pietro.r.wav";
        analizza[16]="04_04_giuseppe.s.wav";
        analizza[17]="04_05_giuseppe.s.wav";
        analizza[18]="04_04_pietro.r.wav";
        analizza[19]="04_05_pietro.r.wav";
        analizza[20]="05_04_giuseppe.s.wav";
        analizza[21]="05_05_giuseppe.s.wav";
        analizza[22]="05_04_pietro.r.wav";
        analizza[23]="05_05_pietro.r.wav";
        analizza[24]="06_04_giuseppe.s.wav";
        analizza[25]="06_05_giuseppe.s.wav";
        analizza[26]="06_04_pietro.r.wav";
        analizza[27]="06_05_pietro.r.wav";
        analizza[28]="07_04_giuseppe.s.wav";
        analizza[29]="07_05_giuseppe.s.wav";
        analizza[30]="07_04_pietro.r.wav";
        analizza[31]="07_05_pietro.r.wav";
        analizza[32]="08_04_giuseppe.s.wav";
        analizza[33]="08_05_giuseppe.s.wav";
        analizza[34]="08_04_pietro.r.wav";
        analizza[35]="08_05_pietro.r.wav";
        analizza[36]="09_04_giuseppe.s.wav";
        analizza[37]="09_05_giuseppe.s.wav";
        analizza[38]="09_04_pietro.r.wav";
        analizza[39]="09_05_pietro.r.wav";
	#elif SET==1
        analizza[0]="00_01_giuseppe.s.wav";
        analizza[1]="00_02_giuseppe.s.wav";
        analizza[2]="00_03_giuseppe.s.wav";
        analizza[3]="00_01_pietro.r.wav";
        analizza[4]="00_02_pietro.r.wav";
        analizza[5]="00_03_pietro.r.wav";
        analizza[6]="01_01_giuseppe.s.wav";
        analizza[7]="01_02_giuseppe.s.wav";
        analizza[8]="01_03_giuseppe.s.wav";
        analizza[9]="01_01_pietro.r.wav";
        analizza[10]="01_02_pietro.r.wav";
        analizza[11]="01_03_pietro.r.wav";
        analizza[12]="02_01_giuseppe.s.wav";
        analizza[13]="02_02_giuseppe.s.wav";
        analizza[14]="02_03_giuseppe.s.wav";
        analizza[15]="02_01_pietro.r.wav";
        analizza[16]="02_02_pietro.r.wav";
        analizza[17]="02_03_pietro.r.wav";
        analizza[18]="03_01_giuseppe.s.wav";
        analizza[19]="03_02_giuseppe.s.wav";
        analizza[20]="03_03_giuseppe.s.wav";
        analizza[21]="03_01_pietro.r.wav";
        analizza[22]="03_02_pietro.r.wav";
        analizza[23]="03_03_pietro.r.wav";
        analizza[24]="04_01_giuseppe.s.wav";
        analizza[25]="04_02_giuseppe.s.wav";
        analizza[26]="04_03_giuseppe.s.wav";
        analizza[27]= "04_01_pietro.r.wav";
        analizza[28]="04_02_pietro.r.wav";
        analizza[29]="04_03_pietro.r.wav";
        analizza[30]="05_01_giuseppe.s.wav";
        analizza[31]="05_02_giuseppe.s.wav";
        analizza[32]="05_03_giuseppe.s.wav";
        analizza[33]="05_01_pietro.r.wav";
        analizza[34]="05_02_pietro.r.wav";
        analizza[35]="05_03_pietro.r.wav";
        analizza[36]="06_01_giuseppe.s.wav";
        analizza[37]="06_02_giuseppe.s.wav";
        analizza[38]="06_03_giuseppe.s.wav";
        analizza[39]="06_01_pietro.r.wav";
        analizza[40]="06_02_pietro.r.wav";
        analizza[41]="06_03_pietro.r.wav";
        analizza[42]="07_01_giuseppe.s.wav";
        analizza[43]="07_02_giuseppe.s.wav";
        analizza[44]="07_03_giuseppe.s.wav";
        analizza[45]="07_01_pietro.r.wav";
        analizza[46]="07_02_pietro.r.wav";
        analizza[47]="07_03_pietro.r.wav";
        analizza[48]="08_01_giuseppe.s.wav";
        analizza[49]="08_02_giuseppe.s.wav";
        analizza[50]="08_03_giuseppe.s.wav";
        analizza[51]="08_01_pietro.r.wav";
        analizza[52]="08_02_pietro.r.wav";
        analizza[53]="08_03_pietro.r.wav";
        analizza[54]="09_01_giuseppe.s.wav";
        analizza[55]="09_02_giuseppe.s.wav";
        analizza[56]="09_03_giuseppe.s.wav";
        analizza[57]="09_01_pietro.r.wav";
        analizza[58]="09_02_pietro.r.wav";
        analizza[59]="09_03_pietro.r.wav";
    #elif SET==2
        analizza[0]="00_04_giuseppe.s_n.wav";
        analizza[1]="00_05_giuseppe.s_n.wav";
        analizza[2]="00_04_pietro.r_n.wav";
        analizza[3]="00_05_pietro.r_n.wav";
        analizza[4]="01_04_giuseppe.s_n.wav";
        analizza[5]="01_05_giuseppe.s_n.wav";
        analizza[6]="01_04_pietro.r_n.wav";
        analizza[7]="01_05_pietro.r_n.wav";
        analizza[8]="02_04_giuseppe.s_n.wav";
        analizza[9]="02_05_giuseppe.s_n.wav";
        analizza[10]="02_04_pietro.r_n.wav";
        analizza[11]="02_05_pietro.r_n.wav";
        analizza[12]="03_04_giuseppe.s_n.wav";
        analizza[13]="03_05_giuseppe.s_n.wav";
        analizza[14]="03_04_pietro.r_n.wav";
        analizza[15]="03_05_pietro.r_n.wav";
        analizza[16]="04_04_giuseppe.s_n.wav";
        analizza[17]="04_05_giuseppe.s_n.wav";
        analizza[18]="04_04_pietro.r_n.wav";
        analizza[19]="04_05_pietro.r_n.wav";
        analizza[20]="05_04_giuseppe.s_n.wav";
        analizza[21]="05_05_giuseppe.s_n.wav";
        analizza[22]="05_04_pietro.r_n.wav";
        analizza[23]="05_05_pietro.r_n.wav";
        analizza[24]="06_04_giuseppe.s_n.wav";
        analizza[25]="06_05_giuseppe.s_n.wav";
        analizza[26]="06_04_pietro.r_n.wav";
        analizza[27]="06_05_pietro.r_n.wav";
        analizza[28]="07_04_giuseppe.s_n.wav";
        analizza[29]="07_05_giuseppe.s_n.wav";
        analizza[30]="07_04_pietro.r_n.wav";
        analizza[31]="07_05_pietro.r_n.wav";
        analizza[32]="08_04_giuseppe.s_n.wav";
        analizza[33]="08_05_giuseppe.s_n.wav";
        analizza[34]="08_04_pietro.r_n.wav";
        analizza[35]="08_05_pietro.r_n.wav";
        analizza[36]="09_04_giuseppe.s_n.wav";
        analizza[37]="09_05_giuseppe.s_n.wav";
        analizza[38]="09_04_pietro.r_n.wav";
        analizza[39]="09_05_pietro.r_n.wav";
    #elif SET==3
         analizza[0]="00_01_giuseppe.s_n.wav";
        analizza[1]="00_02_giuseppe.s_n.wav";
        analizza[2]="00_03_giuseppe.s_n.wav";
        analizza[3]="00_01_pietro.r_n.wav";
        analizza[4]="00_02_pietro.r_n.wav";
        analizza[5]="00_03_pietro.r_n.wav";
        analizza[6]="01_01_giuseppe.s_n.wav";
        analizza[7]="01_02_giuseppe.s_n.wav";
        analizza[8]="01_03_giuseppe.s_n.wav";
        analizza[9]="01_01_pietro.r_n.wav";
        analizza[10]="01_02_pietro.r_n.wav";
        analizza[11]="01_03_pietro.r_n.wav";
        analizza[12]="02_01_giuseppe.s_n.wav";
        analizza[13]="02_02_giuseppe.s_n.wav";
        analizza[14]="02_03_giuseppe.s_n.wav";
        analizza[15]="02_01_pietro.r_n.wav";
        analizza[16]="02_02_pietro.r_n.wav";
        analizza[17]="02_03_pietro.r_n.wav";
        analizza[18]="03_01_giuseppe.s_n.wav";
        analizza[19]="03_02_giuseppe.s_n.wav";
        analizza[20]="03_03_giuseppe.s_n.wav";
        analizza[21]="03_01_pietro.r_n.wav";
        analizza[22]="03_02_pietro.r_n.wav";
        analizza[23]="03_03_pietro.r_n.wav";
        analizza[24]="04_01_giuseppe.s_n.wav";
        analizza[25]="04_02_giuseppe.s_n.wav";
        analizza[26]="04_03_giuseppe.s_n.wav";
        analizza[27]= "04_01_pietro.r_n.wav";
        analizza[28]="04_02_pietro.r_n.wav";
        analizza[29]="04_03_pietro.r_n.wav";
        analizza[30]="05_01_giuseppe.s_n.wav";
        analizza[31]="05_02_giuseppe.s_n.wav";
        analizza[32]="05_03_giuseppe.s_n.wav";
        analizza[33]="05_01_pietro.r_n.wav";
        analizza[34]="05_02_pietro.r_n.wav";
        analizza[35]="05_03_pietro.r_n.wav";
        analizza[36]="06_01_giuseppe.s_n.wav";
        analizza[37]="06_02_giuseppe.s_n.wav";
        analizza[38]="06_03_giuseppe.s_n.wav";
        analizza[39]="06_01_pietro.r_n.wav";
        analizza[40]="06_02_pietro.r_n.wav";
        analizza[41]="06_03_pietro.r_n.wav";
        analizza[42]="07_01_giuseppe.s_n.wav";
        analizza[43]="07_02_giuseppe.s_n.wav";
        analizza[44]="07_03_giuseppe.s_n.wav";
        analizza[45]="07_01_pietro.r_n.wav";
        analizza[46]="07_02_pietro.r_n.wav";
        analizza[47]="07_03_pietro.r_n.wav";
        analizza[48]="08_01_giuseppe.s_n.wav";
        analizza[49]="08_02_giuseppe.s_n.wav";
        analizza[50]="08_03_giuseppe.s_n.wav";
        analizza[51]="08_01_pietro.r_n.wav";
        analizza[52]="08_02_pietro.r_n.wav";
        analizza[53]="08_03_pietro.r_n.wav";
        analizza[54]="09_01_giuseppe.s_n.wav";
        analizza[55]="09_02_giuseppe.s_n.wav";
        analizza[56]="09_03_giuseppe.s_n.wav";
        analizza[57]="09_01_pietro.r_n.wav";
        analizza[58]="09_02_pietro.r_n.wav";
        analizza[59]="09_03_pietro.r_n.wav";
    #endif

    FILE *risultati=NULL;
    risultati=fopen("risultati.txt","w");

    for(j=0; j<parole; j++)
    {
        corrente = lista;

        frequenza_campionamento=campioni_segnale=n_finestre=0;

        segnale=apertura_wav(analizza[j], &frequenza_campionamento, &campioni_segnale);
        caratteristiche=analisi(segnale, frequenza_campionamento, campioni_segnale, FREQUENZA_MAX, DIM_FINESTRA, DIM_PASSO, N_COEFF_CEP, &n_finestre);

        sequenza = vq(codebook, caratteristiche, n_finestre, K, N_vq);

        T = n_finestre;

        q = ivector(1,T);
        delta = dmatrix(1, T, 1, N);
        psi = imatrix(1, T, 1, N);

        i=0;

        double probi[n_voci];

        while(corrente != NULL)
        {
            proba[i]=0;
            Carica_HMM(corrente->nome_file_hmm, &hmm, N, M);

            ViterbiLog_C(&hmm, T, sequenza, delta, psi, q, &proba[i]);

            probi[i]=proba[i];
            i++;
            corrente = corrente->prossimo;
        }

        double alt_proba;
        double basso_proba;
        int val_alto, val_basso;

        alt_proba=basso_proba=proba[0];
        val_alto=val_basso=0;

        int posizione[n_voci];

        for(i=0; i<n_voci; i++)
        {
            posizione[i]=n_voci+1;
            if(proba[i]>alt_proba)
            {
                alt_proba = proba[i];
                val_alto = i;
            }
            if(proba[i]<basso_proba)
            {
                basso_proba=proba[i];
                val_basso = i;
            }
        }

        fprintf(risultati,"%s\tlo associo a\t%s\n", analizza[j], set[val_alto]);

    }
    return 0;
}

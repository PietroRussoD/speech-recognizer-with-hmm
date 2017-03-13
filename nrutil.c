/**     File:       nrutil.c
**      Purpose:    Memory allocation routines borrowed from the
**		            book "Numerical Recipes" by Press, Flannery, Teukolsky,
**		            and Vetterling.
**                  state sequence and probablity of observing a sequence
**                  given the model.
**      Organization: University of Maryland
**
**      $Id: nrutil.c,v 1.2 1998/02/19 16:31:35 kanungo Exp kanungo $
*/

#include <malloc.h>
#include <stdio.h>
static char rcsid[] = "$Id: nrutil.c,v 1.2 1998/02/19 16:31:35 kanungo Exp kanungo $";

void nrerror(char error_text[])
{
	void exit();

	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}

/*  Prepare a memory location that will contain a vector of type int */
int *ivector(int nl, int nh)
{
	int *v;

	v=(int *)calloc((unsigned) (nh-nl+1),sizeof(int));
	if (!v) nrerror("allocation failure in ivector()");
	return v-nl;
}

/*  Prepare a memory location that will contain a vector of type double */
double *dvector(int nl, int nh)
{
	double *v;

	v=(double *)calloc((unsigned) (nh-nl+1),sizeof(double));
	if (!v) nrerror("allocation failure in dvector()");
	return v-nl;
}

/*  Prepare a memory location that will contain a matrix of type double */
double **dmatrix(int nrl, int nrh, int ncl, int nch)
{
	int i;
	double **m;

	m=(double **) calloc((unsigned) (nrh-nrl+1),sizeof(double*));
	if (!m) nrerror("allocation failure 1 in dmatrix()");
	m -= nrl;

	for(i=nrl;i<=nrh;i++)
	{
		m[i]=(double *) calloc((unsigned) (nch-ncl+1),sizeof(double));
		if (!m[i]) nrerror("allocation failure 2 in dmatrix()");
		m[i] -= ncl;
	}
	return m;
}

/*  Prepare a memory location that will contain a matrix of type int */
int **imatrix(int nrl, int nrh, int ncl, int nch)
{
	int i,**m;

	m=(int **)calloc((unsigned) (nrh-nrl+1),sizeof(int*));
	if (!m) nrerror("allocation failure 1 in imatrix()");
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i]=(int *)calloc((unsigned) (nch-ncl+1),sizeof(int));
		if (!m[i]) nrerror("allocation failure 2 in imatrix()");
		m[i] -= ncl;
	}
	return m;
}

/*  free() Free a dynamically allocated memory block */
void free_ivector(int *v, int nl, int nh)
{
	free((char*) (v+nl));
}

void free_dvector(double *v, int nl, int nh)
{
	free((char*) (v+nl));
}

void free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch)
{
	int i;

	for(i=nrh; i>=nrl; i--)
        free((char*) (m[i]+ncl));
	free((char*) (m+nrl));
}

void free_imatrix(int **m, int nrl, int nrh, int ncl, int nch)
{
	int i;

	for(i=nrh; i>=nrl; i--)
        free((char*) (m[i]+ncl));
	free((char*) (m+nrl));
}

/*
**      File:       nrutil.h
**      Purpose:    Memory allocation routines borrowed from the
**                  book "Numerical Recipes" by Press, Flannery, Teukolsky,
**                  and Vetterling.
**                  state sequence and probablity of observing a sequence
**                  given the model.
**      Organization: University of Maryland
**
**      $Id: nrutil.h,v 1.2 1998/02/19 16:32:42 kanungo Exp kanungo $
*/

double *dvector();
double **dmatrix();
int *ivector();
int **imatrix();

void free_dvector();
void free_ivector();
void free_dmatrix();
void free_imatrix();

void nrerror();

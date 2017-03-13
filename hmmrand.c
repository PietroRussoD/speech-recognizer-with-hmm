/*
**	Author:     Tapas Kanungo, kanungo@cfar.umd.edu
**	File:	    hmmrand.c
**	Date:	    4 May 1999
**	Purpose:    To separate out the random number generator
** 		        functions so that the rest of the code can be
**		        platform independent.
*/

/*
        Change
        Author:     Pietro Russo
        Date:       21/02/2011
        File:       hmmrand.c
        Purpose:    Changing the function to generate the seed for the random numbers.
*/

/*  This part of the functions belonging to the random number generator has been divided from the rest of the code,
    in such a way that is independent of the platform */

/*  In stdlib.h there are function srand() e rand() */
#include <stdlib.h>
/*  sys/types.h it is used from the standard POSIX to define all the derived type of data. */
/*   */
//#include <sys/types.h>
#include <time.h>
#include <unistd.h>


/*  Generate a random seed for the random number generator */
unsigned int  hmmgetseed(void)
{
    /*  Changing the code to make the creation of multiplatform random number */
    return ((unsigned)time(NULL));
	/*  Line of code that was there before, linked to a linux system function. */
	//return ((int) getpid());
}

/*  Sets the seed for the random number generator to a specific value */
void hmmsetseed(unsigned int seed)
{   /*  Serve to generate a different sequence of random numbers at each execution,
        if seed is always different. */
	srand(seed);
}

/*  Returns a double value, which is a pseudo random number in the range [0,1). */
double hmmgetrand(void)
{
	return (double) rand()/RAND_MAX;
    /*  RAND_MAX is defined in stdlib.h file and it is equal to 0x7FFF.
        It is the maximum number that can return the rand function. The minimum is zero. */
}

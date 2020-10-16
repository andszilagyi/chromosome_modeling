#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "rand.h"

#define AA 471
#define B 1586
#define CC 6988
#define DD 9689
#define M 16383
#define RIMAX 2147483648.0        /* = 2^31 */
#define RandomInteger (++nd, ra[nd & M] = ra[(nd-AA) & M] ^ ra[(nd-B) & M] ^ ra[(nd-CC) & M] ^ ra[(nd-DD) & M])

static long ra[M+1], nd;

void seed(long seed)
{
 int  i;

 if(seed<=0) { puts("SEED error."); exit(1); }
 ra[0]= (long) fmod(16807.0*(double)seed, 2147483647.0);
 for(i=1; i<=M; i++)
 {
  ra[i] = (long)fmod( 16807.0 * (double) ra[i-1], 2147483647.0);
 }
}

long randl(long num)      /* random number between 0 and num-1 */
{
 return(RandomInteger % num);
}

double randd(void)
{
 return((double) RandomInteger / RIMAX);
}


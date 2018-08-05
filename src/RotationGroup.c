#include <math.h>
#include <gsl/gsl_sf_coupling.h>

double threej_(int *j1, int *j2, int *j3, int *m1, int *m2, int *m3)
{
  return gsl_sf_coupling_3j(*j1, *j2, *j3, *m1, *m2, *m3);
}

double cg_(int *j1, int *j2, int *j3, int *m1, int *m2, int *m3)
{
  int M = -*m3;
  return pow(-1.0, (*j1-*j2+*m3) * 0.5) * sqrt(*j3+1) *
    threej_(j1, j2, j3, m1,m2,&M);
}

double sixj_(int *j1, int *j2, int *j3, int *l1, int *l2, int *l3)
{
  return gsl_sf_coupling_6j(*j1, *j2, *j3, *l1, *l2, *l3);
}

double ninej_(int *j1, int *j2, int *j12, int *j3, int *j4, int *j34, int *j13, int *j24, int *j)
{
  return gsl_sf_coupling_9j(*j1, *j2, *j12, *j3, *j4, *j34, *j13, *j24, *j);
}

#include <stdio.h>
#include <gsl/gsl_sf_coupling.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_laguerre.h>
#include <gsl/gsl_integration.h>

double coupling_3j(int j1, int j2, int j3, int m1, int m2, int m3)
{ return gsl_sf_coupling_3j(j1, j2, j3, m1, m2, m3); }

double coupling_6j(int j1, int j2, int j3, int j4, int j5, int j6)
{ return gsl_sf_coupling_6j(j1, j2, j3, j4, j5, j6); }

double coupling_9j(int j1, int j2, int j3, int j4, int j5, int j6, int j7, int j8, int j9)
{ return gsl_sf_coupling_9j(j1, j2, j3, j4, j5, j6, j7, j8, j9); }

double factorial(unsigned int n)
{
  if(n > 170) {
    printf("Error in factorial, n is too large: n = %d\n",n);
    return 0.0;
  }
  return gsl_sf_fact(n);
}

double double_factorial(unsigned int n)
{
  if(n > 297) {
    printf("Error in double_factorial, n is too large: n = %d\n",n);
    return 0.0;
  }
  return gsl_sf_doublefact(n);
}

double gamma_function(double x)
{ return gsl_sf_gamma(x); }

double spherical_bessel(int l, double x)
{ return gsl_sf_bessel_jl(l, x); }

double legendre_polynomial(int l, double x)
{ return gsl_sf_legendre_Pl(l, x); }

double laguerre(int n, double a, double x)
{ return gsl_sf_laguerre_n(n, a, x); }

// Gauss-Legendre quadrature
gsl_integration_glfixed_table * gauss_legendre_allocate(unsigned int n)
{ return gsl_integration_glfixed_table_alloc(n); }

int gauss_legendre_ith_point_weight(double a, double b, unsigned int i,
    double * xi, double * wi, gsl_integration_glfixed_table * t)
{ return gsl_integration_glfixed_point(a, b, i, xi, wi, t); }

void gauss_legendre_release(gsl_integration_glfixed_table * t)
{ gsl_integration_glfixed_table_free(t); }
// Gauss-Legendre quadrature

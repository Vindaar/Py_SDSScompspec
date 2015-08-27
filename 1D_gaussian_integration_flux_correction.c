// This function is the integrand of the function which is to be calculated
// in the flux correction, if perform_2d is used. That is, both integral are
// calculated numerically.

#include <math.h>
# define M_PIl          3.141592653589793238462643383279502884L /* pi */

double f(int n, double args[n]){
  double r_prime_2 = 0;
  double value = 0;
  double r     = args[0];
  double delta_y = args[1];
  double sigma   = args[2];
  double a       = 1.0;

  value = 1/(2*M_PIl)*(1 - exp(-(a*a - delta_y*delta_y + 2*delta_y*sin(r)*(delta_y*sin(r) + sqrt(a*a + delta_y*delta_y*(sin(r)*sin(r) - 1))))/(2*sigma*sigma)));

  return value;
}

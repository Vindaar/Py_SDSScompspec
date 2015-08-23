// This function is the integrand of the function which is to be calculated
// in the flux correction, if perform_2d is used. That is, both integral are
// calculated numerically.

#include <math.h>
# define M_PIl          3.141592653589793238462643383279502884L /* pi */

double f(int n, double args[n]){
  double r_prime_2 = 0;
  double value = 0;
  double r     = args[0];
  double theta = args[1];
  double delta_y = args[2];
  double sigma   = args[3];

  r_prime_2 = r*r + delta_y*delta_y - 2*r*delta_y*cos(3.0/2.0*M_PIl - theta);

  value =  1/(2*M_PIl*sigma*sigma)*exp(-r_prime_2/(2*sigma*sigma))*r;
  
  return value;
}

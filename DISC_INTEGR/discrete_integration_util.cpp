#include "discrete_integration_util.hpp"

/**
 * @brief Performs a trapezoidal integration step:
 */
DISC_INTEGR_API void trapz_without_dt(
  double val_left,
  double val_right, 
  double *result
)
{
	*result += (val_left + val_right)/2.;
}

/**
 * @brief Integrate a function with a simple trapezoidal method.
 */
DISC_INTEGR_API double trapzIntegral(
  const boost::function<double (double)>& f,
  double start_x,
  double end_x,
  double dx)
{
  double res = 0.0;
  while((start_x+=dx) <= end_x){
    res += (f(start_x)+f(start_x-dx))/2 * dx;
  }
  return res;
}

// Following function is from http://www.bnikolic.co.uk/nqm/1dinteg/gausslobatto.html
/** \brief Perform a single step of the Gauss-Lobatto integration
 
    \param f Function to integrate
    \param a Lower integration limit
    \param b Upper integration limit

    \param fa Value of function at the lower limit (used to save an
    evaluation when refinement is used)
    
    \param fa Value of function at the upper limit (used to save an
    evaluation when refinement is used)
    
    \param neval Number of evaluations made so far
    
    \param maxeval Maximum number of evalutions which should not be
    exceeded

    \param acc Required accuracy expressed in units of
    std::numeric_limits<double>::epsilon(). This allows less-than
    comparison by using addition and equality.
 */
DISC_INTEGR_API double GaussLobattoIntStep(const boost::function<double (double)>& f, 
			   double a, double b,
			   double fa, double fb,
			   size_t &neval,
			   size_t maxeval,
			   double acc,
         double frequency)
{

  // Constants used in the algorithm
  const double alpha = std::sqrt(2.0/3.0); 
  const double beta  = 1.0/std::sqrt(5.0);

  if (neval >= maxeval)
  {
    throw std::runtime_error("Maximum number of evaluations reached in GaussLobatto");
  }

  // Here the abcissa points and function values for both the 4-point
  // and the 7-point rule are calculated (the points at the end of
  // interval come from the function call, i.e., fa and fb. Also note
  // the 7-point rule re-uses all the points of the 4-point rule.)
  const double h=(b-a)/2; 
  const double m=(a+b)/2;
        
  const double mll=m-alpha*h; 
  const double ml =m-beta*h; 
  const double mr =m+beta*h; 
  const double mrr=m+alpha*h;
  
  const double fmll= f(mll * frequency);
  const double fml = f(ml  * frequency);
  const double fm  = f(m   * frequency);
  const double fmr = f(mr  * frequency);
  const double fmrr= f(mrr * frequency);
  neval+=5;
        
  // Both the 4-point and 7-point rule integrals are evaluted
  const double integral2=(h/6)*(fa+fb+5*(fml+fmr));
  const double integral1=(h/1470)*(77*(fa+fb)
				   +432*(fmll+fmrr)+625*(fml+fmr)+672*fm);
  
  // The difference betwen the 4-point and 7-point integrals is the
  // estimate of the accuracy
  const double estacc=(integral1-integral2);

  // The volatile keyword should prevent the floating point
  // destination value from being stored in extended precision
  // registers which actually have a very different
  // std::numeric_limits<double>::epsilon(). 
  volatile double dist = acc + estacc;

  if(dist==acc || mll<=a || b<=mrr) 
  {
    if (not (m>a && b>m))
    {
      throw std::runtime_error("Integration reached an interval with no more machine numbers");
    }
    return integral1;
  }
  else {
    return  GaussLobattoIntStep(f, a, mll, fa, fmll, neval, maxeval, acc,frequency)  
      + GaussLobattoIntStep(f, mll, ml, fmll, fml, neval, maxeval, acc,frequency)
      + GaussLobattoIntStep(f, ml, m, fml, fm, neval, maxeval, acc,frequency)
      + GaussLobattoIntStep(f, m, mr, fm, fmr, neval, maxeval, acc,frequency)
      + GaussLobattoIntStep(f, mr, mrr, fmr, fmrr, neval, maxeval, acc,frequency)
      + GaussLobattoIntStep(f, mrr, b, fmrr, fb, neval, maxeval, acc,frequency);
         
  }
}

// Following function is from http://www.bnikolic.co.uk/nqm/1dinteg/gausslobatto.html
/** \brief Compute the Gauss-Lobatto integral

    \param f The function to be integrated

    \param a The lower integration limit

    \param b The upper integration limit

    \param abstol Absolute tolerance -- integration stops when the
    error estimate is smaller than this

    \param maxeval Maxium of evaluations to make. If this number of
    evalution is made without reaching the requied accuracy, an
    exception of type std::runtime_error is thrown.
 */
DISC_INTEGR_API double GaussLobattoInt(const boost::function<double (double)>& f, 
		       double a, double b,
           double frequency,
		       double abstol, 
		       size_t maxeval) 
{
  const double tol_epsunit = abstol/std::numeric_limits<double>::epsilon();
  size_t neval=0;
  return GaussLobattoIntStep(f, a, b,
			     f(frequency*a), f(frequency*b),
			     neval,
			     maxeval,
			     tol_epsunit,
           frequency);
}
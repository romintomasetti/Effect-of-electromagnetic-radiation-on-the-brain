#ifndef DISCRETE_INTEGRATION_UTIL_HPP
#define DISCRETE_INTEGRATION_UTIL_HPP

#if defined(WIN32)
#ifdef discr_integr_EXPORTS
#define DISC_INTEGR_API __declspec(dllexport)
#else
#define DISC_INTEGR_API __declspec(dllimport)
#endif
#else
#define DISC_INTEGR_API
#endif

#include <cmath>
#include <boost/function.hpp>

#include <iostream>
#include <boost/format.hpp>

/**
 * @brief Integrate a function with a simple trapezoidal method.
 */
DISC_INTEGR_API double trapzIntegral(
  const boost::function<double (double)>& f,
  double start_x,
  double end_x,
  double dx);


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
			   double acc);

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
		       double abstol, 
		       size_t maxeval);

#endif

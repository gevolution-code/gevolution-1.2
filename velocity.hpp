/////////////////////////
// velocity.hpp
//////////////////////////
//
// Utilities for the velocity field
//
// Author: Goran Jelic-Cizmek (Université de Genève)
// Author: Francesca Lepori (SISSA Trieste & INFN Trieste & Université de Genève)
// Author: Julian Adamek (Queen Mary University of London)
//
// Last modified: April 2019
//
//////////////////////////

#ifndef VELOCITY_HEADER
#define VELOCITY_HEADER

#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>

/**
    differential equation for the growth rate D_1
**/

static double g1(double a, void *p)
{
  struct cosmology par = *(struct cosmology *) p;
  double g1 = 1./(pow(a, 3)*pow(1 - par.Omega_m + par.Omega_m/pow(a, 3),1.5));
  return g1;
}

/**
    computes and stores all the background functions
**/

double D1_prime(
    struct cosmology par, double a
)
{
    gsl_function G1;
    G1.function = &g1, G1.params = &par;
    const double delta_a = 0.005;
    double a_minus = a*(1 - delta_a);
    double a_plus  = a*(1 + delta_a);
    gsl_integration_workspace *w;
    double prec = 1E-5;
    double result, error;

    w = gsl_integration_workspace_alloc(10000);
    gsl_integration_qag (&G1, 0, a_plus, 0, prec, 10000,
		         GSL_INTEG_GAUSS61, w, &result, &error);

    double result_plus = (5./2*par.Omega_m*sqrt(1 - par.Omega_m + par.Omega_m/pow(a_plus, 3)))*result;

    gsl_integration_qag (&G1, 0, a_minus, 0, prec, 10000,
                         GSL_INTEG_GAUSS61, w, &result, &error);

    double result_minus = (5./2*par.Omega_m*sqrt(1 - par.Omega_m + par.Omega_m/pow(a_minus, 3)))*result;

    double result_der = (result_plus -  result_minus)/(2*delta_a*a);

    result = result_der*Hconf(a, 1, par)*a;

    gsl_integration_workspace_free(w);

    return result;
}

//////////////////////////
// compute_vi_rescaled
//////////////////////////
// Description:
//   Compute the velocity field as v^i = T^i_0/T^0_0, if a = 1 then vi = a v^i
//   If T^0_0 = 0 the velocity field is set to be the one at the previous time step,
//   rescaled as v^i(a) = v^i(a_past) a*Hconf(a) dD1/da (velocity method = rescaled)
//   [see G. Jelic-Cizmek, F. Lepori, J. Adamek, and R. Durrer, JCAP 1809, 006 (2018)]
//
// Arguments:
//   cosmo      structure containing the cosmological parameters
//   vi         reference to the velocity field, contains (a^3 T^i_0)
//   source     reference to the field source (a^4 T^0_0)
//   Ti0        reference to the field Ti0 (a^3 T^i_0)
//   a          scale factor at current time step
//   a_old      scale factor at previous time step
//
// Returns:
//
//////////////////////////

void compute_vi_rescaled(cosmology & cosmo, Field<Real> * vi, Field<Real> * source, Field<Real> * Ti0, double a = 1., double a_old = 1.)
{
	Site xvi(vi->lattice());

	Real rescale = D1_prime(cosmo, a)/D1_prime(cosmo, a_old)*a/a_old;

	for(xvi.first(); xvi.test(); xvi.next())
	{
		if ((*source)(xvi) < 1.E-30)
		{
			(*vi)(xvi,0) *= rescale;
			(*vi)(xvi,1) *= rescale;
			(*vi)(xvi,2) *= rescale;
		}
		else
		{
			(*vi)(xvi,0) = (*Ti0)(xvi,0) / (*source)(xvi) / a;
			(*vi)(xvi,1) = (*Ti0)(xvi,1) / (*source)(xvi) / a;
			(*vi)(xvi,2) = (*Ti0)(xvi,2) / (*source)(xvi) / a;
		}
	}
}

#endif


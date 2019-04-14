//////////////////////////
// background.hpp
//////////////////////////
// 
// code components related to background evolution
//
// Author: Julian Adamek (Université de Genève & Observatoire de Paris & Queen Mary University of London)
//
// Last modified: September 2018
//
//////////////////////////

#ifndef BACKGROUND_HEADER
#define BACKGROUND_HEADER

#include <gsl/gsl_integration.h>


double FermiDiracIntegrand(double q, void * w)
{
	return q * q * sqrt(q * q + *(double *)w) / (exp(q) + 1.0l);
}

//////////////////////////
// FermiDiracIntegral
//////////////////////////
// Description:
//   computes the integral of the relativistic Fermi-Dirac distribution
// 
// Arguments:
//   w          parameter in the F-D distribution, "(m a / kB T)^2"
//
// Returns: value for the integral
// 
//////////////////////////

double FermiDiracIntegral(double &w)
{
	double result;
	gsl_function f;
	double err;
	size_t n;
	
	f.function = &FermiDiracIntegrand;
	f.params = &w;
		
	gsl_integration_qng(&f, 0.0l, 24.0l, 5.0e-7, 1.0e-7, &result, &err, &n);
	
	return result;
}


//////////////////////////
// bg_ncdm (1)
//////////////////////////
// Description:
//   computes the background model for one ncdm species by integrating the relativistic
//   Fermi-Dirac distribution
// 
// Arguments:
//   a          scale factor at which to compute the background model
//   cosmo      structure containing the cosmological parameters
//   p          index of the ncdm species
//
// Returns: value for the background model
// 
//////////////////////////

double bg_ncdm(const double a, const cosmology cosmo, const int p)
{
	if (p < 0 || p >= cosmo.num_ncdm)
		return 0;
	else
	{
		double w = a * cosmo.m_ncdm[p] / (pow(cosmo.Omega_g * cosmo.h * cosmo.h / C_PLANCK_LAW, 0.25) * cosmo.T_ncdm[p] * C_BOLTZMANN_CST);
		w *= w;
		
		return FermiDiracIntegral(w) * cosmo.Omega_ncdm[p] * pow(cosmo.Omega_g * cosmo.h * cosmo.h / C_PLANCK_LAW, 0.25) * cosmo.T_ncdm[p] * C_BOLTZMANN_CST / cosmo.m_ncdm[p] / C_FD_NORM / a;
	}
}


//////////////////////////
// bg_ncdm (2)
//////////////////////////
// Description:
//   computes the background model for all ncdm species by integrating the relativistic
//   Fermi-Dirac distribution
// 
// Arguments:
//   a          scale factor at which to compute the background model
//   cosmo      structure containing the cosmological parameters
//
// Note:
//   For optimization, the last value of a is stored in a static variable such that
//   multiple calls at the same value of a will not result in multiple integrations
//   being carried out. This assumes that the cosmological model should not change!
//
// Returns: value for the background model
// 
//////////////////////////

double bg_ncdm(const double a, const cosmology cosmo)
{
	double w;
	static double result = -1.0;
	static double a_prev = -1.0;
	
	if (a != a_prev)
	{
		result = 0.0;
		a_prev = a;
		
		for (int p = 0; p < cosmo.num_ncdm; p++)
			result += bg_ncdm(a, cosmo, p);
	}
	
	return result;
}


//////////////////////////
// Hconf
//////////////////////////
// Description:
//   computes the conformal Hubble rate at given scale factor
// 
// Arguments:
//   a          scale factor
//   fourpiG    "4 pi G"
//   cosmo      structure containing the cosmological parameters
//
// Returns: conformal Hubble rate
// 
//////////////////////////

double Hconf(const double a, const double fourpiG, const cosmology cosmo)
{	
	return sqrt((2. * fourpiG / 3.) * (((cosmo.Omega_cdm + cosmo.Omega_b + bg_ncdm(a, cosmo)) / a) + (cosmo.Omega_Lambda * a * a) + (cosmo.Omega_rad / a / a) + (cosmo.Omega_fld * exp(3. * cosmo.wa_fld * (a - 1.)) / pow(a, 1. + 3. * (cosmo.w0_fld + cosmo.wa_fld)))));
}


double Omega_m(const double a, const cosmology cosmo) { return cosmo.Omega_m / (cosmo.Omega_cdm + cosmo.Omega_b + bg_ncdm(a, cosmo) + cosmo.Omega_Lambda * a * a * a + cosmo.Omega_rad / a + cosmo.Omega_fld * exp(3. * cosmo.wa_fld * (a - 1.)) / pow(a, 3. * (cosmo.w0_fld + cosmo.wa_fld))); }

double Omega_rad(const double a, const cosmology cosmo) { return (cosmo.Omega_rad + (bg_ncdm(a, cosmo) + cosmo.Omega_cdm + cosmo.Omega_b - cosmo.Omega_m) * a) / ((cosmo.Omega_cdm + cosmo.Omega_b + bg_ncdm(a, cosmo)) * a + cosmo.Omega_Lambda * a * a * a * a + cosmo.Omega_rad + cosmo.Omega_fld * exp(3. * cosmo.wa_fld * (a - 1.)) / pow(a, 3. * (cosmo.w0_fld + cosmo.wa_fld) - 1.)); }

double Omega_Lambda(const double a, const cosmology cosmo) { return cosmo.Omega_Lambda / ((cosmo.Omega_cdm + cosmo.Omega_b + bg_ncdm(a, cosmo)) / a / a / a + cosmo.Omega_Lambda + cosmo.Omega_rad / a / a / a / a + cosmo.Omega_fld * exp(3. * cosmo.wa_fld * (a - 1.)) / pow(a, 3. + 3. * (cosmo.w0_fld + cosmo.wa_fld))); }


//////////////////////////
// rungekutta4bg
//////////////////////////
// Description:
//   integrates the Friedmann equation for the background model using a fourth-order
//   Runge-Kutta method
// 
// Arguments:
//   a          scale factor (will be advanced by dtau)
//   fourpiG    "4 pi G"
//   cosmo      structure containing the cosmological parameters
//   dtau       time step by which the scale factor should be advanced
//
// Returns:
// 
//////////////////////////

void rungekutta4bg(double &a, const double fourpiG, const cosmology cosmo, const double dtau)
{
	double k1a, k2a, k3a, k4a;

	k1a = a * Hconf(a, fourpiG, cosmo); 
	k2a = (a + k1a * dtau / 2.) * Hconf(a + k1a * dtau / 2., fourpiG, cosmo);
	k3a = (a + k2a * dtau / 2.) * Hconf(a + k2a * dtau / 2., fourpiG, cosmo);
	k4a = (a + k3a * dtau) * Hconf(a + k3a * dtau, fourpiG, cosmo);

	a += dtau * (k1a + 2. * k2a + 2. * k3a + k4a) / 6.;
}


double particleHorizonIntegrand(double sqrta, void * cosmo)
{
	return 2. / (sqrta * Hconf(sqrta*sqrta, 1., *(cosmology *)cosmo));
}

//////////////////////////
// particleHorizon
//////////////////////////
// Description:
//   computes the particle horizon (tau) at given scale factor
// 
// Arguments:
//   a          scale factor
//   fourpiG    "4 pi G"
//   cosmo      structure containing the cosmological parameters
//
// Returns: particle horizon (tau)
// 
//////////////////////////

double particleHorizon(const double a, const double fourpiG, cosmology & cosmo)
{
	double result;
	gsl_function f;
	double err;
	size_t n;
	
	f.function = &particleHorizonIntegrand;
	f.params = &cosmo;
	
	gsl_integration_qng(&f, sqrt(a) * 1.0e-7, sqrt(a), 5.0e-7, 1.0e-7, &result, &err, &n);
	
	return result / sqrt(fourpiG);
}

#endif


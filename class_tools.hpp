//////////////////////////
// class_tools.hpp
//////////////////////////
// 
// interface to linear Boltzmann code CLASS
//
// Author: Julian Adamek (Université de Genève & Observatoire de Paris & Queen Mary University of London)
//
// Last modified: November 2019
//
//////////////////////////

#ifndef CLASS_TOOLS_HEADER
#define CLASS_TOOLS_HEADER

#ifdef HAVE_CLASS

#define thermo thermodynamics
#define perturbs perturbations
#define nonlinear fourier
#define transfers transfer
#define spectra harmonic
#define perturb_init perturbations_init
#define perturb_free perturbations_free
#define perturb_output_titles perturbations_output_titles
#define perturb_output_data perturbations_output_data

#include <gsl/gsl_spline.h>
#include "parser.hpp"

using namespace std;
using namespace LATfield2;

//////////////////////////
// initializeCLASSstructures
//////////////////////////
// Description:
//   initializes CLASS structures containing interpolation tables for various transfer functions
// 
// Arguments:
//   sim               simulation metadata structure
//   ic                settings for IC generation
//   cosmo             cosmological parameter structure
//   class_background  CLASS structure that will contain the background
//   class_thermo      CLASS structure that will contain thermodynamics
//   class_perturbs    CLASS structure that will contain perturbations
//   params            pointer to array of precision settings (optional)
//   numparam          number of precision settings (default 0)
//   output_value      CLASS parameter value specifying the output (optional)
//
// Returns:
// 
//////////////////////////

void initializeCLASSstructures(metadata & sim, icsettings & ic, cosmology & cosmo, background & class_background, thermo & class_thermo, perturbs & class_perturbs, parameter * params = NULL, int numparam = 0, const char * output_value = "dTk, vTk")
{
	precision class_precision;
	transfers class_transfers;
  	primordial class_primordial;
	nonlinear class_nonlinear;
	spectra class_spectra;
  	lensing class_lensing;
	distortions class_distortions;
  	output class_output;
	file_content class_filecontent;
	ErrorMsg class_errmsg;
	char filename[] = "initializeCLASSstructures";
	char tmp[8 * PARAM_MAX_LENGTH];
	double recfast_z_initial;
	double perturb_sampling_stepsize;
	int recfast_Nz0;
	int i;
	int num_entries = 19;
#ifdef CLASS_K_PER_DECADE_FOR_PK
	int k_per_decade_for_pk;
	if (numparam == 0 || !parseParameter(params, numparam, "k_per_decade_for_pk", k_per_decade_for_pk))
	{
		num_entries += 1;
		k_per_decade_for_pk = CLASS_K_PER_DECADE_FOR_PK;
	}
#endif

	if (cosmo.num_ncdm > 0) num_entries += 3;
	if (parallel.isRoot()) num_entries += 7;
	if ((1.015 * ic.z_ic + 0.01) > 9999.)
	{
		if (numparam == 0 || !parseParameter(params, numparam, "recfast_z_initial", recfast_z_initial))
		{
			num_entries += 1;
			recfast_z_initial = ic.z_ic * 1.02;
		}
		if (numparam == 0 || !parseParameter(params, numparam, "recfast_Nz0", recfast_Nz0))
		{
			num_entries += 1;
			recfast_Nz0 = 2 * (int) ceil(ic.z_ic * 1.02);
		}
	}
	if (numparam == 0 || !parseParameter(params, numparam, "perturb_sampling_stepsize", perturb_sampling_stepsize))
	{
		num_entries += 1;
		perturb_sampling_stepsize = 0.01;		
	}
	
	num_entries += numparam;

	parser_init(&class_filecontent, num_entries, filename, class_errmsg);

	for (i = 0; i < num_entries; i++)
		class_filecontent.read[i] = _FALSE_;

	i = 0;

	sprintf(class_filecontent.name[i], "root");
	sprintf(class_filecontent.value[i++], "%s%s_class", sim.output_path, sim.basename_generic);

	sprintf(class_filecontent.name[i], "k_pivot");
	sprintf(class_filecontent.value[i++], "%e", ic.k_pivot);

	sprintf(class_filecontent.name[i], "A_s");
	sprintf(class_filecontent.value[i++], "%e", ic.A_s);

	sprintf(class_filecontent.name[i], "n_s");
	sprintf(class_filecontent.value[i++], "%f", ic.n_s);

	sprintf(class_filecontent.name[i], "z_pk");
	if (ic.z_ic > sim.z_in)
		sprintf(class_filecontent.value[i++], "%f, %f, 0", 1.015 * ic.z_ic + 0.01, sim.z_in);
	else
		sprintf(class_filecontent.value[i++], "%f, 0", 1.015 * ic.z_ic + 0.01);

	sprintf(class_filecontent.name[i], "output");
	sprintf(class_filecontent.value[i++], "%s", output_value);

	sprintf(class_filecontent.name[i], "gauge");
	sprintf(class_filecontent.value[i++], "Newtonian");

	sprintf(class_filecontent.name[i], "P_k_ini type");
	sprintf(class_filecontent.value[i++], "analytic_Pk");

	sprintf(class_filecontent.name[i], "P_k_max_h/Mpc");
	sprintf(class_filecontent.value[i++], "%f", 2. * M_PI * (double) sim.numpts / sim.boxsize);

	sprintf(class_filecontent.name[i], "h");
	sprintf(class_filecontent.value[i++], "%f", cosmo.h);

	sprintf(class_filecontent.name[i], "Omega_cdm");
	sprintf(class_filecontent.value[i++], "%e", cosmo.Omega_cdm);

	sprintf(class_filecontent.name[i], "Omega_b");
	sprintf(class_filecontent.value[i++], "%e", cosmo.Omega_b);

	sprintf(class_filecontent.name[i], "Omega_g");
	sprintf(class_filecontent.value[i++], "%e", cosmo.Omega_g);

	sprintf(class_filecontent.name[i], "Omega_ur");
	sprintf(class_filecontent.value[i++], "%e", cosmo.Omega_ur);

	sprintf(class_filecontent.name[i], "Omega_fld");
	sprintf(class_filecontent.value[i++], "%e", cosmo.Omega_fld);

	sprintf(class_filecontent.name[i], "w0_fld");
	sprintf(class_filecontent.value[i++], "%g", cosmo.w0_fld);

	sprintf(class_filecontent.name[i], "wa_fld");
	sprintf(class_filecontent.value[i++], "%g", cosmo.wa_fld);

	sprintf(class_filecontent.name[i], "cs2_fld");
	sprintf(class_filecontent.value[i++], "%g", cosmo.cs2_fld);

	sprintf(class_filecontent.name[i], "N_ncdm");
	sprintf(class_filecontent.value[i++], "%d", cosmo.num_ncdm);

	sprintf(class_filecontent.name[i], "perturb_sampling_stepsize");
	sprintf(class_filecontent.value[i++], "%g", perturb_sampling_stepsize);

#ifdef CLASS_K_PER_DECADE_FOR_PK
	sprintf(class_filecontent.name[i], "k_per_decade_for_pk");
	sprintf(class_filecontent.value[i++], "%d", k_per_decade_for_pk);
#endif

	if ((1.015 * ic.z_ic + 0.01) > 9999.)
	{
		sprintf(class_filecontent.name[i], "recfast_z_initial");
		sprintf(class_filecontent.value[i++], "%e", recfast_z_initial);

		sprintf(class_filecontent.name[i], "recfast_Nz0");
		sprintf(class_filecontent.value[i++], "%d", recfast_Nz0);
	}

	if (parallel.isRoot())
	{
		sprintf(class_filecontent.name[i], "background_verbose");
		sprintf(class_filecontent.value[i++], "1");

		sprintf(class_filecontent.name[i], "perturbations_verbose");
		sprintf(class_filecontent.value[i++], "1");

		sprintf(class_filecontent.name[i], "spectra_verbose");
		sprintf(class_filecontent.value[i++], "1");

		sprintf(class_filecontent.name[i], "thermodynamics_verbose");
		sprintf(class_filecontent.value[i++], "1");

		sprintf(class_filecontent.name[i], "transfer_verbose");
		sprintf(class_filecontent.value[i++], "1");

		sprintf(class_filecontent.name[i], "primordial_verbose");
		sprintf(class_filecontent.value[i++], "1");

		sprintf(class_filecontent.name[i], "nonlinear_verbose");
		sprintf(class_filecontent.value[i++], "1");
	}
	
	while (numparam > 0)
	{
		numparam--;
		if (!params[numparam].used)
		{
			sprintf(class_filecontent.name[i], "%s", params[numparam].name);
			sprintf(class_filecontent.value[i++], "%s", params[numparam].value);
		}
	}

	if (cosmo.num_ncdm > 0)
	{
		sprintf(class_filecontent.name[num_entries-3], "m_ncdm");
		sprintf(class_filecontent.value[num_entries-3], "%f", cosmo.m_ncdm[0]);
		for (i = 1; i < cosmo.num_ncdm; i++)
		{
			sprintf(tmp, "%s, %f", class_filecontent.value[num_entries-3], cosmo.m_ncdm[i]);
			sprintf(class_filecontent.value[num_entries-3], "%s", tmp);
		}

		sprintf(class_filecontent.name[num_entries-2], "T_ncdm");
		sprintf(class_filecontent.value[num_entries-2], "%f", cosmo.T_ncdm[0]);
		for (i = 1; i < cosmo.num_ncdm; i++)
		{
			sprintf(tmp, "%s, %f", class_filecontent.value[num_entries-2], cosmo.T_ncdm[i]);
			sprintf(class_filecontent.value[num_entries-2], "%s", tmp);
		}

		sprintf(class_filecontent.name[num_entries-1], "deg_ncdm");
		sprintf(class_filecontent.value[num_entries-1], "%f", cosmo.deg_ncdm[0]);
		for (i = 1; i < cosmo.num_ncdm; i++)
		{
			sprintf(tmp, "%s, %f", class_filecontent.value[num_entries-1], cosmo.deg_ncdm[i]);
			sprintf(class_filecontent.value[num_entries-1], "%s", tmp);
		}
	}

	COUT << " gevolution is calling CLASS..." << endl << endl;

	if (input_read_from_file(&class_filecontent, &class_precision, &class_background, &class_thermo, &class_perturbs, &class_transfers, &class_primordial, &class_spectra, &class_nonlinear, &class_lensing, &class_distortions, &class_output, class_errmsg) == _FAILURE_)
	{
		COUT << " error: calling input_init from CLASS library failed!" << endl << " following error message was passed: " << class_errmsg << endl;
		parallel.abortForce();
	}

	parser_free(&class_filecontent);

	if (background_init(&class_precision, &class_background) == _FAILURE_)
	{
		COUT << " error: calling background_init from CLASS library failed!" << endl << " following error message was passed: " << class_background.error_message << endl;
		parallel.abortForce();
	}

	if (thermodynamics_init(&class_precision, &class_background, &class_thermo) == _FAILURE_)
	{
		COUT << " error: calling thermodynamics_init from CLASS library failed!" << endl << " following error message was passed: " << class_thermo.error_message << endl;
		parallel.abortForce();
	}

	if (perturb_init(&class_precision, &class_background, &class_thermo, &class_perturbs) == _FAILURE_)
	{
		COUT << " error: calling perturb_init from CLASS library failed!" << endl << " following error message was passed: " << class_perturbs.error_message << endl;
		parallel.abortForce();
	}

	COUT << endl << " CLASS structures initialized successfully." << endl;
}


//////////////////////////
// freeCLASSstructures
//////////////////////////
// Description:
//   frees CLASS structures containing interpolation tables for various transfer functions
// 
// Arguments:
//   class_background  CLASS structure that contains the background
//   class_thremo      CLASS structure that contains thermodynamics
//   class_perturbs    CLASS structure that contains perturbations
//
// Returns:
// 
//////////////////////////

void freeCLASSstructures(background & class_background, thermo & class_thermo, perturbs & class_perturbs)
{
	if (perturb_free(&class_perturbs) == _FAILURE_)
	{
		COUT << " error: calling perturb_free from CLASS library failed!" << endl << " following error message was passed: " << class_perturbs.error_message << endl;
		parallel.abortForce();
	}
	
	if (thermodynamics_free(&class_thermo) == _FAILURE_)
	{
		COUT << " error: calling thermodynamics_free from CLASS library failed!" << endl << " following error message was passed: " << class_thermo.error_message << endl;
		parallel.abortForce();
	}

	if (background_free(&class_background) == _FAILURE_)
	{
		COUT << " error: calling background_free from CLASS library failed!" << endl << " following error message was passed: " << class_background.error_message << endl;
		parallel.abortForce();
	}
}


//////////////////////////
// loadTransferFunctions (2)
//////////////////////////
// Description:
//   loads a set of tabulated transfer functions from some precomputed CLASS structures
// 
// Arguments:
//   class_background  CLASS structure that contains the background
//   class_perturbs    CLASS structure that contains the perturbations
//   tk_delta          will point to the gsl_spline which holds the tabulated
//                     transfer function for delta (memory will be allocated)
//   tk_theta          will point to the gsl_spline which holds the tabulated
//                     transfer function for theta (memory will be allocated)
//   qname             string containing the name of the component (e.g. "cdm"); if no string is
//                     specified, the transfer functions for phi and psi are returned instead!
//   boxsize           comoving box size (in the same units as used in the CLASS output)
//   z                 redshift at which the transfer functions are to be obtained
//   h                 conversion factor between 1/Mpc and h/Mpc (theta is in units of 1/Mpc)
//
// Returns:
// 
//////////////////////////

void loadTransferFunctions(background & class_background, perturbs & class_perturbs, gsl_spline * & tk_delta, gsl_spline * & tk_theta, const char * qname, const double boxsize, const double z, double h)
{
	int cols = 0, dcol = -1, tcol = -1, kcol = -1;
	double * k;
	double * tk_d;
	double * tk_t;
	double * data;
	char coltitles[_MAXTITLESTRINGLENGTH_] = {0};
	char dname[16];
	char tname[16];
	char kname[8];
	char * ptr;

	perturb_output_titles(&class_background, &class_perturbs, class_format, coltitles);

	if (qname != NULL)
	{
		sprintf(dname, "d_%s", qname);
		sprintf(tname, "t_%s", qname);
		h /= boxsize;
    }
	else
	{
		sprintf(dname, "phi");
		sprintf(tname, "psi");
		h = 1.;
	}
	sprintf(kname, "k");

	ptr = strtok(coltitles, _DELIMITER_);
	while (ptr != NULL)
	{
    	if (strncmp(ptr, dname, strlen(dname)) == 0) dcol = cols;
		else if (strncmp(ptr, tname, strlen(tname)) == 0) tcol = cols;
		else if (strncmp(ptr, kname, strlen(kname)) == 0) kcol = cols;
		cols++;
    	ptr = strtok(NULL, _DELIMITER_);
  	}

	if (dcol < 0 || tcol < 0 || kcol < 0)
	{
		COUT << " error in loadTransferFunctions (HAVE_CLASS)! Unable to identify requested columns!" << endl;
		parallel.abortForce();
	}

	data = (double *) malloc(sizeof(double) * cols*class_perturbs.k_size[class_perturbs.index_md_scalars]);
	k = (double *) malloc(sizeof(double) * class_perturbs.k_size[class_perturbs.index_md_scalars]);
	tk_d = (double *) malloc(sizeof(double) * class_perturbs.k_size[class_perturbs.index_md_scalars]);
	tk_t = (double *) malloc(sizeof(double) * class_perturbs.k_size[class_perturbs.index_md_scalars]);

	perturb_output_data(&class_background, &class_perturbs, class_format, z, cols, data);

	for (int i = 0; i < class_perturbs.k_size[class_perturbs.index_md_scalars]; i++)
	{
		k[i] = data[i*cols + kcol] * boxsize;
		tk_d[i] = data[i*cols + dcol];
		tk_t[i] = data[i*cols + tcol] / h;
		if (i > 0)
		{
			if (k[i] < k[i-1])
			{
				COUT << " error in loadTransferFunctions (HAVE_CLASS)! k-values are not strictly ordered." << endl;
				parallel.abortForce();
			}
		}
	}

	free(data);
	
	tk_delta = gsl_spline_alloc(gsl_interp_cspline, class_perturbs.k_size[class_perturbs.index_md_scalars]);
	tk_theta = gsl_spline_alloc(gsl_interp_cspline, class_perturbs.k_size[class_perturbs.index_md_scalars]);
	
	gsl_spline_init(tk_delta, k, tk_d, class_perturbs.k_size[class_perturbs.index_md_scalars]);
	gsl_spline_init(tk_theta, k, tk_t, class_perturbs.k_size[class_perturbs.index_md_scalars]);
	
	free(k);
	free(tk_d);
	free(tk_t);
}

#endif

#endif


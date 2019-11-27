//////////////////////////
// radiation.hpp
//////////////////////////
// 
// code components related to radiation and linear relativistic species
//
// Author: Julian Adamek (Université de Genève & Observatoire de Paris & Queen Mary University of London)
//
// Last modified: November 2019
//
//////////////////////////

#ifndef RADIATION_HEADER
#define RADIATION_HEADER

#ifdef HAVE_CLASS

//////////////////////////
// projection_T00_project (radiation module)
//////////////////////////
// Description:
//   provides a realization of the linear density field of radiation and
//   non-cold species using linear transfer functions precomputed with CLASS;
//   the contributions for the various species are included only until some
//   individual redshift values are reached (after which no linear treatment
//   is requested)
// 
// Arguments:
//   class_background  CLASS structure that contains the background
//   class_perturbs    CLASS structure that contains the perturbations
//   source            reference to field that will contain the realization
//   scalarFT          reference to Fourier image of that field
//   plan_source       pointer to FFT planner
//   sim               simulation metadata structure
//   ic                settings for IC generation (contains the random seed)
//   cosmo             cosmological parameter structure
//   fourpiG           4 pi G (in code units)
//   a                 scale factor
//   coeff             multiplicative coefficient (default 1)
//
// Returns:
// 
//////////////////////////

void projection_T00_project(background & class_background, perturbs & class_perturbs, Field<Real> & source, Field<Cplx> & scalarFT, PlanFFT<Cplx> * plan_source, metadata & sim, icsettings & ic, cosmology & cosmo, const double fourpiG, double a, double coeff = 1.)
{
	gsl_spline * tk1 = NULL;
	gsl_spline * tk2 = NULL;
	double * delta = NULL;
	double * k = NULL;
	char ncdm_name[8];
	int i, p, n = 0;
	double rescale, Omega_ncdm = 0., Omega_rad = 0., Omega_fld = 0.;
	Site x(source.lattice());
	rKSite kFT(scalarFT.lattice());

	if (a < 1. / (sim.z_switch_deltarad + 1.) && cosmo.Omega_g > 0 && sim.radiation_flag == 1)
	{
		loadTransferFunctions(class_background, class_perturbs, tk1, tk2, "g", sim.boxsize, (1. / a) - 1., cosmo.h);
		Omega_rad += cosmo.Omega_g;

		n = tk1->size;
		delta = (double *) malloc(n * sizeof(double));
		k = (double *) malloc(n * sizeof(double));
		
		for (i = 0; i < n; i++)
		{
			delta[i] = -tk1->y[i] * coeff * cosmo.Omega_g * M_PI * sqrt(Pk_primordial(tk1->x[i] * cosmo.h / sim.boxsize, ic) / tk1->x[i]) / tk1->x[i] / a;
			k[i] = tk1->x[i];
		}

		gsl_spline_free(tk1);
		gsl_spline_free(tk2);
	}

	if (a < 1. / (sim.z_switch_deltarad + 1.) && cosmo.Omega_ur > 0 && sim.radiation_flag == 1)
	{
		loadTransferFunctions(class_background, class_perturbs, tk1, tk2, "ur", sim.boxsize, (1. / a) - 1., cosmo.h);
		Omega_rad += cosmo.Omega_ur;

		if (delta == NULL)
		{
			n = tk1->size;
			delta = (double *) malloc(n * sizeof(double));
			k = (double *) malloc(n * sizeof(double));

			for (i = 0; i < n; i++)
			{
				delta[i] = -tk1->y[i] * coeff * cosmo.Omega_ur * M_PI * sqrt(Pk_primordial(tk1->x[i] * cosmo.h / sim.boxsize, ic) / tk1->x[i]) / tk1->x[i] / a;
				k[i] = tk1->x[i];
			}
		}
		else
		{
			for (i = 0; i < n; i++)
				delta[i] -= tk1->y[i] * coeff * cosmo.Omega_ur * M_PI * sqrt(Pk_primordial(tk1->x[i] * cosmo.h / sim.boxsize, ic) / tk1->x[i]) / tk1->x[i] / a;
		}

		gsl_spline_free(tk1);
		gsl_spline_free(tk2);
	}

	if (a < 1. && cosmo.Omega_fld > 0 && sim.fluid_flag == 1)
	{
		loadTransferFunctions(class_background, class_perturbs, tk1, tk2, "fld", sim.boxsize, (1. / a) - 1., cosmo.h);
		Omega_fld = cosmo.Omega_fld / pow(a, 3. * cosmo.w0_fld);

		if (delta == NULL)
		{
			n = tk1->size;
			delta = (double *) malloc(n * sizeof(double));
			k = (double *) malloc(n * sizeof(double));

			for (i = 0; i < n; i++)
			{
				delta[i] = -tk1->y[i] * coeff * Omega_fld * M_PI * sqrt(Pk_primordial(tk1->x[i] * cosmo.h / sim.boxsize, ic) / tk1->x[i]) / tk1->x[i];
				k[i] = tk1->x[i];
			}
		}
		else
		{
			for (i = 0; i < n; i++)
				delta[i] -= tk1->y[i] * coeff * Omega_fld * M_PI * sqrt(Pk_primordial(tk1->x[i] * cosmo.h / sim.boxsize, ic) / tk1->x[i]) / tk1->x[i];
		}

		gsl_spline_free(tk1);
		gsl_spline_free(tk2);
	}

	for (p = 0; p < cosmo.num_ncdm; p++)
	{
		if (a < 1. / (sim.z_switch_deltancdm[p] + 1.) && cosmo.Omega_ncdm[p] > 0)
		{
			sprintf(ncdm_name, "ncdm[%d]", p);
			loadTransferFunctions(class_background, class_perturbs, tk1, tk2, ncdm_name, sim.boxsize, (1. / a) - 1., cosmo.h);
			rescale = bg_ncdm(a, cosmo, p);
			Omega_ncdm += rescale;

			if (delta == NULL)
			{
				n = tk1->size;
				delta = (double *) malloc(n * sizeof(double));
				k = (double *) malloc(n * sizeof(double));

				for (i = 0; i < n; i++)
				{
					delta[i] = -tk1->y[i] * coeff * rescale * M_PI * sqrt(Pk_primordial(tk1->x[i] * cosmo.h / sim.boxsize, ic) / tk1->x[i]) / tk1->x[i];
					k[i] = tk1->x[i];
				}
			}
			else
			{
				for (i = 0; i < n; i++)
					delta[i] -= tk1->y[i] * coeff * rescale * M_PI * sqrt(Pk_primordial(tk1->x[i] * cosmo.h / sim.boxsize, ic) / tk1->x[i]) / tk1->x[i];
			}

			gsl_spline_free(tk1);
			gsl_spline_free(tk2);
		}
	}

	if (n > 0)
	{
		if (sim.gr_flag == 0) // add gauge correction for N-body gauge
		{
			loadTransferFunctions(class_background, class_perturbs, tk1, tk2, "tot", sim.boxsize, (1. / a) - 1., cosmo.h);
			rescale = Hconf(a, fourpiG, cosmo);

			for (i = 0; i < n; i++)
				delta[i] -= coeff * (4. * Omega_rad / a + 3. * Omega_ncdm + 3. * (1. + cosmo.w0_fld) * Omega_fld) * rescale * M_PI * tk2->y[i] * sqrt(Pk_primordial(tk2->x[i] * cosmo.h / sim.boxsize, ic) / tk2->x[i]) / tk2->x[i] / tk2->x[i] / tk2->x[i];

			gsl_spline_free(tk1);
			gsl_spline_free(tk2);
		}

		tk1 = gsl_spline_alloc(gsl_interp_cspline, n);
		gsl_spline_init(tk1, k, delta, n);		

		generateRealization(scalarFT, 0., tk1, (unsigned int) ic.seed, ic.flags & ICFLAG_KSPHERE);
		plan_source->execute(FFT_BACKWARD);

		gsl_spline_free(tk1);
		free(delta);
		free(k);

		for (x.first(); x.test(); x.next())
			source(x) += Omega_ncdm;
	}
}


//////////////////////////
// prepareFTchiLinear
//////////////////////////
// Description:
//   provides a (Fourier-space) realization of chi (generated by radiation and
//   non-cold species) from the linear transfer functions precomputed with CLASS
// 
// Arguments:
//   class_background  CLASS structure that contains the background
//   class_perturbs    CLASS structure that contains the perturbations
//   scalarFT          reference to Fourier image of field; will contain the
//                     (Fourier image) of the realization
//   sim               simulation metadata structure
//   ic                settings for IC generation (contains the random seed)
//   cosmo             cosmological parameter structure
//   fourpiG           4 pi G (in code units)
//   a                 scale factor
//   coeff             multiplicative coefficient (default 1)
//
// Returns:
// 
//////////////////////////

void prepareFTchiLinear(background & class_background, perturbs & class_perturbs, Field<Cplx> & scalarFT, metadata & sim, icsettings & ic, cosmology & cosmo, const double fourpiG, double a, double coeff = 1.)
{
	gsl_spline * tk1 = NULL;
	gsl_spline * tk2 = NULL;
	double * chi = NULL;
	int i;
	rKSite k(scalarFT.lattice());

	loadTransferFunctions(class_background, class_perturbs, tk1, tk2, NULL, sim.boxsize, (1. / a) - 1., cosmo.h);

	chi = (double *) malloc(tk1->size * sizeof(double));

	for (i = 0; i < tk1->size; i++)
		chi[i] = (tk2->y[i] - tk1->y[i]) * coeff * M_PI * sqrt(Pk_primordial(tk1->x[i] * cosmo.h / sim.boxsize, ic) / tk1->x[i]) / tk1->x[i];

	if (sim.gr_flag == 0) // add gauge correction for N-body gauge
	{
		double * l1 = (double *) malloc(tk1->size * sizeof(double));
		double * l2 = (double *) malloc(tk1->size * sizeof(double));
		double * l3 = (double *) malloc(tk1->size * sizeof(double));
		double * l4 = (double *) malloc(tk1->size * sizeof(double));
		double * l5 = (double *) malloc(tk1->size * sizeof(double));
		double Hconf1 = Hconf(0.99 * a, fourpiG, cosmo);
		double Hconf2 = Hconf(0.995 * a, fourpiG, cosmo);
		double Hconf3 = Hconf(a, fourpiG, cosmo);
		double Hconf4 = Hconf(1.005 * a, fourpiG, cosmo);
		double Hconf5 = Hconf(1.01 * a, fourpiG, cosmo);

		for (i = 0; i < tk1->size; i++)
			l3[i] = -tk1->y[i];

		gsl_spline_free(tk1);
		gsl_spline_free(tk2);

		loadTransferFunctions(class_background, class_perturbs, tk1, tk2, NULL, sim.boxsize, (1. / (1.005 * a)) - 1., cosmo.h);

		for (i = 0; i < tk1->size; i++)
			l4[i] = -tk1->y[i];

		gsl_spline_free(tk1);
		gsl_spline_free(tk2);

		loadTransferFunctions(class_background, class_perturbs, tk1, tk2, "tot", sim.boxsize, (1. / (1.005 * a)) - 1., cosmo.h);

		for (i = 0; i < tk1->size; i++)
			l4[i] -= tk2->y[i] * Hconf4 / tk2->x[i] / tk2->x[i];

		gsl_spline_free(tk1);
		gsl_spline_free(tk2);

		loadTransferFunctions(class_background, class_perturbs, tk1, tk2, NULL, sim.boxsize, (1. / (1.01 * a)) - 1., cosmo.h);

		for (i = 0; i < tk1->size; i++)
			l5[i] = -tk1->y[i];

		gsl_spline_free(tk1);
		gsl_spline_free(tk2);

		loadTransferFunctions(class_background, class_perturbs, tk1, tk2, "tot", sim.boxsize, (1. / (1.01 * a)) - 1., cosmo.h);

		for (i = 0; i < tk1->size; i++)
			l5[i] -= tk2->y[i] * Hconf5 / tk2->x[i] / tk2->x[i];

		gsl_spline_free(tk1);
		gsl_spline_free(tk2);

		loadTransferFunctions(class_background, class_perturbs, tk1, tk2, NULL, sim.boxsize, (1. / (0.995 * a)) - 1., cosmo.h);

		for (i = 0; i < tk1->size; i++)
			l2[i] = -tk1->y[i];

		gsl_spline_free(tk1);
		gsl_spline_free(tk2);

		loadTransferFunctions(class_background, class_perturbs, tk1, tk2, "tot", sim.boxsize, (1. / (0.995 * a)) - 1., cosmo.h);

		for (i = 0; i < tk1->size; i++)
			l2[i] -= tk2->y[i] * Hconf2 / tk2->x[i] / tk2->x[i];

		gsl_spline_free(tk1);
		gsl_spline_free(tk2);

		loadTransferFunctions(class_background, class_perturbs, tk1, tk2, NULL, sim.boxsize, (1. / (0.99 * a)) - 1., cosmo.h);

		for (i = 0; i < tk1->size; i++)
			l1[i] = -tk1->y[i];

		gsl_spline_free(tk1);
		gsl_spline_free(tk2);

		loadTransferFunctions(class_background, class_perturbs, tk1, tk2, "tot", sim.boxsize, (1. / (0.99 * a)) - 1., cosmo.h);

		for (i = 0; i < tk1->size; i++)
			l1[i] -= tk2->y[i] * Hconf1 / tk2->x[i] / tk2->x[i];

		gsl_spline_free(tk1);
		gsl_spline_free(tk2);

		loadTransferFunctions(class_background, class_perturbs, tk1, tk2, "tot", sim.boxsize, (1. / a) - 1., cosmo.h);

		for (i = 0; i < tk1->size; i++)
			l3[i] -= tk2->y[i] * Hconf3 / tk2->x[i] / tk2->x[i];

		Hconf1 = (8.0802 * Hconf4 - 7.9202 * Hconf2 - 1.0201 * Hconf5 + 0.9801 * Hconf1) / 12.;

		for (i = 0; i < tk1->size; i++)
			chi[i] += 10000. * Hconf3 * (Hconf1 * (8. * l4[i] - 8. * l2[i] - l5[i] + l1[i]) + Hconf3 * (16. * l4[i] + 16. * l2[i] - 30. * l3[i] - l5[i] - l1[i])) * M_PI * sqrt(Pk_primordial(tk1->x[i] * cosmo.h / sim.boxsize, ic) / tk1->x[i]) / tk1->x[i] / tk1->x[i] / tk1->x[i];

		free(l1);
		free(l2);
		free(l3);
		free(l4);
		free(l5);
	}

	gsl_spline_free(tk2);
	tk2 = gsl_spline_alloc(gsl_interp_cspline, tk1->size);
	gsl_spline_init(tk2, tk1->x, chi, tk1->size);

	generateRealization(scalarFT, 0., tk2, (unsigned int) ic.seed, ic.flags & ICFLAG_KSPHERE, 0);

	gsl_spline_free(tk1);
	gsl_spline_free(tk2);
	free(chi);
}
#endif

#endif


//////////////////////////
// ic_prevolution.hpp
//////////////////////////
// 
// initial condition generator for gevolution consistently
// generating the full phase space for relativistic particles
// [see C.-P. Ma and E. Bertschinger, Astrophys. J. 429, 22 (1994)]
//
// Author: Julian Adamek (Université de Genève & Observatoire de Paris & Queen Mary University of London)
//
// Last modified: November 2019
//
//////////////////////////

#ifndef IC_PREVOLUTION_HEADER
#define IC_PREVOLUTION_HEADER

#ifndef HAVE_CLASS
#error This version of the prevolution IC generator requires CLASS to be linked as a library! Recompile with -DHAVE_CLASS -lclass
#endif

#include "tools.hpp"

using namespace std;
using namespace LATfield2;


//////////////////////////
// generateIC_prevolution
//////////////////////////
// Description:
//   initial condition generator consistently generating the full phase
//   space for relativistic particles
//
// Arguments: 
//   sim            simulation metadata structure
//   ic             settings for IC generation
//   cosmo          cosmological parameter structure
//   fourpiG        4 pi G (in code units)
//   a              reference to scale factor
//   tau            reference to conformal coordinate time
//   dtau           time step
//   dtau_old       previous time step (will be passed back)
//   pcls_cdm       pointer to (uninitialized) particle handler for CDM
//   pcls_b         pointer to (uninitialized) particle handler for baryons
//   pcls_ncdm      array of (uninitialized) particle handlers for
//                  non-cold DM (may be set to NULL)
//   maxvel         array that will contain the maximum q/m/a (max. velocity)
//   phi            pointer to allocated field
//   chi            pointer to allocated field
//   Bi             pointer to allocated field
//   source         pointer to allocated field
//   Sij            pointer to allocated field
//   scalarFT       pointer to allocated field
//   BiFT           pointer to allocated field
//   SijFT          pointer to allocated field
//   plan_phi       pointer to FFT planner
//   plan_chi       pointer to FFT planner
//   plan_Bi        pointer to FFT planner
//   plan_source    pointer to FFT planner
//   plan_Sij       pointer to FFT planner
//   params         pointer to array of precision settings for CLASS (can be NULL)
//   numparam       number of precision settings for CLASS (can be 0)
//
// Returns:
// 
//////////////////////////

void generateIC_prevolution(metadata & sim, icsettings & ic, cosmology & cosmo, const double fourpiG, double & a, double & tau, double & dtau, double & dtau_old, Particles<part_simple,part_simple_info,part_simple_dataType> * pcls_cdm, Particles<part_simple,part_simple_info,part_simple_dataType> * pcls_b, Particles<part_simple,part_simple_info,part_simple_dataType> * pcls_ncdm, double * maxvel, Field<Real> * phi, Field<Real> * chi, Field<Real> * Bi, Field<Real> * source, Field<Real> * Sij, Field<Cplx> * scalarFT, Field<Cplx> * BiFT, Field<Cplx> * SijFT, PlanFFT<Cplx> * plan_phi, PlanFFT<Cplx> * plan_chi, PlanFFT<Cplx> * plan_Bi, PlanFFT<Cplx> * plan_source, PlanFFT<Cplx> * plan_Sij, parameter * params, int & numparam)
{
#ifdef HAVE_CLASS
	int i, j, p;
	float * pcldata = NULL;
	string h5filename;
	char filename[2*PARAM_MAX_LENGTH+24];
	gsl_spline * phispline = NULL;
	gsl_spline * phispline2 = NULL;
	gsl_spline * chispline = NULL;
	gsl_spline * tk_d1 = NULL;
	gsl_spline * tk_d2 = NULL;
	gsl_spline * tk_t1 = NULL;
	gsl_spline * tk_t2 = NULL;
	double * temp1 = NULL;
	double * temp2 = NULL;
	Site x(phi->lattice());
	double max_displacement;
	double rescale;
	double mean_q;
	double dx = 1.0 / (double) sim.numpts;
	long tmp;
	part_simple_info pcls_cdm_info;
	part_simple_dataType pcls_cdm_dataType;
	part_simple_info pcls_b_info;
	part_simple_dataType pcls_b_dataType;
	part_simple_info pcls_ncdm_info[MAX_PCL_SPECIES];
	part_simple_dataType pcls_ncdm_dataType;
	Real boxSize[3] = {1.,1.,1.};
	char ncdm_name[8];
	Field<Real> * ic_fields[2];
	int cycle = 0;
    int relax_cycles = 0;

	Real * kbin;
	Real * power;
	Real * kscatter;
	Real * pscatter;
	int * occupation;
    long numpts3d = (long) sim.numpts * (long) sim.numpts * (long) sim.numpts;

  	background class_background;
  	thermo class_thermo;
  	perturbs class_perturbs;
	
	ic_fields[0] = phi;
	ic_fields[1] = chi;

	h5filename.reserve(2*PARAM_MAX_LENGTH);
	h5filename.assign(sim.output_path);
	h5filename += sim.basename_generic;

	kbin = (Real *) malloc(sim.numbins * sizeof(Real));
	power = (Real *) malloc(sim.numbins * sizeof(Real));
	kscatter = (Real *) malloc(sim.numbins * sizeof(Real));
	pscatter = (Real *) malloc(sim.numbins * sizeof(Real));
	occupation = (int *) malloc(sim.numbins * sizeof(int));

	a = 1. / (1. + ic.z_ic);
	tau = particleHorizon(a, fourpiG, cosmo);

	COUT << " initial particle horizon tau = " << tau * sim.numpts << " lattice units." << endl;

	initializeCLASSstructures(sim, ic, cosmo, class_background, class_thermo, class_perturbs, params, numparam);

	loadTransferFunctions(class_background, class_perturbs, tk_d1, tk_d2, NULL, sim.boxsize, ic.z_ic, cosmo.h);

	temp1 = (double *) malloc(tk_d1->size * sizeof(double));
	temp2 = (double *) malloc(tk_d1->size * sizeof(double));

	for (i = 0; i < tk_d1->size; i++)
	{
		temp1[i] = -tk_d1->y[i] * M_PI * sqrt(Pk_primordial(tk_d1->x[i] * cosmo.h / sim.boxsize, ic) / tk_d1->x[i]) / tk_d1->x[i];
		temp2[i] = (tk_d2->y[i] - tk_d1->y[i]) * M_PI * sqrt(Pk_primordial(tk_d1->x[i] * cosmo.h / sim.boxsize, ic) / tk_d1->x[i]) / tk_d1->x[i];
	}
	phispline = gsl_spline_alloc(gsl_interp_cspline, tk_d1->size);
	chispline = gsl_spline_alloc(gsl_interp_cspline, tk_d1->size);
	gsl_spline_init(phispline, tk_d1->x, temp1, tk_d1->size);
	gsl_spline_init(chispline, tk_d1->x, temp2, tk_d1->size);
	gsl_spline_free(tk_d1);
	gsl_spline_free(tk_d2);
	
	generateCICKernel(*source);

	for (p = 0; p < cosmo.num_ncdm; p++)
	{
		if (ic.numtile[((sim.baryon_flag == 1) ? 2 : 1)+p] < 1) continue;

		loadHomogeneousTemplate(ic.pclfile[((sim.baryon_flag == 1) ? 2 : 1)+p], sim.numpcl[((sim.baryon_flag == 1) ? 2 : 1)+p], pcldata);
	
		if (pcldata == NULL)
		{
			COUT << " error: particle data was empty!" << endl;
			parallel.abortForce();
		}
		
		sprintf(ncdm_name, "ncdm[%d]", p);
		loadTransferFunctions(class_background, class_perturbs, tk_d1, tk_t1, ncdm_name, sim.boxsize, ic.z_ic, cosmo.h);
	
		if (tk_d1 == NULL || tk_t1 == NULL)
		{
			COUT << " error: ncdm transfer function was empty! (species " << p << ")" << endl;
			parallel.abortForce();
		}

		rescale = pow(cosmo.Omega_g * cosmo.h * cosmo.h / C_PLANCK_LAW, 0.25) * cosmo.T_ncdm[p] * C_BOLTZMANN_CST / cosmo.m_ncdm[p];
		
		for (i = 0; i < tk_d1->size; i++)
		{
			temp1[i] = -3. * phispline->y[i] - 0.75 * tk_d1->y[i] * M_PI * sqrt(Pk_primordial(tk_d1->x[i] * cosmo.h / sim.boxsize, ic) / tk_d1->x[i]) / tk_d1->x[i];
			temp2[i] = -0.25 * tk_d1->y[i] * M_PI * sqrt(Pk_primordial(tk_d1->x[i] * cosmo.h / sim.boxsize, ic) * tk_d1->x[i]);
		}
		gsl_spline_free(tk_d1);
		tk_d1 = gsl_spline_alloc(gsl_interp_cspline, phispline->size);
		tk_d2 = gsl_spline_alloc(gsl_interp_cspline, phispline->size);
		gsl_spline_init(tk_d1, phispline->x, temp1, phispline->size);
		gsl_spline_init(tk_d2, phispline->x, temp2, phispline->size);

		for (i = 0; i < tk_t1->size; i++)
			temp2[i] = -4.2018325 * rescale * tk_t1->y[i] * M_PI * sqrt(Pk_primordial(tk_t1->x[i] * cosmo.h / sim.boxsize, ic) / tk_t1->x[i]) / tk_t1->x[i];
		gsl_spline_free(tk_t1);
		tk_t1 = gsl_spline_alloc(gsl_interp_cspline, phispline->size);
		gsl_spline_init(tk_t1, phispline->x, temp2, phispline->size);
		
		plan_source->execute(FFT_FORWARD);
		generateDisplacementField(*scalarFT, 0., tk_d1, (unsigned int) ic.seed, ic.flags & ICFLAG_KSPHERE);
		plan_chi->execute(FFT_BACKWARD);
		chi->updateHalo();
		gsl_spline_free(tk_d1);
			
		plan_source->execute(FFT_FORWARD);
		generateDisplacementField(*scalarFT, 0., tk_t1, (unsigned int) ic.seed, ic.flags & ICFLAG_KSPHERE, 0);
		plan_phi->execute(FFT_BACKWARD);
		phi->updateHalo();
		gsl_spline_free(tk_t1);
		
		strcpy(pcls_ncdm_info[p].type_name, "part_simple");
		pcls_ncdm_info[p].mass = cosmo.Omega_ncdm[p] / (Real) (sim.numpcl[((sim.baryon_flag == 1) ? 2 : 1)+p]*(long)ic.numtile[((sim.baryon_flag == 1) ? 2 : 1)+p]*(long)ic.numtile[((sim.baryon_flag == 1) ? 2 : 1)+p]*(long)ic.numtile[((sim.baryon_flag == 1) ? 2 : 1)+p]);
		pcls_ncdm_info[p].relativistic = true;
		
		pcls_ncdm[p].initialize(pcls_ncdm_info[p], pcls_ncdm_dataType, &(phi->lattice()), boxSize);
		
		initializeParticlePositions(sim.numpcl[((sim.baryon_flag == 1) ? 2 : 1)+p], pcldata, ic.numtile[((sim.baryon_flag == 1) ? 2 : 1)+p], pcls_ncdm[p]);
		i = MAX;
		pcls_ncdm[p].moveParticles(displace_pcls_ic_basic, 1., &chi, 1, NULL, &max_displacement, &i, 1);
		
		sim.numpcl[((sim.baryon_flag == 1) ? 2 : 1)+p] *= (long) ic.numtile[((sim.baryon_flag == 1) ? 2 : 1)+p] * (long) ic.numtile[((sim.baryon_flag == 1) ? 2 : 1)+p] * (long) ic.numtile[((sim.baryon_flag == 1) ? 2 : 1)+p];
	
		COUT << " " << sim.numpcl[((sim.baryon_flag == 1) ? 2 : 1)+p] << " ncdm particles initialized for species " << p+1 << ": maximum displacement = " << max_displacement * sim.numpts << " lattice units." << endl;
		
		free(pcldata);
		
		pcls_ncdm[p].updateVel(initialize_q_ic_basic, 1., &phi, 1);

		plan_source->execute(FFT_FORWARD);
		generateDisplacementField(*scalarFT, 0., tk_d2, (unsigned int) ic.seed, ic.flags & ICFLAG_KSPHERE);
		plan_chi->execute(FFT_BACKWARD);
		chi->updateHalo();
		gsl_spline_free(tk_d2);
	
		sprintf(filename, "_dT_ncdm%d.h5", p);
		chi->saveHDF5(h5filename + filename);

		mean_q = applyMomentumDistribution(pcls_ncdm+p, (unsigned int) (ic.seed + p), rescale, chi);
		parallel.sum(mean_q);
		COUT << " species " << p+1 << " Fermi-Dirac distribution had mean q/m = " << mean_q / sim.numpcl[((sim.baryon_flag == 1) ? 2 : 1)+p] << endl;
	}

	dtau = ic.Cf * dx;
	    
	dtau_old = 0.;

	COUT << " evolving particle phase space..." << endl;

	do
	{   
	    if (cycle > 0 && (1. / a < ic.z_relax + 1. || 1. / a <= sim.z_in + 1.))
	    {
			if (relax_cycles == 0) // initialize cold species
			{
				loadHomogeneousTemplate(ic.pclfile[0], sim.numpcl[0], pcldata);
	
				if (pcldata == NULL)
				{
					COUT << " error: particle data was empty!" << endl;
					parallel.abortForce();
				}
	
				if (ic.flags & ICFLAG_CORRECT_DISPLACEMENT)
					generateCICKernel(*source, sim.numpcl[0], pcldata, ic.numtile[0]);
				else
					generateCICKernel(*source);
	
				plan_source->execute(FFT_FORWARD);

				loadTransferFunctions(class_background, class_perturbs, tk_d1, tk_d2, NULL, sim.boxsize, (1. / a) - 1., cosmo.h);

				for (i = 0; i < tk_d1->size; i++)
					temp1[i] = -tk_d1->y[i] * M_PI * sqrt(Pk_primordial(tk_d1->x[i] * cosmo.h / sim.boxsize, ic) / tk_d1->x[i]) / tk_d1->x[i];

				phispline2 = gsl_spline_alloc(gsl_interp_cspline, tk_d1->size);
				gsl_spline_init(phispline2, tk_d1->x, temp1, tk_d1->size);
				gsl_spline_free(tk_d1);
				gsl_spline_free(tk_d2);

				loadTransferFunctions(class_background, class_perturbs, tk_d1, tk_t1, "cdm", sim.boxsize, (1. / a) - 1., cosmo.h);
		
				if (sim.baryon_flag > 0)
				{
					loadTransferFunctions(class_background, class_perturbs, tk_d2, tk_t2, "b", sim.boxsize, (1. / a) - 1., cosmo.h);

					if (tk_d2->size != tk_d1->size)
					{
						COUT << " error: baryon transfer function line number mismatch!" << endl;
						parallel.abortForce();
					}
				}
		
				if (sim.baryon_flag == 2)
				{
					for (i = 0; i < tk_d1->size; i++)
					{
						temp1[i] = -3. * phispline2->y[i] - ((cosmo.Omega_cdm * tk_d1->y[i] + cosmo.Omega_b * tk_d2->y[i]) / (cosmo.Omega_cdm + cosmo.Omega_b)) * M_PI * sqrt(Pk_primordial(tk_d1->x[i] * cosmo.h / sim.boxsize, ic) / tk_d1->x[i]) / tk_d1->x[i];
						temp2[i] = -a * ((cosmo.Omega_cdm * tk_t1->y[i] + cosmo.Omega_b * tk_t2->y[i]) / (cosmo.Omega_cdm + cosmo.Omega_b)) * M_PI * sqrt(Pk_primordial(tk_d1->x[i] * cosmo.h / sim.boxsize, ic) / tk_d1->x[i]) / tk_d1->x[i];
					}
					gsl_spline_free(tk_d1);
					gsl_spline_free(tk_t1);
					tk_d1 = gsl_spline_alloc(gsl_interp_cspline, tk_d2->size);
					tk_t1 = gsl_spline_alloc(gsl_interp_cspline, tk_d2->size);
					gsl_spline_init(tk_d1, tk_d2->x, temp1, tk_d2->size);
					gsl_spline_init(tk_t1, tk_d2->x, temp2, tk_d2->size);
					gsl_spline_free(tk_d2);
					gsl_spline_free(tk_t2);
				}
		
				if (sim.baryon_flag == 3)
				{
					if (8. * cosmo.Omega_b / (cosmo.Omega_cdm + cosmo.Omega_b) > 1.)
					{
						for (i = 0; i < tk_d1->size; i++)
						{
							temp1[i] = -3. * phispline2->y[i] - ((8. * cosmo.Omega_cdm * tk_d1->y[i] + (7. * cosmo.Omega_b - cosmo.Omega_cdm) * tk_d2->y[i]) / (cosmo.Omega_cdm + cosmo.Omega_b) / 7.) * M_PI * sqrt(Pk_primordial(tk_d1->x[i] * cosmo.h / sim.boxsize, ic) / tk_d1->x[i]) / tk_d1->x[i];
							temp2[i] = -a * ((8. * cosmo.Omega_cdm * tk_t1->y[i] + (7. * cosmo.Omega_b - cosmo.Omega_cdm) * tk_t2->y[i]) / (cosmo.Omega_cdm + cosmo.Omega_b) / 7.) * M_PI * sqrt(Pk_primordial(tk_d1->x[i] * cosmo.h / sim.boxsize, ic) / tk_d1->x[i]) / tk_d1->x[i];
						}
						gsl_spline_free(tk_d1);
						gsl_spline_free(tk_t1);
						tk_d1 = gsl_spline_alloc(gsl_interp_cspline, tk_d2->size);
						tk_t1 = gsl_spline_alloc(gsl_interp_cspline, tk_d2->size);
						gsl_spline_init(tk_d1, tk_d2->x, temp1, tk_d2->size);
						gsl_spline_init(tk_t1, tk_d2->x, temp2, tk_d2->size);
					}
					else
					{
						for (i = 0; i < tk_d1->size; i++)
						{
							temp1[i] = -3. * phispline2->y[i] - (((cosmo.Omega_cdm - 7. * cosmo.Omega_b) * tk_d1->y[i] + 8. * cosmo.Omega_b * tk_d2->y[i]) / (cosmo.Omega_cdm + cosmo.Omega_b)) * M_PI * sqrt(Pk_primordial(tk_d1->x[i] * cosmo.h / sim.boxsize, ic) / tk_d1->x[i]) / tk_d1->x[i];
							temp2[i] = -a * (((cosmo.Omega_cdm - 7. * cosmo.Omega_b) * tk_t1->y[i] + 8. * cosmo.Omega_b * tk_t2->y[i]) / (cosmo.Omega_cdm + cosmo.Omega_b)) * M_PI * sqrt(Pk_primordial(tk_d1->x[i] * cosmo.h / sim.boxsize, ic) / tk_d1->x[i]) / tk_d1->x[i];
						}
						gsl_spline_free(tk_d2);
						gsl_spline_free(tk_t2);
						tk_d2 = gsl_spline_alloc(gsl_interp_cspline, tk_d1->size);
						tk_t2 = gsl_spline_alloc(gsl_interp_cspline, tk_d1->size);
						gsl_spline_init(tk_d2, tk_d2->x, temp1, tk_d1->size);
						gsl_spline_init(tk_t2, tk_d2->x, temp2, tk_d1->size);
					}
				}
		
				if (sim.baryon_flag == 1 || (sim.baryon_flag == 3 && 8. * cosmo.Omega_b / (cosmo.Omega_cdm + cosmo.Omega_b) > 1.))
				{
					for (i = 0; i < tk_d1->size; i++)
					{
						temp1[i] = -3. * phispline2->y[i] - tk_d2->y[i] * M_PI * sqrt(Pk_primordial(tk_d2->x[i] * cosmo.h / sim.boxsize, ic) / tk_d2->x[i]) / tk_d2->x[i];
						temp2[i] = -a * tk_t2->y[i] * M_PI * sqrt(Pk_primordial(tk_d2->x[i] * cosmo.h / sim.boxsize, ic) / tk_d2->x[i]) / tk_d2->x[i];
					}
					gsl_spline_free(tk_d2);
					gsl_spline_free(tk_t2);
					tk_d2 = gsl_spline_alloc(gsl_interp_cspline, tk_d1->size);
					tk_t2 = gsl_spline_alloc(gsl_interp_cspline, tk_d1->size);
					gsl_spline_init(tk_d2, tk_d1->x, temp1, tk_d1->size);
					gsl_spline_init(tk_t2, tk_d1->x, temp2, tk_d1->size);
				}
		
				if (sim.baryon_flag < 2 || (sim.baryon_flag == 3 && 8. * cosmo.Omega_b / (cosmo.Omega_cdm + cosmo.Omega_b) <= 1.))
				{
					for (i = 0; i < tk_d1->size; i++)
					{
						temp1[i] = -3. * phispline2->y[i] - tk_d1->y[i] * M_PI * sqrt(Pk_primordial(tk_d1->x[i] * cosmo.h / sim.boxsize, ic) / tk_d1->x[i]) / tk_d1->x[i];
						temp2[i] = -a * tk_t1->y[i] * M_PI * sqrt(Pk_primordial(tk_d1->x[i] * cosmo.h / sim.boxsize, ic) / tk_d1->x[i]) / tk_d1->x[i];
					}
					gsl_spline_free(tk_d1);
					gsl_spline_free(tk_t1);
					tk_d1 = gsl_spline_alloc(gsl_interp_cspline, phispline2->size);
					tk_t1 = gsl_spline_alloc(gsl_interp_cspline, phispline2->size);
					gsl_spline_init(tk_d1, phispline2->x, temp1, phispline2->size);
					gsl_spline_init(tk_t1, phispline2->x, temp2, phispline2->size);
				}

				gsl_spline_free(phispline2);
		
				if ((sim.baryon_flag == 1 && !(ic.flags & ICFLAG_CORRECT_DISPLACEMENT)) || sim.baryon_flag == 3)
				{
					generateDisplacementField(*scalarFT, 0., tk_d2, (unsigned int) ic.seed, ic.flags & ICFLAG_KSPHERE);
					gsl_spline_free(tk_d2);
					plan_chi->execute(FFT_BACKWARD);
					chi->updateHalo();
					plan_source->execute(FFT_FORWARD);
				}
		
				generateDisplacementField(*scalarFT, 0., tk_d1, (unsigned int) ic.seed, ic.flags & ICFLAG_KSPHERE);
				gsl_spline_free(tk_d1);
		
				plan_phi->execute(FFT_BACKWARD);
				phi->updateHalo();
	
				strcpy(pcls_cdm_info.type_name, "part_simple");
				if (sim.baryon_flag == 1)
					pcls_cdm_info.mass = cosmo.Omega_cdm / (Real) (sim.numpcl[0]*(long)ic.numtile[0]*(long)ic.numtile[0]*(long)ic.numtile[0]);
				else
					pcls_cdm_info.mass = (cosmo.Omega_cdm + cosmo.Omega_b) / (Real) (sim.numpcl[0]*(long)ic.numtile[0]*(long)ic.numtile[0]*(long)ic.numtile[0]);
				pcls_cdm_info.relativistic = false;
	
				pcls_cdm->initialize(pcls_cdm_info, pcls_cdm_dataType, &(phi->lattice()), boxSize);
	
				initializeParticlePositions(sim.numpcl[0], pcldata, ic.numtile[0], *pcls_cdm);
				i = MAX;
				if (sim.baryon_flag == 3)
					pcls_cdm->moveParticles(displace_pcls_ic_basic, 1., ic_fields, 2, NULL, &max_displacement, &i, 1);
				else
					pcls_cdm->moveParticles(displace_pcls_ic_basic, 1., &phi, 1, NULL, &max_displacement, &i, 1);
	
				sim.numpcl[0] *= (long) ic.numtile[0] * (long) ic.numtile[0] * (long) ic.numtile[0];
	
				COUT << " " << sim.numpcl[0] << " cdm particles initialized: maximum displacement = " << max_displacement * sim.numpts << " lattice units." << endl;
	
				free(pcldata);
	
				if (sim.baryon_flag == 1)
				{
					loadHomogeneousTemplate(ic.pclfile[1], sim.numpcl[1], pcldata);
	
					if (pcldata == NULL)
					{
						COUT << " error: particle data was empty!" << endl;
						parallel.abortForce();
					}
		
					if (ic.flags & ICFLAG_CORRECT_DISPLACEMENT)
					{
						generateCICKernel(*chi, sim.numpcl[1], pcldata, ic.numtile[1]);
						plan_chi->execute(FFT_FORWARD);
						generateDisplacementField(*scalarFT, 0., tk_d2, (unsigned int) ic.seed, ic.flags & ICFLAG_KSPHERE);
						gsl_spline_free(tk_d2);
						plan_chi->execute(FFT_BACKWARD);
						chi->updateHalo();
					}
		
					strcpy(pcls_b_info.type_name, "part_simple");
					pcls_b_info.mass = cosmo.Omega_b / (Real) (sim.numpcl[1]*(long)ic.numtile[1]*(long)ic.numtile[1]*(long)ic.numtile[1]);
					pcls_b_info.relativistic = false;
	
					pcls_b->initialize(pcls_b_info, pcls_b_dataType, &(phi->lattice()), boxSize);
	
					initializeParticlePositions(sim.numpcl[1], pcldata, ic.numtile[1], *pcls_b);
					i = MAX;
					pcls_b->moveParticles(displace_pcls_ic_basic, 1., &chi, 1, NULL, &max_displacement, &i, 1);
	
					sim.numpcl[1] *= (long) ic.numtile[1] * (long) ic.numtile[1] * (long) ic.numtile[1];
	
					COUT << " " << sim.numpcl[1] << " baryon particles initialized: maximum displacement = " << max_displacement * sim.numpts << " lattice units." << endl;
	
					free(pcldata);
				}
	
				if (ic.flags & ICFLAG_CORRECT_DISPLACEMENT)
					generateCICKernel(*source);
		
				plan_source->execute(FFT_FORWARD);
		
				if (sim.baryon_flag == 1 || sim.baryon_flag == 3)
				{
					generateDisplacementField(*scalarFT, 0., tk_t2, (unsigned int) ic.seed, ic.flags & ICFLAG_KSPHERE, 0);
					plan_chi->execute(FFT_BACKWARD);
					chi->updateHalo();
					gsl_spline_free(tk_t2);
					plan_source->execute(FFT_FORWARD);
				}
		
				generateDisplacementField(*scalarFT, 0., tk_t1, (unsigned int) ic.seed, ic.flags & ICFLAG_KSPHERE, 0);
				plan_phi->execute(FFT_BACKWARD);
				phi->updateHalo();
				gsl_spline_free(tk_t1);			
		
				if (sim.baryon_flag == 3)
					maxvel[0] = pcls_cdm->updateVel(initialize_q_ic_basic, 1., ic_fields, 2) / a;
				else
					maxvel[0] = pcls_cdm->updateVel(initialize_q_ic_basic, 1., &phi, 1) / a;
		
				if (sim.baryon_flag == 1)
					maxvel[1] = pcls_b->updateVel(initialize_q_ic_basic, 1., &chi, 1) / a;
	
				if (sim.baryon_flag > 1) sim.baryon_flag = 0;

				generateRealization(*scalarFT, 0., phispline, (unsigned int) ic.seed, ic.flags & ICFLAG_KSPHERE, 0);
				plan_phi->execute(FFT_BACKWARD);
				phi->updateHalo();

				if (1. / a <= sim.z_in + 1.)
				{
					generateRealization(*scalarFT, 0., chispline, (unsigned int) ic.seed, ic.flags & ICFLAG_KSPHERE, 0);
					plan_chi->execute(FFT_BACKWARD);
					chi->updateHalo();

					for (p = 0; p < cosmo.num_ncdm; p++)
					{
						if (ic.numtile[1+sim.baryon_flag+p] < 1) maxvel[1+sim.baryon_flag+p] = 0;
						else maxvel[1+sim.baryon_flag+p] = pcls_ncdm[p].updateVel(update_q, dtau_old / 2., ic_fields, 2, &a);
					}

					dtau_old = 0.;
					break;
				}
			}

			for (i = 0; i < phispline->size; i++)
				temp1[i] = 0.;

			if (cosmo.Omega_g > 0)
			{
				loadTransferFunctions(class_background, class_perturbs, tk_d1, tk_t1, "g", sim.boxsize, ic.z_ic, cosmo.h);

				for (i = 0; i < tk_d1->size; i++)
					temp1[i] -= tk_d1->y[i] * cosmo.Omega_g * M_PI * sqrt(Pk_primordial(tk_d1->x[i] * cosmo.h / sim.boxsize, ic) / tk_d1->x[i]) / tk_d1->x[i] / a;

				gsl_spline_free(tk_d1);
				gsl_spline_free(tk_t1);
			}

			if (cosmo.Omega_ur > 0)
			{
				loadTransferFunctions(class_background, class_perturbs, tk_d1, tk_t1, "ur", sim.boxsize, ic.z_ic, cosmo.h);

				for (i = 0; i < tk_d1->size; i++)
					temp1[i] -= tk_d1->y[i] * cosmo.Omega_ur * M_PI * sqrt(Pk_primordial(tk_d1->x[i] * cosmo.h / sim.boxsize, ic) / tk_d1->x[i]) / tk_d1->x[i] / a;

				gsl_spline_free(tk_d1);
				gsl_spline_free(tk_t1);
			}
			
			for (p = 0; p < cosmo.num_ncdm; p++)
			{
				sprintf(ncdm_name, "ncdm[%d]", p);
				loadTransferFunctions(class_background, class_perturbs, tk_d1, tk_t1, ncdm_name, sim.boxsize, ic.z_ic, cosmo.h);

				rescale = bg_ncdm(a, cosmo, p);

				for (i = 0; i < tk_d1->size; i++)
					temp1[i] -= tk_d1->y[i] * rescale * M_PI * sqrt(Pk_primordial(tk_d1->x[i] * cosmo.h / sim.boxsize, ic) / tk_d1->x[i]) / tk_d1->x[i];

				gsl_spline_free(tk_d1);
				gsl_spline_free(tk_t1);
			}

			if (cosmo.num_ncdm > 0 || cosmo.Omega_rad > 0)
			{
				tk_d1 = gsl_spline_alloc(gsl_interp_cspline, phispline->size);
				gsl_spline_init(tk_d1, phispline->x, temp1, phispline->size);
				
				generateRealization(*scalarFT, 0., tk_d1, (unsigned int) ic.seed, ic.flags & ICFLAG_KSPHERE);
				gsl_spline_free(tk_d1);
			}

			projection_init(source);
			if (cosmo.num_ncdm > 0 || cosmo.Omega_rad > 0)
				plan_source->execute(FFT_BACKWARD);
	        	projection_T00_project(pcls_cdm, source, a, phi);
			if (sim.baryon_flag)
				projection_T00_project(pcls_b, source, a, phi);
			projection_T00_comm(source);
	        
	        prepareFTsource<Real>(*phi, *chi, *source, cosmo.Omega_cdm + cosmo.Omega_b, *source, 3. * Hconf(a, fourpiG, cosmo) * dx * dx / dtau_old, fourpiG * dx * dx / a, 3. * Hconf(a, fourpiG, cosmo) * Hconf(a, fourpiG, cosmo) * dx * dx);
	        plan_source->execute(FFT_FORWARD);
	        solveModifiedPoissonFT(*scalarFT, *scalarFT, 1. / (dx * dx), 3. * Hconf(a, fourpiG, cosmo) / dtau_old);
	        plan_chi->execute(FFT_BACKWARD);     // nonlinear phi temporarily stored in chi
            
			relax_cycles++;
	    }

		loadTransferFunctions(class_background, class_perturbs, tk_d1, tk_d2, NULL, sim.boxsize, (1. / a) - 1., cosmo.h);

		for (i = 0; i < tk_d1->size; i++)
		{
			temp1[i] = -tk_d1->y[i] * M_PI * sqrt(Pk_primordial(tk_d1->x[i] * cosmo.h / sim.boxsize, ic) / tk_d1->x[i]) / tk_d1->x[i];
			temp2[i] = (tk_d2->y[i] - tk_d1->y[i]) * M_PI * sqrt(Pk_primordial(tk_d1->x[i] * cosmo.h / sim.boxsize, ic) / tk_d1->x[i]) / tk_d1->x[i];
		}
		gsl_spline_free(phispline);
		gsl_spline_free(chispline);
		phispline = gsl_spline_alloc(gsl_interp_cspline, tk_d1->size);
		chispline = gsl_spline_alloc(gsl_interp_cspline, tk_d1->size);
		gsl_spline_init(phispline, tk_d1->x, temp1, tk_d1->size);
		gsl_spline_init(chispline, tk_d1->x, temp2, tk_d1->size);
		gsl_spline_free(tk_d1);
		gsl_spline_free(tk_d2);
		
		generateRealization(*scalarFT, 0., phispline, (unsigned int) ic.seed, ic.flags & ICFLAG_KSPHERE, 0);
		plan_phi->execute(FFT_BACKWARD);

		if (relax_cycles == 1)
	    {
	    	plan_phi->execute(FFT_FORWARD);
			extractPowerSpectrum(*scalarFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, false, KTYPE_LINEAR);
			sprintf(filename, "%s%sICtarget_phi.dat", sim.output_path, sim.basename_pk);
			writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI, filename, "power spectrum of phi", a);
			plan_chi->execute(FFT_FORWARD);
			extractPowerSpectrum(*scalarFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, false, KTYPE_LINEAR);
			sprintf(filename, "%s%sICactual_phi.dat", sim.output_path, sim.basename_pk);
			writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI, filename, "power spectrum of phi", a);
		}

	    if (relax_cycles > 0)
	    {	
	        for (x.first(); x.test(); x.next())     // interpolate between linear and nonlinear solution
	            (*phi)(x) = (((1. / (sim.z_in + 1.)) - a) * (*phi)(x) + (a - (1. / (ic.z_relax + 1.))) * (*chi)(x)) / ((1. / (sim.z_in + 1.)) - (1. / (ic.z_relax + 1.)));
		}
		phi->updateHalo();

		if (cycle == 0)
		{
			sprintf(filename, "_ICinit_phi.h5");
			phi->saveHDF5(h5filename + filename);

			plan_phi->execute(FFT_FORWARD);
			extractPowerSpectrum(*scalarFT, kbin, power, kscatter, pscatter, occupation, sim.numbins, false, KTYPE_LINEAR);
			sprintf(filename, "%s%sICinit_phi.dat", sim.output_path, sim.basename_pk);
			writePowerSpectrum(kbin, power, kscatter, pscatter, occupation, sim.numbins, sim.boxsize, (Real) numpts3d * (Real) numpts3d * 2. * M_PI * M_PI, filename, "power spectrum of phi", a);
		}

		generateRealization(*scalarFT, 0., chispline, (unsigned int) ic.seed, ic.flags & ICFLAG_KSPHERE, 0);
		plan_chi->execute(FFT_BACKWARD);

		if (relax_cycles > 0)
	    {
			projection_init(Sij);
	        projection_Tij_project(pcls_cdm, Sij, a, phi);
			if (sim.baryon_flag)
				projection_Tij_project(pcls_b, Sij, a, phi);
	        projection_Tij_comm(Sij);
	
	        prepareFTsource<Real>(*phi, *Sij, *Sij, 2. * fourpiG / a / (double) sim.numpts / (double) sim.numpts);	
	        plan_Sij->execute(FFT_FORWARD);	
            projectFTscalar(*SijFT, *scalarFT);
            plan_source->execute(FFT_BACKWARD);	     // nonlinear chi temporarily stored in source

			for (x.first(); x.test(); x.next())     // interpolate between linear and nonlinear solution
	            (*chi)(x) += (a - (1. / (ic.z_relax + 1.))) * (*source)(x) / ((1. / (sim.z_in + 1.)) - (1. / (ic.z_relax + 1.)));             
	    }
		chi->updateHalo();
	    
		if (relax_cycles > 0)
		{
	    	maxvel[0] = pcls_cdm->updateVel(update_q, (dtau + dtau_old) / 2., ic_fields, 1, &a);
			if (sim.baryon_flag)
				maxvel[1] = pcls_b->updateVel(update_q, (dtau + dtau_old) / 2., ic_fields, 1, &a);
		}

	    for (p = 0; p < cosmo.num_ncdm; p++)
		{
			if (ic.numtile[((sim.baryon_flag == 1) ? 2 : 1)+p] < 1)	maxvel[((sim.baryon_flag == 1) ? 2 : 1)+p] = 0;
			else maxvel[((sim.baryon_flag == 1) ? 2 : 1)+p] = pcls_ncdm[p].updateVel(update_q, (dtau + dtau_old) / 2., ic_fields, 2, &a);
		}

		rungekutta4bg(a, fourpiG, cosmo, 0.5 * dtau);

		if (relax_cycles > 0)
		{
	    	pcls_cdm->moveParticles(update_pos, dtau, ic_fields, 0, &a);
			if (sim.baryon_flag)
				pcls_b->moveParticles(update_pos, dtau, ic_fields, 0, &a);
		}

	    for (p = 0; p < cosmo.num_ncdm; p++)
		{
			if (ic.numtile[((sim.baryon_flag == 1) ? 2 : 1)+p] < 1) continue;
		        pcls_ncdm[p].moveParticles(update_pos, dtau, ic_fields, 2, &a);
		}

		rungekutta4bg(a, fourpiG, cosmo, 0.5 * dtau);
	    tau += dtau;
	    
	    dtau_old = dtau;
	    
		if (relax_cycles > 0 && ic.Cf > sim.Cf)
	        dtau = (((1. / (sim.z_in + 1.)) - a) * ic.Cf * dx + (a - (1. / (ic.z_relax + 1.))) * ((sim.Cf * dx < sim.steplimit / Hconf(a, fourpiG, cosmo)) ? (sim.Cf * dx) : (sim.steplimit / Hconf(a, fourpiG, cosmo)))) / ((1. / (sim.z_in + 1.)) - (1. / (ic.z_relax + 1.)));
		else if (relax_cycles > 0 && ic.Cf * dx > sim.steplimit / Hconf(a, fourpiG, cosmo))
	        dtau = (((1. / (sim.z_in + 1.)) - a) * ic.Cf * dx + (a - (1. / (ic.z_relax + 1.))) * sim.steplimit / Hconf(a, fourpiG, cosmo)) / ((1. / (sim.z_in + 1.)) - (1. / (ic.z_relax + 1.)));
		else
			dtau = ic.Cf * dx;

	    cycle++;
	}
	while (1. / a > sim.z_in + 1. || relax_cycles == 0);

	COUT << " needed " << cycle << " steps and " << relax_cycles << " nonlinear relaxation operations." << endl;
	
	if (sim.gr_flag == 0)
	{
		COUT << " gravity theory = Newton: computing gauge transformation..." << endl;

		loadTransferFunctions(class_background, class_perturbs, tk_d1, tk_t1, "tot", sim.boxsize, (1. / a) - 1., cosmo.h);
		loadTransferFunctions(class_background, class_perturbs, tk_d2, tk_t2, NULL, sim.boxsize, (1. / a) - 1., cosmo.h);

		for (i = 0; i < tk_d2->size; i++)
			temp1[i] = -tk_d2->y[i] * M_PI * sqrt(Pk_primordial(tk_d2->x[i] * cosmo.h / sim.boxsize, ic) / tk_d2->x[i]) / tk_d2->x[i];
		gsl_spline_free(phispline);
		phispline = gsl_spline_alloc(gsl_interp_cspline, tk_d2->size);
		gsl_spline_init(phispline, tk_d2->x, temp1, tk_d2->size);
		gsl_spline_free(tk_d2);
		gsl_spline_free(tk_t2);

		rescale = 3. * Hconf(a, fourpiG, cosmo)  * M_PI;

		for (i = 0; i < tk_d1->size; i++)
			temp1[i] = 3. * phispline->y[i] - rescale * tk_t1->y[i] * sqrt(Pk_primordial(tk_d1->x[i] * cosmo.h / sim.boxsize, ic) / tk_d1->x[i]) / tk_d1->x[i] / tk_d1->x[i] / tk_d1->x[i];

		gsl_spline_free(tk_d1);
		gsl_spline_free(tk_t1);
		tk_d1 = gsl_spline_alloc(gsl_interp_cspline, phispline->size);
		gsl_spline_init(tk_d1, phispline->x, temp1, phispline->size);

		if (ic.flags & ICFLAG_CORRECT_DISPLACEMENT)
		{
			loadHomogeneousTemplate(ic.pclfile[0], tmp, pcldata);
			generateCICKernel(*source, tmp, pcldata, ic.numtile[0]);
			free(pcldata);
		}
		else
			generateCICKernel(*source);
		plan_source->execute(FFT_FORWARD);
		generateDisplacementField(*scalarFT, 0., tk_d1, (unsigned int) ic.seed, ic.flags & ICFLAG_KSPHERE);
		plan_source->execute(FFT_BACKWARD);
		source->updateHalo();

		i = MAX;
		pcls_cdm->moveParticles(displace_pcls_ic_basic, 1., &source, 1, NULL, &max_displacement, &i, 1);
		COUT << " Poisson gauge -> N-body gauge, cdm particles maximum displacement = " << max_displacement * sim.numpts << " lattice units." << endl;
		if (sim.baryon_flag)
		{
			if (ic.flags & ICFLAG_CORRECT_DISPLACEMENT)
			{
				loadHomogeneousTemplate(ic.pclfile[1], tmp, pcldata);
				generateCICKernel(*source, tmp, pcldata, ic.numtile[1]);
				free(pcldata);
				plan_source->execute(FFT_FORWARD);
				generateDisplacementField(*scalarFT, 0., tk_d1, (unsigned int) ic.seed, ic.flags & ICFLAG_KSPHERE);
				plan_source->execute(FFT_BACKWARD);
				source->updateHalo();
			}
			pcls_b->moveParticles(displace_pcls_ic_basic, 1., &source, 1, NULL, &max_displacement, &i, 1);
			COUT << " Poisson gauge -> N-body gauge, baryon particles maximum displacement = " << max_displacement * sim.numpts << " lattice units." << endl;
		}
		if (cosmo.num_ncdm > 0 && ic.flags & ICFLAG_CORRECT_DISPLACEMENT)
		{
			generateCICKernel(*source);
			plan_source->execute(FFT_FORWARD);
			generateDisplacementField(*scalarFT, 0., tk_d1, (unsigned int) ic.seed, ic.flags & ICFLAG_KSPHERE);
			plan_source->execute(FFT_BACKWARD);
			source->updateHalo();
		}
		for (p = 0; p < cosmo.num_ncdm; p++)
		{
			if (ic.numtile[((sim.baryon_flag == 1) ? 2 : 1)+p] < 1) continue;
			pcls_ncdm[p].moveParticles(displace_pcls_ic_basic, 1., &source, 1, NULL, &max_displacement, &i, 1);
			COUT << " Poisson gauge -> N-body gauge, ncdm species " << p+1 << " maximum displacement = " << max_displacement * sim.numpts << " lattice units." << endl;
		}

		loadTransferFunctions(class_background, class_perturbs, tk_d2, tk_t2, NULL, sim.boxsize, (1. / (0.99 * a)) - 1., cosmo.h);

		for (i = 0; i < tk_d2->size; i++)
			temp1[i] = -tk_d2->y[i] * M_PI * sqrt(Pk_primordial(tk_d2->x[i] * cosmo.h / sim.boxsize, ic) / tk_d2->x[i]) / tk_d2->x[i];
		phispline2 = gsl_spline_alloc(gsl_interp_cspline, tk_d2->size);
		gsl_spline_init(phispline2, tk_d2->x, temp1, tk_d2->size);
		gsl_spline_free(tk_d2);
		gsl_spline_free(tk_t2);

		loadTransferFunctions(class_background, class_perturbs, tk_d2, tk_t2, "tot", sim.boxsize, (1. / (0.99 * a)) - 1., cosmo.h);

		rescale = 3. * Hconf(0.99 * a, fourpiG, cosmo) * M_PI;
		mean_q = -99.5 * Hconf(0.995 * a, fourpiG, cosmo);

		for (i = 0; i < tk_d2->size; i++)
			temp1[i] = mean_q * (tk_d1->y[i] - (3. * phispline2->y[i] - rescale * tk_t2->y[i] * sqrt(Pk_primordial(tk_d2->x[i] * cosmo.h / sim.boxsize, ic) / tk_d2->x[i]) / tk_d2->x[i] / tk_d2->x[i] / tk_d2->x[i]));
		gsl_spline_free(tk_d2);
		gsl_spline_free(tk_t2);
		gsl_spline_free(tk_d1);
		gsl_spline_free(phispline2);
		tk_d1 = gsl_spline_alloc(gsl_interp_cspline, phispline->size);
		gsl_spline_init(tk_d1, phispline->x, temp1, phispline->size);

		generateCICKernel(*source);
		plan_source->execute(FFT_FORWARD);
		generateDisplacementField(*scalarFT, 0., tk_d1, (unsigned int) ic.seed, ic.flags & ICFLAG_KSPHERE);
		plan_source->execute(FFT_BACKWARD);
		source->updateHalo();
		gsl_spline_free(tk_d1);

		COUT << " Poisson gauge -> N-body gauge, cdm particles maximum velocity: " << maxvel[0] << " -> ";
		maxvel[0] = pcls_cdm->updateVel(update_q_Newton, 1., &source, 1, &a);
		COUT << maxvel[0] << endl;
		if (sim.baryon_flag)
		{
			COUT << " Poisson gauge -> N-body gauge, baryon particles maximum velocity: " << maxvel[1] << " -> ";
			maxvel[1] = pcls_b->updateVel(update_q_Newton, 1., &source, 1, &a);
			COUT << maxvel[1] << endl;
		}
		for (p = 0; p < cosmo.num_ncdm; p++)
		{
			if (ic.numtile[((sim.baryon_flag == 1) ? 2 : 1)+p] < 1)
			{
				maxvel[p+1+sim.baryon_flag] = 0;
				continue;
			}
			COUT << " Poisson gauge -> N-body gauge, ncdm species " << p+1 << " maximum velocity: " << maxvel[p+1+sim.baryon_flag] << " -> ";
			maxvel[p+1+sim.baryon_flag] = pcls_ncdm[p].updateVel(update_q_Newton, 1., &source, 1, &a);
			COUT << maxvel[p+1+sim.baryon_flag] << endl;
		}
		
		projection_init(source);
		scalarProjectionCIC_project(pcls_cdm, source);
		for (p = 0; p < cosmo.num_ncdm; p++)
		{
			if (ic.numtile[((sim.baryon_flag == 1) ? 2 : 1)+p] < 1) continue;
			scalarProjectionCIC_project(pcls_ncdm+p, source);
		}
		projection_T00_comm(source);
		plan_source->execute(FFT_FORWARD);
		solveModifiedPoissonFT(*scalarFT, *scalarFT, fourpiG / a);
		plan_phi->execute(FFT_BACKWARD);
		phi->updateHalo();
	}


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

	gsl_spline_free(phispline);
	gsl_spline_free(chispline);
	free(temp1);
	free(temp2);

    free(kbin);
	free(power);
	free(kscatter);
	free(pscatter);
	free(occupation);
	
	projection_init(Bi);
	projection_T0i_project(pcls_cdm, Bi, phi);
	if (sim.baryon_flag)
		projection_T0i_project(pcls_b, Bi, phi);
	projection_T0i_comm(Bi);
	plan_Bi->execute(FFT_FORWARD);
	projectFTvector(*BiFT, *BiFT, fourpiG / (double) sim.numpts / (double) sim.numpts);	
	plan_Bi->execute(FFT_BACKWARD);	
	Bi->updateHalo();
#endif
}

#endif


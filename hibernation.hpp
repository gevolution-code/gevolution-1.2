//////////////////////////
// hibernation.hpp
//////////////////////////
// 
// Auxiliary functions for hibernation
//
// Author: Julian Adamek (Université de Genève & Observatoire de Paris & Queen Mary University of London)
//
// Last modified: February 2019
//
//////////////////////////

#ifndef HIBERNATION_HEADER
#define HIBERNATION_HEADER

//////////////////////////
// writeRestartSettings
//////////////////////////
// Description:
//   writes a settings file containing all the relevant metadata for restarting
//   a run from a hibernation point
//
// Arguments: 
//   sim            simulation metadata structure
//   ic             settings for IC generation
//   cosmo          cosmological parameter structure
//   a              scale factor
//   tau            conformal coordinate time
//   dtau           time step
//   cycle          current main control loop cycle count
//   restartcount   restart counter aka number of hibernation point (default -1)
//                  if < 0 no number is associated to the hibernation point
//
// Returns:
// 
//////////////////////////

void writeRestartSettings(metadata & sim, icsettings & ic, cosmology & cosmo, const double a, const double tau, const double dtau, const int cycle, const int restartcount = -1)
{
	char buffer[2*PARAM_MAX_LENGTH+24];
	FILE * outfile;
	int i;
	
	if (!parallel.isRoot()) return;
	
	if (restartcount >= 0)
		sprintf(buffer, "%s%s%03d.ini", sim.restart_path, sim.basename_restart, restartcount);
	else
		sprintf(buffer, "%s%s.ini", sim.restart_path, sim.basename_restart);
	outfile = fopen(buffer, "w");
	if (outfile == NULL)
	{
		cout << " error opening file for restart settings!" << endl;
	}
	else
	{
		fprintf(outfile, "# automatically generated settings for restart after hibernation ");
		if (restartcount < 0)
			fprintf(outfile, "due to wallclock limit ");
		else
			fprintf(outfile, "requested ");
		fprintf(outfile, "at redshift z=%f\n\n", (1./a)-1.);
		
		fprintf(outfile, "# info related to IC generation\n\n");
		fprintf(outfile, "IC generator       = restart\n");
		if (restartcount >= 0)
			sprintf(buffer, "%03d", restartcount);
		else
			buffer[0] = '\0';
		fprintf(outfile, "particle file      = %s%s%s_cdm.h5", sim.restart_path, sim.basename_restart, buffer);
		if (sim.baryon_flag)
			fprintf(outfile, ", %s%s%s_b.h5", sim.restart_path, sim.basename_restart, buffer);
		for (i = 0; i < cosmo.num_ncdm; i++)
		{
			if (sim.numpcl[1+sim.baryon_flag+i] < 1)
				fprintf(outfile, ", /dev/null");
			else
				fprintf(outfile, ", %s%s%s_ncdm%d.h5", sim.restart_path, sim.basename_restart, buffer, i);
		}
		fprintf(outfile, "\n");
		if (sim.gr_flag > 0)
		{
			fprintf(outfile, "metric file        = %s%s%s_phi.h5", sim.restart_path, sim.basename_restart, buffer);
			fprintf(outfile, ", %s%s%s_chi.h5", sim.restart_path, sim.basename_restart, buffer);
			if (sim.vector_flag == VECTOR_PARABOLIC)
				fprintf(outfile, ", %s%s%s_B.h5\n", sim.restart_path, sim.basename_restart, buffer);
			else
#ifdef CHECK_B
				fprintf(outfile, ", %s%s%s_B_check.h5\n", sim.restart_path, sim.basename_restart, buffer);
#else
				fprintf(outfile, "\n");
#endif
		}
		else if (sim.vector_flag == VECTOR_PARABOLIC)
			fprintf(outfile, "metric file        = %s%s%s_B.h5\n", sim.restart_path, sim.basename_restart, buffer);
#ifdef CHECK_B
		else
			fprintf(outfile, "metric file        = %s%s%s_B_check.h5\n", sim.restart_path, sim.basename_restart, buffer);
#endif
			
		fprintf(outfile, "restart redshift   = %.15lf\n", (1./a) - 1.);
		fprintf(outfile, "cycle              = %d\n", cycle);
		fprintf(outfile, "tau                = %.15le\n", tau);
		fprintf(outfile, "dtau               = %.15le\n", dtau);
		fprintf(outfile, "gevolution version = %g\n\n", GEVOLUTION_VERSION);
		fprintf(outfile, "seed               = %d\n", ic.seed);
		if (ic.flags & ICFLAG_KSPHERE)
			fprintf(outfile, "k-domain           = sphere\n");
		else
			fprintf(outfile, "k-domain           = cube\n");
		fprintf(outfile, "\n\n# primordial power spectrum\n\n");
		fprintf(outfile, "k_pivot = %lg\n", ic.k_pivot);
		fprintf(outfile, "A_s     = %lg\n", ic.A_s);
		fprintf(outfile, "n_s     = %lg\n", ic.n_s);
		fprintf(outfile, "\n\n# cosmological parameters\n\n");
		fprintf(outfile, "h         = %lg\n", cosmo.h);
		fprintf(outfile, "Omega_cdm = %.15le\n", cosmo.Omega_cdm);
		fprintf(outfile, "Omega_b   = %.15le\n", cosmo.Omega_b);
		fprintf(outfile, "Omega_g   = %.15le\n", cosmo.Omega_g);
		fprintf(outfile, "Omega_ur  = %.15le\n", cosmo.Omega_ur);
		if (cosmo.Omega_fld > 0.)
		{
			fprintf(outfile, "Omega_fld = %.15le\n", cosmo.Omega_fld);
			fprintf(outfile, "w0_fld    = %lg\n", cosmo.w0_fld);
			fprintf(outfile, "wa_fld    = %lg\n", cosmo.wa_fld);
			if (sim.fluid_flag > 0)
				fprintf(outfile, "cs2_fld   = %lg\n", cosmo.cs2_fld);
		}
		fprintf(outfile, "N_ncdm    = %d\n", cosmo.num_ncdm);
		if (cosmo.num_ncdm > 0)
		{
			fprintf(outfile, "m_cdm     = ");
			for (i = 0; i < cosmo.num_ncdm - 1; i++)
				fprintf(outfile, "%9lf, ", cosmo.m_ncdm[i]);
			fprintf(outfile, "%9lf\n", cosmo.m_ncdm[i]);
			fprintf(outfile, "T_cdm     = ");
			for (i = 0; i < cosmo.num_ncdm - 1; i++)
				fprintf(outfile, "%9lf, ", cosmo.T_ncdm[i]);
			fprintf(outfile, "%9lf\n", cosmo.T_ncdm[i]);
			fprintf(outfile, "deg_cdm   = ");
			for (i = 0; i < cosmo.num_ncdm - 1; i++)
				fprintf(outfile, "%lf, ", cosmo.deg_ncdm[i]);
			fprintf(outfile, "%lg\n", cosmo.deg_ncdm[i]);
		}
		fprintf(outfile, "\n\n# simulation settings\n\n");
		if (sim.baryon_flag > 0)
			fprintf(outfile, "baryon treatment    = sample\n");
		if (sim.radiation_flag > 0)
		{
			fprintf(outfile, "radiation treatment = CLASS\n");
			fprintf(outfile, "switch delta_rad    = %lf\n", sim.z_switch_deltarad);
			if (cosmo.num_ncdm > 0)
			{
				fprintf(outfile, "switch delta_ncdm   = ");
				for (i = 0; i < cosmo.num_ncdm - 1; i++)
					fprintf(outfile, "%lf, ", sim.z_switch_deltancdm[i]);
				fprintf(outfile, "%lf\n", sim.z_switch_deltancdm[i]);
			}
			fprintf(outfile, "switch linear chi   = %lf\n", sim.z_switch_linearchi);
		}
		if (sim.fluid_flag > 0)
			fprintf(outfile, "fluid treatment     = CLASS\n");
		if (sim.gr_flag > 0)
			fprintf(outfile, "gravity theory      = GR\n");
		else
			fprintf(outfile, "gravity theory      = N\n");
		if (sim.vector_flag == VECTOR_ELLIPTIC)
			fprintf(outfile, "vector method       = elliptic\n");
		else
			fprintf(outfile, "vector method       = parabolic\n");
		fprintf(outfile, "\ninitial redshift    = %lg\n", sim.z_in);
		fprintf(outfile, "boxsize             = %lg\n", sim.boxsize);
		fprintf(outfile, "Ngrid               = %d\n", sim.numpts);
		fprintf(outfile, "Courant factor      = %lg\n", sim.Cf);
		fprintf(outfile, "time step limit     = %lg\n", sim.steplimit);
		if (cosmo.num_ncdm > 0)
			fprintf(outfile, "move limit          = %lg\n", sim.movelimit);
		fprintf(outfile, "\n\n# output\n\n");
		fprintf(outfile, "output path         = %s\n", sim.output_path);
		fprintf(outfile, "generic file base   = %s\n", sim.basename_generic);
		fprintf(outfile, "snapshot file base  = %s\n", sim.basename_snapshot);
		fprintf(outfile, "Pk file base        = %s\n", sim.basename_pk);
		fprintf(outfile, "lightcone file base = %s\n", sim.basename_lightcone);
		if (sim.num_snapshot > 0)
		{
			fprintf(outfile, "snapshot redshifts  = ");
			for (i = 0; i < sim.num_snapshot - 1; i++)
				fprintf(outfile, "%lg, ", sim.z_snapshot[i]);
			fprintf(outfile, "%lg\n", sim.z_snapshot[i]);
		}
		if (sim.out_snapshot)
		{
			fprintf(outfile, "snapshot outputs    = ");
			if(sim.out_snapshot & MASK_PHI)
			{
				fprintf(outfile, "phi");
				if (sim.out_snapshot > MASK_CHI)
					fprintf(outfile, ", ");
			}
			if(sim.out_snapshot & MASK_CHI)
			{
				fprintf(outfile, "chi");
				if (sim.out_snapshot > MASK_POT)
					fprintf(outfile, ", ");
			}
			if(sim.out_snapshot & MASK_POT)
			{
				fprintf(outfile, "psiN");
				if (sim.out_snapshot > MASK_B)
					fprintf(outfile, ", ");
			}
			if(sim.out_snapshot & MASK_B)
			{
				fprintf(outfile, "B");
				if (sim.out_snapshot > MASK_T00)
					fprintf(outfile, ", ");
			}
			if(sim.out_snapshot & MASK_T00)
			{
				fprintf(outfile, "T00");
				if (sim.out_snapshot > MASK_TIJ)
					fprintf(outfile, ", ");
			}
			if(sim.out_snapshot & MASK_TIJ)
			{
				fprintf(outfile, "Tij");
				if (sim.out_snapshot > MASK_RBARE)
					fprintf(outfile, ", ");
			}
			if(sim.out_snapshot & MASK_RBARE)
			{
				fprintf(outfile, "rhoN");
				if (sim.out_snapshot > MASK_HIJ)
					fprintf(outfile, ", ");
			}
			if(sim.out_snapshot & MASK_HIJ)
			{
				fprintf(outfile, "hij");
				if (sim.out_snapshot > MASK_P)
					fprintf(outfile, ", ");
			}
			if(sim.out_snapshot & MASK_P)
			{
				fprintf(outfile, "p");
				if (sim.out_snapshot > MASK_GADGET)
					fprintf(outfile, ", ");
			}
			if(sim.out_snapshot & MASK_GADGET)
			{
				if (sim.out_snapshot & MASK_MULTI)
				{
					fprintf(outfile, "multi-Gadget2");
					if (sim.out_snapshot - MASK_MULTI > MASK_PCLS)
						fprintf(outfile, ", ");
				}
				else
				{
					fprintf(outfile, "Gadget2");
					if (sim.out_snapshot > MASK_PCLS)
						fprintf(outfile, ", ");
				}
			}
			if(sim.out_snapshot & MASK_PCLS)
			{
				fprintf(outfile, "particles");
				if ((sim.out_snapshot & MASK_MULTI == 0 && sim.out_snapshot > MASK_DELTA) || sim.out_snapshot - MASK_MULTI > MASK_DELTA)
					fprintf(outfile, ", ");
			}
			if(sim.out_snapshot & MASK_DELTA)
			{
				fprintf(outfile, "delta");
				if ((sim.out_snapshot & MASK_MULTI == 0 && sim.out_snapshot > MASK_DBARE) || sim.out_snapshot - MASK_MULTI > MASK_DBARE)
					fprintf(outfile, ", ");
			}
			if(sim.out_snapshot & MASK_DBARE)
			{
				fprintf(outfile, "deltaN");
			}
			fprintf(outfile, "\n");
		}
		if (sim.out_snapshot & MASK_GADGET)
		{
			fprintf(outfile, "tracer factor       = %d", sim.tracer_factor[0]);
			for (i = 1; i <= sim.baryon_flag + cosmo.num_ncdm; i++)
				fprintf(outfile, ", %d", sim.tracer_factor[i]);
			fprintf(outfile, "\n");
		}
		if (sim.downgrade_factor > 1)
			fprintf(outfile, "downgrade factor    = %d", sim.downgrade_factor);
		if (sim.num_pk > 0)
		{
			fprintf(outfile, "Pk redshifts        = ");
			for (i = 0; i < sim.num_pk - 1; i++)
				fprintf(outfile, "%lg, ", sim.z_pk[i]);
			fprintf(outfile, "%lg\n", sim.z_pk[i]);
		}
		if (sim.out_pk)
		{
			fprintf(outfile, "Pk outputs          = ");
			if(sim.out_pk & MASK_PHI)
			{
				fprintf(outfile, "phi");
				if (sim.out_pk > MASK_CHI)
					fprintf(outfile, ", ");
			}
			if(sim.out_pk & MASK_CHI)
			{
				fprintf(outfile, "chi");
				if (sim.out_pk > MASK_POT)
					fprintf(outfile, ", ");
			}
			if(sim.out_pk & MASK_POT)
			{
				fprintf(outfile, "psiN");
				if (sim.out_pk > MASK_B)
					fprintf(outfile, ", ");
			}
			if(sim.out_pk & MASK_B)
			{
				fprintf(outfile, "B");
				if (sim.out_pk > MASK_T00)
					fprintf(outfile, ", ");
			}
			if(sim.out_pk & MASK_T00)
			{
				fprintf(outfile, "T00");
				if (sim.out_pk > MASK_TIJ)
					fprintf(outfile, ", ");
			}
			if(sim.out_pk & MASK_TIJ)
			{
				fprintf(outfile, "Tij");
				if (sim.out_pk > MASK_RBARE)
					fprintf(outfile, ", ");
			}
			if(sim.out_pk & MASK_RBARE)
			{
				fprintf(outfile, "rhoN");
				if (sim.out_pk > MASK_HIJ)
					fprintf(outfile, ", ");
			}
			if(sim.out_pk & MASK_HIJ)
			{
				fprintf(outfile, "hij");
				if (sim.out_pk > MASK_P)
					fprintf(outfile, ", ");
			}
			if(sim.out_pk & MASK_P)
			{
				fprintf(outfile, "p");
				if (sim.out_pk > MASK_XSPEC)
					fprintf(outfile, ", ");
			}
			if(sim.out_pk & MASK_XSPEC)
			{
				fprintf(outfile, "X-spectra");
				if (sim.out_pk > MASK_DELTA)
					fprintf(outfile, ", ");
			}
			if(sim.out_pk & MASK_DELTA)
			{
				fprintf(outfile, "delta");
				if (sim.out_pk > MASK_DBARE)
					fprintf(outfile, ", ");
			}
			if(sim.out_pk & MASK_DBARE)
			{
				fprintf(outfile, "deltaN");
			}
			fprintf(outfile, "\n");
		}
		fprintf(outfile, "Pk bins             = %d\n", sim.numbins);
		if (sim.num_lightcone == 1)
		{
			fprintf(outfile, "lightcone vertex    = %lg, %lg, %lg\n", sim.lightcone[0].vertex[0], sim.lightcone[0].vertex[1], sim.lightcone[0].vertex[2]);
			fprintf(outfile, "lightcone outputs   = ");
			if(sim.out_lightcone[0] & MASK_PHI)
			{
				fprintf(outfile, "phi");
				if (sim.out_lightcone[0] > MASK_CHI)
					fprintf(outfile, ", ");
			}
			if(sim.out_lightcone[0] & MASK_CHI)
			{
				fprintf(outfile, "chi");
				if (sim.out_lightcone[0] > MASK_POT)
					fprintf(outfile, ", ");
			}
			if(sim.out_lightcone[0] & MASK_B)
			{
				fprintf(outfile, "B");
				if (sim.out_lightcone[0] > MASK_T00)
					fprintf(outfile, ", ");
			}
			if(sim.out_lightcone[0] & MASK_HIJ)
			{
				fprintf(outfile, "hij");
				if (sim.out_lightcone[0] > MASK_P)
					fprintf(outfile, ", ");
			}
			if(sim.out_lightcone[0] & MASK_GADGET)
				fprintf(outfile, "Gadget2");
			fprintf(outfile, "\n");
			if (sim.lightcone[0].opening > -1.)
				fprintf(outfile, "lightcone opening half-angle = %lg\n", acos(sim.lightcone[0].opening) * 180. / M_PI);
			fprintf(outfile, "lightcone distance  = %lg, %lg\n", sim.lightcone[0].distance[1], sim.lightcone[0].distance[0]);
			if (sim.lightcone[0].z != 0)
				fprintf(outfile, "lightcone redshift  = %lg\n", sim.lightcone[0].z);
			fprintf(outfile, "lightcone direction = %.15le, %.15le, %.15le\n", sim.lightcone[0].direction[0], sim.lightcone[0].direction[1], sim.lightcone[0].direction[2]);
			fprintf(outfile, "lightcone covering  = %lg\n", sim.covering[0]);
			if (sim.Nside[0][0] != sim.Nside[0][1])
				fprintf(outfile, "lightcone Nside     = %d, %d\n", sim.Nside[0][0], sim.Nside[0][1]);
			else
				fprintf(outfile, "lightcone Nside     = %d\n", sim.Nside[0][0]);
			fprintf(outfile, "lightcone pixel factor = %lg\n", sim.pixelfactor[0]);
			fprintf(outfile, "lightcone shell factor = %lg\n", sim.shellfactor[0]);	
		}
		else if (sim.num_lightcone > 1)
		{
			for (i = 0; i < sim.num_lightcone; i++)
			{
				fprintf(outfile, "lightcone %d vertex    = %lg, %lg, %lg\n", i, sim.lightcone[0].vertex[0], sim.lightcone[0].vertex[1], sim.lightcone[0].vertex[2]);
				fprintf(outfile, "lightcone %d outputs   = ", i);
				if(sim.out_lightcone[0] & MASK_PHI)
				{
					fprintf(outfile, "phi");
					if (sim.out_lightcone[0] > MASK_CHI)
						fprintf(outfile, ", ");
				}
				if(sim.out_lightcone[0] & MASK_CHI)
				{
					fprintf(outfile, "chi");
					if (sim.out_lightcone[0] > MASK_POT)
						fprintf(outfile, ", ");
				}
				if(sim.out_lightcone[0] & MASK_B)
				{
					fprintf(outfile, "B");
					if (sim.out_lightcone[0] > MASK_T00)
						fprintf(outfile, ", ");
				}
				if(sim.out_lightcone[0] & MASK_HIJ)
				{
					fprintf(outfile, "hij");
					if (sim.out_lightcone[0] > MASK_P)
						fprintf(outfile, ", ");
				}
				if(sim.out_lightcone[0] & MASK_GADGET)
					fprintf(outfile, "Gadget2");
				fprintf(outfile, "\n");
				if (sim.lightcone[0].opening > -1.)
					fprintf(outfile, "lightcone %d opening half-angle = %lg\n", i, acos(sim.lightcone[0].opening) * 180. / M_PI);
					fprintf(outfile, "lightcone %d distance  = %lg, %lg\n", i, sim.lightcone[0].distance[1], sim.lightcone[0].distance[0]);
				if (sim.lightcone[0].z != 0)
					fprintf(outfile, "lightcone %d redshift  = %lg\n", i, sim.lightcone[0].z);
				fprintf(outfile, "lightcone %d direction = %.15le, %.15le, %.15le\n", i, sim.lightcone[0].direction[0], sim.lightcone[0].direction[1], sim.lightcone[0].direction[2]);
				fprintf(outfile, "lightcone %d covering  = %lg\n", i, sim.covering[0]);
				if (sim.Nside[0][0] != sim.Nside[0][1])
					fprintf(outfile, "lightcone %d Nside     = %d, %d\n", i, sim.Nside[0][0], sim.Nside[0][1]);
				else
					fprintf(outfile, "lightcone %d Nside     = %d\n", i, sim.Nside[0][0]);
				fprintf(outfile, "lightcone %d pixel factor = %lg\n", i, sim.pixelfactor[0]);
				fprintf(outfile, "lightcone %d shell factor = %lg\n", i, sim.shellfactor[0]);	
			}
		}
		fprintf(outfile, "\n\n# hibernations\n\n");
		if (sim.num_restart > 0)
		{
			fprintf(outfile, "hibernation redshifts       = ");
			for (i = 0; i < sim.num_restart - 1; i++)
				fprintf(outfile, "%lg, ", sim.z_restart[i]);
			fprintf(outfile, "%lg\n", sim.z_restart[i]);
		}
		if (sim.wallclocklimit > 0.)
			fprintf(outfile, "hibernation wallclock limit = %lg\n", sim.wallclocklimit);
		if (sim.restart_path[0] != '\0')
			fprintf(outfile, "hibernation path            = %s\n", sim.restart_path);
		fprintf(outfile, "hibernation file base       = %s\n", sim.basename_restart);
			
		fclose(outfile);
	}
}


//////////////////////////
// hibernate
//////////////////////////
// Description:
//   creates a hibernation point by writing snapshots of the simulation data and metadata
//
// Arguments: 
//   sim            simulation metadata structure
//   ic             settings for IC generation
//   cosmo          cosmological parameter structure
//   pcls_cdm       pointer to particle handler for CDM
//   pcls_b         pointer to particle handler for baryons
//   pcls_ncdm      array of particle handlers for non-cold DM
//   phi            reference to field containing first Bardeen potential
//   chi            reference to field containing difference of Bardeen potentials
//   Bi             reference to vector field containing frame-dragging potential
//   a              scale factor
//   tau            conformal coordinate time
//   dtau           time step
//   cycle          current main control loop cycle count
//   restartcount   restart counter aka number of hibernation point (default -1)
//                  if < 0 no number is associated to the hibernation point
//
// Returns:
// 
//////////////////////////

void hibernate(metadata & sim, icsettings & ic, cosmology & cosmo, Particles<part_simple,part_simple_info,part_simple_dataType> * pcls_cdm, Particles<part_simple,part_simple_info,part_simple_dataType> * pcls_b, Particles<part_simple,part_simple_info,part_simple_dataType> * pcls_ncdm, Field<Real> & phi, Field<Real> & chi, Field<Real> & Bi, const double a, const double tau, const double dtau, const int cycle, const int restartcount = -1)
{
	string h5filename;
	char buffer[5];
	int i;
	Site x(Bi.lattice());

	h5filename.reserve(2*PARAM_MAX_LENGTH);
	h5filename.assign(sim.restart_path);
	h5filename += sim.basename_restart;
	if (restartcount >= 0)
	{
		sprintf(buffer, "%03d", restartcount);
		h5filename += buffer;
	}

	writeRestartSettings(sim, ic, cosmo, a, tau, dtau, cycle, restartcount);
	
#ifndef CHECK_B
	if (sim.vector_flag == VECTOR_PARABOLIC)
#endif
	for (x.first(); x.test(); x.next())
	{
		Bi(x,0) /= a * a * sim.numpts;
		Bi(x,1) /= a * a * sim.numpts;
		Bi(x,2) /= a * a * sim.numpts;
	}
	
#ifdef EXTERNAL_IO
	while (ioserver.openOstream()== OSTREAM_FAIL);
	
	pcls_cdm->saveHDF5_server_open(h5filename + "_cdm");
	if (sim.baryon_flag)
		pcls_b->saveHDF5_server_open(h5filename + "_b");
	for (i = 0; i < cosmo.num_ncdm; i++)
	{
		if (sim.numpcl[1+sim.baryon_flag+i] < 1) continue;
		sprintf(buffer, "%d", i);
		pcls_ncdm[i].saveHDF5_server_open(h5filename + "_ncdm" + buffer);
	}
		
	if (sim.gr_flag > 0)
	{
		phi.saveHDF5_server_open(h5filename + "_phi");
		chi.saveHDF5_server_open(h5filename + "_chi");
	}
	
	if (sim.vector_flag == VECTOR_PARABOLIC)
		Bi.saveHDF5_server_open(h5filename + "_B");
#ifdef CHECK_B
	else
		Bi.saveHDF5_server_open(h5filename + "_B_check");
#endif
		
	pcls_cdm->saveHDF5_server_write();
	if (sim.baryon_flag)
		pcls_b->saveHDF5_server_write();
	for (i = 0; i < cosmo.num_ncdm; i++)
	{
		if (sim.numpcl[1+sim.baryon_flag+i] < 1) continue;
		pcls_ncdm[i].saveHDF5_server_write();
	}
		
	if (sim.gr_flag > 0)
	{
		phi.saveHDF5_server_write(NUMBER_OF_IO_FILES);
		chi.saveHDF5_server_write(NUMBER_OF_IO_FILES);
	}

#ifndef CHECK_B
	if (sim.vector_flag == VECTOR_PARABOLIC)
#endif
		Bi.saveHDF5_server_write(NUMBER_OF_IO_FILES);
		
	ioserver.closeOstream();
#else
	pcls_cdm->saveHDF5(h5filename + "_cdm", 1);
	if (sim.baryon_flag)
		pcls_b->saveHDF5(h5filename + "_b", 1);
	for (i = 0; i < cosmo.num_ncdm; i++)
	{
		if (sim.numpcl[1+sim.baryon_flag+i] < 1) continue;
		sprintf(buffer, "%d", i);
		pcls_ncdm[i].saveHDF5(h5filename + "_ncdm" + buffer, 1);
	}
		
	if (sim.gr_flag > 0)
	{
		phi.saveHDF5(h5filename + "_phi.h5");
		chi.saveHDF5(h5filename + "_chi.h5");
	}
	
	if (sim.vector_flag == VECTOR_PARABOLIC)
		Bi.saveHDF5(h5filename + "_B.h5");
#ifdef CHECK_B
	else
		Bi.saveHDF5(h5filename + "_B_check.h5");
#endif
#endif
}

#endif


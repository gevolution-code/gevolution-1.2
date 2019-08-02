//////////////////////////
// ic_read.hpp
//////////////////////////
// 
// read initial conditions from disk
//
// Author: Julian Adamek (Université de Genève & Observatoire de Paris & Queen Mary University of London)
//
// Last modified: August 2019
//
//////////////////////////

#ifndef IC_READ_HEADER
#define IC_READ_HEADER

//////////////////////////
// readIC
//////////////////////////
// Description:
//   reads initial conditions from disk
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
//   cycle          reference to cycle counter (for restart after hibernation)
//   snapcount      reference to snapshot counter (for restart after hibernation)
//   pkcount        reference to spectra counter (for restart after hibernation)
//   restartcount   reference to restart counter (for restart after hibernation)
//
// Returns:
// 
//////////////////////////

void readIC(metadata & sim, icsettings & ic, cosmology & cosmo, const double fourpiG, double & a, double & tau, double & dtau, double & dtau_old, Particles_gevolution<part_simple,part_simple_info,part_simple_dataType> * pcls_cdm, Particles_gevolution<part_simple,part_simple_info,part_simple_dataType> * pcls_b, Particles_gevolution<part_simple,part_simple_info,part_simple_dataType> * pcls_ncdm, double * maxvel, Field<Real> * phi, Field<Real> * chi, Field<Real> * Bi, Field<Real> * source, Field<Real> * Sij, Field<Cplx> * scalarFT, Field<Cplx> * BiFT, Field<Cplx> * SijFT, PlanFFT<Cplx> * plan_phi, PlanFFT<Cplx> * plan_chi, PlanFFT<Cplx> * plan_Bi, PlanFFT<Cplx> * plan_source, PlanFFT<Cplx> * plan_Sij, int & cycle, int & snapcount, int & pkcount, int & restartcount, set<long> * IDbacklog)
{
	part_simple_info pcls_cdm_info;
	part_simple_dataType pcls_cdm_dataType;
	part_simple_info pcls_b_info;
	part_simple_dataType pcls_b_dataType;
	part_simple_info pcls_ncdm_info[MAX_PCL_SPECIES];
	part_simple_dataType pcls_ncdm_dataType;
	Real boxSize[3] = {1.,1.,1.};
	string filename;
	string buf;
	int i, p, c;
	char * ext;
	char line[PARAM_MAX_LINESIZE];
	FILE * bgfile;
	FILE * lcfile;
	struct fileDsc fd;
	gadget2_header hdr;
	long * numpcl;
	Real * dummy1;
	Real * dummy2;
	Site x(Bi->lattice());
	Site xPart(pcls_cdm->lattice());
	rKSite kFT(scalarFT->lattice());
	double d;
	long count;
	void * IDbuffer;
	void * buf2;
	set<long> IDlookup;
	
	filename.reserve(PARAM_MAX_LENGTH);
	hdr.npart[1] = 0;
	
	projection_init(phi);
	
	if (ic.z_ic > -1.)
		a = 1. / (1. + ic.z_ic);
	
	strcpy(pcls_cdm_info.type_name, "part_simple");
	pcls_cdm_info.mass = 0.;
	pcls_cdm_info.relativistic = false;
	
	pcls_cdm->initialize(pcls_cdm_info, pcls_cdm_dataType, &(phi->lattice()), boxSize);
	
	if ((ext = strstr(ic.pclfile[0], ".h5")) != NULL)
	{
		filename.assign(ic.pclfile[0], ext-ic.pclfile[0]);
		get_fileDsc_global(filename + ".h5", fd);
		numpcl = (long *) malloc(fd.numProcPerFile * sizeof(long));
		dummy1 = (Real *) malloc(3 * fd.numProcPerFile * sizeof(Real));
		dummy2 = (Real *) malloc(3 * fd.numProcPerFile * sizeof(Real));
		get_fileDsc_local(filename + ".h5", numpcl, dummy1, dummy2, fd.numProcPerFile);
		for (i = 0; i < fd.numProcPerFile; i++)
			sim.numpcl[0] += numpcl[i];
		pcls_cdm->loadHDF5(filename, 1);
		free(numpcl);
		free(dummy1);
		free(dummy2);
	}
	else
	{
		i = 0;
		do
		{
			filename.assign(ic.pclfile[0]);
			pcls_cdm->loadGadget2(filename, hdr);
			if (hdr.npart[1] == 0) break;
			if (hdr.time / a > 1.001 || hdr.time / a < 0.999)
			{
				COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": redshift indicated in Gadget2 header does not match initial redshift of simulation!" << endl;
			}
			sim.numpcl[0] += hdr.npart[1];
			i++;
			if (hdr.num_files > 1)
			{
				ext = ic.pclfile[0];
				while (strchr(ext, (int) '.') != NULL)
					ext = strchr(ext, (int) '.')+1;
				sprintf(ext, "%d", i);
			}
		}
		while (i < hdr.num_files);
		
		if (sim.baryon_flag == 1)
			pcls_cdm->parts_info()->mass = cosmo.Omega_cdm / (Real) sim.numpcl[0];
		else
			pcls_cdm->parts_info()->mass = (cosmo.Omega_cdm + cosmo.Omega_b) / (Real) sim.numpcl[0];
	}
	
	COUT << " " << sim.numpcl[0] << " cdm particles read successfully." << endl;
	maxvel[0] = pcls_cdm->updateVel(update_q, 0., &phi, 1, &a);
	
	if (sim.baryon_flag == 1)
	{
		strcpy(pcls_b_info.type_name, "part_simple");
		pcls_b_info.mass = 0.;
		pcls_b_info.relativistic = false;
	
		pcls_b->initialize(pcls_b_info, pcls_b_dataType, &(phi->lattice()), boxSize);
		
		if ((ext = strstr(ic.pclfile[1], ".h5")) != NULL)
		{
			filename.assign(ic.pclfile[1], ext-ic.pclfile[1]);
			get_fileDsc_global(filename + ".h5", fd);
			numpcl = (long *) malloc(fd.numProcPerFile * sizeof(long));
			dummy1 = (Real *) malloc(3 * fd.numProcPerFile * sizeof(Real));
			dummy2 = (Real *) malloc(3 * fd.numProcPerFile * sizeof(Real));
			get_fileDsc_local(filename + ".h5", numpcl, dummy1, dummy2, fd.numProcPerFile);
			for (i = 0; i < fd.numProcPerFile; i++)
				sim.numpcl[1] += numpcl[i];
			pcls_b->loadHDF5(filename, 1);
			free(numpcl);
			free(dummy1);
			free(dummy2);
		}
		else
		{
			i = 0;
			do
			{
				filename.assign(ic.pclfile[1]);
				pcls_b->loadGadget2(filename, hdr);
				if (hdr.npart[1] == 0) break;
				if (hdr.time / a > 1.001 || hdr.time / a < 0.999)
				{
					COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": redshift indicated in Gadget2 header does not match initial redshift of simulation!" << endl;
				}
				sim.numpcl[1] += hdr.npart[1];
				i++;
				if (hdr.num_files > 1)
				{
					ext = ic.pclfile[1];
					while (strchr(ext, (int) '.') != NULL)
						ext = strchr(ext, (int) '.')+1;
					sprintf(ext, "%d", i);
				}
			}
			while (i < hdr.num_files);
			
			pcls_b->parts_info()->mass = cosmo.Omega_b / (Real) sim.numpcl[1];
		}
		
		COUT << " " << sim.numpcl[1] << " baryon particles read successfully." << endl;
		maxvel[1] = pcls_b->updateVel(update_q, 0., &phi, 1, &a);
	}
	else
		sim.baryon_flag = 0;
		
	for (p = 0; p < cosmo.num_ncdm; p++)
	{
		if (ic.numtile[1+sim.baryon_flag+i] < 1)
		{
			maxvel[sim.baryon_flag+1+p] = 0;
			continue;
		}
		strcpy(pcls_ncdm_info[p].type_name, "part_simple");
		pcls_ncdm_info[p].mass = 0.;
		pcls_ncdm_info[p].relativistic = true;
	
		pcls_ncdm[p].initialize(pcls_ncdm_info[p], pcls_ncdm_dataType, &(phi->lattice()), boxSize);
		
		if ((ext = strstr(ic.pclfile[sim.baryon_flag+1+p], ".h5")) != NULL)
		{
			filename.assign(ic.pclfile[sim.baryon_flag+1+p], ext-ic.pclfile[sim.baryon_flag+1+p]);
			get_fileDsc_global(filename + ".h5", fd);
			numpcl = (long *) malloc(fd.numProcPerFile * sizeof(long));
			dummy1 = (Real *) malloc(3 * fd.numProcPerFile * sizeof(Real));
			dummy2 = (Real *) malloc(3 * fd.numProcPerFile * sizeof(Real));
			get_fileDsc_local(filename + ".h5", numpcl, dummy1, dummy2, fd.numProcPerFile);
			for (i = 0; i < fd.numProcPerFile; i++)
				sim.numpcl[1] += numpcl[i];
			pcls_ncdm[p].loadHDF5(filename, 1);
			free(numpcl);
			free(dummy1);
			free(dummy2);
		}
		else
		{
			i = 0;
			do
			{
				filename.assign(ic.pclfile[sim.baryon_flag+1+p]);
				pcls_ncdm[p].loadGadget2(filename, hdr);
				if (hdr.npart[1] == 0) break;
				if (hdr.time / a > 1.001 || hdr.time / a < 0.999)
				{
					COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": redshift indicated in Gadget2 header does not match initial redshift of simulation!" << endl;
				}
				sim.numpcl[sim.baryon_flag+1+p] += hdr.npart[1];
				i++;
				if (hdr.num_files > 1)
				{
					ext = ic.pclfile[sim.baryon_flag+1+p];
					while (strchr(ext, (int) '.') != NULL)
						ext = strchr(ext, (int) '.')+1;
					sprintf(ext, "%d", i);
				}
			}
			while (i < hdr.num_files);
			
			pcls_ncdm[p].parts_info()->mass = cosmo.Omega_ncdm[p] / (Real) sim.numpcl[sim.baryon_flag+1+p];
		}
		
		COUT << " " << sim.numpcl[sim.baryon_flag+1+p] << " ncdm particles read successfully." << endl;
		maxvel[sim.baryon_flag+1+p] = pcls_ncdm[p].updateVel(update_q, 0., &phi, 1, &a);
	}
	
	if (sim.gr_flag > 0 && ic.metricfile[0][0] != '\0')
	{
		filename.assign(ic.metricfile[0]);
		phi->loadHDF5(filename);
	}
	else
	{
		projection_init(source);
		scalarProjectionCIC_project(pcls_cdm, source);
		if (sim.baryon_flag)
			scalarProjectionCIC_project(pcls_b, source);
		scalarProjectionCIC_comm(source);
	
		plan_source->execute(FFT_FORWARD);
	
		kFT.first();
		if (kFT.coord(0) == 0 && kFT.coord(1) == 0 && kFT.coord(2) == 0)
			(*scalarFT)(kFT) = Cplx(0.,0.);
				
		solveModifiedPoissonFT(*scalarFT, *scalarFT, fourpiG / a, 3. * sim.gr_flag * (Hconf(a, fourpiG, cosmo) * Hconf(a, fourpiG, cosmo) + fourpiG * cosmo.Omega_m / a));
		plan_phi->execute(FFT_BACKWARD);
	}

	phi->updateHalo();
	
	if (ic.restart_tau > 0.)
		tau = ic.restart_tau;
	else
		tau = particleHorizon(a, fourpiG, cosmo);
		
	if (ic.restart_dtau > 0.)
		dtau_old = ic.restart_dtau;
		
	if (sim.Cf / (double) sim.numpts < sim.steplimit / Hconf(a, fourpiG, cosmo))
		dtau = sim.Cf / (double) sim.numpts;
	else
		dtau = sim.steplimit / Hconf(a, fourpiG, cosmo);

	if (ic.restart_cycle >= 0)
	{	
#ifndef CHECK_B	
		if (sim.vector_flag == VECTOR_PARABOLIC)
#endif
		{
			filename.assign(ic.metricfile[2*sim.gr_flag]);
			Bi->loadHDF5(filename);

			for (x.first(); x.test(); x.next())
			{
				(*Bi)(x,0) *= a * a / (sim.numpts * sim.numpts);
				(*Bi)(x,1) *= a * a / (sim.numpts * sim.numpts);
				(*Bi)(x,2) *= a * a / (sim.numpts * sim.numpts);
			}
			plan_Bi->execute(FFT_FORWARD);
		}
		
		if (sim.gr_flag > 0)
		{
			filename.assign(ic.metricfile[1]);
			chi->loadHDF5(filename);
			chi->updateHalo();
		}
		
		if (parallel.isRoot())
		{
			sprintf(line, "%s%s_background.dat", sim.output_path, sim.basename_generic);
			bgfile = fopen(line, "r");
			if (bgfile == NULL)
			{
				COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": unable to locate file for background output! A new file will be created" << endl;
				bgfile = fopen(line, "w");
				if (bgfile == NULL)
				{
					COUT << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to create file for background output!" << endl;
					parallel.abortForce();
				}
				else
				{
					fprintf(bgfile, "# background statistics\n# cycle   tau/boxsize    a             conformal H/H0  phi(k=0)       T00(k=0)\n");
					fclose(bgfile);
				}
			}
			else
			{
				buf.reserve(PARAM_MAX_LINESIZE);
				buf.clear();
				
				if (fgets(line, PARAM_MAX_LINESIZE, bgfile) == 0)
				{
					COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": unable to read file for background output! A new file will be created" << endl;
				}
				else if (line[0] != '#')
				{
					COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": file for background output has unexpected format! Contents will be overwritten!" << endl;
				}
				else
				{
					if (fgets(line, PARAM_MAX_LINESIZE, bgfile) == 0)
					{
						COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": unable to read file for background output! A new file will be created" << endl;
					}
					else if (line[0] != '#')
					{
						COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": file for background output has unexpected format! Contents will be overwritten!" << endl;
					}
				}
				
				while (fgets(line, PARAM_MAX_LINESIZE, bgfile) != 0)
				{
					if (sscanf(line, " %d", &i) != 1)
						break;
					
					if (i > ic.restart_cycle)
						break;
					else
						buf += line;
				}
	
				fclose(bgfile);
				
				sprintf(line, "%s%s_background.dat", sim.output_path, sim.basename_generic);
				bgfile = fopen(line, "w");
				
				if (bgfile == NULL)
				{
					COUT << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to create file for background output!" << endl;
					parallel.abortForce();
				}
				else
				{
					fprintf(bgfile, "# background statistics\n# cycle   tau/boxsize    a             conformal H/H0  phi(k=0)       T00(k=0)\n");
					fwrite((const void *) buf.data(), sizeof(char), buf.length(), bgfile);
					fclose(bgfile);
					buf.clear();
				}
			}
		}
		
		for (i = 0; i < sim.num_lightcone; i++)
		{
			if (parallel.isRoot())
			{
				if (sim.num_lightcone > 1)
					sprintf(line, "%s%s%d_info.dat", sim.output_path, sim.basename_lightcone, i);
				else
					sprintf(line, "%s%s_info.dat", sim.output_path, sim.basename_lightcone);
						
				lcfile = fopen(line, "r");
				
				if (lcfile == NULL)
				{
					COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": unable to locate file for lightcone information! A new file will be created" << endl;
					lcfile = fopen(line, "w");
					if (lcfile == NULL)
					{
						COUT << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to create file for lightcone information!" << endl;
						parallel.abortForce();
					}
					else
					{
						if (sim.num_lightcone > 1)
							fprintf(lcfile, "# information file for lightcone %d\n# geometric parameters:\n# vertex = (%f, %f, %f) Mpc/h\n# redshift = %f\n# distance = (%f - %f) Mpc/h\n# opening half-angle = %f degrees\n# direction = (%f, %f, %f)\n# cycle   tau/boxsize    a              pcl_inner        pcl_outer        metric_inner     metric_outer\n", i, sim.lightcone[i].vertex[0]*sim.boxsize, sim.lightcone[i].vertex[1]*sim.boxsize, sim.lightcone[i].vertex[2]*sim.boxsize, sim.lightcone[i].z, sim.lightcone[i].distance[0]*sim.boxsize, sim.lightcone[i].distance[1]*sim.boxsize, (sim.lightcone[i].opening > -1.) ? acos(sim.lightcone[i].opening) * 180. / M_PI : 180., sim.lightcone[i].direction[0], sim.lightcone[i].direction[1], sim.lightcone[i].direction[2]);
						else
							fprintf(lcfile, "# information file for lightcone\n# geometric parameters:\n# vertex = (%f, %f, %f) Mpc/h\n# redshift = %f\n# distance = (%f - %f) Mpc/h\n# opening half-angle = %f degrees\n# direction = (%f, %f, %f)\n# cycle   tau/boxsize    a              pcl_inner        pcl_outer        metric_inner     metric_outer\n", sim.lightcone[i].vertex[0]*sim.boxsize, sim.lightcone[i].vertex[1]*sim.boxsize, sim.lightcone[i].vertex[2]*sim.boxsize, sim.lightcone[i].z, sim.lightcone[i].distance[0]*sim.boxsize, sim.lightcone[i].distance[1]*sim.boxsize, (sim.lightcone[i].opening > -1.) ? acos(sim.lightcone[i].opening) * 180. / M_PI : 180., sim.lightcone[i].direction[0], sim.lightcone[i].direction[1], sim.lightcone[i].direction[2]);
						fclose(lcfile);
					}
				}
				else
				{
					buf.reserve(PARAM_MAX_LINESIZE);
					buf.clear();
				
					for(p = 0; p < 8; p++)
					{
						if (fgets(line, PARAM_MAX_LINESIZE, lcfile) == 0)
						{
							COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": unable to read file for lightcone information! A new file will be created" << endl;
							break;
						}
						else if (line[0] != '#')
						{
							COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": file for lightcone information has unexpected format! Contents will be overwritten!" << endl;
							break;
						}
					}
					
					p = 0;
					while (fgets(line, PARAM_MAX_LINESIZE, lcfile) != 0)
					{
						if (sscanf(line, " %d", &c) != 1)
							break;
						
						if (c > ic.restart_cycle)
							break;
						else
						{
							buf += line;
							p++;
						}
					}
		
					fclose(lcfile);
					
					if (sim.num_lightcone > 1)
						sprintf(line, "%s%s%d_info.dat", sim.output_path, sim.basename_lightcone, i);
					else
						sprintf(line, "%s%s_info.dat", sim.output_path, sim.basename_lightcone);
						
					lcfile = fopen(line, "w");
					
					if (lcfile == NULL)
					{
						COUT << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to create file for lightcone information!" << endl;
						parallel.abortForce();
					}
					else
					{
						if (sim.num_lightcone > 1)
							fprintf(lcfile, "# information file for lightcone %d\n# geometric parameters:\n# vertex = (%f, %f, %f) Mpc/h\n# redshift = %f\n# distance = (%f - %f) Mpc/h\n# opening half-angle = %f degrees\n# direction = (%f, %f, %f)\n# cycle   tau/boxsize    a              pcl_inner        pcl_outer        metric_inner     metric_outer\n", i, sim.lightcone[i].vertex[0]*sim.boxsize, sim.lightcone[i].vertex[1]*sim.boxsize, sim.lightcone[i].vertex[2]*sim.boxsize, sim.lightcone[i].z, sim.lightcone[i].distance[0]*sim.boxsize, sim.lightcone[i].distance[1]*sim.boxsize, (sim.lightcone[i].opening > -1.) ? acos(sim.lightcone[i].opening) * 180. / M_PI : 180., sim.lightcone[i].direction[0], sim.lightcone[i].direction[1], sim.lightcone[i].direction[2]);
						else
							fprintf(lcfile, "# information file for lightcone\n# geometric parameters:\n# vertex = (%f, %f, %f) Mpc/h\n# redshift = %f\n# distance = (%f - %f) Mpc/h\n# opening half-angle = %f degrees\n# direction = (%f, %f, %f)\n# cycle   tau/boxsize    a              pcl_inner        pcl_outer        metric_inner     metric_outer\n", sim.lightcone[i].vertex[0]*sim.boxsize, sim.lightcone[i].vertex[1]*sim.boxsize, sim.lightcone[i].vertex[2]*sim.boxsize, sim.lightcone[i].z, sim.lightcone[i].distance[0]*sim.boxsize, sim.lightcone[i].distance[1]*sim.boxsize, (sim.lightcone[i].opening > -1.) ? acos(sim.lightcone[i].opening) * 180. / M_PI : 180., sim.lightcone[i].direction[0], sim.lightcone[i].direction[1], sim.lightcone[i].direction[2]);
						fwrite((const void *) buf.data(), sizeof(char), buf.length(), lcfile);
						fclose(lcfile);
						buf.clear();
					}
					
					if (sim.num_lightcone > 1)
						sprintf(line, "%s%s%d_info.bin", sim.output_path, sim.basename_lightcone, i);
					else
						sprintf(line, "%s%s_info.bin", sim.output_path, sim.basename_lightcone);
						
					lcfile = fopen(line, "r");
					
					if (lcfile == NULL)
					{
						if (p > 0)
						{
							COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": unable to read binary file for lightcone information!" << endl;
						}
					}
					else
					{
						if (p > 0)
						{
							buf2 = malloc(p * (sizeof(int) + 6 * sizeof(double)));
							c = fread(buf2, 1, p * (sizeof(int) + 6 * sizeof(double)), lcfile);
							if (c != p * (sizeof(int) + 6 * sizeof(double)))
							{
								COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": unable to read binary file for lightcone information!" << endl;
							}
							
							fclose(lcfile);
							lcfile = fopen(line, "r");
							
							if (lcfile == NULL)
							{
								COUT << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to create file for lightcone information!" << endl;
								parallel.abortForce();
							}
							else
							{
								fwrite(buf2, 1, c, lcfile);
								fclose(lcfile);
							}
							
							free(buf2);
						}
						else
						{
							fclose(lcfile);
							lcfile = fopen(line, "r");
							if (lcfile != NULL)
								fclose(lcfile);
						}
					}
				}
			}
		
			d = particleHorizon(1. / (1. + sim.lightcone[i].z), fourpiG, cosmo);
			if (sim.out_lightcone[i] & MASK_GADGET && sim.lightcone[i].distance[0] > d - tau + 0.5 * dtau_old && sim.lightcone[i].distance[1] <= d - tau + 0.5 * dtau_old && d - tau + 0.5 * dtau_old > 0.)
			{
				for (p = 0; p < 1 + sim.baryon_flag + cosmo.num_ncdm; p++)
				{
					if (sim.numpcl[p] == 0) continue;
					
					if (p == 0)
					{
						for (xPart.first(); xPart.test(); xPart.next())
						{
							for (std::list<part_simple>::iterator it = (pcls_cdm->field())(xPart).parts.begin(); it != (pcls_cdm->field())(xPart).parts.end(); ++it)
								IDlookup.insert((*it).ID);
						}
					}
					else if (p == 1 && sim.baryon_flag > 0)
					{
						for (xPart.first(); xPart.test(); xPart.next())
						{
							for (std::list<part_simple>::iterator it = (pcls_b->field())(xPart).parts.begin(); it != (pcls_b->field())(xPart).parts.end(); ++it)
								IDlookup.insert((*it).ID);
						}
					}
					else
					{
						for (xPart.first(); xPart.test(); xPart.next())
						{
							for (std::list<part_simple>::iterator it = (pcls_ncdm[p-1-sim.baryon_flag].field())(xPart).parts.begin(); it != (pcls_ncdm[p-1-sim.baryon_flag].field())(xPart).parts.end(); ++it)
								IDlookup.insert((*it).ID);
						}
					}
				
					if (parallel.isRoot())
					{
						if (sim.num_lightcone > 1)
						{
							if (p == 0)
								sprintf(line, "%s%s%d_%04d_cdm", sim.output_path, sim.basename_lightcone, i, ic.restart_cycle);
							else if (p == 1 && sim.baryon_flag > 0)
								sprintf(line, "%s%s%d_%04d_b", sim.output_path, sim.basename_lightcone, i, ic.restart_cycle);
							else
								sprintf(line, "%s%s%d_%04d_ncdm%d", sim.output_path, sim.basename_lightcone, i, ic.restart_cycle, p-1-sim.baryon_flag);
						}
						else
						{
							if (p == 0)
								sprintf(line, "%s%s_%04d_cdm", sim.output_path, sim.basename_lightcone, ic.restart_cycle);
							else if (p == 1 && sim.baryon_flag > 0)
								sprintf(line, "%s%s_%04d_b", sim.output_path, sim.basename_lightcone, ic.restart_cycle);
							else
								sprintf(line, "%s%s_%04d_ncdm%d", sim.output_path, sim.basename_lightcone, ic.restart_cycle, p-1-sim.baryon_flag);
						}
						
						lcfile = fopen(line, "r");
					
						if (lcfile == NULL)
						{
							COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": unable to open " << line << " for retrieving particle ID backlog" << endl;
							count = 0;
						}
						else if (fseek(lcfile, 4, SEEK_SET) || fread(&hdr, sizeof(gadget2_header), 1, lcfile) != 1 || fseek(lcfile, 276l + 6l * sizeof(float) * (hdr.npartTotal[1] + ((long) hdr.npartTotalHW[1] << 32)), SEEK_SET))
						{
							COUT << COLORTEXT_YELLOW << " /!\\ warning" << COLORTEXT_RESET << ": unable to read " << line << " for retrieving particle ID backlog" << endl;
							count = 0;
							fclose(lcfile);
							lcfile = NULL;
						}
						else
							count = hdr.npartTotal[1] + ((long) hdr.npartTotalHW[1] << 32);
					}
				
					parallel.broadcast<long>(count, 0);
				
					IDbuffer = malloc(PCLBUFFER * GADGET_ID_BYTES);
				
					while (count > 0)
					{
						if (count > PCLBUFFER)
						{
							if (parallel.isRoot() && fread(IDbuffer, 1, PCLBUFFER * GADGET_ID_BYTES, lcfile) != PCLBUFFER * GADGET_ID_BYTES)
							{
								COUT << COLORTEXT_RED << " /!\\ error" << COLORTEXT_RESET << ": unable to read particle ID block from " << line << endl;
								parallel.abortForce();
							}
						
							parallel.broadcast<char>((char *) IDbuffer, PCLBUFFER * GADGET_ID_BYTES, 0);
						
							for (int j = 0; j < PCLBUFFER; j++)
							{
#if GADGET_ID_BYTES == 8
								if (IDlookup.erase((long) *(((int64_t *) IDbuffer) + j)))
									IDbacklog[p].insert((long) *(((int64_t *) IDbuffer) + j));
#else
								if (IDlookup.erase((long) *(((int32_t *) IDbuffer) + j)))
									IDbacklog[p].insert((long) *(((int32_t *) IDbuffer) + j));
#endif
							}
						
							count -= PCLBUFFER;
						}
						else
						{
							if (parallel.isRoot() && fread(IDbuffer, 1, count * GADGET_ID_BYTES, lcfile) != count * GADGET_ID_BYTES)
							{
								COUT << COLORTEXT_RED << " /!\\ error" << COLORTEXT_RESET << ": unable to read particle ID block from " << line << endl;
								parallel.abortForce();
							}
						
							parallel.broadcast<char>((char *) IDbuffer, count * GADGET_ID_BYTES, 0);
						
							for (int j = 0; j < count; j++)
							{
#if GADGET_ID_BYTES == 8
								if (IDlookup.erase((long) *(((int64_t *) IDbuffer) + j)))
									IDbacklog[p].insert((long) *(((int64_t *) IDbuffer) + j));
#else
								if (IDlookup.erase((long) *(((int32_t *) IDbuffer) + j)))
									IDbacklog[p].insert((long) *(((int32_t *) IDbuffer) + j));
#endif
							}
						
							count = 0;
						}
					}
				
					free(IDbuffer);
				
					if (parallel.isRoot() && lcfile != NULL)
						fclose(lcfile);
					
					IDlookup.clear();
				}
			}
		}
	}
	else
	{
		projection_init(Bi);
		projection_T0i_project(pcls_cdm, Bi, phi);
		if (sim.baryon_flag)
			projection_T0i_project(pcls_b, Bi, phi);
		projection_T0i_comm(Bi);
		plan_Bi->execute(FFT_FORWARD);
		projectFTvector(*BiFT, *BiFT, fourpiG / (double) sim.numpts / (double) sim.numpts);	
		plan_Bi->execute(FFT_BACKWARD);	
		Bi->updateHalo();
		
		projection_init(Sij);
		projection_Tij_project(pcls_cdm, Sij, a, phi);
		if (sim.baryon_flag)
			projection_Tij_project(pcls_b, Sij, a, phi);
		projection_Tij_comm(Sij);
	
		prepareFTsource<Real>(*phi, *Sij, *Sij, 2. * fourpiG / a / (double) sim.numpts / (double) sim.numpts);	
		plan_Sij->execute(FFT_FORWARD);	
		projectFTscalar(*SijFT, *scalarFT);
		plan_chi->execute(FFT_BACKWARD);		
		chi->updateHalo();
	}

	if (ic.restart_cycle >= 0)
		cycle = ic.restart_cycle + 1;
		
	while (snapcount < sim.num_snapshot && 1. / a < sim.z_snapshot[snapcount] + 1.)
		snapcount++;
		
	while (pkcount < sim.num_pk && 1. / a < sim.z_pk[pkcount] + 1.)
		pkcount++;
		
	while (restartcount < sim.num_restart && 1. / a < sim.z_restart[restartcount] + 1.)
		restartcount++;
}

#endif


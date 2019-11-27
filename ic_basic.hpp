//////////////////////////
// ic_basic.hpp
//////////////////////////
// 
// basic initial condition generator for gevolution
//
// Author: Julian Adamek (Université de Genève & Observatoire de Paris & Queen Mary University of London)
//
// Last modified: November 2019
//
//////////////////////////

#ifndef IC_BASIC_HEADER
#define IC_BASIC_HEADER

#include "prng_engine.hpp"
#include <gsl/gsl_spline.h>
#include "parser.hpp"
#include "tools.hpp"

#ifndef Cplx
#define Cplx Imag
#endif  

#define MAX_LINESIZE 2048

using namespace std;
using namespace LATfield2;

// should be larger than maximum Ngrid
#ifndef HUGE_SKIP
#define HUGE_SKIP   65536
#endif


//////////////////////////
// displace_pcls_ic_basic
//////////////////////////
// Description:
//   displaces particles according to gradient of displacement field
//   (accepts two displacement fields for baryon treatment = hybrid)
// 
// Arguments:
//   coeff             coefficient to be applied to displacement
//   lat_resolution    1 / Ngrid
//   part              particle to be displaced
//   ref_dist          distance vector (in lattice units) to reference lattice point
//   partInfo          particle metadata
//   fields            array of pointers to displacement fields
//   sites             array of respective sites
//   nfields           number of fields passed
//   params            additional parameters (ignored)
//   outputs           array for reduction variables; first entry will contain the displacement
//   noutputs          number of reduction variables
//
// Returns:
// 
//////////////////////////

void displace_pcls_ic_basic(double coeff, double lat_resolution, part_simple * part, double * ref_dist, part_simple_info partInfo, Field<Real> ** fields, Site * sites, int nfield, double * params, double * outputs, int noutputs)
{
	int i;
	Real gradxi[3] = {0, 0, 0};
	
	if (nfield > 1 && (*part).ID % 8 == 0)
		i = 1;
	else
		i = 0;
	
	gradxi[0] = (1.-ref_dist[1]) * (1.-ref_dist[2]) * ((*fields[i])(sites[i]+0) - (*fields[i])(sites[i]));
	gradxi[1] = (1.-ref_dist[0]) * (1.-ref_dist[2]) * ((*fields[i])(sites[i]+1) - (*fields[i])(sites[i]));
	gradxi[2] = (1.-ref_dist[0]) * (1.-ref_dist[1]) * ((*fields[i])(sites[i]+2) - (*fields[i])(sites[i]));
	gradxi[0] += ref_dist[1] * (1.-ref_dist[2]) * ((*fields[i])(sites[i]+1+0) - (*fields[i])(sites[i]+1));
	gradxi[1] += ref_dist[0] * (1.-ref_dist[2]) * ((*fields[i])(sites[i]+1+0) - (*fields[i])(sites[i]+0));
	gradxi[2] += ref_dist[0] * (1.-ref_dist[1]) * ((*fields[i])(sites[i]+2+0) - (*fields[i])(sites[i]+0));
	gradxi[0] += (1.-ref_dist[1]) * ref_dist[2] * ((*fields[i])(sites[i]+2+0) - (*fields[i])(sites[i]+2));
	gradxi[1] += (1.-ref_dist[0]) * ref_dist[2] * ((*fields[i])(sites[i]+2+1) - (*fields[i])(sites[i]+2));
	gradxi[2] += (1.-ref_dist[0]) * ref_dist[1] * ((*fields[i])(sites[i]+2+1) - (*fields[i])(sites[i]+1));
	gradxi[0] += ref_dist[1] * ref_dist[2] * ((*fields[i])(sites[i]+2+1+0) - (*fields[i])(sites[i]+2+1));
	gradxi[1] += ref_dist[0] * ref_dist[2] * ((*fields[i])(sites[i]+2+1+0) - (*fields[i])(sites[i]+2+0));
	gradxi[2] += ref_dist[0] * ref_dist[1] * ((*fields[i])(sites[i]+2+1+0) - (*fields[i])(sites[i]+1+0));
	
	gradxi[0] /= lat_resolution;
	gradxi[1] /= lat_resolution;
	gradxi[2] /= lat_resolution;
	
	if (noutputs > 0)
		*outputs = coeff * sqrt(gradxi[0]*gradxi[0] + gradxi[1]*gradxi[1] + gradxi[2]*gradxi[2]);

	for (i = 0; i < 3; i++) (*part).pos[i] += coeff*gradxi[i];
}


//////////////////////////
// initialize_q_ic_basic
//////////////////////////
// Description:
//   initializes velocities proportional to gradient of potential
//   (accepts two potentials for baryon treatment = hybrid)
// 
// Arguments:
//   coeff             coefficient to be applied to velocity
//   lat_resolution    1 / Ngrid
//   part              particle to be manipulated
//   ref_dist          distance vector (in lattice units) to reference lattice point
//   partInfo          particle metadata
//   fields            array of pointers to velocity potentials
//   sites             array of respective sites
//   nfields           number of fields passed
//   params            additional parameters (ignored)
//   outputs           array for reduction variables (ignored)
//   noutputs          number of reduction variables
//
// Returns: square of the velocity, (q/m)^2
// 
//////////////////////////

Real initialize_q_ic_basic(double coeff, double lat_resolution, part_simple * part, double * ref_dist, part_simple_info partInfo, Field<Real> ** fields, Site * sites, int nfield, double * params, double * outputs, int noutputs)
{
	int i;
	Real gradPhi[3] = {0, 0, 0};
	Real v2 = 0.;
	
	if (nfield > 1 && (*part).ID % 8 == 0)
		i = 1;
	else
		i = 0;
	
	gradPhi[0] = (1.-ref_dist[1]) * (1.-ref_dist[2]) * ((*fields[i])(sites[i]+0) - (*fields[i])(sites[i]));
	gradPhi[1] = (1.-ref_dist[0]) * (1.-ref_dist[2]) * ((*fields[i])(sites[i]+1) - (*fields[i])(sites[i]));
	gradPhi[2] = (1.-ref_dist[0]) * (1.-ref_dist[1]) * ((*fields[i])(sites[i]+2) - (*fields[i])(sites[i]));
	gradPhi[0] += ref_dist[1] * (1.-ref_dist[2]) * ((*fields[i])(sites[i]+1+0) - (*fields[i])(sites[i]+1));
	gradPhi[1] += ref_dist[0] * (1.-ref_dist[2]) * ((*fields[i])(sites[i]+1+0) - (*fields[i])(sites[i]+0));
	gradPhi[2] += ref_dist[0] * (1.-ref_dist[1]) * ((*fields[i])(sites[i]+2+0) - (*fields[i])(sites[i]+0));
	gradPhi[0] += (1.-ref_dist[1]) * ref_dist[2] * ((*fields[i])(sites[i]+2+0) - (*fields[i])(sites[i]+2));
	gradPhi[1] += (1.-ref_dist[0]) * ref_dist[2] * ((*fields[i])(sites[i]+2+1) - (*fields[i])(sites[i]+2));
	gradPhi[2] += (1.-ref_dist[0]) * ref_dist[1] * ((*fields[i])(sites[i]+2+1) - (*fields[i])(sites[i]+1));
	gradPhi[0] += ref_dist[1] * ref_dist[2] * ((*fields[i])(sites[i]+2+1+0) - (*fields[i])(sites[i]+2+1));
	gradPhi[1] += ref_dist[0] * ref_dist[2] * ((*fields[i])(sites[i]+2+1+0) - (*fields[i])(sites[i]+2+0));
	gradPhi[2] += ref_dist[0] * ref_dist[1] * ((*fields[i])(sites[i]+2+1+0) - (*fields[i])(sites[i]+1+0));
	
	gradPhi[0] /= lat_resolution;
	gradPhi[1] /= lat_resolution;
	gradPhi[2] /= lat_resolution;  
	
	for (i = 0 ; i < 3; i++)
	{
		(*part).vel[i] = -gradPhi[i] * coeff;
		v2 += (*part).vel[i] * (*part).vel[i];
	}
	
	return v2;
}


//////////////////////////
// Pk_primordial
//////////////////////////
// Description:
//   power spectrum of the primordial curvature perturbation
// 
// Arguments:
//   k          wavenumber in units of inverse Mpc (!)
//   ic         settings for IC generation (datastructure containing parameters)
//
// Returns: amplitude of primordial power spectrum at given wavenumber
// 
//////////////////////////

inline double Pk_primordial(const double k, const icsettings & ic)
{
	return ic.A_s * pow(k / ic.k_pivot, ic.n_s - 1.);  // note that k_pivot is in units of inverse Mpc!
}


//////////////////////////
// loadHomogeneousTemplate
//////////////////////////
// Description:
//   loads a homogeneous template from a GADGET-2 file
// 
// Arguments:
//   filename   string containing the path to the template file
//   numpart    will contain the number of particles of the template
//   partdata   will contain the particle positions (memory will be allocated)
//
// Returns:
// 
//////////////////////////

void loadHomogeneousTemplate(const char * filename, long & numpart, float * & partdata)
{
	int i;

	if (parallel.grid_rank()[0] == 0) // read file
	{
		FILE * templatefile;
		int blocksize1, blocksize2, num_read;
		gadget2_header filehdr;
		
		templatefile = fopen(filename, "r");
		
		if (templatefile == NULL)
		{
			cerr << " proc#" << parallel.rank() << ": error in loadHomogeneousTemplate! Unable to open template file " << filename << "." << endl;
			parallel.abortForce();
		}
		
		fread(&blocksize1, sizeof(int), 1, templatefile);
		if (blocksize1 != sizeof(filehdr))
		{
			cerr << " proc#" << parallel.rank() << ": error in loadHomogeneousTemplate! Unknown template file format - header not recognized." << endl;
			fclose(templatefile);
			parallel.abortForce();
		}		
		fread(&filehdr, sizeof(filehdr), 1, templatefile);
		fread(&blocksize2, sizeof(int), 1, templatefile);
		if (blocksize1 != blocksize2)
		{
			cerr << " proc#" << parallel.rank() << ": error in loadHomogeneousTemplate! Unknown template file format - block size mismatch while reading header." << endl;
			fclose(templatefile);
			parallel.abortForce();
		}
		
		// analyze header for compatibility
		if (filehdr.num_files != 1)
		{
			cerr << " proc#" << parallel.rank() << ": error in loadHomogeneousTemplate! Multiple input files (" << filehdr.num_files << ") currently not supported." << endl;
			fclose(templatefile);
			parallel.abortForce();
		}
		if (filehdr.BoxSize <= 0.)
		{
			cerr << " proc#" << parallel.rank() << ": error in loadHomogeneousTemplate! BoxSize = " << filehdr.BoxSize << " not allowed." << endl;
			fclose(templatefile);
			parallel.abortForce();
		}
		if (filehdr.npart[1] <= 0)
		{
			cerr << " proc#" << parallel.rank() << ": error in loadHomogeneousTemplate! No particles declared." << endl;
			fclose(templatefile);
			parallel.abortForce();
		}
		
		partdata = (float *) malloc(3 * sizeof(float) * filehdr.npart[1]);
		if (partdata == NULL)
		{
			cerr << " proc#" << parallel.rank() << ": error in loadHomogeneousTemplate! Memory error." << endl;
			fclose(templatefile);
			parallel.abortForce();
		}
		
		fread(&blocksize1, sizeof(int), 1, templatefile);
		if (filehdr.npart[0] > 0)
		{
			if (fseek(templatefile, 3 * sizeof(float) * filehdr.npart[0], SEEK_CUR))
			{
				cerr << " proc#" << parallel.rank() << ": error in loadHomogeneousTemplate! Unsuccesful attempt to skip gas particles (" << filehdr.npart[0] << ")." << endl;
				fclose(templatefile);
				parallel.abortForce();
			}
		}
		num_read = fread(partdata, sizeof(float), 3 * filehdr.npart[1], templatefile);
		if (num_read != 3 * filehdr.npart[1])
		{
			cerr << " proc#" << parallel.rank() << ": error in loadHomogeneousTemplate! Unable to read particle data." << endl;
			fclose(templatefile);
			parallel.abortForce();
		}
		for (i = 2; i < 6; i++)
		{
			if (filehdr.npart[i] > 0)
			{
				if (fseek(templatefile, 3 * sizeof(float) * filehdr.npart[i], SEEK_CUR))
				{
					cerr << " proc#" << parallel.rank() << ": error in loadHomogeneousTemplate! Unsuccesful attempt to skip particles (" << filehdr.npart[i] << ")." << endl;
					parallel.abortForce();
				}
			}
		}
		fread(&blocksize2, sizeof(int), 1, templatefile);
		if (blocksize1 != blocksize2)
			{
			cerr << " proc#" << parallel.rank() << ": error in loadHomogeneousTemplate! Unknown template file format - block size mismatch while reading particles." << endl;
			fclose(templatefile);
			parallel.abortForce();
		}
		
		fclose(templatefile);
		
		// reformat and check particle data
		for (i = 0; i < 3 * filehdr.npart[1]; i++)
		{
			partdata[i] /= filehdr.BoxSize;
			if (partdata[i] < 0. || partdata[i] > 1.)
			{
				cerr << " proc#" << parallel.rank() << ": error in loadHomogeneousTemplate! Particle data corrupted." << endl;
				parallel.abortForce();
			}
		}
		numpart = (long) filehdr.npart[1];
		
		parallel.broadcast_dim0<long>(numpart, 0);
	}
	else
	{
		parallel.broadcast_dim0<long>(numpart, 0);
		
		if (numpart <= 0)
		{
			cerr << " proc#" << parallel.rank() << ": error in loadHomogeneousTemplate! Communication error." << endl;
			parallel.abortForce();
		}
		
		partdata = (float *) malloc(3 * sizeof(float) * numpart);
		if (partdata == NULL)
		{
			cerr << " proc#" << parallel.rank() << ": error in loadHomogeneousTemplate! Memory error." << endl;
			parallel.abortForce();
		}
	}
	
	parallel.broadcast_dim0<float>(partdata, 3 * numpart, 0);
}


//////////////////////////
// loadPowerSpectrum
//////////////////////////
// Description:
//   loads a tabulated matter power spectrum from a file
// 
// Arguments:
//   filename   string containing the path to the template file
//   pkspline   will point to the gsl_spline which holds the tabulated
//              power spectrum (memory will be allocated)
//   boxsize    comoving box size (in the same units as used in the file)
//
// Returns:
// 
//////////////////////////

void loadPowerSpectrum(const char * filename, gsl_spline * & pkspline, const double boxsize)
{
	int i = 0, numpoints = 0;
	double * k;
	double * pk;

	if (parallel.grid_rank()[0] == 0) // read file
	{
		FILE * pkfile;
		char line[MAX_LINESIZE];
		double dummy1, dummy2;
		
		line[MAX_LINESIZE-1] = 0;
		
		pkfile = fopen(filename, "r");
		
		if (pkfile == NULL)
		{
			cerr << " proc#" << parallel.rank() << ": error in loadPowerSpectrum! Unable to open file " << filename << "." << endl;
			parallel.abortForce();
		}
		
		while (!feof(pkfile) && !ferror(pkfile))
		{
			fgets(line, MAX_LINESIZE, pkfile);
			if (line[MAX_LINESIZE-1] != 0)
			{
				cerr << " proc#" << parallel.rank() << ": error in loadPowerSpectrum! Character limit (" << (MAX_LINESIZE-1) << "/line) exceeded in file " << filename << "." << endl;
				fclose(pkfile);
				parallel.abortForce();
			}
			
			if (sscanf(line, " %lf %lf", &dummy1, &dummy2) == 2 && !feof(pkfile) && !ferror(pkfile)) numpoints++;
		}
		
		if (numpoints < 2)
		{
			cerr << " proc#" << parallel.rank() << ": error in loadPowerSpectrum! No valid data found in file " << filename << "." << endl;
			fclose(pkfile);
			parallel.abortForce();
		}
		
		k = (double *) malloc(sizeof(double) * numpoints);
		pk = (double *) malloc(sizeof(double) * numpoints);
		
		if (k == NULL || pk == NULL)
		{
			cerr << " proc#" << parallel.rank() << ": error in loadPowerSpectrum! Memory error." << endl;
			fclose(pkfile);
			parallel.abortForce();
		}
		
		rewind(pkfile);
		
		while (!feof(pkfile) && !ferror(pkfile))
		{
			fgets(line, MAX_LINESIZE, pkfile);
			
			if (sscanf(line, " %lf %lf", &dummy1, &dummy2) == 2 && !feof(pkfile) && !ferror(pkfile))
			{
				if (dummy1 < 0. || dummy2 < 0.)
				{
					cerr << " proc#" << parallel.rank() << ": error in loadPowerSpectrum! Negative entry encountered." << endl;
					free(k);
					free(pk);
					fclose(pkfile);
					parallel.abortForce();
				}
				
				if (i > 0)
				{
					if (k[i-1] >= dummy1 * boxsize)
					{
						cerr << " proc#" << parallel.rank() << ": error in loadPowerSpectrum! k-values are not strictly ordered." << endl;
						free(k);
						free(pk);
						fclose(pkfile);
						parallel.abortForce();
					}
				}
				
				k[i] = dummy1 * boxsize;
				pk[i] = sqrt(0.5 * dummy2 * boxsize);
				i++;
			}
		}
		
		fclose(pkfile);
		
		if (i != numpoints)
		{
			cerr << " proc#" << parallel.rank() << ": error in loadPowerSpectrum! File may have changed or file pointer corrupted." << endl;
			free(k);
			free(pk);
			parallel.abortForce();
		}
		
		parallel.broadcast_dim0<int>(numpoints, 0);
	}
	else
	{
		parallel.broadcast_dim0<int>(numpoints, 0);
		
		if (numpoints < 2)
		{
			cerr << " proc#" << parallel.rank() << ": error in loadPowerSpectrum! Communication error." << endl;
			parallel.abortForce();
		}
		
		k = (double *) malloc(sizeof(double) * numpoints);
		pk = (double *) malloc(sizeof(double) * numpoints);
		
		if (k == NULL || pk == NULL)
		{
			cerr << " proc#" << parallel.rank() << ": error in loadPowerSpectrum! Memory error." << endl;
			parallel.abortForce();
		}
	}
	
	parallel.broadcast_dim0<double>(k, numpoints, 0);
	parallel.broadcast_dim0<double>(pk, numpoints, 0);
	
	pkspline = gsl_spline_alloc(gsl_interp_cspline, numpoints);
	
	gsl_spline_init(pkspline, k, pk, numpoints);
	
	free(k);
	free(pk);
}


//////////////////////////
// loadTransferFunctions (1)
//////////////////////////
// Description:
//   loads a set of tabulated transfer functions from a file
// 
// Arguments:
//   filename   string containing the path to the template file
//   tk_delta   will point to the gsl_spline which holds the tabulated
//              transfer function for delta (memory will be allocated)
//   tk_theta   will point to the gsl_spline which holds the tabulated
//              transfer function for theta (memory will be allocated)
//   qname      string containing the name of the component (e.g. "cdm")
//   boxsize    comoving box size (in the same units as used in the file)
//   h          conversion factor between 1/Mpc and h/Mpc (theta is in units of 1/Mpc)
//
// Returns:
// 
//////////////////////////

void loadTransferFunctions(const char * filename, gsl_spline * & tk_delta, gsl_spline * & tk_theta, const char * qname, const double boxsize, const double h)
{
	int i = 0, numpoints = 0;
	double * k;
	double * tk_d;
	double * tk_t;

	if (parallel.grid_rank()[0] == 0) // read file
	{
		FILE * tkfile;
		char line[MAX_LINESIZE];
		char format[MAX_LINESIZE];
		char * ptr;
		double dummy[3];
		int kcol = -1, dcol = -1, tcol = -1, colmax;
		
		line[MAX_LINESIZE-1] = 0;
		
		tkfile = fopen(filename, "r");
		
		if (tkfile == NULL)
		{
			cerr << " proc#" << parallel.rank() << ": error in loadTransferFunctions! Unable to open file " << filename << "." << endl;
			parallel.abortForce();
		}
		
		while (!feof(tkfile) && !ferror(tkfile))
		{
			fgets(line, MAX_LINESIZE, tkfile);
			if (line[MAX_LINESIZE-1] != 0)
			{
				cerr << " proc#" << parallel.rank() << ": error in loadTransferFunctions! Character limit (" << (MAX_LINESIZE-1) << "/line) exceeded in file " << filename << "." << endl;
				fclose(tkfile);
				parallel.abortForce();
			}
			
			if (line[0] != '#' && !feof(tkfile) && !ferror(tkfile)) numpoints++;
		}
		
		if (numpoints < 2)
		{
			cerr << " proc#" << parallel.rank() << ": error in loadTransferFunctions! No valid data found in file " << filename << "." << endl;
			fclose(tkfile);
			parallel.abortForce();
		}
		
		k = (double *) malloc(sizeof(double) * numpoints);
		tk_d = (double *) malloc(sizeof(double) * numpoints);
		tk_t = (double *) malloc(sizeof(double) * numpoints);
		
		if (k == NULL || tk_d == NULL || tk_t == NULL)
		{
			cerr << " proc#" << parallel.rank() << ": error in loadTransferFunctions! Memory error." << endl;
			fclose(tkfile);
			parallel.abortForce();
		}
		
		rewind(tkfile);
		
		while (!feof(tkfile) && !ferror(tkfile))
		{
			fgets(line, MAX_LINESIZE, tkfile);
			for (ptr = line, i = 0; (ptr = strchr(ptr, ':')) != NULL; i++)
			{
				ptr++;
				if (*ptr == 'k') kcol = i;
				else if (*ptr == 'd')
				{
					if (strncmp(ptr+2, qname, strlen(qname)) == 0) dcol = i;
				}
				else if (*ptr == 't')
				{
					if (strncmp(ptr+2, qname, strlen(qname)) == 0) tcol = i;
				}
			}
			
			if (kcol >= 0 && dcol >= 0 && tcol >= 0) break;
		}
		
		if (kcol < 0 || dcol < 0 || tcol < 0)
		{
			cerr << " proc#" << parallel.rank() << ": error in loadTransferFunctions! Unable to identify requested columns!" << endl;
			fclose(tkfile);
			free(k);
			free(tk_d);
			free(tk_t);
			parallel.abortForce();
		}
		
		colmax = i;
		for (i = 0, ptr=format; i < colmax; i++)
		{
		    if (i == kcol || i == dcol || i == tcol)
		    {
		        strncpy(ptr, " %lf", 4);
		        ptr += 4;
		    }
		    else
		    {
		        strncpy(ptr, " %*lf", 5);
		        ptr += 5;
		    }
		}
		*ptr = '\0';
		
		if (kcol < dcol && dcol < tcol)
		{
			kcol = 0; dcol = 1; tcol = 2;
		}
		else if (kcol < tcol && tcol < dcol)
		{
			kcol = 0; dcol = 2; tcol = 1;
		}
		else if (dcol < kcol && kcol < tcol)
		{
			kcol = 1; dcol = 0; tcol = 2;
		}
		else if (dcol < tcol && tcol < kcol)
		{
			kcol = 2; dcol = 0; tcol = 1;
		}
		else if (tcol < kcol && kcol < dcol)
		{
			kcol = 1; dcol = 2; tcol = 0;
		}
		else if (tcol < dcol && dcol < kcol)
		{
			kcol = 2; dcol = 1; tcol = 0;
		}
		else
		{
			cerr << " proc#" << parallel.rank() << ": error in loadTransferFunctions! Inconsistent columns!" << endl;
			fclose(tkfile);
			free(k);
			free(tk_d);
			free(tk_t);
			parallel.abortForce();
		}
		
		i = 0;
		while (!feof(tkfile) && !ferror(tkfile))
		{
			fgets(line, MAX_LINESIZE, tkfile);
			
			if (sscanf(line, format, dummy, dummy+1, dummy+2) == 3 && !feof(tkfile) && !ferror(tkfile))
			{
				if (dummy[kcol] < 0.)
				{
					cerr << " proc#" << parallel.rank() << ": error in loadTransferFunctions! Negative k-value encountered." << endl;
					free(k);
					free(tk_d);
					free(tk_t);
					fclose(tkfile);
					parallel.abortForce();
				}
				
				if (i > 0)
				{
					if (k[i-1] >= dummy[kcol] * boxsize)
					{
						cerr << " proc#" << parallel.rank() << ": error in loadTransferFunctions! k-values are not strictly ordered." << endl;
						free(k);
						free(tk_d);
						free(tk_t);
						fclose(tkfile);
						parallel.abortForce();
					}
				}
				
				k[i] = dummy[kcol] * boxsize;
				tk_d[i] = dummy[dcol];
				tk_t[i] = dummy[tcol] * boxsize / h;
				i++;
			}
		}
		
		fclose(tkfile);
		
		if (i != numpoints)
		{
			cerr << " proc#" << parallel.rank() << ": error in loadTransferFunctions! File may have changed or file pointer corrupted." << endl;
			free(k);
			free(tk_d);
			free(tk_t);
			parallel.abortForce();
		}
		
		parallel.broadcast_dim0<int>(numpoints, 0);
	}
	else
	{
		parallel.broadcast_dim0<int>(numpoints, 0);
		
		if (numpoints < 2)
		{
			cerr << " proc#" << parallel.rank() << ": error in loadTransferFunctions! Communication error." << endl;
			parallel.abortForce();
		}
		
		k = (double *) malloc(sizeof(double) * numpoints);
		tk_d = (double *) malloc(sizeof(double) * numpoints);
		tk_t = (double *) malloc(sizeof(double) * numpoints);
		
		if (k == NULL || tk_d == NULL || tk_t == NULL)
		{
			cerr << " proc#" << parallel.rank() << ": error in loadTransferFunctions! Memory error." << endl;
			parallel.abortForce();
		}
	}
	
	parallel.broadcast_dim0<double>(k, numpoints, 0);
	parallel.broadcast_dim0<double>(tk_d, numpoints, 0);
	parallel.broadcast_dim0<double>(tk_t, numpoints, 0);
	
	tk_delta = gsl_spline_alloc(gsl_interp_cspline, numpoints);
	tk_theta = gsl_spline_alloc(gsl_interp_cspline, numpoints);
	
	gsl_spline_init(tk_delta, k, tk_d, numpoints);
	gsl_spline_init(tk_theta, k, tk_t, numpoints);
	
	free(k);
	free(tk_d);
	free(tk_t);
}


//////////////////////////
// generateCICKernel
//////////////////////////
// Description:
//   generates convolution kernel for CIC projection
// 
// Arguments:
//   ker        reference to allocated field that will contain the convolution kernel
//   numpcl     number of particles in the pcldata pointer
//   pcldata    raw particle data from which the kernel will be constructed;
//              a standard kernel will be provided in case no particles are specified
//   numtile    tiling factor used for the particle template
//
// Returns:
// 
//////////////////////////

void generateCICKernel(Field<Real> & ker, const long numpcl = 0, float * pcldata = NULL, const int numtile = 1)
{
	const long linesize = ker.lattice().sizeLocal(0);
	Real renorm = linesize * linesize;
	long i, oct, sx, sy, sz;
	float wx, wy, wz, q1, q2, q3, q4, ww;
	Site x(ker.lattice());
	
	for (x.first(); x.test(); x.next())
	{
		ker(x) = 0.;
	}
	
	if (numpcl == 0 || pcldata == NULL) // standard kernel
	{
		if (x.setCoord(0, 0, 0))
		{
			ker(x) = 6. * renorm;
			ker(x+0) = -renorm;
			x.setCoord(linesize-1, 0, 0);
			ker(x) = -renorm;
		}
	
		if (x.setCoord(0, 1, 0))
			ker(x) = -renorm;
	
		if (x.setCoord(0, 0, 1))
			ker(x) = -renorm;
	
		if (x.setCoord(0, linesize-1, 0))
			ker(x) = -renorm;
	
		if (x.setCoord(0, 0, linesize-1))
			ker(x) = -renorm;
		
		return;
	}
	
	// compute kernel explicitly
	
	renorm /= (Real) (numpcl * (long) numtile * (long) numtile * (long) numtile) / (Real) (linesize * linesize * linesize);
	
	for (i = 0; i < numpcl; i++)
	{
		for (oct = 0; oct < 8; oct++)
		{
			if (oct == 0)
			{
				if (linesize * pcldata[3*i] / numtile >= 1. || linesize * pcldata[3*i+1] / numtile >= 1. || linesize * pcldata[3*i+2] / numtile >= 1.)
					continue;
				
				// particle is in the first octant   
				wx = 1. - linesize * pcldata[3*i] / numtile;
				wy = 1. - linesize * pcldata[3*i+1] / numtile;
				wz = 1. - linesize * pcldata[3*i+2] / numtile;
				
				sx = 1;
				sy = 1;
				sz = 1;
			}
			else if (oct == 1)
			{
				if (linesize * (1.-pcldata[3*i]) / numtile >= 1. || linesize * pcldata[3*i+1] / numtile >= 1. || linesize * pcldata[3*i+2] / numtile >= 1.)
					continue;
					
				// particle is in the second octant   
				wx = 1. - linesize * (1.-pcldata[3*i]) / numtile;
				wy = 1. - linesize * pcldata[3*i+1] / numtile;
				wz = 1. - linesize * pcldata[3*i+2] / numtile;
				
				sx = -1;
				sy = 1;
				sz = 1;
			}
			else if (oct == 2)
			{
				if (linesize * (1.-pcldata[3*i]) / numtile >= 1. || linesize * (1.-pcldata[3*i+1]) / numtile >= 1. || linesize * pcldata[3*i+2] / numtile >= 1.)
					continue;
					
				// particle is in the third octant   
				wx = 1. - linesize * (1.-pcldata[3*i]) / numtile;
				wy = 1. - linesize * (1.-pcldata[3*i+1]) / numtile;
				wz = 1. - linesize * pcldata[3*i+2] / numtile;
				
				sx = -1;
				sy = -1;
				sz = 1;
			}
			else if (oct == 3)
			{
				if (linesize * pcldata[3*i] / numtile >= 1. || linesize * (1.-pcldata[3*i+1]) / numtile >= 1. || linesize * pcldata[3*i+2] / numtile >= 1.)
					continue;
					
				// particle is in the fourth octant   
				wx = 1. - linesize * pcldata[3*i] / numtile;
				wy = 1. - linesize * (1.-pcldata[3*i+1]) / numtile;
				wz = 1. - linesize * pcldata[3*i+2] / numtile;
				
				sx = 1;
				sy = -1;
				sz = 1;
			}
			else if (oct == 4)
			{
				if (linesize * pcldata[3*i] / numtile >= 1. || linesize * pcldata[3*i+1] / numtile >= 1. || linesize * (1.-pcldata[3*i+2]) / numtile >= 1.)
					continue;
				
				// particle is in the fifth octant   
				wx = 1. - linesize * pcldata[3*i] / numtile;
				wy = 1. - linesize * pcldata[3*i+1] / numtile;
				wz = 1. - linesize * (1.-pcldata[3*i+2]) / numtile;
				
				sx = 1;
				sy = 1;
				sz = -1;
			}
			else if (oct == 5)
			{
				if (linesize * (1.-pcldata[3*i]) / numtile >= 1. || linesize * pcldata[3*i+1] / numtile >= 1. || linesize * (1.-pcldata[3*i+2]) / numtile >= 1.)
					continue;
					
				// particle is in the sixth octant   
				wx = 1. - linesize * (1.-pcldata[3*i]) / numtile;
				wy = 1. - linesize * pcldata[3*i+1] / numtile;
				wz = 1. - linesize * (1.-pcldata[3*i+2]) / numtile;
				
				sx = -1;
				sy = 1;
				sz = -1;
			}
			else if (oct == 6)
			{
				if (linesize * (1.-pcldata[3*i]) / numtile >= 1. || linesize * (1.-pcldata[3*i+1]) / numtile >= 1. || linesize * (1.-pcldata[3*i+2]) / numtile >= 1.)
					continue;
					
				// particle is in the seventh octant   
				wx = 1. - linesize * (1.-pcldata[3*i]) / numtile;
				wy = 1. - linesize * (1.-pcldata[3*i+1]) / numtile;
				wz = 1. - linesize * (1.-pcldata[3*i+2]) / numtile;
				
				sx = -1;
				sy = -1;
				sz = -1;
			}
			else if (oct == 7)
			{
				if (linesize * pcldata[3*i] / numtile >= 1. || linesize * (1.-pcldata[3*i+1]) / numtile >= 1. || linesize * (1.-pcldata[3*i+2]) / numtile >= 1.)
					continue;
					
				// particle is in the eight-th octant   
				wx = 1. - linesize * pcldata[3*i] / numtile;
				wy = 1. - linesize * (1.-pcldata[3*i+1]) / numtile;
				wz = 1. - linesize * (1.-pcldata[3*i+2]) / numtile;
				
				sx = 1;
				sy = -1;
				sz = -1;
			}
			else
				continue;
			
			// 0-direction
			
			ww = wy*wz*renorm;		   
			q1 = (wx > 0.9) ? 2. : 1.;
			
			if (x.setCoord(0, 0, 0))
			{
				ker(x) += ww * wy * wz * q1;
				x.setCoord((linesize+sx)%linesize, 0, 0);
				ker(x) -= ww * wy * wz;
				if (wx > 0.9)
				{
					x.setCoord((linesize-sx)%linesize, 0, 0);
					ker(x) -= ww * wy * wz;
				}
			}
			
			if (x.setCoord(0, 0, (linesize+sz)%linesize))
			{
				ker(x) += ww * wy * (1.-wz) * q1;
				x.setCoord((linesize+sx)%linesize, 0, (linesize+sz)%linesize);
				ker(x) -= ww * wy * (1.-wz);
				if (wx > 0.9)
				{
					x.setCoord((linesize-sx)%linesize, 0, (linesize+sz)%linesize);
					ker(x) -= ww * wy * (1.-wz);
				}
			}
			
			if (x.setCoord(0, (linesize+sy)%linesize, 0))
			{
				ker(x) += ww * (1.-wy) * wz * q1;
				x.setCoord((linesize+sx)%linesize, (linesize+sy)%linesize, 0);
				ker(x) -= ww * (1.-wy) * wz;
				if (wx > 0.9)
				{
					x.setCoord((linesize-sx)%linesize, (linesize+sy)%linesize, 0);
					ker(x) -= ww * (1.-wy) * wz;
				}
			}
			
			if (x.setCoord(0, (linesize+sy)%linesize, (linesize+sz)%linesize))
			{
				ker(x) += ww * (1.-wy) * (1.-wz) * q1;
				x.setCoord((linesize+sx)%linesize, (linesize+sy)%linesize, (linesize+sz)%linesize);
				ker(x) -= ww * (1.-wy) * (1.-wz);
				if (wx > 0.9)
				{
					x.setCoord((linesize-sx)%linesize, (linesize+sy)%linesize, (linesize+sz)%linesize);
					ker(x) -= ww * (1.-wy) * (1.-wz);
				}
			}
			
			// 1-direction
			
			ww = wx*wz*renorm;
			q1 = (wy > 0.9) ? 2. : 1.;
			
			if (x.setCoord(0, 0, 0))
			{
				ker(x) += ww * wx * wz * q1;
				x.setCoord((linesize+sx)%linesize, 0, 0);
				ker(x) += ww * (1.-wx) * wz * q1;
			}
			
			if (x.setCoord(0, 0, (linesize+sz)%linesize))
			{
				ker(x) += ww * wx * (1.-wz) * q1;
				x.setCoord((linesize+sx)%linesize, 0, (linesize+sz)%linesize);
				ker(x) += ww * (1.-wx) * (1.-wz) * q1;
			}
			
			if (x.setCoord(0, (linesize+sy)%linesize, 0))
			{
				ker(x) -= ww * wx * wz;
				x.setCoord((linesize+sx)%linesize, (linesize+sy)%linesize, 0);
				ker(x) -= ww * (1.-wx) * wz;
			}
			
			if (x.setCoord(0, (linesize+sy)%linesize, (linesize+sz)%linesize))
			{
				ker(x) -= ww * wx * (1.-wz);
				x.setCoord((linesize+sx)%linesize, (linesize+sy)%linesize, (linesize+sz)%linesize);
				ker(x) -= ww * (1.-wx) * (1.-wz);
			}
			
			if (wy > 0.9)
			{
				if (x.setCoord(0, (linesize-sy)%linesize, 0))
				{
					ker(x) -= ww * wx * wz;
					x.setCoord((linesize+sx)%linesize, (linesize-sy)%linesize, 0);
					ker(x) -= ww * (1.-wx) * wz;
				}
				
				if (x.setCoord(0, (linesize-sy)%linesize, (linesize+sz)%linesize))
				{
					ker(x) -= ww * wx * (1.-wz);
					x.setCoord((linesize+sx)%linesize, (linesize-sy)%linesize, (linesize+sz)%linesize);
					ker(x) -= ww * (1.-wx) * (1.-wz);
				}
			}
						
			// 2-direction
			
			ww = wx*wy*renorm;
			q1 = (wz > 0.9) ? 2. : 1.;
			
			if (x.setCoord(0, 0, 0))
			{
				ker(x) += ww * wx * wy * q1;
				x.setCoord((linesize+sx)%linesize, 0, 0);
				ker(x) += ww * (1.-wx) * wy * q1;
			}
			
			if (x.setCoord(0, (linesize+sy)%linesize, 0))
			{
				ker(x) += ww * wx * (1.-wy) * q1;
				x.setCoord((linesize+sx)%linesize, (linesize+sy)%linesize, 0);
				ker(x) += ww * (1.-wx) * (1.-wy) * q1;
			}
			
			if (x.setCoord(0, 0, (linesize+sz)%linesize))
			{
				ker(x) -= ww * wx * wy;
				x.setCoord((linesize+sx)%linesize, 0, (linesize+sz)%linesize);
				ker(x) -= ww * (1.-wx) * wy;
			}
			
			if (x.setCoord(0, (linesize+sy)%linesize, (linesize+sz)%linesize))
			{
				ker(x) -= ww * wx * (1.-wy);
				x.setCoord((linesize+sx)%linesize, (linesize+sy)%linesize, (linesize+sz)%linesize);
				ker(x) -= ww * (1.-wx) * (1.-wy);
			}
			
			if (wz > 0.9)
			{
				if (x.setCoord(0, 0, (linesize-sz)%linesize))
				{
					ker(x) -= ww * wx * wy;
					x.setCoord((linesize+sx)%linesize, 0, (linesize-sz)%linesize);
					ker(x) -= ww * (1.-wx) * wy;
				}
				
				if (x.setCoord(0, (linesize+sy)%linesize, (linesize-sz)%linesize))
				{
					ker(x) -= ww * wx * (1.-wy);
					x.setCoord((linesize+sx)%linesize, (linesize+sy)%linesize, (linesize-sz)%linesize);
					ker(x) -= ww * (1.-wx) * (1.-wy);
				}
			}
		}
	}
}


#ifdef FFT3D

//////////////////////////
// generateDisplacementField (generateRealization)
//////////////////////////
// Description:
//   generates particle displacement field
//
// Non-type template parameters:
//   ignorekernel  this is effectively an optimization flag defaulted to 0; instantiating with 1 instead will cause
//                 the function to ignore the convolution kernel, allowing the function to be used for generating
//                 realizations (generateRealization is simply an alias for generateDisplacementField<1>)
// 
// Arguments:
//   potFT         reference to allocated field that contains the convolution kernel relating the potential
//                 (generating the displacement field) with the bare density perturbation; will contain the
//                 Fourier image of the potential generating the displacement field
//   coeff         gauge correction coefficient "H_conformal^2"
//   pkspline      pointer to a gsl_spline which holds a tabulated power spectrum
//   seed          initial seed for random number generator
//   ksphere       flag to indicate that only a sphere in k-space should be initialized
//                 (default = 0: full k-space cube is initialized)
//   deconvolve_f  flag to indicate deconvolution function
//                 0: no deconvolution
//                 1: sinc (default)
//
// Returns:
// 
//////////////////////////

#ifndef generateRealization
#define generateRealization generateDisplacementField<1>
#endif

template<int ignorekernel = 0>
void generateDisplacementField(Field<Cplx> & potFT, const Real coeff, const gsl_spline * pkspline, const unsigned int seed, const int ksphere = 0, const int deconvolve_f = 1)
{
	const int linesize = potFT.lattice().size(1);
	const int kmax = (linesize / 2) - 1;
	rKSite k(potFT.lattice());
	int kx, ky, kz, i, j;
	int kymin, kymax, kzmin, kzmax;
	long jumpy, jumpz;
	float r1, r2, k2, s;
	float * sinc;
	sitmo::prng_engine prng;
	uint64_t huge_skip = HUGE_SKIP;
	gsl_interp_accel * acc = gsl_interp_accel_alloc();
	
	sinc = (float *) malloc(linesize * sizeof(float));
	
	sinc[0] = 1.;
	if (deconvolve_f == 1)
	{
		for (i = 1; i < linesize; i++)
			sinc[i] = sin(M_PI * (float) i / (float) linesize) * (float) linesize / (M_PI * (float) i);
	}
	else
	{
		for (i = 1; i < linesize; i++)
			sinc[i] = 1.;
	}
	
	k.initialize(potFT.lattice(), potFT.lattice().siteLast());
	kymax = k.coord(1);
	kzmax = k.coord(2);
	k.initialize(potFT.lattice(), potFT.lattice().siteFirst());
	kymin = k.coord(1);
	kzmin = k.coord(2);
		
	if (kymin < (linesize / 2) + 1 && kzmin < (linesize / 2) + 1)
	{
		prng.seed(seed);
		   
		if (kymin == 0 && kzmin == 0)
		{
			k.setCoord(0, 0, 0);
			potFT(k) = Cplx(0.,0.);
			kx = 1;
		}
		else
		{
			kx = 0;
			prng.discard(((uint64_t) kzmin * huge_skip + (uint64_t) kymin) * huge_skip); 
		}
		
		for (kz = kzmin; kz < (linesize / 2) + 1 && kz <= kzmax; kz++)
		{
			for (ky = kymin, j = 0; ky < (linesize / 2) + 1 && ky <= kymax; ky++, j++)
			{
				for (i = 0; kx < (linesize / 2) + 1; kx++)
				{
					k.setCoord(kx, ky, kz);
					
					k2 = (float) (kx * kx) + (float) (ky * ky) + (float) (kz * kz);
					
					if (kx >= kmax || ky >= kmax || kz >= kmax || (k2 >= kmax * kmax && ksphere > 0))
					{
						potFT(k) = Cplx(0., 0.);
					}
					else
					{
						s = sinc[kx] * sinc[ky] * sinc[kz];
						k2 *= 4. * M_PI * M_PI;
						do
						{
							r1 = (float) prng() / (float) sitmo::prng_engine::max();
							i++;
						}
						while (r1 == 0.);
						r2 = (float) prng() / (float) sitmo::prng_engine::max();
						i++;
					
						potFT(k) = (ignorekernel ? Cplx(cos(2. * M_PI * r2), sin(2. * M_PI * r2)) : Cplx(cos(2. * M_PI * r2), sin(2. * M_PI * r2)) * (1. + 7.5 * coeff / k2) / potFT(k)) * sqrt(-2. * log(r1)) * gsl_spline_eval(pkspline, sqrt(k2), acc) * s;
					}
				}
				prng.discard(huge_skip - (uint64_t) i);
				kx = 0;
			}
			prng.discard(huge_skip * (huge_skip - (uint64_t) j));
		}
	}
	
	if (kymax >= (linesize / 2) + 1 && kzmin < (linesize / 2) + 1)
	{
		prng.seed(seed);
		prng.discard(((huge_skip + (uint64_t) kzmin) * huge_skip + (uint64_t) (linesize - kymax)) * huge_skip);
		
		for (kz = kzmin; kz < (linesize / 2) + 1 && kz <= kzmax; kz++)
		{
			for (ky = kymax, j = 0; ky >= (linesize / 2) + 1 && ky >= kymin; ky--, j++)
			{
				for (kx = 0, i = 0; kx < (linesize / 2) + 1; kx++)
				{
					k.setCoord(kx, ky, kz);
										
					k2 = (float) (kx * kx) + (float) ((linesize-ky) * (linesize-ky)) + (float) (kz * kz);
					
					if (kx >= kmax || (linesize-ky) >= kmax || kz >= kmax || (k2 >= kmax * kmax && ksphere > 0))
					{
						potFT(k) = Cplx(0., 0.);
					}
					else
					{
						s = sinc[kx] * sinc[linesize-ky] * sinc[kz];
						k2 *= 4. * M_PI * M_PI;
						do
						{
							r1 = (float) prng() / (float) sitmo::prng_engine::max();
							i++;
						}
						while (r1 == 0.);
						r2 = (float) prng() / (float) sitmo::prng_engine::max();
						i++;
						
						potFT(k) = (ignorekernel ? Cplx(cos(2. * M_PI * r2), sin(2. * M_PI * r2)) : Cplx(cos(2. * M_PI * r2), sin(2. * M_PI * r2)) * (1. + 7.5 * coeff / k2) / potFT(k)) * sqrt(-2. * log(r1)) * gsl_spline_eval(pkspline, sqrt(k2), acc) * s;
					}
				}
				prng.discard(huge_skip - (uint64_t) i);
			}
			prng.discard(huge_skip * (huge_skip - (uint64_t) j));
		}
	}
	
	if (kymin < (linesize / 2) + 1 && kzmax >= (linesize / 2) + 1)
	{
		prng.seed(seed);
		prng.discard(((huge_skip + huge_skip + (uint64_t) (linesize - kzmax)) * huge_skip + (uint64_t) kymin) * huge_skip);
		
		for (kz = kzmax; kz >= (linesize / 2) + 1 && kz >= kzmin; kz--)
		{
			for (ky = kymin, j = 0; ky < (linesize / 2) + 1 && ky <= kymax; ky++, j++)
			{
				for (kx = 1, i = 0; kx < (linesize / 2) + 1; kx++)
				{
					k.setCoord(kx, ky, kz);
					
					k2 = (float) (kx * kx) + (float) (ky * ky) + (float) ((linesize-kz) * (linesize-kz));
					
					if (kx >= kmax || ky >= kmax || (linesize-kz) >= kmax || (k2 >= kmax * kmax && ksphere > 0))
					{
						potFT(k) = Cplx(0., 0.);
					}
					else
					{
						s = sinc[kx] * sinc[ky] * sinc[linesize-kz];
						k2 *= 4. * M_PI * M_PI;
						do
						{
							r1 = (float) prng() / (float) sitmo::prng_engine::max();
							i++;
						}
						while (r1 == 0.);
						r2 = (float) prng() / (float) sitmo::prng_engine::max();
						i++;
					
						potFT(k) = (ignorekernel ? Cplx(cos(2. * M_PI * r2), sin(2. * M_PI * r2)) : Cplx(cos(2. * M_PI * r2), sin(2. * M_PI * r2)) * (1. + 7.5 * coeff / k2) / potFT(k)) * sqrt(-2. * log(r1)) * gsl_spline_eval(pkspline, sqrt(k2), acc) * s;
					}
				}
				prng.discard(huge_skip - (uint64_t) i);
			}
			prng.discard(huge_skip * (huge_skip - (uint64_t) j));
		}
		
		prng.seed(seed);
		prng.discard(((uint64_t) (linesize - kzmax) * huge_skip + (uint64_t) kymin) * huge_skip);
		kx = 0;
		
		for (kz = kzmax; kz >= (linesize / 2) + 1 && kz >= kzmin; kz--)
		{
			for (ky = kymin, j = 0; ky < (linesize / 2) + 1 && ky <= kymax; ky++, j++)
			{
				k.setCoord(kx, ky, kz);
					
				k2 = (float) (ky * ky) + (float) ((linesize-kz) * (linesize-kz));
				i = 0;
				
				if (ky >= kmax || (linesize-kz) >= kmax || (k2 >= kmax * kmax && ksphere > 0))
				{
					potFT(k) = Cplx(0., 0.);
				}
				else
				{
					s = sinc[ky] * sinc[linesize-kz];
					k2 *= 4. * M_PI * M_PI;
					do
					{
						r1 = (float) prng() / (float) sitmo::prng_engine::max();
						i++;
					}
					while (r1 == 0.);
					r2 = (float) prng() / (float) sitmo::prng_engine::max();
					i++;
				
					potFT(k) = (ignorekernel? Cplx(cos(2. * M_PI * r2), -sin(2. * M_PI * r2)) : Cplx(cos(2. * M_PI * r2), -sin(2. * M_PI * r2)) * (1. + 7.5 * coeff / k2) / potFT(k)) * sqrt(-2. * log(r1)) * gsl_spline_eval(pkspline, sqrt(k2), acc) * s;
				}
								
				prng.discard(huge_skip - (uint64_t) i);
			}
			prng.discard(huge_skip * (huge_skip - (uint64_t) j));
		}
	}
	
	if (kymax >= (linesize / 2) + 1 && kzmax >= (linesize / 2) + 1)
	{
		prng.seed(seed);
		prng.discard(((huge_skip + huge_skip + huge_skip + (uint64_t) (linesize - kzmax)) * huge_skip + (uint64_t) (linesize - kymax)) * huge_skip);
		
		for (kz = kzmax; kz >= (linesize / 2) + 1 && kz >= kzmin; kz--)
		{
			for (ky = kymax, j = 0; ky >= (linesize / 2) + 1 && ky >= kymin; ky--, j++)
			{
				for (kx = 1, i = 0; kx < (linesize / 2) + 1; kx++)
				{
					k.setCoord(kx, ky, kz);
					
					k2 = (float) (kx * kx) + (float) ((linesize-ky) * (linesize-ky)) + (float) ((linesize-kz) * (linesize-kz));
					
					if (kx >= kmax || (linesize-ky) >= kmax || (linesize-kz) >= kmax || (k2 >= kmax * kmax && ksphere > 0))
					{
						potFT(k) = Cplx(0., 0.);
					}
					else
					{
						s = sinc[kx] * sinc[linesize-ky] * sinc[linesize-kz];
						k2 *= 4. * M_PI * M_PI;
						do
						{
							r1 = (float) prng() / (float) sitmo::prng_engine::max();
							i++;
						}
						while (r1 == 0.);
						r2 = (float) prng() / (float) sitmo::prng_engine::max();
						i++;
					
						potFT(k) = (ignorekernel ? Cplx(cos(2. * M_PI * r2), sin(2. * M_PI * r2)) : Cplx(cos(2. * M_PI * r2), sin(2. * M_PI * r2)) * (1. + 7.5 * coeff / k2) / potFT(k)) * sqrt(-2. * log(r1)) * gsl_spline_eval(pkspline, sqrt(k2), acc) * s;
					}
				}
				prng.discard(huge_skip - (uint64_t) i);
			}
			prng.discard(huge_skip * (huge_skip - (uint64_t) j));
		}
		
		prng.seed(seed);
		prng.discard(((huge_skip + huge_skip + (uint64_t) (linesize - kzmax)) * huge_skip + (uint64_t) (linesize - kymax)) * huge_skip);
		kx = 0;
		
		for (kz = kzmax; kz >= (linesize / 2) + 1 && kz >= kzmin; kz--)
		{
			for (ky = kymax, j = 0; ky >= (linesize / 2) + 1 && ky >= kymin; ky--, j++)
			{
				k.setCoord(kx, ky, kz);
					
				k2 = (float) ((linesize-ky) * (linesize-ky)) + (float) ((linesize-kz) * (linesize-kz));
				i = 0;
				
				if ((linesize-ky) >= kmax || (linesize-kz) >= kmax || (k2 >= kmax * kmax && ksphere > 0))
				{
					potFT(k) = Cplx(0., 0.);
				}
				else
				{
					s = sinc[linesize-ky] * sinc[linesize-kz];
					k2 *= 4. * M_PI * M_PI;
					do
					{
						r1 = (float) prng() / (float) sitmo::prng_engine::max();
						i++;
					}
					while (r1 == 0.);
					r2 = (float) prng() / (float) sitmo::prng_engine::max();
					i++;
				
					potFT(k) = (ignorekernel ? Cplx(cos(2. * M_PI * r2), -sin(2. * M_PI * r2)) : Cplx(cos(2. * M_PI * r2), -sin(2. * M_PI * r2)) * (1. + 7.5 * coeff / k2) / potFT(k)) * sqrt(-2. * log(r1)) * gsl_spline_eval(pkspline, sqrt(k2), acc) * s;
				}
				
				prng.discard(huge_skip - (uint64_t) i);
			}
			prng.discard(huge_skip * (huge_skip - (uint64_t) j));
		}
	}
	
	gsl_interp_accel_free(acc);
	free(sinc);
}
#endif


//////////////////////////
// initializeParticlePositions
//////////////////////////
// Description:
//   initializes particle positions using a homogeneous template
// 
// Arguments:
//   numpart    number of particles of the template
//   partdata   particle positions in the template
//   numtile    tiling factor for homogeneous template - total particle number will be
//              numpart * numtile^3
//   pcls       reference to (empty) particle object which will contain the new particle ensemble
//
// Returns:
// 
//////////////////////////

void initializeParticlePositions(const long numpart, const float * partdata, const int numtile, Particles<part_simple,part_simple_info,part_simple_dataType> & pcls)
{
	long xtile, ytile, ztile, i;
	Site p(pcls.lattice());
	
	part_simple part;
	
	part.vel[0] = 0.;
	part.vel[1] = 0.;
	part.vel[2] = 0.;
	
	for (ztile = (pcls.lattice().coordSkip()[0] * numtile) / pcls.lattice().size(2); ztile <= ((pcls.lattice().coordSkip()[0] + pcls.lattice().sizeLocal(2)) * numtile) / pcls.lattice().size(2); ztile++)
	{
		if (ztile >= numtile) break;
		
		for (ytile = (pcls.lattice().coordSkip()[1] * numtile) / pcls.lattice().size(1); ytile <= ((pcls.lattice().coordSkip()[1] + pcls.lattice().sizeLocal(1)) * numtile) / pcls.lattice().size(1); ytile++)
		{
			if (ytile >= numtile) break;
			
			for (xtile = 0; xtile < numtile; xtile++)
			{
				for (i = 0; i < numpart; i++)
				{
					part.pos[0] = ((Real) xtile + partdata[3*i]) / (Real) numtile;
					part.pos[1] = ((Real) ytile + partdata[3*i+1]) / (Real) numtile;
					part.pos[2] = ((Real) ztile + partdata[3*i+2]) / (Real) numtile;
					
					part.ID = i + numpart * (xtile + (long) numtile * (ytile + (long) numtile * ztile));
					
					pcls.addParticle_global(part);
				}
			}
		}
	}
}


//////////////////////////
// applyMomentumDistribution
//////////////////////////
// Description:
//   adds a random momentum vector drawn from a Fermi-Dirac distribution to
//   each particle. The current implementation uses a simple rejection-sampling
//   method based on the ziggurat algorithm [G. Marsaglia and W.W. Tsang,
//   J. Stat. Softw. 5 (2000) 1]. The "ziggurat" is hard-coded and was precomputed
//   for the ultrarelativistic limit of a Fermi-Dirac distribution with zero
//   chemical potential. A method to construct the ziggurat on-the-fly for more
//   general distribution functions could be implemented in the future.
// 
// Arguments:
//   pcls   pointer to particle handler
//   seed   seed for random number generator
//   T_m    dimensionless parameter in the distribution function
//          in most cases the ratio of temperature and fundamental mass
//   delta  Field containing a local dT/T (optional)
//
// Returns: sum of momenta over all particles (for reduction)
// 
//////////////////////////

double applyMomentumDistribution(Particles<part_simple,part_simple_info,part_simple_dataType> * pcls, unsigned int seed, float T_m = 0., Field<Real> * delta = NULL)
{	
	Site xPart(pcls->lattice());
	Site x;	
	std::list<part_simple>::iterator it;
	sitmo::prng_engine prng;
	float r1, r2, q, dT;
	double l, dummy, d[3];
	uint32_t i, r;
	double sum_q = 0.0;
	
	const float ql[] = {0.0f,       0.0453329523f, 0.0851601009f, 0.115766097f,
	                    0.142169202f, 0.166069623f, 0.188283033f, 0.209275639f,
	                    0.229344099f, 0.248691555f, 0.267464827f, 0.285774552f,
	                    0.303706948f, 0.321331123f, 0.338703809f, 0.355872580f,
	                    0.372878086f, 0.389755674f, 0.406536572f, 0.423248791f,
	                    0.439917819f, 0.456567164f, 0.473218795f, 0.489893502f,
	                    0.506611198f, 0.523391180f, 0.540252353f, 0.557213439f,
	                    0.574293166f, 0.591510448f, 0.608884565f, 0.626435339f,
	                    0.644183319f, 0.662149971f, 0.680357893f, 0.698831041f,
	                    0.717594982f, 0.736677192f, 0.756107381f, 0.775917886f,
	                    0.796144127f, 0.816825143f, 0.838004248f, 0.859729814f,
	                    0.882056241f, 0.905045149f, 0.928766878f, 0.953302387f,
	                    0.978745698f, 1.00520708f,  1.03281729f,  1.06173322f,
	                    1.09214584f,  1.12429121f,  1.15846661f,  1.19505456f,
	                    1.23456031f,  1.27767280f,  1.32536981f,  1.37911302f,
	                    1.44124650f,  1.51592808f,  1.61180199f,  1.75307820f, 29.0f};
	                    
	const float qr[] = {29.0f,       11.8477879f, 10.3339062f, 9.58550750f,
	                    9.08034038f, 8.69584518f, 8.38367575f, 8.11975908f,
	                    7.89035001f, 7.68685171f, 7.50352431f, 7.33634187f,
	                    7.18236952f, 7.03940031f, 6.90573157f, 6.78002122f,
	                    6.66119176f, 6.54836430f, 6.44081188f, 6.33792577f,
	                    6.23919052f, 6.14416529f, 6.05246957f, 5.96377208f,
	                    5.87778212f, 5.79424263f, 5.71292467f, 5.63362289f,
	                    5.55615183f, 5.48034280f, 5.40604138f, 5.33310512f,
	                    5.26140176f, 5.19080749f, 5.12120558f, 5.05248501f,
	                    4.98453932f, 4.91726547f, 4.85056276f, 4.78433177f,
	                    4.71847331f, 4.65288728f, 4.58747144f, 4.52212011f,
	                    4.45672261f, 4.39116154f, 4.32531058f, 4.25903193f,
	                    4.19217309f, 4.12456265f, 4.05600493f, 3.98627278f,
	                    3.91509779f, 3.84215661f, 3.76705129f, 3.68928014f,
	                    3.60819270f, 3.52291720f, 3.43223655f, 3.33436084f,
	                    3.22646839f, 3.10364876f, 2.95592669f, 2.75624893f, 0.0f};
	                    
	const float f[] = {0.0f,          0.0010042516f, 0.0034718142f, 0.00631345904f,
	                  0.00938886471f, 0.0126471707f, 0.0160614805f, 0.0196150986f,
	                   0.0232967063f, 0.0270982040f, 0.0310135938f, 0.0350383391f,
	                   0.0391689707f, 0.0434028310f, 0.0477379005f, 0.0521726764f,
	                   0.0567060859f, 0.0613374223f, 0.0660662983f, 0.0708926105f,
	                   0.0758165136f, 0.0808384013f, 0.0859588926f, 0.0911788222f,
	                   0.0964992355f, 0.101921386f,  0.107446735f,  0.113076958f,
	                   0.118813945f,  0.124659815f,  0.130616922f,  0.136687871f,
	                   0.142875536f,  0.149183078f,  0.155613967f,  0.162172016f,
	                   0.168861408f,  0.175686736f,  0.182653051f,  0.189765914f,
	                   0.197031455f,  0.204456455f,  0.212048433f,  0.219815749f,
	                   0.227767740f,  0.235914877f,  0.244268957f,  0.252843349f,
	                   0.261653294f,  0.270716296f,  0.280052614f,  0.289685921f,
	                   0.299644172f,  0.309960782f,  0.320676286f,  0.331840691f,
	                   0.343516979f,  0.355786485f,  0.368757588f,  0.382580625f,
	                   0.397475563f,  0.413789108f,  0.432131941f,  0.453799050f, 0.482830296f};
	
	if (delta != NULL)
	{
		l = (double) delta->lattice().size(0);
		x.initialize(delta->lattice());
		x.first();
	}
	                    
	for (xPart.first(); xPart.test(); xPart.next())
	{
		if (pcls->field()(xPart).size != 0)
		{
			for (it=(pcls->field())(xPart).parts.begin(); it != (pcls->field())(xPart).parts.end(); ++it)
			{
				prng.seed(seed);
				prng.discard((uint64_t) (7l * (*it).ID));
				
				while (true)
				{
					r = prng();
					i = r % 64;
					r /= 64;
					
					q = ql[i] + 64.0f * (qr[i]-ql[i]) * ((float) r / (float) sitmo::prng_engine::max());
					
					if (q > ql[i+1] && q < qr[i+1]) break;
					
					if (f[i] + (f[i+1]-f[i]) * ((float) prng() / (float) sitmo::prng_engine::max()) < q * q / (exp(q) + 1.0f)) break;
				}
				
				r1 = acos(2. * ((float) prng() / (float) sitmo::prng_engine::max()) - 1.);
				r2 = 2 * M_PI * ((float) prng() / (float) sitmo::prng_engine::max());
				
				if (delta != NULL)
				{
					for (i = 0; i < 3; i++)
						d[i] = modf((*it).pos[i] * l, &dummy);
					
					dT = (1. - d[0]) * (1. - d[1]) * (1. - d[2]) * (*delta)(x);
					dT += d[0] * (1. - d[1]) * (1. - d[2]) * (*delta)(x+0);
					dT += (1. - d[0]) * d[1] * (1. - d[2]) * (*delta)(x+1);
					dT += d[0] * d[1] * (1. - d[2]) * (*delta)(x+0+1);
					dT += (1. - d[0]) * (1. - d[1]) * d[2] * (*delta)(x+2);
					dT += d[0] * (1. - d[1]) * d[2] * (*delta)(x+0+2);
					dT += (1. - d[0]) * d[1] * d[2] * (*delta)(x+1+2);
					dT += d[0] * d[1] * d[2] * (*delta)(x+2);
					
					q *= T_m * (1. + dT);
				}
				else
					q *= T_m;
				
				(*it).vel[0] += cos(r2) * sin(r1) * q;
				(*it).vel[1] += sin(r2) * sin(r1) * q;
				(*it).vel[2] += cos(r1) * q;
				
				sum_q += q;
			}
		}
		
		if (delta != NULL) x.next();
	}
	
	return sum_q;
}


#ifdef FFT3D

//////////////////////////
// generateIC_basic
//////////////////////////
// Description:
//   basic initial condition generator
// 
// Arguments:
//   sim            simulation metadata structure
//   ic             settings for IC generation
//   cosmo          cosmological parameter structure
//   fourpiG        4 pi G (in code units)
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

void generateIC_basic(metadata & sim, icsettings & ic, cosmology & cosmo, const double fourpiG, Particles<part_simple,part_simple_info,part_simple_dataType> * pcls_cdm, Particles<part_simple,part_simple_info,part_simple_dataType> * pcls_b, Particles<part_simple,part_simple_info,part_simple_dataType> * pcls_ncdm, double * maxvel, Field<Real> * phi, Field<Real> * chi, Field<Real> * Bi, Field<Real> * source, Field<Real> * Sij, Field<Cplx> * scalarFT, Field<Cplx> * BiFT, Field<Cplx> * SijFT, PlanFFT<Cplx> * plan_phi, PlanFFT<Cplx> * plan_chi, PlanFFT<Cplx> * plan_Bi, PlanFFT<Cplx> * plan_source, PlanFFT<Cplx> * plan_Sij, parameter * params, int & numparam)
{
	int i, j, p;
	double a = 1. / (1. + sim.z_in);
	float * pcldata = NULL;
	gsl_spline * pkspline = NULL;
	gsl_spline * nbspline = NULL;
	gsl_spline * vnbspline = NULL;
	gsl_spline * tk_d1 = NULL;
	gsl_spline * tk_d2 = NULL;
	gsl_spline * tk_t1 = NULL;
	gsl_spline * tk_t2 = NULL;
	double * temp1 = NULL;
	double * temp2 = NULL;
	Site x(phi->lattice());
	rKSite kFT(scalarFT->lattice());
	double max_displacement;
	double rescale;
	double mean_q;
	part_simple_info pcls_cdm_info;
	part_simple_dataType pcls_cdm_dataType;
	part_simple_info pcls_b_info;
	part_simple_dataType pcls_b_dataType;
	part_simple_info pcls_ncdm_info[MAX_PCL_SPECIES];
	part_simple_dataType pcls_ncdm_dataType;
	Real boxSize[3] = {1.,1.,1.};
	char ncdm_name[8];
	Field<Real> * ic_fields[2];
	
	ic_fields[0] = chi;
	ic_fields[1] = phi;

#ifdef HAVE_CLASS
  	background class_background;
  	thermo class_thermo;
  	perturbs class_perturbs;
#endif
	
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
	
	if (ic.pkfile[0] != '\0')	// initial displacements & velocities are derived from a single power spectrum
	{
		loadPowerSpectrum(ic.pkfile, pkspline, sim.boxsize);
	
		if (pkspline == NULL)
		{
			COUT << " error: power spectrum was empty!" << endl;
			parallel.abortForce();
		}
		
		temp1 = (double *) malloc(pkspline->size * sizeof(double));
		temp2 = (double *) malloc(pkspline->size * sizeof(double));
		
		for (i = 0; i < pkspline->size; i++)
		{
			temp1[i] = pkspline->x[i];
			temp2[i] = pkspline->y[i] / sim.boxsize / sim.boxsize;
		}
		gsl_spline_free(pkspline);
		pkspline = gsl_spline_alloc(gsl_interp_cspline, i);
		gsl_spline_init(pkspline, temp1, temp2, i);
	
		generateDisplacementField(*scalarFT, sim.gr_flag * Hconf(a, fourpiG, cosmo) * Hconf(a, fourpiG, cosmo), pkspline, (unsigned int) ic.seed, ic.flags & ICFLAG_KSPHERE);
	}
	else					// initial displacements and velocities are set by individual transfer functions
	{
#ifdef HAVE_CLASS
		if (ic.tkfile[0] == '\0')
		{
			initializeCLASSstructures(sim, ic, cosmo, class_background, class_thermo, class_perturbs, params, numparam);
			loadTransferFunctions(class_background, class_perturbs, tk_d1, tk_t1, "tot", sim.boxsize, sim.z_in, cosmo.h);
		}
		else
#endif
		loadTransferFunctions(ic.tkfile, tk_d1, tk_t1, "tot", sim.boxsize, cosmo.h);
		
		if (tk_d1 == NULL || tk_t1 == NULL)
		{
			COUT << " error: total transfer function was empty!" << endl;
			parallel.abortForce();
		}
		
		temp1 = (double *) malloc(tk_d1->size * sizeof(double));
		temp2 = (double *) malloc(tk_d1->size * sizeof(double));
		
		rescale = 3. * Hconf(a, fourpiG, cosmo) * Hconf(a, fourpiG, cosmo) * Hconf(a, fourpiG, cosmo) * (1. + 0.5 * Hconf(a, fourpiG, cosmo) * Hconf(a, fourpiG, cosmo) * ((1. / Hconf(0.98 * a, fourpiG, cosmo) / Hconf(0.98 * a, fourpiG, cosmo)) - (8. / Hconf(0.99 * a, fourpiG, cosmo) / Hconf(0.99 * a, fourpiG, cosmo)) + (8. / Hconf(1.01 * a, fourpiG, cosmo) / Hconf(1.01 * a, fourpiG, cosmo)) - (1. / Hconf(1.02 * a, fourpiG, cosmo) / Hconf(1.02 * a, fourpiG, cosmo))) / 0.12);
		for (i = 0; i < tk_d1->size; i++) // construct phi
			temp1[i] = (1.5 * (Hconf(a, fourpiG, cosmo) * Hconf(a, fourpiG, cosmo) - Hconf(1., fourpiG, cosmo) * Hconf(1., fourpiG, cosmo) * a * a * cosmo.Omega_Lambda) * tk_d1->y[i] + rescale * tk_t1->y[i] / tk_d1->x[i] / tk_d1->x[i]) * M_PI * sqrt(Pk_primordial(tk_d1->x[i] * cosmo.h / sim.boxsize, ic) / tk_d1->x[i]) / tk_d1->x[i];

		if (sim.gr_flag == 0)
		{
			for (i = 0; i < tk_t1->size; i++) // construct gauge correction for N-body gauge (3 Hconf theta_tot / k^2)
				temp2[i] = -3. * Hconf(a, fourpiG, cosmo)  * M_PI * tk_t1->y[i] * sqrt(Pk_primordial(tk_d1->x[i] * cosmo.h / sim.boxsize, ic) / tk_d1->x[i]) / tk_d1->x[i] / tk_d1->x[i] / tk_d1->x[i];

			nbspline = gsl_spline_alloc(gsl_interp_cspline, tk_t1->size);
			gsl_spline_init(nbspline, tk_t1->x, temp2, tk_t1->size);
		}
		
		pkspline = gsl_spline_alloc(gsl_interp_cspline, tk_d1->size);
		gsl_spline_init(pkspline, tk_d1->x, temp1, tk_d1->size);
		gsl_spline_free(tk_d1);
		gsl_spline_free(tk_t1);

#ifdef HAVE_CLASS
		if (ic.tkfile[0] == '\0')
		{
			if (sim.gr_flag == 0)
			{
				loadTransferFunctions(class_background, class_perturbs, tk_d1, tk_t1, NULL, sim.boxsize, sim.z_in, cosmo.h);

				for (i = 0; i < tk_d1->size; i++)
					temp1[i] = -tk_d1->y[i];

				gsl_spline_free(tk_d1);
				gsl_spline_free(tk_t1);

				loadTransferFunctions(class_background, class_perturbs, tk_d1, tk_t1, NULL, sim.boxsize, (sim.z_in + 0.01) / 0.99, cosmo.h);

				for (i = 0; i < tk_d1->size; i++)
					temp1[i] += tk_d1->y[i];

				gsl_spline_free(tk_d1);
				gsl_spline_free(tk_t1);

				loadTransferFunctions(class_background, class_perturbs, tk_d1, tk_t1, "tot", sim.boxsize, (sim.z_in + 0.01) / 0.99, cosmo.h);

				for (i = 0; i < tk_d1->size; i++) // construct gauge correction for N-body gauge velocities
					temp1[i] = -99.5 * Hconf(0.995 * a, fourpiG, cosmo) * (3. * temp1[i] * M_PI * sqrt(Pk_primordial(tk_d1->x[i] * cosmo.h / sim.boxsize, ic) / tk_d1->x[i]) / tk_d1->x[i] + (temp2[i] + 3. * Hconf(0.99 * a, fourpiG, cosmo)  * M_PI * tk_t1->y[i] * sqrt(Pk_primordial(tk_d1->x[i] * cosmo.h / sim.boxsize, ic) / tk_d1->x[i]) / tk_d1->x[i] / tk_d1->x[i] / tk_d1->x[i]));

				vnbspline = gsl_spline_alloc(gsl_interp_cspline, tk_t1->size);
				gsl_spline_init(vnbspline, tk_t1->x, temp1, tk_t1->size);

				gsl_spline_free(tk_d1);
				gsl_spline_free(tk_t1);
			}

			loadTransferFunctions(class_background, class_perturbs, tk_d1, tk_t1, "cdm", sim.boxsize, sim.z_in, cosmo.h);
		}
		else
#endif		
		loadTransferFunctions(ic.tkfile, tk_d1, tk_t1, "cdm", sim.boxsize, cosmo.h);	// get transfer functions for CDM
		
		if (tk_d1 == NULL || tk_t1 == NULL)
		{
			COUT << " error: cdm transfer function was empty!" << endl;
			parallel.abortForce();
		}
		
		if (sim.baryon_flag > 0)
		{
#ifdef HAVE_CLASS
			if (ic.tkfile[0] == '\0')
				loadTransferFunctions(class_background, class_perturbs, tk_d2, tk_t2, "b", sim.boxsize, sim.z_in, cosmo.h);
			else
#endif
			loadTransferFunctions(ic.tkfile, tk_d2, tk_t2, "b", sim.boxsize, cosmo.h);	// get transfer functions for baryons
		
			if (tk_d2 == NULL || tk_t2 == NULL)
			{
				COUT << " error: baryon transfer function was empty!" << endl;
				parallel.abortForce();
			}
			if (tk_d2->size != tk_d1->size)
			{
				COUT << " error: baryon transfer function line number mismatch!" << endl;
				parallel.abortForce();
			}
		}
		
		if (sim.baryon_flag == 2)	// baryon treatment = blend; compute displacement & velocity from weighted average
		{
			if (sim.gr_flag > 0)
			{
				for (i = 0; i < tk_d1->size; i++)
					temp1[i] = -3. * pkspline->y[i] / pkspline->x[i] / pkspline->x[i] - ((cosmo.Omega_cdm * tk_d1->y[i] + cosmo.Omega_b * tk_d2->y[i]) / (cosmo.Omega_cdm + cosmo.Omega_b)) * M_PI * sqrt(Pk_primordial(tk_d1->x[i] * cosmo.h / sim.boxsize, ic) / tk_d1->x[i]) / tk_d1->x[i];
			}
			else
			{
				for (i = 0; i < tk_d1->size; i++)
					temp1[i] = nbspline->y[i] - ((cosmo.Omega_cdm * tk_d1->y[i] + cosmo.Omega_b * tk_d2->y[i]) / (cosmo.Omega_cdm + cosmo.Omega_b)) * M_PI * sqrt(Pk_primordial(tk_d1->x[i] * cosmo.h / sim.boxsize, ic) / tk_d1->x[i]) / tk_d1->x[i];
			}
			if (sim.gr_flag > 0 || vnbspline == NULL)
			{
				for (i = 0; i < tk_d1->size; i++)
					temp2[i] = -a * ((cosmo.Omega_cdm * tk_t1->y[i] + cosmo.Omega_b * tk_t2->y[i]) / (cosmo.Omega_cdm + cosmo.Omega_b)) * M_PI * sqrt(Pk_primordial(tk_d1->x[i] * cosmo.h / sim.boxsize, ic) / tk_d1->x[i]) / tk_d1->x[i];
			}
			else
			{
				for (i = 0; i < tk_d1->size; i++)
					temp2[i] = a * vnbspline->y[i] - a * ((cosmo.Omega_cdm * tk_t1->y[i] + cosmo.Omega_b * tk_t2->y[i]) / (cosmo.Omega_cdm + cosmo.Omega_b)) * M_PI * sqrt(Pk_primordial(tk_d1->x[i] * cosmo.h / sim.boxsize, ic) / tk_d1->x[i]) / tk_d1->x[i];
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
		
		if (sim.baryon_flag == 3)	// baryon treatment = hybrid; compute displacement & velocity from weighted average (sub-species)
		{
			if (8. * cosmo.Omega_b / (cosmo.Omega_cdm + cosmo.Omega_b) > 1.)
			{
				if (sim.gr_flag > 0)
				{
					for (i = 0; i < tk_d1->size; i++)
						temp1[i] = -3. * pkspline->y[i] / pkspline->x[i] / pkspline->x[i] - ((8. * cosmo.Omega_cdm * tk_d1->y[i] + (7. * cosmo.Omega_b - cosmo.Omega_cdm) * tk_d2->y[i]) / (cosmo.Omega_cdm + cosmo.Omega_b) / 7.) * M_PI * sqrt(Pk_primordial(tk_d1->x[i] * cosmo.h / sim.boxsize, ic) / tk_d1->x[i]) / tk_d1->x[i];
				}
				else
				{
					for (i = 0; i < tk_d1->size; i++)
						temp1[i] = nbspline->y[i] - ((8. * cosmo.Omega_cdm * tk_d1->y[i] + (7. * cosmo.Omega_b - cosmo.Omega_cdm) * tk_d2->y[i]) / (cosmo.Omega_cdm + cosmo.Omega_b) / 7.) * M_PI * sqrt(Pk_primordial(tk_d1->x[i] * cosmo.h / sim.boxsize, ic) / tk_d1->x[i]) / tk_d1->x[i];
				}
				if (sim.gr_flag > 0 || vnbspline == NULL)
				{
					for (i = 0; i < tk_d1->size; i++)
						temp2[i] = -a * ((8. * cosmo.Omega_cdm * tk_t1->y[i] + (7. * cosmo.Omega_b - cosmo.Omega_cdm) * tk_t2->y[i]) / (cosmo.Omega_cdm + cosmo.Omega_b) / 7.) * M_PI * sqrt(Pk_primordial(tk_d1->x[i] * cosmo.h / sim.boxsize, ic) / tk_d1->x[i]) / tk_d1->x[i];
				}
				else
				{
					for (i = 0; i < tk_d1->size; i++)
						temp2[i] = a * vnbspline->y[i] - a * ((8. * cosmo.Omega_cdm * tk_t1->y[i] + (7. * cosmo.Omega_b - cosmo.Omega_cdm) * tk_t2->y[i]) / (cosmo.Omega_cdm + cosmo.Omega_b) / 7.) * M_PI * sqrt(Pk_primordial(tk_d1->x[i] * cosmo.h / sim.boxsize, ic) / tk_d1->x[i]) / tk_d1->x[i];
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
				if (sim.gr_flag > 0)
				{
					for (i = 0; i < tk_d1->size; i++)
						temp1[i] = -3. * pkspline->y[i] / pkspline->x[i] / pkspline->x[i] - (((cosmo.Omega_cdm - 7. * cosmo.Omega_b) * tk_d1->y[i] + 8. * cosmo.Omega_b * tk_d2->y[i]) / (cosmo.Omega_cdm + cosmo.Omega_b)) * M_PI * sqrt(Pk_primordial(tk_d1->x[i] * cosmo.h / sim.boxsize, ic) / tk_d1->x[i]) / tk_d1->x[i];
				}
				else
				{
					for (i = 0; i < tk_d1->size; i++)
						temp1[i] = nbspline->y[i] - (((cosmo.Omega_cdm - 7. * cosmo.Omega_b) * tk_d1->y[i] + 8. * cosmo.Omega_b * tk_d2->y[i]) / (cosmo.Omega_cdm + cosmo.Omega_b)) * M_PI * sqrt(Pk_primordial(tk_d1->x[i] * cosmo.h / sim.boxsize, ic) / tk_d1->x[i]) / tk_d1->x[i];
				}
				if (sim.gr_flag > 0 || vnbspline == NULL)
				{
					for (i = 0; i < tk_d1->size; i++)
						temp2[i] = -a * (((cosmo.Omega_cdm - 7. * cosmo.Omega_b) * tk_t1->y[i] + 8. * cosmo.Omega_b * tk_t2->y[i]) / (cosmo.Omega_cdm + cosmo.Omega_b)) * M_PI * sqrt(Pk_primordial(tk_d1->x[i] * cosmo.h / sim.boxsize, ic) / tk_d1->x[i]) / tk_d1->x[i];
				}
				else
				{
					for (i = 0; i < tk_d1->size; i++)
						temp2[i] = a * vnbspline->y[i] - a * (((cosmo.Omega_cdm - 7. * cosmo.Omega_b) * tk_t1->y[i] + 8. * cosmo.Omega_b * tk_t2->y[i]) / (cosmo.Omega_cdm + cosmo.Omega_b)) * M_PI * sqrt(Pk_primordial(tk_d1->x[i] * cosmo.h / sim.boxsize, ic) / tk_d1->x[i]) / tk_d1->x[i];
				}

				gsl_spline_free(tk_d2);
				gsl_spline_free(tk_t2);
				tk_d2 = gsl_spline_alloc(gsl_interp_cspline, tk_d1->size);
				tk_t2 = gsl_spline_alloc(gsl_interp_cspline, tk_d1->size);
				gsl_spline_init(tk_d2, tk_d2->x, temp1, tk_d1->size);
				gsl_spline_init(tk_t2, tk_d2->x, temp2, tk_d1->size);
			}
		}
		
		if (sim.baryon_flag == 1 || (sim.baryon_flag == 3 && 8. * cosmo.Omega_b / (cosmo.Omega_cdm + cosmo.Omega_b) > 1.)) // compute baryonic displacement & velocity
		{
			if (sim.gr_flag > 0)
			{
				for (i = 0; i < tk_d1->size; i++)
					temp1[i] = -3. * pkspline->y[i] / pkspline->x[i] / pkspline->x[i] - tk_d2->y[i] * M_PI * sqrt(Pk_primordial(tk_d2->x[i] * cosmo.h / sim.boxsize, ic) / tk_d2->x[i]) / tk_d2->x[i];
			}
			else
			{
				for (i = 0; i < tk_d1->size; i++)
					temp1[i] = nbspline->y[i] - tk_d2->y[i] * M_PI * sqrt(Pk_primordial(tk_d2->x[i] * cosmo.h / sim.boxsize, ic) / tk_d2->x[i]) / tk_d2->x[i];
			}
			if (sim.gr_flag > 0 || vnbspline == NULL)
			{
				for (i = 0; i < tk_d1->size; i++)
					temp2[i] = -a * tk_t2->y[i] * M_PI * sqrt(Pk_primordial(tk_d2->x[i] * cosmo.h / sim.boxsize, ic) / tk_d2->x[i]) / tk_d2->x[i];
			}
			else
			{
				for (i = 0; i < tk_d1->size; i++)
					temp2[i] = a * vnbspline->y[i] - a * tk_t2->y[i] * M_PI * sqrt(Pk_primordial(tk_d2->x[i] * cosmo.h / sim.boxsize, ic) / tk_d2->x[i]) / tk_d2->x[i];
			}

			gsl_spline_free(tk_d2);
			gsl_spline_free(tk_t2);
			tk_d2 = gsl_spline_alloc(gsl_interp_cspline, tk_d1->size);
			tk_t2 = gsl_spline_alloc(gsl_interp_cspline, tk_d1->size);
			gsl_spline_init(tk_d2, tk_d1->x, temp1, tk_d1->size);
			gsl_spline_init(tk_t2, tk_d1->x, temp2, tk_d1->size);
		}
		
		if (sim.baryon_flag < 2 || (sim.baryon_flag == 3 && 8. * cosmo.Omega_b / (cosmo.Omega_cdm + cosmo.Omega_b) <= 1.))	// compute CDM displacement & velocity
		{
			if (sim.gr_flag > 0)
			{
				for (i = 0; i < tk_d1->size; i++)
					temp1[i] = -3. * pkspline->y[i] / pkspline->x[i] / pkspline->x[i] - tk_d1->y[i] * M_PI * sqrt(Pk_primordial(tk_d1->x[i] * cosmo.h / sim.boxsize, ic) / tk_d1->x[i]) / tk_d1->x[i];
			}
			else
			{
				for (i = 0; i < tk_d1->size; i++)
					temp1[i] = nbspline->y[i] - tk_d1->y[i] * M_PI * sqrt(Pk_primordial(tk_d1->x[i] * cosmo.h / sim.boxsize, ic) / tk_d1->x[i]) / tk_d1->x[i];
			}
			if (sim.gr_flag > 0 || vnbspline == NULL)
			{	
				for (i = 0; i < tk_d1->size; i++)
					temp2[i] = -a * tk_t1->y[i] * M_PI * sqrt(Pk_primordial(tk_d1->x[i] * cosmo.h / sim.boxsize, ic) / tk_d1->x[i]) / tk_d1->x[i];
			}
			else
			{
				for (i = 0; i < tk_d1->size; i++)
					temp2[i] = a * vnbspline->y[i] - a * tk_t1->y[i] * M_PI * sqrt(Pk_primordial(tk_d1->x[i] * cosmo.h / sim.boxsize, ic) / tk_d1->x[i]) / tk_d1->x[i];
			}

			gsl_spline_free(tk_d1);
			gsl_spline_free(tk_t1);
			tk_d1 = gsl_spline_alloc(gsl_interp_cspline, pkspline->size);
			tk_t1 = gsl_spline_alloc(gsl_interp_cspline, pkspline->size);
			gsl_spline_init(tk_d1, pkspline->x, temp1, pkspline->size);
			gsl_spline_init(tk_t1, pkspline->x, temp2, pkspline->size);
		}
		
		if ((sim.baryon_flag == 1 && !(ic.flags & ICFLAG_CORRECT_DISPLACEMENT)) || sim.baryon_flag == 3)
		{
			generateDisplacementField(*scalarFT, 0., tk_d2, (unsigned int) ic.seed, ic.flags & ICFLAG_KSPHERE);
			gsl_spline_free(tk_d2);
			plan_phi->execute(FFT_BACKWARD);
			phi->updateHalo();	// phi now contains the baryonic displacement
			plan_source->execute(FFT_FORWARD);
		}
		
		generateDisplacementField(*scalarFT, 0., tk_d1, (unsigned int) ic.seed, ic.flags & ICFLAG_KSPHERE);
		gsl_spline_free(tk_d1);
	}
		
	plan_chi->execute(FFT_BACKWARD);
	chi->updateHalo();	// chi now contains the CDM displacement
	
	strcpy(pcls_cdm_info.type_name, "part_simple");
	if (sim.baryon_flag == 1)
		pcls_cdm_info.mass = cosmo.Omega_cdm / (Real) (sim.numpcl[0]*(long)ic.numtile[0]*(long)ic.numtile[0]*(long)ic.numtile[0]);
	else
		pcls_cdm_info.mass = (cosmo.Omega_cdm + cosmo.Omega_b) / (Real) (sim.numpcl[0]*(long)ic.numtile[0]*(long)ic.numtile[0]*(long)ic.numtile[0]);
	pcls_cdm_info.relativistic = false;
	
	pcls_cdm->initialize(pcls_cdm_info, pcls_cdm_dataType, &(phi->lattice()), boxSize);
	
	initializeParticlePositions(sim.numpcl[0], pcldata, ic.numtile[0], *pcls_cdm);
	i = MAX;
	if (sim.baryon_flag == 3)	// baryon treatment = hybrid; displace particles using both displacement fields
		pcls_cdm->moveParticles(displace_pcls_ic_basic, 1., ic_fields, 2, NULL, &max_displacement, &i, 1);
	else
		pcls_cdm->moveParticles(displace_pcls_ic_basic, 1., &chi, 1, NULL, &max_displacement, &i, 1);	// displace CDM particles
	
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
			generateCICKernel(*phi, sim.numpcl[1], pcldata, ic.numtile[1]);
			plan_phi->execute(FFT_FORWARD);
			generateDisplacementField(*scalarFT, 0., tk_d2, (unsigned int) ic.seed, ic.flags & ICFLAG_KSPHERE);
			gsl_spline_free(tk_d2);
			plan_phi->execute(FFT_BACKWARD);
			phi->updateHalo();
		}
		
		strcpy(pcls_b_info.type_name, "part_simple");
		pcls_b_info.mass = cosmo.Omega_b / (Real) (sim.numpcl[1]*(long)ic.numtile[1]*(long)ic.numtile[1]*(long)ic.numtile[1]);
		pcls_b_info.relativistic = false;
	
		pcls_b->initialize(pcls_b_info, pcls_b_dataType, &(phi->lattice()), boxSize);
	
		initializeParticlePositions(sim.numpcl[1], pcldata, ic.numtile[1], *pcls_b);
		i = MAX;
		pcls_b->moveParticles(displace_pcls_ic_basic, 1., &phi, 1, NULL, &max_displacement, &i, 1);	// displace baryon particles
	
		sim.numpcl[1] *= (long) ic.numtile[1] * (long) ic.numtile[1] * (long) ic.numtile[1];
	
		COUT << " " << sim.numpcl[1] << " baryon particles initialized: maximum displacement = " << max_displacement * sim.numpts << " lattice units." << endl;
	
		free(pcldata);
	}
	
	if (ic.pkfile[0] == '\0')	// set velocities using transfer functions
	{
		if (ic.flags & ICFLAG_CORRECT_DISPLACEMENT)
			generateCICKernel(*source);
		
		plan_source->execute(FFT_FORWARD);
		
		if (sim.baryon_flag == 1 || sim.baryon_flag == 3)
		{
			generateDisplacementField(*scalarFT, 0., tk_t2, (unsigned int) ic.seed, ic.flags & ICFLAG_KSPHERE, 0);
			plan_phi->execute(FFT_BACKWARD);
			phi->updateHalo();	// phi now contains the baryonic velocity potential
			gsl_spline_free(tk_t2);
			plan_source->execute(FFT_FORWARD);
		}
		
		generateDisplacementField(*scalarFT, 0., tk_t1, (unsigned int) ic.seed, ic.flags & ICFLAG_KSPHERE, 0);
		plan_chi->execute(FFT_BACKWARD);
		chi->updateHalo();	// chi now contains the CDM velocity potential
		gsl_spline_free(tk_t1);			
		
		if (sim.baryon_flag == 3)	// baryon treatment = hybrid; set velocities using both velocity potentials
			maxvel[0] = pcls_cdm->updateVel(initialize_q_ic_basic, 1., ic_fields, 2) / a;
		else
			maxvel[0] = pcls_cdm->updateVel(initialize_q_ic_basic, 1., &chi, 1) / a;	// set CDM velocities
		
		if (sim.baryon_flag == 1)
			maxvel[1] = pcls_b->updateVel(initialize_q_ic_basic, 1., &phi, 1) / a;	// set baryon velocities
	}
	
	if (sim.baryon_flag > 1) sim.baryon_flag = 0;
	
	for (p = 0; p < cosmo.num_ncdm; p++)	// initialization of non-CDM species
	{
		if (ic.numtile[1+sim.baryon_flag+p] < 1) continue;

		loadHomogeneousTemplate(ic.pclfile[1+sim.baryon_flag+p], sim.numpcl[1+sim.baryon_flag+p], pcldata);
	
		if (pcldata == NULL)
		{
			COUT << " error: particle data was empty!" << endl;
			parallel.abortForce();
		}
		
		if (ic.pkfile[0] == '\0')
		{
			sprintf(ncdm_name, "ncdm[%d]", p);
#ifdef HAVE_CLASS
			if (ic.tkfile[0] == '\0')
				loadTransferFunctions(class_background, class_perturbs, tk_d1, tk_t1, ncdm_name, sim.boxsize, sim.z_in, cosmo.h);
			else
#endif
			loadTransferFunctions(ic.tkfile, tk_d1, tk_t1, ncdm_name, sim.boxsize, cosmo.h);
		
			if (tk_d1 == NULL || tk_t1 == NULL)
			{
				COUT << " error: ncdm transfer function was empty! (species " << p << ")" << endl;
				parallel.abortForce();
			}
			
			if (sim.gr_flag > 0)
			{
				for (i = 0; i < tk_d1->size; i++)
					temp1[i] = -3. * pkspline->y[i] / pkspline->x[i] / pkspline->x[i] - tk_d1->y[i] * M_PI * sqrt(Pk_primordial(tk_d1->x[i] * cosmo.h / sim.boxsize, ic) / tk_d1->x[i]) / tk_d1->x[i];
			}
			else
			{
				for (i = 0; i < tk_d1->size; i++)
					temp1[i] = nbspline->y[i] - tk_d1->y[i] * M_PI * sqrt(Pk_primordial(tk_d1->x[i] * cosmo.h / sim.boxsize, ic) / tk_d1->x[i]) / tk_d1->x[i];
			}
			if (sim.gr_flag > 0 || vnbspline == NULL)
			{
				for (i = 0; i < tk_d1->size; i++)
					temp2[i] = -a * tk_t1->y[i] * M_PI * sqrt(Pk_primordial(tk_d1->x[i] * cosmo.h / sim.boxsize, ic) / tk_d1->x[i]) / tk_d1->x[i];
			}
			else
			{
				for (i = 0; i < tk_d1->size; i++)
					temp2[i] = a * vnbspline->y[i] - a * tk_t1->y[i] * M_PI * sqrt(Pk_primordial(tk_d1->x[i] * cosmo.h / sim.boxsize, ic) / tk_d1->x[i]) / tk_d1->x[i];
			}

			gsl_spline_free(tk_d1);
			gsl_spline_free(tk_t1);
			tk_d1 = gsl_spline_alloc(gsl_interp_cspline, pkspline->size);
			tk_t1 = gsl_spline_alloc(gsl_interp_cspline, pkspline->size);
			gsl_spline_init(tk_d1, pkspline->x, temp1, pkspline->size);
			gsl_spline_init(tk_t1, pkspline->x, temp2, pkspline->size);
			
			plan_source->execute(FFT_FORWARD);
			generateDisplacementField(*scalarFT, 0., tk_d1, (unsigned int) ic.seed, ic.flags & ICFLAG_KSPHERE);
			plan_chi->execute(FFT_BACKWARD);	// chi now contains the displacement for the non-CDM species
			chi->updateHalo();
			gsl_spline_free(tk_d1);
			
			plan_source->execute(FFT_FORWARD);
			generateDisplacementField(*scalarFT, 0., tk_t1, (unsigned int) ic.seed, ic.flags & ICFLAG_KSPHERE, 0);
			plan_phi->execute(FFT_BACKWARD);	// phi now contains the velocity potential for the non-CDM species
			phi->updateHalo();
			gsl_spline_free(tk_t1);
		}
		
		strcpy(pcls_ncdm_info[p].type_name, "part_simple");
		pcls_ncdm_info[p].mass = cosmo.Omega_ncdm[p] / (Real) (sim.numpcl[1+sim.baryon_flag+p]*(long)ic.numtile[1+sim.baryon_flag+p]*(long)ic.numtile[1+sim.baryon_flag+p]*(long)ic.numtile[1+sim.baryon_flag+p]);
		pcls_ncdm_info[p].relativistic = true;
		
		pcls_ncdm[p].initialize(pcls_ncdm_info[p], pcls_ncdm_dataType, &(phi->lattice()), boxSize);
		
		initializeParticlePositions(sim.numpcl[1+sim.baryon_flag+p], pcldata, ic.numtile[1+sim.baryon_flag+p], pcls_ncdm[p]);
		i = MAX;
		pcls_ncdm[p].moveParticles(displace_pcls_ic_basic, 1., &chi, 1, NULL, &max_displacement, &i, 1);	// displace non-CDM particles
		
		sim.numpcl[1+sim.baryon_flag+p] *= (long) ic.numtile[1+sim.baryon_flag+p] * (long) ic.numtile[1+sim.baryon_flag+p] * (long) ic.numtile[1+sim.baryon_flag+p];
	
		COUT << " " << sim.numpcl[1+sim.baryon_flag+p] << " ncdm particles initialized for species " << p+1 << ": maximum displacement = " << max_displacement * sim.numpts << " lattice units." << endl;
		
		free(pcldata);
		
		if (ic.pkfile[0] == '\0')	// set non-CDM velocities using transfer functions
			pcls_ncdm[p].updateVel(initialize_q_ic_basic, 1., &phi, 1);
	}

	free(temp1);
	free(temp2);
	
	if (ic.pkfile[0] == '\0')
	{
		plan_source->execute(FFT_FORWARD);
		generateDisplacementField(*scalarFT, 0., pkspline, (unsigned int) ic.seed, ic.flags & ICFLAG_KSPHERE, 0);
#ifdef HAVE_CLASS
		if (ic.tkfile[0] == '\0')
			freeCLASSstructures(class_background, class_thermo, class_perturbs);
#endif
	}
	else
	{
		projection_init(source);
		scalarProjectionCIC_project(pcls_cdm, source);
		if (sim.baryon_flag)
			scalarProjectionCIC_project(pcls_b, source);
		for (p = 0; p < cosmo.num_ncdm; p++)
		{
			if (ic.numtile[1+sim.baryon_flag+p] < 1) continue;
			scalarProjectionCIC_project(pcls_ncdm+p, source);
		}
		scalarProjectionCIC_comm(source);
	
		plan_source->execute(FFT_FORWARD);
	
		kFT.first();
		if (kFT.coord(0) == 0 && kFT.coord(1) == 0 && kFT.coord(2) == 0)
			(*scalarFT)(kFT) = Cplx(0.,0.);
				
		solveModifiedPoissonFT(*scalarFT, *scalarFT, fourpiG / a, 3. * sim.gr_flag * (Hconf(a, fourpiG, cosmo) * Hconf(a, fourpiG, cosmo) + fourpiG * cosmo.Omega_m / a));
	}
	
	plan_phi->execute(FFT_BACKWARD);
	phi->updateHalo();	// phi now finally contains phi
	
	if (ic.pkfile[0] != '\0')	// if power spectrum is used instead of transfer functions, set velocities using linear approximation
	{	
		rescale = a / Hconf(a, fourpiG, cosmo) / (1.5 * Omega_m(a, cosmo) + 2. * Omega_rad(a, cosmo));
		maxvel[0] = pcls_cdm->updateVel(initialize_q_ic_basic, rescale, &phi, 1) / a;
		if (sim.baryon_flag)
			maxvel[1] = pcls_b->updateVel(initialize_q_ic_basic, rescale, &phi, 1) / a;
	}
			
	for (p = 0; p < cosmo.num_ncdm; p++)
	{
		if (ic.numtile[1+sim.baryon_flag+p] < 1)
		{
			maxvel[1+sim.baryon_flag+p] = 0;
			continue;
		}

		if (ic.pkfile[0] != '\0') // if power spectrum is used instead of transfer functions, set bulk velocities using linear approximation
		{		
			rescale = a / Hconf(a, fourpiG, cosmo) / (1.5 * Omega_m(a, cosmo) + Omega_rad(a, cosmo));
			pcls_ncdm[p].updateVel(initialize_q_ic_basic, rescale, &phi, 1);
		}
		
		if (cosmo.m_ncdm[p] > 0.) // add velocity dispersion for non-CDM species
		{
			rescale = pow(cosmo.Omega_g * cosmo.h * cosmo.h / C_PLANCK_LAW, 0.25) * cosmo.T_ncdm[p] * C_BOLTZMANN_CST / cosmo.m_ncdm[p];
			mean_q = applyMomentumDistribution(pcls_ncdm+p, (unsigned int) (ic.seed + p), rescale);
			parallel.sum(mean_q);
			COUT << " species " << p+1 << " Fermi-Dirac distribution had mean q/m = " << mean_q / sim.numpcl[1+sim.baryon_flag+p] << endl;
		}
		maxvel[1+sim.baryon_flag+p] = pcls_ncdm[p].updateVel(update_q, 0., &phi, 1, &a);
	}
	
	projection_init(Bi);
	projection_T0i_project(pcls_cdm, Bi, phi);
	if (sim.baryon_flag)
		projection_T0i_project(pcls_b, Bi, phi);
	projection_T0i_comm(Bi);
	plan_Bi->execute(FFT_FORWARD);
	projectFTvector(*BiFT, *BiFT, fourpiG / (double) sim.numpts / (double) sim.numpts);	
	plan_Bi->execute(FFT_BACKWARD);	
	Bi->updateHalo();	// B initialized
	
	projection_init(Sij);
	projection_Tij_project(pcls_cdm, Sij, a, phi);
	if (sim.baryon_flag)
		projection_Tij_project(pcls_b, Sij, a, phi);
	projection_Tij_comm(Sij);
	
	prepareFTsource<Real>(*phi, *Sij, *Sij, 2. * fourpiG / a / (double) sim.numpts / (double) sim.numpts);	
	plan_Sij->execute(FFT_FORWARD);	
	projectFTscalar(*SijFT, *scalarFT);
	plan_chi->execute(FFT_BACKWARD);		
	chi->updateHalo();	// chi now finally contains chi

	gsl_spline_free(pkspline);
	if (sim.gr_flag == 0)
		gsl_spline_free(nbspline);
	if (vnbspline != NULL)
		gsl_spline_free(vnbspline);
}

#endif

#endif

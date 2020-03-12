#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <vector>
#include <map>
#include <array>
#include "chealpix.h"
#include "healpix_base.h"
#include "metadata.hpp"
#include "parser.hpp"
#include "background.hpp"


using namespace std;

struct background_data
{
	int cycle;
	double tau;
	double a;
};

struct metric_data
{
	healpix_header hdr;
	float * pixel;
};

struct metric_container
{
	std::map<int,metric_data> healpix_data;
	background_data cinfo;
	char name[4];
	char * dir;
	char * basename;
	void init(background_data & c, char * d, char * b, const char * n)
	{
		cinfo = c;
		strcpy(name, n);
		dir = d;
		basename = b;
	}
	void truncate(double min_dist)
	{
		std::map<int,metric_data>::iterator it;

		for (it = healpix_data.begin(); it != healpix_data.end(); it++)
		{
			if (it->second.hdr.distance >= min_dist) break;
			else free(it->second.pixel);
		}

		if (it != healpix_data.begin()) healpix_data.erase(healpix_data.begin(), it);
	}
	void clear()
	{
		for (std::map<int,metric_data>::iterator it = healpix_data.begin(); it != healpix_data.end(); it++)
			free(it->second.pixel);
		healpix_data.clear();
	}
};


bool kappa(float * pixel, const int64_t Nside, int64_t ipix, float & result);

int loadHealpixData(metric_container * field, double min_dist, double max_dist, char * lightconeparam = NULL);

float lin_int(float a, float b, float f)
{
	return a + (f*(b-a));
}


int main(int argc, char **argv)
{

	char * settingsfile = NULL;
	char * redshiftparam = NULL;
	char outputfile[1024];
	char * distanceparam = NULL;
	char * lightconeparam = NULL;

	parameter * params = NULL;
	metadata sim;
	cosmology cosmo;
	icsettings ic;
	int numparam = 0;
	int usedparams = 0;

	int numoutputs=0;

	char filename0[1024];
	char filename1[1024];

	char coordsys = 'G';
	Healpix_Base2 helper;
	pointing ptg;
	fix_arr<int64,4> nnpix;
	fix_arr<double,4> nnwgt;

	if (argc < 2)
	{
		cout << COLORTEXT_WHITE << " LCARS tools: lcmap" << COLORTEXT_RESET << endl;
		cout << " extracts local and integrated potential terms from metric light-cone output" << endl << endl;
		
		cout << " List of command-line options:" << endl;
		cout << " -s <filename>       : gevolution settings file of the simulation (mandatory)" << endl;
		cout << " -d <distance>[,...] : conformal distance to source field [Mpc/h] (alternative 1)" << endl;
		cout << " -z <redshift>[,...] : redshift of source field (alternative 2)" << endl;
		cout << " -l <ID1>[,<ID2>,...]: IDs of light cones to be included (optional, must refer" << endl;
		cout << "                       to the same observation event)" << endl;
		cout << " The output will be written to HEALPix FITS files that follow the naming conventions" << endl;
		cout << " specified in the settings file." << endl;

		return 0;
	}


	for (int i = 1 ; i < argc ; i++ )
	{
		if ( argv[i][0] != '-' )
			continue;
		switch(argv[i][1]) {
			case 'z':
				redshiftparam = argv[++i]; //redshifts for calculation 
				break;
			case 's':
				settingsfile = argv[++i]; //settings file name
				break;
			case 'd':
				distanceparam= argv[++i]; //distances for calculation
				break;
			case 'l':
				lightconeparam = argv[++i];
		}
	}

	if (settingsfile == NULL)
	{
		cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": no settings file specified!" << endl;
		return -1;
	}

	cout << COLORTEXT_WHITE << " LCARS tools: lcmap" << endl << endl << " opening settings file of simulation: " << settingsfile << endl << " parser output:" << COLORTEXT_RESET << endl << endl;

	numparam = loadParameterFile(settingsfile, params);
	usedparams = parseMetadata(params, numparam, sim, cosmo, ic);

	cout << endl << " file contains " << numparam << " parameters, " << usedparams << " of which could be parsed." << endl << endl;

	cout << " number of lightcones: " << sim.num_lightcone << endl << endl;

	for (int i = 0; i < sim.num_lightcone; i++)
	{
		cout << " lightcone " << i << " parameters:" << endl << "  vertex = " << sim.lightcone[i].vertex[0] << ", " << sim.lightcone[i].vertex[1] << ", " << sim.lightcone[i].vertex[2] << endl;
		cout << "  redshift of observation = " << sim.lightcone[i].z << endl;
		cout << "  direction = " << sim.lightcone[i].direction[0] << ", " << sim.lightcone[i].direction[1] << ", " << sim.lightcone[i].direction[2] << endl;
		cout << "  opening half-angle = " << ((sim.lightcone[i].opening > -1.) ? acos(sim.lightcone[i].opening) * 180. / M_PI : 180.) << " degrees" << endl;
		cout << "  distance interval = " << sim.lightcone[i].distance[0] << ", " << sim.lightcone[i].distance[1] << endl << endl;
	}
	
	double tauobs = particleHorizon(1, 1.5 * sim.boxsize * sim.boxsize / C_SPEED_OF_LIGHT / C_SPEED_OF_LIGHT, cosmo);
	
	FILE * background_file;
	char * buffer;

	int numlines = -2;

	sprintf(filename0,"%s%s_background.dat",sim.output_path,sim.basename_generic);
	background_file = fopen(filename0, "r");

	if (background_file==NULL)
	{
		cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": background file " << filename0 << " cannot be opened!" << endl;
		return -1;
	}

	while(!feof(background_file))
	{
		if(fscanf(background_file, "%*[^\n]\n") != 0)
		{
			cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": read error! " << endl;
			return 0;
		}
		numlines++;
	}

	rewind(background_file);

	if(fscanf(background_file, "%*[^\n]\n") != 0)
	{
		cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": read error! " << endl;
		return 0;
	}
	if(fscanf(background_file, "%*[^\n]\n") != 0)
	{
		cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": read error! " << endl;
		return 0;
	}

	background_data back[numlines];

	for (int i=0;i<numlines;i++)
	{
		if(fscanf(background_file, "%i %lf %lf %*f %*f %*f \n", &back[numlines-1-i].cycle, &back[numlines-1-i].tau, &back[numlines-1-i].a) != 3)
		{
			cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": read error! " << endl;
			return 0;
		}
	}

	fclose(background_file);
	
	double maxredshift = (1./back[numlines-1].a) - 1.;
	double maxdistance = (tauobs - back[numlines-1].tau) * sim.boxsize;;

	char * token = NULL;

	if(redshiftparam != NULL && distanceparam != NULL)
	{
		cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": please specify EITHER redshift OR distance to source field, not both!" << endl;
		return 0;
	}
	else if(redshiftparam == NULL && distanceparam == NULL)
	{
		cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": please specify redshift (parameter -z) OR distance (parameter -d) to source field!" << endl;
		return 0;		
	}
	else if(redshiftparam == NULL && distanceparam != NULL)
	{
		sprintf(filename0, "%s", distanceparam);

		token = strtok(filename0, ",");

		while (token != NULL)
		{
			numoutputs++;
			token = strtok(NULL, ",");
		}		
	}
	else if(redshiftparam != NULL && distanceparam == NULL)
	{
		sprintf(filename0, "%s", redshiftparam);

		token = strtok(filename0, ",");

		while (token != NULL)
		{
			numoutputs++;
			token = strtok(NULL, ",");
		}		
	}


	double redshifts[numoutputs];
	double distances[numoutputs];

	int n = 0;

	if(redshiftparam == NULL && distanceparam != NULL)
	{

		sprintf(filename0, "%s", distanceparam);

		token = strtok(filename0, ",");

		while (token != NULL)
		{
			distances[n] = atof(token);
			if (distances[n] > maxdistance)
			{
				cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": distance of " << distances[n] << " Mpc/h is too large, maximum distance is " << maxdistance << " Mpc/h" << endl;
				return -1;
			}
			n++;
			token = strtok(NULL, ",");
		}
		
		qsort((void *) distances, numoutputs, sizeof(double), [](const void * a, const void * b)->int{ return (int) (*((double*) a) > *((double*) b)); });

		cout << " distances (Mpc/h) chosen are: ";

		for (int i = 0; i < numoutputs; i++)
		{
			if(i!=0) cout << ", ";
			cout << distances[i];
			distances[i] /= sim.boxsize;
		}
		cout << "." << endl << endl;
	}
	else if(redshiftparam != NULL && distanceparam == NULL)
	{

		sprintf(filename0, "%s", redshiftparam);

		token = strtok(filename0, ",");

		while (token != NULL)
		{
			redshifts[n] = atof(token);
			if (redshifts[n] > maxredshift)
			{
				cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": redshift of " << redshifts[n] << " is too large, maximum redshift is " << maxredshift << endl;
				return -1;
			}
			n++;
			token = strtok(NULL, ",");
		}
		
		qsort((void *) redshifts, numoutputs, sizeof(double), [](const void * a, const void * b)->int{ return (int) (*((double*) a) > *((double*) b)); });

		cout << " redshifts chosen are: ";

		for (int i = 0; i < numoutputs; i++)
		{
			distances[i] = tauobs - particleHorizon(1./(redshifts[i]+1.), 1.5 * sim.boxsize * sim.boxsize / C_SPEED_OF_LIGHT / C_SPEED_OF_LIGHT, cosmo);
			
			if(i!=0) cout << ", ";
			cout << redshifts[i];
		}
		cout << "." << endl << endl;
	}

	metric_container phi0;
	metric_container phi1;

	phi0.init(back[0], sim.output_path, sim.basename_lightcone, "phi");
	phi1.init(back[1], sim.output_path, sim.basename_lightcone, "phi");

	float * map_phi_final[numoutputs];
	float * map_isw_final = NULL;
	float * map_shapiro_final = NULL;
	float * map_kappa_final = NULL;
	float pot_obs = 0.;

	uint32_t Nside_final=2;
	uint32_t Nside_initial=2;
	int64_t Npix_final=48;
	std::vector<int64_t> Npix_interp;
	std::vector<int64_t>::iterator itNpix;

	int64_t p, q, ipix, jpix, ring, pixoffset=0;

	int step=0;
	int outcnt = 0;
	int cnt = 0;
	int thresh = 1;

	double dist = 0;
	double monopole;
	double v1[3];

	std::map<int,metric_data>::iterator it0;
	std::map<int,metric_data>::iterator it1;
 	
	cout << COLORTEXT_YELLOW << " processing light-cone output..." << COLORTEXT_RESET << endl << endl;


	while(outcnt < numoutputs)
	{
		cout << " interpolating cycle number " << COLORTEXT_CYAN << back[step].cycle << COLORTEXT_RESET << " and " << COLORTEXT_CYAN << back[step+1].cycle << COLORTEXT_RESET << "." << endl << " distance interval: " << (tauobs - back[step].tau) * sim.boxsize << " to " << (tauobs - back[step + 1].tau) * sim.boxsize << endl << endl;

		loadHealpixData(&phi0, tauobs - back[step].tau, tauobs - back[step + 1].tau, lightconeparam);
		loadHealpixData(&phi1, tauobs - back[step].tau, tauobs - back[step + 1].tau, lightconeparam);

		for (it0 = phi0.healpix_data.begin(); it0 != phi0.healpix_data.end() && it0->second.hdr.distance < tauobs - back[step].tau; it0++);
		for (it1 = phi1.healpix_data.begin(); it1 != phi1.healpix_data.end() && it1->second.hdr.distance < tauobs - back[step].tau; it1++);
		
		if (it0 == phi0.healpix_data.end() || it1 == phi0.healpix_data.end())
		{
			cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": missing HEALPix data beyond " << dist * sim.boxsize << " Mpc/h!" << endl;
			return -1;
		}

		if (map_isw_final == NULL)
		{
			Nside_final = it0->second.hdr.Nside;
			Nside_initial = Nside_final;
			Npix_final = it0->second.hdr.Npix;
			for (int m = 0; m < numoutputs; m++)
				map_phi_final[m] = (float *) malloc(Nside_final * Nside_final * 12 * sizeof(float));
			map_isw_final = (float *) malloc(Nside_final * Nside_final * 12 * sizeof(float));
			map_shapiro_final = (float *) malloc(Nside_final * Nside_final * 12 * sizeof(float));
			
			monopole = 0.;

#pragma omp parallel for reduction(+:monopole)
			for (int l = 0; l < it0->second.hdr.Npix; l++)
				monopole += lin_int(it0->second.pixel[l], it1->second.pixel[l], (it0->second.hdr.distance - (tauobs - back[step].tau))/((back[step].tau - back[step+1].tau)));
				
			pot_obs = monopole / it0->second.hdr.Npix;
			
#pragma omp parallel for
			for(int i = 0; i < (Nside_final * Nside_final * 12); i++)
			{
				for (int m = 0; m < numoutputs; m++)
					map_phi_final[m][i] = 0;
				map_isw_final[i] = 0;
				map_shapiro_final[i] = 0;
			}
		}

		for(; it0->second.hdr.distance < tauobs - back[step+1].tau && it1->second.hdr.distance < tauobs - back[step+1].tau && dist < distances[numoutputs-1]; it0++, it1++)
		{
			if (it0 == phi0.healpix_data.end() || it1 == phi0.healpix_data.end())
			{
				cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": missing HEALPix data beyond " << dist * sim.boxsize << " Mpc/h!" << endl;
				return -1;
			}

			if(it0->second.hdr.Nside != Nside_final)
			{
				cout << COLORTEXT_CYAN << " map resolution change" << COLORTEXT_RESET;
				cout << ": from Nside of " << COLORTEXT_WHITE << Nside_final << COLORTEXT_RESET << " to " << COLORTEXT_WHITE <<  it0->second.hdr.Nside << COLORTEXT_RESET << "." << endl << endl;
					 
#pragma omp parallel for	
				for (long l = Npix_final + pixoffset; l < pixoffset + 12l * Nside_final * Nside_final; l++)
				{
					for (int m = outcnt; m < numoutputs; m++)
						map_phi_final[m][l] = -1.6375e30;
					map_isw_final[l] = -1.6375e30;
					map_shapiro_final[l] = -1.6375e30;
				}
				
				Npix_interp.push_back(Npix_final);
				
				Nside_final = it0->second.hdr.Nside;
					 
				for (int m = outcnt; m < numoutputs; m++)
					map_phi_final[m] = (float *) realloc(map_phi_final[m], (16l * Nside_final * Nside_final - 4l * Nside_initial * Nside_initial) * sizeof(float));
				
				map_isw_final = (float *) realloc(map_isw_final, (16l * Nside_final * Nside_final - 4l * Nside_initial * Nside_initial) * sizeof(float));
				map_shapiro_final = (float *) realloc(map_shapiro_final, (16l * Nside_final * Nside_final - 4l * Nside_initial * Nside_initial) * sizeof(float));
				
				pixoffset = 4l * (Nside_final * Nside_final - Nside_initial * Nside_initial);
				
#pragma omp parallel for
				for(long l = pixoffset; l < pixoffset + (Nside_final * Nside_final * 12); l++)
				{
					for (int m = outcnt; m < numoutputs; m++)
						map_phi_final[m][l] = 0;
					map_isw_final[l] = 0;
					map_shapiro_final[l] = 0;
				}
			}

			if (it0->second.hdr.distance <= 0)
			{
				continue;
			}
			
			if (it0->second.hdr.distance != it1->second.hdr.distance || it0->second.hdr.Nside != it1->second.hdr.Nside || it0->second.hdr.Npix != it1->second.hdr.Npix)
			{
				cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": map properties do not match on iteration " << cnt << "! distances: " << it0->second.hdr.distance << ", " << it1->second.hdr.distance << "; Nside: " << it0->second.hdr.Nside << ", " << it1->second.hdr.Nside << "; Npix: " << it0->second.hdr.Npix << ", " << it1->second.hdr.Npix << endl << endl;
				return -1;
			}
			
			Npix_final = it0->second.hdr.Npix;

#pragma omp parallel for
			for (int l = 0; l < Npix_final; l++)
			{
				if(map_isw_final[l+pixoffset] < -1.5e29 || map_shapiro_final[l+pixoffset] < -1.5e29)
				{
					continue;
				}
				else if(it0->second.pixel[l] < -1.5e29)
				{
					for (int m = outcnt; m < numoutputs; m++)
						map_phi_final[m][l+pixoffset] = it0->second.pixel[l];
					map_isw_final[l+pixoffset] = it0->second.pixel[l];
					map_shapiro_final[l+pixoffset] = it0->second.pixel[l];
				}
				else if(it1->second.pixel[l] < -1.5e29)
				{
					for (int m = outcnt; m < numoutputs; m++)
						map_phi_final[m][l+pixoffset] = it1->second.pixel[l];
					map_isw_final[l+pixoffset] = it1->second.pixel[l];
					map_shapiro_final[l+pixoffset] = it1->second.pixel[l];
				}
				else
				{
					for (int m = outcnt; m < numoutputs; m++)
						map_phi_final[m][l+pixoffset] -= 2.*(it0->second.hdr.distance - dist)*((distances[m]-it0->second.hdr.distance)/(distances[m]*it0->second.hdr.distance))*lin_int(it0->second.pixel[l],it1->second.pixel[l], (it0->second.hdr.distance - (tauobs - back[step].tau))/((back[step].tau - back[step+1].tau)));
					map_isw_final[l+pixoffset] += 2.*(it0->second.hdr.distance - dist) * (it0->second.pixel[l] - it1->second.pixel[l]) / (back[step].tau - back[step+1].tau);
					map_shapiro_final[l+pixoffset] -= 2.*sim.boxsize*(it0->second.hdr.distance - dist) * lin_int(it0->second.pixel[l],it1->second.pixel[l], (it0->second.hdr.distance - (tauobs - back[step].tau))/((back[step].tau - back[step+1].tau)));
				}
			}
			
			if (it0->second.hdr.distance >= distances[outcnt])
			{
				cout << " source distance reached of " << distances[outcnt]*sim.boxsize << " Mpc/h!" << endl << endl;

#pragma omp parallel for	
				for (long l = Npix_final; l < 12l * Nside_final * Nside_final; l++)
					map_phi_final[outcnt][l+pixoffset] = -1.6375e30;
					
				cout << " computing convergence..." << endl << endl;
				
				map_kappa_final = (float *) malloc((16l * Nside_final * Nside_final - 4l * Nside_initial * Nside_initial) * sizeof(float));

#pragma omp parallel for				
				for (long l = 0; l < Npix_final; l++)
				{	
					if(!kappa(map_phi_final[outcnt]+pixoffset, Nside_final, l, map_kappa_final[l+pixoffset]))
						map_kappa_final[l+pixoffset] = -1.6375e30;
				}
				
#pragma omp parallel for	
				for (long l = Npix_final; l < 12l * Nside_final * Nside_final; l++)
					map_kappa_final[l+pixoffset] = -1.6375e30;
				
				if (Nside_final > Nside_initial)
				{
					itNpix = Npix_interp.begin();
					
					for (uint32_t Nside_interp = Nside_initial; Nside_interp < Nside_final; Nside_interp <<= 1)
					{
						p = 4l * (Nside_interp * Nside_interp - Nside_initial * Nside_initial);
						
#pragma omp parallel for
						for (long l = 0; l < *itNpix; l++)
						{	
							if(!kappa(map_phi_final[outcnt]+p, Nside_interp, l, map_kappa_final[l+p]))
								map_kappa_final[l+p] = -1.6375e30;
						}
						
#pragma omp parallel for	
						for (long l = *itNpix; l < 12l * Nside_interp * Nside_interp; l++)
							map_kappa_final[l+p] = -1.6375e30;
						
						itNpix++;
					}
					
					for (long l = 0; l < Npix_final; l++)
					{
						pix2vec_ring64(Nside_final, l, v1);
						
						helper.SetNside(Nside_final, RING);
						ptg = helper.pix2ang(l);
						
						for (uint32_t Nside_interp = Nside_initial; Nside_interp < Nside_final; Nside_interp <<= 1)
						{
							helper.SetNside(Nside_interp, RING);
							helper.get_interpol(ptg, nnpix, nnwgt);
							
							p = 4l * (Nside_interp * Nside_interp - Nside_initial * Nside_initial);
							
							for (int nn = 0; nn < nnpix.size(); nn++)
							{
								if (map_kappa_final[nnpix[nn]+p] > -1e30)
									map_kappa_final[l+pixoffset] += map_kappa_final[nnpix[nn]+p] * nnwgt[nn];
								else
								{
									map_kappa_final[l+pixoffset] = -1.6375e30;
									break;
								}
							}
							
							for (int nn = 0; nn < nnpix.size(); nn++)
							{
								if (map_phi_final[outcnt][nnpix[nn]+p] > -1e30)
									map_phi_final[outcnt][l+pixoffset] += map_phi_final[outcnt][nnpix[nn]+p] * nnwgt[nn];
								else
								{
									map_phi_final[outcnt][l+pixoffset] = -1.6375e30;
									break;
								}
							}
						}
					}
				}
	
				if (numoutputs > 1)
					sprintf(outputfile, "%s%s_kappa_%d.fits", sim.output_path, sim.basename_lightcone, outcnt);	
				else
					sprintf(outputfile, "%s%s_kappa.fits", sim.output_path, sim.basename_lightcone);	
				write_healpix_map(map_kappa_final+pixoffset, Nside_final, outputfile, 0, &coordsys);
				
				free(map_kappa_final);
				
				if (numoutputs > 1)
					sprintf(outputfile, "%s%s_lensingphi_%d.fits", sim.output_path, sim.basename_lightcone, outcnt);	
				else
					sprintf(outputfile, "%s%s_lensingphi.fits", sim.output_path, sim.basename_lightcone);	
				write_healpix_map(map_phi_final[outcnt]+pixoffset, Nside_final, outputfile, 0, &coordsys);
	
#pragma omp parallel for
				for (long l = 0; l < Npix_final; l++)
					map_phi_final[outcnt][l+pixoffset] = (map_isw_final[l+pixoffset] < -1.5e29) ? map_isw_final[l+pixoffset] : map_isw_final[l+pixoffset] + 2.*(distances[outcnt]+0.5*dist-1.5*it0->second.hdr.distance) * (it0->second.pixel[l] - it1->second.pixel[l]) / (back[step].tau - back[step+1].tau);
				
				if (Nside_final > Nside_initial)
				{
					for (long l = 0; l < Npix_final; l++)
					{
						helper.SetNside(Nside_final, RING);
						ptg = helper.pix2ang(l);
						for (uint32_t Nside_interp = Nside_initial; Nside_interp < Nside_final; Nside_interp <<= 1)
						{
							helper.SetNside(Nside_interp, RING);
							helper.get_interpol(ptg, nnpix, nnwgt);
							
							q = 4l * (Nside_interp * Nside_interp - Nside_initial * Nside_initial);
							
							for (int nn = 0; nn < nnpix.size(); nn++)
							{
								if (map_isw_final[nnpix[nn]+q] > -1e30)
									map_phi_final[outcnt][l+pixoffset] += map_isw_final[nnpix[nn]+q] * nnwgt[nn];
								else
								{
									map_phi_final[outcnt][l+pixoffset] = -1.6375e30;
									break;
								}
							}
						}
					}
				}
	
				if (numoutputs > 1)
					sprintf(outputfile, "%s%s_isw_%d.fits", sim.output_path, sim.basename_lightcone, outcnt);
				else
					sprintf(outputfile, "%s%s_isw.fits", sim.output_path, sim.basename_lightcone);	
				write_healpix_map(map_phi_final[outcnt]+pixoffset, Nside_final, outputfile, 0, &coordsys);

				if (it0 == phi0.healpix_data.begin() || it1 == phi1.healpix_data.begin())
				{
					cout << COLORTEXT_YELLOW << " warning" << COLORTEXT_RESET << ": evaluating potential at boundary of covered range" << endl;
#pragma omp parallel for
					for (long l = 0; l < Npix_final; l++)
						map_phi_final[outcnt][l+pixoffset] = (map_isw_final[l+pixoffset] < -1.5e29) ? map_isw_final[l+pixoffset] : pot_obs - lin_int(it0->second.pixel[l],it1->second.pixel[l], (distances[outcnt] - (tauobs - back[step].tau))/((back[step].tau - back[step+1].tau)));
				}
				else if (std::prev(it0)->second.hdr.Nside != Nside_final)
				{
#pragma omp parallel for private(v1,q)
					for (long l = 0; l < Npix_final; l++)
					{
						pix2vec_ring64(Nside_final, l, v1);
						vec2pix_ring64(std::prev(it0)->second.hdr.Nside, v1, &q);
						map_phi_final[outcnt][l+pixoffset] = (map_isw_final[l+pixoffset] < -1.5e29) ? map_isw_final[l+pixoffset] : pot_obs - ((it0->second.hdr.distance-distances[outcnt]) * lin_int(it0->second.pixel[l],it1->second.pixel[l], (distances[outcnt] - (tauobs - back[step].tau))/((back[step].tau - back[step+1].tau))) + (distances[outcnt]-dist) * lin_int(std::prev(it0)->second.pixel[q],std::prev(it1)->second.pixel[q], (distances[outcnt] - (tauobs - back[step].tau))/((back[step].tau - back[step+1].tau)))) / (it0->second.hdr.distance - dist);
					}
				}
				else
				{
#pragma omp parallel for
					for (long l = 0; l < Npix_final; l++)
						map_phi_final[outcnt][l+pixoffset] = (map_isw_final[l+pixoffset] < -1.5e29) ? map_isw_final[l+pixoffset] : pot_obs - ((it0->second.hdr.distance-distances[outcnt]) * lin_int(it0->second.pixel[l],it1->second.pixel[l], (distances[outcnt] - (tauobs - back[step].tau))/((back[step].tau - back[step+1].tau))) + (distances[outcnt]-dist) * lin_int(std::prev(it0)->second.pixel[l],std::prev(it1)->second.pixel[l], (distances[outcnt] - (tauobs - back[step].tau))/((back[step].tau - back[step+1].tau)))) / (it0->second.hdr.distance - dist);
				}
	
				if (numoutputs > 1)
					sprintf(outputfile, "%s%s_potential_%d.fits", sim.output_path, sim.basename_lightcone, outcnt);
				else
					sprintf(outputfile, "%s%s_potential.fits", sim.output_path, sim.basename_lightcone);	
				write_healpix_map(map_phi_final[outcnt]+pixoffset, Nside_final, outputfile, 0, &coordsys);
	
				#pragma omp parallel for
				for (long l = 0; l < Npix_final; l++)
					map_phi_final[outcnt][l+pixoffset] = (map_shapiro_final[l+pixoffset] < -1.5e29) ? map_shapiro_final[l+pixoffset] : map_shapiro_final[l+pixoffset] - 2.*sim.boxsize*(distances[outcnt]+0.5*dist-1.5*it0->second.hdr.distance) * lin_int(it0->second.pixel[l],it1->second.pixel[l], (it0->second.hdr.distance - (tauobs - back[step].tau))/((back[step].tau - back[step+1].tau)));
					
				if (Nside_final > Nside_initial)
				{
					for (long l = 0; l < Npix_final; l++)
					{
						helper.SetNside(Nside_final, RING);
						ptg = helper.pix2ang(l);
						for (uint32_t Nside_interp = Nside_initial; Nside_interp < Nside_final; Nside_interp <<= 1)
						{
							helper.SetNside(Nside_interp, RING);
							helper.get_interpol(ptg, nnpix, nnwgt);
							
							q = 4l * (Nside_interp * Nside_interp - Nside_initial * Nside_initial);
							
							for (int nn = 0; nn < nnpix.size(); nn++)
							{
								if (map_shapiro_final[nnpix[nn]+q] > -1e30)
									map_phi_final[outcnt][l+pixoffset] += map_shapiro_final[nnpix[nn]+q] * nnwgt[nn];
								else
								{
									map_phi_final[outcnt][l+pixoffset] = -1.6375e30;
									break;
								}
							}
						}
					}
				}
				
				if (numoutputs > 1)
					sprintf(outputfile, "%s%s_shapiro_%d.fits", sim.output_path, sim.basename_lightcone, outcnt);	
				else
					sprintf(outputfile, "%s%s_shapiro.fits", sim.output_path, sim.basename_lightcone);
				write_healpix_map(map_phi_final[outcnt]+pixoffset, Nside_final, outputfile, 0, &coordsys);
				free(map_phi_final[outcnt]);
			
				outcnt++;
			}
			
			dist = it0->second.hdr.distance;
			
			if (cnt == thresh && it0->second.hdr.Npix == 12 * it0->second.hdr.Nside * it0->second.hdr.Nside)
			{
				for (int m = outcnt; m < numoutputs; m++)
				{
					monopole = 0;
				
#pragma omp parallel for reduction(+:monopole)			
					for (int l = pixoffset; l < pixoffset + it0->second.hdr.Npix; l++)
					{
						if (map_phi_final[m][l] > -1.5e29)
							monopole += map_phi_final[m][l];
					}
					
					monopole /= (double) it0->second.hdr.Npix;

#pragma omp parallel for				
					for (int l = pixoffset; l < pixoffset + it0->second.hdr.Npix; l++)
						map_phi_final[m][l] -= monopole;
				}
				
				thresh *= 2;
			}
			
			cnt++;
		}

		phi0.clear();
		phi1.clear();

		step++;
		
		if (step > numlines-2)
		{
			cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": reached cycle " << back[step].cycle << ", no further integration steps possible!" << endl;
			return -1;
		}
		
		phi0.cinfo = back[step];
		phi1.cinfo = back[step+1];
	}
	
	free(map_isw_final);
	free(map_shapiro_final);
	
	cout << COLORTEXT_GREEN << " normal completion." << COLORTEXT_RESET << endl;

	return 0;
}


int loadHealpixData(metric_container * field, double min_dist, double max_dist, char * lightconeparam)
{
	metric_data metric;
	char filename[1024];
	char tokens[128];
	char * token = NULL;
	uint32_t blocksize[2];
	FILE * infile = NULL;
	FILE * infile2 = NULL;
	int count;
	long backtrack;
	double vec[3];
	int64_t j, q;
	int ring;
	int pixbatch_delim[3];
	int pixbatch_size[3] = {0, 0, 0};


	metric.pixel = NULL;

	if (!field->healpix_data.empty())
	{
		if (field->healpix_data.begin()->second.hdr.distance < min_dist)
		{
			if ((min_dist = field->healpix_data.rbegin()->second.hdr.distance) > max_dist)
				return field->healpix_data.rbegin()->first;
		}	
		else if (field->healpix_data.begin()->second.hdr.distance > max_dist)
			max_dist = field->healpix_data.begin()->second.hdr.distance;
	}
	
	if (lightconeparam != NULL)
	{
		strcpy(tokens, lightconeparam);
		token = strtok(tokens, ",");
		do
		{
			sprintf(filename, "%s%s%d_%04d_%s.map", field->dir, field->basename, atoi(token), field->cinfo.cycle, field->name);
			infile = fopen(filename, "rb");
			token = strtok(NULL, ",");
		}
		while (infile == NULL && token != NULL);
	}
	else
	{
		sprintf(filename, "%s%s_%04d_%s.map", field->dir, field->basename, field->cinfo.cycle, field->name);
		infile = fopen(filename, "rb");
	}

	if (infile == NULL)
	{
		cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": map file " << filename << " could not be opened for reading!" << endl;
		return -1;
	}

	if (fread(blocksize+1, sizeof(uint32_t), 1, infile) != 1)
	{
		fclose(infile);
		return -1;
	}

	for (count = 0; !feof(infile) && !ferror(infile); count++)
	{
		if (blocksize[1] != 256)
		{
			fclose(infile);
			cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": invalid header block size in map file " << filename << "!" << endl;
			return -1;
		}

		if (fread(&metric.hdr, sizeof(metric.hdr), 1, infile) != 1)
		{
			fclose(infile);
			cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to read header block in map file " << filename << "!" << endl;
			return -1;
		}
		
		if (fread(blocksize, sizeof(uint32_t), 2, infile) != 2)
		{
			fclose(infile);
			cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": invalid block structure in map file " << filename << "!" << endl;
			return -1;
		}

		if (blocksize[0] != 256)
		{
			fclose(infile);
			cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": invalid header block size in map file " << filename << "!" << endl;
			return -1;
		}

		if (blocksize[1] != metric.hdr.precision * metric.hdr.Npix)
		{
			fclose(infile);
			cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": invalid data block size in map file " << filename << "!" << endl;
			return -1;
		}

		if (metric.hdr.distance > min_dist) break;

		backtrack = ftell(infile) - (256 + 2 * sizeof(uint32_t));
		
		if (infile2 != NULL)
		{
			fclose(infile2);
			infile2 = NULL;
		}

		if (fseek(infile, blocksize[1], SEEK_CUR))
		{
			cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to skip data block in map file " << filename << "!" << endl;
			fclose(infile);
			return -1;
		}

		if (fread(blocksize, sizeof(uint32_t), 2, infile) != 2)
		{	
			if (token != NULL && feof(infile))
			{
				infile2 = infile;
				sprintf(filename, "%s%s%d_%04d_%s.map", field->dir, field->basename, atoi(token), field->cinfo.cycle, field->name);
				infile = fopen(filename, "rb");
				token = strtok(NULL, ",");
				if (infile == NULL)
				{
					cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": map file " << filename << " could not be opened for reading!" << endl;
					return -1;
				}
				if (fread(blocksize+1, sizeof(uint32_t), 1, infile) != 1)
				{
					fclose(infile);
					return -1;
				}
			}
			else
			{
				fclose(infile);
				return -1;
			}
		}
	}

	if (feof(infile) || ferror(infile))
	{
		cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": read error occured in map file " << filename << "!" << endl;
		fclose(infile);
		return -1;
	}

	if (count > 0)
	{
		if (infile2 != NULL)
		{
			fclose(infile);
			infile = infile2;
		}
		
		if (fseek(infile, backtrack, SEEK_SET))
		{
			cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to backtrack in map file " << filename << "!" << endl;
			fclose(infile);
			return -1;
		}

		if (fread(&metric.hdr, sizeof(metric.hdr), 1, infile) != 1)
		{
			fclose(infile);
			cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to read header block in map file " << filename << "!" << endl;
			return -1;
		}
		
		if (fread(blocksize, sizeof(uint32_t), 2, infile) != 2)
		{
			fclose(infile);
			cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": invalid block structure in map file " << filename << "!" << endl;
			return -1;
		}

		if (blocksize[0] != 256)
		{
			fclose(infile);
			cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": invalid header block size in map file " << filename << "!" << endl;
			return -1;
		}

		if (blocksize[1] != metric.hdr.precision * metric.hdr.Npix)
		{
			fclose(infile);
			cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": invalid data block size in map file " << filename << "!" << endl;
			return -1;
		}

		count--;
	}

	while (true)
	{
		if (field->healpix_data.find(count) == field->healpix_data.end()) // data not present
		{
			if (metric.hdr.Nside < 2)
			{
				cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": invalid Nside = " << metric.hdr.Nside << " in map file " << filename << "!" << endl;
				fclose(infile);
				return -1;
			}

			if (metric.hdr.distance < 0)
			{
				cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": invalid distance = " << metric.hdr.distance << " in map file " << filename << "!" << endl;
				fclose(infile);
				return -1;
			}

			metric.pixel = (float *) malloc(metric.hdr.Npix * sizeof(float));

			if (metric.hdr.Nside_ring > 0 && metric.hdr.Nside_ring < metric.hdr.Nside)
			{	
				if ((long) metric.hdr.Npix <= 2 * (long) metric.hdr.Nside * (metric.hdr.Nside + 1))
					ring = (int) floor((sqrt(2. * metric.hdr.Npix + 1.01) - 1.) / 2.);
				else if ((long) metric.hdr.Npix <= 2 * (long) metric.hdr.Nside * (metric.hdr.Nside + 1) + 4 * (2 * metric.hdr.Nside - 1) * (long) metric.hdr.Nside)
					ring = ((metric.hdr.Npix - 2 * metric.hdr.Nside * (metric.hdr.Nside + 1)) / 4 / metric.hdr.Nside) + metric.hdr.Nside;
				else if ((long) metric.hdr.Npix < 12 * (long) metric.hdr.Nside * metric.hdr.Nside)
				{
					ring = 12 * (long) metric.hdr.Nside * metric.hdr.Nside - (long) metric.hdr.Npix;
					ring = (int) floor((sqrt(2. * ring + 1.01) - 1.) / 2.);
					ring = 4 * metric.hdr.Nside - 1 - ring;
				}
				else
					ring = 4 * metric.hdr.Nside - 1;
					
				pixbatch_size[0] = (metric.hdr.Nside / metric.hdr.Nside_ring);
						
				pixbatch_delim[1] = ring / pixbatch_size[0];
				pixbatch_delim[0] = (pixbatch_delim[1] > 0) ? pixbatch_delim[1]-1 : 0;
				pixbatch_delim[2] = pixbatch_delim[1]+1;
				pixbatch_size[1] = (pixbatch_size[0] * (pixbatch_size[0]+1) + (2*pixbatch_size[0] - 1 - ring%pixbatch_size[0]) * (ring%pixbatch_size[0])) / 2;
				pixbatch_size[2] = ((ring%pixbatch_size[0] + 1) * (ring%pixbatch_size[0])) / 2;
				pixbatch_size[0] *= pixbatch_size[0];
				for (int p = 0; p < 3; p++)
				{
					if (pixbatch_delim[p] <= metric.hdr.Nside_ring)
						pixbatch_delim[p] = 2 * pixbatch_delim[p] * (pixbatch_delim[p]+1);
					else if (pixbatch_delim[p] <= 3 * metric.hdr.Nside_ring)
						pixbatch_delim[p] = 2 * metric.hdr.Nside_ring * (metric.hdr.Nside_ring+1) + (pixbatch_delim[p]-metric.hdr.Nside_ring) * 4 * metric.hdr.Nside_ring;
					else if (pixbatch_delim[p] < 4 * metric.hdr.Nside_ring)
						pixbatch_delim[p] = 12 * metric.hdr.Nside_ring * metric.hdr.Nside_ring - 2 * (4 * metric.hdr.Nside_ring - 1 - pixbatch_delim[p]) * (4 * metric.hdr.Nside_ring - pixbatch_delim[p]);
					else
						pixbatch_delim[p] = 12 * metric.hdr.Nside_ring * metric.hdr.Nside_ring;
				}
					
				if (metric.hdr.precision == 4)
				{
					float * fpix = (float *) malloc (metric.hdr.Npix * sizeof(float));
					if (fread(fpix, sizeof(float), metric.hdr.Npix, infile) != metric.hdr.Npix)
					{
						cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to read data block in map file " << filename << "!" << endl;
						fclose(infile);
						free(metric.pixel);
						return -1;
					}

#pragma omp parallel for private(j) collapse(2)
					for (int p = 0; p < pixbatch_delim[0]; p++)
					{
						for (int i = 0; i < pixbatch_size[0]; i++)
						{
							ring2nest64(metric.hdr.Nside_ring, p, &j);
							j = j*pixbatch_size[0] + i;
							nest2ring64(metric.hdr.Nside, j, &j);
							metric.pixel[j] = fpix[pixbatch_size[0]*p + i];
						}
					}
#pragma omp parallel for private(q,j)
					for (int p = pixbatch_delim[0]; p < pixbatch_delim[1]; p++)
					{
						q = 0;
						for (int i = 0; i < pixbatch_size[0]; i++)
						{
							ring2nest64(metric.hdr.Nside_ring, p, &j);
							j = j*pixbatch_size[0] + i;
							nest2ring64(metric.hdr.Nside, j, &j);
							if (j < metric.hdr.Npix)
							{
								metric.pixel[j] = fpix[pixbatch_size[0]*pixbatch_delim[0] + pixbatch_size[1]*(p-pixbatch_delim[0]) + q];
								q++;
							}
						}
					}
#pragma omp parallel for private(q,j)
					for (int p = pixbatch_delim[1]; p < pixbatch_delim[2]; p++)
					{
						q = 0;
						for (int i = 0; i < pixbatch_size[0]; i++)
						{
							ring2nest64(metric.hdr.Nside_ring, p, &j);
							j = j*pixbatch_size[0] + i;
							nest2ring64(metric.hdr.Nside, j, &j);
							if (j < metric.hdr.Npix)
							{
								metric.pixel[j] = fpix[pixbatch_size[0]*pixbatch_delim[0] + pixbatch_size[1]*(pixbatch_delim[1]-pixbatch_delim[0]) + pixbatch_size[2]*(p-pixbatch_delim[1]) + q];
								q++;
							}
						}
					}
	
					free(fpix);
				}
				else if (metric.hdr.precision == 8)
				{
					double * dpix = (double *) malloc (metric.hdr.Npix * sizeof(double));
					if (fread(dpix, sizeof(double), metric.hdr.Npix, infile) != metric.hdr.Npix)
					{
						cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to read data block in map file " << filename << "!" << endl;
						fclose(infile);
						free(metric.pixel);
						return -1;
					}
	
#pragma omp parallel for private(j) collapse(2)
					for (int p = 0; p < pixbatch_delim[0]; p++)
					{
						for (int i = 0; i < pixbatch_size[0]; i++)
						{
							ring2nest64(metric.hdr.Nside_ring, p, &j);
							j = j*pixbatch_size[0] + i;
							nest2ring64(metric.hdr.Nside, j, &j);
							metric.pixel[j] = dpix[pixbatch_size[0]*p + i];
						}
					}
#pragma omp parallel for private(q,j)
					for (int p = pixbatch_delim[0]; p < pixbatch_delim[1]; p++)
					{
						q = 0;
						for (int i = 0; i < pixbatch_size[0]; i++)
						{
							ring2nest64(metric.hdr.Nside_ring, p, &j);
							j = j*pixbatch_size[0] + i;
							nest2ring64(metric.hdr.Nside, j, &j);
							if (j < metric.hdr.Npix)
							{
								metric.pixel[j] = dpix[pixbatch_size[0]*pixbatch_delim[0] + pixbatch_size[1]*(p-pixbatch_delim[0]) + q];
								q++;
							}
						}
					}
#pragma omp parallel for private(q,j)
					for (int p = pixbatch_delim[1]; p < pixbatch_delim[2]; p++)
					{
						q = 0;
						for (int i = 0; i < pixbatch_size[0]; i++)
						{
							ring2nest64(metric.hdr.Nside_ring, p, &j);
							j = j*pixbatch_size[0] + i;
							nest2ring64(metric.hdr.Nside, j, &j);
							if (j < metric.hdr.Npix)
							{
								metric.pixel[j] = dpix[pixbatch_size[0]*pixbatch_delim[0] + pixbatch_size[1]*(pixbatch_delim[1]-pixbatch_delim[0]) + pixbatch_size[2]*(p-pixbatch_delim[1]) + q];
								q++;
							}
						}
					}
	
					free(dpix);
				}
				else
				{
					cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": precision " << metric.hdr.precision << " bytes not supported for map files!" << endl;
					free(metric.pixel);
				}
			}
			else
			{
				if (metric.hdr.precision == 4)
				{
					if (fread(metric.pixel, sizeof(float), metric.hdr.Npix, infile) != metric.hdr.Npix)
					{
						cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to read data block in map file " << filename << "!" << endl;
						fclose(infile);
						free(metric.pixel);
						return -1;
					}
				}
				else if (metric.hdr.precision == 8)
				{
					double * dpix = (double *) malloc (metric.hdr.Npix * sizeof(double));
					if (fread(dpix, sizeof(double), metric.hdr.Npix, infile) != metric.hdr.Npix)
					{
						cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to read data block in map file " << filename << "!" << endl;
						fclose(infile);
						free(metric.pixel);
						return -1;
					}
					for (int p = 0; p < metric.hdr.Npix; p++)
						metric.pixel[p] = dpix[p];
					free(dpix);
				}
				else
				{
					cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": precision " << metric.hdr.precision << " bytes not supported for map files!" << endl;
					free(metric.pixel);
				}
			}
			
			while (metric.hdr.Nside > 8192)
			{
				float * fpix = (float *) malloc (metric.hdr.Npix * sizeof(float) / 4);
#pragma omp parallel for
				for (int p = 0; p < metric.hdr.Npix/4; p++)
					fpix[p] = 0.;
					
#pragma omp parallel for private(j)
				for (int p = 0; p < metric.hdr.Npix; p++)
				{
					ring2nest64(metric.hdr.Nside, p, &j);
					j /= 4;
					nest2ring64(metric.hdr.Nside/2, j, &j);
					if (j < metric.hdr.Npix/4)
					{
						if (metric.pixel[p] > -1.e30)
							fpix[j] += metric.pixel[p] / 4.;
						else
							fpix[j] = -1.6375e30;
					}
				}
				
				free(metric.pixel);
				metric.pixel = fpix;
				metric.hdr.Nside /= 2;
				metric.hdr.Npix /= 4;
			}

			field->healpix_data.insert(std::pair<int,metric_data>(count, metric));
		}
		else
		{
			if (fseek(infile, blocksize[1], SEEK_CUR))
			{
				cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to skip data block in map file " << filename << "!" << endl;
				fclose(infile);
				return -1;
			}
		}

		if (metric.hdr.distance > max_dist) break;

		if (fread(blocksize, sizeof(uint32_t), 2, infile) != 2)
		{
			if (infile2 != NULL)
			{
				fclose(infile2);
				infile2 = NULL;
				infile = fopen(filename, "rb");
				if (infile == NULL)
				{
					cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": map file " << filename << " could not be opened for reading!" << endl;
					return -1;
				}
				if (fread(blocksize+1, sizeof(uint32_t), 1, infile) != 1)
				{
					fclose(infile);
					return -1;
				}
			}
			else if (token != NULL && feof(infile))
			{
				fclose(infile);
				sprintf(filename, "%s%s%d_%04d_%s.map", field->dir, field->basename, atoi(token), field->cinfo.cycle, field->name);
				infile = fopen(filename, "rb");
				token = strtok(NULL, ",");
				if (infile == NULL)
				{
					cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": map file " << filename << " could not be opened for reading!" << endl;
					return -1;
				}
				if (fread(blocksize+1, sizeof(uint32_t), 1, infile) != 1)
				{
					fclose(infile);
					return -1;
				}
			}
			else
			{
				fclose(infile);
				return -2;
			}
		}

		if (blocksize[1] != 256)
		{
			fclose(infile);
			cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": invalid header block size in map file " << filename << "!" << endl;
			return -1;
		}

		if (fread(&metric.hdr, sizeof(metric.hdr), 1, infile) != 1)
		{
			fclose(infile);
			cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to read header block in map file " << filename << "!" << endl;
			return -1;
		}
		
		if (fread(blocksize, sizeof(uint32_t), 2, infile) != 2)
		{
			fclose(infile);
			cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": invalid block structure in map file " << filename << "!" << endl;
			return -1;
		}

		if (blocksize[0] != 256)
		{
			fclose(infile);
			cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": invalid header block size in map file " << filename << "!" << endl;
			return -1;
		}

		if (blocksize[1] != metric.hdr.precision * metric.hdr.Npix)
		{
			fclose(infile);
			cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": invalid data block size in map file " << filename << "!" << endl;
			return -1;
		}

		count++;
	}

	fclose(infile);

	return count;
}

bool kappa(float * pixel, const int64_t Nside, int64_t ipix, float & result)
{
	int64_t j, k, l, q, ring;
	float temp, temp2, w1, w2;

	if (pixel[ipix] < -1e30) return false;

	if (ipix < Nside * (Nside + 1) * 2l) // north polar cap
	{
		ring = (1 + (int64_t) sqrt(1.5 + 2 * (double) ipix)) / 2;
		j = ipix - 2 * ring * (ring-1);
		q = j / ring;
		j %= ring;
		
		// phi-derivative
		k = ipix+1;
		l = ipix-1;
		if (q == 3 && j == ring-1)
			k -= 4*ring;
		if (q == 0 && j == 0)
			l += 4*ring;
		
		if (pixel[k] < -1e30 || pixel[l] < -1e30) return false;
	
		temp2 = (pixel[k] + pixel[l] - 2. * pixel[ipix]);
		
		// ring derivative
		if (ring == Nside)
		{
			k = ring * (ring+1) * 2l + q * ring + j;
			l = k+1;
			
			if (q == 3 && j == Nside-1)
				l -= 4*Nside;
		
			if (pixel[k] < -1e30 || pixel[l] < -1e30) return false;
		
			result = (0.5 * (pixel[k] + pixel[l]) - pixel[ipix]) * 8. / 3.;
			temp = (0.5 * (pixel[k] + pixel[l]));
			w1 = 0.125;
		}
		else
		{
			k = ring * (ring+1) * 2l + q * (ring+1) + j;
			l = k+1;
		
			if (pixel[k] < -1e30 || pixel[l] < -1e30) return false;
		
			result = (1. - (j+0.5)/ring) * pixel[k] + ((j+0.5)/ring) * pixel[l];
			temp = result;
			w1 = (ring-j-0.5) * (j+0.5) * 0.5 / (ring+1) / (ring+1);
		}
	
		if (ring == 1)
		{
			if (pixel[0] < -1e30 || pixel[1] < -1e30 || pixel[2] < -1e30 || pixel[3] < -1e30) return false;
				
			result -= 0.25 * (pixel[0] + pixel[1] + pixel[2] + pixel[3]);
			temp += 0.25 * (pixel[0] + pixel[1] + pixel[2] + pixel[3]);
			w2 = 0;
		}
		else
		{
			k = (ring-2) * (ring-1) * 2l + q * (ring-1) + j;
			l = k-1;
			if (q == 0 && j == 0)
				l += 4*(ring-1);
			if (q == 3 && j == ring-1)
				k -= 4*(ring-1);
				
			if (pixel[k] < -1e30 || pixel[l] < -1e30) return false;
			
			if (ring == Nside)
				result += pixel[ipix] - ((1. - (j+0.5)/ring) * pixel[k] + ((j+0.5)/ring) * pixel[l]);
			else
				result -= (1. - (j+0.5)/ring) * pixel[k] + ((j+0.5)/ring) * pixel[l];
			temp += (1. - (j+0.5)/ring) * pixel[k] + ((j+0.5)/ring) * pixel[l];
			
			w2 = (ring-j-0.5) * (j+0.5) * 0.5 / (ring-1) / (ring-1);
		}
		
		result -= (w1 - w2) * temp2;
	
		result *= (6 * Nside * Nside - 3 * ring * ring) / 8. / ring;
		result += (6 * Nside * Nside - ring * ring) * (temp - 2. * pixel[ipix] - (w1 + w2) * temp2) / 4.;
	
		result += 36. * Nside * Nside * Nside * Nside * temp2 / M_PI / M_PI / (6 * Nside * Nside - ring * ring);
	}
	else if (ipix < 2l * Nside * (5l * Nside - 1l)) // equatorial region
	{
		ring = (ipix - 2l * Nside * (Nside-1)) / (4l * Nside); // + Nside
		j = ipix - 2l * Nside * (Nside-1) - 4l * Nside * ring;
		
		k = (j == 4l*Nside-1) ? ipix+1-4l*Nside : ipix+1;
		l = (j == 0) ? ipix+4l*Nside-1 : ipix-1;
		
		if (pixel[k] < -1e30 || pixel[l] < -1e30) return false;
		
		temp2 = (pixel[k] + pixel[l] - 2.*pixel[ipix]);
		
		k = ipix + 4l * Nside;
		
		if (ring % 2)
		{
			l = (j == 0) ? k+(4l*Nside-1) : k-1;
			
			if (pixel[k] < -1e30 || pixel[l] < -1e30) return false;
			
			result = 0.5 * (pixel[k] + pixel[l]);
			temp = result;
			
			k = ipix - 4l * Nside;
			l = (j == 0) ? ipix-1 : k-1;
		}
		else
		{
			l = (j == 4l*Nside-1) ? ipix+1 : k+1;
			
			if (pixel[k] < -1e30 || pixel[l] < -1e30) return false;
			
			result = 0.5 * (pixel[k] + pixel[l]);
			temp = result;
			
			k = ipix - 4l * Nside;
			l = (j == 4l*Nside-1) ? k+1-4l*Nside : k+1;
		}
		
		if (pixel[k] < -1e30 || pixel[l] < -1e30) return false;
			
		result -= 0.5 * (pixel[k] + pixel[l]);
		temp += 0.5 * (pixel[k] + pixel[l]);
		
		result *= Nside-ring;
		result += (temp - 2.*pixel[ipix] - 0.25 * temp2) * (2.25*Nside*Nside - (Nside-ring)*(Nside-ring));
		
		result += temp2 * 4. * Nside * Nside / M_PI / M_PI / (1. - (Nside-ring)*(Nside-ring)/2.25/Nside/Nside);
	}
	else // south polar cap
	{
		ring = (1 + (int64_t) sqrt(1.5 + 2 * (double) (12l*Nside*Nside-1 - ipix))) / 2;
		j = 12l*Nside*Nside-1 - ipix - 2 * ring * (ring-1);
		q = j / ring;
		j %= ring;
		
		// phi-derivative
		k = ipix+1;
		l = ipix-1;
		if (q == 3 && j == ring-1)
			l += 4*ring;
		if (q == 0 && j == 0)
			k -= 4*ring;
		
		if (pixel[k] < -1e30 || pixel[l] < -1e30) return false;
	
		temp2 = (pixel[k] + pixel[l] - 2. * pixel[ipix]);
		
		// ring derivative
		if (ring == Nside)
		{
			k = ring * (ring+1) * 2l + q * ring + j;
			l = k-1;
			if (q == 0 && j == 0)
				l += 4*Nside;
		
			if (pixel[12l*Nside*Nside-1-k] < -1e30 || pixel[12l*Nside*Nside-1-l] < -1e30) return false;
		
			result = (0.5 * (pixel[12l*Nside*Nside-1-k] + pixel[12l*Nside*Nside-1-l]) - pixel[ipix]) * 8. / 3.;
			temp = (0.5 * (pixel[12l*Nside*Nside-1-k] + pixel[12l*Nside*Nside-1-l]));
			w1 = 0.125;
		}
		else
		{
			k = ring * (ring+1) * 2l + q * (ring+1) + j;
			l = k+1;
		
			if (pixel[12l*Nside*Nside-1-k] < -1e30 || pixel[12l*Nside*Nside-1-l] < -1e30) return false;
		
			result = (1. - (j+0.5)/ring) * pixel[12l*Nside*Nside-1-k] + ((j+0.5)/ring) * pixel[12l*Nside*Nside-1-l];
			temp = result;
			w1 = (ring-j-0.5) * (j+0.5) * 0.5 / (ring+1) / (ring+1);
		}
	
		if (ring == 1)
		{
			if (pixel[12l*Nside*Nside-1] < -1e30 || pixel[12l*Nside*Nside-2] < -1e30 || pixel[12l*Nside*Nside-3] < -1e30 || pixel[12l*Nside*Nside-4] < -1e30) return false;
				
			result -= 0.25 * (pixel[12l*Nside*Nside-1] + pixel[12l*Nside*Nside-2] + pixel[12l*Nside*Nside-3] + pixel[12l*Nside*Nside-4]);
			temp += 0.25 * (pixel[12l*Nside*Nside-1] + pixel[12l*Nside*Nside-2] + pixel[12l*Nside*Nside-3] + pixel[12l*Nside*Nside-4]);
			w2 = 0;
		}
		else
		{
			k = (ring-2) * (ring-1) * 2l + q * (ring-1) + j;
			l = k-1;
			if (q == 0 && j == 0)
				l += 4*(ring-1);
			if (q == 3 && j == ring-1)
				k -= 4*(ring-1);
				
			if (pixel[12l*Nside*Nside-1-k] < -1e30 || pixel[12l*Nside*Nside-1-l] < -1e30) return false;
			
			if (ring == Nside)
				result += pixel[ipix] - ((1. - (j+0.5)/ring) * pixel[12l*Nside*Nside-1-k] + ((j+0.5)/ring) * pixel[12l*Nside*Nside-1-l]);
			else
				result -= (1. - (j+0.5)/ring) * pixel[12l*Nside*Nside-1-k] + ((j+0.5)/ring) * pixel[12l*Nside*Nside-1-l];
			temp += (1. - (j+0.5)/ring) * pixel[12l*Nside*Nside-1-k] + ((j+0.5)/ring) * pixel[12l*Nside*Nside-1-l];
			w2 = (ring-j-0.5) * (j+0.5) * 0.5 / (ring-1) / (ring-1);
		}

		result -= (w1 - w2) * temp2;
		result *= (6 * Nside * Nside - 3 * ring * ring) / 8. / ring;
		result += (6 * Nside * Nside - ring * ring) * (temp - 2. * pixel[ipix] - (w1 + w2) * temp2) / 4.;
		
		result += 36. * Nside * Nside * Nside * Nside * temp2  / M_PI / M_PI / (6 * Nside * Nside - ring * ring);
	}
	
	result /= -2.;
	
	return true;
}




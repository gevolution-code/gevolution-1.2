#include <stdlib.h>
#include <iostream>
#include <cmath>
#include "metadata.hpp"
#include "parser.hpp"

using namespace std;

int main(int argc, char **argv)
{
	char * settingsfile = NULL;
	char * cycleparam = NULL;
	char * lightconeparam = NULL;
	parameter * params = NULL;
	metadata sim;
	cosmology cosmo;
	icsettings ic;
	int numparam = 0;
	int usedparams = 0;

	bool use_lightcone[MAX_OUTPUTS];
	int min_cycle = 0;
	int max_cycle = 0;

	char filename[1024];
	char ofilename[1024];
	FILE * infile;
	FILE * outfile;

	uint64_t numpart_tot = 0;
	uint64_t numpart_write = 0;
	uint32_t blocksize = 0;
	int numfiles = 1;
	long numread = 0;
	long numwrite = 0;
	gadget2_header hdr;
	gadget2_header outhdr;

	double * vertex = NULL;
	double z_obs = -2.;

	long backtrack;
	long fastforward;
	float * posbatch = NULL;
	float * velbatch = NULL;
	uint32_t batch;
#if GADGET_ID_BYTES == 8
	uint64_t * IDbatch = NULL;
#else
	uint32_t * IDbatch = NULL;
#endif

	double offset = 0.;
	
	if (argc < 2)
	{
		cout << COLORTEXT_WHITE << " LCARS tools: lccat" << COLORTEXT_RESET << endl;
		cout << " catenates particle light-cone output and generates Gadget-2 binaries" << endl << endl;
		
		cout << " List of command-line options:" << endl;
		cout << " -s <filename>       : gevolution settings file of the simulation (mandatory)" << endl;
		cout << " -c <min>-<max>      : range of simulation cycles for which output should be" << endl;
		cout << "                       catenated (mandatory)" << endl;
		cout << " -l <ID1>[,<ID2>,...]: IDs of light cones to be included (optional, must refer" << endl;
		cout << "                       to the same observation event)" << endl;
		cout << " -n <numfiles>       : number of Gadget-2 binaries to distribute the catenated" << endl;
		cout << "                       output over (optional, default 1)" << endl;
		cout << " -o <offset>         : constant offset added to all particle coordinates, e.g." << endl;
		cout << "                       in order to avoid negative values (optional, default 0)" << endl << endl;
		cout << " The output will be written to <numfiles> approximately equal-sized Gadget-2" << endl;
		cout << " binaries that follow the naming conventions specified in the settings file." << endl;
		return 0;
	}

	for (int i = 1 ; i < argc ; i++ ){
		if ( argv[i][0] != '-' )
			continue;
		switch(argv[i][1]) {
			case 'c':
				cycleparam = argv[++i]; // cycle range selector
				break;
			case 'l':
				lightconeparam = argv[++i]; // light cone selector
				break;
			case 'n':
				numfiles = atoi(argv[++i]); // number of output files
				break;
			case 'o':
				offset = atof(argv[++i]); // position offset (to make resulting contiguous data region fit into the final cube)
				break;
			case 's':
				settingsfile = argv[++i]; //settings file name
				break;
			default:
				cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unknown command-line parameter " << argv[i] << endl << " call lccat without arguments to display help" << endl;
				return -1;
		}
	}

	if (settingsfile == NULL)
	{
		cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": no settings file specified!" << endl;
		return -1;
	}
	else if (numfiles < 1 || !isfinite(numfiles))
	{
		cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": number of output files not recognized!" << endl;
		return -1;
	}
	
	cout << COLORTEXT_WHITE << " LCARS tools: lccat" << COLORTEXT_RESET << endl << endl << " opening settings file of simulation: " << settingsfile << endl << " parser output:" << endl << endl;

	numparam = loadParameterFile(settingsfile, params);
	usedparams = parseMetadata(params, numparam, sim, cosmo, ic);
	free(params);

	cout << endl << " file contains " << numparam << " parameters, " << usedparams << " of which could be parsed." << endl << endl;

	offset /= GADGET_LENGTH_CONVERSION;

	cout << " number of lightcones: " << sim.num_lightcone << endl << endl;

	for (int i = 0; i < sim.num_lightcone; i++)
	{
		cout << " lightcone " << i << " parameters:" << endl << "  vertex = " << sim.lightcone[i].vertex[0] << ", " << sim.lightcone[i].vertex[1] << ", " << sim.lightcone[i].vertex[2] << endl;
		cout << "  redshift of observation = " << sim.lightcone[i].z << endl;
		cout << "  direction = " << sim.lightcone[i].direction[0] << ", " << sim.lightcone[i].direction[1] << ", " << sim.lightcone[i].direction[2] << endl;
		cout << "  opening half-angle = " << ((sim.lightcone[i].opening > -1.) ? acos(sim.lightcone[i].opening) * 180. / M_PI : 180.) << " degrees" << endl;
		cout << "  distance interval = " << sim.lightcone[i].distance[0] << ", " << sim.lightcone[i].distance[1] << endl << endl;
	}

	if (lightconeparam == NULL)
	{
		cout << " no light cones selected (parameter -l), using all light cones with particles: ";
		for (int i = 0; i < sim.num_lightcone; i++)
		{
			if (sim.out_lightcone[i] & MASK_GADGET)
			{
				cout << i << " ";
				use_lightcone[i] = true;

				if (z_obs < -1.)
				{
					vertex = sim.lightcone[i].vertex;
					z_obs = sim.lightcone[i].z;
				}
				else if (sim.lightcone[i].z != z_obs || sim.lightcone[i].vertex[0] != vertex[0] || sim.lightcone[i].vertex[1] != vertex[1] || sim.lightcone[i].vertex[2] != vertex[2])
				{
					cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": not all selected light cones refer to the same observation event!" << endl << " incompatible lightcone " << i << " will be omitted!" << endl << endl;
					use_lightcone[i] = false;
				}
			}
			else use_lightcone[i] = false;
		}
		cout << endl << endl;
	}
	else
	{
		for (int i = 0; i < sim.num_lightcone; i++) use_lightcone[i] = false;
		
		char * token = NULL;

		sprintf(filename, "%s", lightconeparam);

		token = strtok(filename, ",");

		while (token != NULL)
		{
			int i = atoi(token);

			if (i < 0 || i >= sim.num_lightcone)
			{
				cout << " light cone " << i << " (parameter -l) does not exist!" << endl;
			}
			else if (!(sim.out_lightcone[i] & MASK_GADGET))
			{
				cout << " light cone " << i << " contains no particles and will be omitted!" << endl;
			}
			else
			{
				use_lightcone[i] = true;

				if (z_obs < -1.)
				{
					vertex = sim.lightcone[i].vertex;
					z_obs = sim.lightcone[i].z;
				}
				else if (sim.lightcone[i].z != z_obs || sim.lightcone[i].vertex[0] != vertex[0] || sim.lightcone[i].vertex[1] != vertex[1] || sim.lightcone[i].vertex[2] != vertex[2])
				{
					cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": not all selected light cones refer to the same observation event!" << endl << " incompatible lightcone " << i << " will be omitted!" << endl << endl;
					use_lightcone[i] = false;
				}
			}

			token = strtok(NULL, ",");
		}

		cout << " using light cone(s) ";
		for (int i = 0; i < sim.num_lightcone; i++)
		{
			if (use_lightcone[i]) cout << i << " ";
		}
		cout << endl << endl;
	}

	if (cycleparam == NULL)
	{
		cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": no range of cycles selected (parameter -c)" << endl;
		return -1;
	}
	else
	{
		sprintf(filename, "%s", cycleparam);
		char * dash = strchr(filename, '-');
		int tmp;

		if (dash == NULL)
		{
			cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": range of cycle (parameter -c) could not be interpreted!" << endl;
			return -1;
		}

		max_cycle = atoi(dash+1);

		dash[0] = '\0';
		min_cycle = atoi(filename);

		cout << " range of cycles set to: " << min_cycle << "-" << max_cycle << endl;
	}

	cout << " reading particle headers..." << endl << endl;

	for (int cycle = min_cycle; cycle <= max_cycle; cycle++)
	{
		for (int i = 0; i < sim.num_lightcone; i++)
		{
			if (!use_lightcone[i]) continue;

			if (sim.num_lightcone > 1)
				sprintf(filename, "%s%s%d_%04d_cdm", sim.output_path, sim.basename_lightcone, i, cycle);
			else
				sprintf(filename, "%s%s_%04d_cdm", sim.output_path, sim.basename_lightcone, cycle);

			infile = fopen(filename, "rb");

			if (infile != NULL)
			{
				fread(&blocksize, sizeof(uint32_t), 1, infile);

				if (blocksize != sizeof(hdr))
				{
					cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unknown file format " << filename << "!" << endl;
					fclose(infile);
					continue;
				}

				if(fread(&hdr, sizeof(hdr), 1, infile) != 1)
				{
					cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to read header block from " << filename << "!" << endl;
					fclose(infile);
					continue;
				}

				numpart_tot += (uint64_t) hdr.npartTotal[1] + ((uint64_t) hdr.npartTotalHW[1] << 32);
				numread++;
				fclose(infile);
			}
		}
	}

	cout << " " << numread << " particle headers read successfully. Total number of particles = " << numpart_tot << endl << endl;

	for (int i = 0; i < 6; i++)
	{
		outhdr.npart[i] = 0;
		outhdr.mass[i] = 0.;
		outhdr.npartTotal[i] = 0;
		outhdr.npartTotalHW[i] = 0;
	}
	for (int i = 0; i < 256 - 6 * 4 - 6 * 8 - 2 * 8 - 2 * 4 - 6 * 4 - 2 * 4 - 4 * 8 - 2 * 4 - 6 * 4; i++)
		outhdr.fill[i] = 0;

	outhdr.Omega0 = cosmo.Omega_m;
	outhdr.OmegaLambda = 1. - cosmo.Omega_m;
	outhdr.HubbleParam = cosmo.h;
	outhdr.BoxSize = sim.boxsize / GADGET_LENGTH_CONVERSION;
	outhdr.flag_sfr = 0;
	outhdr.flag_cooling = 0;
	outhdr.flag_feedback = 0;
	outhdr.flag_age = 0;
	outhdr.flag_metals = 0;
	outhdr.time = 1. / (z_obs + 1.);
	outhdr.redshift = z_obs;

	outhdr.npart[1] = (uint32_t) ((((int64_t) numpart_tot) / numfiles) % (1ll << 32));
	outhdr.npartTotal[1] = (uint32_t) (((int64_t) numpart_tot) % (1ll << 32));
	outhdr.npartTotalHW[1] = (uint32_t) (((int64_t) numpart_tot) >> 32);
	outhdr.mass[1] = hdr.mass[1];
	outhdr.num_files = numfiles;
	
	posbatch = (float *) malloc(3 * (outhdr.npart[1] + numfiles) * sizeof(float));
	velbatch = (float *) malloc(3 * (outhdr.npart[1] + numfiles) * sizeof(float));
#if GADGET_ID_BYTES == 8
	IDbatch = (uint64_t *) malloc((outhdr.npart[1] + numfiles) * sizeof(uint64_t));
#else
	IDbatch = (uint32_t *) malloc((outhdr.npart[1] + numfiles) * sizeof(uint32_t));
#endif

	if (posbatch == NULL || velbatch == NULL || IDbatch == NULL)
	{
		cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to allocate memory for data read!" << endl;
		return -1;
	}

	cout << " building up particle light cone..." << endl << endl;

	numread = 0;

	for (int cycle = min_cycle; cycle <= max_cycle; cycle++)
	{
		cout << " cycle " << cycle << " ..." << endl;
		
		for (int i = 0; i < sim.num_lightcone; i++)
		{
			if (!use_lightcone[i]) continue;

			if (sim.num_lightcone > 1)
				sprintf(filename, "%s%s%d_%04d_cdm", sim.output_path, sim.basename_lightcone, i, cycle);
			else
				sprintf(filename, "%s%s_%04d_cdm", sim.output_path, sim.basename_lightcone, cycle);

			infile = fopen(filename, "rb");

			if (infile != NULL)
			{
				fread(&blocksize, sizeof(uint32_t), 1, infile);

				if (blocksize != sizeof(hdr))
				{
					cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unknown file format " << filename << "!" << endl;
					fclose(infile);
					continue;
				}

				if(fread(&hdr, sizeof(hdr), 1, infile) != 1)
				{
					cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to read header block from " << filename << "!" << endl;
					fclose(infile);
					continue;
				}
				
				fread(&blocksize, sizeof(uint32_t), 1, infile);
				fread(&blocksize, sizeof(uint32_t), 1, infile);

				long blockoffset = 3l * sizeof(float) * (long) hdr.npart[1] + 2l * sizeof(uint32_t);

				backtrack = ftell(infile);

				if (fseek(infile, 2 * blockoffset, SEEK_CUR))
				{
					cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to fast forward to ID block in " << filename << "!" << endl;
					fclose(infile);
					continue;
				}

				for (int64_t p = 0; p < (int64_t) hdr.npartTotal[1] + ((int64_t) hdr.npartTotalHW[1] << 32); p += batch)
				{
					batch = ((int64_t) hdr.npartTotal[1] + ((int64_t) hdr.npartTotalHW[1] << 32) - p >= (int64_t) outhdr.npart[1] - numread) ? (outhdr.npart[1] - (uint32_t) numread) : (uint32_t) ((int64_t) hdr.npartTotal[1] + ((int64_t) hdr.npartTotalHW[1] << 32) - p);

					if (
#if GADGET_ID_BYTES == 8
					fread(IDbatch+numread, sizeof(uint64_t), batch, infile)
#else
					fread(IDbatch+numread, sizeof(uint32_t), batch, infile)
#endif
					!= batch)
					{
						cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to read ID batch from " << filename << "!" << endl;
						fclose(infile);
						return -1;
					}

					fastforward = ftell(infile);

					if (fseek(infile, backtrack, SEEK_SET))
					{
						cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to rewind to positions block in " << filename << "!" << endl;
						fclose(infile);
						return -1;
					}

					if (fread(posbatch+3l*numread, sizeof(float), 3l*batch, infile) != 3l*batch)
					{
						cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to read position data from " << filename << "!" << endl;
						fclose(infile);
						return -1;
					}

					backtrack = ftell(infile);

					if (fseek(infile, blockoffset - 3l * batch * sizeof(float), SEEK_CUR))
					{
						cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to fast forward to velocities block in " << filename << "!" << endl;
						fclose(infile);
						return -1;
					}

					if (fread(velbatch+3l*numread, sizeof(float), 3l*batch, infile) != 3l*batch)
					{
						cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to read velocity data from " << filename << "!" << endl;
						fclose(infile);
						return -1;
					}

					if (fseek(infile, fastforward, SEEK_SET))
					{
						cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to fast forward to ID block in " << filename << "!" << endl;
						fclose(infile);
						return -1;
					}

					numread += batch;

					if ((numread == outhdr.npart[1] && numwrite < numfiles-1) || (uint64_t) numread + numpart_write == numpart_tot)
					{
						for (long q = 0; q < 3l * numread; q++)
						{
							posbatch[q] += offset;
							if (posbatch[q] > outhdr.BoxSize) outhdr.BoxSize = posbatch[q];
						}

						if (numfiles > 1)
							sprintf(ofilename, "%s%s_cdm.%ld", sim.output_path, sim.basename_lightcone, numwrite);
						else
							sprintf(ofilename, "%s%s_cdm", sim.output_path, sim.basename_lightcone);

						outfile = fopen(ofilename, "wb");
	
						if (outfile == NULL)
						{
							fclose(infile);
							cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to open file " << ofilename << " for output!" << endl;
							return -1;
						}
	
						blocksize = sizeof(outhdr);

						cout << COLORTEXT_CYAN << " writing" << COLORTEXT_RESET << " output file " << ofilename << " ..." << endl;

						fwrite(&blocksize, sizeof(uint32_t), 1, outfile);
						fwrite(&outhdr, sizeof(outhdr), 1, outfile);
						fwrite(&blocksize, sizeof(uint32_t), 1, outfile);

						blocksize = 3l * outhdr.npart[1] * sizeof(float);

						fwrite(&blocksize, sizeof(uint32_t), 1, outfile);

						if (fwrite(posbatch, sizeof(float), 3l * outhdr.npart[1], outfile) != 3l * outhdr.npart[1])
						{
							fclose(infile);
							fclose(outfile);
							cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to write position block!" << endl;
							return -1;
						}

						fwrite(&blocksize, sizeof(uint32_t), 1, outfile);
						fwrite(&blocksize, sizeof(uint32_t), 1, outfile);

						if (fwrite(velbatch, sizeof(float), 3l * outhdr.npart[1], outfile) != 3l * outhdr.npart[1])
						{
							fclose(infile);
							fclose(outfile);
							cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to write velocity block!" << endl;
							return -1;
						}

						fwrite(&blocksize, sizeof(uint32_t), 1, outfile);

						blocksize = outhdr.npart[1] * GADGET_ID_BYTES;

						fwrite(&blocksize, sizeof(uint32_t), 1, outfile);

						if (fwrite(IDbatch, GADGET_ID_BYTES, outhdr.npart[1], outfile) != outhdr.npart[1])
						{
							fclose(infile);
							fclose(outfile);
							cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to write ID block!" << endl;
							return -1;
						}

						fwrite(&blocksize, sizeof(uint32_t), 1, outfile);

						fclose(outfile);

						cout << " output file written, contains " << outhdr.npart[1] << " particles." << endl;
	
						numpart_write += outhdr.npart[1];
						numwrite++;

						if (numwrite == numfiles-1)
							outhdr.npart[1] += (uint32_t) (((long) numpart_tot) % numfiles);

						numread = 0;
					}
				}
				
				fclose(infile);
			}
		}
	}

	cout << endl << COLORTEXT_GREEN << " particle light cone complete." << endl << COLORTEXT_RESET << endl;

	free(IDbatch);
	free(posbatch);
	free(velbatch);

	if (outhdr.BoxSize != sim.boxsize / GADGET_LENGTH_CONVERSION)
	{
		cout << " correcting header information (BoxSize = " << outhdr.BoxSize << ") ..." << endl;

		outhdr.npart[1] = (uint32_t) ((((long) numpart_tot) / numfiles) % (1ll << 32));
		blocksize = sizeof(outhdr);

		for (int i = 0; i < numwrite-1; i++)
		{
			sprintf(ofilename, "%s%s_cdm.%d", sim.output_path, sim.basename_lightcone, i);

			outfile = fopen(ofilename, "r+b");

			if (outfile == NULL)
			{
				cout << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": unable to open file " << ofilename << " for output!" << endl;
				return -1;
			}

			fwrite(&blocksize, sizeof(uint32_t), 1, outfile);
			fwrite(&outhdr, sizeof(outhdr), 1, outfile);

			fclose(outfile);
		}
	}

	cout << endl << COLORTEXT_GREEN << " normal completion." << COLORTEXT_RESET << endl;

	return 0;
}


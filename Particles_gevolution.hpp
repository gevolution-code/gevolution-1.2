//////////////////////////
// Particles_gevolution.hpp
//////////////////////////
//
// Author: Julian Adamek (Université de Genève & Observatoire de Paris & Queen Mary University of London)
//
// Last modified: March 2020
//
//////////////////////////

#ifndef PARTICLES_GEVOLUTION_HEADER
#define PARTICLES_GEVOLUTION_HEADER

#ifndef PCLBUFFER
#define PCLBUFFER 1048576
#endif

using namespace LATfield2;

template <typename part, typename part_info, typename part_dataType>
class Particles_gevolution: public Particles<part, part_info, part_dataType>
{
	public:
		void saveGadget2(string filename, gadget2_header & hdr, const int tracer_factor = 1, double dtau_pos = 0., double dtau_vel = 0., Field<Real> * phi = NULL);
		void saveGadget2(string filename, gadget2_header & hdr, lightcone_geometry & lightcone, double dist, double dtau, double dtau_old, double dadtau, double vertex[MAX_INTERSECTS][3], const int vertexcount, set<long> & IDbacklog, set<long> & IDprelog, Field<Real> * phi, const int tracer_factor = 1);
		void loadGadget2(string filename, gadget2_header & hdr);
};

template <typename part, typename part_info, typename part_dataType>
void Particles_gevolution<part,part_info,part_dataType>::saveGadget2(string filename, gadget2_header & hdr, const int tracer_factor, double dtau_pos, double dtau_vel, Field<Real> * phi)
{
	float * posdata;
	float * veldata;
	void * IDs;
	MPI_File outfile;
	long count, npart;
	MPI_Offset offset_pos, offset_vel, offset_ID;
	MPI_Status status;
	uint32_t blocksize;
	uint32_t i;
	char fname[filename.length()+8];
	double rescale_vel = 1. / sqrt(hdr.time) / GADGET_VELOCITY_CONVERSION;
#ifdef EXACT_OUTPUT_REDSHIFTS
	Real phip = 0.;
	Real gradphi[3] = {0., 0., 0.};
	double ref_dist[3];
	LATfield2::Site xField;
#endif
	
	filename.copy(fname, filename.length());
	fname[filename.length()] = '\0';
	
	LATfield2::Site xPart(this->lat_part_);
	typename std::list<part>::iterator it;
	
	if (hdr.num_files != 1 && hdr.num_files != parallel.grid_size()[1])
	{
		COUT << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": number of Gadget2 files does not match the number of processes in dim-1!" << endl;
		return;
	}
	
	posdata = (float *) malloc(3 * sizeof(float) * PCLBUFFER);
	veldata = (float *) malloc(3 * sizeof(float) * PCLBUFFER);

#if GADGET_ID_BYTES == 8
	IDs = malloc(sizeof(int64_t) * PCLBUFFER);
#else
	IDs = malloc(sizeof(int32_t) * PCLBUFFER);
#endif
	
	npart = 0;
	for(xPart.first(); xPart.test(); xPart.next())
	{
		if(this->field_part_(xPart).size!=0)
		{
			for (it=(this->field_part_)(xPart).parts.begin(); it != (this->field_part_)(xPart).parts.end(); ++it)
			{
				if ((*it).ID % tracer_factor == 0)
					npart++;
			}
		}
	}

	if (hdr.num_files == 1)
	{	
		if (parallel.rank() == 0)
		{
			parallel.send<long>(npart, 1);
			parallel.receive<long>(count, parallel.size()-1);
			if (count != hdr.npart[1]) cout << " error: number of particles in saveGadget2 does not match request!" << endl;
			count = 0;
		}
		else
		{
			parallel.receive<long>(count, parallel.rank()-1);
			npart += count;
			parallel.send<long>(npart, (parallel.rank()+1)%parallel.size());
		}
	
		MPI_File_open(parallel.lat_world_comm(), fname, MPI_MODE_WRONLY | MPI_MODE_CREATE,  MPI_INFO_NULL, &outfile);
	}
	else
	{
		if (parallel.grid_rank()[0] == 0)
		{
			parallel.send_dim0<long>(npart, 1);
			parallel.receive_dim0<long>(count, parallel.grid_size()[0]-1);
			hdr.npart[1] = (uint32_t) count;
			count = 0;
		}
		else
		{
			parallel.receive_dim0<long>(count, parallel.grid_rank()[0]-1);
			npart += count;
			parallel.send_dim0<long>(npart, (parallel.grid_rank()[0]+1)%parallel.grid_size()[0]);
		}

		parallel.broadcast_dim0<uint32_t>(hdr.npart[1], 0);
		
		sprintf(fname+filename.length(), ".%d", parallel.grid_rank()[1]);

		MPI_File_open(parallel.dim0_comm()[parallel.grid_rank()[1]], fname, MPI_MODE_WRONLY | MPI_MODE_CREATE,  MPI_INFO_NULL, &outfile);
	}
	
	offset_pos = (MPI_Offset) hdr.npart[1];
	offset_pos *= (MPI_Offset) (6 * sizeof(float) + ((GADGET_ID_BYTES == 8) ? sizeof(int64_t) : sizeof(int32_t)));
	offset_pos += (MPI_Offset) (8 * sizeof(uint32_t) + sizeof(hdr));
	MPI_File_set_size(outfile, offset_pos);
	
	offset_pos = (MPI_Offset) (3 * sizeof(uint32_t) + sizeof(hdr)) + ((MPI_Offset) count) * ((MPI_Offset) (3 * sizeof(float)));
	offset_vel = offset_pos + (MPI_Offset) (2 * sizeof(uint32_t)) + ((MPI_Offset) hdr.npart[1]) * ((MPI_Offset) (3 * sizeof(float)));
	offset_ID = offset_vel + (MPI_Offset) (2 * sizeof(uint32_t)) + ((MPI_Offset) hdr.npart[1] - (MPI_Offset) count) * ((MPI_Offset) (3 * sizeof(float))) + ((MPI_Offset) count) * ((MPI_Offset) ((GADGET_ID_BYTES == 8) ? sizeof(int64_t) : sizeof(int32_t)));
	
	if ((hdr.num_files == 1 && parallel.rank() == 0) || (hdr.num_files > 1 && parallel.grid_rank()[0] == 0))
	{
		blocksize = sizeof(hdr);		
		MPI_File_write_at(outfile, 0, &blocksize, 1, MPI_UNSIGNED, &status);
		MPI_File_write_at(outfile, sizeof(uint32_t), &hdr, sizeof(hdr), MPI_BYTE, &status);
		MPI_File_write_at(outfile, sizeof(hdr) + sizeof(uint32_t), &blocksize, 1, MPI_UNSIGNED, &status);
		blocksize = 3 * sizeof(float) * hdr.npart[1];
		MPI_File_write_at(outfile, sizeof(hdr) + 2*sizeof(uint32_t), &blocksize, 1, MPI_UNSIGNED, &status);
		MPI_File_write_at(outfile, offset_vel - 2*sizeof(uint32_t), &blocksize, 1, MPI_UNSIGNED, &status);
		MPI_File_write_at(outfile, offset_vel - sizeof(uint32_t), &blocksize, 1, MPI_UNSIGNED, &status);
		MPI_File_write_at(outfile, offset_ID - 2*sizeof(uint32_t), &blocksize, 1, MPI_UNSIGNED, &status);
		blocksize = ((GADGET_ID_BYTES == 8) ? sizeof(int64_t) : sizeof(int32_t)) * hdr.npart[1];
		MPI_File_write_at(outfile, offset_ID - sizeof(uint32_t), &blocksize, 1, MPI_UNSIGNED, &status);
		MPI_File_write_at(outfile, offset_ID + blocksize, &blocksize, 1, MPI_UNSIGNED, &status);
	}
	
	count = 0;
#ifdef EXACT_OUTPUT_REDSHIFTS
	if (phi != NULL)
		xField.initialize(phi->lattice());
	else
		xField.initialize(this->lat_part_);
	xField.first();
#endif
	for(xPart.first(); xPart.test(); xPart.next())
	{
		if(this->field_part_(xPart).size!=0)
		{
			for (it=(this->field_part_)(xPart).parts.begin(); it != (this->field_part_)(xPart).parts.end(); ++it)
			{
				if ((*it).ID % tracer_factor == 0)
				{
#ifdef EXACT_OUTPUT_REDSHIFTS
					if (phi != NULL)
					{
						for (i = 0; i < 3; i++)
							ref_dist[i] = modf((*it).pos[i] / this->lat_resolution_, &phip);
							
						phip = (*phi)(xField) * (1.-ref_dist[0]) * (1.-ref_dist[1]) * (1.-ref_dist[2]);
						phip += (*phi)(xField+0) * ref_dist[0] * (1.-ref_dist[1]) * (1.-ref_dist[2]);
						phip += (*phi)(xField+1) * (1.-ref_dist[0]) * ref_dist[1] * (1.-ref_dist[2]);
						phip += (*phi)(xField+0+1) * ref_dist[0] * ref_dist[1] * (1.-ref_dist[2]);
						phip += (*phi)(xField+2) * (1.-ref_dist[0]) * (1.-ref_dist[1]) * ref_dist[2];
						phip += (*phi)(xField+0+2) * ref_dist[0] * (1.-ref_dist[1]) * ref_dist[2];
						phip += (*phi)(xField+1+2) * (1.-ref_dist[0]) * ref_dist[1] * ref_dist[2];
						phip += (*phi)(xField+0+1+2) * ref_dist[0] * ref_dist[1] * ref_dist[2];
						
						gradphi[0] = (1.-ref_dist[1]) * (1.-ref_dist[2]) * ((*phi)(xField+0) - (*phi)(xField));
						gradphi[1] = (1.-ref_dist[0]) * (1.-ref_dist[2]) * ((*phi)(xField+1) - (*phi)(xField));
						gradphi[2] = (1.-ref_dist[0]) * (1.-ref_dist[1]) * ((*phi)(xField+2) - (*phi)(xField));
						gradphi[0] += ref_dist[1] * (1.-ref_dist[2]) * ((*phi)(xField+1+0) - (*phi)(xField+1));
						gradphi[1] += ref_dist[0] * (1.-ref_dist[2]) * ((*phi)(xField+1+0) - (*phi)(xField+0));
						gradphi[2] += ref_dist[0] * (1.-ref_dist[1]) * ((*phi)(xField+2+0) - (*phi)(xField+0));
						gradphi[0] += (1.-ref_dist[1]) * ref_dist[2] * ((*phi)(xField+2+0) - (*phi)(xField+2));
						gradphi[1] += (1.-ref_dist[0]) * ref_dist[2] * ((*phi)(xField+2+1) - (*phi)(xField+2));
						gradphi[2] += (1.-ref_dist[0]) * ref_dist[1] * ((*phi)(xField+2+1) - (*phi)(xField+1));
						gradphi[0] += ref_dist[1] * ref_dist[2] * ((*phi)(xField+2+1+0) - (*phi)(xField+2+1));
						gradphi[1] += ref_dist[0] * ref_dist[2] * ((*phi)(xField+2+1+0) - (*phi)(xField+2+0));
						gradphi[2] += ref_dist[0] * ref_dist[1] * ((*phi)(xField+2+1+0) - (*phi)(xField+1+0));
					}
					
					ref_dist[0] = (*it).vel[0]*(*it).vel[0] + (*it).vel[1]*(*it).vel[1] + (*it).vel[2]*(*it).vel[2];
					ref_dist[1] = ref_dist[0] + hdr.time * hdr.time;
					ref_dist[2] = sqrt(ref_dist[1]);
					ref_dist[0] += ref_dist[1];
					ref_dist[1] = 1. + (4. - (ref_dist[0] / ref_dist[1])) * phip;
					
					for (i = 0; i < 3; i++)
						posdata[3*count+i] = modf((*it).pos[i] + dtau_pos * (*it).vel[i] * ref_dist[1] / ref_dist[2], &phip) * hdr.BoxSize;
					
					for (i = 0; i < 3; i++)
						veldata[3*count+i] = ((*it).vel[i] - dtau_vel * ref_dist[0] * gradphi[i] / this->lat_resolution_ / ref_dist[2]) * rescale_vel / hdr.time;
#else						
					for (i = 0; i < 3; i++)
						posdata[3*count+i] = (*it).pos[i] * hdr.BoxSize;
					
					for (i = 0; i < 3; i++)
						veldata[3*count+i] = (*it).vel[i] * rescale_vel / hdr.time;
#endif
					
#if GADGET_ID_BYTES == 8
					*((int64_t *) IDs + count) = (int64_t) (*it).ID;
#else	
					*((int32_t *) IDs + count) = (int32_t) (*it).ID;
#endif
					
					count++;
						
					if (count == PCLBUFFER)
					{
						MPI_File_write_at(outfile, offset_pos, posdata, 3 * count, MPI_FLOAT, &status);
						offset_pos += 3 * PCLBUFFER * sizeof(float);
						MPI_File_write_at(outfile, offset_vel, veldata, 3 * count, MPI_FLOAT, &status);
						offset_vel += 3 * PCLBUFFER * sizeof(float);
						count *= (GADGET_ID_BYTES == 8) ? sizeof(int64_t) : sizeof(int32_t);
						MPI_File_write_at(outfile, offset_ID, IDs, count, MPI_BYTE, &status);
						offset_ID += count;
						count = 0;
					}
				}
			}
		}
#ifdef EXACT_OUTPUT_REDSHIFTS
		xField.next();
#endif
	}
	
	MPI_File_write_at_all(outfile, offset_pos, posdata, 3 * count, MPI_FLOAT, &status);
	MPI_File_write_at_all(outfile, offset_vel, veldata, 3 * count, MPI_FLOAT, &status);
	count *= (GADGET_ID_BYTES == 8) ? sizeof(int64_t) : sizeof(int32_t);
	MPI_File_write_at_all(outfile, offset_ID, IDs, count, MPI_BYTE, &status);
	
	MPI_File_close(&outfile);
	
	free(posdata);
	free(veldata);
	free(IDs);
}


template <typename part, typename part_info, typename part_dataType>
void Particles_gevolution<part,part_info,part_dataType>::saveGadget2(string filename, gadget2_header & hdr, lightcone_geometry & lightcone, double dist, double dtau, double dtau_old, double dadtau, double vertex[MAX_INTERSECTS][3], const int vertexcount, set<long> & IDbacklog, set<long> & IDprelog, Field<Real> * phi, const int tracer_factor)
{
	float * posdata;
	float * veldata;
	void * IDs;
	MPI_File outfile;
	long count, npart;
	MPI_Offset offset_pos, offset_vel, offset_ID;
	MPI_Status status;
	uint32_t blocksize;
	uint32_t i;
	char fname[filename.length()+1];
	double rescale_vel = 1. / GADGET_VELOCITY_CONVERSION;
	double inner = dist - 0.5 * dtau;
	double outer = dist + (0.5 + LIGHTCONE_IDCHECK_ZONE) * dtau_old;
	double d, v2, e2;
	double ref_dist[3];
	Real gradphi[3];
	
	filename.copy(fname, filename.length());
	fname[filename.length()] = '\0';
	
	LATfield2::Site xPart(this->lat_part_);
	LATfield2::Site xField(phi->lattice());
	typename std::list<part>::iterator it;
	
	if (hdr.num_files != 1)
	{
		COUT << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": writing multiple Gadget2 files not currently supported!" << endl;
		return;
	}
	
	posdata = (float *) malloc(3 * sizeof(float) * PCLBUFFER);
	veldata = (float *) malloc(3 * sizeof(float) * PCLBUFFER);

#if GADGET_ID_BYTES == 8
	IDs = malloc(sizeof(int64_t) * PCLBUFFER);
#else
	IDs = malloc(sizeof(int32_t) * PCLBUFFER);
#endif
	
	npart = 0;
	if (vertexcount > 0)
	{
		for(xPart.first(), xField.first(); xPart.test(); xPart.next(), xField.next())
		{
			if(this->field_part_(xPart).size!=0)
			{
				for (it=(this->field_part_)(xPart).parts.begin(); it != (this->field_part_)(xPart).parts.end(); ++it)
				{
					if ((*it).ID % tracer_factor == 0)
					{
						for (i = 0; i < (uint32_t) vertexcount; i++)
						{
							d = sqrt(((*it).pos[0]-vertex[i][0])*((*it).pos[0]-vertex[i][0]) + ((*it).pos[1]-vertex[i][1])*((*it).pos[1]-vertex[i][1]) + ((*it).pos[2]-vertex[i][2])*((*it).pos[2]-vertex[i][2]));

							if (d < inner || d >= outer) continue;

							if (lightcone.opening == -1. || (((*it).pos[0]-vertex[i][0])*lightcone.direction[0] + ((*it).pos[1]-vertex[i][1])*lightcone.direction[1] + ((*it).pos[2]-vertex[i][2])*lightcone.direction[2]) / d > lightcone.opening)
							{
								if (outer - d > LIGHTCONE_IDCHECK_ZONE * dtau_old || IDbacklog.find((*it).ID) == IDbacklog.end())
								{
									if (d - inner < 2. * LIGHTCONE_IDCHECK_ZONE * dtau)
										IDprelog.insert((*it).ID);

									for (int j = 0; j < 3; j++)
										ref_dist[j] = modf((*it).pos[j] / this->lat_resolution_, &v2);
									
									v2 = (*it).vel[0] * (*it).vel[0] + (*it).vel[1] * (*it).vel[1] + (*it).vel[2] * (*it).vel[2];
									e2 = v2 + hdr.time * (hdr.time + (dist - d - 0.5 * dtau_old) * dadtau);
	
									gradphi[0] = (1.-ref_dist[1]) * (1.-ref_dist[2]) * ((*phi)(xField+0) - (*phi)(xField));
									gradphi[1] = (1.-ref_dist[0]) * (1.-ref_dist[2]) * ((*phi)(xField+1) - (*phi)(xField));
									gradphi[2] = (1.-ref_dist[0]) * (1.-ref_dist[1]) * ((*phi)(xField+2) - (*phi)(xField));
									gradphi[0] += ref_dist[1] * (1.-ref_dist[2]) * ((*phi)(xField+1+0) - (*phi)(xField+1));
									gradphi[1] += ref_dist[0] * (1.-ref_dist[2]) * ((*phi)(xField+1+0) - (*phi)(xField+0));
									gradphi[2] += ref_dist[0] * (1.-ref_dist[1]) * ((*phi)(xField+2+0) - (*phi)(xField+0));
									gradphi[0] += (1.-ref_dist[1]) * ref_dist[2] * ((*phi)(xField+2+0) - (*phi)(xField+2));
									gradphi[1] += (1.-ref_dist[0]) * ref_dist[2] * ((*phi)(xField+2+1) - (*phi)(xField+2));
									gradphi[2] += (1.-ref_dist[0]) * ref_dist[1] * ((*phi)(xField+2+1) - (*phi)(xField+1));
									gradphi[0] += ref_dist[1] * ref_dist[2] * ((*phi)(xField+2+1+0) - (*phi)(xField+2+1));
									gradphi[1] += ref_dist[0] * ref_dist[2] * ((*phi)(xField+2+1+0) - (*phi)(xField+2+0));
									gradphi[2] += ref_dist[0] * ref_dist[1] * ((*phi)(xField+2+1+0) - (*phi)(xField+1+0));

									gradphi[0] *= (v2 + e2) / e2 / this->lat_resolution_;
									gradphi[1] *= (v2 + e2) / e2 / this->lat_resolution_;
									gradphi[2] *= (v2 + e2) / e2 / this->lat_resolution_;
						
									e2 = sqrt(e2);

									for (uint32_t j = 0; j < 3; j++)
										veldata[3*(npart%PCLBUFFER)+j] = ((*it).vel[j] - (dist - d + 0.5 * dtau_old) * e2 * gradphi[j]) * rescale_vel / (hdr.time + (dist - d) * dadtau);
										
									if (d >= dist)
									{
										e2 = sqrt(v2 + hdr.time * (hdr.time - dtau_old * dadtau));
										
										for (uint32_t j = 0; j < 3; j++)
											posdata[3*(npart%PCLBUFFER)+j] = ((*it).pos[j] - vertex[i][j] + lightcone.vertex[j] + (dist - d) * (*it).vel[j] / e2) * hdr.BoxSize;
									}
									else
									{
										e2 = sqrt(v2 + hdr.time * (hdr.time + dtau * dadtau));
										v2 = sqrt(v2 + hdr.time * hdr.time);
										
										for (uint32_t j = 0; j < 3; j++)
											posdata[3*(npart%PCLBUFFER)+j] = ((*it).pos[j] - vertex[i][j] + lightcone.vertex[j] + (dist - d) * ((*it).vel[j] - dtau * v2 * gradphi[j]) / e2) * hdr.BoxSize;
									}
					
#if GADGET_ID_BYTES == 8
									*((int64_t *) IDs + (npart%PCLBUFFER)) = (int64_t) (*it).ID;
#else	
									*((int32_t *) IDs + (npart%PCLBUFFER)) = (int32_t) (*it).ID;
#endif

									npart++;
								}

								break;
							}
						}
					}
				}
			}
		}
	}
	
	if (parallel.rank() == 0)
	{
		parallel.send<long>(npart, 1);
		parallel.receive<long>(count, parallel.size()-1);
		hdr.npart[1] = (uint32_t) (count % (1ll << 32));
		hdr.npartTotal[1] = (uint32_t) (count % (1ll << 32));
		hdr.npartTotalHW[1] = (uint32_t) (count / (1ll << 32));
		count = 0;
	}
	else
	{
		parallel.receive<long>(count, parallel.rank()-1);
		count += npart;
		parallel.send<long>(count, (parallel.rank()+1)%parallel.size());
		count -= npart;
	}

	parallel.broadcast<uint32_t>(hdr.npartTotal[1], 0);
	parallel.broadcast<uint32_t>(hdr.npartTotalHW[1], 0);

	if (hdr.npartTotal[1] + ((int64_t) hdr.npartTotalHW[1] << 32) > 0)
	{
		MPI_File_open(parallel.lat_world_comm(), fname, MPI_MODE_WRONLY | MPI_MODE_CREATE,  MPI_INFO_NULL, &outfile);
	
		offset_pos = (MPI_Offset) ((int64_t) hdr.npartTotal[1] + ((int64_t) hdr.npartTotalHW[1] << 32));
		offset_pos *= (MPI_Offset) (6 * sizeof(float) + ((GADGET_ID_BYTES == 8) ? sizeof(int64_t) : sizeof(int32_t)));
		offset_pos += (MPI_Offset) (8 * sizeof(uint32_t) + sizeof(hdr));
		MPI_File_set_size(outfile, offset_pos);
	
		offset_pos = (MPI_Offset) (3 * sizeof(uint32_t) + sizeof(hdr)) + ((MPI_Offset) count) * ((MPI_Offset) (3 * sizeof(float)));
		offset_vel = offset_pos + (MPI_Offset) (2 * sizeof(uint32_t)) + ((MPI_Offset) ((int64_t) hdr.npartTotal[1] + ((int64_t) hdr.npartTotalHW[1] << 32))) * ((MPI_Offset) (3 * sizeof(float)));
		offset_ID = offset_vel + (MPI_Offset) (2 * sizeof(uint32_t)) + ((MPI_Offset) ((int64_t) hdr.npartTotal[1] + ((int64_t) hdr.npartTotalHW[1] << 32)) - (MPI_Offset) count) * ((MPI_Offset) (3 * sizeof(float))) + ((MPI_Offset) count) * ((MPI_Offset) ((GADGET_ID_BYTES == 8) ? sizeof(int64_t) : sizeof(int32_t)));
	
		if (parallel.rank() == 0)
		{
			blocksize = sizeof(hdr);		
			MPI_File_write_at(outfile, 0, &blocksize, 1, MPI_UNSIGNED, &status);
			MPI_File_write_at(outfile, sizeof(uint32_t), &hdr, sizeof(hdr), MPI_BYTE, &status);
			MPI_File_write_at(outfile, sizeof(hdr) + sizeof(uint32_t), &blocksize, 1, MPI_UNSIGNED, &status);
			blocksize = 3 * sizeof(float) * hdr.npart[1];
			MPI_File_write_at(outfile, sizeof(hdr) + 2*sizeof(uint32_t), &blocksize, 1, MPI_UNSIGNED, &status);
			MPI_File_write_at(outfile, offset_vel - 2*sizeof(uint32_t), &blocksize, 1, MPI_UNSIGNED, &status);
			MPI_File_write_at(outfile, offset_vel - sizeof(uint32_t), &blocksize, 1, MPI_UNSIGNED, &status);
			MPI_File_write_at(outfile, offset_ID - 2*sizeof(uint32_t), &blocksize, 1, MPI_UNSIGNED, &status);
			blocksize = ((GADGET_ID_BYTES == 8) ? sizeof(int64_t) : sizeof(int32_t)) * hdr.npart[1];
			MPI_File_write_at(outfile, offset_ID - sizeof(uint32_t), &blocksize, 1, MPI_UNSIGNED, &status);
			MPI_File_write_at(outfile, offset_ID + blocksize, &blocksize, 1, MPI_UNSIGNED, &status);
		}
		
		count = (npart < PCLBUFFER) ? npart : PCLBUFFER;
		npart -= count;
		MPI_File_write_at_all(outfile, offset_pos, posdata, 3 * count, MPI_FLOAT, &status);
		offset_pos += 3 * count * sizeof(float);
		MPI_File_write_at_all(outfile, offset_vel, veldata, 3 * count, MPI_FLOAT, &status);
		offset_vel += 3 * count * sizeof(float);
		count *= (GADGET_ID_BYTES == 8) ? sizeof(int64_t) : sizeof(int32_t);
		MPI_File_write_at_all(outfile, offset_ID, IDs, count, MPI_BYTE, &status);
		offset_ID += count;
		count = 0;

		if (npart > 0)
		{
			for(xPart.first(); xPart.test() && npart > 0; xPart.next())
			{
				if(this->field_part_(xPart).size!=0)
				{
					for (it=(this->field_part_)(xPart).parts.begin(); it != (this->field_part_)(xPart).parts.end(); ++it)
					{
						if ((*it).ID % tracer_factor == 0)
						{
							for (i = 0; i < (uint32_t) vertexcount; i++)
							{
								d = sqrt(((*it).pos[0]-vertex[i][0])*((*it).pos[0]-vertex[i][0]) + ((*it).pos[1]-vertex[i][1])*((*it).pos[1]-vertex[i][1]) + ((*it).pos[2]-vertex[i][2])*((*it).pos[2]-vertex[i][2]));

								if (d < inner || d >= outer) continue;

								if (lightcone.opening == -1. || (((*it).pos[0]-vertex[i][0])*lightcone.direction[0] + ((*it).pos[1]-vertex[i][1])*lightcone.direction[1] + ((*it).pos[2]-vertex[i][2])*lightcone.direction[2]) / d > lightcone.opening)
								{
									if (outer - d > LIGHTCONE_IDCHECK_ZONE * dtau_old || IDbacklog.find((*it).ID) == IDbacklog.end())
									{
										for (int j = 0; j < 3; j++)
											ref_dist[j] = modf((*it).pos[j] / this->lat_resolution_, &v2);
									
										v2 = (*it).vel[0] * (*it).vel[0] + (*it).vel[1] * (*it).vel[1] + (*it).vel[2] * (*it).vel[2];
										e2 = v2 + hdr.time * (hdr.time + (dist - d - 0.5 * dtau_old) * dadtau);
	
										gradphi[0] = (1.-ref_dist[1]) * (1.-ref_dist[2]) * ((*phi)(xField+0) - (*phi)(xField));
										gradphi[1] = (1.-ref_dist[0]) * (1.-ref_dist[2]) * ((*phi)(xField+1) - (*phi)(xField));
										gradphi[2] = (1.-ref_dist[0]) * (1.-ref_dist[1]) * ((*phi)(xField+2) - (*phi)(xField));
										gradphi[0] += ref_dist[1] * (1.-ref_dist[2]) * ((*phi)(xField+1+0) - (*phi)(xField+1));
										gradphi[1] += ref_dist[0] * (1.-ref_dist[2]) * ((*phi)(xField+1+0) - (*phi)(xField+0));
										gradphi[2] += ref_dist[0] * (1.-ref_dist[1]) * ((*phi)(xField+2+0) - (*phi)(xField+0));
										gradphi[0] += (1.-ref_dist[1]) * ref_dist[2] * ((*phi)(xField+2+0) - (*phi)(xField+2));
										gradphi[1] += (1.-ref_dist[0]) * ref_dist[2] * ((*phi)(xField+2+1) - (*phi)(xField+2));
										gradphi[2] += (1.-ref_dist[0]) * ref_dist[1] * ((*phi)(xField+2+1) - (*phi)(xField+1));
										gradphi[0] += ref_dist[1] * ref_dist[2] * ((*phi)(xField+2+1+0) - (*phi)(xField+2+1));
										gradphi[1] += ref_dist[0] * ref_dist[2] * ((*phi)(xField+2+1+0) - (*phi)(xField+2+0));
										gradphi[2] += ref_dist[0] * ref_dist[1] * ((*phi)(xField+2+1+0) - (*phi)(xField+1+0));

										gradphi[0] *= (v2 + e2) / e2 / this->lat_resolution_;
										gradphi[1] *= (v2 + e2) / e2 / this->lat_resolution_;
										gradphi[2] *= (v2 + e2) / e2 / this->lat_resolution_;
						
										e2 = sqrt(e2);

										for (uint32_t j = 0; j < 3; j++)
											veldata[3*(npart%PCLBUFFER)+j] = ((*it).vel[j] - (dist - d + 0.5 * dtau_old) * e2 * gradphi[j]) * rescale_vel / (hdr.time + (dist - d) * dadtau);
											
										if (d >= dist)
										{
											e2 = sqrt(v2 + hdr.time * (hdr.time - dtau_old * dadtau));
										
											for (uint32_t j = 0; j < 3; j++)
												posdata[3*(npart%PCLBUFFER)+j] = ((*it).pos[j] - vertex[i][j] + lightcone.vertex[j] + (dist - d) * (*it).vel[j] / e2) * hdr.BoxSize;
										}
										else
										{
											e2 = sqrt(v2 + hdr.time * (hdr.time + dtau * dadtau));
											v2 = sqrt(v2 + hdr.time * hdr.time);
										
											for (uint32_t j = 0; j < 3; j++)
												posdata[3*(npart%PCLBUFFER)+j] = ((*it).pos[j] - vertex[i][j] + lightcone.vertex[j] + (dist - d) * ((*it).vel[j] - dtau * v2 * gradphi[j]) / e2) * hdr.BoxSize;
										}
						
#if GADGET_ID_BYTES == 8
										*((int64_t *) IDs + count) = (int64_t) (*it).ID;
#else	
										*((int32_t *) IDs + count) = (int32_t) (*it).ID;
#endif
	
										npart--;
										count++;
									}
									break;
								}
							}
							
							if (count == PCLBUFFER)
							{
								MPI_File_write_at(outfile, offset_pos, posdata, 3 * count, MPI_FLOAT, &status);
								offset_pos += 3 * PCLBUFFER * sizeof(float);
								MPI_File_write_at(outfile, offset_vel, veldata, 3 * count, MPI_FLOAT, &status);
								offset_vel += 3 * PCLBUFFER * sizeof(float);
								count *= (GADGET_ID_BYTES == 8) ? sizeof(int64_t) : sizeof(int32_t);
								MPI_File_write_at(outfile, offset_ID, IDs, count, MPI_BYTE, &status);
								offset_ID += count;
								count = 0;
							}

							if (npart <= 0) break;
						}
					}
				}
			}

			if (count > 0)
			{
					MPI_File_write_at(outfile, offset_pos, posdata, 3 * count, MPI_FLOAT, &status);
					MPI_File_write_at(outfile, offset_vel, veldata, 3 * count, MPI_FLOAT, &status);
					count *= (GADGET_ID_BYTES == 8) ? sizeof(int64_t) : sizeof(int32_t);
					MPI_File_write_at(outfile, offset_ID, IDs, count, MPI_BYTE, &status);
			}
		}	
	
		MPI_File_close(&outfile);
	}
	
	free(posdata);
	free(veldata);
	free(IDs);
}


template <typename part, typename part_info, typename part_dataType>
void Particles_gevolution<part,part_info,part_dataType>::loadGadget2(string filename, gadget2_header & hdr)
{
	float * posdata;
	float * veldata;
	void * IDs;
	part pcl;
	MPI_File infile;
	uint32_t i, count, npart = 0;
	MPI_Offset offset_pos, offset_vel, offset_ID;
	MPI_Status status;
	uint32_t blocksize;
	char fname[filename.length()+1];
	double rescale_vel = 1. / GADGET_VELOCITY_CONVERSION;
	
	filename.copy(fname, filename.length());
	fname[filename.length()] = '\0';
	
	posdata = (float *) malloc(3 * sizeof(float) * PCLBUFFER);
	veldata = (float *) malloc(3 * sizeof(float) * PCLBUFFER);

#if GADGET_ID_BYTES == 8
	IDs = malloc(sizeof(int64_t) * PCLBUFFER);
#else
	IDs = malloc(sizeof(int32_t) * PCLBUFFER);
#endif

	MPI_File_open(parallel.lat_world_comm(), fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &infile);

	MPI_File_read_all(infile, &blocksize, 1, MPI_UNSIGNED, &status);

	if (blocksize != sizeof(hdr))
	{
		COUT << COLORTEXT_RED << " error" << COLORTEXT_RESET << ": file type not recognized when reading Gadget2 file!" << endl;
		return;
	}

	MPI_File_read_all(infile, &hdr, sizeof(hdr), MPI_BYTE, &status);

	rescale_vel /= sqrt(hdr.time);
	offset_pos = (MPI_Offset) sizeof(hdr) + (MPI_Offset) (3 * sizeof(uint32_t));
	offset_vel = offset_pos + ((MPI_Offset) hdr.npart[1]) * ((MPI_Offset) (3 * sizeof(float))) + (MPI_Offset) (2 * sizeof(uint32_t));
	offset_ID = offset_vel + offset_vel - offset_pos;

	MPI_File_seek(infile, offset_pos, MPI_SEEK_SET);
	while (npart < hdr.npart[1])
	{
		count = (hdr.npart[1] - npart > PCLBUFFER) ? PCLBUFFER : (hdr.npart[1] - npart);

		MPI_File_read_all(infile, posdata, 3 * count, MPI_FLOAT, &status);
		offset_pos += (MPI_Offset) (3 * count * sizeof(float));
		MPI_File_seek(infile, offset_vel, MPI_SEEK_SET);
		MPI_File_read_all(infile, veldata, 3 * count, MPI_FLOAT, &status);
		offset_vel += (MPI_Offset) (3 * count * sizeof(float));
		MPI_File_seek(infile, offset_ID, MPI_SEEK_SET);
#if GADGET_ID_BYTES == 8
		MPI_File_read_all(infile, IDs, count * sizeof(int64_t), MPI_BYTE, &status);
		offset_ID += (MPI_Offset) (count * sizeof(int64_t));
#else
		MPI_File_read_all(infile, IDs, count * sizeof(int32_t), MPI_BYTE, &status);
		offset_ID += (MPI_Offset) (count * sizeof(int32_t));
#endif
		MPI_File_seek(infile, offset_pos, MPI_SEEK_SET);

		for (i = 0; i < 3 * count; i++)
		{
			posdata[i] /= hdr.BoxSize;
			if (posdata[i] >= 1.) posdata[i] -= 1.;
			veldata[i] *= hdr.time / rescale_vel;
		}

		for (i = 0; i < count; i++)
		{
#if GADGET_ID_BYTES == 8
			pcl.ID = *((int64_t *) IDs + i);
#else
			pcl.ID = *((int32_t *) IDs + i);
#endif
			pcl.pos[0] = posdata[3*i];
			pcl.pos[1] = posdata[3*i+1];
			pcl.pos[2] = posdata[3*i+2];
			pcl.vel[0] = veldata[3*i];
			pcl.vel[1] = veldata[3*i+1];
			pcl.vel[2] = veldata[3*i+2];
			this->addParticle_global(pcl);
		}
		
		npart += count;
	}

	MPI_File_close(&infile);
	
	free(posdata);
	free(veldata);
	free(IDs);
}

#endif

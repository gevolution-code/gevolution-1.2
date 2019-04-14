//////////////////////////
// metadata.hpp
//////////////////////////
// 
// Constants and metadata structures
//
// Author: Julian Adamek (Université de Genève & Observatoire de Paris & Queen Mary University of London)
//
// Last modified: April 2019
//
//////////////////////////

#ifndef METADATA_HEADER
#define METADATA_HEADER

#define GEVOLUTION_VERSION 1.2

#ifndef MAX_OUTPUTS
#define MAX_OUTPUTS 32
#endif

#ifndef PARAM_MAX_LENGTH
#define PARAM_MAX_LENGTH 256
#endif

#ifndef PARAM_MAX_LINESIZE
#define PARAM_MAX_LINESIZE 1024
#endif

#ifndef MAX_INTERSECTS
#define MAX_INTERSECTS 12
#endif

#ifndef LIGHTCONE_IDCHECK_ZONE
#define LIGHTCONE_IDCHECK_ZONE 0.05
#endif

#define LIGHTCONE_PHI_OFFSET 0
#define LIGHTCONE_CHI_OFFSET 1
#define LIGHTCONE_B_OFFSET   2
#define LIGHTCONE_HIJ_OFFSET 5
#define LIGHTCONE_MAX_FIELDS 10

#ifndef MAX_PCL_SPECIES
#define MAX_PCL_SPECIES 6
#endif

#ifndef CYCLE_INFO_INTERVAL
#define CYCLE_INFO_INTERVAL 10
#endif

#define MASK_PHI    1
#define MASK_CHI    2
#define MASK_POT    4
#define MASK_B      8
#define MASK_T00    16
#define MASK_TIJ    32
#define MASK_RBARE  64
#define MASK_HIJ    128
#define MASK_P      256
#define MASK_GADGET 512
#define MASK_PCLS   1024
#define MASK_XSPEC  2048
#define MASK_DELTA  4096
#define MASK_DBARE  8192
#define MASK_MULTI  16384
#define MASK_VEL    32768

#define ICFLAG_CORRECT_DISPLACEMENT 1
#define ICFLAG_KSPHERE              2

// Identifiers for IC generator modules
#define ICGEN_BASIC                 0
#define ICGEN_READ_FROM_DISK        1
#ifdef ICGEN_PREVOLUTION
#undef ICGEN_PREVOLUTION
#define ICGEN_PREVOLUTION           2
#endif
#ifdef ICGEN_SONG
#undef ICGEN_SONG
#define ICGEN_SONG                  3
#endif
#ifdef ICGEN_FALCONIC
#undef ICGEN_FALCONIC
#define ICGEN_FALCONIC              4
#endif

#define VECTOR_PARABOLIC            0
#define VECTOR_ELLIPTIC             1

// Physical constants
#define C_PLANCK_LAW      4.48147e-7    // omega_g / (T_cmb [K])^4
#define C_BOLTZMANN_CST   8.61733e-5    // Boltzmann constant [eV/K]
#define C_SPEED_OF_LIGHT  2997.92458    // speed of light [100 km/s]
#define C_RHO_CRIT        2.77459457e11 // critical density [M_sun h^2 / Mpc^3]
#define C_FD_NORM         1.80308535    // Integral[q*q/(exp(q)+1), 0, infinity]

// default physical parameters (used in parser.hpp)
#define P_HUBBLE          0.67556       // default value for h
#define P_T_NCDM          0.71611       // default value for T_ncdm
#define P_NCDM_MASS_OMEGA 93.14         // m_ncdm / omega_ncdm [eV]
#define P_N_UR            3.046         // default value for N_ur (= N_eff)
#define P_SPECTRAL_AMP    2.215e-9      // default value for A_s
#define P_SPECTRAL_INDEX  0.9619        // default value for n_s
#define P_PIVOT_SCALE     0.05          // default pivot scale [Mpc^-1]

#ifndef GADGET_LENGTH_CONVERSION
#define GADGET_LENGTH_CONVERSION 0.001  // Gadget length unit in Mpc / h
#endif
#ifndef GADGET_MASS_CONVERSION
#define GADGET_MASS_CONVERSION 1.0e10   // Gadget mass unit in M_sun / h
#endif
#ifndef GADGET_VELOCITY_CONVERSION
#define GADGET_VELOCITY_CONVERSION 3.335640952e-6  // Gadget velocity unit / speed of light
#endif
#ifndef GADGET_ID_BYTES
#define GADGET_ID_BYTES 8
#endif

#ifdef EXTERNAL_IO
#ifndef NUMBER_OF_IO_FILES
#define NUMBER_OF_IO_FILES 4
#endif
#endif

// color escape sequences for terminal highlighting (enable with -DCOLORTERMINAL)
#ifdef COLORTERMINAL
#define COLORTEXT_WHITE     "\033[37;1m"
#define COLORTEXT_CYAN      "\033[36;1m"
#define COLORTEXT_GREEN     "\033[32;1m"
#define COLORTEXT_RED       "\033[31;1m"
#define COLORTEXT_YELLOW    "\033[33;1m"
#define COLORTEXT_RESET     "\033[0m"
#else
#define COLORTEXT_WHITE     ""
#define COLORTEXT_CYAN      ""
#define COLORTEXT_GREEN     ""
#define COLORTEXT_RED       ""
#define COLORTEXT_YELLOW    ""
#define COLORTEXT_RESET     ""
#endif

// header structure for GADGET-2 files [V. Springel, N. Yoshida, and S.D. White, New Astron. 6 (2001) 79
// and V. Springel, Mon. Not. R. Astron. Soc. 364 (2005) 1105]

#ifndef GADGET2_HEADER
#define GADGET2_HEADER
struct gadget2_header
{
	uint32_t npart[6];
	double mass[6];
	double time;
	double redshift;
	int32_t flag_sfr;
	int32_t flag_feedback;
	uint32_t npartTotal[6];
	int32_t flag_cooling;
	int32_t num_files;
	double BoxSize;
	double Omega0;
	double OmegaLambda;
	double HubbleParam;
	int32_t flag_age;
	int32_t flag_metals;
	uint32_t npartTotalHW[6];
	char fill[256 - 6 * 4 - 6 * 8 - 2 * 8 - 2 * 4 - 6 * 4 - 2 * 4 - 4 * 8 - 2 * 4 - 6 * 4]; /* fills to 256 Bytes */
};
#endif

#ifdef HAVE_HEALPIX
#include "chealpix.h"
#ifndef PIXBUFFER
#define PIXBUFFER 1048576
#endif

struct healpix_header
{
	uint32_t Nside;
	uint32_t Npix;
	uint32_t precision;
	uint32_t Ngrid;
	double direction[3];
	double distance;
	double boxsize;
	uint32_t Nside_ring;
	char fill[256 - 5 * 4 - 5 * 8]; /* fills to 256 Bytes */
};
#endif

struct lightcone_geometry
{
	double vertex[3];
	double z;
	double direction[3];
	double opening;
	double distance[2];
};

struct metadata
{
	int numpts;
	int downgrade_factor;
	long numpcl[MAX_PCL_SPECIES];
	int tracer_factor[MAX_PCL_SPECIES];
	int baryon_flag;
	int gr_flag;
	int vector_flag;
	int radiation_flag;
	int fluid_flag;
	int out_pk;
	int out_snapshot;
	int out_lightcone[MAX_OUTPUTS];
	int num_pk;
	int numbins;
	int num_snapshot;
	int num_lightcone;
	int num_restart;
	int Nside[MAX_OUTPUTS][2];
	double Cf;
	double movelimit;
	double steplimit;
	double boxsize;
	double wallclocklimit;
	double pixelfactor[MAX_OUTPUTS];
	double shellfactor[MAX_OUTPUTS];
	double covering[MAX_OUTPUTS];
	double z_in;
	double z_snapshot[MAX_OUTPUTS];
	double z_pk[MAX_OUTPUTS];
	double z_restart[MAX_OUTPUTS];
	double z_switch_deltarad;
	double z_switch_linearchi;
	double z_switch_deltancdm[MAX_PCL_SPECIES-2];
	double z_switch_Bncdm[MAX_PCL_SPECIES-2];
	lightcone_geometry lightcone[MAX_OUTPUTS];
	char basename_lightcone[PARAM_MAX_LENGTH];
	char basename_snapshot[PARAM_MAX_LENGTH];
	char basename_pk[PARAM_MAX_LENGTH];
	char basename_generic[PARAM_MAX_LENGTH];
	char output_path[PARAM_MAX_LENGTH];
	char restart_path[PARAM_MAX_LENGTH];
	char basename_restart[PARAM_MAX_LENGTH];
};

struct icsettings
{
	int numtile[MAX_PCL_SPECIES];
	int seed;
	int flags;
	int generator;
	int restart_cycle;
	char pclfile[MAX_PCL_SPECIES][PARAM_MAX_LENGTH];
	char pkfile[PARAM_MAX_LENGTH];
	char tkfile[PARAM_MAX_LENGTH];
	char metricfile[3][PARAM_MAX_LENGTH];
	double restart_tau;
	double restart_dtau;
	double restart_version;
	double z_ic;
	double z_relax;
	double Cf;
	double A_s;
	double n_s;
	double k_pivot;
};

struct cosmology
{
	double Omega_cdm;
	double Omega_b;
	double Omega_m;
	double Omega_Lambda;
	double Omega_fld;
	double w0_fld;
	double wa_fld;
	double cs2_fld;
	double Omega_g;
	double Omega_ur;
	double Omega_rad;
	double Omega_ncdm[MAX_PCL_SPECIES-2];
	double h;
	double m_ncdm[MAX_PCL_SPECIES-2];
	double T_ncdm[MAX_PCL_SPECIES-2];
	double deg_ncdm[MAX_PCL_SPECIES-2];
	int num_ncdm;
};

#endif

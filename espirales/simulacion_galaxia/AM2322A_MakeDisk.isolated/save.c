#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "nrsrc/nr.h"
#include "nrsrc/nrutil.h"
#include "prototypes.h"
#include "globvars.h"


#ifdef T3E
typedef short int int4byte;	/* Note: int has 8 Bytes on the T3E ! */
#else
typedef int int4byte;
#endif



// header for gadget basic
/*
  struct io_header_1
  {
  int4byte npart[6];
  double mass[6];
  double time;
  double redshift;
  int4byte flag_sfr;
  int4byte flag_feedback;
  int4byte npartTotal[6];
  int4byte flag_cooling;
  int4byte num_files;
  double BoxSize;
  double Omega0;
  double OmegaLambda;
  double HubbleParam;
  char fill[256 - 6 * 4 - 6 * 8 - 2 * 8 - 2 * 4 - 6 * 4 - 2 * 4 - 4 * 8];	// fills to 256 Bytes
  }
  header1;
*/

// header for P-gadget3
/*
  struct io_header_1
  {
  int npart[6];                 //!< number of particles of each type in this file 
  double mass[6];               //!< mass of particles of each type. If 0, then the masses are explicitly                                                                                                                       //stored in the mass-block of the snapshot file, otherwise they are omitted   
  double time;                  //!< time of snapshot file
  double redshift;              //!< redshift of snapshot file 
  int flag_sfr;                 //!< flags whether the simulation was including star formation 
  int flag_feedback;            //!< flags whether feedback was included (obsolete) 
  unsigned int npartTotal[6];   //!< total number of particles of each type in this snapshot. This can be                                                                                        
				//  different from npart if one is dealing with a multi-file snapshot. 
  int flag_cooling;             //!< flags whether cooling was included  
  int num_files;                //!< number of files in multi-file snapshot 
  double BoxSize;               //!< box-size of simulation in case periodic boundaries were used 
  double Omega0;                //!< matter density in units of critical density 
  double OmegaLambda;           //!< cosmological constant parameter 
  double HubbleParam;           //!< Hubble parameter in units of 100 km/sec/Mpc 
  int flag_stellarage;          //!< flags whether the file contains formation times of star particles 
  int flag_metals;              //!< flags whether the file contains metallicity values for gas and star                                                                               				  //particles 
  unsigned int npartTotalHighWord[6];   //!< High word of the total number of particles of each type 
  int flag_entropy_instead_u;   //!< flags that IC-file contains entropy instead of u 
  int flag_doubleprecision;     //!< flags that snapshot contains double-precision instead of single precision 
  
  int flag_ic_info;   
  float lpt_scalingfactor;      //!< scaling factor for 2lpt initial conditions 
  char fill[18];                //!< fills to 256 Bytes 
  
  
  char names[15][2];
  }
  header1;                         //!< holds header for snapshot files 
*/


struct io_header_1
{
  int npart[6];
  double mass[6];
		
  double time;
  double redshift;
  int flag_sfr;
  int flag_feedback;
  unsigned int npartTotal[6];
			
  int flag_cooling;
  int num_files;
  double BoxSize;
  double Omega0;
  double OmegaLambda;
  double HubbleParam;
  int flag_stellarage;
  int flag_metals;
  unsigned int npartTotalHighWord[6];
  int flag_doubleprecision;

  int flag_ic_info;
  
  float lpt_scalingfactor;

  char fill[18];
  char names[15][2];
} header1;       



/* Count the number of lines in the file infile */

int countlines(char *infile)
{
  
  int NDAT,c;  
  static FILE *pf;
  
  if((pf=fopen(infile,"r"))==NULL)
    {
      printf(" ERROR: I can not open file %s\n",infile);
    }
  
  NDAT=0;
  
  while((c=fgetc(pf))!= EOF)
    {
      if(c == '\n')
	{
	  ++NDAT;
	}
    }
  
  fclose(pf);
  
  return NDAT;
}




void save_particles(char *fname)
{
  FILE *fd, *pf_subhalos;
  int i, k;
  float xyz[3];
  int4byte blklen;
  char subhalofile[500];
  
#define BLKLEN fwrite(&blklen, sizeof(blklen), 1, fd);


  if(!(fd = fopen(fname, "w")))
    {
      printf("error opening file %s\n", fname);
      exit(0);
    }

  printf("saveing initial conditions to file `%s'\n\n", fname);

  N_HALO = 0;

  /*  Include the subhalos in sec. 5 BLACKHOLE */

  sprintf(subhalofile,"%s","data_subhaloes.csv");
  N_BLACKHOLE = countlines(subhalofile);
  N_BLACKHOLE = 0;
  
  header1.npart[0] = header1.npartTotal[0] = N_GAS;
  header1.npart[1] = header1.npartTotal[1] = N_HALO;
  header1.npart[2] = header1.npartTotal[2] = N_DISK;
  header1.npart[3] = header1.npartTotal[3] = N_BULGE;
  header1.npart[4] = 0;
  header1.npart[5] = header1.npartTotal[5] = N_BLACKHOLE;

  
  pf_subhalos = fopen(subhalofile, "r");
  double *x_subhalo, *y_subhalo, *z_subhalo;
  double *vx_subhalo, *vy_subhalo, *vz_subhalo, *Msubhalo;

  x_subhalo = (double *) malloc((size_t ) sizeof(double)*N_BLACKHOLE);
  y_subhalo = (double *) malloc((size_t ) sizeof(double)*N_BLACKHOLE);
  z_subhalo = (double *) malloc((size_t ) sizeof(double)*N_BLACKHOLE);
    
  vx_subhalo = (double *) malloc((size_t ) sizeof(double)*N_BLACKHOLE);
  vy_subhalo = (double *) malloc((size_t ) sizeof(double)*N_BLACKHOLE);
  vz_subhalo = (double *) malloc((size_t ) sizeof(double)*N_BLACKHOLE);
  Msubhalo = (double *) malloc((size_t ) sizeof(double)*N_BLACKHOLE);

  
  for(i=0; i<N_BLACKHOLE; i++)
    {
      fscanf(pf_subhalos, "%lf %lf %lf %lf %lf %lf %lf",
	     &x_subhalo[i], &y_subhalo[i], &z_subhalo[i],
	     &vx_subhalo[i], &vy_subhalo[i], &vz_subhalo[i], &Msubhalo[i]);

      printf("%lf %lf %lf %lf %lf %lf %lf\n",
	     x_subhalo[i], y_subhalo[i], z_subhalo[i],
	     vx_subhalo[i], vy_subhalo[i], vz_subhalo[i], Msubhalo[i]);
      
    }
  fclose(pf_subhalos);

  getchar();
  
  header1.num_files = 0;

  for(i = 0; i < 6; i++)
    header1.mass[i] = 0;

  if(N_GAS)
    header1.mass[0] = mp_gas[1];

  if(N_HALO)
    header1.mass[1] = mp_halo[1];

  if(N_DISK)
    header1.mass[2] = mp_disk[1];

  if(N_BULGE)
    header1.mass[3] = mp_bulge[1];

  if(N_BLACKHOLE)
    header1.mass[5] = 0.0; /// Individual masses!!!

  header1.flag_sfr = 0;
  header1.flag_feedback = 0;
  header1.flag_cooling = 0;


  blklen = sizeof(header1);
  BLKLEN;
  fwrite(&header1, sizeof(header1), 1, fd);
  BLKLEN;


  blklen = 3 * (N_GAS + N_HALO + N_DISK + N_BULGE + N_BLACKHOLE) * sizeof(float);
  BLKLEN;
  for(i = 1; i <= N_GAS; i++)
    {
      xyz[0] = xp_gas[i];
      xyz[1] = yp_gas[i];
      xyz[2] = zp_gas[i];
      fwrite(xyz, sizeof(float), 3, fd);
    }
  for(i = 1; i <= N_HALO; i++)
    {
      xyz[0] = xp_halo[i];
      xyz[1] = yp_halo[i];
      xyz[2] = zp_halo[i];
      fwrite(xyz, sizeof(float), 3, fd);
    }
  FILE *PF;
  PF = fopen("disco.dat","w");
  for(i = 1; i <= N_DISK; i++)
    {
      xyz[0] = xp_disk[i];
      xyz[1] = yp_disk[i];
      xyz[2] = zp_disk[i];
      fwrite(xyz, sizeof(float), 3, fd);
      fprintf(PF,"%lf %lf %lf\n",xp_disk[i],yp_disk[i],zp_disk[i]);
    }
  fclose(PF);
  for(i = 1; i <= N_BULGE; i++)
    {
      xyz[0] = xp_bulge[i];
      xyz[1] = yp_bulge[i];
      xyz[2] = zp_bulge[i];
      fwrite(xyz, sizeof(float), 3, fd);
    }

  for(i = 0; i < N_BLACKHOLE; i++)
    {
      xyz[0] = x_subhalo[i];
      xyz[1] = y_subhalo[i];
      xyz[2] = z_subhalo[i];
      fwrite(xyz, sizeof(float), 3, fd);
    }

  BLKLEN;



  blklen = 3 * (N_GAS + N_HALO + N_DISK + N_BULGE + N_BLACKHOLE) * sizeof(float);
  BLKLEN;
  for(i = 1; i <= N_GAS; i++)
    {
      xyz[0] = vxp_gas[i];
      xyz[1] = vyp_gas[i];
      xyz[2] = vzp_gas[i];
      fwrite(xyz, sizeof(float), 3, fd);
    }
  for(i = 1; i <= N_HALO; i++)
    {
      xyz[0] = vxp_halo[i];
      xyz[1] = vyp_halo[i];
      xyz[2] = vzp_halo[i];
      fwrite(xyz, sizeof(float), 3, fd);
    }
  for(i = 1; i <= N_DISK; i++)
    {
      xyz[0] = vxp_disk[i];
      xyz[1] = vyp_disk[i];
      xyz[2] = vzp_disk[i];
      fwrite(xyz, sizeof(float), 3, fd);
    }
  for(i = 1; i <= N_BULGE; i++)
    {
      xyz[0] = vxp_bulge[i];
      xyz[1] = vyp_bulge[i];
      xyz[2] = vzp_bulge[i];
      fwrite(xyz, sizeof(float), 3, fd);
    }
  for(i = 0; i < N_BLACKHOLE; i++)
    {

      xyz[0] = vx_subhalo[i];
      xyz[1] = vy_subhalo[i];
      xyz[2] = vz_subhalo[i];
      
      fwrite(xyz, sizeof(float), 3, fd);
    }
  BLKLEN;


  blklen = (N_GAS + N_HALO + N_DISK + N_BULGE + N_BLACKHOLE) * sizeof(int);
  BLKLEN;
  k=0;
  for(i = 1; i <= (N_GAS + N_HALO + N_DISK + N_BULGE + N_BLACKHOLE); i++)
    {
      fwrite(&k, sizeof(int), 1, fd);	/* ID */
      k++;
    }
  BLKLEN;

  ////////////////////////////////////
  /* Individual masses for subhalos */

  if(N_BLACKHOLE)
    {
      blklen = (N_BLACKHOLE) * sizeof(float);
      BLKLEN;
      for(i = 0; i < N_BLACKHOLE; i++)
	{
	  xyz[0] = Msubhalo[i];
	  fwrite(xyz, sizeof(float), 1, fd);
	}
      BLKLEN;
    }

  /* Individual masses for subhalos */
  ////////////////////////////////////
  
  if(N_GAS)
    {
      blklen = (N_GAS) * sizeof(float);
      BLKLEN;
      for(i = 1; i <= N_GAS; i++)
	{
	  xyz[0] = u_gas[i];
	  fwrite(xyz, sizeof(float), 1, fd);
	}
      BLKLEN;
    }

  fclose(fd);
}

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "nrsrc/nr.h"
#include "nrsrc/nrutil.h"

#include "prototypes.h"
#include "globvars.h"



int main(int argc, char *argv[])
{
  char parameter_filename[100];
  char buf[1000];
  FILE *fd;
  int rep, i;

  if(argc != 2)
    {
      fprintf(stderr, "\n\nwrong argument(s).  Specify a parameterfile.\n\n");
      exit(0);
    }
  strcpy(parameter_filename, argv[1]);

  read_parameterfile(parameter_filename);

  init_units();			/* set system of units */

  structure();			/* determine structure of halo, disk, and bulge */

  init();			/* allocate arrays */

  init_clouds();


  if(M_DISK > 0)
    {

      init_central_densities();
      
      integrate_gasdensity();
      
      determine_cumulative_gasdensity();
      

      for(rep = 0; rep < 5; rep++)
	{
	  set_dummy_particles();
	  
	  compute_vertical_force_field();
	  
	  integrate_and_adjust();
	}
      
      dump_gas_density();
    }

  compute_radial_force_field();

  compute_vstream_gas();
  compute_velocity_dispersions_disk();
  compute_velocity_dispersions_halo();
  compute_velocity_dispersions_bulge();

  set_gas_positions();
  set_halo_positions();
  set_disk_positions();
  set_bulge_positions();

  set_halo_velocities();
  set_disk_velocities();
  set_gas_velocities();
  set_bulge_velocities();
  
  
  dump_veldisp_field();
  mkdir(OutputDir, 02755);

  sprintf(buf, "%s/%s", OutputDir, OutputFile);
  


  ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////  
  /* corrections to CM */ //QUIROGA
  
  double xcm, ycm, zcm;
  double vxcm, vycm, vzcm;
  double mtot;

  double galxcm, galycm, galzcm;
  double galvxcm, galvycm, galvzcm;
  double galmtot;

  double galxcm1, galycm1, galzcm1;
  double galvxcm1, galvycm1, galvzcm1;


  // Centering Halo

  if( N_HALO > 0 )
    {
      
      xcm = ycm = zcm = 0.0;
      vxcm = vycm = vzcm = 0.0;
      mtot = 0;
      
      for(i=1; i<=N_HALO; i++)
	{
	  xcm += mp_halo[i]*xp_halo[i];
	  ycm += mp_halo[i]*yp_halo[i];
	  zcm += mp_halo[i]*zp_halo[i];
	  
	  vxcm += mp_halo[i]*vxp_halo[i];
	  vycm += mp_halo[i]*vyp_halo[i];
	  vzcm += mp_halo[i]*vzp_halo[i];
	  
	  mtot += mp_halo[i];
	}
      
      xcm = xcm/mtot;
      ycm = ycm/mtot;
      zcm = zcm/mtot;
      
      vxcm = vxcm/mtot;
      vycm = vycm/mtot;
      vzcm = vzcm/mtot;
      
      printf("\nHalo berofe:\n");
      printf("Mtot=%g xcm=%g ycm=%g zcm=%g\n",mtot, xcm, ycm, zcm);
      printf("Mtot=%g vxcm=%g vycm=%g vzcm=%g\n",mtot, vxcm, vycm, vzcm);
      
      galxcm = xcm;
      galycm = ycm;
      galzcm = zcm;
      galmtot = mtot;
      
      for(i=1; i<=N_HALO; i++)
	{
	  xp_halo[i] = xp_halo[i] - xcm;
	  yp_halo[i] = yp_halo[i] - ycm;
	  zp_halo[i] = zp_halo[i] - zcm;
	  
	  vxp_halo[i] = vxp_halo[i] - vxcm;
	  vyp_halo[i] = vyp_halo[i] - vycm;
	  vzp_halo[i] = vzp_halo[i] - vzcm;
	}  
      
      xcm = ycm = zcm = 0.0;
      vxcm = vycm = vzcm = 0.0;
      mtot = 0;
      
      for(i=1; i<=N_HALO; i++)
	{
	  xcm += mp_halo[i]*xp_halo[i];
	  ycm += mp_halo[i]*yp_halo[i];
	  zcm += mp_halo[i]*zp_halo[i];
	  
	  vxcm += mp_halo[i]*vxp_halo[i];
	  vycm += mp_halo[i]*vyp_halo[i];
	  vzcm += mp_halo[i]*vzp_halo[i];
	  
	  mtot += mp_halo[i];
	}
      
      xcm = xcm/mtot;
      ycm = ycm/mtot;
      zcm = zcm/mtot;
      
      vxcm = vxcm/mtot;
      vycm = vycm/mtot;
      vzcm = vzcm/mtot;
      
      printf("Halo after:\n");
      printf("Mtot=%g xcm=%g ycm=%g zcm=%g\n",mtot, xcm, ycm, zcm);
      printf("Mtot=%g vxcm=%g vycm=%g vzcm=%g\n",mtot, vxcm, vycm, vzcm);
      
      galxcm1 = xcm;
      galycm1 = ycm;
      galzcm1 = zcm;
      
    }

  // Centering disk

  if( N_DISK>0 )
    {
      
      xcm = ycm = zcm = 0.0;
      vxcm = vycm = vzcm = 0.0;
      mtot = 0;
      
      for(i=1; i<=N_DISK; i++)
	{
	  xcm += mp_disk[i]*xp_disk[i];
	  ycm += mp_disk[i]*yp_disk[i];
	  zcm += mp_disk[i]*zp_disk[i];
	  
	  vxcm += mp_disk[i]*vxp_disk[i];
	  vycm += mp_disk[i]*vyp_disk[i];
	  vzcm += mp_disk[i]*vzp_disk[i];
	  
	  mtot += mp_disk[i];
	}
      
      xcm = xcm/mtot;
      ycm = ycm/mtot;
      zcm = zcm/mtot;
      
      vxcm = vxcm/mtot;
      vycm = vycm/mtot;
      vzcm = vzcm/mtot;
      
      printf("\nDisk berofe:\n");
      printf("Mtot=%g xcm=%g ycm=%g zcm=%g\n",mtot, xcm, ycm, zcm);
      printf("Mtot=%g vxcm=%g vycm=%g vzcm=%g\n",mtot, vxcm, vycm, vzcm);
      
      galxcm = galxcm + xcm;
      galycm = galxcm + ycm;
      galzcm = galxcm + zcm;
      galmtot = galmtot + mtot;
      
      for(i=1; i<=N_DISK; i++)
	{
	  xp_disk[i] = xp_disk[i] - xcm;
	  yp_disk[i] = yp_disk[i] - ycm;
	  zp_disk[i] = zp_disk[i] - zcm;
	  
	  vxp_disk[i] = vxp_disk[i] - vxcm;
	  vyp_disk[i] = vyp_disk[i] - vycm;
	  vzp_disk[i] = vzp_disk[i] - vzcm;
	}
      
      xcm = ycm = zcm = 0.0;
      vxcm = vycm = vzcm = 0.0;
      mtot = 0;
      
      for(i=1; i<=N_DISK; i++)
	{
	  xcm += mp_disk[i]*xp_disk[i];
	  ycm += mp_disk[i]*yp_disk[i];
	  zcm += mp_disk[i]*zp_disk[i];
	  
	  vxcm += mp_disk[i]*vxp_disk[i];
	  vycm += mp_disk[i]*vyp_disk[i];
	  vzcm += mp_disk[i]*vzp_disk[i];
	  
	  mtot += mp_disk[i];
	}
      
      xcm = xcm/mtot;
      ycm = ycm/mtot;
      zcm = zcm/mtot;
      
      vxcm = vxcm/mtot;
      vycm = vycm/mtot;
      vzcm = vzcm/mtot;
      
      printf("Disk after:\n");
      printf("Mtot=%g xcm=%g ycm=%g zcm=%g\n",mtot, xcm, ycm, zcm);
      printf("Mtot=%g vxcm=%g vycm=%g vzcm=%g\n",mtot, vxcm, vycm, vzcm);
      
      galxcm1 = galxcm1 + xcm;
      galycm1 = galxcm1 + ycm;
      galzcm1 = galxcm1 + zcm;
      
    }

  // Centering Bulge

  if( N_BULGE>0 )
    {

      xcm = ycm = zcm = 0.0;
      vxcm = vycm = vzcm = 0.0;
      mtot = 0;
      
      for(i=1; i<=N_BULGE; i++)
	{
	  xcm += mp_bulge[i]*xp_bulge[i];
	  ycm += mp_bulge[i]*yp_bulge[i];
	  zcm += mp_bulge[i]*zp_bulge[i];
	  
	  vxcm += mp_bulge[i]*vxp_bulge[i];
	  vycm += mp_bulge[i]*vyp_bulge[i];
	  vzcm += mp_bulge[i]*vzp_bulge[i];
	  
	  mtot += mp_bulge[i];
	}
      
      xcm = xcm/mtot;
      ycm = ycm/mtot;
      zcm = zcm/mtot;
      
      vxcm = vxcm/mtot;
      vycm = vycm/mtot;
      vzcm = vzcm/mtot;
      
      printf("\nBulge before:\n");
      printf("Mtot=%g xcm=%g ycm=%g zcm=%g\n",mtot, xcm, ycm, zcm);
      printf("Mtot=%g vxcm=%g vycm=%g vzcm=%g\n",mtot, vxcm, vycm, vzcm);
      
      galxcm = galxcm + xcm;
      galycm = galxcm + ycm;
      galzcm = galxcm + zcm;
      galmtot = galmtot + mtot;
      
      for(i=1; i<=N_BULGE; i++)
	{
	  xp_bulge[i] = xp_bulge[i] - xcm;
	  yp_bulge[i] = yp_bulge[i] - ycm;
	  zp_bulge[i] = zp_bulge[i] - zcm;
	  
	  vxp_bulge[i] = vxp_bulge[i] - vxcm;
	  vyp_bulge[i] = vyp_bulge[i] - vycm;
	  vzp_bulge[i] = vzp_bulge[i] - vzcm;
	}
      
      xcm = ycm = zcm = 0.0;
      vxcm = vycm = vzcm = 0.0;
      mtot = 0;
      
      for(i=1; i<=N_BULGE; i++)
	{
	  xcm += mp_bulge[i]*xp_bulge[i];
	  ycm += mp_bulge[i]*yp_bulge[i];
	  zcm += mp_bulge[i]*zp_bulge[i];
	  
	  vxcm += mp_bulge[i]*vxp_bulge[i];
	  vycm += mp_bulge[i]*vyp_bulge[i];
	  vzcm += mp_bulge[i]*vzp_bulge[i];
	  
	  mtot += mp_bulge[i];
	}
      
      xcm = xcm/mtot;
      ycm = ycm/mtot;
      zcm = zcm/mtot;
      
      vxcm = vxcm/mtot;
      vycm = vycm/mtot;
      vzcm = vzcm/mtot;
      
      printf("Bulge after:\n");
      printf("Mtot=%g xcm=%g ycm=%g zcm=%g\n",mtot, xcm, ycm, zcm);
      printf("Mtot=%g vxcm=%g vycm=%g vzcm=%g\n",mtot, vxcm, vycm, vzcm);
      
      galxcm1 = galxcm1 + xcm;
      galycm1 = galxcm1 + ycm;
      galzcm1 = galxcm1 + zcm;
      
    }

  // Centering gas disk
  
  if( N_GAS>0 )  
    {
      xcm = ycm = zcm = 0.0;
      vxcm = vycm = vzcm = 0.0;
      mtot = 0;
      
      for(i=1; i<=N_GAS; i++)
	{
	  xcm += mp_gas[i]*xp_gas[i];
	  ycm += mp_gas[i]*yp_gas[i];
	  zcm += mp_gas[i]*zp_gas[i];
	  
	  vxcm += mp_gas[i]*vxp_gas[i];
	  vycm += mp_gas[i]*vyp_gas[i];
	  vzcm += mp_gas[i]*vzp_gas[i];
	  
	  mtot += mp_gas[i];
	}
      
      xcm = xcm/mtot;
      ycm = ycm/mtot;
      zcm = zcm/mtot;
      
      vxcm = vxcm/mtot;
      vycm = vycm/mtot;
      vzcm = vzcm/mtot;
      
      printf("\nGas berofe:\n");
      printf("Mtot=%g xcm=%g ycm=%g zcm=%g\n",mtot, xcm, ycm, zcm);
      printf("Mtot=%g vxcm=%g vycm=%g vzcm=%g\n",mtot, vxcm, vycm, vzcm);
      
      galxcm = galxcm + xcm;
      galycm = galxcm + ycm;
      galzcm = galxcm + zcm;
      galmtot = galmtot + mtot;
      
      for(i=1; i<=N_GAS; i++)
	{
	  xp_gas[i] = xp_gas[i] - xcm;
	  yp_gas[i] = yp_gas[i] - ycm;
	  zp_gas[i] = zp_gas[i] - zcm;
	  
	  vxp_gas[i] = vxp_gas[i] - vxcm;
	  vyp_gas[i] = vyp_gas[i] - vycm;
	  vzp_gas[i] = vzp_gas[i] - vzcm;
	}
      
      xcm = ycm = zcm = 0.0;
      vxcm = vycm = vzcm = 0.0;
      mtot = 0;
      
      for(i=1; i<=N_GAS; i++)
	{
	  xcm += mp_gas[i]*xp_gas[i];
	  ycm += mp_gas[i]*yp_gas[i];
	  zcm += mp_gas[i]*zp_gas[i];
	  
	  vxcm += mp_gas[i]*vxp_gas[i];
	  vycm += mp_gas[i]*vyp_gas[i];
	  vzcm += mp_gas[i]*vzp_gas[i];
	  
	  mtot += mp_gas[i];
	}
      
      xcm = xcm/mtot;
      ycm = ycm/mtot;
      zcm = zcm/mtot;
      
      vxcm = vxcm/mtot;
      vycm = vycm/mtot;
      vzcm = vzcm/mtot;
      
      printf("Gas disk after:\n");
      printf("Mtot=%g xcm=%g ycm=%g zcm=%g\n",mtot, xcm, ycm, zcm);
      printf("Mtot=%g vxcm=%g vycm=%g vzcm=%g\n",mtot, vxcm, vycm, vzcm);
      
      galxcm1 = galxcm1 + xcm;
      galycm1 = galxcm1 + ycm;
      galzcm1 = galxcm1 + zcm;

    }

  printf("\nCm galaxy before:\n");
  printf("Mtot=%g xcm=%g ycm=%g zcm=%g\n",galmtot, galxcm, galycm, galzcm);
  printf("Mtot=%g vxcm=%g vycm=%g vzcm=%g\n",galmtot, galvxcm, galvycm, galvzcm);
  
  printf("Cm galaxy after:\n");
  printf("Mtot=%g xcm=%g ycm=%g zcm=%g\n",galmtot, galxcm1, galycm1, galzcm1);
  printf("Mtot=%g vxcm=%g vycm=%g vzcm=%g\n",galmtot, galvxcm1, galvycm1, galvzcm1);
  
  
  ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////

  

  save_particles(buf);
  sprintf(buf, "%s/%s_curve.txt", OutputDir,OutputFile);

  if((fd = fopen(buf, "w")))
    {
      printf("writing circular velocity curve + Toomre's Q\n");
      write_header(fd);
      //plot_circular_speeds(fd);
      plot_circular_speeds();
      //plot_toomre_stability(fd);
      plot_toomre_stability();
      fclose(fd);
    }
  else
    {
      fprintf(stderr, "Can't open file '%s'.\n", buf);
      exit(0);
    }
  printf("done.\n");


  printf("Disk scale length: %g\n", H);



  return 0;
}



void write_header(FILE * fd)
{
  fprintf(fd, "%g\n", CC);
  fprintf(fd, "%g\n", V200);
  fprintf(fd, "%g\n", LAMBDA);
  fprintf(fd, "%g\n", MD);
  fprintf(fd, "%g\n", JD);
  fprintf(fd, "%g\n", MB);
  fprintf(fd, "%g\n", DiskHeight);
  fprintf(fd, "%g\n", BulgeSize);
  fprintf(fd, "\n%g\n", R200);
  fprintf(fd, "%g\n\n", H);
}

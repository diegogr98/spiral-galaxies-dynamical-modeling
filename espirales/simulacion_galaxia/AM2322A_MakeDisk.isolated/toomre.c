#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "nrsrc/nr.h"
#include "nrsrc/nrutil.h"


#include "prototypes.h"
#include "globvars.h"



double epicyclic_kappa2(double R)
{
  double dR, dphi, dphi_;

  if(R > 0)
    {
      dR = R * 0.01;

      dphi = comp_Dphi_R(R, 0);

      /*comp_Dphi_R_halo(R, 0) + comp_Dphi_R_disk(R, 0);
       */

      dphi_ = comp_Dphi_R(R + dR, 0);

      /*
         comp_Dphi_R_halo(R + dR, 0) + comp_Dphi_R_disk(R + dR, 0);
       */

      return 3 * dphi / R + (dphi_ - dphi) / dR;

    }
  else
    return 0;
}



//void plot_toomre_stability(FILE * fd)
void plot_toomre_stability(void)
{
  double *Q;
  int i;
  double Sigma0;
  int count;
  FILE *pf;
  char buff[1000];

  for(i = 0, count = 0; i <= RSIZE; i++)
    if(list_R[i] <= 6 * H)
      count++;

  Sigma0 = (M_DISK) / (2 * PI * H * H);
  Q = dvector(0, RSIZE);


  for(i = 0; i < count; i++)
    Q[i] = RadialDispersionFactor * sqrt(VelDispRz_disk[i][0]) * sqrt(epi_kappa2[i]) / (3.36 * G * Sigma0 * exp(-list_R[i] / H));


  /*
  fprintf(fd, "\n%d\n", count);

  for(i = 0; i < count; i++)
    fprintf(fd, "%g\n", list_R[i]);

    for(i = 0; i < count; i++)
    fprintf(fd, "%g\n", Q[i]);
  */
  
  sprintf(buff,"%s/%s_MyToomre.dat",OutputDir, OutputFile);
  pf = fopen(buff,"w");
  for(i = 0; i < count; i++)
    fprintf(pf, "%g %g\n", list_R[i], Q[i]);
  fclose(pf);
}







void plot_circular_speeds(void)
//void plot_circular_speeds(FILE * fd)
{
  int i;
  double R;
  double RMAX;
  double vc2_1, vc2_2, vc2_3, vc2_4;
  FILE *pf;
  char buff[1000];
  

#define POINTS 1000

  RMAX = R200;

  //fprintf(fd, "%d\n", POINTS);

  /*
    for(i = 1; i <= POINTS; i++)
    {
      R = (RMAX / POINTS) * i;
      fprintf(fd, "%f\n", R);
    }
  for(i = 1; i <= POINTS; i++)
    {
      R = (RMAX / POINTS) * i;
      fprintf(fd, "%f\n", vc2 = R * comp_Dphi_R(R, 0));
    }

  for(i = 1; i <= POINTS; i++)
    {
      R = (RMAX / POINTS) * i;
      fprintf(fd, "%f\n", vc2 = R * comp_Dphi_R_disk_tree(R, 0));
    }
  for(i = 1; i <= POINTS; i++)
    {
      R = (RMAX / POINTS) * i;
      fprintf(fd, "%f\n", vc2 = R * comp_Dphi_R_halo(R, 0));
    }

  for(i = 1; i <= POINTS; i++)
    {
      R = (RMAX / POINTS) * i;
      fprintf(fd, "%f\n", vc2 = R * comp_Dphi_R_bulge(R, 0));
    }
  */

  sprintf(buff,"%s/%s_MyVelCurve.dat",OutputDir, OutputFile);
  pf = fopen(buff, "w");
  for(i = 1; i <= POINTS; i++)
    {
      R = (RMAX / POINTS) * i;
      vc2_1 = R * comp_Dphi_R(R, 0);
      vc2_2 = R * comp_Dphi_R_disk_tree(R, 0);
      vc2_3 = R * comp_Dphi_R_halo(R, 0);
      vc2_4 = R * comp_Dphi_R_bulge(R, 0);
      
      fprintf(pf, "%16.8f %16.8f %16.8f %16.8f %16.8f\n", R, vc2_1, vc2_2, vc2_3, vc2_4);
    }
  fclose(pf);

#undef POINTS
}

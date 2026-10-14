#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "nrsrc/nr.h"
#include "nrsrc/nrutil.h"

#include "prototypes.h"
#include "globvars.h"

void structure(void)
{
  double jhalo, jdisk, jd;
  double hnew, dh;

   //R200 = V200 / (10 * H0);
  //M200 = pow(V200, 3) / (10 * G * H0);
  
  M_DISK = MD * M200; 
  M_BULGE = MB * M200;

  M_BLACKHOLE = MBH * M200;
  
  RS = R200 / CC;

  if(MBH>0)
    N_BLACKHOLE = 1;
  else
    N_BLACKHOLE = 0;

  M_GAS = M_DISK * GasFraction;
  M_HALO = M200 - M_DISK - M_BULGE - M_BLACKHOLE;
  
  RH = RS * sqrt(2 * (log(1 + CC) - CC / (1 + CC)));	/* scale length of Hernquist profile */
  
  M_HALO = (M200*(R200 + RH)*(R200 + RH)/(R200*R200)) - M_DISK - M_BULGE - M_BLACKHOLE;
    
  jhalo = LAMBDA * sqrt(G) * pow(M200, 1.5) * sqrt(2 *  R200 / fc(CC));
  jdisk = JD * jhalo;
  
  printf("fc=%g\n", fc(CC));
  
  halo_spinfactor = 1.5 * LAMBDA * sqrt(2 * CC / fc(CC)) * pow(log(1 + CC) - CC / (1 + CC), 1.5) / gc(CC);
  
  printf("halo_spinfactor %g\n", halo_spinfactor);
  printf("R200: %g\n", R200);
  printf("M200: %g\n", M200);
  printf("RH: %g\n", RH);
  printf("Rs: %g\n", RS);
  printf("Total mass is: %g\n", M_HALO); 
  
  H = sqrt(2.0) / 2.0 * LAMBDA / fc(CC) * R200;	/* first guess for disk scale length */

  printf("H = %g\n",H);

  Z0 = DiskHeight * H;		/* sets disk thickness for stars */
  A = BulgeSize * H;		/* sets bulge size */
  
  if(M_DISK > 0)
    {
      do
	{
	  jd = disk_angmomentum();	/* computes disk momentum */
	  
	  hnew = jdisk / jd * H;
	  
	  dh = hnew - H;
	  
	  if(fabs(dh) > 0.5 * H)
	    {
	      dh = 0.5 * H * dh / fabs(dh);
	    }
	  else
	    dh = dh * 0.1;
	  
	  H = H + dh;
	  
	  printf("Jd/J=%g   hnew: %g  \n", jd / jhalo, H);
	  
	  
	  Z0 = DiskHeight * H;	/* sets disk thickness */
	  A = BulgeSize * H;	/* sets bulge size */
	}
      while(fabs(dh) / H > 1e-5);
    }
  
  
  printf(" after iteration disk scale length is %f\n",H);

}








double halo_mass(double r)
{
  return M_HALO * pow(r / (r + RH), 2);
}


double mass_total_dark(double radius)
{
  return halo_mass(radius);
}


double halo_q_to_r(double q)
{
  if(q>0)
    return RH * (q + sqrt(q)) / (1 - q);
  else
    return 0;
}


double halo_rho(double r)
{
  if(r > 0)
    return M_HALO / (2 * M_PI) * RH / r / pow(r + RH, 3);
  else
    return 0;
}









double disk_angmomentum(void)
{
  double jdisk_int(double);

  return M_DISK * qromb(jdisk_int, 0, dmin(30 * H, R200));
}




double jdisk_int(double x)
{
  double vc2, Sigma0, vc, y;
  double mass_cumulative_disk(double);
  double mass_total_dark(double radius);

  if(x > 1e-20)
    {
      vc2 = G * (mass_total_dark(x) + mass_cumulative_bulge(x)) / x;
    }
  else
    vc2 = 0;


  if(vc2 < 0)
    {
      vc2 = 0;
      printf("wwww\n");
      exit(0);
    }

  Sigma0 = (M_DISK) / (2 * PI * H * H);
  y = x / (2 * H);

  if(y > 1e-4)
    vc2 += x * 2 * PI * G * Sigma0 * y * (bessi0(y) * bessk0(y) - bessi1(y) * bessk1(y));

  vc = sqrt(vc2);

  return pow(x / H, 2) * vc * exp(-x / H);
}


double fc(double c)
{
  return c * (0.5 - 0.5 / pow(1 + c, 2) - log(1 + c) / (1 + c)) / pow(log(1 + c) - c / (1 + c), 2);
}

double gc(double c)
{
  double gc_int(double);

  return qromb(gc_int, 0, c);
}

double gc_int(double x)
{
  return pow(log(1 + x) - x / (1 + x), 0.5) * pow(x, 1.5) / pow(1 + x, 2);
}

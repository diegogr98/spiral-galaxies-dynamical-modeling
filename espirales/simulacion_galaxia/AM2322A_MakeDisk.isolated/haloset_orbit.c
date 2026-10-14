#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "nrsrc/nr.h"
#include "nrsrc/nrutil.h"

#include "prototypes.h"
#include "globvars.h"


#if defined(EXACT_HERNQUIST) || defined(OUTPUT_DF)
/*E to q conversion*/
double halo_E_to_q(double E)
{
 return sqrt(-E*RH / (G*M_HALO));
}


/*Hernquist potential*/
double halo_pot(double rad)
{
 return - G*M_HALO / (rad + RH);
}


/*Hernquist density of states (as a function of q)*/
double halo_g_q(double q)
{
 double pre = 2*sqrt(2)*PI*PI*RH*RH*RH*sqrt(G*M_HALO/RH);

 return pre * (3*(8*q*q*q*q - 4*q*q + 1)*acos(q) - q*sqrt(1 - q*q)*(4*q*q- 1 )*(2*q*q + 3)) / (3*q*q*q*q*q);
}


/*Hernquist distribution function (as a function of q)*/
double halo_f_q(double q)
{
 double pre = (M_HALO / (RH*RH*RH)) / (4*PI*PI*PI*pow(2*G*M_HALO/RH,1.5));
 
 return pre*(3*asin(q)+q*sqrt(1-q*q)*(1-2*q*q)*(8*q*q*q*q-8*q*q-3))/pow(1-q*q,2.5);
}


/*Hernquist distribution function (as a function of radius and velocity)*/
double halo_f(double rad, double vel)
{
 double E = 0.5 * vel * vel + halo_pot(rad);
 double q = halo_E_to_q(E);

 return halo_f_q(q);
}


/*generate velocities with von Neumann rejection technique*/
double generate_v(double rad)
{
 double pot = halo_pot(rad);
 double v_max = sqrt(-2*pot); // escape velocity
 double v_guess, x_aux;
 double f_max = v_max * v_max * halo_f(rad,0);  

 v_guess = drand48() * v_max;
 x_aux = drand48() * f_max;

 while (x_aux >= v_guess * v_guess * halo_f(rad, v_guess))
 {
   v_guess = drand48() * v_max;
   x_aux = drand48() * f_max;
  }  
 return v_guess;
}


/*function that writes out energy distribution dM/dE (both theoretical and the one given by MakeNewDisk)*/
void write_f()
{
 double E, q;
 int i, j;
 int bins=500, bin[bins];
 double fac = G * M_HALO / RH;
 double Emin = 0;
 double Emax = 1.0 * fac;
 double dE = (Emax-Emin) / bins;
 double mass = M_HALO / N_HALO;

 FILE *fd = fopen("DF.txt","w");
 
 for (j = 0; j < bins; j++)
  bin[j] = 0;
  
 for (i = 1; i <= N_HALO; i++)
 {
   E = -(0.5 * (vxp_halo[i]*vxp_halo[i] + vyp_halo[i]*vyp_halo[i] + vzp_halo[i]*vzp_halo[i]) + 
         halo_pot(sqrt(xp_halo[i]*xp_halo[i] + yp_halo[i]*yp_halo[i] + zp_halo[i]*zp_halo[i])));

   for (j = 0; j < bins; j++) 
   {
    if (E > 0.0)
     if ((Emin + j*dE < E) && (Emin + (j+1)*dE > E)) {bin[j]++; j = bins+1;}
   }
 }
 for (j = 0; j < bins; j++)
 {
  E = (Emin + (2*j+1)*dE / 2.0);
  q = halo_E_to_q(-E); 
  fprintf(fd, "%e %e %e %e %e\n", 
               E / fac, 
               ((bin[j]*mass) / dE) / halo_g_q(q) , 
               halo_f_q(q), 
               (bin[j]*mass) / dE, 
               halo_g_q(q) * halo_f_q(q));
 } 
 fclose(fd);
}
#endif

#ifdef EXACT_HERNQUIST
double set_halo_velocities(void) 
{
 double vr, r, phi, theta;
 int i;
 
 
 printf("set halo velocities...\n");
   
 srand48(22); 
 
 for(i = 1; i <= N_HALO; i++)
   {
    phi = drand48() * PI * 2;
    theta = acos(drand48() * 2 - 1);

    r = sqrt(xp_halo[i]*xp_halo[i] + yp_halo[i]*yp_halo[i] + zp_halo[i]*zp_halo[i]);

    vr = generate_v(r);
    
    /*isotropic velocity distribution*/
    vxp_halo[i] = vr*sin(theta)*cos(phi);
    vyp_halo[i] = vr*sin(theta)*sin(phi);
    vzp_halo[i] = vr*cos(theta);
   }

 printf("done.\n");
 fflush(stdout);

#ifdef OUTPUT_DF
 write_f();
#endif 
 return 0;
}
#endif


#ifndef EXACT_HERNQUIST
static double R, z;

double set_halo_velocities(void)
{
  int i, rangecount=0;
  long dum;
  int iz, ir;
  double ur, uz;
  double vdisp_rz, vdisp_phi, vstream_phi;
  double vr, vphi;
  double vx, vy, vz;


  if(N_HALO == 0)
    return 0;

  dum = drand48() * 1e8;

  printf("set halo velocities...\n");
  fflush(stdout);

  for(i = 1; i <= N_HALO; i++)
    {
      R = sqrt(xp_halo[i] * xp_halo[i] + yp_halo[i] * yp_halo[i]);
      z = fabs(zp_halo[i]);

      if(R < Baselen)
	ir = 0;
      else
	ir = (log(R) - log(Baselen)) / (log(LL) - log(Baselen)) * (RSIZE - 1) + 1;

      ur = (R - list_R[ir]) / (list_R[ir + 1] - list_R[ir]);

      if(z < Baselen)
	iz = 0;
      else
	iz = (log(z) - log(Baselen)) / (log(LL) - log(Baselen)) * (ZSIZE - 1) + 1;

      uz = (z - list_z[iz]) / (list_z[iz + 1] - list_z[iz]);

      if(ir < 0 || ir >= RSIZE || ur < 0 || ur > 1)
	{
	  rangecount++;
	  /* 
	     printf("ir=%d out of range R=%g\n", ir, R);
	  */
	  vxp_halo[i] = 0;
	  vyp_halo[i] = 0;
	  vzp_halo[i] = 0;
	  continue;
	}

      if(iz < 0 || iz >= ZSIZE || uz < 0 || uz > 1)
	{
	  rangecount++;
	  /*
	  printf("iz=%d out of range Z=%g\n", iz, z);
	  */
	  vxp_halo[i] = 0;
	  vyp_halo[i] = 0;
	  vzp_halo[i] = 0;
	  continue;
	}


      vdisp_rz = VelDispRz_halo[ir][iz] * (1 - ur) * (1 - uz)
	+ VelDispRz_halo[ir + 1][iz] * (ur) * (1 - uz)
	+ VelDispRz_halo[ir][iz + 1] * (1 - ur) * (uz) + VelDispRz_halo[ir + 1][iz + 1] * (ur) * (uz);

      vdisp_phi = VelDispPhi_halo[ir][iz] * (1 - ur) * (1 - uz)
	+ VelDispPhi_halo[ir + 1][iz] * (ur) * (1 - uz)
	+ VelDispPhi_halo[ir][iz + 1] * (1 - ur) * (uz) + VelDispPhi_halo[ir + 1][iz + 1] * (ur) * (uz);

      vstream_phi = VelStreamPhi_halo[ir][iz] * (1 - ur) * (1 - uz)
	+ VelStreamPhi_halo[ir + 1][iz] * (ur) * (1 - uz)
	+ VelStreamPhi_halo[ir][iz + 1] * (1 - ur) * (uz) + VelStreamPhi_halo[ir + 1][iz + 1] * (ur) * (uz);

      if(vdisp_rz < 0)
	{
	  printf("in halo: vdisp_rz:%g   %g %g %d %d \n", vdisp_rz, ur, uz, ir, iz);
	  vdisp_rz = -vdisp_rz;
	}
      if(vdisp_phi < 0)
	{
	  printf("in halo: vdisp_phi:%g  %g %g %d %d\n", vdisp_phi, ur, uz, ir, iz);

	  vdisp_phi = -vdisp_phi;
	}

      vr = gasdev(&dum) * sqrt(vdisp_rz);
      vz = gasdev(&dum) * sqrt(vdisp_rz);

      vphi = vstream_phi + gasdev(&dum) * sqrt(vdisp_phi);


      vx = vr * xp_halo[i] / R - vphi * yp_halo[i] / R;
      vy = vr * yp_halo[i] / R + vphi * xp_halo[i] / R;

      vxp_halo[i] = vx;
      vyp_halo[i] = vy;
      vzp_halo[i] = vz;
    }

  printf("rangecount=%d\n", rangecount);
  printf("done.\n");
  fflush(stdout);
#ifdef OUTPUT_DF
  write_f(); 
#endif  
  return 0;
}
#endif





double set_halo_positions(void)
{
  int i, countr, countz;
  double q, R, phi, theta;
  double halo_q_to_r(double q);

  if(N_HALO == 0)
    return 0;

  srand48(22);

  printf("set halo positions...\n");

  for(i = 1, countr = countz = 0; i <= N_HALO;)
    {

      do
	{
	  q = drand48();

	  R = halo_q_to_r(q);
	}
      while((R > 2 * R200)); 
           
            
      phi = drand48() * PI * 2;
      theta = acos(drand48() * 2 - 1);

      xp_halo[i] = R * sin(theta) * cos(phi);
      yp_halo[i] = R * sin(theta) * sin(phi);
      zp_halo[i] = R * cos(theta);

      i++;
    }

  for(i = 1; i <= N_HALO; i++)
    mp_halo[i] = M_HALO / N_HALO;

  return 0;
}



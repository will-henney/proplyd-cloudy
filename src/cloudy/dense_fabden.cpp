/* This file is part of Cloudy and is copyright (C) 1978-2005 by Gary J. Ferland.
 * For conditions of distribution and use, see copyright notice in license.txt */
/*fabden called by dlaw command, returns density for any density law */
#include "cddefines.h"
#include "physconst.h"
#include "dense.h"
#include "rfield.h"
#include "phycon.h"
#include "geometry.h"
#include "struc.h"
#include "timesc.h"
#include "thirdparty.h"

/* Table of dimensionless radius from center of globule: R = r/r0 */ 
static double Rtab [50] = {1.0,           1.0018393337, 1.00724258596, 1.01608739516, 1.02830925054, 
						   1.04389360183, 1.06287048426, 1.08531103017, 1.11132544558, 1.1410621687, 
						   1.1747080229,  1.21248924461, 1.25467331637, 1.30157157394, 1.35354258677,
						   1.41099633775, 1.47439925166, 1.54428014425, 1.62123718677, 1.70594600442,
						   1.79916905279, 1.90176644544, 2.01470843744, 2.1390898074, 2.27614642265, 
						   2.42727432236, 2.59405171128, 2.7782643246, 2.98193470474, 3.20735602471, 
						   3.45713120398, 3.73421819406, 4.04198246663, 4.38425792141, 4.7654176498,
						   5.19045625101, 5.66508570687, 6.19584719123, 6.79024163074, 7.45688236054,
						   8.20567384917, 9.04802122222, 9.99707622054, 11.0680263174, 12.2784350283,
						   13.6486430234, 15.2022415548, 16.9666320047, 18.9736881379, 21.260541};
/* Table of dimensionless velocities: U = vc0 
   Note that U = 1 at R = 1
*/ 
static double Utab [50] = {1.0, 1.0612244898, 1.12244897959, 1.18367346939, 1.24489795918, 
						   1.30612244898, 1.36734693878, 1.42857142857, 1.48979591837, 1.55102040816,
						   1.61224489796, 1.67346938776, 1.73469387755, 1.79591836735, 1.85714285714,
						   1.91836734694, 1.97959183673, 2.04081632653, 2.10204081633, 2.16326530612,
						   2.22448979592, 2.28571428571, 2.34693877551, 2.40816326531, 2.4693877551, 
						   2.5306122449, 2.59183673469, 2.65306122449, 2.71428571429, 2.77551020408,
						   2.83673469388, 2.89795918367, 2.95918367347, 3.02040816327, 3.08163265306,
						   3.14285714286, 3.20408163265, 3.26530612245, 3.32653061224, 3.38775510204,
						   3.44897959184, 3.51020408163, 3.57142857143, 3.63265306122, 3.69387755102,
						   3.75510204082, 3.81632653061, 3.87755102041, 3.9387755102, 4.0};


double dense_proplyd(double depth, double r0, double Rmax, double rho_Rmax)
{
  
  /* New variables for interpolating velocity of globule flow 14 Jun 2011 */ 
  double R; /* Dimensionless globule radius */
  double U; /* Dimensionless velocity */
  static double rho_m; /* Density at the sonic point, R=1 */
  static double c_m; /* Isothermal sound speed at sonic point, R=1*/
  double fabden_v;
  static double U_Rmax;
  static double Uspline [50] ; 		// Coefficients for spline interpolation
  

  /* 
	 14 Jun 2011 - Nahiely y Will
	 Accelerating divergent photoevaporation flow (applied to proplyds)
	 
	 1. Fully ionized region: Density follows from mass conservation
	 in spherical geometry Velocity law is prescribed: U(R)
	 
	 2. Partially ionized region (TODO): Density is function of sound speed
	 
  */

  /* Calculamos R en funcion de la profundidad z del paso en el que vaya el modelo*/
  R = Rmax-((depth)/r0);

  /* First time called, so set up spline interpolation for U vs R using the tables Utab and Rtab */
  if (nzone == 0) spline(Rtab, Utab, 50, 2e31, 2e31, Uspline);

  // Find the velocity by spline interpolation
  splint(Rtab, Utab, Uspline, 50, R, &U);
  
  if (nzone == 0) {
	U_Rmax = U;
	fabden_v = rho_Rmax;
  }	else {
	if (R >= 1) {
	  // Calculate density from continuity in spherical symmetry
	  fabden_v = rho_Rmax * (U_Rmax/U) * POW2(Rmax/R) ;
	  c_m = timesc.sound_speed_isothermal;
	  rho_m = fabden_v;
	}	else {
	  fabden_v = (rho_m)/(1-sqrt(1.2-POW2(timesc.sound_speed_isothermal/c_m)));
	}
  }	
  fprintf(ioQQQ,"%.4i %.2e  %.2e %.4f \n", nzone, c_m,rho_m, R);

  return(fabden_v);

}


double dense_proplyd_LinealInterpol(double depth, double r0, double Rmax, double rho_Rmax)
{
  
  /* New variables for interpolating velocity of globule flow 14 Jun 2011 */ 
  double R; /* Dimensionless globule radius */
  double U=3.5; /* Dimensionless velocity */
  static double U_Rmax;
  double a;
  int j; /* index into Rtab and Utab */
  bool lgHit; /* true when we have found the right place in the table */
  double fabden_v;
  
  /* 
	 14 Jun 2011 - Nahiely y Will
	 Accelerating divergent photoevaporation flow (applied to proplyds)
	 
	 1. Fully ionized region: Density follows from mass conservation in spherical geometry
	 Velocity law is prescribed: U(R)
	 
	 2. Partially ionized region (TODO): Density is function of sound speed
	 
  */

  /* Calculamos R en funcion de la profundidad z del paso en el que vaya el modelo*/
  R = Rmax-((depth)/r0);

  if (nzone == 0) {
    U_Rmax = U;
    fabden_v = rho_Rmax;
  }
  else {
	/* Buscamos si R esta dentro del ango de valores de la tabla que hemso dado */
	if (R < Rtab[0] || R >= Rtab[49])
	  {		
		fprintf( ioQQQ, " requested radius outside range of Rtab\n" );
		fprintf( ioQQQ, " radius was%10.2e min, max=%10.2e%10.2e\n", 
				 R, Rtab[0], Rtab[49]);
		cdEXIT(EXIT_FAILURE);
	  }
	/* Si esta, entonces interpolamos la tabla de U(R) para el R que hemos calculado y sacamos el valor de	 U(R)*/
	else
	  {
		lgHit = false;
		j = 1;
		while( !lgHit && j < 49 )
		  {
			if( Rtab[j] <= (realnum)R && Rtab[j+1] > (realnum)R )
			  {
				a = (R - Rtab[j])/(Rtab[j+1] - Rtab[j]);
				U = (1-a) * Utab[j] + a * Utab[j+1];
				lgHit = true;
			  }
			j += 1;
		  }
	  }	
	fprintf(ioQQQ," %.2f  %.2f \n",dense.xIonDense[ipHYDROGEN][1],dense.xIonDense[ipHYDROGEN][0]);
	/* Finalmente calculamos la densidad a partir de todo lo demas*/
	fabden_v =  rho_Rmax * (U_Rmax/U) * POW2(Rmax/R);
  }	

  return(fabden_v);

}


double dense_fabden(double radius, 
					double depth)
{
  double fabden_v;
  double fluc;

#	ifdef DEBUG_FUN
  fputs( "<+>fabden()\n", debug_fp );
#	endif

  /* this routine is called when the DLAW command is given,
   * and must return the hydrogen density at the current radius or depth
   * RADIUS is the radius, the distance from the center of symmetry.
   * DEPTH is the depth into the cloud, from the illuminated face
   * both are in cm
   *
   * radius, depth, and the array DensityLaw are double precision, although
   * FABDEN itself is not
   *
   * this is one way to generate a density
   * fabden = radius * depth
   *
   * this is how to use the parameters on the dlaw command
   * fabden = DensityLaw(1) + radius * DensityLaw(2)
   *
   * following must be removed if this sub is to be used */
  fabden_v = 0.0;
  /*    fprintf(ioQQQ,"TTT 2%.2f\n",fabden_v);*/

  if( depth > radius )
	TotalInsanity();


  switch((int)(dense.DensityLaw[0]))
	{
	case 0:
	  /* dens = P1 */
	  fprintf(ioQQQ," %.2f  %.2f \n",timesc.sound_speed_adiabatic,timesc.sound_speed_isothermal);
	  fabden_v = dense.DensityLaw[1];
	  break;
	case 1:
	  /* dens = P1 * exp(-(Depth/P2)^2) */
	  fabden_v = dense.DensityLaw[1]*
		exp(-(pow(depth/pow(10.,dense.DensityLaw[2]),2)));
	  break;
	case 2:
	  /* dens = 10^P1 * exp((R*sin(P2)-P3)/P4) */
	  fabden_v = pow(10.,dense.DensityLaw[1])*
		exp((radius*sin(dense.DensityLaw[2])-pow(10.,dense.DensityLaw[3]))/
			pow(10.,dense.DensityLaw[4]));
	  break;
	case 21:
	  /* dens = 10^P1 * MIN[exp((R*sin(P2)-P3)/P4),1.0] */
	  fabden_v = pow(10.,dense.DensityLaw[1])*
		MIN2(exp((radius*sin(dense.DensityLaw[2])-pow(10.,dense.DensityLaw[3]))/
				 pow(10.,dense.DensityLaw[4])),1.0);
	  break;
	case 22:
	  fabden_v = MAX2(pow(10.,dense.DensityLaw[1])*
					  MIN2(exp((radius*sin(dense.DensityLaw[2])-pow(10.,dense.DensityLaw[3]))/
							   pow(10.,dense.DensityLaw[4])),1.0),pow(10.,dense.DensityLaw[5]));
	  break;
	case 3:
	  /* dens = P1 + P2*exp(-((depth-P3)/P4)^2) + P5*exp(-((depth-P6)/P7)^2   */
	  fabden_v = dense.DensityLaw[1]+dense.DensityLaw[2]*
		exp(-(pow((depth-pow(10.,dense.DensityLaw[3]))/pow(10.,dense.DensityLaw[4]),2)))+
		dense.DensityLaw[5]*
		exp(-(pow((depth-pow(10.,dense.DensityLaw[6]))/pow(10.,dense.DensityLaw[7]),2)));
	  break;
	case 31:
	  /* dens = P1 + P2*exp(-((depth-P3)/P4)^2) + P5*exp(-((depth-P6)/P7)^2   */
	  /* if depth LT P3 then dens = dens + P8 */
	  fabden_v = dense.DensityLaw[1]+dense.DensityLaw[2]*
		exp(-(pow((depth-pow(10.,dense.DensityLaw[3]))/pow(10.,dense.DensityLaw[4]),2)))+
		dense.DensityLaw[5]*
		exp(-(pow((depth-pow(10.,dense.DensityLaw[6]))/pow(10.,dense.DensityLaw[7]),2)));
	  if (depth < pow(10.,dense.DensityLaw[3])) fabden_v = fabden_v + dense.DensityLaw[8];
	  break;
	case 32:
	  /* dens = P2*exp(-((depth-P3)/P4)^2) + P5*exp(-((depth-P6)/P7)^2   */
	  /* dens = P1+dens*(1.-exp(P9*(sin(10.^depth/10.^p.P8)))/exp(P9)) */
	  fabden_v = dense.DensityLaw[1]+dense.DensityLaw[2]*
		exp(-(pow((depth-pow(10.,dense.DensityLaw[3]))/pow(10.,dense.DensityLaw[4]),2)))+
		dense.DensityLaw[5]*
		exp(-(pow((depth-pow(10.,dense.DensityLaw[6]))/pow(10.,dense.DensityLaw[7]),2)));
	  fabden_v = fabden_v * exp(dense.DensityLaw[9]*sin((depth/pow(10.,dense.DensityLaw[8]))-pow(10.,dense.DensityLaw[10])))/
		exp(dense.DensityLaw[9]);
	  fabden_v = fabden_v + dense.DensityLaw[1]/4.;
	  break;
	case 33:
	  /* fluc = 1/(1+exp(10*sin(depth/P8)+P9))*/
	  /* dens = P1/4.+fluc*(P1+P2*exp(-((depth-P3)/P4)^2) + P5*exp(-((depth-P6)/P7)^2))   */
	  fabden_v = dense.DensityLaw[1]+dense.DensityLaw[2]*
		exp(-(pow((depth-pow(10.,dense.DensityLaw[3]))/pow(10.,dense.DensityLaw[4]),2)))+
		dense.DensityLaw[5]*
		exp(-(pow((depth-pow(10.,dense.DensityLaw[6]))/pow(10.,dense.DensityLaw[7]),2)));
	  fluc = 1./(1.+exp(10.*sin(depth/pow(10.,dense.DensityLaw[8]))+dense.DensityLaw[9]));
	  fabden_v = dense.DensityLaw[1]/4. + fabden_v * fluc;
	  break;
	case 34:
	  /* fluc = 1/(1+exp(P3*sin(depth/P4)+P5))*/
	  /* dens = P1/2.+fluc*(P1*2.+P2*6.3*exp(-((depth-P2)/P2)^2) + P2*1000.*exp(-((depth-P2)/P2)^2))   */
	  fabden_v = dense.DensityLaw[1]*2.0+dense.DensityLaw[1]*6.3*
		exp(-(pow((depth-pow(10.,dense.DensityLaw[2]+0.2))/pow(10.,dense.DensityLaw[2]+0.05),2)))+
		dense.DensityLaw[1]*10.*
		exp(-(pow((depth-pow(10.,dense.DensityLaw[2]+1.))/pow(10.,dense.DensityLaw[2]+0.55),2)));
	  fluc = 1./(1.+exp(dense.DensityLaw[3]*sin(depth/pow(10.,dense.DensityLaw[4])+dense.DensityLaw[7])
						+dense.DensityLaw[5]));
	  fabden_v = dense.DensityLaw[6] + fabden_v * fluc;
	  break;
	case 35:
	  /* fluc = 1/(1+exp(P3*sin(depth/P4+P6)+P5))*/
	  /* dens = P1*fluc + P2   */
	  fabden_v = dense.DensityLaw[1];
	  fluc = 1./(1.+exp(dense.DensityLaw[3]*sin(depth/pow(10.,dense.DensityLaw[4])+dense.DensityLaw[6])
						+dense.DensityLaw[5]));
	  fabden_v = dense.DensityLaw[2] + fabden_v * fluc;
	  break;
	case 302:
	  /* dens = 10^P1 + 10^P2*(10^P3-depth)^(-2) + 10^P4/T4^P5 ut	*/
	  fabden_v = pow(10.,dense.DensityLaw[1])+
		pow(10.,dense.DensityLaw[2])*pow((pow(10.,dense.DensityLaw[3])-depth)/pow(10.,dense.DensityLaw[4]),-2.)+
		pow(10.,dense.DensityLaw[5])/pow(phycon.te/1e4,dense.DensityLaw[6]);
	  break;
	case 4:
	  /* dens = 10.^P1 * exp(-depth / 10.^P2) */ 
	  fabden_v = dense.DensityLaw[1] * exp(-radius/pow(10.,dense.DensityLaw[2]));
	  break;
	case 41:
	  /* dens = 10.^P1 + 10^P2*((10^P3-depth)/10^P4)^(P5) */ 
	  fabden_v = pow(10.,dense.DensityLaw[1])+ 
		pow(10.,dense.DensityLaw[2])*pow((pow(10.,dense.DensityLaw[3])-depth)/pow(10.,dense.DensityLaw[4]),-dense.DensityLaw[5]);
	  break;	
	case 5:
	  fprintf(ioQQQ,"1 y 0 %.2f  %.2f \n",dense.xIonDense[ipHYDROGEN][1],dense.xIonDense[ipHYDROGEN][0]);
	  /* dens = 10.^3 * x */ 
	  if (nzone > 2) fabden_v = (pow(10.,dense.DensityLaw[1]) * dense.xIonDense[ipHYDROGEN][1]/(dense.xIonDense[ipHYDROGEN][1]+dense.xIonDense[ipHYDROGEN][0]));
	  else fabden_v = (pow(10.,dense.DensityLaw[1]));
	  break;
	case 51:
	  fprintf(ioQQQ,"1 y 0 %.2f  %.2f \n",dense.xIonDense[ipHYDROGEN][1],dense.xIonDense[ipHYDROGEN][0]);
	  if (nzone > 2)
		if ((dense.xIonDense[ipHYDROGEN][1]/(dense.xIonDense[ipHYDROGEN][1]+dense.xIonDense[ipHYDROGEN][0])) < 1) fabden_v = pow(phycon.te,2.);
		else fabden_v = pow(10.,dense.DensityLaw[1])+ 
			   pow(10.,dense.DensityLaw[2])*pow((pow(10.,dense.DensityLaw[3])-depth)/pow(10.,dense.DensityLaw[4]),-dense.DensityLaw[5]);
	  else fabden_v = (pow(10.,dense.DensityLaw[1]));
	  break;	  
	case 52:
	  fabden_v = dense_proplyd_LinealInterpol(depth,  dense.DensityLaw[1], dense.DensityLaw[2], dense.DensityLaw[3]); /* refactored into a separate function */
	  break;		  
	case 53:
	  fabden_v = dense_proplyd(depth,  dense.DensityLaw[1], dense.DensityLaw[2], dense.DensityLaw[3]); /* refactored into a separate function */
	  break;	  
	}

#	ifdef DEBUG_FUN
  fputs( " <->fabden()\n", debug_fp );
#	endif
  return( fabden_v );

}





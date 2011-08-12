#include <math.h>
#include <stdio.h>


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



double dense_proplyd(double depth, double r0, double Rmax, double rho_Rmax, double x0, double A, double ifrac_m)
{
  
  /* New variables for interpolating velocity of globule flow 14 Jun 2011 */ 
  double R; /* Dimensionless globule radius */
  double U; /* Dimensionless velocity */
  static double U_Rmax; /* Dimensionless velocity at the beginnig of model */
  static double rho_m; /* Density at the sonic point, R=1 */
  static double c_m; /* Isothermal sound speed at sonic point, R=1*/
  static double dR; /* Dimension of the Ioniation Front*/
  static double sigma=6.3e-18; 
  double fabden_v; /* The density */
  double c_s; /* Calculated sound speed */
  double ifrac; /* Ionization fraction X(R) */
  double x; /* ifrac argument x(R)*/
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
  dR = A/(r0*sigma*rho_Rmax);

  /* First time called, so set up spline interpolation for U vs R using the tables Utab and Rtab */
  if (nzone == 0) spline(Rtab, Utab, 50, 2e31, 2e31, Uspline);

  if (nzone == 0) {
	U_Rmax = U;
	fabden_v = rho_Rmax;
  }	else {
	if (R > 1) {
	  // In outer, supersonic part
	  // Find the velocity by spline interpolation
	  splint(Rtab, Utab, Uspline, 50, R, &U);
	} else {
	  // In inner, subsonic part
	  x = x0 * (R - dR -1.0)/dR;
	  ifrac = 0.5 * (tanh(x) + 1.0);
	  ifrac_m = 0.5 * (tanh(x0) + 1.0);
	  U = (1.0 - sqrt(1.0-((1.0+ifrac)/(1.0+ifrac_m))));
	} 
	/* Calculate density by continuity */
	fabden_v = rho_Rmax * (U_Rmax/U) * POW2(Rmax/R) ;
  }	

  return(fabden_v);

}


/* This has to be outside so that all the routines can see it */
static int nzone;
const int NZONES = 300;


int main()
{
  double densidad[NZONES], depth[NZONES], r0, Rmax, rho_Rmax, x0, A, ifrac_m;

  depth = 1.;
  r0 = 8.e14;
  Rmax=9.;
  rho_Rmax=3.8e6;
  x0=2.3;
  A=10.;
  ifrac_m=0.9;

  
  for (nzone=0; nzone < NZONES; nzone++) {	
	depth[nzone] = depth_min + (depth_max - depth_min)*double(nzone)/double(NZONES - 1)
	densidad[nzone] = dense_proplyd(depth[nzone], r0, Rmax, rho_Rmax, x0, A, ifrac_m);
  }	
  
  printf("%f\n",densidad);
}


#ifndef _INCLUDE_AIRTOICERAYTRACING_H_
#define _INCLUDE_AIRTOICERAYTRACING_H_

#include <string> 
#include <sstream>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <vector>
#include <cstdlib>
#include <chrono>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_spline.h>

#include <sys/time.h>

namespace AirToIceRayTracing{
  
  //static constexpr double pi=4.0*atan(1.0); /**< Gives back value of Pi */
  static const double pi=3.1415927; /**< Gives back value of Pi */
  static const double spedc=299792458.0; /**< Speed of Light in m/s */

  ////Define std::vectors to store data from the file
  static std::vector <std::vector <double> > nh_data;////n(h) refractive index profile of the atmosphere as a function of height
  static std::vector <std::vector <double> > lognh_data;////log(n(h)-1) log of the refractive index profile of the atmosphere as a function of height subtracted by 1
  static std::vector <std::vector <double> > h_data;////height data
  
  ////Define Arrays for storing values of ATMLAY and a,b and c parameters taken from the Atmosphere.dat file
  static double ATMLAY[5];
  static double abc[5][3];

  ////define dummy variables which will be filled in later after fitting
  static double C_air[5];
  static double B_air[5];

  ////define variables which are going to be used by GSL for linear interpolation
  static gsl_interp_accel * accelerator;
  static gsl_spline *spline;
  
  ////The variable which will store the max layers available in an atmosphere model
  static int MaxLayers=0;
  
  ////This Function reads in the values of ATMLAY and a,b and c parameters taken from the GDAS Atmosphere.dat file. The a,b and c values are mass overburden values and are not required in this code.
  int readATMpar(std::string atmosFileName);

  ////This Function reads in the tavulated refractive index profile from the GDAS Atmosphere.dat file and fills in the nh_data, lognh_data and h_data vectors
  int readnhFromFile(std::string atmosFileName);

  static bool UseConstantRefractiveIndex=false;
  
  ////Set the value of the asymptotic parameter of the ice refractive index model
  static const double A_ice=1.78;

  ////Get the value of the B parameter for the ice refractive index model
  double GetB_ice(double z);

  ////Get the value of the C parameter for the ice refractive index model
  double GetC_ice(double z);

  ////Get the value of refractive index model for a given depth in ice
  double Getnz_ice(double z);

  ////Set the value of the asymptotic parameter of the air refractive index model
  static double A_air=1.00;

  ////Set the value of the asymptotic parameter of the air refractive index model
  static double A_const=1.00;

  int FillInAirRefractiveIndex();
  
  ////Get the value of the B parameter for the air refractive index model
  double GetB_air(double z);

  ////Get the value of the C parameter for the air refractive index model
  double GetC_air(double z);

  ////Get the value of refractive index model for a given depth in air
  double Getnz_air(double z);

  /* E-feild Fresnel coefficient for S-polarised wave which is perpendicular to the plane of propogation/incidence. This function gives you back the reflection coefficient. The transmittance is t=1+r */
  double Refl_S(double thetai, double IceLayerHeight);
  double Trans_S(double thetai, double IceLayerHeight);
  
  /* E-feild Fresnel coefficient for P-polarised wave which is parallel to the plane of propogation/incidence. This function gives you back the reflection coefficient. The transmittance is t=(n_1/n_2)*(1+R) */
  double Refl_P(double thetai, double IceLayerHeight);
  double Trans_P(double thetai, double IceLayerHeight);
  
  ////Use GSL minimiser which uses Brent's Method to find root for a given function
  double FindFunctionRoot(gsl_function F,double x_lo, double x_hi,const gsl_root_fsolver_type *T,double tolerance, int iterations);

  /////Functions used for Raytracing in Ice using the analytical solution

  ////Analytical solution describing the ray path in ice
  struct fDnfR_params { double a, b, c, l; };
  double fDnfR(double x,void *params);

  // ////Define the function that will be minimised to calculate the angle of reciept (from the vertical) on the antenna and the hit point of the ray on the ice surface given a ray incident angle
  // struct fdxdz_params { double lang, z0,z1; int airorice;};
  // double fdxdz(double x,void *params);

  ////The function used to calculate ray propogation time in ice
  struct ftimeD_params { double a, b, c, speedc,l; int airorice; };
  double ftimeD(double x,void *params);
  double fpathD(double x,void *params);
  
  ////Get the distance on the 2ndLayer at which the ray should hit given an incident angle such that it hits an target at depth of z0 m in the second layer.
  //// n_layer1 is the refractive index value of the previous layer at the boundary of the two mediums
  //// A,B and C values are the values required for n(z)=A+B*exp(Cz) for the second layer
  //// TxDepth is the starting height or depth
  //// RxDepth is the final height or depth
  //// AirOrIce variable is used to determine whether we are working in air or ice as that sets the range for the GSL root finder.

  ////These two functions are used in GetLayerHitPointPar() which is the main function for doing the minimisation
  double GetRayHorizontalPath(double A, double RxDepth, double TxDepth, double Lvalue, int AirOrIce);
  double GetRayPropagationTime(double A, double RxDepth, double TxDepth, double Lvalue, int AirOrIce);
  double GetRayGeometricPath(double A, double RxDepth, double TxDepth, double Lvalue, int AirOrIce);

  double *GetLayerHitPointPar(double n_layer1, double RxDepth,double TxDepth, double IncidentAng, int AirOrIce);
  
  ////This function flattens out 2d std::vectors into 1d std::vectors
  std::vector<double> flatten(const std::vector<std::vector<double> >& v);

  ////This function loads in the GDAS atmosphere file. It calls the other functions to load in the tabulated refractive index values and the sea level refractive index value from the file. It also reads the mass overburden A,B and C values from the file
  int MakeAtmosphere(std::string atmosFileName);
  
  ////Get Propogation parameters for ray propagating in air
  double * GetAirPropagationPar(double LaunchAngle, double AirTxHeight, double IceLayerHeight);
  
  ////Get Propogation parameters for ray propagating in ice
  double * GetIcePropagationPar(double IncidentAngleonIce, double IceLayerHeight, double AntennaDepth, double Lvalue);
  
  ////This function will be minimised to get the launch angle in air. It looks at the total horizontal distance and tries to find a ray in air and ice which can cover that horizontal distance perfectly. It takes into account the discontinuity in refractive index at the air ice boundary
  struct MinforLAng_params { double airtxheight, icelayerheight, antennadepth, horizontaldistance; };
  double MinimizeforLaunchAngle(double x, void *params); 
 
  ////This function uses my raw code to calculate values for CoREAS. Since its directly using the minimiser to calculate launch angles and distances it is slightly slower than its _Table version.  
  bool GetRayTracingSolution(double SrcHeightASL, double HorizontalDistanceToRx,double RxDepthBelowIceBoundary,double IceLayerHeight, double& opticalPathLengthInIce, double& opticalPathLengthInAir, double& geometricalPathLengthInIce, double& geometricalPathLengthInAir, double& launchAngle, double& horizontalDistanceToIntersectionPoint,  double& AngleOfIncidenceOnIce, double &RecievedAngleInIce);
  
  void Air2IceRayTracing(double AirTxHeight, double HorizontalDistance, double IceLayerHeight,double AntennaDepth, double StraightAngle, double dummy[20]);
  
}
#endif

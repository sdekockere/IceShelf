/*/////////////////////////////////////////////////////////////////////////////////////
This header defines some ice density models
*//////////////////////////////////////////////////////////////////////////////////////

#ifndef ICEDENSITYMODELS_h
#define ICEDENSITYMODELS_h

#include "G4Types.hh"
#include "IceRayTracing.hh"

namespace IceDensityModels{

// The exponential model based on the paper "In situ index-of-refraction measurements of the South Polar firn with the RICE detector" at South Pole
// NOTE: depth here is in units length and is a positive value (axis points downwards)

inline G4double rhoSouthPole(G4double depth){
    G4double rho_ice = 917.*kg/m3;
    G4double rho_surface = 359.*kg/m3;
    G4double t_firn = 100.*m;
    return rho_ice - (rho_ice - rho_surface)*std::exp(-1.9*depth/t_firn);
}

// NOTE: depth here is in units length and is a positive value (axis points downwards)
inline G4double rhoTaylorDome(G4double depth){
    return 0.46*g/cm3 + 0.468*g/cm3*(1.-std::exp(-0.02*(1./m)*depth));
}

inline G4double rhoConst(){
    //return 400.*kg/m3;
    return 359.*kg/m3;
}

// The double exponential profile for Greenland (https://arxiv.org/pdf/1805.12576.pdf)
// NOTE: depth here is in units length and is a positive value (axis points downwards)
inline G4double rhoGreenland(G4double depth){
    G4double rho_val = -999.;
    if(depth <= IceRayTracing::TransitionBoundary*m){
        rho_val = (0.917*g/cm3) - (0.594*g/cm3)*std::exp(-1.*depth/(30.8*m));
    } else {
        rho_val = (0.917*g/cm3) - (0.367*g/cm3)*std::exp(-1.*(depth - IceRayTracing::TransitionBoundary*m)/(40.5*m));
    }
    return rho_val;
}

// Here you can choose the model you want to use
// NOTE: depth here is in units length and is a positive value (axis points downwards)
inline G4double rho(G4double depth){
    //return rhoTaylorDome(depth);    
    return rhoSouthPole(depth);    
    //return rhoConst();
    //return rhoGreenland(depth);
}

}

#endif

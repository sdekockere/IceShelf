/*//////////////////////////////////////////////////////////////////////////////////////////////////
This header file contains some global variables.
*///////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef ICESHELFGLOBALVARIABLES_h
#define ICESHELFGLOBALVARIABLES_h

#include "G4Types.hh"
#include "G4SystemOfUnits.hh"
#include "IceRayTracing.hh"

// Set on true if you want to include ray tracing (highly recommended, but slows down simulation)
static const G4bool rayTracingOn = true;
// The offset used in case a particle has a depth of 0, which the ray tracing cannot handle
static const G4double zeroDepthOffset = 0.00001*m;
// The expected depth of the center of the shower, < 0 if under the ice surface
// Needed to create the tables for the ray tracing
static const G4double centerShowerDepth = -10.*m;
// Set true if you want to include the direct rays
static const G4bool includeDirectRays = true;
// Set true if you want to include the indirect rays
static const G4bool includeIndirectRays = true;

// Set true if you want to include air-to-ice ray tracing, meaning that points with depth > 0 will be 
// treated as being in air during ray tracing. If not enabled, ray tracing will fail for points with 
// depth > 0, meaning the corresponding step will not be taken into account for emission calculation.
// If ray tracing is off, this parameter does not influence the emission calculation: particles with
// depth > 0 will automatically be treated as if IN ICE.
//
// This feature was originally introduced to calculate transition radiation explicitly, but in the end
// we decided this was not needed. Switching this on will enable you to use the existing functions in 
// IceShelfSteppingAction to calculate air-to-air step contributions, but note that this also means
// air-to-ice and ice-to-air steps will start to contribute to the emission, in case a particle goes 
// from depth  > 0 to depth < 0 (or vice versa) in a single step. This is NOT how the end-point 
// formalism should be used, a single step should always be in the same medium. It is therefor adviced
// that when setting includeAirToIceRT to true, the code is adapted to make sure boundary-crossing steps
// (depth < 0 to depth > 0 in a single step, or vice versa) are ruled out.
static const G4bool includeAirToIceRT = false;

// Set true if you want to have a separate trace file for the sudden appearance radiation. This only
// means the sudden appearance radiation will be stored in a separate file too, it will still be 
// included in the direct and indirect emission trace files!
// Setting to false means this extra trace file will not be generated, but as stated above the sudden 
// appearance radiation will still be included in the direct and indirect emission trace files.
static const G4bool writeSARadSeparateToo = true;

// The limits on the focusing factor, calculated during ray tracing
static const G4double focFactUnderlim = 0.5;
static const G4double focFactUpperlim = 2.0;

// Dimensions of the simulated ice shelf
static const G4double shelfSizeX = 40.*m;
static const G4double shelfSizeY = 20.*m;
static const G4double shelfSizeZ = 40.*m;

// Number of layers of ice, used for the density gradient.
static const G4int numberOfLayers = 2000;

// Index of refraction of ice, in case ray tracing is disabled
//static const G4double nIce = 1.52; // For the PRD paper
//static const G4double nIce = 1.78; // For deep south pole ice
//static const G4double nIce = 1.4;
static const G4double nIce = 1.2; // To mimic the cherenkov ring we get from ray tracing using Uzair's
                                  // analytical model, at 150m down in the ice.

// Sampling rate of antenna
static const G4double sPeriod = 0.2*ns;

// The length of the recorded antenna trace in time calculated during the simulation.
// The length of the calculated trace to save is determined by traceLengthToSave
static const G4double traceLength = 40000.*ns;

// The time length of the part of the traces we want to save, 
// beginning from the point of time the first non-zero value corresponds to minus startBuffer
static const G4double startBuffer = 20.*ns;
static const G4double traceLengthToSave = 400.*ns;

// The energy threshold above which the kinetic energy of a charged particle should be 
// to take part in the electrical field calculations
static const G4double kinEnergyThreshold = 0.1*MeV;

// The size of the sides of the blocks for the calculation of energy density slice histos
static const G4double blockSize = 1.*cm;

// The number of bins for the r-axis of the radial energy histograms
// and the shifting value over which the range [0, r_max[ will be shifted in order to create 
// a logarithmic binning (log(0) is not defined...). After creating the binning all the bin limits will
// be shifted over the same value back again, so that the original [0, r_max[ range is restored.
// Note this value determines the binwidth dX0 of the first bin.
static const G4int nBinsRadialHist = 1000;
//static const G4double shiftingValRadialHist = 37.5115*cm; // Gives dX0 = 0.1 cm
                                                        // if number of bins is 1000
                                                        // and r_max is 500 cm.
static const G4double shiftingValRadialHist = 22.089741*cm; // Gives dX0 = 0.1 cm
                                                        // if number of bins is 1000
                                                        // and r_max is 2000 cm.

// The maximum depth value for the depth histograms and the number of bins
// as well as the kinetic energy cuts applied and the number of bins and under- and upper limits
// for the energy distribution on the y-axis
// Particles falling below the kinetic energy cut will not contribute to the histograms.
// NOTE: The bin width represents the depth intervals between two consecutive depth points 
// for which the total amount of particles will be given. More number of bins means more detail, 
// but also longer computation time.
static const G4double maxDepthValue = 20.*m;
static const G4int depthNumberOfBins = 2000;
static const G4double depthKinELimVals[4] = {0.299*GeV, 0.299*GeV, 0.00299*GeV, 0.00299*GeV};  
                                                       // Kinetic energy cut 
                                                       // for hadrons (without pi0's) (i = 0), 
                                                       // muons (i = 1), 
                                                       // electrons (i = 2) and 
                                                       // photons (including pi0's) (i = 3)
static const G4double depthEnergyNumberOfBins = 50;
static const G4double depthEnergyUnderLim = 0.003*GeV;
static const G4double depthEnergyUpperLim = 1*GeV; 

// The maximum time value for the charge in function of time histograms, the number of bins,
// the energy limit above which the kinetic energy of the particle needs to be to contribute 
// to the charge in function of time histograms
// and the number of bins and under- and upper limits
// for the energy distribution on the y-axis
// NOTE: The bin width represents the time intervals between two consecutive time points 
// on which the total amount of charge will be given. More number of bins means more detail, 
// but also longer computation time.
static const G4double maxTimeValue = 200.*ns;
static const G4int timeNumberOfBins = 2000;
static const G4double chargeKinELim = kinEnergyThreshold; // Don't forget to add the unit if you're
                                                          // putting an other value than
                                                          // kinEnergyThreshold here!
static const G4double chargeEnergyNumberOfBins = 50;
static const G4double chargeEnergyUnderLim = 0.003*GeV;
static const G4double chargeEnergyUpperLim = 1*GeV; 

// The cascade front depth values for which a snapshot will be made during the simulation
// in units mass/area,
// together with a kinetic energy cut above which the particles need to be to contribute.
//static const G4double snapshotDepthsMassArea[0] = {};
static const G4double snapshotDepthsMassArea[1] = {350.*g/cm2};
//static const G4double snapshotDepthsMassArea[18] 
//                    = {50.*g/cm2, 100.*g/cm2, 150.*g/cm2, 200.*g/cm2, 250.*g/cm2, 300.*g/cm2, 350.*g/cm2, 400.*g/cm2, 450.*g/cm2, 500.*g/cm2, 550.*g/cm2, 600.*g/cm2, 650.*g/cm2, 700.*g/cm2, 750.*g/cm2, 800.*g/cm2, 850.*g/cm2, 900.*g/cm2};
static const G4int numberOfSnapshotDepthsMassArea = sizeof(snapshotDepthsMassArea)/sizeof(*snapshotDepthsMassArea);
static const G4double snapshotKinELim = 0.; // Don't forget to add the unit if you're
                                                            // putting an other value than
                                                            // kinEnergyThreshold here!

// The number pi
static const G4double pi = 3.14159265359;

// The speed of light in vacuum
static const G4double c0 = 299792458.*m/s;

// The charge of a positron in franklin (unit charge of the electrostatic CGS unit system)
// As the unit of charge in geant4 is the charge of the positron (charge e+ = 1.),
// you should multiply all charges by this number to get it in franklin
static const G4double e_CGS = 4.80319946105e-10;

// The unit of charge in the electrostatic CGS unit system, expressed in geant4 internal units
static const G4double franklin = eplus/e_CGS;

// The conversion factor for getting electrical field in V/m instead of electrostatic CGS units.
// Unit of voltage in electrostatic CGS unit system is statvolt. 1 statvolt = 299.792458 volts.
// So 1 statvolt/cm = 29979.2458 volt/m. So to go from electrostatic CGS units to SI units you should 
// multiply by:
static const G4double CGSToSIFactor = 29979.2458;

// The threshold for the doppler factor. If the factor is smaller than this value we need to use
// an approximation method in the endpoint formalism
static const G4double approxThreshold = 0.9e-1;
// Threshold used when calculating the transition radiation (not using the interpoint fallback)
//static const G4double approxThresholdTrans = 1e-2;

#endif

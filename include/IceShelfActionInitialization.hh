/*///////////////////////////////////////////////////////////////////////////////////////////////////
The ActionInitialization class initializes the user action classes.
In the case of this project there are 4 action classes:
    - IceShelfPrimaryGeneratorAction (mandatory): deals with the primary particles
    - IceShelfRunAction: tells the program what to do before and after a single run
    - IceShelfEventAction: tells the program what to do before and after a single event
    - IceShelfSteppingAction: tells the program what to do at each end of stepping process
      -> In a step all the transient information is being calculated (e.g. energy deposit in the
         detector). If you want to access this information you have implement a SteppingAction class.
The Build() function is where all the initialization happens.
The private field fOutoutFileName contains the name of the output file. The private field fPrimaryList contains a pointer to the primary particle list, which is passed on to the IceShelfPrimaryGeneratorAction object. The private field fDetConstruction holds a pointer to the DetectorConstruction instance, which is passed on to the SteppingAction instance.
*////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef ICESHELFACTIONINITIALIZATION_h
#define ICESHELFACTIONINITIALIZATION_h

#include "G4VUserActionInitialization.hh"
#include "globals.hh"
#include <vector>

class IceShelfDetectorConstruction;

// Action initialization class

class IceShelfActionInitialization : public G4VUserActionInitialization
{
    public:
        IceShelfActionInitialization(G4String, IceShelfDetectorConstruction*, std::vector<std::vector<G4double>>*, G4double*, G4double, G4double, G4String, G4String, G4String);
        virtual ~IceShelfActionInitialization();

        virtual void Build() const;

    private:
        G4String fOutputFileName;
        IceShelfDetectorConstruction* fDetConstruction;
        std::vector<std::vector<G4double>>* fPrimaryList;
        G4double* fSnapshotTimes;
        G4double fZenithAngleRad;
        G4double fAzimuthAngleRad;
        G4String fReasFileName;
        G4String fListFileName;
        G4String fAtmosphereFileName;
};

#endif

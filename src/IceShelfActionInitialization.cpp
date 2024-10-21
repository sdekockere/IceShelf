#include "IceShelfActionInitialization.hh"
#include "IceShelfPrimaryGeneratorAction.hh"
#include "IceShelfRunAction.hh"
#include "IceShelfEventAction.hh"
#include "IceShelfSteppingAction.hh"
#include "IceShelfDetectorConstruction.hh"

IceShelfActionInitialization::IceShelfActionInitialization
     (G4String outputFileName, IceShelfDetectorConstruction* detConstruction, std::vector<std::vector<G4double>>* primaryList, G4double* snapshotTimes, G4double zenithAngleRad, G4double azimuthAngleRad,
                         G4String reasFileName, G4String listFileName, G4String atmosphereFileName)
 : G4VUserActionInitialization(),
   fOutputFileName(outputFileName),
   fDetConstruction(detConstruction),
   fPrimaryList(primaryList),
   fSnapshotTimes(snapshotTimes),
   fZenithAngleRad(zenithAngleRad),
   fAzimuthAngleRad(azimuthAngleRad),
   fReasFileName(reasFileName),
   fListFileName(listFileName),
   fAtmosphereFileName(atmosphereFileName)
{}

IceShelfActionInitialization::~IceShelfActionInitialization()
{}

void IceShelfActionInitialization::Build() const
{
    // Setting the RunAction
    auto runAction = new IceShelfRunAction(fOutputFileName, fSnapshotTimes, fZenithAngleRad,
                                    fAzimuthAngleRad, fReasFileName, fListFileName, fAtmosphereFileName);
    SetUserAction(runAction);

    // Setting the EventAction
    auto eventAction = new IceShelfEventAction(runAction);
    SetUserAction(eventAction);

    // Setting the SteppingAction
    auto steppingAction = new IceShelfSteppingAction(fDetConstruction, runAction, eventAction, 
                                                  fSnapshotTimes, fZenithAngleRad);
    SetUserAction(steppingAction);

    // Linking the SteppinAction to the EventAction
    eventAction->SetSteppingAction(steppingAction);

    // Setting the PrimaryGeneratorAction. It needs a pointer to the EventAction instance so the 
    // PrimaryGeneratorAction object can tell it the correct weight of the primary.
    SetUserAction(new IceShelfPrimaryGeneratorAction(fPrimaryList, eventAction));
}

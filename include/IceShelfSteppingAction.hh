/*///////////////////////////////////////////////////////////////////////////////////////////////////
This is the SteppingAction class. It tells the program what to do at each end of a stepping process.
In a step all the transient information is being calculated (e.g. energy deposited in the detector).
If you want to access this information you have to implement a SteppingAction class.
The things that need to be done after each stepping process is implemented in the UserSteppingAction() 
function. 
The fDetConstruction field contains a pointer to the DetectorConstruction instance. The fEventAction field contains the pointer to the EventAction object. This is needed to pass down the quantities of interest to the EventAction object.
The fTimeBar field represents a time bar. It has the same resolution as the histogram keeping track of
charge in function of time (see the RunAction class). It will represent through 0's and 1's at which point in time the current particle is exactly. This infomration is needed in order ot correctly fill up the charge in function of time histogram. fTrackIDs will keep track of the particles that already have
contributed to the charge-time histogram.
*////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef ICESHELFSTEPPINGACTION_h
#define ICESHELFSTEPPINGACTION_h

#include "G4UserSteppingAction.hh"
#include "G4Types.hh"
#include "G4ThreeVector.hh"

#include <vector>
#include "globals.hh"

#include "IceShelfGlobalVariables.hh"

class IceShelfDetectorConstruction;
class IceShelfEventAction;
class IceShelfRunAction;

class IceShelfSteppingAction : public G4UserSteppingAction
{
    public:
        IceShelfSteppingAction(IceShelfDetectorConstruction*, IceShelfRunAction*, 
                                    IceShelfEventAction*, G4double*, G4double);
        virtual ~IceShelfSteppingAction();
        
        // Calculates the radio emission emitted by energetic charged particles
        // using the endpoint formalism
        void CalcRadioEmissionEndpoint(const G4Step* step);

        virtual void UserSteppingAction(const G4Step* step);

        void ResetCurrentTrackID();
        void ResetCurrentParticleDepthIndex();
        void ResetCurrentParticleChargeIndex();
        void ResetCurrentParticleSnapshotTimeIndex();
        // Transforms Geant4 coordinates (xCoord, yCoord, zCoord) to the coordinates
        // of the shower axis system (y-axis aligned with shower axis pointing to where 
        // the shower came from, origin on surface of the ice block)
        void GetShowerAxisCoordinates(G4double xCoord, G4double yCoord, G4double zCoord,
                                            G4double* axisCoords);
        // ************** DEPRECATED ***************
        // Rotates vect counter clockwise over angle theta around axis.
        // Counter clockwise is defines wrt observer which has axis pointing towards it
        // G4ThreeVector RotateVector(G4ThreeVector vect, G4ThreeVector axis, G4double theta);
        // *****************************************
        // Faster version for rotations (does not evaluate cos and sin function)
        // Works perfect if you e.g. calculate the cos(rot_angle) through scalar product
        // Rotates vect counter clockwise over angle theta around axis.
        // Counter clockwise is defines wrt observer which has axis pointing towards it
        // WARNING: assumes rot_angle <= pi, so that sin(rot_angle) >= 0
        G4ThreeVector RotateVectorCosPi(G4ThreeVector vect, G4ThreeVector axis, G4double cosVal);
        // Returns the values needed to evaluate the end-point formula, NOT using ray tracing
        // See the cpp file for more information
        G4int GetEPVals(G4ThreeVector emitterObsVect,
                        G4double emitterTime, G4ThreeVector beta,
                        G4double &obsTime, G4ThreeVector &launchDir, G4ThreeVector &receiveDir,
                        G4double &doppler, G4double &R);
        // Returns the values needed to evaluate the end-point formula using ray tracing
        // See the cpp file for more information
        G4int GetEPValsRayTracing(unsigned int k,
                G4ThreeVector emitterPos, G4ThreeVector emitterObsVect,
                G4double emitterTime, G4ThreeVector beta, G4int dirOrIndirRay,
                G4double &obsTime, G4ThreeVector &launchDir, G4ThreeVector &receiveDir,
                G4double &doppler, G4double &R, G4double &incidenceAngle, G4double &focusingFact);
        // Calculates the electric field associated with a start point, end point
        // and receiver antenna using the End-Point formalism
        // See the cpp file for more information
        G4int GetEfieldEP(unsigned int k, G4ThreeVector observerPos,
              G4ThreeVector startPoint, G4ThreeVector endPoint, G4double startTime, G4double endTime,
              G4ThreeVector beta, G4double particleCharge, G4int dirOrIndirRay,
              G4ThreeVector &startE, G4ThreeVector &endE, G4int &iStart, G4int &iEnd);
        // Calculates the reflection coeffs for the s-polarisation (phi-polarisation)
        // and p-polarisation (theta-polarisation)
        // See the cpp file for more information
        void GetReflectionCoeffs(G4double incidenceAngle, G4double nIceAB, G4double nAirAB,
                                        G4double &rs, G4double &rp);
        // Calculates the transition coeffs for the s-polarisation (phi-polarisation)
        // and p-polarisation (theta-polarisation)
        // See the cpp file for more information
        void GetTransitionCoeffs(G4double indicenceAngle, G4double nIceAB, G4double nAirAB,
                                        G4double &ts, G4double &tp);
        // Applies the Fresnel coefficients for the given Efield
        // NOTE: Eflied should be the electric field AT THE RECEIVER
        // See cpp file for more information
        void ApplyFresnelCoeff(G4double incidenceAngle, G4ThreeVector receiveDir,
                                        G4ThreeVector &Efield, G4int reflOrTrans);

    private:
        IceShelfDetectorConstruction* fDetConstruction;
        IceShelfRunAction* fRunAction;
        IceShelfEventAction* fEventAction;
        G4int fCurrentTrackID; // The ID of the particle currently being simulated
        G4int fCurrentParticleDepthIndex; // Gives the index of the nex depth histogram bin
                                               // that should be filled if the current particle
                                               // ever reaches the corresponding depth
        G4int fCurrentParticleChargeIndex; // Gives the index of the nex charge histogram bin
                                                // that should be filled if the current particle
                                                // ever reaches the corresponding time
        G4int fCurrentParticleSnapshotTimeIndex; // Gives the index of the next snapshot time 
                                         // for which a
                                         // Ntuple should be filled if the current particle ever
                                         // reaches the corresponding time
        G4double* fSnapshotTimes; // List holding the time values for the snapshots to make
        G4double fZenithAngleRad;
};

// inline functions

inline void IceShelfSteppingAction::ResetCurrentTrackID(){
    fCurrentTrackID = 0;
}

inline void IceShelfSteppingAction::ResetCurrentParticleDepthIndex(){
    fCurrentParticleDepthIndex = 0;
}

inline void IceShelfSteppingAction::ResetCurrentParticleChargeIndex(){
    fCurrentParticleChargeIndex = 0;
}

inline void IceShelfSteppingAction::ResetCurrentParticleSnapshotTimeIndex(){
    fCurrentParticleSnapshotTimeIndex = 0;
}

inline void IceShelfSteppingAction::GetShowerAxisCoordinates(G4double xCoord, G4double yCoord,
                                                 G4double zCoord, G4double* axisCoords){
    axisCoords[0] = xCoord*cos(fZenithAngleRad) + (yCoord - shelfSizeY/2.)*sin(fZenithAngleRad);
    axisCoords[1] = (yCoord - shelfSizeY/2.)*cos(fZenithAngleRad) - xCoord*sin(fZenithAngleRad);
    axisCoords[2] = zCoord;
}

#endif

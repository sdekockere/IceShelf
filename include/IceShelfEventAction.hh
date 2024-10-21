/*/////////////////////////////////////////////////////////////////////////////////////////////////////
This is the EventAction class. It takes care of everything that needs to be done before and after the
simulation of a single event. This is done using the functions BeginOfEventAction() and 
EndOfEventAction().
In this project the EventAction class first private field contains the particle id of the primary of the event, the second the weight of the primary of the event. This particle id and weight are set during the selection of the primary in IceShelfPrimaryGeneratorAction. During construction both values are 0, and the values are reset to 0 by the EndOfEventAction function. The third private field contains the amount of energy that leaked out of the simulated volume into the world volume during the simulation of this event.
*//////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef ICESHELFEVENTACTION_h
#define ICESHELFEVENTACTION_h

#include "G4UserEventAction.hh"
#include "IceShelfSteppingAction.hh"
#include "IceShelfRunAction.hh"
#include "globals.hh"

#include "G4ThreeVector.hh"

class IceShelfEventAction : public G4UserEventAction
{
    public:
        IceShelfEventAction(IceShelfRunAction* runAction);
        virtual ~IceShelfEventAction();

        virtual void BeginOfEventAction(const G4Event* event);
        virtual void EndOfEventAction(const G4Event* event);

        G4int GetEventID();

        void SetPrimaryPartName(G4String primaryPartName);
        G4String GetPrimaryPartName();

        void SetStartingTime(G4double startingTime);
        G4double GetStartingTime();

        void SetPrimaryWeight(G4double primaryWeight);
        G4double GetPrimaryWeight();

        void SetPrimaryCharge(G4double primaryCharge);
        G4double GetPrimaryCharge();

        void SetPrimaryStartPos(G4ThreeVector primaryPos);
        G4ThreeVector GetPrimaryStartPos();

        void SetPrimaryStartMomentum(G4ThreeVector primaryMom);
        G4ThreeVector GetPrimaryStartMomentum();

        void SetPrimaryMass(G4double primaryMass);
        G4double GetPrimaryMass();

        void AddLeakedEnergy(G4double leakedEnergy);
        G4double GetLeakedEnergy();

        void SetSteppingAction(IceShelfSteppingAction* steppingAction);


    private:
        G4int fEventID; // The id of the event
        G4String fPrimaryPartName; // The particle name of the primary particle
        G4double fStartingTime; // The time at which the primary entered the ice
        G4double fPrimaryWeight; // The weight of the primary particle
        G4double fPrimaryCharge; // The charge of the primary particle
        G4ThreeVector fPrimaryStartPos; // The starting position of the primary particle
        G4ThreeVector fPrimaryStartMomentum; // The starting momentum (p*c) of the primary particle
        G4double fPrimaryMass; // The mass (m*c^2) of the primary particle
        G4double fLeakedEnergy; // The amount of energy that leaked out of the simulated volume
        IceShelfSteppingAction* fSteppingAction; // Pointer to the IceShelfSteppingAction object
        IceShelfRunAction* fRunAction; // Pointer to the IceShelfRunAction object
        clock_t startOfEventTime;

        // Old codes used for transition radiation calculation.
        // In the end, we concluded that this should not be taken into account explicitly, and comes
        // out naturally by combining the CoREAS and Geant4 framework.
        //
        // Used to calculate the transition radiation for the primary at the start of the event
        // This is the old implementation, which uses the interpoint when the doppler term falls below
        // the given threshold
        //void CalcTransitionRadiationEndpointOld();
        // Used to calculate the transition radiation for the primary at the start of the event
        // This is the new implementation, which always uses the start- and endpoints.
        // It does not add a contribution if electric fields diverge too much (doppler term too small)
        //void CalcTransitionRadiationEndpoint();

};

// inline functions

inline G4int IceShelfEventAction::GetEventID(){
    return fEventID;
}

inline void IceShelfEventAction::SetPrimaryPartName(G4String primaryPartName){
    fPrimaryPartName = primaryPartName;
}

inline G4String IceShelfEventAction::GetPrimaryPartName(){
    return fPrimaryPartName;
}

inline void IceShelfEventAction::SetStartingTime(G4double startingTime){
    fStartingTime = startingTime;
}

inline G4double IceShelfEventAction::GetStartingTime(){
    return fStartingTime;
}

inline void IceShelfEventAction::SetPrimaryWeight(G4double primaryWeight){
    fPrimaryWeight = primaryWeight;
}

inline G4double IceShelfEventAction::GetPrimaryWeight(){
    return fPrimaryWeight;
}

inline void IceShelfEventAction::SetPrimaryCharge(G4double primaryCharge){
    fPrimaryCharge = primaryCharge;
}

inline G4double IceShelfEventAction::GetPrimaryCharge(){
    return fPrimaryCharge;
}

inline void IceShelfEventAction::SetPrimaryStartPos(G4ThreeVector primaryPos){
    fPrimaryStartPos = primaryPos;
}

inline G4ThreeVector IceShelfEventAction::GetPrimaryStartPos(){
    return fPrimaryStartPos;
}

inline void IceShelfEventAction::SetPrimaryStartMomentum(G4ThreeVector primaryMom){
    fPrimaryStartMomentum = primaryMom;
}

inline G4ThreeVector IceShelfEventAction::GetPrimaryStartMomentum(){
    return fPrimaryStartMomentum;
}

inline void IceShelfEventAction::SetPrimaryMass(G4double primaryMass){
    fPrimaryMass = primaryMass;
}

inline G4double IceShelfEventAction::GetPrimaryMass(){
    return fPrimaryMass;
}

inline void IceShelfEventAction::AddLeakedEnergy(G4double leakedEnergy){
    fLeakedEnergy += leakedEnergy;
}

inline G4double IceShelfEventAction::GetLeakedEnergy(){
    return fLeakedEnergy;
}

inline void IceShelfEventAction::SetSteppingAction(IceShelfSteppingAction* steppingAction){
    fSteppingAction = steppingAction;
}

#endif

/*///////////////////////////////////////////////////////////////////////////////////////////////
This is the physics list, which determines the physical processes to take into account during the 
simulation. It is a copy of the SimplePhysicsList found in the slac_rf project send to me by Steven.
In the SimplePhysicsList.cpp file you can see that the folowing physics builders are included:
    - G4DecayPhysics
    - G4RadioactiveDecayPhysics
    - G4EmStandardPhysics
*////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef SIMPLEPHYSICSLIST_h
#define SIMPLEPHYSICSLIST_h

#include "G4VModularPhysicsList.hh"

class SimplePhysicsList: public G4VModularPhysicsList
{
    public:
        SimplePhysicsList();
        virtual ~SimplePhysicsList();

        virtual void SetCuts();
};

#endif

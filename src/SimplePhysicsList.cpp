/*///////////////////////////////////////////////////////////////////////////////////////////////
This is the SimplePhysicsList.cc file from the slac_rf project send to me by Steven.
*////////////////////////////////////////////////////////////////////////////////////////////////

#include "SimplePhysicsList.hh"

#include "G4DecayPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4EmStandardPhysics.hh"

SimplePhysicsList::SimplePhysicsList() : G4VModularPhysicsList(){
  SetVerboseLevel(1);

  // Default physics
  RegisterPhysics(new G4DecayPhysics());

  // Radioactive decay
  RegisterPhysics(new G4RadioactiveDecayPhysics());

  // EM physics
  RegisterPhysics(new G4EmStandardPhysics());
}

SimplePhysicsList::~SimplePhysicsList()
{ 
}

void SimplePhysicsList::SetCuts()
{
  G4VUserPhysicsList::SetCuts();
} 

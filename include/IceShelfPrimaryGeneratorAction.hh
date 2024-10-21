/*//////////////////////////////////////////////////////////////////////////////////////////////////
This is the PrimaryGeneratorAction class. It deals with generating the primary particles. This is
implemented in the function GeneratePrimaries().
In this project the primaries are generated using a G4ParticleGun. A pointer to the ParticleGun instance is stored in the first field. The second field is used to store a pointer to a vector of all primary particles that have not been shot in the ice yet by the particle gun. It is initialized by reading an inputfile containing all primaries the gun can choose from. Its initialization is done in the main function of this project. The third field contains a pointer to the EventAction object, which is needed to be able to tell the EventAction object the weight of the primary that is being propagated through the ice.
*///////////////////////////////////////////////////////////////////////////////////////////////////


#ifndef ICESHELFPRIMARYGENERATORACTION_h
#define ICESHELFPRIMARYGENERATORACTION_h

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "globals.hh"
#include "IceShelfEventAction.hh"

#include <vector>
#include <map>

class G4ParticleGun;
class G4Event;
class G4Box;

class IceShelfPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
    public:
        IceShelfPrimaryGeneratorAction(std::vector<std::vector<G4double>>*, IceShelfEventAction*);
        virtual ~IceShelfPrimaryGeneratorAction();

        // method from the base class
        virtual void GeneratePrimaries(G4Event*);

        // method to access particle gun
        const G4ParticleGun* GetParticleGun() const { return fParticleGun; }

    private:
        G4ParticleGun* fParticleGun; // pointer to a G4 gun class
        std::vector<std::vector<G4double>>* fPrimaryList; // list containing all the primaries not used yet
        std::map<G4int, G4String> fParticleList; // map linking Corsika ID's to particle names
        IceShelfEventAction* fEventAction; // pointer to the EventAction object
};

#endif

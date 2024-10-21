#include "IceShelfPrimaryGeneratorAction.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include <stdexcept>

IceShelfPrimaryGeneratorAction::IceShelfPrimaryGeneratorAction
                                                   (std::vector<std::vector<G4double>>* primaryList,
                                                    IceShelfEventAction* eventAction)
: G4VUserPrimaryGeneratorAction(),
  fParticleGun(0),
  fPrimaryList(primaryList),
  fEventAction(eventAction)
{
    // Initializing the particle gun
    G4int n_particle = 1; // If this number is increased to n, the particle gun will generate n primaries
                          // instead of only one. Note however that the initial conditions of these n 
                          // primaries (position, energy, ...) will be identical! The particle gun does
                          // not introduce any randomness. If you want to achieve randomness you have 
                          // to make sure the input you give to the particle gun (position, energy, ...)
                          // is in itself a random variable!
    fParticleGun = new G4ParticleGun(n_particle);

    // Constructing the particle list, which maps Corsika ID's to particle names
    fParticleList[1] = "gamma";
    fParticleList[2] = "e+";
    fParticleList[3] = "e-";
    fParticleList[5] = "mu+";
    fParticleList[6] = "mu-";
    fParticleList[7] = "pi0";
    fParticleList[8] = "pi+";
    fParticleList[9] = "pi-";
    fParticleList[10] = "kaon0L";
    fParticleList[11] = "kaon+";
    fParticleList[12] = "kaon-";
    fParticleList[13] = "neutron";
    fParticleList[14] = "proton";
    fParticleList[15] = "anti_proton";
    fParticleList[16] = "kaon0S";
    fParticleList[17] = "eta";
    fParticleList[18] = "lambda";
    fParticleList[19] = "sigma+";
    fParticleList[20] = "sigma0";
    fParticleList[21] = "sigma-";
    fParticleList[22] = "xi0";
    fParticleList[23] = "xi-";
    fParticleList[24] = "omega-";
    fParticleList[25] = "anti_neutron";
    fParticleList[26] = "anti_lambda";
    fParticleList[27] = "anti_sigma+";
    fParticleList[28] = "anti_sigma0";
    fParticleList[29] = "anti_sigma-";
    fParticleList[30] = "anti_xi0";
    fParticleList[31] = "anti_xi-";
    fParticleList[32] = "anti_omega-";
    fParticleList[201] = "deuteron";
    fParticleList[301] = "triton";
    fParticleList[302] = "He3";
    fParticleList[402] = "alpha";
}

IceShelfPrimaryGeneratorAction::~IceShelfPrimaryGeneratorAction()
{
    delete fParticleGun;
}

void IceShelfPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
    // This function is called at the beginning of each event

    // Choosing a random particle from the primary list
    
    //G4cout << "Choosing a random particle. There are " << fPrimaryList->size() << " left." << G4endl;
    G4double numberOfEntries = fPrimaryList->size();
    G4int randomEntryNumber = (G4int)(floor(G4UniformRand()*numberOfEntries));
    //G4cout << "Found a random number: " << randomEntryNumber << G4endl;
    std::vector<G4double> randomEntry = fPrimaryList->at(randomEntryNumber);
    //G4cout << "Found a random particle." << G4endl;

    // Getting its values
    G4int part_id = (G4int)randomEntry.at(0);
    G4double px = randomEntry.at(1);
    G4double py = randomEntry.at(2);
    G4double pz = randomEntry.at(3);
    G4double x = randomEntry.at(4);
    G4double y = randomEntry.at(5);
    G4double z = randomEntry.at(6);
    G4double t = randomEntry.at(7);
    G4double weight = randomEntry.at(8);
    /*
    G4cout << "Read all particle properties." << G4endl;
    G4cout << "part_id: " << part_id << G4endl;
    G4cout << "px: " << px/GeV << " GeV/c" << G4endl;
    G4cout << "py: " << py/GeV << " GeV/c" << G4endl;
    G4cout << "pz: " << pz/GeV << " GeV/c" << G4endl;
    G4cout << "x: " << x/cm << " cm"  << G4endl;
    G4cout << "y: " << y/cm << " cm" << G4endl;
    G4cout << "z: " << z/cm << " cm" << G4endl;
    G4cout << "t: " << t/ns << " ns" << G4endl;
    G4cout << "weight: " << weight << G4endl;
    */

    // Configuring the particle gun

    if(fParticleList[part_id] == ""){
        G4cerr << "Encountered particle id in primary generator action which is not implemented. Particle id is " << part_id << G4endl;
        exit(EXIT_FAILURE);
    }

    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4String particleName;
    G4ParticleDefinition* particle;

    particle = particleTable->FindParticle(particleName=fParticleList[part_id]);

    G4cout << "Atomic number Z of primary is " << particle->GetAtomicNumber() << G4endl;
    G4cout << "Atomic mass number A of primary is " << particle->GetAtomicMass() << G4endl;
    G4cout << "Mass of the particle is " << particle->GetPDGMass()/MeV << " MeV" << G4endl;

    fParticleGun->SetParticleDefinition(particle);
    G4double p = (G4double)(sqrt(px*px + py*py + pz*pz));
    fParticleGun->SetParticleMomentum(p);
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(px/p, py/p, pz/p));
    fParticleGun->SetParticlePosition(G4ThreeVector(x, y, z));
    //G4cout << "Configured the particle gun." << G4endl;
    
    // Setting the particle id of the primary particle
    fEventAction->SetPrimaryPartName(particleName);
    //G4cout << "Set the particle id of the primary particle." << G4endl;

    // Setting the starting time at which the primary enters the ice in the event object
    fEventAction->SetStartingTime(t);
    //G4cout << "Set the starting time at which the primary enters the ice." << G4endl;

    // Setting the weight of the primary particle in the event object
    fEventAction->SetPrimaryWeight(weight);
    //G4cout << "Set the weight of the primary particle." << G4endl;

    // Setting the charge of the primary particle in the event object
    fEventAction->SetPrimaryCharge(particle->GetPDGCharge());

    // Setting the starting position of the primary particle in the event object
    fEventAction->SetPrimaryStartPos(G4ThreeVector(x, y, z));

    // Setting the starting momentum (p*c) of the primary particle in the event object
    fEventAction->SetPrimaryStartMomentum(G4ThreeVector(px, py, pz));

    // Setting the mass (m*c^2) of the primary particle in the event object
    fEventAction->SetPrimaryMass(particle->GetPDGMass());

    // Removing the random particle from the primary list, so that it will not be choosen again
    fPrimaryList->erase(fPrimaryList->begin() + randomEntryNumber);
    //G4cout << "Removed the primary from the list of primaries still to simulate." << G4endl;

    // Generating a primary
    G4cout << "Activating the particle gun." << G4endl;
    fParticleGun->GeneratePrimaryVertex(anEvent);

}

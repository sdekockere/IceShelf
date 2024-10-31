#include "IceShelfSteppingAction.hh"
#include "IceShelfEventAction.hh"
#include "IceShelfRunAction.hh"
#include "IceShelfDetectorConstruction.hh"
#include "IceShelfAnalysis.hh"
#include "IceShelfAntennas.hh"

#include "G4Step.hh"
#include "IceShelfGlobalVariables.hh"

#include <math.h>

#include "G4VProcess.hh"
#include "G4TrackStatus.hh"

#include "IceDensityModels.hh"

#include "IceRayTracing.hh"
#include "AirToIceRayTracing.h"

IceShelfSteppingAction::IceShelfSteppingAction(IceShelfDetectorConstruction* detConstruction, IceShelfRunAction* runAction, IceShelfEventAction* eventAction, G4double* snapshotTimes, G4double zenithAngleRad)
 : G4UserSteppingAction(),
   fDetConstruction(detConstruction),
   fRunAction(runAction),
   fEventAction(eventAction),
   fCurrentTrackID(0),
   fCurrentParticleDepthIndex(0),
   fCurrentParticleChargeIndex(0),
   fCurrentParticleSnapshotTimeIndex(0),
   fSnapshotTimes(snapshotTimes),
   fZenithAngleRad(zenithAngleRad)
{}

IceShelfSteppingAction::~IceShelfSteppingAction()
{}

// Based on Anne Zilles' work
void IceShelfSteppingAction::CalcRadioEmissionEndpoint(const G4Step* step)
{

    IceShelfAntennas* iceShelfAntennas = fRunAction->GetIceShelfAntennas();

    G4Track* track = step->GetTrack();

    // Unit of charge in geant4 is the charge of the positron (charge e+ = 1.)
    G4int pdgCharge = track->GetParticleDefinition()->GetPDGCharge();
    if(pdgCharge == 0){
        return; // No charge, so no radiation...
    }

    G4double kinEnergy = track->GetKineticEnergy(); // The kinetic energy at end of step
    if(kinEnergy < kinEnergyThreshold){
        return; // Too low in energy to be interesting
    }

    // PreStepPoint, which marks the beginning of the step
    // Remember to add starting time to global time!
    G4ThreeVector preStepPos = step->GetPreStepPoint()->GetPosition();
    G4double preTime = step->GetPreStepPoint()->GetGlobalTime() + fEventAction->GetStartingTime();

    // PostStepPoint, which marks the end of the step
    // Remember to add starting time to global time!
    G4ThreeVector postStepPos = step->GetPostStepPoint()->GetPosition();
    G4double postTime = step->GetPostStepPoint()->GetGlobalTime() + fEventAction->GetStartingTime();

    // The direction and length of the step
    // NOTE: use the step->GetStepLength() function instead? Gives length corrected for multiple scat
    // NOTE: or only use the corrected one once to check if ALSO not 0 (is what Anne does...)
    G4ThreeVector stepDir = postStepPos - preStepPos;
    G4double stepLength = stepDir.mag();

    if(stepLength == 0.){
        //G4cout << "Found a stepLength = 0..." << G4endl;
        return; // Make sure we're not going to divide by 0 later on
    }
    stepDir = stepDir/stepLength;
    
    // Beta during this step
    G4double betaVal = (1./c0)*(stepLength/(postTime - preTime));
    G4ThreeVector beta = betaVal*stepDir;

    // The weight of the primary particle
    G4double primaryWeight = fEventAction->GetPrimaryWeight();

    // From here on we need the positions of the antennas. This means we'll have to do this part
    // for every single antenna in the antennas list.
    for(G4int k = 0; k < iceShelfAntennas->GetNumberOfAntennas(); k++){

        // The observer position
        G4ThreeVector observerPos = G4ThreeVector(iceShelfAntennas->GetAntennaCoord(k, 0),
                                                  iceShelfAntennas->GetAntennaCoord(k, 1),
                                                  iceShelfAntennas->GetAntennaCoord(k, 2));

        // Getting the contributions to the electric field a the antenna from start and end point
        G4ThreeVector startE;
        G4ThreeVector endE;
        G4int iStart;
        G4int iEnd;
        G4int gotContribution;

        if(includeDirectRays){ // The direct ray

            gotContribution = GetEfieldEP(k, observerPos, preStepPos, postStepPos, preTime, postTime,
                                                    beta, pdgCharge, 0, startE, endE, iStart, iEnd);

            if(gotContribution){

                // Taking the weight of the primary into account
                startE = primaryWeight*startE;
                endE = primaryWeight*endE;

                // Write the electrical field contributions to the array

                // Checking if the time window for registration is large enough
                if(iStart < 0 || iStart >= static_cast<int>(ceil(traceLength/sPeriod)) || iEnd < 0 
                        || iEnd >= static_cast<int>(ceil(traceLength/sPeriod))){

                    G4cerr << "Index out of bounds: " << G4endl;
                    G4cerr << "iStart = " << iStart << G4endl;
                    G4cerr << "iEnd = " << iEnd << G4endl;
                    G4cerr << "Size of antenna trace: "
                                            << static_cast<int>(ceil(traceLength/sPeriod)) << G4endl;

                    G4cerr << "startE with weight included: " << CGSToSIFactor*startE
                                                                                  << " V/m" << G4endl;
                    G4cerr << "endE with weight included: " << CGSToSIFactor*endE << " V/m" << G4endl;
                    
                    exit(EXIT_FAILURE);

                }

                fRunAction->AddToAntennaTrace(k, startE, iStart, 0);
                fRunAction->AddToAntennaTrace(k, endE, iEnd, 0);

                if(writeSARadSeparateToo){
                    // Checking if it is the primary particle doing its first step
                    if(track->GetTrackID() == 1 && track->GetCurrentStepNumber() == 1){
                        fRunAction->AddToAntennaTrace(k, startE, iStart, 2);
                    }
                }

            }

        }

        if(includeIndirectRays){ // The indirect ray

            gotContribution = GetEfieldEP(k, observerPos, preStepPos, postStepPos, preTime, postTime,
                                                    beta, pdgCharge, 1, startE, endE, iStart, iEnd);

            if(gotContribution){

                // Taking the weight of the primary into account
                startE = primaryWeight*startE;
                endE = primaryWeight*endE;

                // Write the electrical field contributions to the array

                // Checking if the time window for registration is large enough
                if(iStart < 0 || iStart >= static_cast<int>(ceil(traceLength/sPeriod)) || iEnd < 0 
                        || iEnd >= static_cast<int>(ceil(traceLength/sPeriod))){

                    G4cerr << "Index out of bounds: " << G4endl;
                    G4cerr << "iStart = " << iStart << G4endl;
                    G4cerr << "iEnd = " << iEnd << G4endl;
                    G4cerr << "Size of antenna trace: "
                                            << static_cast<int>(ceil(traceLength/sPeriod)) << G4endl;

                    G4cerr << "startE with weight included: " << CGSToSIFactor*startE
                                                                                  << " V/m" << G4endl;
                    G4cerr << "endE with weight included: " << CGSToSIFactor*endE << " V/m" << G4endl;
                    
                    exit(EXIT_FAILURE);

                }

                fRunAction->AddToAntennaTrace(k, startE, iStart, 1);
                fRunAction->AddToAntennaTrace(k, endE, iEnd, 1);

                if(writeSARadSeparateToo){
                    // Checking if it is the primary particle doing its first step
                    if(track->GetTrackID() == 1 && track->GetCurrentStepNumber() == 1){
                        fRunAction->AddToAntennaTrace(k, startE, iStart, 2);
                    }
                }

            }

        }

    }
}

void IceShelfSteppingAction::UserSteppingAction(const G4Step* step)
{
    /*
    // For keeping track of a particle. To track primary particle processes use if(_parentID == 1).
    G4Track* _track = step->GetTrack();
    G4int _trackID = _track->GetTrackID(); 
    G4int _parentID = _track->GetParentID();
    G4String _name = _track->GetParticleDefinition()->GetParticleName();
    G4double _kinEnergy = _track->GetKineticEnergy();
    //if(_parentID == 1){
    //if(_trackID == 1){
        G4String _process = "Not found...?";
        if(_trackID == 1){
            _process = "primary";
        } else {
            _process = _track->GetCreatorProcess()->GetProcessName();
        }
        auto _preStep = step->GetPreStepPoint()->GetPosition();
        auto _postStep = step->GetPostStepPoint()->GetPosition();
        G4double _preStepTime = step->GetPreStepPoint()->GetGlobalTime() 
                                    + fEventAction->GetStartingTime(); 
        G4double _postStepTime = step->GetPostStepPoint()->GetGlobalTime()
                                    + fEventAction->GetStartingTime(); 
        G4cout << "Track ID, parent ID, particle name, process, post kin energy (MeV), preStep (cm), postStep (cm), preStepTime (ns), postStepTime (ns): "; 
        G4cout << _trackID << "," << _parentID << "," << _name << "," << _process << "," << _kinEnergy/MeV << "," << _preStep/cm << "," << _postStep/cm << "," << _preStepTime/ns << "," << _postStepTime/ns << G4endl;
        G4cout << "Is first step in volume: " << step->IsFirstStepInVolume() << ", and step number"
            << " is " << _track->GetCurrentStepNumber() << G4endl;
    //}
    */

    // Get the energy deposit of this step
    // This is the sum of the energy deposited by the energy loss process and
    // the energy lost by secondaries which have not been generated because each of their energies 
    // was below the cut threshold

    G4double energyDep = step->GetTotalEnergyDeposit();

    auto volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();

    //////////////////////////////////////////////////////////////////////////////////
    ////// Stuff that needs to be done for particles outside the ice block only //////
    //////////////////////////////////////////////////////////////////////////////////

    if(volume->GetName() == "World"){

        fEventAction->AddLeakedEnergy(energyDep);

    } 

    /////////////////////////////////////////////////////////////////////////////////
    ////// Stuff that needs to be done for particles inside the ice block only //////
    /////////////////////////////////////////////////////////////////////////////////

    /*
    else {


    }
    */

    ////////////////////////////////////////////////////////////////////////////////////////////
    ////// Stuff that needs to be done for particles both inside or outside the ice block //////
    ////////////////////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////
    ////// FILLING THE 2D ENERGY DEPOSIT HISTOGRAMS //////
    //////////////////////////////////////////////////////

    // Find the point in the middle between PreStepPoint and PostStepPoint

    G4ThreeVector preStepPointPos = step->GetPreStepPoint()->GetPosition();
    G4double preStepPointPosX = preStepPointPos.getX();
    G4double preStepPointPosY = preStepPointPos.getY();
    G4double preStepPointPosZ = preStepPointPos.getZ();

    G4ThreeVector postStepPointPos = step->GetPostStepPoint()->GetPosition();
    G4double postStepPointPosX = postStepPointPos.getX();
    G4double postStepPointPosY = postStepPointPos.getY();
    G4double postStepPointPosZ = postStepPointPos.getZ();
    
    G4double middlePointX = (preStepPointPosX + postStepPointPosX)/2.;
    G4double middlePointY = (preStepPointPosY + postStepPointPosY)/2.;
    G4double middlePointZ = (preStepPointPosZ + postStepPointPosZ)/2.;

    // Making the coordinate transformation needed to go from the geant4 axis system to 
    // the shower axis system (y-axis aligned with the shower) on the surface of the ice block.

    G4double preAxisCoords[3];
    GetShowerAxisCoordinates(preStepPointPosX, preStepPointPosY, preStepPointPosZ, preAxisCoords);

    G4double postAxisCoords[3];
    GetShowerAxisCoordinates(postStepPointPosX, postStepPointPosY, postStepPointPosZ,
                                                                                    postAxisCoords);
    G4double middlePointAxisCoords[3];
    GetShowerAxisCoordinates(middlePointX, middlePointY, middlePointZ, middlePointAxisCoords);

    // Store the energy density in the 2D histogram if z coordinate falls within the slice thickness.
    // The weight of each bin represents the total energy deposited in this bin.
    // ATTENTION: each primary particle read from the corsika output file 
    // has its own weight, which we need to take into account here as well! 
    // The EventAction object knows the weight of the current primary particle
    // (i.e. the current corsika output particle we're propagating through the ice)

    G4String primaryPartName = fEventAction->GetPrimaryPartName();
    G4double primaryWeight = fEventAction->GetPrimaryWeight();

    // The energy density slice histograms
    if(middlePointAxisCoords[2] < blockSize/2. && middlePointAxisCoords[2] > -1.*blockSize/2.){
        auto analysisManager = G4AnalysisManager::Instance();
        if(primaryPartName == "gamma"){
            analysisManager->FillH2(0, middlePointAxisCoords[0]/cm, middlePointAxisCoords[1]/cm, primaryWeight*(energyDep/MeV)/((blockSize/cm)*(blockSize/cm)*(blockSize/cm)));
        } else if(primaryPartName == "e+" || primaryPartName == "e-"){
            analysisManager->FillH2(1, middlePointAxisCoords[0]/cm, middlePointAxisCoords[1]/cm, primaryWeight*(energyDep/MeV)/((blockSize/cm)*(blockSize/cm)*(blockSize/cm)));
        } else if(primaryPartName == "mu+" || primaryPartName == "mu-"){
            analysisManager->FillH2(2, middlePointAxisCoords[0]/cm, middlePointAxisCoords[1]/cm, primaryWeight*(energyDep/MeV)/((blockSize/cm)*(blockSize/cm)*(blockSize/cm)));
        } else {
            analysisManager->FillH2(3, middlePointAxisCoords[0]/cm, middlePointAxisCoords[1]/cm, primaryWeight*(energyDep/MeV)/((blockSize/cm)*(blockSize/cm)*(blockSize/cm)));
        }
    }

    // The radial energy profile histograms
    G4double radiusPart = sqrt(middlePointAxisCoords[0]*middlePointAxisCoords[0] + middlePointAxisCoords[2]*middlePointAxisCoords[2]);
    if(primaryPartName == "gamma"){
        fRunAction->FillH2(0, radiusPart/cm, middlePointAxisCoords[1]/cm, primaryWeight*(energyDep/MeV));
    } else if(primaryPartName == "e+" || primaryPartName == "e-"){
        fRunAction->FillH2(1, radiusPart/cm, middlePointAxisCoords[1]/cm, primaryWeight*(energyDep/MeV));
    } else if(primaryPartName == "mu+" || primaryPartName == "mu-"){
        fRunAction->FillH2(2, radiusPart/cm, middlePointAxisCoords[1]/cm, primaryWeight*(energyDep/MeV));
    } else {
        fRunAction->FillH2(3, radiusPart/cm, middlePointAxisCoords[1]/cm, primaryWeight*(energyDep/MeV));
    }

    /////////////////////////////////////////////////
    ////// FILLING DEPTH AND CHARGE HISTOGRAMS //////
    ////// FILLING THE SNAPSHOT HISTOS         //////
    /////////////////////////////////////////////////

    G4int trackID = step->GetTrack()->GetTrackID();
    G4double depth = -1.*postAxisCoords[1];
    G4double time = step->GetPostStepPoint()->GetGlobalTime() + fEventAction->GetStartingTime();

    G4double depthBinWidth = maxDepthValue/depthNumberOfBins;       
    G4double timeBinWidth = maxTimeValue/timeNumberOfBins;

    if(trackID != fCurrentTrackID){ // Means Geant4 started calculations for a new particle

        fCurrentTrackID = trackID;

        G4double startDepth = -1.*preAxisCoords[1];

        // The index of the first bin of the depth histogram 
        // we should fill if the current particle ever reaches the corresponding depth
        fCurrentParticleDepthIndex = static_cast<int>(ceil(startDepth/depthBinWidth));

        G4double startTime = step->GetPreStepPoint()->GetGlobalTime()
                                    + fEventAction->GetStartingTime();

        // The index of the first bin of the charge histogram
        // we should fill if the current particle ever reaches the corresponding time
        fCurrentParticleChargeIndex = static_cast<int>(ceil(startTime/timeBinWidth));

        G4int snapshotStartIndex = 0;
        while(snapshotStartIndex < numberOfSnapshotDepthsMassArea && 
                fSnapshotTimes[snapshotStartIndex] < startTime){
            snapshotStartIndex++;
        }
        fCurrentParticleSnapshotTimeIndex = snapshotStartIndex;

        /*
        G4cout << "Started new particle. Starting depth is " << startDepth/m << " m" << G4endl;
        G4cout << "Starting kin energy is " << step->GetPreStepPoint()->GetKineticEnergy()/GeV << 
               " GeV" << G4endl;
        G4cout << "Starting time is " << startTime/ns << " ns" << G4endl;
        G4cout << "Corresponding depth index is " << fCurrentParticleDepthIndex << G4endl;
        G4cout << "Corresponding charge index is " << fCurrentParticleChargeIndex << G4endl;
        G4cout << "Corresponding snapshot index is " << fCurrentParticleSnapshotTimeIndex << G4endl;
        */

    }

    G4double preStepKinEnergy = step->GetPreStepPoint()->GetKineticEnergy();
    G4String particleName = step->GetTrack()->GetParticleDefinition()->GetParticleName();

    // Checking whether we need to fill the depth histogram

    G4double depthKinELim = 0.;
    if(particleName == "gamma" || particleName == "pi0"){
        depthKinELim = depthKinELimVals[3];
    } else if(particleName == "e+" || particleName == "e-"){
        depthKinELim = depthKinELimVals[2];
    } else if(particleName == "mu+" || particleName == "mu-"){
        depthKinELim = depthKinELimVals[1];
    } else {
        depthKinELim = depthKinELimVals[0];
    }

    if(preStepKinEnergy > depthKinELim && depth >= fCurrentParticleDepthIndex*depthBinWidth){

        /*
        G4cout << "Depth is now " << depth/m << " m, going to fill the histogram at " <<
        (fCurrentParticleDepthIndex + 0.5)*depthBinWidth/m << " m" << G4endl;
        G4cout << "Energy during this step was " << step->GetPreStepPoint()->GetKineticEnergy()/GeV
         << " GeV" << G4endl;
        */

        G4int depthHistoNumber;

        if(particleName == "gamma"){
            //G4cout << "Filling gamma histo..." << G4endl;
            depthHistoNumber = 4;
        } else if(particleName == "e+"){
            //G4cout << "Filling e+ histo..." << G4endl;
            depthHistoNumber = 5;
        } else if(particleName == "e-"){
            //G4cout << "Filling e- histo..." << G4endl;
            depthHistoNumber = 6;
        } else if(particleName == "mu+"){
            //G4cout << "Filling mu+ histo..." << G4endl;
            depthHistoNumber = 7;
        } else if(particleName == "mu-"){
            //G4cout << "Filling mu- histo..." << G4endl;
            depthHistoNumber = 8;
        } else {
            //G4cout << "Filling hadrons histo..." << G4endl;
            depthHistoNumber = 9;
        }

        // Filling the depth histogram
        auto analysisManager = G4AnalysisManager::Instance();
        analysisManager->FillH2(depthHistoNumber, 
                                (fCurrentParticleDepthIndex + 0.5)*depthBinWidth/m,
                                preStepKinEnergy,
                                primaryWeight);
        
        // Incrementing fCurrentParticleDepthIndex
        fCurrentParticleDepthIndex += 1;

        // Making sure we fill all the bins we can
        while(depth >= fCurrentParticleDepthIndex*depthBinWidth){

            analysisManager->FillH2(depthHistoNumber, 
                                    (fCurrentParticleDepthIndex + 0.5)*depthBinWidth/m,
                                    preStepKinEnergy,
                                    primaryWeight);

            /*
            G4cout << "Also filled histogram at " << 
            (fCurrentParticleDepthIndex + 0.5)*depthBinWidth/m << " m" << G4endl;
            */

            fCurrentParticleDepthIndex += 1;
            
        }

    }

    // Checking whether we need to fill the charge histogram
    if(preStepKinEnergy > chargeKinELim && time >= fCurrentParticleChargeIndex*timeBinWidth){

        /*
        G4cout << "Time is now " << time/ns << " ns, going to fill the histogram at " <<
        (fCurrentParticleChargeIndex + 0.5)*timeBinWidth/ns << " ns" << G4endl;
        G4cout << "Energy during this step was " << step->GetPreStepPoint()->GetKineticEnergy()/GeV
         << " GeV" << G4endl;
        */

        G4double charge = step->GetTrack()->GetParticleDefinition()->GetPDGCharge();

        G4int timeHistoNumber;

        if(primaryPartName == "gamma"){
            timeHistoNumber = 10;
        } else if(primaryPartName == "e+" || primaryPartName == "e-"){
            timeHistoNumber = 11;
        } else if(primaryPartName == "mu+" || primaryPartName == "mu-"){
            timeHistoNumber = 12;
        } else {
            timeHistoNumber = 13;
        }

        // Filling the charge histogram
        auto analysisManager = G4AnalysisManager::Instance();
        analysisManager->FillH2(timeHistoNumber, 
                                (fCurrentParticleChargeIndex + 0.5)*timeBinWidth/ns,
                                preStepKinEnergy,
                                primaryWeight*charge/eplus);

        // Incrementing fCurrentParticleChargeIndex
        fCurrentParticleChargeIndex += 1;

        // Making sure we fill all the bins we can
        while(time >= fCurrentParticleChargeIndex*timeBinWidth){

            analysisManager->FillH2(timeHistoNumber, 
                                    (fCurrentParticleChargeIndex + 0.5)*timeBinWidth/ns,
                                    preStepKinEnergy,
                                    primaryWeight*charge/eplus);
            
            /*
            G4cout << "Also filled histogram at " <<
            (fCurrentParticleChargeIndex + 0.5)*timeBinWidth/ns << " ns" << G4endl;
            */

            fCurrentParticleChargeIndex += 1;

        }

    }

    // Checking whether we need to fill one of the snapshot ntuples

    while(preStepKinEnergy > snapshotKinELim
            && fCurrentParticleSnapshotTimeIndex < numberOfSnapshotDepthsMassArea
            && time >= fSnapshotTimes[fCurrentParticleSnapshotTimeIndex]){

        // Get a better estimate of where the particle was exactly at the snapshot time

        // Estimation of the velocity during the last step
        G4double timeOfPreStep = step->GetPreStepPoint()->GetGlobalTime()
                                    + fEventAction->GetStartingTime();
        G4double velLastStep[3];
        velLastStep[0] = (postAxisCoords[0] - preAxisCoords[0])
                         /(time - timeOfPreStep);
        velLastStep[1] = (postAxisCoords[1] - preAxisCoords[1])
                         /(time - timeOfPreStep);
        velLastStep[2] = (postAxisCoords[2] - preAxisCoords[2])
                         /(time - timeOfPreStep);

        // Estimation of the point where the particle was at time of the snapshot
        G4double snapshotPosAxisCoords[3]; 
        snapshotPosAxisCoords[0] = postAxisCoords[0]
                     - velLastStep[0]*(time - fSnapshotTimes[fCurrentParticleSnapshotTimeIndex]);
        snapshotPosAxisCoords[1] = postAxisCoords[1]
                     - velLastStep[1]*(time - fSnapshotTimes[fCurrentParticleSnapshotTimeIndex]);
        snapshotPosAxisCoords[2] = postAxisCoords[2]
                     - velLastStep[2]*(time - fSnapshotTimes[fCurrentParticleSnapshotTimeIndex]);

        G4double PDGCharge = step->GetTrack()->GetParticleDefinition()->GetPDGCharge();

        auto analysisManager = G4AnalysisManager::Instance();

        // Filling the snapshot Ntuple
        if(PDGCharge != 0.){
            analysisManager->FillNtupleDColumn(fCurrentParticleSnapshotTimeIndex, 0, 
                                                    snapshotPosAxisCoords[0]/m);
            analysisManager->FillNtupleDColumn(fCurrentParticleSnapshotTimeIndex, 1, 
                                                    snapshotPosAxisCoords[1]/m);
            analysisManager->FillNtupleDColumn(fCurrentParticleSnapshotTimeIndex, 2, 
                                                    snapshotPosAxisCoords[2]/m);
            analysisManager->FillNtupleDColumn(fCurrentParticleSnapshotTimeIndex, 3, 
                                                    PDGCharge/eplus);
            analysisManager->FillNtupleDColumn(fCurrentParticleSnapshotTimeIndex, 4, 
                                                    primaryWeight);
            analysisManager->AddNtupleRow(fCurrentParticleSnapshotTimeIndex);
        }

        fCurrentParticleSnapshotTimeIndex++;

    }

    /////////////////////////////////////////////////////////////////////
    ////// Calculating the radio emission using endpoint formalism //////
    /////////////////////////////////////////////////////////////////////
    
    CalcRadioEmissionEndpoint(step);

}

/*
G4ThreeVector IceShelfSteppingAction::RotateVector(G4ThreeVector vect, G4ThreeVector axis,
                                                                                G4double theta){

    // The vect components
    G4double x = vect.getX();
    G4double y = vect.getY();
    G4double z = vect.getZ();

    // The axis components
    G4double ux = axis.getX();
    G4double uy = axis.getY();
    G4double uz = axis.getZ();

    // The cos and sin values
    G4double cosVal = cos(theta);
    G4double sinVal = sin(theta);

    // The rotated vector components
    G4double xr = ( cosVal + ux*ux*(1.-cosVal) ) * x
                    + ( ux*uy*(1.-cosVal) - uz*sinVal ) * y
                    + ( ux*uz*(1.-cosVal) + uy*sinVal ) * z;
    G4double yr = ( uy*ux*(1.-cosVal) + uz*sinVal ) * x
                    + ( cosVal + uy*uy*(1.-cosVal) ) * y
                    + ( uy*uz*(1.-cosVal) - ux*sinVal ) * z;
    G4double zr = ( uz*ux*(1.-cosVal) - uy*sinVal ) * x
                    + ( uz*uy*(1.-cosVal) + ux*sinVal ) * y
                    + ( cosVal + uz*uz*(1.-cosVal) ) * z;

    // The rotated vector
    G4ThreeVector rotVect = G4ThreeVector(xr, yr, zr);
    return rotVect;

}
*/

// Faster version for rotations (does not evaluate cos and sin function)
// Works perfect if you e.g. calculate the cos(rot_angle) through scalar product
// Rotates according to right-hand rule
// WARNING: assumes rot_angle <= pi, so that sin(rot_angle) >= 0
G4ThreeVector IceShelfSteppingAction::RotateVectorCosPi(G4ThreeVector vect, G4ThreeVector axis,
                                                                                G4double cosVal){

    // The vect components
    G4double x = vect.getX();
    G4double y = vect.getY();
    G4double z = vect.getZ();

    // The axis components
    G4double ux = axis.getX();
    G4double uy = axis.getY();
    G4double uz = axis.getZ();

    // The cos and sin values
    G4double sinVal = sqrt(1-cosVal*cosVal);

    // The rotated vector components
    G4double xr = ( cosVal + ux*ux*(1.-cosVal) ) * x
                    + ( ux*uy*(1.-cosVal) - uz*sinVal ) * y
                    + ( ux*uz*(1.-cosVal) + uy*sinVal ) * z;
    G4double yr = ( uy*ux*(1.-cosVal) + uz*sinVal ) * x
                    + ( cosVal + uy*uy*(1.-cosVal) ) * y
                    + ( uy*uz*(1.-cosVal) - ux*sinVal ) * z;
    G4double zr = ( uz*ux*(1.-cosVal) - uy*sinVal ) * x
                    + ( uz*uy*(1.-cosVal) + ux*sinVal ) * y
                    + ( cosVal + uz*uz*(1.-cosVal) ) * z;

    // The rotated vector
    G4ThreeVector rotVect = G4ThreeVector(xr, yr, zr);
    return rotVect;

}

// Returns the values needed to evaluate the end-point formalism, NOT using ray tracing
// Input:
//  - emitterObsVect, the G4ThreeVector giving the distance from emitter to observer
//  - emitterTime, G4double giving the time of emission
//  - beta, the G4ThreeVector giving the beta vector (v/c) of the emitter
//  #########
//  - obsTime, G4double where the observation time will be stored
//  - launchDir, G4ThreeVector where the launch direction will be stored
//  - receiveDir, G4ThreeVector where the receive direction will be stored
//  - doppler, G4double where the Doppler term will be stored (1 - n*beta . rhat)
//  - R, G4double where the value for R will be stored
// Output:
// - 1 if function ran normally
// NOTE:
// direction of launchDir and receiveDir are defined as if the emitter is shooting an arrow to
// the receiver, i.e. the launch dir points from emitter to receiver, while the receive dir 
// points from receiver AWAY from emitter
G4int IceShelfSteppingAction::GetEPVals(G4ThreeVector emitterObsVect,
                    G4double emitterTime, G4ThreeVector beta,
                    G4double &obsTime, G4ThreeVector &launchDir, G4ThreeVector &receiveDir,
                    G4double &doppler, G4double &R){

    // Setting the value for R
    R = emitterObsVect.mag();

    // Calculating the time of arrival of the radiation
    obsTime = emitterTime + nIce*R/c0;

    // Setting the launch direction and receive direction
    launchDir = emitterObsVect/R; 
    receiveDir = launchDir;    

    // Calculating the Doppler term
    doppler = 1. - nIce * (beta.dot(launchDir));
    if(doppler == 0){
        G4cerr << "WARNING: encountered a doppler term = 0" << G4endl;
    }

    return 1;

}

// Returns the values needed to evaluate the end-point formalism using ray tracing
// Input:
//  - k, unsigned int giving the antenna number
//  - emitterPos, the G4ThreeVector for the position of the emitter
//  - emitterObsVect, the G4ThreeVector giving the distance from emitter to observer
//  - emitterTime, G4double giving the time of emission
//  - beta, the G4ThreeVector giving the beta vector (v/c) of the emitter
//  - dirOrIndirRay, should be 0 if you want to calculate the values for the direct ray,
//          should be 1 if you want to calculate the values for the indirect ray
//  #########
//  - obsTime, G4double where the observation time will be stored
//  - launchDir, G4ThreeVector where the launch direction will be stored
//  - receiveDir, G4ThreeVector where the receive direction will be stored
//  - doppler, G4double where the Doppler term will be stored (1 - n*beta . rhat)
//  - R, G4double where the value for R will be stored
//  - incidenceAngle, G4double where the value for the incidence angle will be stored,
//      in case of the calculation for an indirect ray that reflects on the ice-air boundary
//      or in case of air-to-ice ray tracing
//      -> 90 degrees corresponds to a horzontally going ray
//      FOR ICE-ICE RAY TRACING:
//      -> 0 degrees corresponds to a vertically upgoing ray, coming from the ice and hitting the 
//           ice-air boundary perpendicular to the boundary
//      FOR AIR-ICE RAY TRACING:
//      -> 0 degrees corresponds to a vertically downgoing ray, coming from the air and hitting the 
//           ice-air boundary perpendicular to the boundary
//  - focusingFact, G4double where the focusing factor will be stored,
//       which is limited between focFactUnderlim and focFactUpperlim
// Output:
// - 1 if ray was found, 0 if ray was not found
// NOTE:
// direction of launchDir and receiveDir are defined as if the emitter is shooting an arrow to
// the receiver, i.e. the launch dir points from emitter to receiver, while the receive dir 
// points from receiver AWAY from emitter
G4int IceShelfSteppingAction::GetEPValsRayTracing(unsigned int k, 
    G4ThreeVector emitterPos, G4ThreeVector emitterObsVect, G4double emitterTime, G4ThreeVector beta,
    G4int dirOrIndirRay, G4double &obsTime, G4ThreeVector &launchDir, G4ThreeVector &receiveDir,
    G4double &doppler, G4double &R, G4double &incidenceAngle, G4double &focusingFact){

    // Depth of the particle wrt the surface, < 0 when below surface
    G4double depthPart = emitterPos.getY() - shelfSizeY/2.;

    // Applying the offset in case depth of the particle is exactly 0, which the ray tracer
    // cannot handle
    if(depthPart == 0){
        depthPart = depthPart - zeroDepthOffset;
    }

    // Horizontal distance from particle to receiver
    G4double horDist = sqrt( emitterObsVect.getX()*emitterObsVect.getX()
                                                    + emitterObsVect.getZ()*emitterObsVect.getZ() );

    G4double travelTime;
    G4double geoPath;
    G4double receiveAngle;
    G4double launchAngle;
    G4double nAtEmitter;

    /////////////////////////////////////////////////////////////////////////////////////
    // IN CASE OF AIR TO ICE RAY TRACING (DOES NOT USE TABLES, AND NO FOCUSING FACTOR) //
    /////////////////////////////////////////////////////////////////////////////////////

    if(includeAirToIceRT && depthPart > 0){

        if(dirOrIndirRay == 1){
            return 0; // No indirect ray in case of air-ice ray tracing...
        }


        G4double depthAntenna = fRunAction->GetIceShelfAntennas()->GetAntennaCoord(k, 1) - shelfSizeY/2.;
        G4double iceShelfAltitude = fRunAction->GetIceShelfAntennas()->GetIceShelfAltitude();

        G4double optPathInIce;
        G4double optPathInAir;
        G4double geoPathInIce;
        G4double geoPathInAir;
        G4double horDistToInterPoint;

        // Doing the ray tracing
        G4bool foundSolution = false;
        foundSolution = AirToIceRayTracing::GetRayTracingSolution(
                     (iceShelfAltitude + depthPart)/cm, horDist/cm, depthAntenna/cm, iceShelfAltitude/cm,
                     optPathInIce, optPathInAir, geoPathInIce, geoPathInAir,
                     launchAngle, horDistToInterPoint, incidenceAngle, receiveAngle); 
        if(!foundSolution){
            return 0;
        }

        // Switching to internal units
        optPathInIce = optPathInIce*cm;
        optPathInAir = optPathInAir*cm;
        geoPathInIce = geoPathInIce*cm;
        geoPathInAir = geoPathInAir*cm;
        launchAngle = launchAngle*degree;
        horDistToInterPoint = horDistToInterPoint*cm;
        incidenceAngle = incidenceAngle*degree;
        receiveAngle = receiveAngle*degree;

        // Receive angle should be the angle between receive direction (as defined above) and the
        // vertical.
        receiveAngle = 180.*degree - receiveAngle;

        // Calculating the travel time
        travelTime = optPathInIce/c0 + optPathInAir/c0;

        // Calculating the geometric path length
        geoPath = geoPathInIce + geoPathInAir;

        // Calculating the index of refraction at the emitter
        G4double emitterAltitude = iceShelfAltitude + depthPart;
        nAtEmitter = AirToIceRayTracing::Getnz_air(emitterAltitude/m);

        // Setting the focusing factor to 1, which means applying it later will do nothing...
        focusingFact = 1.;

    }

    /////////////////////////////////////////////////
    // IN CASE OF ICE TO ICE TRACING (USES TABLES) //
    /////////////////////////////////////////////////

    else {

        // We use the same definitions for launch and receive angle as Uzair, where the emitter
        // is the charged particle and the receiver is the antenna:
        // - LAUNCH ANGLE: angle between y-direction (vertical upwards) and emitter direction as 
        //      defined above
        // - RECEIVE ANGLE: angle between y-direction (vertical upwards) and receiver direction as 
        //      defined above

        // Getting the ray tracing values
        // travelTime: time it takes to go from emitter to observer (internal units)
        // geoPath: geometric path length (internal units)
        // receiveAngle: receive angle as defined above (internal units)
        // launchAngle: launch angle as defined above (internal units)
        // optPath: optical path length (internal units)

        if(dirOrIndirRay == 0){

            travelTime = IceRayTracing::GetInterpolatedValue(horDist/m, depthPart/m, 0, k)*s;
            geoPath = IceRayTracing::GetInterpolatedValue(horDist/m, depthPart/m, 1, k)*m;
            launchAngle = IceRayTracing::GetInterpolatedValue(horDist/m, depthPart/m, 2, k)*degree;
            receiveAngle = IceRayTracing::GetInterpolatedValue(horDist/m, depthPart/m, 3, k)*degree;
            incidenceAngle = -1000.*degree; // No incidence angle
            focusingFact = IceRayTracing::GetInterpolatedValue(horDist/m, depthPart/m, 5, k);

        } else if(dirOrIndirRay == 1){

            travelTime = IceRayTracing::GetInterpolatedValue(horDist/m, depthPart/m, 6, k)*s;
            geoPath = IceRayTracing::GetInterpolatedValue(horDist/m, depthPart/m, 7, k)*m;
            launchAngle = IceRayTracing::GetInterpolatedValue(horDist/m, depthPart/m, 8, k)*degree;
            receiveAngle = IceRayTracing::GetInterpolatedValue(horDist/m, depthPart/m, 9, k)*degree;
            incidenceAngle = IceRayTracing::GetInterpolatedValue(horDist/m, depthPart/m, 12, k)*degree;
            focusingFact = IceRayTracing::GetInterpolatedValue(horDist/m, depthPart/m, 11, k);

        } else{

            G4cerr << "The value of dirOrIndirRay should be 0 or 1, instead it is "
                << dirOrIndirRay << G4endl;
            exit(EXIT_FAILURE);

        }

        //G4double optPath = travelTime*c0;

        // Checking whether a ray was found through the default fail values
        if(travelTime == -1000*s || geoPath == -1000*m || receiveAngle == -1000.*degree
                || launchAngle == -1000.*degree){
            return 0;
        }

        // Checking whether a ray was found by checking for nans
        if(std::isnan(travelTime) || std::isnan(geoPath) || std::isnan(receiveAngle)
                                                                || std::isnan(launchAngle)){
            G4cerr << "WARNING: ecountered a 'nan' value for one of the ray tracing values" << G4endl;
            G4cerr << "travelTime = " << travelTime/ns << " ns, " << "geoPath = " << geoPath/m << " m, "
                << "receiveAngle = " << receiveAngle/degree << " degr, "
                << "launchAngle = " << launchAngle/degree << " degr" << G4endl;
            G4cerr << "horDist = " << horDist/m << " m, " << "depthPart = " << depthPart/m << " m, "
                << "antenna number is " << k << G4endl;
            G4cerr << "Probably no suitable ray was found." << G4endl;
            return 0;
        }

        // Checking if the values of the launch angle make sense
        if(launchAngle/degree < 0 || launchAngle/degree > 180){
            G4cerr << "ERROR: encountered a launchAngle of " << launchAngle/degree << " degr, "
                << "which is not possible..." << G4endl;
            exit(EXIT_FAILURE);
            return 0;
        }

        // Checking if the values for the travel time and geometric path length make sense
        if(travelTime < 0. || geoPath < 0.){
            G4cerr << "ERROR: Found value for travelTime of " << travelTime/ns << " ns, and a value for "
                << "geoPath of " << geoPath/m << " m. Not sure what this means..." << G4endl;
            exit(EXIT_FAILURE);
            return 0;
        }

        // Applying the focusing factor limits
        if(focusingFact < focFactUnderlim){
            focusingFact = focFactUnderlim;
        } else if(focusingFact > focFactUpperlim){
            focusingFact = focFactUpperlim;
        }

        // The index of refraction at the emitter
        nAtEmitter = IceRayTracing::Getnz(depthPart/m);

    }

    // Calculating the observation time
    obsTime = emitterTime + travelTime;

    // Calculating the launch and receive direction of the ray
    G4ThreeVector yDir = G4ThreeVector(0., 1., 0.); // Unit vector in y-dir (internal units)
    G4double emitterObsVectMag = emitterObsVect.mag();
    G4ThreeVector rotAxis = yDir.cross(emitterObsVect/emitterObsVectMag); // Rotation axis to get
                                                                          // from yDir to launchDir
    G4double rotAxisMag = rotAxis.mag();
    // Sanity check, if emitterObsVect and yDir are near parallel, the direction of the
    // rotation axis is maybe dominated by numerical noise
    // If so, we set the directions manually
    if(rotAxisMag < 1e-12){

        G4cerr << "WARNING: encoutered rotation axis with a magnitude of "
                << rotAxisMag
                << " during calculation of launch and receive direction of a ray, " 
                << " which means emitter and observer are on a near-vertical line."
                << " Rotation axis is normalized, and its direction might"
                << " be dominated by numerical noise." << G4endl;
        //G4cerr << "launchAngle and receiveAngle used during rotation with this axis were " 
        //        << launchAngle/degree << " degr and " << receiveAngle/degree << " degr." << G4endl;
        G4cerr << "Setting directions manually (vert. up if angle < 2 dgr, vert. down if"
                << " angle > 178 dgr)..." << G4endl;
        G4cerr << "launchAngle is " << launchAngle/degree << " degr, receiveAngle is "
                << receiveAngle/degree << " degr" << G4endl;

        if(launchAngle/degree < 2.){
            launchDir = yDir;
        } else if(launchAngle/degree > 178.){
            launchDir = -1.*yDir;
        } else {
            G4cerr << "ERROR: launchAngle is " << launchAngle/degree << " degr, which does not"
            << " satisfy any of the two conditions given above. Not sure why rotation axis"
            << " has a magnitude close to 0. Aborting..." << G4endl;
            exit(EXIT_FAILURE);
        }

        if(receiveAngle/degree < 2.){
            receiveDir = yDir;
        } else if(receiveAngle/degree > 178.){
            receiveDir = -1.*yDir;
        } else {
            G4cerr << "ERROR: receiveAngle is " << receiveAngle/degree << " degr, which does not"
            << " satisfy any of the two conditions given above. Not sure why rotation axis"
            << " has a magnitude close to 0. Aborting..." << G4endl;
            exit(EXIT_FAILURE);
        }

        G4cerr << "launchDir has been set to " << launchDir 
               << ", receiveDir has been set to " << receiveDir << G4endl;

    } else {
        rotAxis = rotAxis/rotAxisMag;
        launchDir = RotateVectorCosPi(yDir, rotAxis, cos(launchAngle/radian));
        receiveDir = RotateVectorCosPi(yDir, rotAxis, cos(receiveAngle/radian));
    }

    // Calculating the doppler term (1 - n*beta . rhat)
    // with n the index of refraction at the emission point,
    // beta the velocity vector of the emitter divided by c0
    // and rhat the launch direction of the ray
    doppler = 1. - nAtEmitter * (beta.dot(launchDir));
    if(doppler == 0){
        G4cerr << "WARNING: encountered a doppler term = 0" << G4endl;
    }

    // Setting the value for R
    R = geoPath;

    return 1; 

}

// Calculates the electric field associated with a start point, end point and receiver antenna
// using the End-Point formalism.
// Input:
//  - k, unsigned int giving the antenna number
//  - observerPos, the G4ThreeVector giving the position of the antenna
//  - startPoint, the G4ThreeVector giving the start point of the step
//  - endPoint, the G4ThreeVector giving the end point of the step
//  - startTime, the G4double giving the time the emitter is in startPoint
//  - endTime, the G4double giving the time the emitter is in endPoint
//  - beta, the G4ThreeVector giving the beta vector (v/c) of the emitter
//  - particleCharge, a G4double giving the charge of the emitter (internal units)
//  - dirOrIndirRay, should be 0 if you want to calculate the values for the direct ray,
//          should be 1 if you want to calculate the values for the indirect ray,
//          only applicable in case of ray tracing
//  #########
//  - startE, G4ThreeVector where the contribution to the trace of the start point will be stored
//  - endE, G4ThreeVector where the contribution to the trace of the end point will be stored
//  - iStart, G4int where the time index for the trace for the start point will be stored in
//  - iEnd, G4int where the time index for the trace for the end point will be stored in
// Output:
// - 1 if electric field was calculated, 0 if not
//      (e.g. because no suitable ray was found in case of ray tracing)
G4int IceShelfSteppingAction::GetEfieldEP(unsigned int k, G4ThreeVector observerPos,
    G4ThreeVector startPoint, G4ThreeVector endPoint, G4double startTime, G4double endTime,
    G4ThreeVector beta, G4double particleCharge, G4int dirOrIndirRay,
    G4ThreeVector &startE, G4ThreeVector &endE, G4int &iStart, G4int &iEnd){

    // Calculating the vectors pointing from Start/End-point to observer
    G4ThreeVector startObsVect = observerPos - startPoint;
    G4ThreeVector endObsVect = observerPos - endPoint;

    // Getting the values needed for the End-Point formula, both for start and end point

    G4double startObsTime;
    G4ThreeVector startLaunchDir;
    G4ThreeVector startReceiveDir;
    G4double startDoppler;
    G4double startR;
    G4double startIncidenceAngle;
    G4double startFocusingFact;

    G4double endObsTime;
    G4ThreeVector endLaunchDir;
    G4ThreeVector endReceiveDir;
    G4double endDoppler;
    G4double endR;
    G4double endIncidenceAngle;
    G4double endFocusingFact;

    G4int startValsFound;
    G4int endValsFound;

    if(rayTracingOn){
        startValsFound = GetEPValsRayTracing(k, startPoint, startObsVect, startTime, beta, dirOrIndirRay,
               startObsTime, startLaunchDir, startReceiveDir, startDoppler, startR, startIncidenceAngle,
                                                                                    startFocusingFact);
        endValsFound = GetEPValsRayTracing(k, endPoint, endObsVect, endTime, beta, dirOrIndirRay,
                           endObsTime, endLaunchDir, endReceiveDir, endDoppler, endR, endIncidenceAngle,
                                                                                      endFocusingFact);
    } else if(dirOrIndirRay == 0){ // if no ray tracing, so don't do this for indirect ray again
        startValsFound = GetEPVals(startObsVect, startTime, beta,
                                    startObsTime, startLaunchDir, startReceiveDir, startDoppler, startR);
        endValsFound = GetEPVals(endObsVect, endTime, beta,
                                    endObsTime, endLaunchDir, endReceiveDir, endDoppler, endR);
    } else if(dirOrIndirRay == 1){ // If looking for indirect ray without ray tracing, do nothing
        startValsFound = 0;
        endValsFound = 0;
    } else {
        G4cerr << "ERROR: dirOrIndirRay is " << dirOrIndirRay << ", but should be 0 or 1." << G4endl;
        exit(EXIT_FAILURE);
    }
    if(!startValsFound || !endValsFound){
        return 0;
    }

    // We will calculate the emission due to creation and destruction of the charge
    // using formula (8) of arXiv:1007.4146v3 (Physical Review E 84, 056602 (2011))
    // We explicitly take care of the units here, because this formula is only valid
    // for the electrostatic CGS unit system

    G4double sPeriodCGS = sPeriod/s;
    G4double chargeCGS = particleCharge/franklin;
    G4double c0CGS = c0/(cm/s);

    G4double startRCGS = startR/cm;
    G4double endRCGS = endR/cm;

    // Initiating the trace indices we'll use to store the signal in the array.
    // We add the 0.5 so that e.g. both signal received after 5.5*sPeriod and before
    // 6.5*sPeriod will contribute to the value of the electrical field at timepoint 'array[6]'

    iStart = static_cast<int>(floor(startObsTime/sPeriod + 0.5l));
    iEnd = static_cast<int>(floor(endObsTime/sPeriod + 0.5l));

    // Declaring some variables we'll need later on.
    G4double deltaT = 0.;

    // First we deal with the case where we need to use the approximation
    // Defining these variables here, in case we need them later
    G4double interObsTime;
    G4ThreeVector interLaunchDir;
    G4ThreeVector interReceiveDir;
    G4double interDoppler;
    G4double interR;
    G4double interIncidenceAngle;
    G4double interFocusingFact;

    if(fabs(startDoppler) < approxThreshold || fabs(endDoppler) < approxThreshold){

        // The point between startPoint and endPoint
        G4ThreeVector interPoint = 0.5*(startPoint + endPoint);
        G4double interTime = 0.5*(startTime + endTime);

        // Calculating the vector pointing from inter-point to observer
        G4ThreeVector interObsVect = observerPos - interPoint;

        // Getting the values needed for the End-Point formula, now for the inter-point

        G4int interValsFound;

        if(rayTracingOn){
            interValsFound = GetEPValsRayTracing(k, interPoint, interObsVect, interTime, beta,
                    dirOrIndirRay, interObsTime, interLaunchDir, interReceiveDir,
                    interDoppler, interR, interIncidenceAngle, interFocusingFact);
        } else {
            interValsFound = GetEPVals(interObsVect, interTime, beta,
                                    interObsTime, interLaunchDir, interReceiveDir, interDoppler, interR);
        }
        if(!interValsFound){
            return 0;
        }

        // Changing to CGS system
        G4double interRCGS = interR/cm;

        // Calculating the electric field values at the antenna postion
        startE = (1./sPeriodCGS) * ((chargeCGS)/(c0CGS)) *
         ((interLaunchDir.cross(interLaunchDir.cross(beta))) / 
                                                        (fabs(interDoppler)*interRCGS));
        endE = -1.*startE;

        // Applying the focusing factor
        if(rayTracingOn){
            startE = interFocusingFact*startE;
            endE = interFocusingFact*endE;
        }

        // Make sure the time interval deltaT in which the signal 
        // is received scales with the Doppler factor,
        // i.e. the formula dt = dt'(1 - n beta.r) is valid
        // This will effect the values of startObsTime and endObsTime
        // NOTE: why the fabs?
        deltaT = (endTime - startTime)*fabs(interDoppler);
        if(startObsTime < endObsTime){ // startE arrives first
            startObsTime = interObsTime - 0.5*deltaT;
            endObsTime = interObsTime + 0.5*deltaT;
        } else { // endE arrives first
            startObsTime = interObsTime + 0.5*deltaT;
            endObsTime = interObsTime - 0.5*deltaT;
        }
        deltaT = endObsTime - startObsTime; // Can be positive of negative

        // Redistribute contributions over time scale defined by sampling rate of observer
        if(fabs(deltaT) < sPeriod){

            startE = startE*fabs(deltaT/sPeriod);
            endE = endE*fabs(deltaT/sPeriod);

            long double startBinFraction = (startObsTime/sPeriod) - floor(startObsTime/sPeriod);
            long double endBinFraction = (endObsTime/sPeriod) - floor(endObsTime/sPeriod);

            // Recalculate, because the times have changed a little
            iStart = static_cast<int>(floor(startObsTime/sPeriod + 0.5l));
            iEnd = static_cast<int>(floor(endObsTime/sPeriod + 0.5l));

            // Do timing modification if both contributions go to same point in time, 
            // which would mean they would cancel eachother...
            // Splitting them up follows more realistically the actual expected signal 
            // close to the Cherenkov angle

            if(iStart == iEnd){

                // If startE arrives before endE
                if(deltaT >= 0){
                    // Both fall left of the point in time
                    if((startBinFraction >= 0.5) && (endBinFraction >= 0.5)){
                        startObsTime = startObsTime - sPeriod; // Moving startE to previous point in time
                    }
                    // Both fall right of the point in time
                    else if((startBinFraction < 0.5) && (endBinFraction < 0.5)){
                        endObsTime = endObsTime + sPeriod; // Moving endE to the next point in time
                    }
                    // startE is left of point in time, endE is right of point in time
                    else{
                        // If endE is closer to the next point in time than startE is to previous
                        // point in time, we move endE to the next point in time
                        if(endBinFraction >= (1.0 - startBinFraction)){
                            endObsTime = endObsTime + sPeriod;
                        // Otherwise we move startE to the previous point in time
                        } else {
                            startObsTime = startObsTime - sPeriod;
                        }
                    }
                }
                // If endE arrives before startE
                else{
                    // Both fall left of the point in time
                    if((startBinFraction >= 0.5) && (endBinFraction >= 0.5)){
                        endObsTime = endObsTime - sPeriod; // Moving endE to previous point in time
                    }
                    // Both fall right of the point in time
                    else if((startBinFraction < 0.5) && (endBinFraction < 0.5)){
                        startObsTime = startObsTime + sPeriod; // Moving startE to next point in time
                    }
                    // startE is right of point in time, endE is left of point in time
                    else{
                        // If endE is closest to the point in time, we move startE
                        // to the next point in time
                        if(startBinFraction >= (1.0 - endBinFraction)){
                            startObsTime = startObsTime + sPeriod;
                        // Otherwise we move endE to the previous point in time
                        } else{
                            endObsTime = endObsTime - sPeriod;
                        }
                    }
                }

                iStart = static_cast<int>(floor(startObsTime/sPeriod + 0.5l));
                iEnd = static_cast<int>(floor(endObsTime/sPeriod + 0.5l));

            }
        }
    }

    // Now dealing with the case where we do not need to use the approximation
    else{

        startE = (1./sPeriodCGS) * ((chargeCGS)/(c0CGS)) *
         ((startLaunchDir.cross(startLaunchDir.cross(beta))) / 
                                                    (fabs(startDoppler)*startRCGS));
        endE = (-1.0) * (1./sPeriodCGS) * ((chargeCGS)/(c0CGS)) *
         ((endLaunchDir.cross(endLaunchDir.cross(beta))) / 
                                                    (fabs(endDoppler)*endRCGS));

        // Applying the focusing factors
        if(rayTracingOn){
            startE = startFocusingFact*startE;
            endE = endFocusingFact*endE;
        }

    }

    // Doing a rotation, so the electric field is perpendicular to the ray
    // Only needed in case of ray tracing, startE and endE are by construction already
    // perpendicular to the straight path
    // In case of ray tracing, it will be perpendicular to the launch direction, and this has
    // to be changed to be perpendicular to the receive direction
    if(rayTracingOn){

        G4ThreeVector yDir = G4ThreeVector(0., 1., 0.); // Unit vector in y-dir (internal units)
        G4ThreeVector rotAxis;
        G4double cosVal;

        // Checking if we did the approximation

        if(fabs(startDoppler) < approxThreshold || fabs(endDoppler) < approxThreshold){

            cosVal = interLaunchDir.dot(interReceiveDir);
            if(cosVal > 1.){// Can happen due to numerical issues
                cosVal = 1.;
            }

            rotAxis = yDir.cross(interLaunchDir);
            G4double rotAxisMag = rotAxis.mag();

            // Sanity check, if rotAxisMag was very small the direction of rotation might be
            // dominated by numerical noise...
            // In this case, launch dir and normally also receive dir should be close to vertical
            if(rotAxisMag < 1e-12){

                G4cerr << "WARNING: encoutered rotation axis with a magnitude of " << rotAxisMag
                        << " during rotation of electric field with interPoint,"
                        << " which means launch direction is"
                        << " near vertical. Rotation axis is normalized, and its direction might"
                        << " be dominated by numerical noise." << G4endl;
                //G4cerr << "cos(theta) during rotation with this axis was " << cosVal << G4endl;
                G4cerr << "Doing rotation of electric fields manually: no rotation if "
                << "launchDir.dot(receiveDir) > .998, sign flip of field if < -0.998" << G4endl;

                if(cosVal < -0.998){
                    G4cerr << "cosVal is " << cosVal << ", so flipping the signs" << G4endl;
                    startE = -1.*startE;
                    endE = -1.*endE;
                } else if(cosVal > 0.998){
                    G4cerr << "cosVal is " << cosVal << ", so no rotations" << G4endl;
                } else {
                    G4cerr << "ERROR: cosVal is " << cosVal << ", which does not fall "
                    << "under any of the two conditions given above. Not sure why the magnitude"
                    << " of the rotation axis is so small. Aborting..." << G4endl;
                    exit(EXIT_FAILURE);
                }

            } else {

                rotAxis = rotAxis/rotAxisMag;
                // Sanity check, cos(theta_launch) < cos(theta_receive) shouldn't happen,
                // but maybe it does for double exp profile
                // If it does, the sign of the rotation axis should be inverted
                if(yDir.dot(interLaunchDir) < yDir.dot(interReceiveDir)){
                    //G4cerr << "WARNING: encountered launch angle larger than receive angle"
                    //<< G4endl;
                    //G4cerr << "cos(theta_launch) is " << yDir.dot(interLaunchDir)
                    //       << " and cos(theta_receive) is "
                    //       << yDir.dot(interReceiveDir) << G4endl;
                    rotAxis = -1.*rotAxis; // Changing the direction of rotation
                }
                startE = RotateVectorCosPi(startE, rotAxis, cosVal);
                endE = RotateVectorCosPi(endE, rotAxis, cosVal);

            }

        } else {

            // START POINT

            cosVal = startLaunchDir.dot(startReceiveDir);
            if(cosVal > 1.){// Can happen due to numerical issues
                cosVal = 1.;
            }

            rotAxis = yDir.cross(startLaunchDir);
            G4double rotAxisMag = rotAxis.mag();
            // Sanity check, if rotAxisMag was very small the direction of rotation might be
            // dominated by numerical noise...
            // In this case, launch dir and normally also receive dir should be close to vertical
            if(rotAxisMag < 1e-12){

                G4cerr << "WARNING: encoutered rotation axis with a magnitude of " << rotAxisMag
                        << " during rotation of electric field with startPoint,"
                        << " which means launch direction is"
                        << " near vertical. Rotation axis is normalized, and its direction might"
                        << " be dominated by numerical noise." << G4endl;
                //G4cerr << "cos(theta) during rotation with this axis was " << cosVal << G4endl;
                G4cerr << "Doing rotation of electric fields manually: no rotation if "
                << "launchDir.dot(receiveDir) > .998, sign flip of field if < -0.998" << G4endl;

                if(cosVal < -0.998){
                    G4cerr << "cosVal is " << cosVal << ", so flipping the sign" << G4endl;
                    startE = -1.*startE;
                } else if(cosVal > 0.998){
                    G4cerr << "cosVal is " << cosVal << ", so no rotation" << G4endl;
                } else {
                    G4cerr << "ERROR: cosVal is " << cosVal << ", which does not fall "
                    << "under any of the two conditions given above. Not sure why the magnitude"
                    << " of the rotation axis is so small. Aborting..." << G4endl;
                    exit(EXIT_FAILURE);
                }

            } else {

                rotAxis = rotAxis/rotAxisMag;
                // Sanity check, cos(theta_launch) < cos(theta_receive) shouldn't happen,
                // but maybe it does for double exp profile
                // If it does, the sign of the rotation axis should be inverted
                if(yDir.dot(startLaunchDir) < yDir.dot(startReceiveDir)){
                    //G4cerr << "WARNING: encountered launch angle larger than receive angle"
                    //<< G4endl;
                    //G4cerr << "cos(theta_launch) is " << yDir.dot(startLaunchDir)
                    //       << " and cos(theta_receive) is "
                    //       << yDir.dot(startReceiveDir) << G4endl;
                    rotAxis = -1.*rotAxis; // Changing the direction of rotation
                }
                startE = RotateVectorCosPi(startE, rotAxis, cosVal);

            }

            // END POINT

            cosVal = endLaunchDir.dot(endReceiveDir);
            if(cosVal > 1.){// Can happen due to numerical issues
                cosVal = 1.;
            }

            rotAxis = yDir.cross(endLaunchDir);
            rotAxisMag = rotAxis.mag();
            // Sanity check, if rotAxisMag was very small the direction of rotation might be
            // dominated by numerical noise...
            // In this case, launch dir and normally also receive dir should be close to vertical
            if(rotAxisMag < 1e-12){

                G4cerr << "WARNING: encoutered rotation axis with a magnitude of " << rotAxisMag
                        << " during rotation of electric field with endPoint,"
                        << " which means launch direction is"
                        << " near vertical. Rotation axis is normalized, and its direction might"
                        << " be dominated by numerical noise." << G4endl;
                //G4cerr << "cos(theta) during rotation with this axis was " << cosVal << G4endl;
                G4cerr << "Doing rotation of electric fields manually: no rotation if "
                << "launchDir.dot(receiveDir) > .998, sign flip of field if < -0.998" << G4endl;

                if(cosVal < -0.998){
                    G4cerr << "cosVal is " << cosVal << ", so flipping the sign" << G4endl;
                    endE = -1.*endE;
                } else if(cosVal > 0.998){
                    G4cerr << "cosVal is " << cosVal << ", so no rotation" << G4endl;
                } else {
                    G4cerr << "ERROR: cosVal is " << cosVal << ", which does not fall "
                    << "under any of the two conditions given above. Not sure why the magnitude"
                    << " of the rotation axis is so small. Aborting..." << G4endl;
                    exit(EXIT_FAILURE);
                }

            } else {

                rotAxis = rotAxis/rotAxisMag;
                // Sanity check, cos(theta_launch) < cos(theta_receive) shouldn't happen,
                // but maybe it does for double exp profile
                // If it does, the sign of the rotation axis should be inverted
                if(yDir.dot(endLaunchDir) < yDir.dot(endReceiveDir)){
                    //G4cerr << "WARNING: encountered launch angle larger than receive angle"
                    //<< G4endl;
                    //G4cerr << "cos(theta_launch) is " << yDir.dot(endLaunchDir)
                    //       << " and cos(theta_receive) is "
                    //       << yDir.dot(endReceiveDir) << G4endl;
                    rotAxis = -1.*rotAxis; // Changing the direction of rotation
                }
                endE = RotateVectorCosPi(endE, rotAxis, cosVal);

            }

        } 
        
    }

    // Taking into account the Fresnel coefficients
    if(rayTracingOn){

        G4double startDepth = startPoint.getY() - shelfSizeY/2.; // > 0 if above the ice
        G4double endDepth = endPoint.getY() - shelfSizeY/2.; // > 0 if above the ice

        // Checking if we did the approximation
        if(fabs(startDoppler) < approxThreshold || fabs(endDoppler) < approxThreshold){

            G4double interDepth = (startDepth + endDepth)/2.;

            // Checking if we're dealing with a reflected ray, coming from ice
            if(interDepth <= 0 && dirOrIndirRay == 1){
                // Checking if we're dealing with a reflected ray, not a refracted
                if(!std::isnan(interIncidenceAngle) && interIncidenceAngle != -1000.*degree){
                    ApplyFresnelCoeff(interIncidenceAngle, interReceiveDir, startE, 0);
                    ApplyFresnelCoeff(interIncidenceAngle, interReceiveDir, endE, 0);
                }
            }

            // Checking if we're dealing with a transmitted ray, coming from air
            if(includeAirToIceRT && interDepth > 0 && dirOrIndirRay == 0){
                ApplyFresnelCoeff(interIncidenceAngle, interReceiveDir, startE, 1);
                ApplyFresnelCoeff(interIncidenceAngle, interReceiveDir, endE, 1);
            }

        } else { // No approximation

            // Checking if we're dealing with reflected rays, coming from the ice
            if(dirOrIndirRay == 1){
                if(startDepth <= 0){
                    // Reflected, or refracted?
                    if(!std::isnan(startIncidenceAngle) && startIncidenceAngle != -1000.*degree){
                        ApplyFresnelCoeff(startIncidenceAngle, startReceiveDir, startE, 0);
                    }
                }
                if(endDepth <= 0){
                    // Reflected, or refracted?
                    if(!std::isnan(endIncidenceAngle) && endIncidenceAngle != -1000.*degree){
                        ApplyFresnelCoeff(endIncidenceAngle, endReceiveDir, endE, 0);
                    }
                }
            }

            // Checking if we're dealing with a transmitted ray, coming from air
            if(includeAirToIceRT){
                if(dirOrIndirRay == 0){
                    if(startDepth > 0){
                        ApplyFresnelCoeff(startIncidenceAngle, startReceiveDir, startE, 1);
                    }
                    if(endDepth > 0){
                        ApplyFresnelCoeff(endIncidenceAngle, endReceiveDir, endE, 1);
                    }
                }
            }

        }

    }

    return 1;

}

// Calculates the reflection coeffs for the s-polarisation (phi-polarisation)
// and p-polarisation (theta-polarisation), for a ray emitted IN ICE hitting the air-ice boundary
// Input:
//  - incidenceAngle, G4double giving the incidence angle on the ice-air boundary
//          -> 0 degrees corresponds to a vertically upgoing ray, coming from the ice and hitting the 
//             ice-air boundary perpendicular to the boundary
//          -> 90 degrees corresponds to a horzontally going ray
//  - nIceAB, G4double giving the index of refraction of the ice at the ice-air boundary
//  - nAirAB, G4double giving the index of refraction of the air at the ice-air boundary
//  #########
//  - rs, G4double where the reflection coeffs for the s-polarisation will be stored
//  - rp, G4double where the reflection coeffs for the p-polarisation will be stored
void IceShelfSteppingAction::GetReflectionCoeffs(G4double incidenceAngle, G4double nIceAB,
                                                 G4double nAirAB, G4double &rs, G4double &rp){

    G4double cosIncidenceAngle = cos(incidenceAngle);
    G4double cosTransAngle = 0.;

    // Checking whether the incidence angle is smaller than the critical angle at which we get
    // total reflection. If not, the coeff values become complex, with modulus 1 and a certain
    // phase shift. We will not take the phase shift into account, instead using rs = 1 and rp = 1.
    G4double totReflectionAngle = asin(nAirAB/nIceAB)*radian;
    if(incidenceAngle < totReflectionAngle){
        cosTransAngle = sqrt(1. - pow((nIceAB/nAirAB)*sin(incidenceAngle),2)); // Fresnel's law
    }

    // The s-polarisation (phi-polarisation)
    rs = (nIceAB*cosIncidenceAngle - nAirAB*cosTransAngle)/
         (nIceAB*cosIncidenceAngle + nAirAB*cosTransAngle);

    // The p-polarisation (theta-polarisation)
    rp = (nAirAB*cosIncidenceAngle - nIceAB*cosTransAngle)/
         (nAirAB*cosIncidenceAngle + nIceAB*cosTransAngle);

}

// Calculates the transition coeffs for the s-polarisation (phi-polarisation)
// and p-polarisation (theta-polarisation), for a ray emitted IN AIR hitting the air-ice boundary
// Input:
//  - incidenceAngle, G4double giving the incidence angle on the ice-air boundary
//              -> 0 degrees corresponds to a vertically downgoing ray, coming from the air and hitting 
//                 the ice-air boundary perpendicular to the boundary
//              -> 90 degrees corresponds to a horizontally going ray
//  - nIceAB, G4double giving the index of refraction of the ice at the ice-air boundary
//  - nAirAB, G4double giving the index of refraction of the air at the ice-air boundary
//  #########
//  - ts, G4double where the transition coeffs for the s-polarisation will be stored
//  - ts, G4double where the transition coeffs for the p-polarisation will be stored
void IceShelfSteppingAction::GetTransitionCoeffs(G4double incidenceAngle, G4double nIceAB,
                                                 G4double nAirAB, G4double &ts, G4double &tp){

    G4double cosIncidenceAngle = cos(incidenceAngle);

    // Getting the transition angle, which should always be defined since nAirAB < nIceAB
    if(nIceAB <= nAirAB){
        G4cerr << "ERROR: seems like nIceAB (" << nIceAB << ") <= nAirAB (" << nAirAB
               << ") which means you are using the GetTransitionCoeffs function wrongly."
               << " Please read it's description in the cpp file." << G4endl;
        exit(EXIT_FAILURE);;
    }
    G4double cosTransAngle = sqrt(1. - pow((nAirAB/nIceAB)*sin(incidenceAngle),2)); // Fresnel's law

    // The s-polarisation (phi-polarisation)
    ts = 2*nAirAB*cosIncidenceAngle/(nAirAB*cosIncidenceAngle + nIceAB*cosTransAngle);

    // The p-polarisation (theta-polarisation)
    tp = 2*nAirAB*cosIncidenceAngle/(nIceAB*cosIncidenceAngle + nAirAB*cosTransAngle);

}

// Applies the Fresnel coefficients for the given Efield
// - for a ray emitted IN THE ICE and REFLECTING on the ice-air boundary
// - or for a ray emitted IN THE AIR and TRANSITIONING through the air-ice boundary
// Input:
//  - incidenceAngle, G4double giving the incidence angle on the ice-air boundary
//          FOR REFLECTION (EMISSION IN ICE):
//          -> 0 degrees corresponds to a vertically upgoing ray, coming from the ice and hitting the 
//             ice-air boundary perpendicular to the boundary
//          -> 90 degrees corresponds to a horzontally going ray
//          FOR TRANSITION (EMISSION IN AIR):
//          -> 0 degrees corresponds to a vertically downgoing ray, coming from the air and hitting the
//             air-ice boudanry perpendicular to the boundary
//          -> 90 degrees corresponds to a horzontally going ray
//  - receiveDir, G4ThreeVector giving the receive direction at the receiver
//  - Efield, G4ThreeVector to which the Fresnel coefficients will be applied
//    NOTE: Efield should be the electric field AT THE RECEIVER
//  - reflOrTrans, G4int that should be 0 for reflection and 1 for transition
void IceShelfSteppingAction::ApplyFresnelCoeff(G4double incidenceAngle, G4ThreeVector receiveDir,
                                                    G4ThreeVector &Efield, G4int reflOrTrans){

    G4double nIceAtBoundary = IceRayTracing::Getnz(0./m);

    G4double iceShelfAltitude = fRunAction->GetIceShelfAntennas()->GetIceShelfAltitude();
    G4double nAirAtBoundary = AirToIceRayTracing::Getnz_air(iceShelfAltitude/m);

    // Getting the values for s (phi) polarisation and p (theta) polarisation
    G4double sval;
    G4double pval;

    if(reflOrTrans == 0){
        GetReflectionCoeffs(incidenceAngle, nIceAtBoundary, nAirAtBoundary, sval, pval);
    } else if(reflOrTrans == 1){
        GetTransitionCoeffs(incidenceAngle, nIceAtBoundary, nAirAtBoundary, sval, pval);
    } else {
        G4cerr << "Got a reflOrTrans variable in ApplyFresnelCoeff function that has the value "
                << reflOrTrans << ", but it should be 0 or 1." << G4endl;
        exit(EXIT_FAILURE);
    }

    // Constructing the r, theta and phi unit vectors at the receiver
    G4ThreeVector r = receiveDir;
    G4double rx = r.getX();
    G4double rz = r.getZ();

    G4double phix = -1.*rz/sqrt(rz*rz + rx*rx);
    G4double phiy = 0.;
    G4double phiz = rx/sqrt(rz*rz + rx*rx);
    G4ThreeVector phi = G4ThreeVector(phix, phiy, phiz);

    G4ThreeVector theta = phi.cross(r);

    // Getting the r, theta and phi components of the Efield
    G4double Er = Efield.dot(r);
    G4double Ephi = Efield.dot(phi);
    G4double Etheta = Efield.dot(theta);

    // Applying the Rs and Rp factors
    Ephi = sval*Ephi;
    Etheta = pval*Etheta;

    // Reconstructing the Efield
    Efield = Er*r + Ephi*phi + Etheta*theta;

}

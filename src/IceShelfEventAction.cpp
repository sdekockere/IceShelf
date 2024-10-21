#include "IceShelfEventAction.hh"
#include "IceShelfAnalysis.hh"
#include "IceShelfAntennas.hh"

#include "G4Event.hh"
#include "G4UnitsTable.hh"

#include <chrono>

IceShelfEventAction::IceShelfEventAction(IceShelfRunAction* runAction)
 : G4UserEventAction(),
   fEventID(0),
   fStartingTime(0.),
   fPrimaryWeight(0.),
   fPrimaryCharge(0.),
   fPrimaryStartPos(G4ThreeVector(0., 0., 0.)),
   fPrimaryStartMomentum(G4ThreeVector(0., 0., 0.)),
   fPrimaryMass(0.),
   fLeakedEnergy(0.),
   fSteppingAction(0),
   fRunAction(runAction),
   startOfEventTime(0)
{
    fPrimaryPartName = "NONE";
}

IceShelfEventAction::~IceShelfEventAction()
{}

void IceShelfEventAction::BeginOfEventAction(const G4Event* event)
{
    // Do stuff. The weight, charge, time and name of the primary (i.e. the value of the fPrimaryWeight, 
    // fStartingTime, fPrimaryCharge, fPrimaryStartPos, fPrimaryStartMomentum, fPrimaryMass
    // and fPrimaryEventName field) 
    // is set by the IceShelfPrimaryGeneratorAction object, 
    // so don't change it here because that would mean you would lose it.

    fEventID = event->GetEventID();
    startOfEventTime = std::clock();

}

void IceShelfEventAction::EndOfEventAction(const G4Event* event)
{

    // Print per event
    auto eventID = event->GetEventID();
    G4cout << "----> End of event: " << eventID << G4endl;
    G4cout << "Particle name of this event: " << fPrimaryPartName << G4endl;
    G4cout << "Charge of the particle (unit e+): " << fPrimaryCharge/eplus << G4endl;
    G4cout << "Weight of this event: " << fPrimaryWeight << G4endl;
    G4cout << "Energy leaked out of the simulated volume into the world volume: " << G4BestUnit(fLeakedEnergy, "Energy") << G4endl;

    // Set event id, particle id, ... to 0 to prepare for the next event.
    // Not really necessary, but it's cleaner this way...
    fEventID = 0;
    fPrimaryPartName = "NONE";
    fStartingTime = 0.;
    fPrimaryWeight = 0.;
    fPrimaryCharge = 0.;
    fPrimaryStartPos = G4ThreeVector(0., 0., 0.);
    fPrimaryStartMomentum = G4ThreeVector(0., 0., 0.);
    fPrimaryMass = 0.;
    fLeakedEnergy = 0.;

    // Resetting fCurrentTrackID etc
    fSteppingAction->ResetCurrentTrackID();
    fSteppingAction->ResetCurrentParticleDepthIndex();
    fSteppingAction->ResetCurrentParticleChargeIndex();
    fSteppingAction->ResetCurrentParticleSnapshotTimeIndex();

    clock_t endOfEventTime = std::clock(); //
    G4double eventTime = double(endOfEventTime - startOfEventTime)/double(CLOCKS_PER_SEC); //
    G4cout << "Time taken for event to complete: " << eventTime << " s" << G4endl; //

    startOfEventTime = 0;
}

/*
void IceShelfEventAction::CalcTransitionRadiationEndpointOld(){

    if(fPrimaryCharge == 0){
        return;
    }

    // The magnitude of the momentum (p*c) of the primary
    G4double primaryMomMag = fPrimaryStartMomentum.mag();

    // Getting the arrival direction and velocity (beta) of the primary particle
    G4ThreeVector arrivalDir = fPrimaryStartMomentum/primaryMomMag;
    G4double betaMag = primaryMomMag/sqrt(fPrimaryMass*fPrimaryMass + primaryMomMag*primaryMomMag);
    G4ThreeVector beta = betaMag*arrivalDir;

    // Getting the start and end point for the air->boundary step, with corresponding times
    G4ThreeVector airStartPoint = fPrimaryStartPos - transRadStepsize*arrivalDir;
    G4ThreeVector airEndPoint = fPrimaryStartPos -  transRadStepFrac*transRadStepsize*arrivalDir;
    G4double airStartTime = fStartingTime - transRadStepsize/(betaMag*c0);
    G4double airEndTime = fStartingTime - transRadStepFrac*transRadStepsize/(betaMag*c0);

    // Getting the start and end point for the boundary->ice step
    G4ThreeVector iceStartPoint = fPrimaryStartPos + transRadStepFrac*transRadStepsize*arrivalDir;
    G4ThreeVector iceEndPoint = fPrimaryStartPos + transRadStepsize*arrivalDir;
    G4double iceStartTime = fStartingTime + transRadStepFrac*transRadStepsize/(betaMag*c0);
    G4double iceEndTime = fStartingTime + transRadStepsize/(betaMag*c0);

    //G4cout << "This is the event spreaking: " << G4endl;
    //G4cout << "The momentum (p*c) of the primary is " << primaryMomMag/GeV << " GeV" << G4endl;
    //G4cout << "The mass (m*c^2) of the primary is " << fPrimaryMass/MeV << " MeV" << G4endl;
    //G4cout << "arrivalDir of the primary is " << arrivalDir << G4endl;
    //G4cout << "The beta value of the primary is " << beta << G4endl;
    //G4cout << "The airStartPoint is at " << airStartPoint/m << " m" << G4endl;
    //G4cout << "The airEndPoint is at " << airEndPoint/m << " m" << G4endl;
    //G4cout << "The airStartTime is " << airStartTime/ns << " ns" << G4endl;
    //G4cout << "The airEndTime is " << airEndTime/ns << " ns" << G4endl;
    //G4cout << "The iceStartPoint is at " << iceStartPoint/m << " m" << G4endl;
    //G4cout << "The iceEndPoint is at " << iceEndPoint/m << " m" << G4endl;
    //G4cout << "The iceStartTime is " << iceStartTime/ns << " ns" << G4endl;
    //G4cout << "The iceEndTime is " << iceEndTime/ns << " ns" << G4endl;

    for(G4int k = 0; k < fRunAction->GetIceShelfAntennas()->GetNumberOfAntennas(); k++){

        // Air emission
        G4ThreeVector airStartE;
        G4ThreeVector airEndE;
        G4int iStartAir;
        G4int iEndAir;
        G4int gotContributionAir;

        // Direct ice emission
        G4ThreeVector iceStartEDir;
        G4ThreeVector iceEndEDir;
        G4int iStartIceDir;
        G4int iEndIceDir;
        G4int gotContributionIceDir;

        // Indirect ice emission
        G4ThreeVector iceStartEIndir;
        G4ThreeVector iceEndEIndir;
        G4int iStartIceIndir;
        G4int iEndIceIndir;
        G4int gotContributionIceIndir;

        // The antenna position
        G4ThreeVector antennaPosition = G4ThreeVector(
                                          fRunAction->GetIceShelfAntennas()->GetAntennaCoord(k, 0),
                                          fRunAction->GetIceShelfAntennas()->GetAntennaCoord(k, 1),
                                          fRunAction->GetIceShelfAntennas()->GetAntennaCoord(k, 2));

        // First getting the contribution for the step air -> boundary, of which we will only
        // keep the start point emission. Direct ray only.

        gotContributionAir = fSteppingAction->GetEfieldEP(k, antennaPosition,
                                         airStartPoint, airEndPoint, airStartTime, airEndTime,
                                         beta, fPrimaryCharge, 0,
                                         airStartE, airEndE, iStartAir, iEndAir);

        if(gotContributionAir){

            // Taking the weight of the primary into account
            airStartE = fPrimaryWeight*airStartE;
            airEndE = fPrimaryWeight*airEndE;

            // Checking if the time window for registration is large enough
            if(iStartAir < 0 || iStartAir >= static_cast<int>(ceil(traceLength/sPeriod))
                    || iEndAir < 0 || iEndAir >= static_cast<int>(ceil(traceLength/sPeriod))){
                G4cerr << "Index out of bounds during transition radiation calculation"
                        << " (air contribution): " << G4endl;
                G4cerr << "iStartAir = " << iStartAir << G4endl;
                G4cerr << "iEndAir = " << iEndAir << G4endl;
                G4cerr << "Size of antenna trace: "
                            << static_cast<int>(ceil(traceLength/sPeriod)) << G4endl;
                G4cerr << "airStartE (weighted): " << CGSToSIFactor*airStartE << " V/m" << G4endl;
                G4cerr << "airEndE (weighted): " << CGSToSIFactor*airEndE << " V/m" << G4endl;
                exit(EXIT_FAILURE);
            }

            // Writing the contribution to the trace (only for start point!)
            fRunAction->AddToAntennaTrace(k, airStartE, iStartAir, 2);

        }

        // Next we get the contributions for the step boundary -> ice, of which we will only
        // keep the end point emission. Both direct and indirect rays.

        gotContributionIceDir = fSteppingAction->GetEfieldEP(k, antennaPosition,
                                         iceStartPoint, iceEndPoint, iceStartTime, iceEndTime,
                                         beta, fPrimaryCharge, 0,
                                         iceStartEDir, iceEndEDir, iStartIceDir, iEndIceDir);

        if(gotContributionIceDir){

            // Taking the weight of the primary into account
            iceStartEDir = fPrimaryCharge*iceStartEDir;
            iceEndEDir = fPrimaryWeight*iceEndEDir;

            // Checking if the time window for registration is large enough
            if(iStartIceDir < 0 || iStartIceDir >= static_cast<int>(ceil(traceLength/sPeriod))
                || iEndIceDir < 0 || iEndIceDir >= static_cast<int>(ceil(traceLength/sPeriod))){
                G4cerr << "Index out of bounds during transition radiation calculation"
                    << " (direct ice contribution): " << G4endl;
                G4cerr << "iStartIceDir = " << iStartIceDir << G4endl;
                G4cerr << "iEndIceDir = " << iEndIceDir << G4endl;
                G4cerr << "Size of antenna trace: "
                        << static_cast<int>(ceil(traceLength/sPeriod)) << G4endl;
                G4cerr << "iceStartEDir (weighted): " << CGSToSIFactor*iceStartEDir << " V/m"
                        << G4endl;
                G4cerr << "iceEndEDir (weighted): " << CGSToSIFactor*iceEndEDir << " V/m" << G4endl;
                exit(EXIT_FAILURE);
            }

            // Writing the contribution to the trace (only for end point!)
            fRunAction->AddToAntennaTrace(k, iceEndEDir, iEndIceDir, 2);

        }

        gotContributionIceIndir = fSteppingAction->GetEfieldEP(k, antennaPosition,
                                         iceStartPoint, iceEndPoint, iceStartTime, iceEndTime,
                                         beta, fPrimaryCharge, 1,
                                         iceStartEIndir, iceEndEIndir, iStartIceIndir, iEndIceIndir);

        if(gotContributionIceIndir){

            // Taking the weight of the primary into account
            iceStartEIndir = fPrimaryCharge*iceStartEIndir;
            iceEndEIndir = fPrimaryCharge*iceEndEIndir;

            // Checking if the time window for registration is large enough
            if(iStartIceIndir < 0 || iStartIceIndir >= static_cast<int>(ceil(traceLength/sPeriod))
                || iEndIceIndir < 0 || iEndIceIndir >= static_cast<int>(ceil(traceLength/sPeriod))){
                G4cerr << "Index out of bounds during transition radiation calculation"
                        << " (indirect ice contribution): " << G4endl;
                G4cerr << "iStartIceIndir = " << iStartIceIndir << G4endl;
                G4cerr << "iEndIceIndir = " << iEndIceIndir << G4endl;
                G4cerr << "Size of antenna trace: "
                        << static_cast<int>(ceil(traceLength/sPeriod)) << G4endl;
                G4cerr << "iceStartEIndir (weighted): " << CGSToSIFactor*iceStartEIndir << " V/m"
                            << G4endl;
                G4cerr << "iceEndEIndir (weighted): " << CGSToSIFactor*iceEndEIndir << " V/m"
                            << G4endl;
                exit(EXIT_FAILURE);
            }

            // Writing the contribution to the trace (only for end point!)
            fRunAction->AddToAntennaTrace(k, iceEndEIndir, iEndIceIndir, 2);

        }

    }

}

void IceShelfEventAction::CalcTransitionRadiationEndpoint(){

    if(fPrimaryCharge == 0){
        return;
    }

    // The magnitude of the momentum (p*c) of the primary
    G4double primaryMomMag = fPrimaryStartMomentum.mag();

    // Getting the arrival direction and velocity (beta) of the primary particle
    G4ThreeVector arrivalDir = fPrimaryStartMomentum/primaryMomMag;
    G4double betaMag = primaryMomMag/sqrt(fPrimaryMass*fPrimaryMass + primaryMomMag*primaryMomMag);
    G4ThreeVector beta = betaMag*arrivalDir;

    // Getting the start point for the air->boundary step, with corresponding time
    G4ThreeVector airStartPoint = fPrimaryStartPos - transRadStepsize*arrivalDir;
    G4double airStartTime = fStartingTime - transRadStepsize/(betaMag*c0);

    // Getting the end point for the boundary->ice step, with corresponding time
    G4ThreeVector iceEndPoint = fPrimaryStartPos + transRadStepsize*arrivalDir;
    G4double iceEndTime = fStartingTime + transRadStepsize/(betaMag*c0);

    //G4cout << "This is the event spreaking: " << G4endl;
    //G4cout << "The momentum (p*c) of the primary is " << primaryMomMag/GeV << " GeV" << G4endl;
    //G4cout << "The mass (m*c^2) of the primary is " << fPrimaryMass/MeV << " MeV" << G4endl;
    //G4cout << "arrivalDir of the primary is " << arrivalDir << G4endl;
    //G4cout << "The beta value of the primary is " << beta << G4endl;
    //G4cout << "The airStartPoint is at " << airStartPoint/m << " m" << G4endl;
    //G4cout << "The airStartTime is " << airStartTime/ns << " ns" << G4endl;
    //G4cout << "The iceEndPoint is at " << iceEndPoint/m << " m" << G4endl;
    //G4cout << "The iceEndTime is " << iceEndTime/ns << " ns" << G4endl;

    for(G4int k = 0; k < fRunAction->GetIceShelfAntennas()->GetNumberOfAntennas(); k++){

        // The antenna position
        G4ThreeVector antennaPosition = G4ThreeVector(
                                          fRunAction->GetIceShelfAntennas()->GetAntennaCoord(k, 0),
                                          fRunAction->GetIceShelfAntennas()->GetAntennaCoord(k, 1),
                                          fRunAction->GetIceShelfAntennas()->GetAntennaCoord(k, 2));

        // Doing the ray tracing for the air (direct) and ice (direct + indirect) points

        // Air
        G4double airObsTime;
        G4ThreeVector airLaunchDir;
        G4ThreeVector airReceiveDir;
        G4double airDoppler;
        G4double airR;
        G4double airIncidenceAngle;

        // Ice, direct
        G4double iceObsTimeD;
        G4ThreeVector iceLaunchDirD;
        G4ThreeVector iceReceiveDirD;
        G4double iceDopplerD;
        G4double iceRD;
        G4double iceIncidenceAngleD;

        // Ice, indirect
        G4double iceObsTimeI;
        G4ThreeVector iceLaunchDirI;
        G4ThreeVector iceReceiveDirI;
        G4double iceDopplerI;
        G4double iceRI;
        G4double iceIncidenceAngleI;

        G4ThreeVector airObsVect = antennaPosition - airStartPoint;
        G4ThreeVector iceObsVect = antennaPosition - iceEndPoint;

        // Ray tracing
        G4int gotContributionAir = fSteppingAction->GetEPValsRayTracing(k, airStartPoint, airObsVect,
                                                                               airStartTime, beta, 0,
                                                             airObsTime, airLaunchDir, airReceiveDir,
                                                                 airDoppler, airR, airIncidenceAngle);

        G4int gotContributionIceD = fSteppingAction->GetEPValsRayTracing(k, iceEndPoint, iceObsVect,
                                                                                iceEndTime, beta, 0,
                                                         iceObsTimeD, iceLaunchDirD, iceReceiveDirD,
                                                             iceDopplerD, iceRD, iceIncidenceAngleD);

        G4int gotContributionIceI = fSteppingAction->GetEPValsRayTracing(k, iceEndPoint, iceObsVect,
                                                                                iceEndTime, beta, 1,
                                                         iceObsTimeI, iceLaunchDirI, iceReceiveDirI,
                                                             iceDopplerI, iceRI, iceIncidenceAngleI);

        // Skipping this antenna if one of the doppler terms is too small
        if(gotContributionAir && airDoppler < approxThresholdTrans){
            continue;
        }
        if(gotContributionIceD && iceDopplerD < approxThresholdTrans){
            continue;
        }
        if(gotContributionIceI && iceDopplerI < approxThresholdTrans){
            continue;
        }

        // Applying the end-point formula

        // First switching to CGS unit system

        G4double sPeriodCGS = sPeriod/s;
        G4double chargeCGS = fPrimaryCharge/franklin;
        G4double c0CGS = c0/(cm/s);

        G4double airRCGS = airR/cm;
        G4double iceRDCGS = iceRD/cm;
        G4double iceRICGS = iceRI/cm;

        // Getting the time indices for the antenna trace array

        G4int iAir = static_cast<int>(floor(airObsTime/sPeriod + 0.5l));
        G4int iIceD = static_cast<int>(floor(iceObsTimeD/sPeriod + 0.5l));
        G4int iIceI = static_cast<int>(floor(iceObsTimeI/sPeriod + 0.5l));

        // Getting the electric fields

        G4ThreeVector airE = G4ThreeVector(0., 0., 0.);
        if(gotContributionAir){
            airE = (1./sPeriodCGS) * ((chargeCGS)/(c0CGS)) *
             ((airLaunchDir.cross(airLaunchDir.cross(beta))) /
                                                        (fabs(airDoppler)*airRCGS));
        }

        G4ThreeVector iceED = G4ThreeVector(0., 0., 0.);
        if(gotContributionIceD){
            iceED = (-1.0) * (1./sPeriodCGS) * ((chargeCGS)/(c0CGS)) *
             ((iceLaunchDirD.cross(iceLaunchDirD.cross(beta))) /
                                                        (fabs(iceDopplerD)*iceRDCGS));
        }

        G4ThreeVector iceEI = G4ThreeVector(0., 0., 0.);
        if(gotContributionIceI){
            iceEI = (-1.0) * (1./sPeriodCGS) * ((chargeCGS)/(c0CGS)) *
             ((iceLaunchDirI.cross(iceLaunchDirI.cross(beta))) /
                                                        (fabs(iceDopplerI)*iceRICGS));
        }

        // Doing the rotation...

        G4ThreeVector yDir = G4ThreeVector(0., 1., 0.); // Unit vector in y-dir (internal units)
        G4ThreeVector rotAxis;
        G4double cosVal;

        // ... for air
        rotAxis = yDir.cross(airLaunchDir);
        rotAxis = rotAxis/rotAxis.mag();
        cosVal = airLaunchDir.dot(airReceiveDir);
        airE = fSteppingAction->RotateVectorCosPi(airE, rotAxis, cosVal);

        // ... for ice (direct)
        rotAxis = yDir.cross(iceLaunchDirD);
        rotAxis = rotAxis/rotAxis.mag();
        cosVal = iceLaunchDirD.dot(iceReceiveDirD);
        iceED = fSteppingAction->RotateVectorCosPi(iceED, rotAxis, cosVal);

        // ... for ice (indirect)
        rotAxis = yDir.cross(iceLaunchDirI);
        rotAxis = rotAxis/rotAxis.mag();
        cosVal = iceLaunchDirI.dot(iceReceiveDirI);
        iceEI = fSteppingAction->RotateVectorCosPi(iceEI, rotAxis, cosVal);

        // Taking the Fresnel coefficients into account...

        // ... for air (transition)
        fSteppingAction->ApplyFresnelCoeff(airIncidenceAngle, airReceiveDir, airE, 1);

        // ... for ice (indirect), if reflection on the surface
        if(!std::isnan(iceIncidenceAngleI) && iceIncidenceAngleI != -1000.*degree){
            fSteppingAction->ApplyFresnelCoeff(iceIncidenceAngleI, iceReceiveDirI, iceEI, 0);
        }

        // Taking the weight of the primary into account
        airE = fPrimaryWeight*airE;
        iceED = fPrimaryWeight*iceED;
        iceEI = fPrimaryWeight*iceEI;

        // Writing the electric field contributions to the antenna array...

        // ... for air
        if(gotContributionAir){
            if(iAir < 0 || iAir >= static_cast<int>(ceil(traceLength/sPeriod))){
                G4cerr << "Index out of bounds during transition radiation calculation (air): "
                       << G4endl;
                G4cerr << "iAir = " << iAir << G4endl;
                G4cerr << "Size of antenna trace: " << static_cast<int>(ceil(traceLength/sPeriod))
                       << G4endl;
                G4cerr << "airE (weighted): " << CGSToSIFactor*airE << " V/m" << G4endl;
                exit(EXIT_FAILURE);
            }
            fRunAction->AddToAntennaTrace(k, airE, iAir, 2);
        }

        // ... for ice (direct)
        if(gotContributionIceD){
            if(iIceD < 0 || iIceD >= static_cast<int>(ceil(traceLength/sPeriod))){
                G4cerr << "Index out of bounds during transition radiation calculation (ice, direct): "
                       << G4endl;
                G4cerr << "iIceD = " << iIceD << G4endl;
                G4cerr << "Size of antenna trace: " << static_cast<int>(ceil(traceLength/sPeriod))
                       << G4endl;
                G4cerr << "iceED (weighted): " << CGSToSIFactor*iceED << " V/m" << G4endl;
                exit(EXIT_FAILURE);
            }
            fRunAction->AddToAntennaTrace(k, iceED, iIceD, 2);
        }

        // ... for ice (indirect)
        if(gotContributionIceI){
            if(iIceI < 0 || iIceI >= static_cast<int>(ceil(traceLength/sPeriod))){
                G4cerr << "Index out of bounds during transition radiation calculation (ice, indirect): "
                       << G4endl;
                G4cerr << "iIceI = " << iIceI << G4endl;
                G4cerr << "Size of antenna trace: " << static_cast<int>(ceil(traceLength/sPeriod))
                       << G4endl;
                G4cerr << "iceEI (weighted): " << CGSToSIFactor*iceEI << " V/m" << G4endl;
                exit(EXIT_FAILURE);
            }
            fRunAction->AddToAntennaTrace(k, iceEI, iIceI, 2);
        }
        
    }
}
*/

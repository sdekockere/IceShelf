#include "IceShelfRunAction.hh"
#include "IceShelfAnalysis.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include <fstream>

#include "IceShelfGlobalVariables.hh"

#include "IceShelfAntennas.hh"

#include "TFile.h"

#include "IceRayTracing.hh"
#include "AirToIceRayTracing.h"

IceShelfRunAction::IceShelfRunAction(G4String outputFileName, G4double* snapshotTimes, 
         G4double zenithAngleRad, G4double azimuthAngleRad, G4String reasFileName, G4String listFileName,
                                                                                G4String atmosFileName)
 : G4UserRunAction(),
   fOutputFileName(outputFileName),
   startOfRunTime(0),
   fSnapshotTimes(snapshotTimes)
{
    // The antenna coordinates
    if(reasFileName != "" && listFileName != ""){
        fIceShelfAntennas = new IceShelfAntennas(reasFileName, listFileName, azimuthAngleRad);
    } else {
        fIceShelfAntennas = new IceShelfAntennas(zenithAngleRad);
    }
    // For debugging
    G4cout << "### PRINTING THE ANTENNA LOCATIONS ###" << G4endl;
    for(G4int i = 0; i < fIceShelfAntennas->GetNumberOfAntennas(); i++){
        G4cout << "Antenna " << i << ": (" << fIceShelfAntennas->GetAntennaCoord(i, 0)/m << " m, "
            << fIceShelfAntennas->GetAntennaCoord(i, 1)/m << " m, "
            << fIceShelfAntennas->GetAntennaCoord(i, 2)/m << " m)" << G4endl;
    }

    if(rayTracingOn){

        // Setting the values for A, B and C
        //G4double newA = 1.78;
        //G4double newB = -0.454;
        //G4double newC = 0.0202;
        //IceRayTracing::SetA(newA);
        //IceRayTracing::SetB(newB);
        //IceRayTracing::SetC(newC);

        G4cout << G4endl;
        G4cout << "#############################" << G4endl;
        G4cout << "# RAYTRACING ON" << G4endl;
        if(includeIndirectRays){
            G4cout << "# INDIRECT RAY INCLUDED" << G4endl;
        } else {
            G4cout << "# INDIRECT RAY NOT INCLUDED" << G4endl;
        }
        if(includeAirToIceRT){
            G4cout << "# AIR TO ICE RT ENABLED" << G4endl;
        } else {
            G4cout << "# AIR TO ICE RT DISABLED" << G4endl;
        }
        if(IceRayTracing::TransitionBoundary == 0){
            G4cout << "# USING EXP. PROFILE FOR N: N(Z) = A + B*exp(-C*|Z|)" << G4endl;
            G4cout << "# WITH A = " << IceRayTracing::A_ice << ", B = " << IceRayTracing::GetB(0)
                       << ", C = " << IceRayTracing::GetC(0) << G4endl;
        } else {
            G4cout << "# USING GREENLAND DOUBLE EXP. PROFILE FOR N(Z)" << G4endl;
        }
        G4cout << "#############################" << G4endl;
        G4cout << G4endl;

        // Creating the interpolation tables for raytracing, one for each antenna
        IceRayTracing::SetNumberOfAntennas(fIceShelfAntennas->GetNumberOfAntennas());
        for(G4int i = 0; i < fIceShelfAntennas->GetNumberOfAntennas(); i++){

            G4double xPos = fIceShelfAntennas->GetAntennaCoord(i, 0);
            G4double yPos = fIceShelfAntennas->GetAntennaCoord(i, 1);
            G4double zPos = fIceShelfAntennas->GetAntennaCoord(i, 2);

            G4double showerHitDistance = sqrt(xPos*xPos + zPos*zPos);
            G4double antennaDepth = yPos - shelfSizeY/2.;

            IceRayTracing::MakeTable(showerHitDistance/m, centerShowerDepth/m, antennaDepth/m, i);

            // For debugging
            /*
            G4cout << "Log for antenna " << i << " out of " << fIceShelfAntennas->GetNumberOfAntennas()
                 << G4endl;
            G4cout << "Antenna position (x, y z): (" << xPos/m << " m, " <<
                yPos/m << " m, " << zPos/m << "m)" << G4endl;
            G4cout << "Made table with showerHitDistance = " << showerHitDistance/m
                << " m and antennaDepth = " << antennaDepth/m << " m" << G4endl;
            //// Use only when using the non-silent IceRayTracing.cc and IceRayTracing.hh
            //IceRayTracing::printInterpolTable(i);
            */

        }

        // Creating the atmosphere for indirect radiation calculations or in case
        // air to ice RT is enabled

        if(includeIndirectRays || includeAirToIceRT){
            G4cout << "ATMOSPHERE INFORMATION: " << G4endl;
            AirToIceRayTracing::MakeAtmosphere(atmosFileName);
            G4cout << "* Made the atmosphere using " << atmosFileName << G4endl;

            AirToIceRayTracing::A_const =
                              AirToIceRayTracing::Getnz_air(fIceShelfAntennas->GetIceShelfAltitude()/m);
            G4cout << "* Set A_const using an ice layer height of "
                    << fIceShelfAntennas->GetIceShelfAltitude()/m << " m" << G4endl;

            AirToIceRayTracing::UseConstantRefractiveIndex = true;
            AirToIceRayTracing::A_air = AirToIceRayTracing::A_const;

            G4cout << G4endl;
        }

    }

    // set printing run number for each run and event number for each event
    G4RunManager::GetRunManager()->SetPrintProgress(1);

    // Create vector that will contain additional histograms that can't be made with the 
    // limited G4AnalysisManager
    extra2DHistos = new std::vector<TH2*>;

    // Creat analysis manager
    // The choice of analysis technology is done via IceShelfAnalysis.hh
    auto analysisManager = G4AnalysisManager::Instance();
    G4cout << "Using " << analysisManager->GetType() << G4endl;

    // Creating the slice histograms
    // The histograms represent a slice of the simulated volume made out of blocks of which the 
    // dimension can be found in IceShelfGlobalVariables.hh. The slice follows the x-axis. 
    // The slice of blocks is centered around the x-axis, meaning each block goes from z = -0.5*blocksize
    // to z = 0.5*blocksize (remember: in Geant4 the z axis lies in the horizontal plane together with
    // the x-axis). The slice covers the entire simuated volume in x- and y-direction.
    //
    analysisManager->CreateH2("EnergyDensityProfile_gammas", "Energy density profile gammas", shelfSizeX/blockSize, (-1*shelfSizeX/2.)/cm, (shelfSizeX/2.)/cm, shelfSizeY/blockSize, (-1.*shelfSizeY)/cm, (0.*m)/cm);
    analysisManager->SetH2Title(0, "Deposited energy density in 1-cm-wide slice along y-axis (#gamma primaries only)");
    analysisManager->SetH2XAxisTitle(0, "x (cm)");
    analysisManager->SetH2YAxisTitle(0, "z (cm)");
    analysisManager->SetH2ZAxisTitle(0, "Deposited energy density (MeV/cm^{3})");

    analysisManager->CreateH2("EnergyDensityProfile_electrons", "Energy density profile e+ and e-", shelfSizeX/blockSize, (-1*shelfSizeX/2.)/cm, (shelfSizeX/2.)/cm, shelfSizeY/blockSize, (-1.*shelfSizeY)/cm, (0.*m)/cm);
    analysisManager->SetH2Title(1, "Deposited energy density in 1-cm-wide slice along y-axis (e^{#plus} and e^{#minus} primaries only)");
    analysisManager->SetH2XAxisTitle(1, "x (cm)");
    analysisManager->SetH2YAxisTitle(1, "z (cm)");
    analysisManager->SetH2ZAxisTitle(1, "Deposited energy density (MeV/cm^{3})");

    analysisManager->CreateH2("EnergyDensityProfile_muons", "Energy density profile mu+ and mu-", shelfSizeX/blockSize, (-1*shelfSizeX/2.)/cm, (shelfSizeX/2.)/cm, shelfSizeY/blockSize, (-1.*shelfSizeY)/cm, (0.*m)/cm);
    analysisManager->SetH2Title(2, "Deposited energy density in 1-cm-wide slice along y-axis (#mu^{#plus} and #mu^{#minus} primaries only)");
    analysisManager->SetH2XAxisTitle(2, "x (cm)");
    analysisManager->SetH2YAxisTitle(2, "z (cm)");
    analysisManager->SetH2ZAxisTitle(2, "Deposited energy density (MeV/cm^{3})");

    analysisManager->CreateH2("EnergyDensityProfile_hadrons", "Energy density profile hadrons", shelfSizeX/blockSize, (-1*shelfSizeX/2.)/cm, (shelfSizeX/2.)/cm, shelfSizeY/blockSize, (-1.*shelfSizeY)/cm, (0.*m)/cm);
    analysisManager->SetH2Title(3, "Deposited energy density in 1-cm-wide slice along y-axis (hadron primaries only)");
    analysisManager->SetH2XAxisTitle(3, "x (cm)");
    analysisManager->SetH2YAxisTitle(3, "z (cm)");
    analysisManager->SetH2ZAxisTitle(3, "Deposited energy density (MeV/cm^{3})");

    // Creating the radial histograms
    // These are 2D histograms, with the first dimension being radius and the second one depth
    // They will contain the energy left behind in the radius-depth bins, NOT the energy density
    // Note you have to take into account the volume represented by the bins when interpreting
    // these histograms.
    // The max value of the radius will be the smallest of shelfSizeX and shelfSizeZ
    // The binning is logarithmic to achieve better accuracy at smaller radii
    G4double histoMaxRadius = shelfSizeX/2.;
    if(shelfSizeZ < shelfSizeX){
        histoMaxRadius = shelfSizeZ/2.;
    }

    // The start of the first bin
    G4double x0 = 0.;
    x0 += shiftingValRadialHist;

    // The end of the first bin
    G4double xn = histoMaxRadius;
    xn += shiftingValRadialHist;

    // The bin width of the first bin
    G4double dX0 = x0*(pow(xn/x0, 1./nBinsRadialHist) - 1.);

    // The vector containing all the bin widths. Carefull, these are in internal units!
    G4double binWidths[nBinsRadialHist];
    binWidths[0] = dX0;
    for(G4int i = 1; i < nBinsRadialHist; i++){
        G4double binWidth = (dX0 + x0)*binWidths[i-1]/x0;
        binWidths[i] = binWidth;
    }
    
    // The vector containing the limits of the bins
    G4double binLimits[nBinsRadialHist + 1];
    binLimits[0] = x0;
    for(G4int i = 1; i < nBinsRadialHist + 1; i++){
        binLimits[i] = (binLimits[i-1] + binWidths[i-1]);
    }
    for(G4int i = 0; i < nBinsRadialHist + 1; i++){
        binLimits[i] -= shiftingValRadialHist;
    }
    for(G4int i = 0; i < nBinsRadialHist + 1; i++){
        binLimits[i] = binLimits[i]/cm;
    }

    // Creating the histogram
    extra2DHistos->push_back(new TH2D("RadialEnergyProfile_gammas", "Radial energy profile (#gamma primaries only)", nBinsRadialHist, 0., histoMaxRadius/cm, shelfSizeY/blockSize, (-1.*shelfSizeY)/cm, (0.*m)/cm));
    extra2DHistos->push_back(new TH2D("RadialEnergyProfile_electrons", "Radial energy profile (e^{#plus} and e^{#minus} primaries only)", nBinsRadialHist, 0., histoMaxRadius/cm, shelfSizeY/blockSize, (-1.*shelfSizeY)/cm, (0.*m)/cm));
    extra2DHistos->push_back(new TH2D("RadialEnergyProfile_muons", "Radial energy profile (#mu^{#plus} and #mu^{#minus} primaries only)", nBinsRadialHist, 0., histoMaxRadius/cm, shelfSizeY/blockSize, (-1.*shelfSizeY)/cm, (0.*m)/cm));
    extra2DHistos->push_back(new TH2D("RadialEnergyProfile_hadrons", "Radial energy profile (hadron primaries only)", nBinsRadialHist, 0., histoMaxRadius/cm, shelfSizeY/blockSize, (-1.*shelfSizeY)/cm, (0.*m)/cm));
    extra2DHistos->at(0)->GetXaxis()->SetTitle("Radius (cm)");
    extra2DHistos->at(1)->GetXaxis()->SetTitle("Radius (cm)");
    extra2DHistos->at(2)->GetXaxis()->SetTitle("Radius (cm)");
    extra2DHistos->at(3)->GetXaxis()->SetTitle("Radius (cm)");
    extra2DHistos->at(0)->GetYaxis()->SetTitle("Depth (cm)");
    extra2DHistos->at(1)->GetYaxis()->SetTitle("Depth (cm)");
    extra2DHistos->at(2)->GetYaxis()->SetTitle("Depth (cm)");
    extra2DHistos->at(3)->GetYaxis()->SetTitle("Depth (cm)");
    extra2DHistos->at(0)->GetZaxis()->SetTitle("Energy in bin (MeV)");
    extra2DHistos->at(1)->GetZaxis()->SetTitle("Energy in bin (MeV)");
    extra2DHistos->at(2)->GetZaxis()->SetTitle("Energy in bin (MeV)");
    extra2DHistos->at(3)->GetZaxis()->SetTitle("Energy in bin (MeV)");
    
    // Adjusting the bin limits of the histograms, changing it to a logarithmic binning
    extra2DHistos->at(0)->GetXaxis()->Set(nBinsRadialHist, binLimits);
    extra2DHistos->at(1)->GetXaxis()->Set(nBinsRadialHist, binLimits);
    extra2DHistos->at(2)->GetXaxis()->Set(nBinsRadialHist, binLimits);
    extra2DHistos->at(3)->GetXaxis()->Set(nBinsRadialHist, binLimits);

    // Creating the number of particles in function of depth histograms
    // Each bin represents the number of particles that reached a depth of lower bin limit
    // Depth here is in the shower axis frame, so it is the distance along the shower axis
    // The y-axis has energy distributions
    analysisManager->CreateH2("NumberOfParticlesDepth_gammas", "Number of #gamma in function of depth", depthNumberOfBins, 0., maxDepthValue/m, depthEnergyNumberOfBins, depthEnergyUnderLim/MeV, depthEnergyUpperLim/MeV);
    analysisManager->SetH2Title(4, "Number of #gamma in function of depth");
    analysisManager->SetH2XAxisTitle(4, "Depth (m)");
    analysisManager->SetH2YAxisTitle(4, "Kinetic energy (MeV)");
    analysisManager->SetH2ZAxisTitle(4, "Number of #gamma");

    analysisManager->CreateH2("NumberOfParticlesDepth_positrons", "Number of e^{#plus} in function of depth", depthNumberOfBins, 0., maxDepthValue/m, depthEnergyNumberOfBins, depthEnergyUnderLim/MeV, depthEnergyUpperLim/MeV);
    analysisManager->SetH2Title(5, "Number of e^{#plus} in function of depth");
    analysisManager->SetH2XAxisTitle(5, "Depth (m)");
    analysisManager->SetH2YAxisTitle(5, "Kinetic energy (MeV)");
    analysisManager->SetH2ZAxisTitle(5, "Number of e^{#plus}");

    analysisManager->CreateH2("NumberOfParticlesDepth_electrons", "Number of e^{#minus} in function of depth", depthNumberOfBins, 0., maxDepthValue/m, depthEnergyNumberOfBins, depthEnergyUnderLim/MeV, depthEnergyUpperLim/MeV);
    analysisManager->SetH2Title(6, "Number of e^{#minus} in function of depth");
    analysisManager->SetH2XAxisTitle(6, "Depth (m)");
    analysisManager->SetH2YAxisTitle(6, "Kinetic energy (MeV)");
    analysisManager->SetH2ZAxisTitle(6, "Number of e^{#minus}");

    analysisManager->CreateH2("NumberOfParticlesDepth_antimuons", "Number of #mu^{#plus} in function of depth", depthNumberOfBins, 0., maxDepthValue/m, depthEnergyNumberOfBins, depthEnergyUnderLim/MeV, depthEnergyUpperLim/MeV);
    analysisManager->SetH2Title(7, "Number of #mu^{#plus} in function of depth");
    analysisManager->SetH2XAxisTitle(7, "Depth (m)");
    analysisManager->SetH2YAxisTitle(7, "Kinetic energy (MeV)");
    analysisManager->SetH2ZAxisTitle(7, "Number of #mu^{#plus}");

    analysisManager->CreateH2("NumberOfParticlesDepth_muons", "Number of #mu^{#minus} in function of depth", depthNumberOfBins, 0., maxDepthValue/m, depthEnergyNumberOfBins, depthEnergyUnderLim/MeV, depthEnergyUpperLim/MeV);
    analysisManager->SetH2Title(8, "Number of #mu^{#minus} in function of depth");
    analysisManager->SetH2XAxisTitle(8, "Depth (m)");
    analysisManager->SetH2YAxisTitle(8, "Kinetic energy (MeV)");
    analysisManager->SetH2ZAxisTitle(8, "Number of #mu^{#minus}");

    analysisManager->CreateH2("NumberOfParticlesDepth_hadrons", "Number of hadrons in function of depth", depthNumberOfBins, 0., maxDepthValue/m, depthEnergyNumberOfBins, depthEnergyUnderLim/MeV, depthEnergyUpperLim/MeV);
    analysisManager->SetH2Title(9, "Number of hadrons in function of depth");
    analysisManager->SetH2XAxisTitle(9, "Depth (m)");
    analysisManager->SetH2YAxisTitle(9, "Kinetic energy (MeV)");
    analysisManager->SetH2ZAxisTitle(9, "Number of hadrons");

    // Creating the charge in function of time histograms
    // Each bin represents the amount of charges at time t = lower bin limit
    // The y-axis has energy distributions
    analysisManager->CreateH2("ChargeTimeProfile_gammas", "Total charge in function of time (#gamma primaries only)", timeNumberOfBins, 0., maxTimeValue/ns, chargeEnergyNumberOfBins, chargeEnergyUnderLim/MeV, chargeEnergyUpperLim/MeV);
    analysisManager->SetH2Title(10, "Total charge in function of time (#gamma primaries only)");
    analysisManager->SetH2XAxisTitle(10, "Time (ns)");
    analysisManager->SetH2YAxisTitle(10, "Kinetic energy (MeV)");
    analysisManager->SetH2ZAxisTitle(10, "Charge (e^{#plus})");

    analysisManager->CreateH2("ChargeTimeProfile_electrons", "Total charge in function of time (e^{#plus} and e^{#minus} primaries only)", timeNumberOfBins, 0., maxTimeValue/ns, chargeEnergyNumberOfBins, chargeEnergyUnderLim/MeV, chargeEnergyUpperLim/MeV);
    analysisManager->SetH2Title(11, "Total charge in function of time (e^{#plus} and e^{#minus} primaries only)");
    analysisManager->SetH2XAxisTitle(11, "Time (ns)");
    analysisManager->SetH2YAxisTitle(11, "Kinetic energy (MeV)");
    analysisManager->SetH2ZAxisTitle(11, "Charge (e^{#plus})");

    analysisManager->CreateH2("ChargeTimeProfile_muons", "Total charge in function of time (#mu^{#plus} and #mu^{#minus} primaries only)", timeNumberOfBins, 0., maxTimeValue/ns, chargeEnergyNumberOfBins, chargeEnergyUnderLim/MeV, chargeEnergyUpperLim/MeV);
    analysisManager->SetH2Title(12, "Total charge in function of time (#mu^{#plus} and #mu^{#minus} primaries only)");
    analysisManager->SetH2XAxisTitle(12, "Time (ns)");
    analysisManager->SetH2YAxisTitle(12, "Kinetic energy (MeV)");
    analysisManager->SetH2ZAxisTitle(12, "Charge (e^{#plus})");

    analysisManager->CreateH2("ChargeTimeProfile_hadrons", "Total charge in function of time (hadron primaries only)", timeNumberOfBins, 0., maxTimeValue/ns, chargeEnergyNumberOfBins, chargeEnergyUnderLim/MeV, chargeEnergyUpperLim/MeV);
    analysisManager->SetH2Title(13, "Total charge in function of time (hadron primaries only)");
    analysisManager->SetH2XAxisTitle(13, "Time (ns)");
    analysisManager->SetH2YAxisTitle(13, "Kinetic energy (MeV)");
    analysisManager->SetH2ZAxisTitle(13, "Charge (e^{#plus})");

    // Creating the snapshot ntuples
    // x, y and z will be in the shower axis coord system as defined by the 
    // GetShowerAxisCoordinates function in IceShelfSteppingAction

    analysisManager->SetNtupleMerging(true);

    for(G4int i = 0; i < numberOfSnapshotDepthsMassArea; i++){
        G4String tupleName = "SnapshotNtuple_" + std::to_string(i);
        G4String tupleTitle 
            = "x (m), y (m), z (m), charge (e+) & weight for all charges at cascade front depth "
            + std::to_string(round(snapshotDepthsMassArea[i]/(g/cm2))) + " g/cm2 (t = "
            + std::to_string(fSnapshotTimes[i]/ns) + " ns)";
        analysisManager->CreateNtuple(tupleName, tupleTitle);
        analysisManager->CreateNtupleDColumn("x");
        analysisManager->CreateNtupleDColumn("y");
        analysisManager->CreateNtupleDColumn("z");
        analysisManager->CreateNtupleDColumn("charge");
        analysisManager->CreateNtupleDColumn("weight");
        analysisManager->FinishNtuple();
    }

    // Creating the vector array representing the trace of the antennas for the direct emission

    if(includeDirectRays){
        fAntennaTracesDir = new std::vector<std::vector<G4ThreeVector>>
                                                (fIceShelfAntennas->GetNumberOfAntennas());
        for(G4int k = 0; k < fIceShelfAntennas->GetNumberOfAntennas(); k++){
            for(int j = 0; j < static_cast<int>(ceil(traceLength/sPeriod)); j++){
                fAntennaTracesDir->at(k).push_back(G4ThreeVector(0.,0.,0.));
            }
        }
    }

    // Creating the vector array representing the trace of the antennas for the direct emission

    if(includeIndirectRays){
        fAntennaTracesIndir = new std::vector<std::vector<G4ThreeVector>>
                                                (fIceShelfAntennas->GetNumberOfAntennas());
        for(G4int k = 0; k < fIceShelfAntennas->GetNumberOfAntennas(); k++){
            for(int j = 0; j < static_cast<int>(ceil(traceLength/sPeriod)); j++){
                fAntennaTracesIndir->at(k).push_back(G4ThreeVector(0.,0.,0.));
            }
        }
    }

    // Creating the vector array representing the trace of the antennas for
    // the sudden appearance emission

    if(writeSARadSeparateToo){
        fAntennaTracesSA = new std::vector<std::vector<G4ThreeVector>>
                                                (fIceShelfAntennas->GetNumberOfAntennas());
        for(G4int k = 0; k < fIceShelfAntennas->GetNumberOfAntennas(); k++){
            for(int j = 0; j < static_cast<int>(ceil(traceLength/sPeriod)); j++){
                fAntennaTracesSA->at(k).push_back(G4ThreeVector(0.,0.,0.));
            }
        }
    }

    /*
    // For testing
    AddToAntennaTrace(0, G4ThreeVector(1., 2., 3.) , 2);
    AddToAntennaTrace(1, G4ThreeVector(3., 2., 1.) , 2);
    AddToAntennaTrace(1, G4ThreeVector(1., 1., 1.) , 2);
    G4cout << "###########################" << G4endl;
    for(unsigned int k = 0; k < fAntennaTracesDir->size(); k++){
        G4cout << "This is antenna " << k << G4endl;
        G4cout << "Its position: (" << fIceShelfAntennas->GetAntennaCoord(k, 0)/m << "m, " 
            << fIceShelfAntennas->GetAntennaCoord(k, 1)/m << "m, "
            << fIceShelfAntennas->GetAntennaCoord(k, 2)/m << "m)" << G4endl;
        G4cout << "Its trace size is " << fAntennaTracesDir->at(k).size() << G4endl;
        G4cout << "Its first trace element: " << fAntennaTracesDir->at(k).at(0) << G4endl;
        G4cout << "Its last trace element: " << fAntennaTracesDir->at(k).at(fAntennaTracesDir->at(k).size()-1) << G4endl;
        G4cout << "Its third trace element: " << fAntennaTracesDir->at(k).at(2) << G4endl;
    }
    */
}

IceShelfRunAction::~IceShelfRunAction()
{
    if(includeDirectRays){
        delete fAntennaTracesDir;
    }
    if(includeIndirectRays){
        delete fAntennaTracesIndir;
    }
    if(writeSARadSeparateToo){
        delete fAntennaTracesSA;
    }
    delete G4AnalysisManager::Instance();
    for(unsigned int i = 0; i < extra2DHistos->size(); i++){
        delete extra2DHistos->at(i);
    }
    delete extra2DHistos;
    delete fIceShelfAntennas;
}

void IceShelfRunAction::BeginOfRunAction(const G4Run*)
{
    // Get analysis manager
    auto analysisManager = G4AnalysisManager::Instance();
    
    // Open an output file for the histogram and Ntuple
    analysisManager->OpenFile(fOutputFileName + ".root");

    // Set startOfRunTime
    startOfRunTime = std::clock();
}

void IceShelfRunAction::EndOfRunAction(const G4Run*)
{

    // Save the root file
    G4cout << "Saving the root file..." << G4endl;
    auto analysisManager = G4AnalysisManager::Instance();
    analysisManager->Write();
    analysisManager->CloseFile();

    // Save the antenna traces: writing it to a txt file
    G4cout << "Saving the antenna traces..." << G4endl;
    if(includeDirectRays){
        SaveAntennaTraces(fAntennaTracesDir, fOutputFileName + "_dir.txt");
    }
    if(includeIndirectRays){
        SaveAntennaTraces(fAntennaTracesIndir, fOutputFileName + "_indir.txt");
    }
    if(writeSARadSeparateToo){
        SaveAntennaTraces(fAntennaTracesSA, fOutputFileName + "_SA.txt");
    }

    // Save the extra 2D histos
    G4cout << "Saving the extra root file..." << G4endl;
    G4String fnameExtra2DHistos = fOutputFileName + "_extra_2D_histos.root";
    TFile* fExtra2DHistos = new TFile(fnameExtra2DHistos, "RECREATE");
    fExtra2DHistos->cd();
    for(unsigned int i = 0; i < extra2DHistos->size(); i++){
        extra2DHistos->at(i)->Write();
    }
    fExtra2DHistos->Close();
    delete fExtra2DHistos;

    // Wrapping up
    clock_t endOfRunTime = std::clock();
    G4double runTime = double(endOfRunTime - startOfRunTime)/double(CLOCKS_PER_SEC);
    G4cout << "Time taken for run to complete: " << runTime << " s" << G4endl;
}

// traceType should be 0 if adding to direct emission traces,
// traceType should be 1 if adding to indirect emission traces,
// traceType should be 2 if adding to SA radiation traces
void IceShelfRunAction::AddToAntennaTrace(G4int antennaNumber, G4ThreeVector elecField,
                                          G4int index, G4int traceType){
    if(traceType == 0){
        if(includeDirectRays){
            fAntennaTracesDir->at(antennaNumber).at(index)
                          = fAntennaTracesDir->at(antennaNumber).at(index) + elecField;
        } else {
            G4cerr << "Trying to add to direct emission traces, but includeDirectRays is false."
                   << G4endl;
            exit(EXIT_FAILURE);
        }
    }
    else if(traceType == 1){
        if(includeIndirectRays){
            fAntennaTracesIndir->at(antennaNumber).at(index)
                          = fAntennaTracesIndir->at(antennaNumber).at(index) + elecField;
        } else {
            G4cerr << "Trying to add to indirect emission traces, but includeIndirectRays is false."
                   << G4endl;
            exit(EXIT_FAILURE);
        }
    }
    else if(traceType == 2){
        if(writeSARadSeparateToo){
            fAntennaTracesSA->at(antennaNumber).at(index)
                          = fAntennaTracesSA->at(antennaNumber).at(index) + elecField;
        } else {
            G4cerr << "Trying to add to sudden appearance radiation traces, but writeSARadSeparateToo" 
                   << " is false." << G4endl;
            exit(EXIT_FAILURE);
        }
    }
    else{
        G4cerr << "ERROR in AddToAntennaTrace from the IceShelfRunAction class: traceType is "
               << traceType << ", but should be 0, 1 or 2." << G4endl;
        exit(EXIT_FAILURE);
    }
}

void IceShelfRunAction::SaveAntennaTraces(std::vector<std::vector<G4ThreeVector>>* antennaTraces, 
                                                                    G4String outputFileName){

    std::ofstream traceFile;
    traceFile.precision(15);
    traceFile.open(outputFileName);
    traceFile << "# The energy threshold above which the kinetic energy of a charged particle should be to take part in the electrical field calculations is " << kinEnergyThreshold/MeV << " MeV" << G4endl;
    traceFile << "# time (ns)\telectric field (V/m)\n";
    for(G4int k = 0; k < fIceShelfAntennas->GetNumberOfAntennas(); k++){
        traceFile << "# Antenna " << k << " at (" << fIceShelfAntennas->GetAntennaCoord(k, 0)/m << "m, " 
            << fIceShelfAntennas->GetAntennaCoord(k, 1)/m << "m, "
            << fIceShelfAntennas->GetAntennaCoord(k, 2)/m << "m) (Geant4 coordinates)\n";
        bool foundFirstNonZero = false;
        int i = 0;
        while(!foundFirstNonZero && i < static_cast<int>(ceil(traceLength/sPeriod))){
            if(antennaTraces->at(k).at(i).x() != 0. || antennaTraces->at(k).at(i).y() != 0. 
                                                        || antennaTraces->at(k).at(i).z() != 0.){
                foundFirstNonZero = true;
            } else {
                i++;
            }
        }
        i = i - static_cast<int>(ceil(startBuffer/sPeriod));
        if(i < 0){ // Not enough buffer before the first non-zero value
            i = 0;
        }
        G4int iEnd = i + static_cast<int>(ceil(traceLengthToSave/sPeriod));
        while(i <= iEnd && i < static_cast<int>(ceil(traceLength/sPeriod))){
            G4double tPoint = i*sPeriod;
            traceFile << tPoint/ns << " " << CGSToSIFactor*antennaTraces->at(k).at(i) << "\n";
            i++;
        }
    }
    traceFile.close();

}

G4bool IceShelfRunAction::FillH2(G4int id, G4double xvalue, G4double yvalue, G4double weight){

    extra2DHistos->at(id)->Fill(xvalue, yvalue, weight); 
    return true;

}

IceShelfAntennas* IceShelfRunAction::GetIceShelfAntennas(){
    return fIceShelfAntennas;
}

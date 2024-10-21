/*///////////////////////////////////////////////////////////////////////////////////////////////
This is the main code of the ice_shelf project.
Launching the executable without any arguments will start up the interactive GUI, showing the ice shelf
and allowing the user to shoot gammas through it while visualizing the result.
If the user does not want to use the interactive GUI, a run macro file (like run.mac) should be provided
as argument which will then control the simulation without activating the interactive GUI.

The program expects at least one argument, being the full path to the input file containing a list of
primary particles. Use the python script make_geant4_input_file.py for creating the right input files.
*////////////////////////////////////////////////////////////////////////////////////////////////

#include "IceShelfDetectorConstruction.hh"
#include "IceShelfActionInitialization.hh"

#include "G4RunManager.hh"
#include "G4UImanager.hh"

//#include "FTFP_BERT.hh"
#include "SimplePhysicsList.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "Randomize.hh"

#include "G4SystemOfUnits.hh"

#include <vector>
#include <stdexcept>
#include <sstream>
#include <fstream>
#include <iostream>

#include "IceShelfGlobalVariables.hh"
#include "IceDensityModels.hh"

#include "IceRayTracing.hh"

int main(int argc, char** argv)
{

    G4String reasFileName = "";
    G4String listFileName = "";
    G4String atmosphereFileName = "";

    if(includeIndirectRays || includeAirToIceRT){
        if(argc != 9){
            G4cerr << "To run this program use " << argv[0] << " <primary particles input file> <name for output file> <random number used for engine seed (index of seedtable will be this number modulo 215, so using 0 or 215 e.g. will result in same engine seed)> <zenith angle of primary particle in degrees> <azimuth angle of primary particle in degrees> <reas file from the CoREAS simulation> <list file from the CoREAS simulation> <atmosphere dat file used in the corresponding CoREAS sim>" << G4endl;
            exit(EXIT_FAILURE);
        }
        reasFileName = argv[6];
        listFileName = argv[7];
        atmosphereFileName = argv[8];
    } else {
        if(argc != 6 && argc != 8){
            G4cerr << "To run this program use " << argv[0] << " <primary particles input file> <name for output file> <random number used for engine seed (index of seedtable will be this number modulo 215, so using 0 or 215 e.g. will result in same engine seed)> <zenith angle of primary particle in degrees> <azimuth angle of primary particle in degrees> \n+ OPTIONAL: <reas file from the CoREAS simulation> <list file from the CoREAS simulation>" << G4endl;
            exit(EXIT_FAILURE);
        }
        if(argc == 8){
            reasFileName = argv[6];
            listFileName = argv[7];
        }
    }

    clock_t start, end;
    start = std::clock();
    G4cout << "Starting the timing for total time calculation here" << G4endl;

    //// Getting the angles in internal units (radians) ////

    G4double zenithAngleRad = std::stof(argv[4])*degree;
    G4double azimuthAngleRad = std::stof(argv[5])*degree;

    //// Choose the Random engine ////
    
    G4Random::setTheEngine(new CLHEP::RanecuEngine);
    std::stringstream charToInt;
    charToInt << argv[3];
    G4int randNumber;
    charToInt >> randNumber;
    G4Random::setTheSeed(randNumber);
    G4cout << "Random number received as argument: " << randNumber << G4endl;
    G4cout << "Zenith angle received as argument (transformed into radians): " << zenithAngleRad 
           << G4endl;
    G4cout << "Azimuth angle received as argument (transformed into radians): " << azimuthAngleRad
           << G4endl;
    G4cout << "Index of the SeedTable: " << G4Random::getTheSeed() << G4endl;
    G4cout << "Dimension of the simulated ice shelf: (" << shelfSizeX/m << "m, " 
           << shelfSizeY/m << "m, " << shelfSizeZ/m << "m)" << G4endl;
    G4cout << "Number of ice layers used for the density gradient: " << numberOfLayers << G4endl; 
    G4cout << "Breaking index of the ice: " << nIce << G4endl;
    G4cout << "sPeriod: " << sPeriod/ns << " ns" << G4endl;
    G4cout << "traceLength: " << traceLength/ns << " ns" << G4endl;
    G4cout << "startBuffer: " << startBuffer/ns << " ns" << G4endl;
    G4cout << "traceLengthToSave: " << traceLengthToSave/ns << " ns" << G4endl;
    G4cout << "kinEnergyThreshold: " << kinEnergyThreshold/MeV << " MeV" << G4endl;
    G4cout << "blockSize: " << blockSize/cm << " cm" << G4endl;
    G4cout << "nBinsRadialHist: " << nBinsRadialHist << G4endl;
    G4cout << "shiftingValRadialHist: " << shiftingValRadialHist/cm << " cm" << G4endl;
    G4cout << "maxDepthValue: " << maxDepthValue/m << " m" << G4endl;
    G4cout << "depthNumberOfBins: " << depthNumberOfBins << G4endl;
    G4cout << "depthKinELimVals: [" << depthKinELimVals[0]/GeV << " GeV, "
           << depthKinELimVals[1]/GeV << " GeV, " << depthKinELimVals[2]/GeV << " GeV, "
           << depthKinELimVals[3]/GeV << " GeV]" << G4endl;
    G4cout << "depthEnergyNumberOfBins: " << depthEnergyNumberOfBins << G4endl;
    G4cout << "depthEnergyUnderLim: " << depthEnergyUnderLim/GeV << " GeV" << G4endl;
    G4cout << "depthEnergyUpperLim: " << depthEnergyUpperLim/GeV << " GeV" << G4endl;
    G4cout << "maxTimeValue: " << maxTimeValue/ns << " ns" << G4endl;
    G4cout << "timeNumberOfBis: " << timeNumberOfBins << G4endl;
    G4cout << "chargeKinELim: " << chargeKinELim/MeV << " MeV" << G4endl;
    G4cout << "chargeEnergyNumberOfBins: " << chargeEnergyNumberOfBins << G4endl;
    G4cout << "chargeEnergyUnderLim: " << chargeEnergyUnderLim/GeV << " GeV" << G4endl;
    G4cout << "chargeEnergyUpperLim: " << chargeEnergyUpperLim/GeV << " GeV" << G4endl;
    G4cout << "numberOfSnapshotDepthsMassArea: " << numberOfSnapshotDepthsMassArea << G4endl;
    G4cout << "snapshotDepthsMassArea: [";
    for(int i = 0; i < numberOfSnapshotDepthsMassArea; i++){
        if(i == numberOfSnapshotDepthsMassArea - 1){
            G4cout << snapshotDepthsMassArea[i]/(g/cm2) << " g/cm2";
        } else {
            G4cout << snapshotDepthsMassArea[i]/(g/cm2) << " g/cm2, ";
        }
    }
    G4cout << "]" << G4endl;
    G4cout << "snapshotKinELim: " << snapshotKinELim/GeV << " GeV" << G4endl;
    G4cout << "approxThreshold: " << approxThreshold << G4endl;
    G4cout << "zeroDepthOffset: " << zeroDepthOffset/m << " m" << G4endl;

    //// Construct the default run manager ////

    G4RunManager* runManager = new G4RunManager;

    //// Read the input file and store the particles in a list
    G4String inputFile = argv[1];
    std::vector<std::vector<G4double>>* primaryList = new std::vector<std::vector<G4double>>;
    std::ifstream file(inputFile);
    //// The units used in the input file
    G4String energyUnit;
    G4String lengthUnit;
    G4String timeUnit;
    //// Reading line per line, searching for the units used (should be mentioned in the comments block in the beginning)
    G4String line = ""; 
    while(std::getline(file, line)){
        G4String lineCheck = line.substr(0, 7);
        if(lineCheck == "# Units"){
            std::istringstream iss(line);
            for(G4int i = 0; i < 3; i++){
                std::getline(iss, energyUnit, ' '); // Skipping the first 3 words
            }
            std::getline(iss, energyUnit, ' '); 
            std::getline(iss, lengthUnit, ' '); 
            std::getline(iss, timeUnit, ' '); 
            std::cout << "Energy unit: " << energyUnit << std::endl;
            std::cout << "Length unit: " << lengthUnit << std::endl;
            std::cout << "Time unit: " << timeUnit << std::endl;
            break;
        }
    }
    //// Checking if we found the units used
    if(energyUnit.size() == 0 or lengthUnit.size() == 0 or timeUnit.size() == 0){
        throw std::invalid_argument("Could not determine all the units used in the input file.");
    }
    if(energyUnit == lengthUnit or energyUnit == timeUnit or lengthUnit == timeUnit){
        throw std::invalid_argument("Something wrong with the units. Please check.");
    }   
    //// Going back to the beginning of the file
    file.clear();
    file.seekg(0, std::ios::beg);
    //// Checking the number of particles in the file
    G4int numberOfParticles = 0;
    while(std::getline(file, line)){
        if(line.at(0) != '#'){
            numberOfParticles++;
        }
    }
    //// Ending if number of particles is 0
    if(numberOfParticles == 0){
        G4cout << "No particles in input file. Aborting..." << G4endl;
        return 0;
    }
    //// Going back to the beginning of the file
    file.clear();
    file.seekg(0, std::ios::beg);
    //// Reading the file, now storing the primary particles info
    while(std::getline(file, line)){
        if(line.at(0) != '#'){
            std::vector<G4double> entryList;
            //// Splitting the line on spaces into entries
            std::stringstream ss(line);
            G4String entry = "";
            G4int counter = 0;
            while(std::getline(ss, entry, ' ')){
                G4double entryD = (G4double)std::stod(entry);
                if(counter == 1 or counter == 2 or counter == 3){
                    //// Reading momenta
                    if(energyUnit == "GeV"){
                        entryD = entryD*GeV;
                    } else if (energyUnit == "MeV"){
                        entryD = entryD*MeV;
                    } else if (energyUnit == "TeV"){
                        entryD = entryD*TeV;
                    } else {
                        throw std::invalid_argument("Used energy unit is not yet supported!");
                    }
                }
                if(counter == 4 or counter == 5 or counter == 6){
                    //// Reading positions
                    if(lengthUnit == "cm"){
                        entryD = entryD*cm;
                    } else if(lengthUnit == "m"){
                        entryD = entryD*m;
                    } else {
                        throw std::invalid_argument("Used length unit is not yet supported!");
                    }
                }
                if(counter == 7){
                    //// Reading time
                    if(timeUnit == "ns"){
                        entryD = entryD*ns;
                    } else {
                        throw std::invalid_argument("Used time unit is not yet supported!");
                    }
                }
                entryList.push_back(entryD);
                counter++;
            }
            primaryList->push_back(entryList);
        }
    }

    //// Converting the snapshotDepthsMassArea values to time values ////
    //// We start at depth = 0 and go down the shower axis in tiny steps to find the length ////
    //// along the shower axis corresponding to the given depth (in mass/area) values ////
    //// Using velocity shower = c0, we then convert this to a time value ////

    G4double snapshotTimes[numberOfSnapshotDepthsMassArea];

    G4double integrationStepSize = (shelfSizeY/numberOfLayers)/1000.;
    G4double depthMassAreaVal = 0.;
    G4double pathAlongShowerAxis = 0.;
    G4int depthMassAreaIndex = 0;

    while(depthMassAreaIndex < numberOfSnapshotDepthsMassArea){
         
        if(depthMassAreaVal > snapshotDepthsMassArea[depthMassAreaIndex]){
            snapshotTimes[depthMassAreaIndex] = pathAlongShowerAxis/c0;
            depthMassAreaIndex++;
        }

        // Incrementing the path along shower axis
        pathAlongShowerAxis += integrationStepSize;

        // Incrementing the shower depth (units mass/area)
        G4double meanVerticalDepth 
                            = (pathAlongShowerAxis-0.5*integrationStepSize)*cos(zenithAngleRad);
        G4double meanDensity = IceDensityModels::rho(meanVerticalDepth);
        depthMassAreaVal += meanDensity*integrationStepSize;
        
    }

    /*
    G4cout << "The time values corresponding to the given depth values are:" << G4endl;
    for(G4int i = 0; i < numberOfSnapshotDepthsMassArea; i++){
        G4cout << snapshotTimes[i]/ns << " ns, ";
    }
    G4cout << G4endl;
    */

    //// Set mandatory initialization classes ////

    // Detector construction
    auto detConstruction = new IceShelfDetectorConstruction();
    runManager->SetUserInitialization(detConstruction);

    // Physics list (contains the physics processes)
    //auto physicsList = new FTFP_BERT;
    auto physicsList = new SimplePhysicsList;
    physicsList->SetVerboseLevel(1);
    runManager->SetUserInitialization(physicsList);

    // User action initialization (includes primary particle generation)
    G4String outputFileName = argv[2];
    auto actionInitialization = new IceShelfActionInitialization(outputFileName, 
             detConstruction, primaryList, snapshotTimes, zenithAngleRad, azimuthAngleRad,
                                                      reasFileName, listFileName, atmosphereFileName);
    runManager->SetUserInitialization(actionInitialization);

    //// Initialize visualization ////
    
    G4VisManager* visManager = new G4VisExecutive;
    // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance
    // G4VisManager* visManager = new G4VisExecutive("Quiet");
    visManager->Initialize();

    //// Get the pointer to the User Interface mamanger ////

    // The UI manager is used to call set fuctions of certain class objects 
    // of which we do not know the pointers
    // Other commands which allow for detailed control of the simulation are also available
    
    G4UImanager* UImanager = G4UImanager::GetUIpointer();

    //// Process macro or start UI session

    // batch mode
    G4String command = "/run/initialize";
    UImanager->ApplyCommand(command);
    command = "/control/verbose 0";
    UImanager->ApplyCommand(command);
    command = "/run/verbose 0";
    UImanager->ApplyCommand(command);
    command = "/event/verbose 0";
    UImanager->ApplyCommand(command);
    command = "/tracking/verbose 0";
    UImanager->ApplyCommand(command);
    runManager->BeamOn(numberOfParticles);

    //// job termination ////

    // user actions, physics_list and detector_description are owned and deleted by the run manager,
    // so they should not be deleted in the main program

    delete visManager;
    delete runManager;
    delete primaryList;

    end = std::clock();
    G4double time_taken = double(end - start)/double(CLOCKS_PER_SEC);
    G4cout << "Total time taken: " << time_taken << " s" << G4endl;

    return 0;
}

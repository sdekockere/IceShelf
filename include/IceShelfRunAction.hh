/*/////////////////////////////////////////////////////////////////////////////////////////////////////
This is the RunAction class. It takes care of things that should happen right before the start of a run
or at the end of a run. This is done using the functions BeginOfRunAction() and EndOfRunAction().
In this project the constructor of the RunAction class makes several histograms. Each histogram represents a slice of the simulated volume made out of blocks of which the dimension can be found in IceShelfGlobalVariables.hh. The slice follows the x-axis. The slice of blocks is centered around the x-axis, meaning each block goes from z = -0.5*blocksize to z = 0.5*blocksize (remember: in Geant4 the z axis lies in the horizontal plane together with the x-axis). The slice covers the entire simuated volume in x- and y-direction. The histograms will contain the energy density in each block, one for gamma primaries, one for electrons and positron primaries and one for muon and antimuon primaries.
The histograms are filled by the IceShelfSteppingAction object.
The RunAction class creates an output file at the beginning of a run using the BeginOfRunAction()
function and saves the histograms at the end of the run using the EndOfRunAction() function.
The private field fOutputFileName contains the name of the outputfile.
*//////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef ICESHELFRUNACTION_h
#define ICESHELFRUNACTION_h

#include "G4UserRunAction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"

#include <fstream>

#include "TH2.h"

class G4Run;
class IceShelfAntennas;

class IceShelfRunAction : public G4UserRunAction
{
    public:
        IceShelfRunAction(G4String, G4double*, G4double, G4double, G4String, G4String, G4String);
        virtual ~IceShelfRunAction();

        virtual void BeginOfRunAction(const G4Run*);
        virtual void EndOfRunAction(const G4Run*);

        // Adds elecField to the G4ThreeVector at index for antenna antennaNumber
        void AddToAntennaTrace(G4int antennaNumber, G4ThreeVector elecField, G4int index,
                                                                                    G4int traceType);

        // Fills the 2D histogram at index id in extra2DHistos
        G4bool FillH2(G4int id, G4double xvalue, G4double yvalue, G4double weight=1.0);

        IceShelfAntennas* GetIceShelfAntennas();

    private:
        G4String fOutputFileName;
        // The antenna traces for the direct emission
        std::vector<std::vector<G4ThreeVector>>* fAntennaTracesDir;
        // The antenna traces for the indirect emission
        std::vector<std::vector<G4ThreeVector>>* fAntennaTracesIndir;
        // The antenna traces for the sudden appearance emission
        std::vector<std::vector<G4ThreeVector>>* fAntennaTracesSA;
        // 2D histos that can't be made with limited G4AnalysisManager
        std::vector<TH2*>* extra2DHistos;
        clock_t startOfRunTime;
        G4double* fSnapshotTimes;
        IceShelfAntennas* fIceShelfAntennas;
        // Saves the antenna traces to a txt file
        void SaveAntennaTraces(std::vector<std::vector<G4ThreeVector>>*, G4String);
};

#endif

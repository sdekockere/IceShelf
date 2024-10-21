/*/////////////////////////////////////////////////////////////////////////////////////////////////////
* This is the IceShelfAntennas class, which will represent the antennas in the simulaion.
*//////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef ICESHELFANTENNAS_H
#define ICESHELFANTENNAS_H

#include "IceShelfGlobalVariables.hh"
#include <vector>
#include <G4String.hh>

class IceShelfAntennas
{

    public:
        IceShelfAntennas(G4double zenithAngleRad);
        IceShelfAntennas(G4String coreasReasFile, G4String coreasListFile, G4double azimuthAngleRad);
        virtual ~IceShelfAntennas();

        G4int GetNumberOfAntennas();
        // Returns coordinate k (0 for x, 1 for y and 2 for z) of antenna antennaNumber
        G4double GetAntennaCoord(G4int antennaNumber, G4int i);
        G4double GetIceShelfAltitude();

    private:

        std::vector<std::vector<G4double>> fAntennas; // The positions of the antennas (G4 axis system)
        G4double fIceShelfAltitude; // The height wrt sea level of the start of the ice shelf

        // Functions that can be used for the constructor without reas and list file
        void starGrid();
        void armGrid(G4double xStart, G4double xEnd, G4double dx);
        void cherenkovAngleGridFixedRad(G4double zenithAngleRad, G4double fixedRad);
        void cherenkovAngleGridFixedDepth(G4double zenithAngleRad, G4double fixedDepth);

};

inline G4int IceShelfAntennas::GetNumberOfAntennas(){
    return fAntennas.size();
}

inline G4double IceShelfAntennas::GetAntennaCoord(G4int antennaNumber, G4int i){
    return fAntennas.at(antennaNumber).at(i);
}

inline G4double IceShelfAntennas::GetIceShelfAltitude(){
    return fIceShelfAltitude;
}

#endif

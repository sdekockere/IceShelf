#include "IceShelfAntennas.hh"
#include "IceRayTracing.hh"
#include <math.h>
#include <iostream>
#include <sstream>
#include <fstream>

#include "G4UImanager.hh"
#include <G4String.hh>

// This constructor defines the antenna setup using the function specified in its body
IceShelfAntennas::IceShelfAntennas(G4double zenithAngleRad)
{

    armGrid(-200*m, 200*m, 2*m); 
    //starGrid();
    //cherenkovAngleGridFixedRad(zenithAngleRad, 100.*m);
    //cherenkovAngleGridFixedDepth(zenithAngleRad, 100.*m);

    // As we're not using the CoREAS input files, we will not be connecting the simulation to
    // a CoREAS simulation, and therefor the height of the ice shelf does not matter
    fIceShelfAltitude = -9999*m;

}

// This constructor defines the antenna setup using the CoREAS input files.
// Arguments are the .reas file and the .list file, as specified by the CoREAS manual.
IceShelfAntennas::IceShelfAntennas(G4String coreasReasFile, G4String coreasListFile, 
                                                                            G4double azimuthAngleRad)
{

    std::vector<std::vector<G4double>> antennas = {};

    // Reading the reas file, getting the spatial observer configuration
    G4String coreCoordinateNorth = "";
    G4String coreCoordinateWest = "";
    G4String coreCoordinateVertical = "";
    G4String iceBoundaryAltitude = "";
    std::ifstream reasFile(coreasReasFile);
    G4String line = "";
    while(std::getline(reasFile, line)){

        std::istringstream iss(line); // Needed to go through the line

        // Looking for the right lines
        G4String lineCheck = "";
        std::getline(iss, lineCheck, ' '); // Getting the first word
        if(lineCheck == "CoreCoordinateNorth"){
            for(G4int i = 0; i < 2; i++){
                std::getline(iss, coreCoordinateNorth, ' '); // Getting the third word
            }
            std::istringstream cleaner(coreCoordinateNorth); // To get rid of the extra tabs in the end
            std::getline(cleaner, coreCoordinateNorth, '\t');
        } else if(lineCheck == "CoreCoordinateWest"){
            for(G4int i = 0; i < 2; i++){
                std::getline(iss, coreCoordinateWest, ' '); // Getting the third word
            }
            std::istringstream cleaner(coreCoordinateWest); // To get rid of the extra tabs in the end
            std::getline(cleaner, coreCoordinateWest, '\t');
        } else if(lineCheck == "CoreCoordinateVertical"){
            for(G4int i = 0; i < 2; i++){
                std::getline(iss, coreCoordinateVertical, ' '); // Getting the third word
            }
            std::istringstream cleaner(coreCoordinateVertical); // To get rid of extra tabs in the end
            std::getline(cleaner, coreCoordinateVertical, '\t');
        } else if(lineCheck == "IceBoundaryAltitude"){
            for(G4int i = 0; i < 2; i++){
                std::getline(iss, iceBoundaryAltitude, ' '); // Getting the third word
            }
            std::istringstream cleaner(iceBoundaryAltitude); // To get rid of extra tabs in the end
            std::getline(cleaner, iceBoundaryAltitude, '\t');
        }

    }

    // Checking whether the spatial observer was placed at the air-ice boundary
    if(coreCoordinateVertical != iceBoundaryAltitude){
        G4cerr << "According to the .reas file, the spatial observer was not placed at the air-ice"
        << " boundary altitude." << G4endl;
        exit(EXIT_FAILURE);
    }

    // Changing to Geant4 internal units
    G4double coreNorth = (G4double)std::stod(coreCoordinateNorth);
    coreNorth = coreNorth*cm;
    G4double coreWest = (G4double)std::stod(coreCoordinateWest);
    coreWest = coreWest*cm;
    G4double coreVertical = (G4double)std::stod(coreCoordinateVertical);
    coreVertical = coreVertical*cm;

    // Setting the correct ice shelf altitude
    fIceShelfAltitude = coreVertical;

    // Checking whether coreCoordNorth and coreCoordWest are 0.
    // If not, you will probably miss the block of ice to start with, so we want to avoid this.
    if(coreNorth != 0 || coreWest != 0){
        G4cerr << "The core of the footprint should be at CoreCoordinateNorth = 0 and CoreCoordinateWest = 0, which is not the case according to the reas file. CoreCoordinateNorth = " << coreNorth/cm << " cm and CoreCoordinateWest = " << coreWest/cm << " cm." << G4endl;
        exit(EXIT_FAILURE);
    }

    // Reading the list file, getting the antenna coordinates
    std::ifstream listFile(coreasListFile);
    line = "";
    while(std::getline(listFile, line)){

        // Reading an antenna line
        G4String northCoordStr = "";
        G4String westCoordStr = "";
        G4String heightASLStr = "";
        std::istringstream iss(line);
        for(G4int i = 0; i < 2; i++){
            std::getline(iss, northCoordStr, ' '); // Skipping the first two words
        }
        std::getline(iss, northCoordStr, ' ');
        std::getline(iss, westCoordStr, ' ');
        std::getline(iss, heightASLStr, ' ');

        // Changing to Geant4 internal units
        G4double northCoord = (G4double)std::stod(northCoordStr);
        northCoord = northCoord*cm;
        G4double heightASL = (G4double)std::stod(heightASLStr);
        heightASL = heightASL*cm;
        G4double westCoord = (G4double)std::stod(westCoordStr);
        westCoord = westCoord*cm;
        // and now changing to Geant4 axis system:
        // * x pointing to the north, but rotated ccw over the azimuth angle
        //      (so in this frame, azimuth angle = 0)
        // * y points upwards, starting in the middle of the block of ice
        // * z pointing to the east, but rotated ccw over the azimuth angle
        G4double xCoord = northCoord*cos(azimuthAngleRad) + westCoord*sin(azimuthAngleRad);
        G4double yCoord = heightASL - (coreVertical - shelfSizeY/2.);
        G4double zCoord = -1.*(westCoord*cos(azimuthAngleRad) - northCoord*sin(azimuthAngleRad));

        // Adding it to the antenna list
        antennas.push_back({xCoord, yCoord, zCoord});

    }

    fAntennas = antennas;

}

IceShelfAntennas::~IceShelfAntennas()
{}

// This setup corresponds to a radial symmetric start shaped grid (121 antennas,
// 15 in each arm + 1 in the center, spacing of 20 m)
void IceShelfAntennas::starGrid()
{

    fAntennas = 
        {
            {-300.0*m, -140.*m, 0.*m},
            {-280.0*m, -140.*m, 0.*m},
            {-260.0*m, -140.*m, 0.*m},
            {-240.0*m, -140.*m, 0.*m},
            {-220.0*m, -140.*m, 0.*m},
            {-200.0*m, -140.*m, 0.*m},
            {-180.0*m, -140.*m, 0.*m},
            {-160.0*m, -140.*m, 0.*m},
            {-140.0*m, -140.*m, 0.*m},
            {-120.0*m, -140.*m, 0.*m},
            {-100.0*m, -140.*m, 0.*m},
            {-80.0*m, -140.*m, 0.*m},
            {-60.0*m, -140.*m, 0.*m},
            {-40.0*m, -140.*m, 0.*m},
            {-20.0*m, -140.*m, 0.*m},
            {20.0*m, -140.*m, 0.*m},
            {40.0*m, -140.*m, 0.*m},
            {60.0*m, -140.*m, 0.*m},
            {80.0*m, -140.*m, 0.*m},
            {100.0*m, -140.*m, 0.*m},
            {120.0*m, -140.*m, 0.*m},
            {140.0*m, -140.*m, 0.*m},
            {160.0*m, -140.*m, 0.*m},
            {180.0*m, -140.*m, 0.*m},
            {200.0*m, -140.*m, 0.*m},
            {220.0*m, -140.*m, 0.*m},
            {240.0*m, -140.*m, 0.*m},
            {260.0*m, -140.*m, 0.*m},
            {280.0*m, -140.*m, 0.*m},
            {300.0*m, -140.*m, 0.*m},
            {0.*m, -140.*m, -300.0*m},
            {0.*m, -140.*m, -280.0*m},
            {0.*m, -140.*m, -260.0*m},
            {0.*m, -140.*m, -240.0*m},
            {0.*m, -140.*m, -220.0*m},
            {0.*m, -140.*m, -200.0*m},
            {0.*m, -140.*m, -180.0*m},
            {0.*m, -140.*m, -160.0*m},
            {0.*m, -140.*m, -140.0*m},
            {0.*m, -140.*m, -120.0*m},
            {0.*m, -140.*m, -100.0*m},
            {0.*m, -140.*m, -80.0*m},
            {0.*m, -140.*m, -60.0*m},
            {0.*m, -140.*m, -40.0*m},
            {0.*m, -140.*m, -20.0*m},
            {0.*m, -140.*m, 20.0*m},
            {0.*m, -140.*m, 40.0*m},
            {0.*m, -140.*m, 60.0*m},
            {0.*m, -140.*m, 80.0*m},
            {0.*m, -140.*m, 100.0*m},
            {0.*m, -140.*m, 120.0*m},
            {0.*m, -140.*m, 140.0*m},
            {0.*m, -140.*m, 160.0*m},
            {0.*m, -140.*m, 180.0*m},
            {0.*m, -140.*m, 200.0*m},
            {0.*m, -140.*m, 220.0*m},
            {0.*m, -140.*m, 240.0*m},
            {0.*m, -140.*m, 260.0*m},
            {0.*m, -140.*m, 280.0*m},
            {0.*m, -140.*m, 300.0*m},
            {-212.132034356*m, -140.*m, -212.132034356*m},
            {-197.989898732*m, -140.*m, -197.989898732*m},
            {-183.847763109*m, -140.*m, -183.847763109*m},
            {-169.705627485*m, -140.*m, -169.705627485*m},
            {-155.563491861*m, -140.*m, -155.563491861*m},
            {-141.421356237*m, -140.*m, -141.421356237*m},
            {-127.279220614*m, -140.*m, -127.279220614*m},
            {-113.13708499*m, -140.*m, -113.13708499*m},
            {-98.9949493661*m, -140.*m, -98.9949493661*m},
            {-84.8528137424*m, -140.*m, -84.8528137424*m},
            {-70.7106781187*m, -140.*m, -70.7106781187*m},
            {-56.5685424949*m, -140.*m, -56.5685424949*m},
            {-42.4264068712*m, -140.*m, -42.4264068712*m},
            {-28.2842712475*m, -140.*m, -28.2842712475*m},
            {-14.1421356237*m, -140.*m, -14.1421356237*m},
            {14.1421356237*m, -140.*m, 14.1421356237*m},
            {28.2842712475*m, -140.*m, 28.2842712475*m},
            {42.4264068712*m, -140.*m, 42.4264068712*m},
            {56.5685424949*m, -140.*m, 56.5685424949*m},
            {70.7106781187*m, -140.*m, 70.7106781187*m},
            {84.8528137424*m, -140.*m, 84.8528137424*m},
            {98.9949493661*m, -140.*m, 98.9949493661*m},
            {113.13708499*m, -140.*m, 113.13708499*m},
            {127.279220614*m, -140.*m, 127.279220614*m},
            {141.421356237*m, -140.*m, 141.421356237*m},
            {155.563491861*m, -140.*m, 155.563491861*m},
            {169.705627485*m, -140.*m, 169.705627485*m},
            {183.847763109*m, -140.*m, 183.847763109*m},
            {197.989898732*m, -140.*m, 197.989898732*m},
            {212.132034356*m, -140.*m, 212.132034356*m},
            {-212.132034356*m, -140.*m, 212.132034356*m},
            {-197.989898732*m, -140.*m, 197.989898732*m},
            {-183.847763109*m, -140.*m, 183.847763109*m},
            {-169.705627485*m, -140.*m, 169.705627485*m},
            {-155.563491861*m, -140.*m, 155.563491861*m},
            {-141.421356237*m, -140.*m, 141.421356237*m},
            {-127.279220614*m, -140.*m, 127.279220614*m},
            {-113.13708499*m, -140.*m, 113.13708499*m},
            {-98.9949493661*m, -140.*m, 98.9949493661*m},
            {-84.8528137424*m, -140.*m, 84.8528137424*m},
            {-70.7106781187*m, -140.*m, 70.7106781187*m},
            {-56.5685424949*m, -140.*m, 56.5685424949*m},
            {-42.4264068712*m, -140.*m, 42.4264068712*m},
            {-28.2842712475*m, -140.*m, 28.2842712475*m},
            {-14.1421356237*m, -140.*m, 14.1421356237*m},
            {14.1421356237*m, -140.*m, -14.1421356237*m},
            {28.2842712475*m, -140.*m, -28.2842712475*m},
            {42.4264068712*m, -140.*m, -42.4264068712*m},
            {56.5685424949*m, -140.*m, -56.5685424949*m},
            {70.7106781187*m, -140.*m, -70.7106781187*m},
            {84.8528137424*m, -140.*m, -84.8528137424*m},
            {98.9949493661*m, -140.*m, -98.9949493661*m},
            {113.13708499*m, -140.*m, -113.13708499*m},
            {127.279220614*m, -140.*m, -127.279220614*m},
            {141.421356237*m, -140.*m, -141.421356237*m},
            {155.563491861*m, -140.*m, -155.563491861*m},
            {169.705627485*m, -140.*m, -169.705627485*m},
            {183.847763109*m, -140.*m, -183.847763109*m},
            {197.989898732*m, -140.*m, -197.989898732*m},
            {212.132034356*m, -140.*m, -212.132034356*m},
            {0.*m, -140.*m, 0.*m}
        };

}

// This function puts antennas at a depth of 150 m (y = -140m in Geant4 axis system)  along the x-axis.
// It starts at xStart, ends at xEnd and using a spacing distance of dx.
void IceShelfAntennas::armGrid(G4double xStart, G4double xEnd, G4double dx)
{

    std::vector<std::vector<G4double>> antennas = {};

    G4double xPos = xStart;
    while(xPos <= xEnd){
        antennas.push_back({xPos, -140.*m, 0.});
        xPos += dx;
    }


    fAntennas = antennas;

}

// This function creates 81 antennas close to the expected Cherenkov angle in the ice,
// at a distance fixedRad from the point of impact on the surface.
// It will place them on the left side of the shower axis (negative x direction).
// The antennas will be placed on a grid, with angular separation of 0.1 degrees.
// Antenna 30 will correspond to the antenna on the Cherenkov angle (wrt the ice surface)
// Don't forget to sepcify the unit of fixedDepth!
void IceShelfAntennas::cherenkovAngleGridFixedRad(G4double zenithAngleRad, G4double fixedRad)
{

    G4double thetaCh;

    if(rayTracingOn){
        // Getting the value for n around emDep, which should be the depth in the ice
        // where you expect most of the radiation to come from
        G4double emDep = -5*m;
        G4double nDeep = IceRayTracing::Getnz(emDep/m);
        thetaCh = acos(1./nDeep);
    }
    else{
        thetaCh = acos(1./nIce);
    }

    // The cherenkov angle wrt the y-axis (negative value if left of y-axis)
    G4double thetaY = zenithAngleRad - thetaCh;

    // The z coordinate of the antenna at the Cherenkov angle
    G4double zCoord = 0.;

    std::vector<std::vector<G4double>> antennas = {};

    G4double distance = fixedRad;

    G4double deltaTheta = 0.1*pi/180.;
    for(G4int i = -30; i < 51; i++){
        G4double xCoord = distance*sin(thetaY + i*deltaTheta);
        G4double yCoord = -1.*distance*cos(thetaY + i*deltaTheta) + shelfSizeY/2.;
        antennas.push_back({xCoord, yCoord, zCoord});
    }

    fAntennas = antennas;

}

// This function creates 81 antennas close to the expected Cherenkov angle.
// It will place them on the left side of the shower axis (negative x direction).
// The antenna on the Cherenkov angle (wrt the surface) will be at fixedDepth, 
// the others will be spread out around this antenna in certain angular steps,
// keeping the total distance to the point where the shower hit the surface the same.
// Antenna 30 will correspond to the antenna on the Cherenkov angle (wrt the ice surface)
// Don't forget to specify the unit of fixedDepth!
void IceShelfAntennas::cherenkovAngleGridFixedDepth(G4double zenithAngleRad, G4double fixedDepth)
{

    G4double thetaCh;

    if(rayTracingOn){
        // Getting the value for n around emDep, which should be the depth in the ice
        // where you expect most of the radiation to come from
        G4double emDep = -5*m;
        G4double nDeep = IceRayTracing::Getnz(emDep/m);
        thetaCh = acos(1./nDeep);
    }
    else{
        thetaCh = acos(1./nIce);
    }

    // The cherenkov angle wrt the y-axis (negative value if left of axis)
    G4double thetaY = zenithAngleRad - thetaCh;

    // The z coordinate of the antenna at the Cherenkov angle
    G4double zCoord = 0.;

    std::vector<std::vector<G4double>> antennas = {};

    G4double distance = fixedDepth/cos(thetaY); // The fixed distance of the antennas to the
                                                // point of impact on the ice

    G4double deltaTheta = 0.1*pi/180.;
    for(G4int i = -30; i < 51; i++){
        G4double xCoord = distance*sin(thetaY + i*deltaTheta);
        G4double yCoord = -1.*distance*cos(thetaY + i*deltaTheta) + shelfSizeY/2.;
        antennas.push_back({xCoord, yCoord, zCoord});
    }

    fAntennas = antennas;

}

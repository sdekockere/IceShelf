#include "IceShelfDetectorConstruction.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"

#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"

#include "G4NistManager.hh"

#include "IceShelfGlobalVariables.hh"

#include <vector>

#include "IceDensityModels.hh"

IceShelfDetectorConstruction::IceShelfDetectorConstruction() : G4VUserDetectorConstruction()
{ 
    fIceLayers = new std::vector<G4VPhysicalVolume*>;
}

IceShelfDetectorConstruction::~IceShelfDetectorConstruction()
{ 
    delete fIceLayers;
}

G4VPhysicalVolume* IceShelfDetectorConstruction::Construct()
{
    //// Creating the material Air for the World volume ////
    G4NistManager* man = G4NistManager::Instance();
    G4Material* Air = man->FindOrBuildMaterial("G4_AIR");

    //// Creating the material Ice for the Ice Shelf volume ////

    G4double z, a; // The atomic number and mass per mole
    G4String name, symbol;
    G4int ncomponents, natoms;
    G4double temp, pres; // The temperature and pressure
    
    a = 1.00794*g/mole;
    G4Element* elH = new G4Element(name="Hydrogen", symbol="H", z=1., a);

    a = 15.9994*g/mole;
    G4Element* elO = new G4Element(name="Oxygen", symbol="O", z=8., a);

    // Introducing a density gradient. For this we divide the ice shelf into layers, each layer
    // having a different density. We follow the density profile 
    // rho(z) = rho_ice - (rho_ice - rho_surface)*exp(-1.9*z/t_firn)).
    // rho_ice is the asymptotic density of polar ice, 
    // rho_surface is the density of the snow at the surface for South Pole, 
    // t_firn is the thickness of the firn layer at South Pole in m.

    G4double dy = shelfSizeY/numberOfLayers; // The width of each layer
    std::vector<G4double> layerMiddlePoints; // This vector stores the middle points of the layer, relative to the middle of the world volume. Layer 0 is at the surface.
    std::vector<G4Material*> iceLayerMaterials;

    for(int i = 0; i < numberOfLayers; i++){
        G4String materialName = "IceLayerMaterial" + std::to_string(i);
        G4double mPoint = -9999.9; // The middle point of layer i
        mPoint = (i + 0.5)*dy;
        layerMiddlePoints.push_back(shelfSizeY/2. - mPoint);
        G4double density = -9999.9; // The density of layer i
        density = IceDensityModels::rho(mPoint);
        G4Material* iceLayerMaterial = new G4Material(name=materialName, density, ncomponents=2, kStateSolid, temp=258.15*kelvin, pres=1.25*atmosphere);
        iceLayerMaterial->AddElement(elH, natoms=2);
        iceLayerMaterial->AddElement(elO, natoms=1);
        iceLayerMaterials.push_back(iceLayerMaterial);    
        /*
        G4cout << "Material name: " << materialName << G4endl;
        G4cout << "Layer width: " << dy/m << " m" << G4endl;
        G4cout << "Layer middle point: " << mPoint/m << " m" << G4endl;
        G4cout << "Density: " << density/(kg/m3) << " kg/m3" << G4endl;
        */
    }
    
    //// Creating the solids //// 

    // Option to switch on/off checking of volumes overlaps
    G4bool checkOverlaps = false; // Put this on true if you change the layout to check for no overlaps 
                                  // once. Checking for overlaps for all these layers takes a while so
                                  // standard value is false.

    // World. Dimensions are based on the size of the ice shelf,
    // which is fixed in IceShelfGlobalVariables.hh
    G4double world_sizeX = 1.5*shelfSizeX;
    G4double world_sizeY = 1.5*shelfSizeY;
    G4double world_sizeZ = 1.5*shelfSizeZ;

    G4Box* solidWorld = new G4Box("World", 0.5*world_sizeX, 0.5*world_sizeY, 0.5*world_sizeZ);
    G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld, Air, "World");
    G4VPhysicalVolume* physWorld = new G4PVPlacement(0,                 // no rotation
                                                     G4ThreeVector(),   // at (0, 0, 0)
                                                     logicWorld,        // its logical volume
                                                     "World",           // its name
                                                     0,                 // its mother volume
                                                     false,             // no boolean operation
                                                     0,                 // its copy number
                                                     checkOverlaps);    // overlaps checking

    // Ice Shelf
    std::vector<G4Box*> solidIceLayers;
    std::vector<G4LogicalVolume*> logicIceLayers;
    for (int i = 0; i < numberOfLayers; i++){
        G4String volumeName = "IceLayer" + std::to_string(i);
        solidIceLayers.push_back(new G4Box(volumeName, 0.5*shelfSizeX, 0.5*dy, 0.5*shelfSizeZ));
        logicIceLayers.push_back(new G4LogicalVolume(solidIceLayers[i], iceLayerMaterials[i], volumeName));
        (*fIceLayers).push_back(new G4PVPlacement(
                 0,                                           // no rotation
                 G4ThreeVector(0, layerMiddlePoints[i] , 0),  // translation to this point
                 logicIceLayers[i],                           // its logical volume
                 volumeName,                                  // its name
                 logicWorld,                                  // its mother volume
                 false,                                       // no boolean operation
                 0,                                           // its copy number
                 checkOverlaps));                              // overlaps checking
    }

    /*
    for(int i = 0; i < numberOfLayers; i++){
        G4cout << "#####################################################" << G4endl;
        G4cout << (*fIceLayers)[i]->GetLogicalVolume()->GetMaterial() << G4endl;
        G4cout << "Position of this layer with respect to center of World (in m): " << (*fIceLayers)[i]->GetTranslation()/m << G4endl;
        G4cout << "#####################################################" << G4endl;
    }
    */

    //// Returning the physical world ////

    return physWorld;
}

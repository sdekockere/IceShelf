/*//////////////////////////////////////////////////////////////////////////////////////////////////
This is the DetectorConstruction class. It takes care of the construction of the detector. It is 
obviously a mandatory class.
The construction of the detector is implemented in the Construct() function. The Construct() function 
creates two volumes: the World volume and the Ice Shelf volume. The Ice Shelf volume is placed inside 
the surrounding World volume. Both are made of ice. The Ice Shelf volume is made up of differebt layers each having a different density.
The private field fIceLayers contains a pointer to the list containing the physical volume instances of the ice shelf layers. This isn't really usefull right now but might be in the future.
*///////////////////////////////////////////////////////////////////////////////////////////////////


#ifndef ICESHELFDETECTORCONSTRUCTION_h
#define ICESHELFDETECTORCONSTRUCTION_h

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

#include <vector>

class G4VPhysicalVolume;
class G4LogicalVolume;

class IceShelfDetectorConstruction : public G4VUserDetectorConstruction
{
    public:
        IceShelfDetectorConstruction();
        virtual ~IceShelfDetectorConstruction();

        virtual G4VPhysicalVolume* Construct();

        const std::vector<G4VPhysicalVolume*>* GetIceShelfPV() const;

    private:
        std::vector<G4VPhysicalVolume*>* fIceLayers; // The ice shelf volume
};

// inline functions

inline const std::vector<G4VPhysicalVolume*>* IceShelfDetectorConstruction::GetIceShelfPV() const{
    return fIceLayers;
}

#endif

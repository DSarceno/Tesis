#include "MagneticField.hh"

#include "G4GenericMessenger.hh"
#include "G4SystemOfUnits.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MagneticField::MagneticField()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MagneticField::~MagneticField()
{ 
}

void MagneticField::GetFieldValue(const G4double Point[4],double *bField) const
{


  if (Point[2]< -20*cm) {
  bField[0] = 0.;
  bField[1] = 0.;
  bField[2] = 0.;
  }

  if (Point[2]> -20*cm) {


  bField[0] = 0.;
  bField[1] = -1.*tesla;
  bField[2] = 0;
  }

  if (Point[2]> 20*cm) {


  bField[0] = 0.0;
  bField[1] = 0.0;
  bField[2] = 0.0;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: TrackingAction.cc 81501 2014-05-31 13:51:51Z ldesorgh $
//
/// \file eventgenerator/exgps/src/TrackingAction.cc
/// \brief Implementation of the TrackingAction class
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "TrackingAction.hh"
#include "EventAction.hh"
#include "RunAction.hh"
#include "Analysis.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4UnitsTable.hh"
#include "G4Track.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

TrackingAction::TrackingAction()
:G4UserTrackingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void TrackingAction::PreUserTrackingAction(const G4Track* track)
{

  G4double time           = track->GetGlobalTime();
  G4int pid               = track->GetDynamicParticle()->GetPDGcode();
  G4double ekin           = track->GetKineticEnergy();
  G4double etot           = track->GetTotalEnergy();
  G4ThreeVector vertex    = track->GetPosition();
  G4ThreeVector momentum  = track->GetMomentum();
  
  auto analysisManager = G4AnalysisManager::Instance();
  
  //Primary Particles
  if (time == 0){
  G4double x = vertex.x(), y = vertex.y(), z = vertex.z();
  
  // fill ntuple  
  analysisManager->FillNtupleDColumn(3,0,pid);
  analysisManager->FillNtupleDColumn(3,1,ekin/(1*MeV));
  analysisManager->FillNtupleDColumn(3,2,etot/(1*MeV));
  analysisManager->FillNtupleDColumn(3,3,x);
  analysisManager->FillNtupleDColumn(3,4,y);
  analysisManager->FillNtupleDColumn(3,5,z);
  analysisManager->FillNtupleDColumn(3,6,momentum.x());
  analysisManager->FillNtupleDColumn(3,7,momentum.y());
  analysisManager->FillNtupleDColumn(3,8,momentum.z());

  analysisManager->AddNtupleRow(3);

  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....







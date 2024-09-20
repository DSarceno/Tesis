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
// $Id: SteppingAction.cc 100946 2016-11-03 11:28:08Z gcosmo $
// 
/// \file SteppingAction.cc
/// \brief Implementation of the SteppingAction class

#include "SteppingAction.hh"
#include "EventAction.hh"
#include "RunAction.hh"
#include "Analysis.hh"

#include "DetectorConstruction.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4SDManager.hh"
#include "G4Track.hh"
#include "G4SystemOfUnits.hh"

#include "G4Step.hh"
#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(
                      const DetectorConstruction* detectorConstruction,
                      EventAction* eventAction)
  : G4UserSteppingAction(),
    fDetConstruction(detectorConstruction),
    fEventAction(eventAction)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* step)
{
// Collect energy and track length step by step

  // get volume of the current step
  auto volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();

  G4String IP_NameCheck = "IP_Sensitive";
  G4int copyNumber;
  G4int partID;
  
  // energy deposit
  auto edep = step->GetTotalEnergyDeposit();

  auto analysisManager = G4AnalysisManager::Instance();

  G4Track* track;

 
  if (volume == fDetConstruction->GetETT_FrontPV()) {
    track = step->GetTrack();

    auto pos = track->GetPosition();
    auto momentum = track->GetMomentum();

    //Debugging

    //G4cout << track->GetDefinition()->GetParticleName() << "\n";
    //G4cout << track->GetDynamicParticle()->GetPDGcode() << "\n";
    //G4cout << track->GetKineticEnergy() << "\n";
    //G4cout << track->GetTotalEnergy() << "\n";
    //G4cout << track->GetPosition() << "\n";
    //G4cout << track->GetMomentum() << "\n";
    //G4cout << track->GetMomentumDirection() << "\n";
    if (momentum.z()>0){
      partID = track->GetDynamicParticle()->GetPDGcode();
      if (partID == 22){
        auto px = momentum.x();
        auto py = momentum.y();
        auto pz = momentum.z();
        auto energy = sqrt(px*px + py*py + pz*pz);

        analysisManager->FillH1(analysisManager->GetH1Id("GammaSpecTarget",true),energy/(1.0*keV));
      }
    }


    if ( (momentum.x() && momentum.y()) != 0) {
    /*
    analysisManager->FillNtupleDColumn(1,0,track->GetDynamicParticle()->GetPDGcode());
    analysisManager->FillNtupleDColumn(1,1,track->GetKineticEnergy());
    analysisManager->FillNtupleDColumn(1,2,track->GetTotalEnergy());
    analysisManager->FillNtupleDColumn(1,3,pos.x());
    analysisManager->FillNtupleDColumn(1,4,pos.y());
    analysisManager->FillNtupleDColumn(1,5,pos.z());
    analysisManager->FillNtupleDColumn(1,6,momentum.x());
    analysisManager->FillNtupleDColumn(1,7,momentum.y());
    analysisManager->FillNtupleDColumn(1,8,momentum.z());

    analysisManager->AddNtupleRow(1);
    */
  }


  }

  if (volume == fDetConstruction->GetETT_BackPV()) {

    track = step->GetTrack();
    auto pos = track->GetPosition();
    auto momentum = track->GetMomentum();

    
    //Debugging

    //G4cout << track->GetDefinition()->GetParticleName() << "\n";
    //G4cout << track->GetDynamicParticle()->GetPDGcode() << "\n";
    //G4cout << track->GetKineticEnergy() << "\n";
    //G4cout << track->GetTotalEnergy() << "\n";
    //G4cout << track->GetPosition() << "\n";
    //G4cout << track->GetMomentum() << "\n";
    //G4cout << track->GetMomentumDirection() << "\n";
    

    if (momentum.z() > 0){

      //Ensures particles going into stack, not back reflecting
      partID = track->GetDynamicParticle()->GetPDGcode();
      if (partID == 22){
        auto px = momentum.x();
        auto py = momentum.y();
        auto pz = momentum.z();
        auto energy = sqrt(px*px + py*py + pz*pz);

        analysisManager->FillH1(analysisManager->GetH1Id("GammaSpecStack",true),energy/(1.0*keV));
      }
    /*
    analysisManager->FillNtupleDColumn(2,0,track->GetDynamicParticle()->GetPDGcode());
    analysisManager->FillNtupleDColumn(2,1,track->GetKineticEnergy());
    analysisManager->FillNtupleDColumn(2,2,track->GetTotalEnergy());
    analysisManager->FillNtupleDColumn(2,3,pos.x());
    analysisManager->FillNtupleDColumn(2,4,pos.y());
    analysisManager->FillNtupleDColumn(2,5,pos.z());
    analysisManager->FillNtupleDColumn(2,6,momentum.x());
    analysisManager->FillNtupleDColumn(2,7,momentum.y());
    analysisManager->FillNtupleDColumn(2,8,momentum.z());

    analysisManager->AddNtupleRow(2);
    */

  }
 
  
  }
  
  if (volume->GetName() == IP_NameCheck){
    if(edep != 0){
      track = step->GetTrack();
      auto pos = track->GetPosition();
      auto momentum = track->GetMomentum();
      copyNumber  = volume->GetCopyNo();
      
      auto px = momentum.x();
      auto py = momentum.y();
      auto pz = momentum.z();
      auto energy = sqrt(px*px + py*py + pz*pz);
      partID = track->GetDynamicParticle()->GetPDGcode();
      
      if (partID==11){
        analysisManager->FillH1(analysisManager->GetH1Id("ElStackDist",true),energy/(1.0*keV));
        analysisManager->FillH1(analysisManager->GetH1Id("ProfileElectron",true),copyNumber,edep/(1.0*keV));
      }

      if (partID==-11){
        analysisManager->FillH1(analysisManager->GetH1Id("PosStackDist",true),energy/(1.0*keV));
        analysisManager->FillH1(analysisManager->GetH1Id("ProfilePositron",true),copyNumber,edep/(1.0*keV));
      }
      /*
      analysisManager->FillNtupleDColumn(4,0,pos.x());
      analysisManager->FillNtupleDColumn(4,1,pos.y());
      analysisManager->FillNtupleDColumn(4,2,pos.z());
      analysisManager->FillNtupleDColumn(4,3,edep/(1.0*keV));
      analysisManager->FillNtupleDColumn(4,4,track->GetTotalEnergy()/(1.0*keV));
      analysisManager->FillNtupleDColumn(4,5,track->GetDynamicParticle()->GetPDGcode());
      analysisManager->FillNtupleDColumn(4,6,copyNumber);
      analysisManager->AddNtupleRow(4);
      */


      analysisManager->FillH1(analysisManager->GetH1Id("Profile",true),copyNumber,edep/(1.0*keV));
  }

  }

  //const G4VProcess* process = step->GetPostStepPoint()->GetProcessDefinedStep();
  //G4cout << process << "\n";
  //if (process) G4cout << "Process = " << process->GetProcessName() << "\n";

  // step length
  G4double stepLength = 0.;
  if ( step->GetTrack()->GetDefinition()->GetPDGCharge() != 0. ) {
    stepLength = step->GetStepLength();
  }
      
  if ( volume == fDetConstruction->GetAbsorberPV() ) {
    fEventAction->AddAbs(edep,stepLength);
  }
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

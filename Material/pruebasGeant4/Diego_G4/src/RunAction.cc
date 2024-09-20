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
// $Id: RunAction.cc 100946 2016-11-03 11:28:08Z gcosmo $
//
/// \file RunAction.cc
/// \brief Implementation of the RunAction class

#include "RunAction.hh"
#include "Analysis.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction()
 : G4UserRunAction()
{ 
  // set printing event number per each event
  G4RunManager::GetRunManager()->SetPrintProgress(10000);     

  // Create analysis manager
  // The choice of analysis technology is done via selectin of a namespace
  // in Analysis.hh
  auto analysisManager = G4AnalysisManager::Instance();
  G4cout << "Using " << analysisManager->GetType() << G4endl;

  // Create directories 
  //analysisManager->SetHistoDirectoryName("histograms");
  //analysisManager->SetNtupleDirectoryName("ntuple");
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetNtupleMerging(true);
    // Note: merging ntuples is available only with Root output

  // Book histograms, ntuple
  //
  
  // Creating histograms
  /*
  analysisManager->CreateH1("Eabs","Edep in absorber", 2500, 0., 2500*MeV); //0
  analysisManager->CreateH2("eKin","eKin in pixels", 90,-22.5*cm,22.5*cm,25000,0,2500*MeV);    // 0
  */

  // Create Histograms
  analysisManager->CreateH1("Profile","Total Energy Deposited [keV] in Each Layer",25,-0.5,24.5);
  analysisManager->CreateH1("ProfileElectron","Total Energy Deposited [keV] in Each Layer",25,-0.5,24.5);
  analysisManager->CreateH1("ProfilePositron","Total Energy Deposited [keV] in Each Layer",25,-0.5,24.5);
  analysisManager->CreateH1("GammaSpecTarget","Number of photons per bin",1600,0,800.5);
  analysisManager->CreateH1("GammaSpecStack","Number of photons per bin",1600,0,800.5);
  analysisManager->CreateH1("ElStackDist","Electron distribution in stack",1600,0,800.5);
  analysisManager->CreateH1("PosStackDist","Positron distribution in stack",1600,0,800.5);
  

  // Creating ntuple
  //
  analysisManager->SetFirstNtupleId(0);

  analysisManager->CreateNtuple("AllParticleInfo", "AllTracks"); //0

  analysisManager->CreateNtupleDColumn("Eabs"); //0
  analysisManager->CreateNtupleDColumn("Labs"); //1
  analysisManager->CreateNtupleDColumn("eKin"); //2
  analysisManager->CreateNtupleDColumn("eDep"); //3
  analysisManager->CreateNtupleDColumn("xPos"); //4
  analysisManager->CreateNtupleDColumn("yPos"); //5
  analysisManager->CreateNtupleDColumn("zPos"); //6
  analysisManager->FinishNtuple();

  analysisManager->CreateNtuple("ETT_FrontHits", "HitsAfterTarget"); //1

  analysisManager->CreateNtupleDColumn("pdgID"); //0
  analysisManager->CreateNtupleDColumn("eKin");  //1
  analysisManager->CreateNtupleDColumn("eTot");  //2
  analysisManager->CreateNtupleDColumn("xPos");  //3
  analysisManager->CreateNtupleDColumn("yPos");  //4
  analysisManager->CreateNtupleDColumn("zPos");  //5
  analysisManager->CreateNtupleDColumn("Px");    //6
  analysisManager->CreateNtupleDColumn("Py");    //7
  analysisManager->CreateNtupleDColumn("Pz");    //8
  analysisManager->FinishNtuple();

  analysisManager->CreateNtuple("ETT_BackHits", "HitsBeforeTarget"); //2

  analysisManager->CreateNtupleDColumn("pdgID"); //0
  analysisManager->CreateNtupleDColumn("eKin");  //1
  analysisManager->CreateNtupleDColumn("eTot");  //2
  analysisManager->CreateNtupleDColumn("xPos");  //3
  analysisManager->CreateNtupleDColumn("yPos");  //4
  analysisManager->CreateNtupleDColumn("zPos");  //5
  analysisManager->CreateNtupleDColumn("Px");    //6
  analysisManager->CreateNtupleDColumn("Py");    //7
  analysisManager->CreateNtupleDColumn("Pz");    //8
  analysisManager->FinishNtuple();

  analysisManager->CreateNtuple("PrimaryParticles", "PrimaryParticleInfo"); //3

  analysisManager->CreateNtupleDColumn("pdgID"); //0
  analysisManager->CreateNtupleDColumn("eKin");  //1
  analysisManager->CreateNtupleDColumn("eTot");  //2
  analysisManager->CreateNtupleDColumn("xPos");  //3
  analysisManager->CreateNtupleDColumn("yPos");  //4
  analysisManager->CreateNtupleDColumn("zPos");  //5
  analysisManager->CreateNtupleDColumn("Px");    //6
  analysisManager->CreateNtupleDColumn("Py");    //7
  analysisManager->CreateNtupleDColumn("Pz");    //8
  analysisManager->FinishNtuple();

  analysisManager->CreateNtuple("IP_Hits", "IP_HITS"); //4

  analysisManager->CreateNtupleDColumn("xPos");  //0
  analysisManager->CreateNtupleDColumn("yPos");  //1
  analysisManager->CreateNtupleDColumn("zPos");  //2
  analysisManager->CreateNtupleDColumn("eDep");  //3
  analysisManager->CreateNtupleDColumn("eTotal");//4
  analysisManager->CreateNtupleDColumn("pdgID"); //5
  analysisManager->CreateNtupleDColumn("copyNo");//6
  analysisManager->FinishNtuple();


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
  delete G4AnalysisManager::Instance();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* /*run*/)
{ 
  //inform the runManager to save random number seed
  //G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  
  // Get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();

  // Open an output file
  //
  const G4String& fileName = "";
  analysisManager->OpenFile(fileName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* /*run*/)
{
  // print histogram statistics
  //
  auto analysisManager = G4AnalysisManager::Instance();

  // save histograms & ntuple
  //
  analysisManager->Write();
  analysisManager->CloseFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

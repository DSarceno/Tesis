from ROOT import *
import math
import numpy as np
import csv

def CreateCSV(Hist):
	nBins = Hist.GetBin
def PrintHitInfo(hitInfo):

	keys = hitInfo.keys()
	for key in keys:

		noOfCounts = hitInfo[key]["Counts"]

		if noOfCounts != 0:
			print "~~~~~~~~~~~~"
			print "Particle: ", key
			print "Counts: ", noOfCounts
			print "MeanEnergy: ", hitInfo[key]["MeanEnergy"]/noOfCounts, "[MeV]"
			print ""

ShotNo = raw_input("What file name? ")
#filePath="/home/ut3group/Desktop/Andrea/Single_Stack_Sim_HZDR/Output/Electron response/%s.root" %ShotNo
filePath="/home/ut3group/Desktop/Andrea/Single_Stack_Sim_HZDR/Output/%s.root" %ShotNo

rFile = TFile(filePath)

# Get Trees

BackHits = rFile.Get('ETT_BackHits')


# How many events in file (1 event = 1 incedent photon)

nBackHits = BackHits.GetEntries()

minE = -1
maxE = 600
binSize = .5
nBins = (maxE-minE)/binSize



#GammaDist_2D = TH2F("Gamma_Dist2D","Gamma Radial Distribution - %s; Gamma Energy [MeV] ; Radial Dist [mm]" %ShotNo, nBins,minE,maxE,200,0,1000) 
#GammaDist_Theta = TH2F("Gamma_Dist_T","Gamma Radial Distribution - %s; Gamma Energy [MeV] ; Angle [rad]" %ShotNo, nBins,minE,maxE,400,0,3.142) 

GammaDist = TH1F("GammaDist","Gamma Energy Distribution - %s; Energy [MeV]" %ShotNo, int(nBins),minE,maxE)
#ElDist = TH1F("ElDist","Electron Energy Distribution - %s; Energy [MeV]" %ShotNo, int(nBins),minE,maxE)
#PosDist = TH1F("PosDist","Positron Energy Distribution - %s; Energy [MeV]" %ShotNo, int(nBins),minE,maxE)
#OtherDist = TH1F("OtherDist","Other Energy Distribution - %s; Energy [MeV]" %ShotNo, int(nBins),minE,maxE)


HitInfo = {"e-":{"Counts":0,"MeanEnergy":0}, "e+":{"Counts":0,"MeanEnergy":0}, "gamma":{"Counts":0,"MeanEnergy":0}, "N":{"Counts":0,"MeanEnergy":0}, "P":{"Counts":0,"MeanEnergy":0}, "Mu-":{"Counts":0,"MeanEnergy":0}, "Mu+":{"Counts":0,"MeanEnergy":0}, "Pi0":{"Counts":0,"MeanEnergy":0}, "Pi+":{"Counts":0,"MeanEnergy":0}, "Pi-":{"Counts":0,"MeanEnergy":0}}
NeutrinoHits = 0
OtherHits = 0
OtherHitIDs = []

ParticleIDs = {11:"e-", -11:"e+", 22:"gamma", 2112:"N", 2212:"P", 13:"Mu-", -13:"Mu+",111:"Pi0",211:"Pi+",-211:"Pi-"}
NeutrinoIDs = [12,14]
Cut = 3.2

for hitNo in range(nBackHits):

	if hitNo%100000 == 0:
		print "Processing Front Hit: ", hitNo, "/", nBackHits

	#Get the Hit
	BackHits.GetEntry(hitNo)

	px = BackHits.GetLeaf('Px').GetValue()
	py = BackHits.GetLeaf('Py').GetValue()
	pz = BackHits.GetLeaf('Pz').GetValue()

	xPos = BackHits.GetLeaf('xPos').GetValue()
	yPos = BackHits.GetLeaf('yPos').GetValue()

	r = math.sqrt(xPos**2 + yPos**2)

	Energy = np.sqrt(px**2 + py**2 + pz**2)

	pdgID = BackHits.GetLeaf('pdgID').GetValue()
	if pdgID in ParticleIDs.keys():
		HitInfo[ParticleIDs[pdgID]]["Counts"] += 1
		HitInfo[ParticleIDs[pdgID]]["MeanEnergy"] += Energy

	elif abs(pdgID) in NeutrinoIDs:
		NeutrinoHits += 1

	else:
		OtherHitIDs.append(pdgID)
		OtherHits += 1

	if pdgID == 22:
		#GammaDist_2D.Fill(Energy,r)
		#GammaDist_Theta.Fill(Energy,np.arctan(r/200))
		GammaDist.Fill(Energy)

	#elif pdgID == 11:
		#ElDist.Fill(Energy)

	#elif pdgID == -11:
		#PosDist.Fill(Energy)

	#else:
		#OtherDist.Fill(Energy)

PrintHitInfo(HitInfo)

print "Neutrino Hits: ", NeutrinoHits
print "Unknown Hits: ", OtherHits, " with IDs: ", OtherHitIDs

#ElectronDist_Theta.Draw('colz')

#ElProfile.Draw('same')

#raw_input('Press enter to continue....')

#PositronDist_Theta.Draw('colz')

#PosProfile.Draw('same')

#raw_input('Press enter to quit....')
#GammaEnergyDist = GammaDist_Theta.ProjectionX()

#GammaDist_2D.SaveAs("Plots/GammaDist_2D_%s.root" %ShotNo )
#GammaDist_Theta.SaveAs("Plots/GammaDistTheta_%s.root"%ShotNo )
#GammaEnergyDist.SaveAs("Plots/GammaEnergyDist_%s.root"%ShotNo)
#GammaDist.SaveAs("Plots/Electron response/GammaDist_%s.root" %ShotNo )
GammaDist.SaveAs("Plots/GammaDist_%s.root" %ShotNo )
#ElDist.SaveAs("Plots/ElDist_%s.root" %ShotNo )
#PosDist.SaveAs("Plots/PosDist_%s.root" %ShotNo )
#OtherDist.SaveAs("Plots/OtherDist_%s.root" %ShotNo )



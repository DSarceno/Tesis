from ROOT import *
import math
import csv
import numpy as numpy


ShotNo   = raw_input("What Shot Number? ")
filePath = "/home/ut3group/Desktop/Andrea/Single_Stack_Sim_HZDR/Output/%s.root" %ShotNo


rFile = TFile(filePath)

# Get Trees

IPHits = rFile.Get('IP_Hits')
nIPHits = IPHits.GetEntries()

nLayers = 25;

histos = [];
sMax = 40;
binSize = 0.2;

nBins = int((sMax*2)/binSize)

profile_h = TH1F("StackProfile", "; LayerNo; Total Energy Deposited [keV]", 25, -.5, 24.5)

profile2D_h = TH2F("StackProfile_2D","; LayerNo; Radial Pos [mm]", 25,-.5,24.5,350,0,35)

for layer in range(nLayers):
	histos.append(TH2F("StackLayer_%s"%layer,"; xPos [mm] ; yPos [mm]", nBins,-1*sMax,sMax,nBins,-1*sMax,sMax) 
)

for hitNo in range(nIPHits):

	if hitNo%10000 == 0:
		print "Processing Front Hit: ", hitNo, "/", nIPHits

	IPHits.GetEntry(hitNo)

	copyNo = IPHits.GetLeaf('copyNo').GetValue()
	eDep   = IPHits.GetLeaf('eDep').GetValue()
	xPos   = IPHits.GetLeaf('yPos').GetValue()
	yPos   = IPHits.GetLeaf('xPos').GetValue()

	r = math.sqrt(xPos**2+yPos**2)

	r_set = 35
	if r < r_set: #1D edep only in area of image plates 
		profile_h.Fill(copyNo,eDep)

	profile2D_h.Fill(copyNo,r,eDep)

	histos[int(copyNo)].Fill(xPos,yPos,eDep)




#Adding from root_to_csv.py to do all in one
xLabel = "Layer"
yLabel = "Energy deposited [keV]"
nBins  = profile_h.GetNbinsX()
print("number of bins = ",str(nBins))


CSV_Dir  = "/home/ut3group/Desktop/Andrea/Single_Stack_Sim_HZDR/Plots/Profile_"+ShotNo+"_r"+str(r_set)+".csv"
with open(CSV_Dir, 'wb') as csvfile:

	writer = csv.writer(csvfile, delimiter=',', quotechar="|")

	for binNo in range(nBins+1):
		if binNo == 0:
			writer.writerow([xLabel,yLabel])
			continue 

		binCenter = profile_h.GetBinCenter(binNo)
		binHeight = profile_h.GetBinContent(binNo)
		writer.writerow([binCenter, binHeight])

#End of added

n = 0
for histo in histos:
	n = n + 1
	#histo.Draw('colz')
	#histo.SaveAs('/Users/Luc/Documents/Research/Geant4SimulationStudies/Geant4Souce/Geant4Sim/SourceCode/Stack_Sim/Plots/Shot_%s-Layer_%s.root' %(ShotNo,n))
	#raw_input("Stack Layer %s..." %n)

#print(profile_h.shape)
profile_h.SaveAs('/home/ut3group/Desktop/Andrea/Single_Stack_Sim_HZDR/Plots/Profile_'+ShotNo+'_r'+str(r_set)+'.root' )

profile2D_h.SaveAs('/home/ut3group/Desktop/Andrea/Single_Stack_Sim_HZDR/Plots/Profile2D_'+ShotNo+'.root' )



import csv
import numpy as np
from ROOT import *

ShotNo = raw_input("What file? ")

#Directory names
# For gamma ray distributions

"""
CSV_Dir = "/home/ut3group/Desktop/Andrea/Single_Stack_Sim_HZDR/Plots/GammaDist_%s.csv" %ShotNo
RootDir = "/home/ut3group/Desktop/Andrea/Single_Stack_Sim_HZDR/Plots/GammaDist_%s.root" %ShotNo

histName = "GammaDist"
xLable = "Energy [MeV]"
yLable = "dNdE [#/MeV]"
"""

# For stack profile distributions

histName = "ProfilePositron"
xLable = "Layer"
yLable = "Energy deposited [keV]"

CSV_Dir = "/home/ut3group/Desktop/Andrea/Single_Stack_Sim_HZDR/Output/"+histName+"_%s.csv" %ShotNo
RootDir = "/home/ut3group/Desktop/Andrea/Single_Stack_Sim_HZDR/Output/%s.root" %ShotNo




print ""
print "Starting..."
print ""

rFile = TFile(RootDir)
myHist = rFile.Get(histName)

#print myHist.GetBinContent(600)
#print myHist.GetBinCenter(600)
nBins = myHist.GetNbinsX()

with open(CSV_Dir, 'wb') as csvfile:

	writer = csv.writer(csvfile, delimiter=',', quotechar="|")

	for binNo in range(0,nBins):
		if binNo == 0:
			writer.writerow([xLable,yLable])
			continue 

		binCenter = myHist.GetBinCenter(binNo)
		binHeight = myHist.GetBinContent(binNo)
		writer.writerow([binCenter+1, binHeight])



print "ROOT to CSV conversion succesful!"
print ""
print "File saved to: %s" %CSV_Dir
print ""

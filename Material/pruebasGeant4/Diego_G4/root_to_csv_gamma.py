import csv
import numpy as np
from ROOT import *

ShotNo = raw_input("What file? ")

#Directory names
# For gamma ray distributions


distname = "GammaSpecStack"

CSV_Dir = "/home/ut3group/Desktop/Andrea/Single_Stack_Sim_HZDR/Plots/"+distname+"_%s.csv" %ShotNo
RootDir = "/home/ut3group/Desktop/Andrea/Single_Stack_Sim_HZDR/Output/%s.root" %ShotNo

histName = distname
xLable = "Energy [MeV]"
yLable = "dNdE [#/MeV]"
"""

CSV_Dir = "/home/ut3group/Desktop/Andrea/Stack_Sim_LWFA8/Plots/"+distname+"Dist_%s.csv" %ShotNo
RootDir = "/home/ut3group/Desktop/Andrea/Stack_Sim_LWFA8/Plots/"+distname+"Dist_%s.root" %ShotNo

histName = distname+"Dist"
xLable = "Energy [MeV]"
yLable = "dNdE [#/MeV]"
"""
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
		writer.writerow([binCenter, binHeight])



print "ROOT to CSV conversion succesful!"
print ""
print "File saved to: %s" %CSV_Dir
print ""

import csv
import numpy as np
from ROOT import *

ShotNo = raw_input("What file? ")

#Directory names
# For gamma ray distributions

CSV_Dir = "/home/ut3group/Desktop/Andrea/Single_Stack_Sim_HZDR/Plots/Profile2D_%s.csv" %ShotNo
RootDir = "/home/ut3group/Desktop/Andrea/Single_Stack_Sim_HZDR/Plots/Profile2D_%s.root" %ShotNo

histName = "StackProfile_2D"
xLabel = "Layer"
yLabel = "Radius [mm]"
zlabel = "Energy deposited [keV]"

"""
# For stack profile distributions

CSV_Dir = "/home/ut3group/Desktop/Andrea/Single_Stack_Sim_HZDR/Output/%s.csv" %ShotNo
RootDir = "/home/ut3group/Desktop/Andrea/Single_Stack_Sim_HZDR/Output/%s.root" %ShotNo

histName = "Profile"
xLabel = "Layer"
yLabel = "Radius [mm]"
zlabel = "Energy deposited [keV]"
"""

print ""
print "Starting..."
print ""

rFile = TFile(RootDir)
myHist = rFile.Get(histName)

#print myHist.GetBinContent(600)
#print myHist.GetBinCenter(600)
nBinsX = myHist.GetNbinsX()
nBinsY = myHist.GetNbinsY()

print("nBinsx",nBinsX)
print("nBinsy",nBinsY)
print(myHist)


with open(CSV_Dir, 'wb') as csvfile:

	writer = csv.writer(csvfile, delimiter=',', quotechar="|")

	writer.writerow([xLabel,yLabel,zlabel])

	for binNoX in range(1,nBinsX):
		for binNoY in range(1,nBinsY):

			del_r = .1 #in mm
			area = 2*3.14*del_r*binNoY/10.
			
			binCenter = myHist.GetBin(binNoX,binNoY)
			#print(binNoX,binNoY,binCenter)
			binHeight  = myHist.GetBinContent(binCenter)
			writer.writerow([binNoX,binNoY/10.,binHeight])



print "ROOT to CSV conversion succesful!"
print ""
print "File saved to: %s" %CSV_Dir
print ""

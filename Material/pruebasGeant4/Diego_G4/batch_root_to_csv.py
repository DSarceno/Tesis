import csv
import numpy as np
from ROOT import *
import fnmatch
import glob, os

root_path = "/home/ut3group/Desktop/Andrea/Stack_Sim_HZDR/Output/Response_matrix_raw_5e5/"
Pnum = 5e5    # If it didnt take it into account first

Estart   = 5
Eend     = 200
spacing  = 1 # energy spacing in keV


Ebins    = np.arange(Estart,Eend+spacing,spacing)

for file in os.listdir(root_path):
	if fnmatch.fnmatch(file,"*.root"):
		energy = file[8:-8]
		#Ebins[np.round(int(energy)/spacing)-1] = energy


		#Directory names
		#CSV_Dir = "/home/ut3group/Desktop/Andrea/Stack_Sim_LWFA8/Output/Response_matrix_raw/Test%s.csv" %ShotNo
		#RootDir = "/home/ut3group/Desktop/Andrea/Stack_Sim_LWFA8/Output/Response_matrix_raw/Test%s.root" %ShotNo

		histName = "Profile"
		xLable = "Layer"
		yLable = "Energy deposited [keV]"

		"""
		print ""
		print "Converting = "
		print ""
		"""

		rFile = TFile(root_path+file)
		myHist = rFile.Get(histName)

		#print myHist.GetBinContent(600)
		#print myHist.GetBinCenter(600)
		nBins = myHist.GetNbinsX()
		CSV_Dir = root_path+file[:-4]+"csv"
		
		with open(CSV_Dir, 'wb') as csvfile:

			writer = csv.writer(csvfile, delimiter=',', quotechar="|")

			for binNo in range(0,nBins):
				if binNo == 0:
					writer.writerow([xLable,yLable])
					continue 

				binCenter = myHist.GetBinCenter(binNo)
				binHeight = myHist.GetBinContent(binNo)
				writer.writerow([binCenter+1, binHeight/Pnum])


np.savetxt("/home/ut3group/Desktop/Andrea/Stack_Sim_HZDR/Ebins_200keV.txt",Ebins)
print "ROOT to CSV conversion succesful!"
print ""
#print "File saved to: %s" %CSV_Dir
#print ""

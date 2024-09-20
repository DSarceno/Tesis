import csv
import numpy as np
from ROOT import *
import fnmatch
import glob, os

root_path = "/home/ut3group/Desktop/Andrea/Stack_Sim_HZDR/Output/Response_matrix_raw_5e5/"
spacing  = 1 #energy spacing in keV

Ebins    = np.loadtxt("/home/ut3group/Desktop/Andrea/Stack_Sim_HZDR/Ebins_200keV.txt")
response = np.zeros((24,len(Ebins)))
#response = np.loadtxt("/home/ut3group/Desktop/Andrea/Stack_Sim_HZDR/response_100000keV_mix.txt")


for ei in range(len(Ebins)):
	for file in os.listdir(root_path):
		if fnmatch.fnmatch(file,"Profile_"+str(int(Ebins[ei]))+"keV.csv"):
			#print(Ebins[ei])
			temp = open(root_path+file,"r")
			contents = temp.readlines()[1:]
			temp.close()
			for li in range(len(contents)):
				
				line = contents[li].split(",")
				response[li,ei] = float(line[1][:-1])


np.savetxt("/home/ut3group/Desktop/Andrea/Stack_Sim_HZDR/response_200keV_5e5.txt",response)


		#Directory names
		#CSV_Dir = "/home/ut3group/Desktop/Andrea/Stack_Sim_LWFA8/Output/Response_matrix_raw/Test%s.csv" %ShotNo
		#RootDir = "/home/ut3group/Desktop/Andrea/Stack_Sim_LWFA8/Output/Response_matrix_raw/Test%s.root" %ShotNo

		



print "Response matrix successfully built!"
print ""
#print "File saved to: %s" %CSV_Dir
#print ""

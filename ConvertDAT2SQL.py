#myFile = open("Run_171900.dat","rb")
from struct import *
import re

with open("Run_171900.dat","rb") as fin:
	header = fin.read(4)

	for j in range(0,3):
		brHeader= fin.read(4)
		if not re.search('B',brHeader):
			fin.seek(-4,1)
			break
		brSerial = unpack('H',brHeader[2:])
		print "Reading Board...",brHeader[:2],brSerial[0]

		for i in range(0,4):
			chHeader = fin.read(4)
			if not re.search('C',chHeader):
				fin.seek(-4,1)
				break

			print "Read channel...",chHeader
			binData = fin.read(4096)
			binWidth = unpack('1024f',binData)
	
				


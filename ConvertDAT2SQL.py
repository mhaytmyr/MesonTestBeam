#myFile = open("Run_171900.dat","rb")
from struct import *
import re

with open("Run_171900.dat","rb") as fin:
	header = fin.read(4)
	print header
	brHeader= fin.read(4)
	brSerial = unpack('H',brHeader[2:])
	print brHeader[:2],brSerial[0]

	#for i in range(0,4):
	chHeader = fin.read(4)
	binWidth = unpack('1024f',chHeader[4:])
	#print chHeader[:2],binWidth
	
	if re.search('C',chHeader):
		print chHeader[:4]	
		#if chHeader[1] not 'C':
				


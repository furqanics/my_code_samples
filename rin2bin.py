#!/usr/bin/python

import os
#os.system("module load python"); #Uncomment/Enable this line when running on TACC

import numpy
import sys
import commands
import datetime
from scipy.io import FortranFile

#rin2bin.py
#put descriptive notes here

#parse input arguments

#This program requires that FreqToRead.txt file exists in the same directory as it is in

def getSysToRead():
	#This function reads the FreqToRead.txt file and returns the list of GNSS to be read
	FreqFile = open('FreqToRead.txt','r').readlines();
	SysToRead = [];
	for line in FreqFile:
		if(line[0:1]!="#"):
			SysToRead.append(line.split()[0]);
	return SysToRead;

def getFreqToRead(sysFlag):
	#This function reads the FreqToRead.txt file and returns the observables to be read for a specified GNSS
	FreqFile = open('FreqToRead.txt','r').readlines();
	FreqToRead = [];
	for line in FreqFile:
		if(line[0:1]!="#"):
			if line.split()[0]==sysFlag:
				FreqToRead = line.split()[1:]
	return FreqToRead;

def getRinVer(RinFile):
	#This function reads the RINEX version from the header
	inFile = open(RinFile,'r').readlines();
	ver = 0;
	for line in inFile:
		if "RINEX VERSION / TYPE" in line[60:]:
			ver = float(line[0:9]);
	return ver;

def getRinFileType(RinFile):
	#This function reads the type of the RINEX file (e.g. O:observation, N:navigation, etc.) from the header
	inFile = open(RinFile,'r').readlines();
	type = 'O'; #default: Observation file
	for line in inFile:
		if "RINEX VERSION / TYPE" in line[60:]:
			type = line[20:21];
	return type;

def getRinSysType(RinFile):
	#This function reads the flag for type of observations (e.g. M:mixed, G:GPS, etc.)  from the header
	inFile = open(RinFile,'r').readlines();
	type = 'G'; #default: GPS only
	for line in inFile:
		if "RINEX VERSION / TYPE" in line[60:]:
			type = line[40:41];
	return type;
	
def getRinPGM(RinFile):
	#This function reads the PGM / RUN BY / DATE label value from the header
	inFile = open(RinFile,'r').readlines();
	pgm = 'XXXX';
	for line in inFile:
		if "PGM / RUN BY / DATE" in line[60:]:
			pgm = line[0:60];
	return pgm;

def getRinComments(RinFile):
	#This function reads all the COMMENTS from RINEX header and puts them in an array
	inFile = open(RinFile,'r').readlines();
	commentsArray = [];
	for line in inFile:
		if "COMMENT" in line[60:]:
			commentsArray.append(line[0:60]);
	return commentsArray;

def getRinMarkerName(RinFile):
        #This function reads the station/marker name from the RINEX file
        inFile = open(RinFile,'r').readlines();
        markerName = 'XXXX';
        for line in inFile:
                if "MARKER NAME" in line[60:]:
                        markerName = line[0:60];
        return markerName;

def getRinMarkerNumber(RinFile):
	#This function reads the marker number from the RINEX header
	inFile = open(RinFile,'r').readlines();
	markerNumber = 'XXXX';
	for line in inFile:
		if "MARKER NUMBER" in line[60:]:
			markerNumber = line[0:20];
	return markerNumber;
		
def getRinMarkerType(RinFile):
	#This function reads the marker type from the RINEX file
	inFile = open(RinFile,'r').readlines();
	markerType = 'XXXX';
	for line in inFile:
		if "MARKER TYPE" in line[60:]:
			markerType = line[0:20];
	return markerType;

def getObsAgency(RinFile):
	#This function reads the Observer / Agency from the RINEX file
	inFile = open(RinFile,'r').readlines();
	Observer = 'XXXX';Agency = 'XXXX';
	for line in inFile:
		if "OBSERVER / AGENCY" in line[60:]:
			Observer = line[0:20];
			Agency = line[20:60];
	return Observer, Agency;

def getRcvInfo(RinFile):
	#This function reads the Receiver Information from the RINEX file
	inFile = open(RinFile,'r').readlines();
	RcvSN = '0000'; #Receiver Serial Number
	RcvType = 'XXXX'; #Receiver Type
	RcvVer = '0000'; #Receiver Version
	for line in inFile:
		if "REC # / TYPE / VERS" in line[60:]:
			RcvSN = line[0:20];
			RcvType = line[20:40];
			RcvVer = line[40:60];
	return RcvSN, RcvType, RcvVer;

def getAntInfo(RinFile):
	#This function reads the Antenna Information from the RINEX file
	inFile = open(RinFile,'r').readlines();
	AntSN = '0000'; #Antenna Serial Number
	AntType = 'XXXX'; #Antenna Type
	for line in inFile:
		if "ANT # / TYPE" in line[60:]:
			AntSN = line[0:20];
			AntType = line[20:40];
	return AntSN, AntType;
	
def getRinXYZ(RinFile):
	#This function reads the Station XYZ from the RINEX file
	inFile = open(RinFile,'r').readlines();
	staX = 0; #Approx X
	staY = 0; #Approx Y
	staZ = 0; #Approx Z
	for line in inFile:
		if "APPROX POSITION XYZ" in line[60:]:
			staX = float(line[0:14]);
			staY = float(line[14:28]);
			staZ = float(line[28:42]);
	return staX,staY,staZ;
	
def getRinEcc(RinFile):
	#This function reads the ARP eccentricity from the RINEX file
	inFile = open(RinFile,'r').readlines();
	ARPh = 0; #ARP deltaH
	ARPe = 0; #ARP deltaE
	ARPn = 0; #ARP deltaN
	for line in inFile:
		if "ANTENNA: DELTA H/E/N" in line[60:]:
			ARPh = float(line[0:14]);
			ARPe = float(line[14:28]);
			ARPn = float(line[28:42]);
	return ARPh,ARPe,ARPn;
	

def getRinSysNumTypes(RinFile):
	#This number reads the number and types of observables for each available GNSS in the RINEX file
	inFile = open(RinFile,'r').readlines();
	availableSys = []; #Array containing flags for all available GNSS in the RINEX file
	availableObsvs = []; #Array of arrays of available observation types for each GNSS
	tempIndices = [];
	for line in inFile:
		if "SYS / # / OBS TYPES" in line[60:]:
			tempIndices.append(inFile.index(line));
			tempIndices.sort();
	startIndex = numpy.min(tempIndices); #startIndex of SYS / # / OBS TYPES block
	endIndex = numpy.max(tempIndices); #endIndex of SYS / # / OBS TYPES block
	for rec in inFile[startIndex:endIndex+1]:
		if "SYS / # / OBS TYPES" in rec[60:]:
			if rec[0:1]!=' ':
				availableSys.append(rec[0:1]);
	for GNSS in availableSys:
		thisObsvs = []; #Types of observables available for this GNSS
		thisBlock = inFile[startIndex:endIndex+1];
		for rec in thisBlock:
			if rec[0:1] == GNSS:
				numTypes = int(rec[3:6]); #Number of types of observables for a GNSS
				numLines = (numTypes/13)+1; #Number of lines for a GNSS in header (1 line can carry maximum 13)
				for l in range(thisBlock.index(rec),thisBlock.index(rec)+numLines):
					for obsv in thisBlock[l][7:60].split():
						thisObsvs.append(obsv);
		availableObsvs.append(thisObsvs);
	return availableSys, availableObsvs; #These two will be useful for reading epoch by epoch data for each system i.e. the order of data in a row
	
def getRinSysObs(GNSS, RinFile):
	#This function returns the list of available observables for a specific GNSS (e.g. 'G', 'R' or 'E') in the RINEX file
	#print getRinSysObs.__name__;
	SysObs = [];
	systems, obslist = getRinSysNumTypes(RinFile);
	if (GNSS in systems):
		SysObs = obslist[systems.index(GNSS)];
	return SysObs;

def getNumOfEpochs(RinFile):
	#This function returns the number of data epochs in a given RINEX obs. file
	inFile = open(RinFile,'r').readlines();
	NumOfEpochs = 0;
	for line in inFile:
		if(line[0:1]==">"):NumOfEpochs+=1;
	return NumOfEpochs;

def getTimeOfFirstObs(RinFile):
	#This function reads the time of first observation from the RINEX file
	inFile = open(RinFile,'r').readlines();
	yyyy = 0;
	mm = 0;
	dd = 0;
	HH = 0;
	MM = 0;
	SS = 0;
	timeSys = "GPS";
	for line in inFile:
		if "TIME OF FIRST OBS" in line[60:]:
			yyyy = int(line[0:6]);
			mm = int(line[6:12]);
			dd = int(line[12:18]);
			HH = int(line[18:24]);
			MM = int(line[24:30]);
			SS = float(line[30:43]);
			timeSys = line[48:51];
	return yyyy,mm,dd,HH,MM,SS,timeSys;
	
def getSysPhaseShift(RinFile):
	#This function reads the "SYS / PHASE SHIFT" records from the Rinex header
	inFile = open(RinFile,'r').readlines();
	sysPhaseShifts = []; #Array of arrays containing these tuples: [system, signalID, corrValue, satsToApply]
	tempIndices = [];
	for line in inFile:
		if "SYS / PHASE SHIFT" in line[60:]:
			tempIndices.append(inFile.index(line));
			tempIndices.sort();
	startIndex = numpy.min(tempIndices); #startIndex of SYS / PHASE SHIFT block
	endIndex = numpy.max(tempIndices); #endIndex of SYS / PHASE SHIFT block
	thisBlock = inFile[startIndex:endIndex+1];
	for rec in thisBlock:
		if "SYS / PHASE SHIFT" in rec[60:]:
			if rec[0:1]!=' ':
				tempphase = [];
				tempphase.append(rec[0:1]);
				tempphase.append(rec[2:5]);
				corrVal = rec[6:14];
				tempphase.append(corrVal);
				satsToApply = "ALLSAT"; #By default, suppose that the phase correction is to be applied to the observables of all the satellites in this system
				#Or if the satellite numbers are given, satsToApply will be the array containing satellite names to be corrected for phase
				if(rec[16:18]!=" " or "0" not in rec[16:18]):
					satsToApply = [];
					numSat = int(rec[16:18]);
					numLines = (numSat/10)+1;
					for l in range(thisBlock.index(rec),thisBlock.index(rec)+numLines):
						for each in thisBlock[l][18:58].split():
							satsToApply.append(each);
				for each in satsToApply:
					tempphase.append(each);
				sysPhaseShifts.append(tempphase);
	return sysPhaseShifts;
	
def getGloSlotFreq(RinFile):
	GloSlotFreq = []; #Array of arrays i.e. [["R",satNum,freqNum],...]
	inFile = open(RinFile,'r').readlines();
	for line in inFile:
		if("GLONASS SLOT / FRQ #" in line[60:]):
			if(int(line[0:3])==0):
				return GloSlotFreq;
	numSat = int(line[0:3]);
	if(numSat!=0):
		numLines = (numSat/8)+1;
		tempIndices = [];
		for line in inFile:
			if "GLONASS SLOT / FRQ #" in line[60:]:
				tempIndices.append(inFile.index(line));
				tempIndices.sort();
		startIndex = numpy.min(tempIndices); #startIndex of GLONASS SLOT / FRQ # block
		endIndex = numpy.max(tempIndices); #endIndex of GLONASS SLOT / FRQ # block
		thisBlock = inFile[startIndex:endIndex+1];
		for i in range(startIndex,endIndex+1):
			for j in range(0,8):
				tempsat = [];
				tempsat.append(inFile[i][0+(j*7):3+(j*7)]);
				tempsat.append(inFile[i][4+(j*7):6+(j*7)]);
				GloSlotFreq.append(tempsat);
	return GloSlotFreq;

def getGloCPB(RinFile):
	#This function reads the GLONASS COD/PHS/BIS record from RINEX header
	GloCPB = []; #Array of C1c/code-phase, C1P/code-phase, C2C/code-phase, and C2P/code-phase bias correction values [meters]
	inFile = open(RinFile,'r').readlines();
	for line in inFile:
		if("GLONASS COD/PHS/BIS" in line[60:]):
			for i in range(0,4):
				temp = [];
				temp.append(line[1+(i*13):4+(i*13)]);
				temp.append(line[5+(i*13):13+(i*13)]);
				GloCPB.append(temp);
	return GloCPB;

def headerEnds(RinFile):
	EOHfound = False;
	inFile = open(RinFile,'r').readlines();
	for line in inFile:
		if("END OF HEADER" in line[60:]):
			EOHfound = True;
			break;
	return EOHfound;

def readRinData(RinFile):
	SysToRead = getSysToRead();
	thisFileData = []; #List of epoch-wise lists of dictionaries to store this file's data where each dictionary contains data from one observation line
	#This function reads epoch by epoch observations and returns each record as a data structure
	inFile = open(RinFile,'r').readlines();
	startIndexes = []; #List of indexes where an epoch starts
	for line in inFile:
		if(line[0:1]==">"):startIndexes.append(inFile.index(line));
	for epoch in startIndexes: #Read the obs. data epoch-by-epoch
		#Parse the time information line, read its values:
		YYYY = int(inFile[epoch][2:6]);#Year (4-digit)
		MM = int(inFile[epoch][7:9]);#Month (2-digit)
		DD = int(inFile[epoch][10:12]);#Day (2-digit)
		hh = int(inFile[epoch][13:15]);#Hour (2-digit)
		mm = int(inFile[epoch][16:18]);#Minute (2-digit)
		ss = float(inFile[epoch][18:29]);#Second (F11.7)
		epochFlag = int(inFile[epoch][29:32]);#Epoch Flag (I1)
		NumSat = int(inFile[epoch][32:35]);#Number of satellites covered in this epoch
		reserved = inFile[epoch][35:41];
		RcvClkOffset = inFile[epoch][41:56];
		print "Reading epoch",YYYY,MM,DD,hh,mm,ss;
		
		if(epochFlag == 0 or epochFlag == 1): #If the epoch record contains observations
			thisEpochDicts = []; #List of line-wise dictionaries for this epoch
			indDataStart = epoch+1;
			indDataEnd = epoch+NumSat;
			thisEpochRawData = inFile[indDataStart:indDataEnd+1];#Array to store raw observation lines for this epoch
		
			for line in thisEpochRawData:
				IGNSSF = line[0]; #GNSS Flag / Identifier
				if(IGNSSF in SysToRead):
					ObsInThisLine = 0;#Number of observables in this line
					lineTemp = line[3:];
					lineTempInd = 0;
					lineStartInd = 0;
					while(lineTemp[lineTempInd].isspace()):
						lineTempInd+=1;
						if(lineTemp[lineTempInd+1].isspace()==False):
							lineStartInd = lineTempInd+1;
							break;
					#print lineTemp,"len:",len(lineTemp)%16,"lineTempInd",;
					lineTemp = " "*((lineTempInd+1)%16)+lineTemp[lineStartInd:];
					ObsCountStart = (((lineTempInd+1)/16)+1)-1;
					print "ObsCountStart: ",ObsCountStart;
					if (((len(lineTemp))%16)==0):ObsInThisLine = (len(lineTemp))/16;
					if (((len(lineTemp))%16)>0):ObsInThisLine = ((len(lineTemp))/16)+1;
					SatNum = line[1:3]; #Satellite number
					ObsTypes = getRinSysObs(IGNSSF, RinFile); #Obs. types available for this GNSS
					print IGNSSF,SatNum,ObsInThisLine,ObsTypes;
					#print "len(line)-3= ",len(line)-3,(len(line)-3)/16,len(ObsTypes);
					thisLineData = dict(); #Dictionary to store this line data
					thisLineData['YYYY'] = YYYY;
					thisLineData['MM'] = MM;
					thisLineData['DD'] = DD;
					thisLineData['hh'] = hh;
					thisLineData['mm'] = mm;
					thisLineData['ss'] = ss;
					thisLineData['epochFlag'] = epochFlag;
					thisLineData['NumSat'] = NumSat;
					thisLineData['reserved'] = reserved;
					thisLineData['IGNSSF'] = IGNSSF;
					thisLineData['SatNum'] = SatNum;
					thisLineData['RCBIAS'] = RcvClkOffset;
					
					for i in range(ObsCountStart,ObsInThisLine): #stuck here, august 15 5:10pm, changing from 0 to ObsCountStart causing error
						ObsType = ObsTypes[i];
						print "ObsType",ObsType;
						if (lineTemp[(16*ObsTypes.index(ObsType)):14+(16*ObsTypes.index(ObsType))]).isspace()==False:
							thisLineData[ObsType] = float(lineTemp[(16*ObsTypes.index(ObsType)):14+(16*ObsTypes.index(ObsType))]);
						if (lineTemp[14+(16*ObsTypes.index(ObsType))]).isspace()==False:
							thisLineData[ObsType+"_LLI"] = float(lineTemp[14+(16*ObsTypes.index(ObsType))]);
						#print "INDEX",17+(16*ObsTypes.index(ObsType))+1,"LEN",len(line);
						try:
							if (lineTemp[14+(16*ObsTypes.index(ObsType))+1]).isspace()==False:
								thisLineData[ObsType+"_SS"] = float(lineTemp[14+(16*ObsTypes.index(ObsType))+1]);
						except IndexError:
							continue;
					#print thisLineData;
					#print " ";
					if(ObsType in getFreqToRead(IGNSSF)):
						print thisLineData;
						thisEpochDicts.append(thisLineData);
			thisFileData.append(thisEpochDicts);
		#Next: Get the GNSS-wise order of observable [using getRinSysObs(GNSS, RinFile)], read these observables and put into a data-cube
		#Make a data structure for each line, put variables with observable names e.g. C1X, etc.
		#Tip: define a global dictionary data type with all possible observable names and from each line, add values to those keys for which data exists
	return thisFileData;
	
def getGPSWeek(yyyy,mm,dd,hh,min,sec):
	sec = str(int(sec));
	min = str(min);
	hh = str(hh);
	dd = str(dd);
	mm = str(mm);
	yyyy = str(yyyy);
	timestamp = datetime.datetime.strptime(mm+" "+dd+" "+yyyy+" "+hh+":"+min+":"+sec,"%m %d %Y %H:%M:%S");
	timestamp = timestamp.strftime("%m %d %Y %H:%M:%S");
	tc_out = commands.getoutput("timeconvert -c \""+timestamp+"\"");
	GPSWeekTS = tc_out.split('\n')[4].split()[2];
	return int(GPSWeekTS);
	DoW = tc_out.split('\n')[3].split()[4];
	SoWTS = tc_out.split('\n')[3].split()[5];
	DoY = tc_out.split('\n')[5].split()[4];

def getDOY(yyyy,mm,dd,hh,min,sec):
	sec = str(int(sec));
	min = str(min);
	hh = str(hh);
	dd = str(dd);
	mm = str(mm);
	yyyy = str(yyyy);
	timestamp = datetime.datetime.strptime(mm+" "+dd+" "+yyyy+" "+hh+":"+min+":"+sec,"%m %d %Y %H:%M:%S");
	timestamp = timestamp.strftime("%m %d %Y %H:%M:%S");
	tc_out = commands.getoutput("timeconvert -c \""+timestamp+"\"");
	DoY = tc_out.split('\n')[5].split()[4];
	return int(DoY);
		
def getDOW(yyyy,mm,dd,hh,min,sec):
	sec = str(int(sec));
	min = str(min);
	hh = str(hh);
	dd = str(dd);
	mm = str(mm);
	yyyy = str(yyyy);
	timestamp = datetime.datetime.strptime(mm+" "+dd+" "+yyyy+" "+hh+":"+min+":"+sec,"%m %d %Y %H:%M:%S");
	timestamp = timestamp.strftime("%m %d %Y %H:%M:%S");
	tc_out = commands.getoutput("timeconvert -c \""+timestamp+"\"");
	DoW = tc_out.split('\n')[3].split()[4];
	return int(DoW);
	
def getSOW(yyyy,mm,dd,hh,min,sec):
	sec = str(int(sec));
	min = str(min);
	hh = str(hh);
	dd = str(dd);
	mm = str(mm);
	yyyy = str(yyyy);
	timestamp = datetime.datetime.strptime(mm+" "+dd+" "+yyyy+" "+hh+":"+min+":"+sec,"%m %d %Y %H:%M:%S");
	timestamp = timestamp.strftime("%m %d %Y %H:%M:%S");
	tc_out = commands.getoutput("timeconvert -c \""+timestamp+"\"");
	SoWTS = float(tc_out.split('\n')[3].split()[5]);
	return int(SoWTS);

def epoch_NumSat_GNSS(epoch,FGNSS):
	numSats = 0;
	for line in epoch:
		if line['IGNSSF'] == FGNSS:
			numSats = numSats + 1;
	return numSats;
	
def writeBinData(InFile,RinData,FGNSS,OutFile):
	#binOut = open(OutFile,'wb');
	binOut = FortranFile(OutFile, 'w');

	IDSTAC = getRinMarkerNumber(InFile)[0:5];
	for epoch in RinData: #epoch is list of dictionaries
		Y4 = epoch[0]['YYYY'];
		MM = epoch[0]['MM'];
		DD = epoch[0]['DD'];
		hh = epoch[0]['hh'];
		mm = epoch[0]['mm'];
		ss = epoch[0]['ss'];
		RCB = epoch[0]['RCBIAS'];
		if RCB=='':
			RCB = 0;
		JWK = getGPSWeek(Y4,MM,DD,hh,mm,ss);
		WSEC = getSOW(Y4,MM,DD,hh,mm,ss);
		JTTAG = int(WSEC)*int(1000);#not rounding yet + float(0.5);
		ITPRN = 0;
		ITPHS = 0;
		TSCALE = 1e9;
		ITCUS = float(RCB)*float(TSCALE);
		ISIGMA = 0;
		NSAT = 0;
		MSAT = epoch_NumSat_GNSS(epoch,FGNSS);#epoch[0]['NumSat'];
		TCORB = False;
		epoch_JSV = [];
		epoch_PCL1 = [];
		epoch_PRL1 = [];
		epoch_PCL2 = [];
		epoch_PRL2 = [];
		epoch_PRL3 = [];
		epoch_PHL1 = [];
		epoch_PHL2 = [];
		epoch_PHL3 = [];
		epoch_NRL1 = [];
		epoch_NRL2 = [];
		epoch_NRL3 = [];
		epoch_IEDPH1 = [];
		epoch_IEDPH2 = [];
		epoch_IEDPH3 = [];
		epoch_ICYFX1 = [];
		epoch_ICYFX2 = [];
		epoch_ICYFX3 = [];
		epoch_RORB = [];
		epoch_IAZ = [];
		epoch_IEL = [];
		epoch_JPASS = [];
		for i in range(0,int(MSAT)):
			epoch_ICYFX1.append(False); #Default value for now (False);
			epoch_ICYFX2.append(False); #Default value for now (False);
			epoch_ICYFX3.append(False); #Default value for now (False);
			epoch_RORB.append(float(24000000)); #Default value for now (24000km)
			epoch_IAZ.append(0); #Default value for now (0)
			epoch_IEL.append(0); #Default value for now (0)
		
		#print "new epoch starts...";
		for line in epoch:
			if line['IGNSSF'] == FGNSS:epoch_JSV.append(line['SatNum']); #Sat. PRN's
			
			#Initialize observables for this line (i.e. this satellite in current epoch)
			PCL1 = 0;
			PRL1 = 0;
			PHL1 = 0;
			NRL1 = 0;
			IEDPH1 = 0;
			PCL2 = 0;
			PRL2 = 0;
			PHL2 = 0;
			NRL2 = 0;
			IEDPH2 = 0;
			PRL3 = 0;
			PHL3 = 0;
			NRL3 = 0;
			IEDPH3 = 0;
			DT = 0;
			
			#Fill in the GPS Observables:
			if (FGNSS == 'G'):
				if line['IGNSSF'] == FGNSS: 
				#For PCL1, assumes that only ONE of the C1C, C1S, C1L or C1X exists in observations
					if 'C1C' in line.keys():
						PCL1 = line['C1C'];
					elif 'C1S' in line.keys():
						PCL1 = line['C1S'];
					elif 'C1L' in line.keys():
						PCL1 = line['C1L'];
					elif 'C1X' in line.keys():
						PCL1 = line['C1X'];
				#For PRL1, assumes that only ONE of the C1P, C1W, C1Y or C1M exists in observations	
					if 'C1P' in line.keys():
						PRL1 = line['C1P'];
					elif 'C1W' in line.keys():
						PRL1 = line['C1W'];
					elif 'C1Y' in line.keys():
						PRL1 = line['C1Y'];
					elif 'C1M' in line.keys():
						PRL1 = line['C1M'];
				#For PHL1, assumes that only ONE of the L1C, L1P, L1W, L1Y or L1M exists in observations
					if 'L1C' in line.keys():
						if(PHL1==0):
							PHL1 = line['L1C'];
							NRL1 = line['S1C'];
					elif 'L1P' in line.keys():
						if(PHL1==0):
							PHL1 = line['L1P'];
							NRL1 = line['S1P'];
					elif 'L1W' in line.keys():
						if(PHL1==0):
							PHL1 = line['L1W'];
							NRL1 = line['S1W'];
					elif 'L1Y' in line.keys():
						if(PHL1==0):
							PHL1 = line['L1Y'];
							NRL1 = line['S1Y'];
					elif 'L1M' in line.keys():
						if(PHL1==0):
							PHL1 = line['L1M'];
							NRL1 = line['S1M'];
				#For PCL2, assumes that only ONE of the C2C, C2S, C2L or C2X exists in observations
					if 'C2C' in line.keys():
						PCL2 = line['C2C'];
					elif 'C2S' in line.keys():
						PCL2 = line['C2S'];
					elif 'C2L' in line.keys():
						PCL2 = line['C2L'];
					elif 'C2X' in line.keys():
						PCL2 = line['C2X'];
				#For PRL2, assumes that only ONE of the C2P, C2W, C2Y or C2M exists in observations	
					if 'C2P' in line.keys():
						PRL2 = line['C2P'];
					elif 'C2W' in line.keys():
						PRL2 = line['C2W'];
					elif 'C2Y' in line.keys():
						PRL2 = line['C2Y'];
					elif 'C2M' in line.keys():
						PRL2 = line['C2M'];
				#For PHL2, assumes that only ONE of the L2C, L2P, L2W, L2Y or L2M exists in observations
					if 'L2C' in line.keys():
						if(PHL2==0):
							PHL2 = line['L2C'];
							NRL2 = line['S2C'];
					elif 'L2P' in line.keys():
						if(PHL2==0):
							PHL2 = line['L2P'];
							NRL2 = line['S2P'];
					elif 'L2W' in line.keys():
						if(PHL2==0):
							PHL2 = line['L2W'];
							NRL2 = line['S2W'];
					elif 'L2Y' in line.keys():
						if(PHL2==0):
							PHL2 = line['L2Y'];
							NRL2 = line['S2Y'];
					elif 'L2M' in line.keys():
						if(PHL2==0):
							PHL2 = line['L2M'];
							NRL2 = line['S2M'];
				#For PRL3, assumes that only ONE of the C5Q, C5X or C5I exists in observations	
					if 'C5Q' in line.keys():
						PRL3 = line['C5Q'];
					elif 'C5X' in line.keys():
						PRL3 = line['C5X'];
					elif 'C5I' in line.keys():
						PRL3 = line['C5I'];
				#For PHL3, assumes that only ONE of the L5Q, L5X or L5I exists in observations	
					if 'L5Q' in line.keys():
						PHL3 = line['L5Q'];
						NRL3 = line['S5Q'];
					elif 'L5X' in line.keys():
						PHL3 = line['L5X'];
						NRL3 = line['S5X'];
					elif 'L5I' in line.keys():
						PHL3 = line['L5I'];
						NRL3 = line['S5I'];
					if PHL1 == 0:IEDPH1=1;
					if PHL2 == 0:IEDPH2=1;
					if PHL3 == 0:IEDPH3=1;
					if PRL1 == 0:PRL1 = PCL1;
					if PRL2 == 0:PRL2 = PCL1;
					SOLR = float(1)/float(299792458.00);
					DT = float(PRL1) * float(SOLR);
					epoch_PCL1.append(PCL1);
					epoch_PRL1.append(PRL1);
					epoch_PHL1.append(PHL1);
					epoch_NRL1.append(NRL1);
					epoch_IEDPH1.append(IEDPH1);
					epoch_IEDPH2.append(IEDPH2);
					epoch_IEDPH3.append(IEDPH3);
					epoch_PCL2.append(PCL2);
					epoch_PRL2.append(PRL2);
					epoch_PHL2.append(PHL2);
					epoch_NRL2.append(NRL2);
					epoch_PRL3.append(PRL3);
					epoch_PHL3.append(PHL3);		
					epoch_NRL3.append(NRL3);
			
			#Fill in the GALILEO Observables:
			if (FGNSS == 'E'):
				if line['IGNSSF'] == FGNSS: 
				#For PCL1, assumes that only ONE of the C1A, C1B, C1C, C1X or C1Z exists in observations
					if 'C1A' in line.keys():
						PCL1 = line['C1A'];
					elif 'C1B' in line.keys():
						PCL1 = line['C1B'];
					elif 'C1C' in line.keys():
						PCL1 = line['C1C'];
					elif 'C1X' in line.keys():
						PCL1 = line['C1X'];
					elif 'C1Z' in line.keys():
						PCL1 = line['C1Z'];
				#For PRL1, assumes same as PCL1
					if 'C1A' in line.keys():
						PRL1 = line['C1A'];
					elif 'C1B' in line.keys():
						PRL1 = line['C1B'];
					elif 'C1C' in line.keys():
						PRL1 = line['C1C'];
					elif 'C1X' in line.keys():
						PRL1 = line['C1X'];
					elif 'C1Z' in line.keys():
						PRL1 = line['C1Z'];
				#For PHL1, assumes that only ONE of the L1A, L1B, L1C, L1X or L1Z exists in observations
					if 'L1A' in line.keys():
						if(PHL1==0):
							PHL1 = line['L1A'];
							NRL1 = line['S1A'];
					elif 'L1B' in line.keys():
						if(PHL1==0):
							PHL1 = line['L1B'];
							NRL1 = line['S1B'];
					elif 'L1C' in line.keys():
						if(PHL1==0):
							PHL1 = line['L1C'];
							NRL1 = line['S1C'];
					elif 'L1X' in line.keys():
						if(PHL1==0):
							PHL1 = line['L1X'];
							NRL1 = line['S1X'];
					elif 'L1Z' in line.keys():
						if(PHL1==0):
							PHL1 = line['L1Z'];
							NRL1 = line['S1Z'];
				#For PCL2, assumes that only ONE of the C7I, C7Q or C7X exists in observations
					if 'C7I' in line.keys():
						PCL2 = line['C7I'];
					elif 'C7Q' in line.keys():
						PCL2 = line['C7Q'];
					elif 'C7X' in line.keys():
						PCL2 = line['C7X'];
				#For PRL2, assumes same as PCL2
					if 'C7I' in line.keys():
						PRL2 = line['C7I'];
					elif 'C7Q' in line.keys():
						PRL2 = line['C7Q'];
					elif 'C7X' in line.keys():
						PRL2 = line['C7X'];
				#For PHL2, assumes that only ONE of the L7I, L7Q or L7X exists in observations
					if 'L7I' in line.keys():
						if(PHL2==0):
							PHL2 = line['L7I'];
							NRL2 = line['S7I'];
					elif 'L7Q' in line.keys():
						if(PHL2==0):
							PHL2 = line['L7Q'];
							NRL2 = line['S7Q'];
					elif 'L7X' in line.keys():
						if(PHL2==0):
							PHL2 = line['L7X'];
							NRL2 = line['S7X'];
				#For PRL3, assumes that only ONE of the C5Q, C5X or C5I exists in observations	
					if 'C5Q' in line.keys():
						PRL3 = line['C5Q'];
					elif 'C5X' in line.keys():
						PRL3 = line['C5X'];
					elif 'C5I' in line.keys():
						PRL3 = line['C5I'];
				#For PHL3, assumes that only ONE of the L5Q, L5X or L5I exists in observations	
					if 'L5Q' in line.keys():
						PHL3 = line['L5Q'];
						NRL3 = line['S5Q'];
					elif 'L5X' in line.keys():
						PHL3 = line['L5X'];
						NRL3 = line['S5X'];
					elif 'L5I' in line.keys():
						PHL3 = line['L5I'];
						NRL3 = line['S5I'];
					if PHL1 == 0:IEDPH1=1;
					if PHL2 == 0:IEDPH2=1;
					if PHL3 == 0:IEDPH3=1;
					if PRL1 == 0:PRL1 = PCL1;
					if PRL2 == 0:PRL2 = PCL1;
					SOLR = float(1)/float(299792458.00);
					DT = float(PRL1) * float(SOLR);
					epoch_PCL1.append(PCL1);
					epoch_PRL1.append(PRL1);
					epoch_PHL1.append(PHL1);
					epoch_NRL1.append(NRL1);
					epoch_IEDPH1.append(IEDPH1);
					epoch_IEDPH2.append(IEDPH2);
					epoch_IEDPH3.append(IEDPH3);
					epoch_PCL2.append(PCL2);
					epoch_PRL2.append(PRL2);
					epoch_PHL2.append(PHL2);
					epoch_NRL2.append(NRL2);
					epoch_PRL3.append(PRL3);
					epoch_PHL3.append(PHL3);		
					epoch_NRL3.append(NRL3);
			
			
			#Fill in the BeiDou Observables:
			if (FGNSS == 'C'):
				if line['IGNSSF'] == FGNSS: 
				#For PCL1, assumes that only ONE of the C2I, C2Q or C2X exists in observations
					if 'C2I' in line.keys():
						PCL1 = line['C2I'];
					elif 'C2Q' in line.keys():
						PCL1 = line['C2Q'];
					elif 'C2X' in line.keys():
						PCL1 = line['C2X'];
				#For PRL1, assumes same as PCL1
					if 'C2I' in line.keys():
						PRL1 = line['C2I'];
					elif 'C2Q' in line.keys():
						PRL1 = line['C2Q'];
					elif 'C2X' in line.keys():
						PRL1 = line['C2X'];
				#For PHL1, assumes that only ONE of the L21, L2Q or L2X exists in observations
					if 'L2I' in line.keys():
						if(PHL1==0):
							PHL1 = line['L2I'];
							NRL1 = line['S2I'];
					elif 'L2Q' in line.keys():
						if(PHL1==0):
							PHL1 = line['L2Q'];
							NRL1 = line['S2Q'];
					elif 'L2X' in line.keys():
						if(PHL1==0):
							PHL1 = line['L2X'];
							NRL1 = line['S2X'];
				#For PCL2, assumes that only ONE of the C7I, C7Q or C7X exists in observations
					if 'C7I' in line.keys():
						PCL2 = line['C7I'];
					elif 'C7Q' in line.keys():
						PCL2 = line['C7Q'];
					elif 'C7X' in line.keys():
						PCL2 = line['C7X'];
				#For PRL2, assumes same as PCL2
					if 'C7I' in line.keys():
						PRL2 = line['C7I'];
					elif 'C7Q' in line.keys():
						PRL2 = line['C7Q'];
					elif 'C7X' in line.keys():
						PRL2 = line['C7X'];
				#For PHL2, assumes that only ONE of the L71, L7Q or L7X exists in observations
					if 'L7I' in line.keys():
						if(PHL2==0):
							PHL2 = line['L7I'];
							NRL2 = line['S7I'];
					elif 'L7Q' in line.keys():
						if(PHL2==0):
							PHL2 = line['L7Q'];
							NRL2 = line['S7Q'];
					elif 'L7X' in line.keys():
						if(PHL2==0):
							PHL2 = line['L7X'];
							NRL2 = line['S7X'];
				#For PRL3, assumes that only ONE of the C6I, C6Q or C6X exists in observations
					if 'C6I' in line.keys():
						PRL3 = line['C6I'];
					elif 'C6Q' in line.keys():
						PRL3 = line['C6Q'];
					elif 'C6X' in line.keys():
						PRL3 = line['C6X'];
				#For PHL3, assumes that only ONE of the L61, L6Q or L6X exists in observations
					if 'L6I' in line.keys():
						if(PHL3==0):
							PHL3 = line['L6I'];
							NRL3 = line['S6I'];
					elif 'L6Q' in line.keys():
						if(PHL3==0):
							PHL3 = line['L6Q'];
							NRL3 = line['S6Q'];
					elif 'L6X' in line.keys():
						if(PHL3==0):
							PHL3 = line['L6X'];
							NRL3 = line['S6X'];
					if PHL1 == 0:IEDPH1=1;
					if PHL2 == 0:IEDPH2=1;
					if PHL3 == 0:IEDPH3=1;
					if PRL1 == 0:PRL1 = PCL1;
					if PRL2 == 0:PRL2 = PCL1;
					SOLR = float(1)/float(299792458.00);
					DT = float(PRL1) * float(SOLR);
					epoch_PCL1.append(PCL1);
					epoch_PRL1.append(PRL1);
					epoch_PHL1.append(PHL1);
					epoch_NRL1.append(NRL1);
					epoch_IEDPH1.append(IEDPH1);
					epoch_PCL2.append(PCL2);
					epoch_PRL2.append(PRL2);
					epoch_PHL2.append(PHL2);
					epoch_NRL2.append(NRL2);
					epoch_IEDPH2.append(IEDPH2);
					epoch_PRL3.append(PRL3);
					epoch_PHL3.append(PHL3);		
					epoch_NRL3.append(NRL3);			
					epoch_IEDPH3.append(IEDPH3);
			
			epoch_JPASS.append(0); #Default value for now (0)
			
			#PCL1 = 
			
		if(MSAT>0):
			#print IDSTAC,JWK,JTTAG,ITPRN,ITPHS,ITCUS,ISIGMA,NSAT,MSAT,TCORB,epoch_JSV,epoch_PCL1,epoch_PCL1,epoch_PRL1,epoch_PRL2,epoch_PRL3,epoch_PHL1,epoch_PHL2,epoch_PHL3,epoch_NRL1,epoch_NRL2,epoch_NRL3,epoch_IEDPH3,epoch_IEDPH2,epoch_IEDPH3,epoch_ICYFX1,epoch_ICYFX2,epoch_ICYFX3,epoch_RORB,epoch_IAZ,epoch_IEL;
			print "Writing epoch",JWK,JTTAG;
			binOut.write_record(numpy.int32(IDSTAC));
			binOut.write_record(numpy.int32(JWK));
			binOut.write_record(numpy.int32(JTTAG));
			binOut.write_record(numpy.int32(ITPRN));
			binOut.write_record(numpy.int32(ITPHS));
			binOut.write_record(numpy.int32(ITCUS));
			binOut.write_record(numpy.int32(ISIGMA));
			binOut.write_record(numpy.int32(MSAT));
			binOut.write_record(numpy.bool_(TCORB));
			for item in epoch_JSV:
				binOut.write_record(numpy.int32(item));
			for item in epoch_PCL1:
				binOut.write_record(numpy.float64(item));
			for item in epoch_PRL1:
				binOut.write_record(numpy.float64(item));
			for item in epoch_PRL2:
				binOut.write_record(numpy.float64(item));
			for item in epoch_PRL3:
				binOut.write_record(numpy.float64(item));
			for item in epoch_PHL1:
				binOut.write_record(numpy.float64(item));
			for item in epoch_PHL2:
				binOut.write_record(numpy.float64(item));
			for item in epoch_PHL3:
				binOut.write_record(numpy.float64(item));
			for item in epoch_NRL1:
				binOut.write_record(numpy.int32(item));
			for item in epoch_NRL2:
				binOut.write_record(numpy.int32(item));
			for item in epoch_NRL3:
				binOut.write_record(numpy.int32(item));
			for item in epoch_IEDPH1:
				binOut.write_record(numpy.int32(item));
			for item in epoch_IEDPH2:
				binOut.write_record(numpy.int32(item));
			for item in epoch_IEDPH3:
				binOut.write_record(numpy.int32(item));
			for item in epoch_ICYFX1:
				binOut.write_record(numpy.bool_(item));
			for item in epoch_ICYFX2:
				binOut.write_record(numpy.bool_(item));
			for item in epoch_ICYFX3:
				binOut.write_record(numpy.bool_(item));
			for item in epoch_RORB:
				binOut.write_record(numpy.float64(item));
			for item in epoch_IAZ:
				binOut.write_record(numpy.int32(item));
			for item in epoch_IEL:
				binOut.write_record(numpy.int32(item));
			for item in epoch_JPASS:
				binOut.write_record(numpy.int32(item));
			
		#print "new epoch ends\n"
		
	# JTTAG
	# ITPRN
	# ITPHS
	# ITCUS
	# ISIGMA
	# NSAT
	# MSAT
	# TCORB
	
	binOut.close();
	return 0;

#Testing / execution below this line:

systems = getSysToRead();

# print getRinVer('SAMPLE.rnx');
# print getRinFileType('SAMPLE.rnx');
# print getRinSysType('SAMPLE.rnx');
# print getRinComments('SAMPLE.rnx');
# print getRcvInfo('SAMPLE.rnx');
# print getRinXYZ('SAMPLE.rnx');
# print getRinEcc('SAMPLE.rnx');
# print getRinSysNumTypes('SAMPLE.rnx');
# print getRinSysObs('G','SAMPLE.rnx');
# print getRinSysObs('R','SAMPLE.rnx');
# print getRinSysObs('E','SAMPLE.rnx');
# print getRinSysObs('S','SAMPLE.rnx');
# print getNumOfEpochs('SAMPLE.rnx');
########

InFile = sys.argv[1]+".rnx";
OutFile = sys.argv[1]+".bin";
FGNSS = 'G';

RinData = readRinData(InFile);
#writeBinData(InFile, RinData, FGNSS, OutFile);
writeBinData(InFile, RinData, 'G', OutFile);


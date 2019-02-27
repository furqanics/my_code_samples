#!/usr/bin/python

import os
import sys
import numpy
from datetime import datetime, timedelta
import time
import commands
import math
import jdcal
from scipy.integrate import ode

#from MyFunctions import *



#satPRN, M0, rootA, deltaN, ecc, omega0, Cuc, Cus, Crc, Crs, Cic, Cis, i0, idot, OMEGA0, OMEGADOT, ToE
def getGEPos(Tk,af0,af1,af2,sysType,satPRN,M0,rootA,deltaN,ecc,omega0,Cuc,Cus,Crc,Crs,Cic,Cis,i0,idot,OMEGA0,OMEGADOT,ToE):
	#Get positions at Tk for G and E satellites
	if(float(Tk)>float(302400)):
		Tk = float(Tk)-float(604800);
	elif(float(Tk)<float(-302400)):
		Tk = float(Tk)+float(604800);
	if debug==1:showVar(Tk);
	
	svClkBias = float(af0) + float(af1)*Tk + float(af2)*pow(Tk,2);
	if debug==1:showVar(svClkBias);
			
	#Set value of GE [m^3/s^2]
	GM = float(3.986004418e14); #For GPS
		
	if(sysType=='E'):GM = float(3.986004418e14); #For Galileo
	
	if debug==1:showVar(GM);
	
	#Get Mean anomaly at SoWTS
	A = float(rootA)*float(rootA);
	A3 = A*A*A;
	#A3 = float(pow(float(rootA),6))
	
	if debug==1:showVar(A3);
	
	M = float(M0) + ((float(numpy.sqrt(GM/A3))+float(deltaN))*Tk);
	if debug==1:showVar(M);
	#Iterative solution of Eccentric Anomaly (E) (we do 20 iterations here)
	ecc=float(ecc);
	E0 = M;
	dE = float(1);
	while(dE > float(1e-11)):
		E = M+ecc*numpy.sin(E0)
		dE = abs(float(E)-float(E0))
		E0 = E;
		
	omegae = float(7.2921151467e-5) #this value is for GPS and Galileo systems
	if debug==1:showVar(E);
	#True anomaly (v):
	v = math.atan2(numpy.sqrt(1-pow(ecc,2))*numpy.sin(E),numpy.cos(E)-ecc)
	if debug==1:showVar(v);
	
	#Corrections for orbital perturbations:
	
	OMEGA = float(OMEGA0) + ((float(OMEGADOT)-float(omegae))*Tk) - (omegae * float(ToE))
	if debug==1:showVar(OMEGA);
	
	
	if debug==1:showVar(omegae);
	#Argument of perigee:
	phi = float(omega0)+float(v)
	phi = float(phi)%float(2*3.1416);
	omega = float(omega0) + float(Cuc)*numpy.cos(2*phi) + float(Cus)*numpy.sin(2*phi);
	if debug==1:showVar(omega);
	
	#Radial distance:
	r = A*(1-ecc*numpy.cos(E)) + float(Crc)*numpy.cos(2*phi) + float(Crs)*numpy.sin(2*phi);
	
	#Inclination:
	#if(float(Cis)<0):Cis=float(Cis)*(-1)
	if debug==1:showVar(Cis);
	i = float(i0) + float(idot)*Tk + float(Cic)*numpy.cos(2*phi) + float(Cis)*numpy.sin(2*phi);
	if debug==1:showVar(i);
	
	
	if debug==1:showVar(r);
	
	#Argument of latitude:
	uk = float(omega)+float(v)+(float(Cuc)*numpy.cos(2*(float(omega)+float(v))))+(float(Cus)*numpy.sin(2*(float(omega)+float(v))));
	if debug==1:showVar(uk);
	
	#Computation of the right ascension, accounting for Earth's rotation:
	
	
	
	#Rotation matrices
	R3a = numpy.zeros([3,3],float);
	R1 = numpy.zeros([3,3],float);
	R3b = numpy.zeros([3,3],float);
	RV = numpy.zeros([3,1],float);
	
	theta = OMEGA*(-1)
	R3a[0,0] = numpy.cos(theta);
	R3a[0,1] = numpy.sin(theta);
	R3a[0,2] = 0;
	R3a[1,0] = (-1)*numpy.sin(theta);
	R3a[1,1] = numpy.cos(theta);
	R3a[1,2] = 0;
	R3a[2,0] = 0;
	R3a[2,1] = 0;
	R3a[2,2] = 1;
	if debug==1:showVar(R3a);
	
	theta = uk*(-1)
	R3b[0,0] = numpy.cos(theta);
	R3b[0,1] = numpy.sin(theta);
	R3b[0,2] = 0;
	R3b[1,0] = (-1)*numpy.sin(theta);
	R3b[1,1] = numpy.cos(theta);
	R3b[1,2] = 0;
	R3b[2,0] = 0;
	R3b[2,1] = 0;
	R3b[2,2] = 1;
	if debug==1:showVar(R3b);
	
	theta = i*(-1)
	R1[0,0] = 1;
	R1[0,1] = 0;
	R1[0,2] = 0;
	R1[1,0] = 0;
	R1[1,1] = numpy.cos(theta);
	R1[1,2] = numpy.sin(theta);
	R1[2,0] = 0;
	R1[2,1] = (-1)*numpy.sin(theta);
	R1[2,2] = numpy.cos(theta);
	if debug==1:showVar(R1);
	
	RV[0,0] = r;
	if debug==1:showVar(r);
	
	temp1 = numpy.dot(R3a,R1)
	if debug==1:showVar(temp1);
	
	temp2 = numpy.dot(R3b,RV)
	if debug==1:showVar(temp2);
	
	trs_pos = numpy.dot(temp1,temp2)/1000;
	if debug==1:showVar(trs_pos);
	
	#Position of the satellite in orbital frame (r_vec):
	r_vec = numpy.zeros([3,1],float);
	r_vec[0,0]=float(r)*float(numpy.cos(v))
	r_vec[1,0]=float(r)*float(numpy.sin(v))
	if debug==1:showVar(r_vec);
	#Rotation matrix (R_ecef) to translate to ECEF position:
	R_ecef = numpy.zeros([3,3],float);
	R_ecef[0,0] = (numpy.cos(OMEGA)*numpy.cos(omega)) - (numpy.sin(OMEGA)*numpy.sin(omega)*numpy.cos(i))
	R_ecef[0,1] = -(numpy.cos(OMEGA)*numpy.sin(omega)) - (numpy.sin(OMEGA)*numpy.cos(omega)*numpy.cos(i))
	R_ecef[0,2] = numpy.sin(OMEGA)*numpy.sin(i)
	R_ecef[1,0] = (numpy.sin(OMEGA)*numpy.cos(omega)) + (numpy.cos(OMEGA)*numpy.sin(omega)*numpy.cos(i))
	R_ecef[1,1] = -(numpy.sin(OMEGA)*numpy.sin(omega)) + (numpy.cos(OMEGA)*numpy.cos(omega)*numpy.cos(i))
	R_ecef[1,2] = -numpy.cos(OMEGA)*numpy.sin(i)
	R_ecef[2,0] = numpy.sin(omega)*numpy.sin(i)
	R_ecef[2,1] = numpy.cos(omega)*numpy.sin(i)
	R_ecef[2,2] = numpy.cos(i)
	if debug==1:showVar(R_ecef);
	#Apply rotation:
	sat_pos_ecef = numpy.dot(R_ecef,r_vec); #X,Y,Z in meters
	#sat_pos_ecef = sat_pos_ecef/float(1000); #X, Y, Z in Kilometers

    #Get GPS/Gal. Satellite Velocities using program "velo.exe roota toe m0 e delta_n smallomega cus cuc crs crc cis cic idot i0 bigomega0 bigomegadot t":
#	print "getting velo for this epoch...";
#	print "./velo.exe "+str(rootA)+" "+str(ToE)+" "+str(M0)+" "+str(ecc)+" "+str(deltaN)+" "+str(omega0)+" "+str(Cus)+" "+str(Cuc)+" "+str(Crs)+" "+str(Crc)+" "+str(Cis)+" "+str(Cic)+" "+str(idot)+" "+str(i0)+" "+str(OMEGA0)+" "+str(OMEGADOT)+" "+str(SoWTS)+" | grep \"BCvel\"";
	vel_temp=commands.getoutput("./velo.exe "+str(rootA)+" "+str(ToE)+" "+str(M0)+" "+str(ecc)+" "+str(deltaN)+" "+str(omega0)+" "+str(Cus)+" "+str(Cuc)+" "+str(Crs)+" "+str(Crc)+" "+str(Cis)+" "+str(Cic)+" "+str(idot)+" "+str(i0)+" "+str(OMEGA0)+" "+str(OMEGADOT)+" "+str(SoWTS)+" | grep \"BCvel\"").split();
#	print vel_temp;
	Vx = float(vel_temp[6]);
	Vy = float(vel_temp[7]);
	Vz = float(vel_temp[8]);
#	print "Vx, Vy, Vz:",Vx,Vy,Vz;
	R,T,N = XYZ2RTN([sat_pos_ecef[0]*1000,sat_pos_ecef[1]*1000,sat_pos_ecef[2]*1000],[Vx,Vy,Vz]); # [m] and [m/s]
#	print "RTN: ",R,T,N;
	
	
	#sat_pos_ecef = trs_pos; #Update to navipedia method
	
	if debug==1:showVar(sat_pos_ecef);	
	return sat_pos_ecef,Vx,Vy,Vz,R,T,N,svClkBias;


def XYZ2RTN(P,V):
	#XYZ2RTN: Convert the position vector P from XYZ to RTN
	#Input vector V is velocity vector
	X = float(P[0]);
	Y = float(P[1]);
	Z = float(P[2]);
	Vx = float(V[0]);
	Vy = float(V[1]);
	Vz = float(V[2]);
	R = 0;
	T = 0;
	N = 0;
	RTNrot = numpy.zeros((3,3));
	temp = 1 / (numpy.sqrt(numpy.square(X)+numpy.square(Y)+numpy.square(Z)));
	
	RTNrot[0,0] = X * temp;
	RTNrot[0,1] = Y * temp;
	RTNrot[0,2] = Z * temp;
	
	RTNrot[2,0] = (Y * Vz) - (Z * Vy);
	RTNrot[2,1] = (Z * Vx) - (X * Vz);
	RTNrot[2,2] = (X * Vy) - (Y * Vx);
	temp = 1 / (numpy.sqrt(numpy.square(RTNrot[2,0])+numpy.square(RTNrot[2,1])+numpy.square(RTNrot[2,2])));
	RTNrot[2,0] = RTNrot[2,0]*temp;
	RTNrot[2,1] = RTNrot[2,1]*temp;
	RTNrot[2,2] = RTNrot[2,2]*temp;
	
	RTNrot[1,0] = RTNrot[2,1]*RTNrot[0,2] - RTNrot[2,2]*RTNrot[0,1];
	RTNrot[1,1] = RTNrot[2,2]*RTNrot[0,0] - RTNrot[2,0]*RTNrot[0,2];
	RTNrot[1,2] = RTNrot[2,0]*RTNrot[0,1] - RTNrot[2,1]*RTNrot[0,0];

	XYZ = numpy.zeros((3,1));
	XYZ[0,0] = X;
	XYZ[1,0] = Y;
	XYZ[2,0] = Z;
	
	RTN = numpy.dot(RTNrot,XYZ);
	R = RTN[0,0];
	T = RTN[1,0];
	N = RTN[2,0];
	return R,T,N;

def glonassfunc(t,y):
	u = float(398600.44); #km^3/s^2
	ae= float(6378.136);
	C20= float(-1082.63e-6);
	w3= float(0.7292115e-4);
	xa = y[0];
	xdot = y[1];
	ya = y[2];
	ydot = y[3];
	za = y[4];
	zdot = y[5];
	r = numpy.sqrt(xa*xa+ya*ya+za*za);
	ubar = u/pow(r,2);
	xdot2 = -ubar*(xa/r) + (3/2)*C20*ubar*(xa/r)*pow(ae/r,2)*(1-5*(pow(za/r,2))) + ddx_te*numpy.cos(theta_Ge) - ddy_te*numpy.sin(theta_Ge);
	ydot2 = -ubar*(ya/r) + (3/2)*C20*ubar*(ya/r)*pow(ae/r,2)*(1-5*(pow(za/r,2))) + ddx_te*numpy.sin(theta_Ge) + ddy_te*numpy.cos(theta_Ge);
	zdot2 = -ubar*(za/r) + (3/2)*C20*ubar*(za/r)*pow(ae/r,2)*(3-5*(pow(za/r,2))) + ddz_te;
	xyz_pv = numpy.zeros((6,1)); #Column vector
	xyz_pv[0] = xdot;
	xyz_pv[1] = xdot2;
	xyz_pv[2] = ydot;
	xyz_pv[3] = ydot2;
	xyz_pv[4] = zdot;
	xyz_pv[5] = zdot2;
	return xyz_pv;

			#xdot_out[k] = integ_gln.y[0]
			#ydot_out[k] = integ_gln.y[1]
			#zdot_out[k] = integ_gln.y[2]
			#xdot2_out[k] = integ_gln.y[3]
			#ydot2_out[k] = integ_gln.y[4]
			#zdot2_out[k] = integ_gln.y[5]
	
	
def namestr(obj, namespace):
	return [name for name in namespace if namespace[name] is obj];

def showVar(var):
	varName = namestr(var, globals());
	print str(varName[0])+" = "+str(var);
	return;

print "Script start time: ",os.system("date");
	
debug = 0; #1 to display variable values using showVar();

inp_brdc = sys.argv[1]

print inp_brdc;

brdc_path = "/home/fahmed/JTP/PROD/EPH/"

if(os.path.exists(brdc_path+inp_brdc)):
	brdc = open(brdc_path+inp_brdc,'r').readlines();
elif(not os.path.exists(brdc_path+inp_brdc)):
	sys.exit("File "+inp_brdc+" not found. Please retry with an existing file.")

#print len(brdc)

startIndex = 0;
lineCount = 8; #variable to identify number of lines in one record of each GNSS type
satsRead = [];

EOheaderFound = 0;

LeapSecNow = 18;

for line in brdc:
	if "LEAP SECONDS" in line:
		LeapSecNow = int(line.split()[0])
		break;

for line in brdc:
	if "END OF HEADER" in line:
		EOheaderFound=1;
		startIndex = brdc.index(line)+1;
		break;
if (EOheaderFound==0):
	sys.exit("Header does not end in file "+inp_brdc+".")

#Open output file:
eceffile = open(brdc_path+inp_brdc.split('.')[0]+inp_brdc.split('.')[1]+".ecef",'w');
#print startIndex;

brdc = brdc[startIndex:];
#print len(brdc)

pi = float(3.14159);
		
for line in brdc:
	sysType = '';
	if (line[0] == 'G' or line[0] == 'E'):# or line[0] == 'J'):# or line[0] == 'R' or  or line[0] == 'C' or  or line[0] == 'I' or line[0] == 'S'):
		thisIndex = brdc.index(line);
		sysType = line[0];
		satsRead.append(line[0:3])
		#print brdc[brdc.index(line)].split();
       
	   # if (sysType == 'G' or sysType == 'C' or sysType == 'E' or sysType == 'I' or sysType == 'J'):
			# lineCount=8;
        # elif (sysType == 'R' or sysType == 'S'):
			# lineCount=4;
			
		satPRN = line[0:3];
		
		yyyy = line[4:8] #GPS time in case of GPS and QZSS, UTC in case of GLONASS, GAL time in case of GALILEO
		mm = line[9:11]
		dd = line[12:14]
		hour = line[15:17]
		min = line[18:20]
		sec = line[21:23]
		sysWeek = brdc[thisIndex+5][42:61];
		M0 = brdc[thisIndex+1][61:80]; #[rad]
		rootA = brdc[thisIndex+2][61:80]; #[sqrt(m)]
		deltaN = brdc[thisIndex+1][42:61]; #[rad/s]
		ecc = brdc[thisIndex+2][23:42]; #[unitless]
		omega0 = brdc[thisIndex+4][42:61]; #[rad]
		Cuc = brdc[thisIndex+2][4:23]; #[rad]
		Cus = brdc[thisIndex+2][42:61]; #[rad]
		Crc = brdc[thisIndex+4][23:42]; #[m]
		Crs = brdc[thisIndex+1][23:42]; #[m]
		Cic = brdc[thisIndex+3][23:42]; #[rad]
		Cis = brdc[thisIndex+3][61:80]; #[rad]
		i0 = brdc[thisIndex+4][4:23]; #[rad]
		idot = brdc[thisIndex+5][4:23]; #[rad/s]
		OMEGA0 = brdc[thisIndex+3][42:61]; #[rad]
		OMEGADOT = brdc[thisIndex+4][61:80]; #[rad/s]
		ToE = float(brdc[thisIndex+3][4:23]); #[sec of GPS/Sys Week]
		ToT = float(brdc[thisIndex+7][4:23]); #[sec of GPS/Sys Week] : Transmission time
		timestamp = datetime.strptime(mm+" "+dd+" "+yyyy+" "+hour+":"+min+":"+sec,"%m %d %Y %H:%M:%S");
		timestamp = timestamp.strftime("%m %d %Y %H:%M:%S");
		tc_out = commands.getoutput("timeconvert -c \""+timestamp+"\"")
		#print tc_out;
		GPSWeekTS = tc_out.split('\n')[4].split()[2]
		DoW = tc_out.split('\n')[3].split()[4]
		#print tc_out.split('\n')[2].split()
		SoWTS = tc_out.split('\n')[3].split()[5]
		DoY = tc_out.split('\n')[5].split()[4]
		#Clock info:
		af0 = brdc[thisIndex+0][23:42]; #[sec]
		af1 = brdc[thisIndex+0][42:61]; #[sec/sec]
		af2 = brdc[thisIndex+0][61:80]; #[sec/sec^2]
		
		if debug==1:
			showVar(satPRN);
			showVar(M0);
			showVar(rootA);
			showVar(deltaN);
			showVar(ecc);
			showVar(omega0);
			showVar(Cuc);
			showVar(Cus);
			showVar(Crc);
			showVar(Crs);
			showVar(Cic);
			showVar(Cis);
			showVar(i0);
			showVar(idot);
			showVar(OMEGA0);
			showVar(OMEGADOT);
			showVar(ToE);
		
		
		#Now, compute the satellite positions based on the data extracted above:
		
		Tk = float(SoWTS)-float(ToE); #Time elapsed since ToE [sec]
        
		sat_pos_ecef,Vx,Vy,Vz,R,T,N,svClkBias=getGEPos(Tk,af0,af1,af2,sysType,satPRN,M0,rootA,deltaN,ecc,omega0,Cuc,Cus,Crc,Crs,Cic,Cis,i0,idot,OMEGA0,OMEGADOT,ToE);
		TkP = Tk+float(0.005);
		TkM = Tk-float(0.005);
		sat_pos_ecefP,VxP,VyP,VzP,RP,TP,NP,svClkBiasP=getGEPos(TkP,af0,af1,af2,sysType,satPRN,M0,rootA,deltaN,ecc,omega0,Cuc,Cus,Crc,Crs,Cic,Cis,i0,idot,OMEGA0,OMEGADOT,ToE);
		sat_pos_ecefM,VxM,VyM,VzM,RM,TM,NM,svClkBiasM=getGEPos(TkM,af0,af1,af2,sysType,satPRN,M0,rootA,deltaN,ecc,omega0,Cuc,Cus,Crc,Crs,Cic,Cis,i0,idot,OMEGA0,OMEGADOT,ToE);
		vx_num = (sat_pos_ecefP[0] - sat_pos_ecefM[0])/float(0.01);
		vy_num = (sat_pos_ecefP[1] - sat_pos_ecefM[1])/float(0.01);
		vz_num = (sat_pos_ecefP[2] - sat_pos_ecefM[2])/float(0.01);
		#print "vx_num,vy_num,vz_num",vx_num,vy_num,vz_num;
#		print "P"+satPRN+" "+yyyy+" "+mm+" "+dd+" "+hour+" "+min+" "+sec+" "+str(float(sat_pos_ecef[0]))+" "+str(float(sat_pos_ecef[1]))+" "+str(float(sat_pos_ecef[2]))+" "+str(float(svClkBias)*float(1000000));
		eceffile.write("P"+satPRN+" "+yyyy+" "+mm+" "+dd+" "+hour+" "+min+" "+sec+" "+str(float(sat_pos_ecef[0]))+" "+str(float(sat_pos_ecef[1]))+" "+str(float(sat_pos_ecef[2]))+" "+str(vx_num[0])+" "+str(vy_num[0])+" "+str(vz_num[0])+" "+str(float(svClkBias)*float(1000000))+"\n");
		#eceffile.write("P"+satPRN+" "+GLOtime.strftime("%Y")+" "+GLOtime.strftime("%m")+" "+GLOtime.strftime("%d")+" "+GLOtime.strftime("%H")+" "+GLOtime.strftime("%M")+" "+GLOtime.strftime("%S")+" "+str(float(x_te*1000))+" "+str(float(y_te*1000))+" "+str(float(z_te*1000))+" "+str(vx_te*1000)+" "+str(vy_te*1000)+" "+str(vz_te*1000)+" "+str(float(svClkBias)*float(1000000))+"\n");
		#print satPRN+" "+yyyy+" "+mm+" "+dd+" "+hour+" "+min+" "+sec+" "+sysWeek+" "+M0+" "+rootA+" "+deltaN+" "+ecc+" "+omega0+" "+Cuc+" "+Cus+" "+Crc+" "+Crs+" "+Cic+" "+Cis+" "+i0+" "+idot+" "+OMEGA0+" "+OMEGADOT+" "+ToE;
	
	elif (line[0] == 'R'):# or line[0] == 'E' or line[0] == 'J'):# or line[0] == 'R' or  or line[0] == 'C' or  or line[0] == 'I' or line[0] == 'S'):
		thisIndex = brdc.index(line);
		sysType = line[0];
		satsRead.append(line[0:3])
		#print brdc[brdc.index(line)].split();
       
	   # if (sysType == 'G' or sysType == 'C' or sysType == 'E' or sysType == 'I' or sysType == 'J'):
			# lineCount=8;
        # elif (sysType == 'R' or sysType == 'S'):
			# lineCount=4;
			
		satPRN = line[0:3];
		yyyy = line[4:8] #GPS time in case of GPS and QZSS, UTC+3 in case of GLONASS, GAL time in case of GALILEO
		mm = line[9:11]
		dd = line[12:14]
		hour = line[15:17]
		min = line[18:20]
		sec = line[21:23]
		
		#Convert GLONASS time (UTC+3) to GPS time (UTC-LeapSecNow)
		GLOtime = datetime.strptime(yyyy+" "+mm+" "+dd+" "+hour+" "+min+" "+sec,"%Y %m %d %H %M %S");
		J2000 = datetime.strptime("2000 01 01 12 00 00","%Y %m %d %H %M %S");
		if(debug==1):showVar(GLOtime);
		GPStime = GLOtime + timedelta(seconds=LeapSecNow)
		#GPStime = GPStime + timedelta(hours=-3)
		#print "just read GLOtime:",GLOtime,"equivalentGPStime is:",GPStime;
		#gps - utc = 17s, glonass = utc + 3;
		if(debug==1):showVar(GPStime);
		yyyy = GPStime.strftime("%Y")
		mm = GPStime.strftime("%m")
		dd = GPStime.strftime("%d")
		hour = GPStime.strftime("%H")
		min = GPStime.strftime("%M")
		sec = GPStime.strftime("%S")
		#############
		
		svClkBias = brdc[thisIndex+0][23:42]; #[sec]
		
		pz90_pos = numpy.zeros([3,1],float);
		pz90_pos[0,0] = brdc[thisIndex+1][4:23];#[km]
		pz90_pos[1,0] = brdc[thisIndex+2][4:23];#[km]
		pz90_pos[2,0] = brdc[thisIndex+3][4:23];#[km]
		
		####Begin block for orbit interpolation####
		te = GLOtime;
		x_te = pz90_pos[0,0];
		y_te = pz90_pos[1,0];
		z_te = pz90_pos[2,0];
		vx_te = float(brdc[thisIndex+1][23:42]);
		vy_te = float(brdc[thisIndex+2][23:42]);
		vz_te = float(brdc[thisIndex+3][23:42]);
		ddx_te = float(brdc[thisIndex+1][42:61]);
		ddy_te = float(brdc[thisIndex+2][42:61]);
		ddz_te = float(brdc[thisIndex+3][42:61]);
		
		te_minus_3h = te + timedelta(hours=-3); #UTC
		sec_of_day = (te_minus_3h - te_minus_3h.replace(hour=0, minute=0, second=0, microsecond=0)).total_seconds(); #elapsed sec
		utc_midnight = te_minus_3h + timedelta(seconds=-sec_of_day); #UTC at midnight
		JD0 = jdcal.gcal2jd(utc_midnight.year,utc_midnight.month,utc_midnight.day);
		#print JD0, str(JD0);
		JD0 = float(JD0[0])+float(JD0[1]);
		#print JD0, str(JD0);
		UT = te_minus_3h;
		Tu = (float(JD0)-float(2451542.0))/float(36525);
		theta_G0 = float(24110.54841) + float(8640184.812866)*Tu + float(0.093104)*(Tu*Tu) - float(6.2e-6)*(Tu*Tu*Tu); #theta_G0 (radians) is the greenwich mean sidereal time at midnight of the day specified by te
		theta_G0 = 86400 - (abs(theta_G0)%86400);
		theta_G0 = 2*pi * (theta_G0/86400); #radian
		omegae = float(7.2921151467e-5);
		theta_Ge = theta_G0 + omegae*(float(sec_of_day)); #radian
		
		#Coordinates transformation to an inertial reference frame:
		xa_te = x_te*numpy.cos(theta_Ge) - y_te*numpy.sin(theta_Ge);
		ya_te = x_te*numpy.sin(theta_Ge) + y_te*numpy.cos(theta_Ge);
		za_te = z_te;
		
		vxa_te = vx_te*numpy.cos(theta_Ge) - vy_te*numpy.sin(theta_Ge) - omegae*ya_te;
		vya_te = vx_te*numpy.sin(theta_Ge) + vy_te*numpy.cos(theta_Ge) + omegae*xa_te;
		vza_te = vz_te;
		
		u = float(398600.44); #km^3/s^2
		ae= float(6378.136);
		C20= float(-1082.63e-6);
		w3= float(0.7292115e-4);
		r = numpy.sqrt(xa_te*xa_te + ya_te*ya_te + za_te*za_te);
		ubar = u/pow(r,2);
		xdot2 = -ubar*(xa_te/r) + (3/2)*C20*ubar*(xa_te/r)*pow(ae/r,2)*(1-5*(pow(za_te/r,2))) + ddx_te*numpy.cos(theta_Ge) - ddy_te*numpy.sin(theta_Ge);
		ydot2 = -ubar*(ya_te/r) + (3/2)*C20*ubar*(ya_te/r)*pow(ae/r,2)*(1-5*(pow(za_te/r,2))) + ddx_te*numpy.sin(theta_Ge) + ddy_te*numpy.cos(theta_Ge);
		zdot2 = -ubar*(za_te/r) + (3/2)*C20*ubar*(za_te/r)*pow(ae/r,2)*(3-5*(pow(za_te/r,2))) + ddz_te;
		
		xyz = [xa_te, vxa_te, ya_te, vya_te, za_te, vza_te]; #Original
		xyz2 = [x_te*1000, vx_te*1000, y_te*1000, vy_te*1000, z_te*1000, vz_te*1000]; #Test

		#print "xyz_xayaza: ";
		#for i in range(0,6):
	#		print float(xyz2[i])-float(xyz[i]);
	#	print "***\n";
		
		#########1. Forward integration, step:30s, duration: 900s
		t_start1 = (GPStime-J2000).total_seconds();
		t_final1 = ((GPStime + timedelta(seconds=882))-J2000).total_seconds();
		delta_t = 1; #sec
		integ_glnF = ode(glonassfunc).set_integrator('dopri5'); #F for forward
		# Number of time steps: 1 extra for initial condition
		num_steps = numpy.floor((t_final1 - t_start1)/delta_t) + 1;
		# Set initial condition(s): for integrating variables and time!
		integ_glnF.set_initial_value(xyz, t_start1)
		# Additional Python step: create vectors to store trajectories
		t_out = numpy.zeros((num_steps, 1));
		t_out[0] = t_start1;
#		showVar(t_start1);
#		showVar(xa_te);
#		showVar(ya_te);
#		showVar(za_te);
#		showVar(vxa_te);
#		showVar(vya_te);
#		showVar(vza_te);
		x_out = numpy.zeros((num_steps, 1));
		x_out[0] = xa_te; #p
		y_out = numpy.zeros((num_steps, 1));
		y_out[0] = ya_te; #p
		z_out = numpy.zeros((num_steps, 1));
		z_out[0] = za_te; #p
		xdot_out = numpy.zeros((num_steps, 1));
		xdot_out[0] = vxa_te; #v
		ydot_out = numpy.zeros((num_steps, 1));
		ydot_out[0] = vya_te; #v
		zdot_out = numpy.zeros((num_steps, 1));
		zdot_out[0] = vza_te; #v
		# Integrate the ODE(s) across each delta_t timestep
		k = 1
		while integ_glnF.successful() and k < num_steps:
			integ_glnF.integrate(integ_glnF.t + delta_t)
 			# Store the results to plot later
			t_out[k] = integ_glnF.t
			x_out[k] = integ_glnF.y[0]
			y_out[k] = integ_glnF.y[2] #changed
			z_out[k] = integ_glnF.y[4]
			xdot_out[k] = integ_glnF.y[1]
			ydot_out[k] = integ_glnF.y[3]
			zdot_out[k] = integ_glnF.y[5]
			k += 1
		forwardT = GPStime + timedelta(seconds=882);
		#print "forwardT",forwardT,"equivalentGPSTime",GPStime;
		forwardX = float(x_out[-1:]);
		forwardY = float(y_out[-1:]);
		forwardZ = float(z_out[-1:]);
		#########################################################################

		#########2. Backwards integration, step:30s, duration: 900s
		t_start2 = (GPStime-J2000).total_seconds();
		t_final2 = ((GPStime + timedelta(seconds=-18))-J2000).total_seconds();
		delta_t = -1; #sec
		integ_glnB = ode(glonassfunc).set_integrator('dopri5'); #B for backwards
		# Number of time steps: 1 extra for initial condition
		num_steps = numpy.floor((t_final2 - t_start2)/delta_t) + 1;
		# Set initial condition(s): for integrating variables and time!
		integ_glnB.set_initial_value(xyz, t_start1)
		# Additional Python step: create vectors to store trajectories
		t_out = numpy.zeros((num_steps, 1));
		t_out[0] = t_start2;
		x_out = numpy.zeros((num_steps, 1));
		x_out[0] = xa_te; #p
		y_out = numpy.zeros((num_steps, 1));
		y_out[0] = ya_te; #p
		z_out = numpy.zeros((num_steps, 1));
		z_out[0] = za_te; #p
		xdot_out = numpy.zeros((num_steps, 1));
		xdot_out[0] = vxa_te; #v
		ydot_out = numpy.zeros((num_steps, 1));
		ydot_out[0] = vya_te; #v
		zdot_out = numpy.zeros((num_steps, 1));
		zdot_out[0] = vza_te; #v
		# Integrate the ODE(s) across each delta_t timestep
		k = 1
		while integ_glnB.successful() and k < num_steps:
			integ_glnB.integrate(integ_glnB.t + delta_t)
 			# Store the results to plot later
			t_out[k] = integ_glnB.t
			x_out[k] = integ_glnB.y[0]
			y_out[k] = integ_glnB.y[2] #changed
			z_out[k] = integ_glnB.y[4]
			xdot_out[k] = integ_glnB.y[1]
			ydot_out[k] = integ_glnB.y[3]
			zdot_out[k] = integ_glnB.y[5]
			k += 1
		backwardT = GPStime + timedelta(seconds=-18);
		backwardX = float(x_out[-1:]);
		backwardY = float(y_out[-1:]);
		backwardZ = float(z_out[-1:]);
		#########################################################################

		#print "back",sysType, satPRN, backwardT, backwardX, backwardY, backwardZ;
		#print sysType, satPRN, GPStime, xa_te, ya_te, za_te;
		#print "forward",sysType, satPRN, forwardT, forwardX, forwardY, forwardZ;
		#print "#######";
		
		
		#theta_G computation for forwardT (theta_Gf) and backwardT (theta_Gb):
		#1. for forwardT:
		tempT = te + timedelta(hours=-3) + timedelta(seconds=883); #UTC
		sec_of_day = (tempT - tempT.replace(hour=0, minute=0, second=0, microsecond=0)).total_seconds(); #elapsed sec
		utc_midnight = tempT + timedelta(seconds=-sec_of_day); #UTC at midnight
		JD0 = jdcal.gcal2jd(utc_midnight.year,utc_midnight.month,utc_midnight.day);
		#print JD0, str(JD0);
		JD0 = float(JD0[0])+float(JD0[1]);
		#print JD0, str(JD0);
		UT = te + timedelta(hours=-3) + timedelta(seconds=883);
		Tu = (float(JD0)-float(2451542.0))/float(36525);
		theta_G0 = float(24110.54841) + float(8640184.812866)*Tu + float(0.093104)*(Tu*Tu) - float(6.2e-6)*(Tu*Tu*Tu); #theta_G0 (radians) is the greenwich mean sidereal time at midnight of the day specified by te
		theta_G0 = 86400 - (abs(theta_G0)%86400);
		theta_G0 = 2*pi * (theta_G0/86400); #radian
		omegae = float(7.2921151467e-5);
		theta_Gf = theta_G0 + omegae*(float(sec_of_day)); #radian
		#2. for backwardT:
		tempT = te + timedelta(hours=-3) + timedelta(seconds=-17); #UTC
		sec_of_day = (tempT - tempT.replace(hour=0, minute=0, second=0, microsecond=0)).total_seconds(); #elapsed sec
		utc_midnight = tempT + timedelta(seconds=-sec_of_day); #UTC at midnight
		JD0 = jdcal.gcal2jd(utc_midnight.year,utc_midnight.month,utc_midnight.day);
		#print JD0, str(JD0);
		JD0 = float(JD0[0])+float(JD0[1]);
		#print JD0, str(JD0);
		UT = te + timedelta(hours=-3) + timedelta(seconds=-17);
		Tu = (float(JD0)-float(2451542.0))/float(36525);
		theta_G0 = float(24110.54841) + float(8640184.812866)*Tu + float(0.093104)*(Tu*Tu) - float(6.2e-6)*(Tu*Tu*Tu); #theta_G0 (radians) is the greenwich mean sidereal time at midnight of the day specified by te
		theta_G0 = 86400 - (abs(theta_G0)%86400);
		theta_G0 = 2*pi * (theta_G0/86400); #radian
		omegae = float(7.2921151467e-5);
		theta_Gb = theta_G0 + omegae*(float(sec_of_day)); #radian
		
		#Change the forwards and backwards back to PZ-90:
		forwardX_pz90 = forwardX*numpy.cos(theta_Gf) + forwardY*numpy.sin(theta_Gf);
		forwardY_pz90 = -forwardX*numpy.sin(theta_Gf) + forwardY*numpy.cos(theta_Gf);
		forwardZ_pz90 = forwardZ;
		backwardX_pz90 = backwardX*numpy.cos(theta_Gb) + backwardY*numpy.sin(theta_Gb);
		backwardY_pz90 = -backwardX*numpy.sin(theta_Gb) + backwardY*numpy.cos(theta_Gb);
		backwardZ_pz90 = backwardZ;
		#print "in pz-90:";
		#print sysType, satPRN, backwardT, backwardX_pz90, backwardY_pz90, backwardZ_pz90;
		#print sysType, satPRN, GPStime, x_te, y_te, z_te;
		#print sysType, satPRN, forwardT, forwardX_pz90, forwardY_pz90, forwardZ_pz90;
		#print "####\n"; 	

		pz90_pos_forward = numpy.zeros([3,1],float);
		pz90_pos_forward[0,0] = forwardX_pz90*1000;#[m]
		pz90_pos_forward[1,0] = forwardY_pz90*1000;#[m]
		pz90_pos_forward[2,0] = forwardZ_pz90*1000;#[m]
		
		pz90_pos_backward = numpy.zeros([3,1],float);
		pz90_pos_backward[0,0] = backwardX_pz90*1000;#[m]
		pz90_pos_backward[1,0] = backwardY_pz90*1000;#[m]
		pz90_pos_backward[2,0] = backwardZ_pz90*1000;#[m]

		
		#Bring the forwards and backwards in coincidence with WGS-84 by rotation:
		tra_pz92_wgs84 = numpy.zeros([3,1],float);
		tra_pz92_wgs84[0,0] = float(-0.36); #m
		tra_pz92_wgs84[1,0] = float(0.08); #m
		tra_pz92_wgs84[2,0] = float(0.18); #m
		
		rot_pz90 = numpy.zeros([3,3],float);
		rot_pz90[0,0] = float(-3e-9);
		rot_pz90[0,1] = float(-353 * (4.84813681e-9));
		rot_pz90[0,2] = float(-4 * (4.84813681e-9));
		rot_pz90[1,0] = float(353 * (4.84813681e-9));
		rot_pz90[1,1] = float(-3e-9);
		rot_pz90[1,2] = float(19 * (4.84813681e-9));
		rot_pz90[2,0] = float(4 * (4.84813681e-9));
		rot_pz90[2,1] = float(-19 * (4.84813681e-9));
		rot_pz90[2,2] = float(-3e-9);
		
		tra_pz90 = numpy.zeros([3,1],float);
		tra_pz90[0,0] = float(0.07) #[m]
		tra_pz90[1,0] = float(-0.0) #[m]
		tra_pz90[2,0] = float(-0.77) #[m]

		forward_ecef = numpy.zeros([3,1],float); #ecef coincident with WGS84
		forward_ecef = (pz90_pos_forward + numpy.dot(rot_pz90,pz90_pos_forward) + tra_pz90)+tra_pz92_wgs84;
		forward_ecef = forward_ecef/1000; #convert to [km]

		backward_ecef = numpy.zeros([3,1],float); #ecef coincident with WGS84
		backward_ecef = (pz90_pos_backward + numpy.dot(rot_pz90,pz90_pos_backward) + tra_pz90)+tra_pz92_wgs84;
		backward_ecef = backward_ecef/1000; #convert to [km]
				
		#print satPRN+" "+yyyy+" "+mm+" "+dd+" "+hour+" "+min+" "+sec+" "+sysWeek
		#print sat_pos_ecef;
		#print svClkBias;
	#	print "P"+satPRN+" "+yyyy+" "+mm+" "+dd+" "+hour+" "+min+" "+sec+" "+str(float(sat_pos_ecef[0]))+" "+str(float(sat_pos_ecef[1]))+" "+str(float(sat_pos_ecef[2]))+" "+str(float(svClkBias)*float(1000000));
		
			
		#just multiply with velocity (wrong way):
		backwardX = (x_te*1000) + (vx_te*1000)*(-18);
		backwardY = (y_te*1000) + (vy_te*1000)*(-18);
		backwardZ = (z_te*1000) + (vz_te*1000)*(-18);
		
		#print "./glo_rk4.py "+str(x_te*1000)+" "+str(y_te*1000)+" "+str(z_te*1000)+" "+str(vx_te*1000)+" "+str(vy_te*1000)+" "+str(vz_te*1000)+" "+str(ddx_te*1000)+" "+str(ddy_te*1000)+" "+str(ddz_te*1000);
		backwardPOS = commands.getoutput("./glo_rk4.py "+str(x_te*1000)+" "+str(y_te*1000)+" "+str(z_te*1000)+" "+str(vx_te*1000)+" "+str(vy_te*1000)+" "+str(vz_te*1000)+" "+str(ddx_te*1000)+" "+str(ddy_te*1000)+" "+str(ddz_te*1000));
		#print "backwardPOS: ", backwardPOS;
		backwardX = float(backwardPOS.split()[0]);
		backwardY = float(backwardPOS.split()[1]);
		backwardZ = float(backwardPOS.split()[2]);
		
		#backwardX = backward_ecef[0,0]; #even worse than the above
		#backwardY = backward_ecef[1,0];
		#backwardZ = backward_ecef[2,0];
		
			
		#print "GLO RTN: ",XYZ2RTN([x_te*1000,y_te*1000,z_te*1000],[vx_te*1000,vy_te*1000,vz_te*1000]); # [m] and [m/s]
		#print "backwardT:",backwardT;
		eceffile.write("P"+satPRN+" "+GLOtime.strftime("%Y")+" "+GLOtime.strftime("%m")+" "+GLOtime.strftime("%d")+" "+GLOtime.strftime("%H")+" "+GLOtime.strftime("%M")+" "+GLOtime.strftime("%S")+" "+str(float(backwardX))+" "+str(float(backwardY))+" "+str(float(backwardZ))+" "+str(vx_te)+" "+str(vy_te)+" "+str(vz_te)+" "+str(float(svClkBias)*float(1000000))+"\n");
		#eceffile.write("P"+satPRN+" "+backwardT.strftime("%Y")+" "+backwardT.strftime("%m")+" "+backwardT.strftime("%d")+" "+backwardT.strftime("%H")+" "+backwardT.strftime("%M")+" "+backwardT.strftime("%S")+" "+str(float(backwardX))+" "+str(float(backwardY))+" "+str(float(backwardZ))+" "+str(float(svClkBias)*float(1000000))+"\n");
		#eceffile.write("P"+satPRN+" "+forwardT.strftime("%Y")+" "+forwardT.strftime("%m")+" "+forwardT.strftime("%d")+" "+forwardT.strftime("%H")+" "+forwardT.strftime("%M")+" "+forwardT.strftime("%S")+" "+str(float(forwardX))+" "+str(float(forwardY))+" "+str(float(forwardZ))+" "+str(float(svClkBias)*float(1000000))+"\n");
#print set(satsRead);

eceffile.close();



print "Script end time: ",os.system("date");

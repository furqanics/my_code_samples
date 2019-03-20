#!/usr/bin/python2.7

import commands;
import sys;
import os;
import numpy as np;
from datetime import datetime as dt;
from datetime import timedelta;
import time

#First change made

def toYearFraction(date):
    def sinceEpoch(date): # returns seconds since epoch
        return time.mktime(date.timetuple())
    s = sinceEpoch

    year = date.year
    startOfThisYear = dt(year=year, month=1, day=1)
    startOfNextYear = dt(year=year+1, month=1, day=1)

    yearElapsed = s(date) - s(startOfThisYear)
    yearDuration = s(startOfNextYear) - s(startOfThisYear)
    fraction = yearElapsed/yearDuration

    return date.year + fraction


dateProc = sys.argv[1]; #should be "YYYY-MM-DD" format
tmpSubDir = sys.argv[2]+"/";

Y4 = dt.strptime(dateProc,"%Y-%m-%d").strftime('%Y');
Y2 = Y4[2:];
m  = dt.strptime(dateProc,"%Y-%m-%d").strftime('%m');
d  = dt.strptime(dateProc,"%Y-%m-%d").strftime('%d');

fC = "GPS1A_"+Y4+"-"+m+"-"+d+"_C_00";
fD = "GPS1A_"+Y4+"-"+m+"-"+d+"_D_00";
fCD = [fC,fD];
print "fCD: ",fCD;

extensions = [".azi.gz",".ele.gz",".d12.gz",".i12.gz",".m12.gz",".m21.gz",".sn1.gz",".sn2.gz"]; #"."+Y2+"S.gz",
extensions = [".m12",".m21",".sn1",".sn2"];#,".azi",".ele"]; #"."+Y2+"S.gz",

#Desired output format: GPS_YYYYMMDDhhmmss, GRC-id, GPS-id, GPS-ele, GPS-azi, GPS-mp1, GPS-mp2, GPS-sn1, GPS-sn2
#ref_time = 12:00:00 noon of January 1, 2000 in GPS Time
ref_time = dt.strptime('2000-01-01 12:00:00','%Y-%m-%d %H:%M:%S');
gps_start_time = "";

for f in fCD:
	print "f: ",tmpSubDir,f;
	for ext in extensions:
		if (os.path.exists(tmpSubDir+f+ext)==False):
                	print f+ext+" missing!";
                	continue;
		ex = ext;
		qty = '';
		if ex[1:] == 'm12': qty = 'MP1';
                if ex[1:] == 'm21': qty = 'MP2';
                if ex[1:] == 'sn1': qty = 'SN1';
                if ex[1:] == 'sn2': qty = 'SN2';
                if ex[1:] == 'azi': qty = 'AZI';
                if ex[1:] == 'ele': qty = 'ELE';

		print "Going to create ",tmpSubDir+f+"_"+qty+".txt",commands.getoutput('date');
		outfile = open(tmpSubDir+f+"_"+qty+".txt",'w');
		#print rPath+f+ext;
		#os.system("cp "+rPath+f+ext+" .")
		#os.system("gunzip "+f+ext)
		#Get start time and 10-sec epoch timestamps from m12 file;
		if (ext):#== ".m12.gz"):
			gps_start_time = ' '.join(commands.getoutput("grep GPS_START_TIME "+f+ex).split()[1:]);
			startH = gps_start_time.split()[3];
			startM = gps_start_time.split()[4];
			startS = gps_start_time.split()[5].split('.')[0];
			com3file = open(tmpSubDir+f+ext,'r').readlines();
			epoch_stamps = [];
			for line in com3file[2:]:
				if("G" in line[15:] or line[12:15]=='-1\n'):
					epoch_stamps.append(line);
					#print line;
			#print "epoch read complete";
			for epoch in epoch_stamps:
				#print "epoch: ",epoch[:-1];
				Y4_epoch = '00';
				Y2_epoch = '00';
				m_epoch = '00';
				d_epoch = '00';
				s_epoch = '00';
				start_time = dt.strptime(Y4+m+d+" "+startH+":"+startM+":"+startS,'%Y%m%d %H:%M:%S');
                                sec = epoch.split()[0];
				epoch_time = start_time + timedelta(seconds=float(sec));
				Y4_epoch = epoch_time.strftime('%Y');
                                Y2_epoch = Y4[0:2];
                                m_epoch = epoch_time.strftime('%m');
                                d_epoch = epoch_time.strftime('%d');
				h_epoch = epoch_time.strftime('%H');
				min_epoch = epoch_time.strftime('%M');
				epoch_time_full = epoch_time.strftime('%Y %m %d %H %M %S');
                                s_epoch = epoch_time.strftime('%S');
				sec_int = int(sec.split('.')[0]);
				sec_frac = float(sec.split('.')[1])*float(1e-6); #microsec
				this_time = dt.strptime(gps_start_time.split()[0]+"-"+gps_start_time.split()[1].zfill(2)+"-"+gps_start_time.split()[2].zfill(2)+" "+startH+":"+startM+":"+startS,"%Y-%m-%d %H:%M:%S") + timedelta(seconds = sec_int);
				this_sec_int = int((this_time - ref_time).total_seconds());
				#print this_time,this_sec_int,sec_frac;
				#print gps_start_time.split()[0]+"-"+gps_start_time.split()[1].zfill(2)+"-"+gps_start_time.split()[2].zfill(2)+" "+gps_start_time.split()[3].zfill(2)+":"+gps_start_time.split()[4].zfill(2)+":"+gps_start_time.split()[5].split('.')[0].zfill(2);
				if(epoch[12:15]!='-1\n'):nsat = int(epoch.split()[1]);
				if(epoch[12:15]!='-1\n'):sats = epoch.split()[2:];
				grace_id = f.split("_")[2];
				values = com3file[com3file.index(epoch)+1].split();
		                thisDate = dt.strptime(Y4+m+d,"%Y%m%d");
				for i in range(nsat):
					this_val = values[i];
					val = values[i];
					if "S" in val: this_val = val.replace("S", "");
					#print gps_start_time+" "+sec+" "+str(this_sec_int)+" "+str(sec_frac)+" "+grace_id+" "+sats[i]+" "+this_val+" "+ex[1:];
					if(float(this_val)!=0):
						this_az = float(0);
						this_el = float(0);
						grepOut = commands.getoutput("grep \'"+str(this_sec_int)+" "+grace_id+" "+sats[i]+"\' "+sats[i]+"_"+grace_id+"_"+dateProc+".azelsnr");
						if "No such file or directory" not in grepOut and grepOut != "":
							#print grepOut;
							this_az = float(grepOut.split()[3]);
							this_el = float(grepOut.split()[4]);
						outfile.write(str(toYearFraction(epoch_time))+" "+epoch_time_full+" "+str(this_sec_int)+" "+sec+" "+str(sec_frac)+" "+grace_id+" "+sats[i]+" "+this_val+" "+qty+" "+str(this_az)+" "+str(this_el)+"\n");
		outfile.close();
		
		print "Created ",tmpSubDir+f+"_"+qty+".txt",commands.getoutput('date');
		#os.system("rm "+f+ex);
'''
580563610 C G02 221.628570557 3.82090640068 19 19 94
580563620 C G02 221.434051514 3.29394173622 19 18 92
[ahmed@fred tmp_2018-05-25_151009]$ more G02_C_2018-05-25.azelsnr
'''

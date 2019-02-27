#!/usr/bin/python2.7

import os
import sys
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import commands
from datetime import date, timedelta
from datetime import datetime as dt

dateProc = sys.argv[1];

QTYS = ["MP1","MP2","SN1","SN2"];
UNTS = ["m","m","dB","dB"];

gfoC_file = open("GPS1A_C_"+dateProc+".txt",'w');
gfoD_file = open("GPS1A_D_"+dateProc+".txt",'w');

daily_metrics_file = open("GPS1A_"+dateProc+".txt",'w');

pObsMissCD = [];
nObsPossibleCD = [];
nObsAvailableCD = [];

#pObsMissTotal = np.sum(nObsPossibleCD)/np.sum(nObsAvailableCD)*100;

elevCutOff = 15;

for GFOsat in ['C','D']:
	this_gfo_file = open("GPS1A_"+GFOsat+"_"+dateProc+".txt",'w');
	RepFile = "GPS1A_"+dateProc+"_"+GFOsat+"_00."+dateProc[2:4]+"S";
	#print commands.getoutput("grep \"Poss. # of obs epochs\" "+RepFile).split();
        #print commands.getoutput("grep \"Epochs w/ observations\" "+RepFile).split();
	nPossibleObs = float(commands.getoutput("grep \"Poss. # of obs epochs\" "+RepFile).split()[6]);
	print GFOsat,"nPossibleObs",nPossibleObs;
	nAvailableObs = float(commands.getoutput("grep \"Epochs w/ observations\" "+RepFile).split()[4]);
        print GFOsat,"nAvailableObs",nAvailableObs;
	pObsMiss = (nPossibleObs - nAvailableObs)/nPossibleObs*100;
	pObsMissCD.append(pObsMiss);
	nObsPossibleCD.append(nPossibleObs);
	nObsAvailableCD.append(nAvailableObs);
	all_prn_mp1 = [];
        all_prn_mp2 = [];
        all_prn_sn1 = [];
        all_prn_sn2 = [];
	for PRN in range(1,33):
		this_prn_mp1 = [];
		this_prn_mp2 = [];
		this_prn_sn1 = [];
		this_prn_sn2 = [];
		for QTY in range(0,4):
			#qty_mean_file = open("mean_"+dateProc+"_"+GFOsat+"_"+QTYS[QTY]+".txt",'w');

			print "Plotting GPS1A_"+dateProc+"_"+GFOsat+"_"+'G'+str(PRN).zfill(2)+"_"+QTYS[QTY]+".png ...";
			prn_id = 'G'+str(PRN).zfill(2);

			#qfname = "GPS1A_"+dateProc+"_"+GFOsat+"_00_"+QTYS[QTY]+".txt";
			####Accumulate begin###
			datesTillNow = [];
			filesTillNow = [];
			startdate = date(2018,5,22);
			endyear = int(dateProc.split('-')[0]);
			endmonth = int(dateProc.split('-')[1]);
			endday = int(dateProc.split('-')[2]);
			#for startdate < date(int(endyear), int(endmonth), int(endday)):
			while startdate < date(int(endyear), int(endmonth), int(endday)):
			        doy = startdate.strftime('%j');
			        dd = startdate.strftime('%d');
			        mm = startdate.strftime('%m');
			        y4 = startdate.strftime('%Y');
			        y2 = startdate.strftime('%y');
			        #print startdate, y4, mm, dd, y2;
			        thisDate = dt.strptime(y4+mm+dd,"%Y%m%d");
			        datesTillNow.append(dt.strftime(thisDate,'%Y-%m-%d'));
			        startdate += timedelta(days=1)
#			print datesTillNow;
#			print "cat "+' '.join(dateTillNow);
			for d in datesTillNow:
				#original: if(os.path.exists("/proj/L1A_QC/GPS/"+d+"/GPS1A_"+d+"_"+GFOsat+"_00_"+QTYS[QTY]+".txt")==True):
				if(os.path.exists("/home/ahmed/GPS2/"+d+"/GPS1A_"+d+"_"+GFOsat+"_00_"+QTYS[QTY]+".txt")==True):
					filesTillNow.append("/home/ahmed/GPS2/"+d+"/GPS1A_"+d+"_"+GFOsat+"_00_"+QTYS[QTY]+".txt");
			filesTillNow.append("GPS1A_"+dateProc+"_"+GFOsat+"_00_"+QTYS[QTY]+".txt");
                        print "cat "+' '.join(filesTillNow)+" > accumulated_"+GFOsat+"_"+QTYS[QTY]+"_"+dateProc+".txt";
                        os.system("cat "+' '.join(filesTillNow)+" > accumulated_"+GFOsat+"_"+QTYS[QTY]+"_"+dateProc+".txt");
			qfname = "accumulated_"+GFOsat+"_"+QTYS[QTY]+"_"+dateProc+".txt";
			####Accumulate end###

			
                        #print "awk \'{print $9}\' "+qfname;
			#epochs = commands.getoutput("awk \'{print $9}\' "+qfname).split();
			#print len(epochs);
			#for all in open(qfname,'r').readlines():
			#	epochs.append(all.split()[8]);
			
			#epochs = list(set(epochs));
			#print len(epochs);

			#for epoch in epochs:
			#	x = commands.getoutput("awk \'{if($9=="+epoch+") print $13}\' "+qfname).split();
			#	print x;
			#	qty_epoch_val = [];
			#	for i in x:
			#		qty_epoch_val.append(float(i));
			#	qty_mean_file.write(epoch+" "+str(np.average(qty_epoch_val))+"\n");
			#qty_mean_file.close();

			tfname = prn_id+"_"+GFOsat+"_"+QTYS[QTY]+".tmp";
			os.system("grep "+prn_id+" "+qfname+" > "+tfname);
			if(os.stat(tfname).st_size != 0):
				QTY_time = [];
				QTY_val = [];

				QTY_time_above_eo = [];
                                QTY_time_below_eo = [];
				QTY_time_no_eo = [];
				QTY_val_above_eo = [];
				QTY_val_below_eo = [];
				QTY_val_no_eo = [];

				#GNV1A_2018-05-28_C_00.azel
				#580737600 C G01 142.019470215 1.64449393749 19 17 80
				#580737610 C G01 142.091964722 1.07756340504 20 16 75
				#2018.39169616 2018 05 24 00 15 30 580306530 930.0000 0.0 C G29 4.923 MP1
				#2018.39169647 2018 05 24 00 15 40 580306540 940.0000 0.0 C G29 2.919 MP1
				#2018.39169679 2018 05 24 00 15 50 580306550 950.0000 0.0 C G29 2.722 MP1
				for line in open(tfname,'r').readlines():
					azi = float(line.split()[14]);
					ele = float(line.split()[15]);
					y_4 = line.split()[1];
                                        m_2 = line.split()[2];
                                        d_2 = line.split()[3];
					h_2 = line.split()[4];
                                        min_2 = line.split()[5];
                                        s_2 = line.split()[6];
					currDOY = dt.strptime(y_4+" "+m_2+" "+d_2+" "+h_2+" "+min_2+" "+s_2,'%Y %m %d %H %M %S').strftime('%j');
					j2000 = line.split()[7];
					to_grep = j2000;
					to_grep = j2000+" "+GFOsat+" "+str(PRN).zfill(2);
					#580321330 C 31
					#print "to_grep: ",to_grep;
					azelgrep = "";
#					if(os.path.exists("GNV1A_"+dateProc+"_"+GFOsat+"_00.azel")==True):
						#print "exists";
#						azelgrep = commands.getoutput("grep \'"+to_grep+"\' GNV1A_"+dateProc+"_"+GFOsat+"_00.azel");
					#if(os.path.exists("G"+prn_id+"_"+GFOsat+"_"+dateProc+".azelsnr")==True):
						#print "exists";
						#azelgrep = commands.getoutput("grep \'"+to_grep+"\' G"+prn_id+"_"+GFOsat+"_"+dateProc+".azelsnr");
						#azelgrep = commands.getoutput("grep \'"+to_grep+"\' G"+prn_id+"_"+GFOsat+"_"+dateProc+".azelsnr");
					#print "azelgrep: ",azelgrep;
					#azi = -999;
					#ele = -999;
					#if len(azelgrep)!=0:azi=float(azelgrep.split()[3]);
                    #                    if len(azelgrep)!=0:ele=float(azelgrep.split()[4]);
					#print "azel: ",azi,ele,"\n";
					if (ele!=-999 and ele>=elevCutOff):
						QTY_time_above_eo.append(float(((dt.strptime(y_4+" "+m_2+" "+d_2+" "+h_2+" "+min_2+" "+s_2,'%Y %m %d %H %M %S') - dt.strptime('2018 05 22 00 00 00','%Y %m %d %H %M %S')).total_seconds()))/float(86400));
                                        	QTY_val_above_eo.append(float(line.split()[12]));
					elif (ele!=-999 and ele<elevCutOff and ele!=0):
                                                QTY_time_below_eo.append(float(((dt.strptime(y_4+" "+m_2+" "+d_2+" "+h_2+" "+min_2+" "+s_2,'%Y %m %d %H %M %S') - dt.strptime('2018 05 22 00 00 00','%Y %m %d %H %M %S')).total_seconds()))/float(86400));
                                                QTY_val_below_eo.append(float(line.split()[12]));
					elif (ele!=-999 and ele==0):
                                                QTY_time_no_eo.append(float(((dt.strptime(y_4+" "+m_2+" "+d_2+" "+h_2+" "+min_2+" "+s_2,'%Y %m %d %H %M %S') - dt.strptime('2018 05 22 00 00 00','%Y %m %d %H %M %S')).total_seconds()))/float(86400));
                                                QTY_val_no_eo.append(float(line.split()[12]));


	                                if QTY == 0:
						all_prn_mp1.append(float(line.split()[12]));
                                                this_prn_mp1.append(float(line.split()[12]));
						qty_ymin = -25;
                                                qty_ymax = 25;
                                                qty_offset = 20;
                                        if QTY == 1:
						all_prn_mp2.append(float(line.split()[12]));
                                                this_prn_mp2.append(float(line.split()[12]));
						qty_ymin = -25;
                                                qty_ymax = 25;
                                                qty_offset = 20;
                                        if QTY == 2:
						all_prn_sn1.append(float(line.split()[12]));
                                                this_prn_sn1.append(float(line.split()[12]));
						qty_ymin = 0;
                                                qty_ymax = 1600;
                                                qty_offset = 800;
                                        if QTY == 3:
						all_prn_sn2.append(float(line.split()[12]));
                                                this_prn_sn2.append(float(line.split()[12]));
						qty_ymin = 0;
						qty_ymax = 1600;
						qty_offset = 800;

				xmin = int(currDOY);
				QTY_val = QTY_val_above_eo + QTY_val_below_eo + QTY_val_no_eo;
				qty_mean = "%.2f"%np.average(QTY_val);
                                qty_std = "%.2f"%np.std(QTY_val);
				qty_rms = "%.2f"%np.sqrt(np.average(np.square(QTY_val)));
				plotname = "GPS1A_"+dateProc+"_"+GFOsat+"_"+'G'+str(PRN).zfill(2)+"_"+QTYS[QTY]+".png";
				fig = plt.figure();
				ax = fig.add_subplot(111);
		                print "len above: ",len(QTY_time_above_eo),len(QTY_val_above_eo);
                                print "len below: ",len(QTY_time_below_eo),len(QTY_val_below_eo);
                                print "len unknown: ",len(QTY_time_no_eo),len(QTY_val_no_eo);
				ax.plot(QTY_time_above_eo,QTY_val_above_eo, color='green', marker='o', linestyle='none', markersize=0.25, label='Above '+str(elevCutOff)+" deg");
                                ax.plot(QTY_time_below_eo,QTY_val_below_eo, color='red', marker='o', linestyle='none', markersize=0.25, label='Below '+str(elevCutOff)+" deg");
				ax.plot(QTY_time_no_eo,QTY_val_no_eo, color='blue', marker='o', linestyle='none', markersize=0.25, label='Elevation unknown');
				ax.grid(True);
				ax.legend(loc='best', shadow=False)
				if(QTY==2 or QTY==3):plt.ylim(qty_ymin, qty_ymax);
				stamp = str(commands.getoutput('date'));
				plt.xlabel('Days since 2018-05-22 00:00:00');
                                plt.ylabel(QTYS[QTY]+" ["+UNTS[QTY]+"]");
				units = UNTS[QTY];
				plt.title('GFO: '+GFOsat+", "+'GPS: '+prn_id+" (Mean: "+str(qty_mean)+" $\pm$ "+str(qty_std)+" "+units+", RMS: "+str(qty_rms)+" "+units+")");
				ax.annotate ( 'Created: '+str ( commands.getoutput ( 'date' ) ) ,xy= ( 5, 5 ) , xycoords='figure points',fontsize=6);
				#plt.text(26,20,'Created on '+str(commands.getoutput('date')),rotation='vertical');
				fig.savefig(plotname, dpi=300)
				plt.close(fig);
				
	


		if(GFOsat == 'C' and len(this_prn_mp1)!=0):gfoC_file.write(dateProc+" "+GFOsat+" G"+str(PRN).zfill(2)+" "+"%.4f"%(np.average(this_prn_mp1))+" "+"%.4f"%(np.std(this_prn_mp1))+" "+"%.4f"%(np.sqrt(np.average(np.square(this_prn_mp1))))+" "+"%.4f"%(np.average(this_prn_mp2))+" "+"%.4f"%(np.std(this_prn_mp2))+" "+"%.4f"%(np.sqrt(np.average(np.square(this_prn_mp2))))+" "+"%.4f"%(np.average(this_prn_sn1))+" "+"%.4f"%(np.std(this_prn_sn1))+" "+"%.4f"%(np.sqrt(np.average(np.square(this_prn_sn1))))+" "+"%.4f"%(np.average(this_prn_sn2))+" "+"%.4f"%(np.std(this_prn_sn2))+" "+"%.4f"%(np.sqrt(np.average(np.square(this_prn_sn2))))+" "+"%.4f"%(pObsMiss)+"\n");
		if(GFOsat == 'D' and len(this_prn_mp1)!=0):gfoD_file.write(dateProc+" "+GFOsat+" G"+str(PRN).zfill(2)+" "+"%.4f"%(np.average(this_prn_mp1))+" "+"%.4f"%(np.std(this_prn_mp1))+" "+"%.4f"%(np.sqrt(np.average(np.square(this_prn_mp1))))+" "+"%.4f"%(np.average(this_prn_mp2))+" "+"%.4f"%(np.std(this_prn_mp2))+" "+"%.4f"%(np.sqrt(np.average(np.square(this_prn_mp2))))+" "+"%.4f"%(np.average(this_prn_sn1))+" "+"%.4f"%(np.std(this_prn_sn1))+" "+"%.4f"%(np.sqrt(np.average(np.square(this_prn_sn1))))+" "+"%.4f"%(np.average(this_prn_sn2))+" "+"%.4f"%(np.std(this_prn_sn2))+" "+"%.4f"%(np.sqrt(np.average(np.square(this_prn_sn2))))+" "+"%.4f"%(pObsMiss)+"\n");


#		aa
	pObsMissTotal = (float(np.sum(nObsPossibleCD))-float(np.sum(nObsAvailableCD)))/float(np.sum(nObsPossibleCD))*100;
	daily_metrics_file.write(dateProc+" "+GFOsat+" "+"%.4f"%(np.average(all_prn_mp1))+" "+"%.4f"%(np.std(all_prn_mp1))+" "+"%.4f"%(np.sqrt(np.average(np.square(all_prn_mp1))))+" "+"%.4f"%(np.average(all_prn_mp2))+" "+"%.4f"%(np.std(all_prn_mp2))+" "+"%.4f"%(np.sqrt(np.average(np.square(all_prn_mp2))))+" "+"%.4f"%(np.average(all_prn_sn1))+" "+"%.4f"%(np.std(all_prn_sn1))+" "+"%.4f"%(np.sqrt(np.average(np.square(all_prn_sn1))))+" "+"%.4f"%(np.average(all_prn_sn2))+" "+"%.4f"%(np.std(all_prn_sn2))+" "+"%.4f"%(np.sqrt(np.average(np.square(all_prn_sn2))))+" "+"%.4f"%(pObsMissTotal)+"\n");

	this_gfo_mp1_hist ="GPS1A_"+dateProc+"_"+GFOsat+"_MP1.png";
	plotname = this_gfo_mp1_hist;
	fig = plt.figure();
        ax = fig.add_subplot(111);
        ax.hist(all_prn_mp1, bins=np.arange(-10, 10, 0.1));
        ax.grid(True);
        plt.xlabel('Value of MP1 [m]');
        plt.ylabel('# of epochs');
        units = "m";
	mp1_mean = "%.2f"%np.average(all_prn_mp1);
        mp1_std = "%.2f"%np.std(all_prn_mp1);
        mp1_rms = "%.2f"%np.sqrt(np.average(np.square(all_prn_mp1)));
        plt.xlim(-10,10);
	#plt.ylim(0,30000);
	plt.title('GFO_'+GFOsat+" (Mean: "+str(mp1_mean)+" $\pm$ "+str(mp1_std)+" "+units+", RMS: "+str(mp1_rms)+" "+units+")");
	ax.annotate ( 'Created: '+str ( commands.getoutput ( 'date' ) ) ,xy= ( 5, 5 ) , xycoords='figure points',fontsize=6);
        #plt.text(27,0,'Created on '+str(commands.getoutput('date')),horizontalalignment='right',verticalalignment='bottom',rotation='vertical');
        fig.savefig(plotname, dpi=300)
        plt.close(fig);

	this_gfo_mp2_hist ="GPS1A_"+dateProc+"_"+GFOsat+"_MP2.png";
        plotname = this_gfo_mp2_hist;
        fig = plt.figure();
        ax = fig.add_subplot(111);
        ax.hist(all_prn_mp2, bins=np.arange(-10, 10, 0.1));
        ax.grid(True);
        plt.xlabel('Value of MP2 [m]');
        plt.ylabel('# of epochs');
        units = "m";
        mp2_mean = "%.2f"%np.average(all_prn_mp2);
        mp2_std = "%.2f"%np.std(all_prn_mp2);
        mp2_rms = "%.2f"%np.sqrt(np.average(np.square(all_prn_mp2)));
        plt.xlim(-10,10);
        #plt.ylim(0,30000);
        plt.title('GFO_'+GFOsat+" (Mean: "+str(mp2_mean)+" $\pm$ "+str(mp2_std)+" "+units+", RMS: "+str(mp2_rms)+" "+units+")");
        ax.annotate ( 'Created: '+str ( commands.getoutput ( 'date' ) ) ,xy= ( 5, 5 ) , xycoords='figure points',fontsize=6);
	#plt.text(27,0,'Created on '+str(commands.getoutput('date')),horizontalalignment='right',verticalalignment='bottom',rotation='vertical');
        fig.savefig(plotname, dpi=300)
        plt.close(fig);

	this_gfo_sn1_hist ="GPS1A_"+dateProc+"_"+GFOsat+"_SN1.png";
        plotname = this_gfo_sn1_hist;
        fig = plt.figure();
        ax = fig.add_subplot(111);
        ax.hist(all_prn_sn1, bins=np.arange(0, 1000, 10));
        ax.grid(True);
        plt.xlabel('Value of SN1 [dB]');
        plt.ylabel('# of epochs');
        units = "dB";
        sn1_mean = "%.2f"%np.average(all_prn_sn1);
        sn1_std = "%.2f"%np.std(all_prn_sn1);
        sn1_rms = "%.2f"%np.sqrt(np.average(np.square(all_prn_sn1)));
        plt.xlim(0,1000);
        #plt.ylim(0,30000);
        plt.title('GFO_'+GFOsat+" (Mean: "+str(sn1_mean)+" $\pm$ "+str(sn1_std)+" "+units+", RMS: "+str(sn1_rms)+" "+units+")");
        ax.annotate ( 'Created: '+str ( commands.getoutput ( 'date' ) ) ,xy= ( 5, 5 ) , xycoords='figure points',fontsize=6);
	#plt.text(1010,0,'Created on '+str(commands.getoutput('date')),horizontalalignment='right',verticalalignment='bottom',rotation='vertical');
        fig.savefig(plotname, dpi=300)
        plt.close(fig);

	this_gfo_sn2_hist ="GPS1A_"+dateProc+"_"+GFOsat+"_SN2.png";
        plotname = this_gfo_sn2_hist;
        fig = plt.figure();
        ax = fig.add_subplot(111);
        ax.hist(all_prn_sn2, bins=np.arange(0, 1000, 10));
        ax.grid(True);
        plt.xlabel('Value of SN2 [dB]');
        plt.ylabel('# of epochs');
        units = "dB";
        sn2_mean = "%.2f"%np.average(all_prn_sn2);
        sn2_std = "%.2f"%np.std(all_prn_sn2);
        sn2_rms = "%.2f"%np.sqrt(np.average(np.square(all_prn_sn2)));
        plt.xlim(0,1000);
        #plt.ylim(0,30000);
        plt.title('GFO_'+GFOsat+" (Mean: "+str(sn2_mean)+" $\pm$ "+str(sn2_std)+" "+units+", RMS: "+str(sn2_rms)+" "+units+")");
        ax.annotate ( 'Created: '+str ( commands.getoutput ( 'date' ) ) ,xy= ( 5, 5 ) , xycoords='figure points',fontsize=6);        
	#plt.text(1010,0,'Created on '+str(commands.getoutput('date')),horizontalalignment='right',verticalalignment='bottom',rotation='vertical');
        fig.savefig(plotname, dpi=300)
        plt.close(fig);




gfoC_file.close();
gfoD_file.close();


daily_metrics_file.close();

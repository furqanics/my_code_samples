#!/usr/bin/python2.7

import os
import matplotlib
import numpy as np
matplotlib.use('Agg')
from matplotlib.pyplot import figure, show, rc, grid
import sys
import commands

satProc = sys.argv[1];
dateProc = sys.argv[2];

#fileProc = "G02_C_2018-04-15.azelsnr";
#plotname = "GPS1A_"+2018-04-15_C_G01_CASNR.png";
this_gfo_files = commands.getoutput("ls G??_C_"+dateProc+".azelsnr").split();
os.system("cat G??_C_"+dateProc+".azelsnr > tempthisGFOazel");

gfoid = satProc;
#gpsid = fileProc.split('_')[0];

#plotname = "GPS1A_"+fileProc.split('_')[2]+"_"+gfoid+"_"+gpsid 2018-04-15_C_G01_CASNR.png";

QTYS = ["L1SNR","L2SNR","CASNR"];
UNTS = ["dB","dB","dB"];

datafile = open("tempthisGFOazel",'r').readlines();
os.system("rm tempthisGFOazel");

az = [];
el = [];
l1snr = [];
l2snr = [];
casnr = [];

for l in datafile:
        az.append(float(l.split()[3]));
        el.append(float(l.split()[4]));
        l1snr.append(float(l.split()[5]));
        l2snr.append(float(l.split()[6]));
        casnr.append(float(l.split()[7]));
		
elneg = []
for i in el:
	elneg.append(90-i);
	
#fig = plt.figure();
#fig.set_size_inches(20, 10, forward=True)
plotname = "GPS1A_"+dateProc+"_"+gfoid+"_"+"L1SNR.png";
cmap = matplotlib.pyplot.cm.rainbow
norm = matplotlib.colors.Normalize(vmin=0, vmax=1600)
# radar green, solid grid lines
rc('grid', color='gray', linewidth=1, linestyle='-')
rc('xtick', labelsize=20)
rc('ytick', labelsize=20)
# force square figure and square axes looks better for polar, IMO
width, height = matplotlib.rcParams['figure.figsize']
size = min(width, height)
# make a square figure
fig = figure(figsize=(size*2.5, size*2.5))
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], polar=True)
ax.set_theta_zero_location("N")
ax.set_theta_direction(-1)
#ax.set_rlim(90, 0, 10)
ax.set_rmin(0);
ax.set_rmax(90.0);
#ax.set_yticks(np.arange(-45, 91, 15))
ax.scatter(az, elneg, s=3, color=cmap(norm(l1snr)))
ax.set_yticks(range              (0, 90, 10))    # (min int, max int, increment) 
ax.set_yticklabels(map(str, range(90, 0, -10)))
ax.set_rmax(90.0);
sm = matplotlib.pyplot.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
matplotlib.pyplot.colorbar(sm);
#ax.set_rmax(2.0)
grid(True)
#ax.text(225,120,'Created: '+str(commands.getoutput('date')),rotation="vertical",horizontalalignment="left",fontsize=6);
ax.annotate ( '[dB]' , xy= ( 665, 55 ) , xycoords='figure points',fontsize=25);
ax.annotate ( 'Created: '+str ( commands.getoutput ( 'date' ) ) ,xy= ( 80, 80 ) , xycoords='figure points',fontsize=8);
ax.set_title('GFO_'+gfoid+' ALL_PRN L1_SNR '+dateProc, fontsize=25)
fig.savefig(plotname);
matplotlib.pyplot.close(fig);

plotname = "GPS1A_"+dateProc+"_"+gfoid+"_"+"L2SNR.png";
cmap = matplotlib.pyplot.cm.rainbow
norm = matplotlib.colors.Normalize(vmin=0, vmax=1600)
# radar green, solid grid lines
rc('grid', color='gray', linewidth=1, linestyle='-')
rc('xtick', labelsize=20)
rc('ytick', labelsize=20)
# force square figure and square axes looks better for polar, IMO
width, height = matplotlib.rcParams['figure.figsize']
size = min(width, height)
# make a square figure
fig = figure(figsize=(size*2.5, size*2.5))
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], polar=True)
ax.set_theta_zero_location("N")
ax.set_theta_direction(-1)
#ax.set_rlim(90, 0, 10)
ax.set_rmin(0);
ax.set_rmax(90.0);
#ax.set_yticks(np.arange(-45, 91, 15))
ax.scatter(az, elneg, s=3, color=cmap(norm(l2snr)))
ax.set_yticks(range              (0, 90, 10))    # (min int, max int, increment) 
ax.set_yticklabels(map(str, range(90, 0, -10)))
ax.set_rmax(90.0);
sm = matplotlib.pyplot.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
matplotlib.pyplot.colorbar(sm);
#ax.set_rmax(2.0)
grid(True)
ax.annotate ( '[dB]' , xy= ( 665, 55 ) , xycoords='figure points',fontsize=25);
ax.annotate ( 'Created: '+str ( commands.getoutput ( 'date' ) ) ,xy= ( 80, 80 ) , xycoords='figure points',fontsize=8);
#ax.text(225,120,'Created: '+str(commands.getoutput('date')),rotation="vertical",horizontalalignment="left");
ax.set_title('GFO_'+gfoid+' ALL_PRN L2_SNR '+dateProc, fontsize=25)
fig.savefig(plotname);
matplotlib.pyplot.close(fig);

plotname = "GPS1A_"+dateProc+"_"+gfoid+"_"+"CASNR.png";
cmap = matplotlib.pyplot.cm.rainbow
norm = matplotlib.colors.Normalize(vmin=0, vmax=1600)
# radar green, solid grid lines
rc('grid', color='gray', linewidth=1, linestyle='-')
rc('xtick', labelsize=20)
rc('ytick', labelsize=20)
# force square figure and square axes looks better for polar, IMO
width, height = matplotlib.rcParams['figure.figsize']
size = min(width, height)
# make a square figure
fig = figure(figsize=(size*2.5, size*2.5))
ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], polar=True)
ax.set_theta_zero_location("N")
ax.set_theta_direction(-1)
#ax.set_rlim(90, 0, 10)
ax.set_rmin(0);
ax.set_rmax(90.0);
#ax.set_yticks(np.arange(-45, 91, 15))
ax.scatter(az, elneg, s=3, color=cmap(norm(casnr)))
ax.set_yticks(range              (0, 90, 10))    # (min int, max int, increment) 
ax.set_yticklabels(map(str, range(90, 0, -10)))
ax.set_rmax(90.0);
sm = matplotlib.pyplot.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
matplotlib.pyplot.colorbar(sm);
#ax.set_rmax(2.0)
grid(True)
#ax.text(225,120,'Created: '+str(commands.getoutput('date')),rotation="vertical",horizontalalignment="left");
ax.annotate ( '[dB]' , xy= ( 665, 55 ) , xycoords='figure points',fontsize=25);
ax.annotate ( 'Created: '+str ( commands.getoutput ( 'date' ) ) ,xy= ( 80, 80 ) , xycoords='figure points',fontsize=8);
ax.set_title('GFO_'+gfoid+' ALL_PRN CA_SNR '+dateProc, fontsize=25)
fig.savefig(plotname);
matplotlib.pyplot.close(fig);


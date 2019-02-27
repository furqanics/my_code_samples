#!/usr/bin/python

import os
import numpy
import sys

def glo_der(t,y, acc):
	gmWGS = float(398.60044e12); #m^3/s^2
	AE = float(6378136.0);
	OMEGA = float(7292115.e-11);
	C20= float(-1082.6257e-6);
	
	#Initial positions:
	rr1 = float(y[0]);
	rr2 = float(y[1]);
	rr3 = float(y[2]);
	
	#Initial velocities:
	vv1 = float(y[3]);
	vv2 = float(y[4]);
	vv3 = float(y[5]);
	
	rho = numpy.sqrt(rr1*rr1+rr2*rr2+rr3*rr3);
	
	t1 = -gmWGS/(rho*rho*rho);
	t2 = float(3.0)/float(2.0) * C20 * (gmWGS*AE*AE)/(rho*rho*rho*rho*rho);
	t3 = OMEGA * OMEGA;
	t4 = float(2.0) * OMEGA;
	z2 = rr3 * rr3;
	
	#Vector of derivatives:
	xyz_pv = numpy.zeros((6,1)); #Column vector
	xyz_pv[0] = vv1;
	xyz_pv[1] = vv2;
	xyz_pv[2] = vv3;
	xyz_pv[3] = (t1 + t2*(float(1.0)-float(5.0)*z2/(rho*rho)) + t3) * rr1 + t4*vv2 + float(acc[0]);
	xyz_pv[4] = (t1 + t2*(float(1.0)-float(5.0)*z2/(rho*rho)) + t3) * rr2 - t4*vv1 + float(acc[1]);
	xyz_pv[5] = (t1 + t2*(float(3.0)-float(5.0)*z2/(rho*rho))     ) * rr3          + float(acc[2]);
	return xyz_pv;

def rk4(xi, yi, dx, acc, glo_der):

	k1 = numpy.zeros((6,1));
	k2 = numpy.zeros((6,1));
	k3 = numpy.zeros((6,1));
	k4 = numpy.zeros((6,1));
	yf = numpy.zeros((6,1));
	
	k1 = glo_der(xi,                      yi,              acc) * float(dx);
	k2 = glo_der(xi+float(dx)/float(2.0), yi+k1/float(2.0),acc) * float(dx);
	k3 = glo_der(xi+float(dx)/float(2.0), yi+k2/float(2.0),acc) * float(dx);
	k4 = glo_der(xi+float(dx)           , yi+k3           ,acc) * float(dx);
	
	yf = yi + k1/float(6.0) + k2/float(3.0) + k3/float(3.0) + k4/float(6.0);
	
	return yf;
	
#Program inputs: satellite initial X, Y, Z [m], satellite initial Vx, Vy, Vz [m/s]


if __name__ == '__main__':

	X = float(sys.argv[1]);
	Y = float(sys.argv[2]);
	Z = float(sys.argv[3]);
	Vx = float(sys.argv[4]);
	Vy = float(sys.argv[5]);
	Vz = float(sys.argv[6]);
	Ax = float(sys.argv[7]);
	Ay = float(sys.argv[8]);
	Az = float(sys.argv[9]);
	

#R01 2017 07 02 00 15 00 1.189019531012e-05 0.000000000000e+00 0.000000000000e+00
#    -6.754042968750e+03-4.355239868164e-01 0.000000000000e+00 0.000000000000e+00
#     1.494270947266e+04-2.640743255615e+00 9.313225746155e-10 1.000000000000e+00
#    -1.953489208984e+04-1.867003440857e+00 1.862645149231e-09 0.000000000000e+00


#	X = float(-6.754042968750e+03)*1000; Vx = float(-4.355239868164e-01)*1000; Ax = float(0.000000000000e+00)*1000;
#	Y = float(1.494270947266e+04)*1000; Vy = float(-2.640743255615e+00)*1000; Ay = float(9.313225746155e-10)*1000;
#	Z = float(-1.953489208984e+04)*1000; Vz = float(-1.867003440857e+00)*1000; Az = float(1.862645149231e-09)*1000;
	
#	print "Welcome to GLONASS RK4 Interpolator";
	
	#Initial state is in yi
	yi = numpy.zeros((6,1));
	acc = [];
	
	yi[0,0] = X;
	yi[1,0] = Y;
	yi[2,0] = Z;
	yi[3,0] = Vx;
	yi[4,0] = Vy;
	yi[5,0] = Vz;
	
	acc.append(Ax);
	acc.append(Ay);
	acc.append(Az);
	
	#Final state is in yf
	yf = numpy.zeros((6,1));
	
	#By default, run the interpolator to bring back the positions -18 seconds
	dx = float(-1);
	xi = float(18);
	xf = float(0);
	nSteps = 18;
	
	for i in range(0,nSteps):
#		print "Iter.: "+str(i+1);
		#print "xi: ", xi-i;
		#print "yi: ", yi;
		yf = rk4(xi, yi, dx, acc, glo_der)
#		print "yf: ", yf/1000;
		yi = yf;
		xi = xi - 1;
		
	
	
	print yf[0,0], yf[1,0], yf[2,0];
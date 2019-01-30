from __future__ import division
import numpy as np
from math import *
import matplotlib.pyplot as plt

tmp = 0

def interpcurve(N,pX,pY):
	#equally spaced in arclength
	N=np.transpose(np.linspace(0,1,N))
	#print(N)

	#how many points will be uniformly interpolated?
	nt=N.size

	#number of points on the curve
	n=pX.size
	pxy=np.array((pX,pY)).T
	p1=pxy[0,:]
	pend=pxy[-1,:]
	last_segment= np.linalg.norm(np.subtract(p1,pend))
	epsilon= 10*np.finfo(float).eps

	#IF the two end points are not close enough lets close the curve
	if last_segment > epsilon*np.linalg.norm(np.amax(abs(pxy),axis=0)):
		pxy=np.vstack((pxy,p1))
		nt = nt + 1
	else:
		# pass
		print('Contour already closed')

	pt=np.zeros((nt,2))

	#Compute the chordal arclength of each segment.
	chordlen = (np.sum(np.diff(pxy,axis=0)**2,axis=1))**(1/2)
	#Normalize the arclengths to a unit total
	chordlen = chordlen/np.sum(chordlen)
	#cumulative arclength
	cumarc = np.append(0,np.cumsum(chordlen))

	tbins= np.digitize(N,cumarc) # bin index in which each N is in

	#catch any problems at the ends
	tbins[np.where(tbins<=0 | (N<=0))]=1
	tbins[np.where(tbins >= n | (N >= 1))] = n - 1      

	s = np.divide((N - cumarc[tbins]),chordlen[tbins-1])
	pt = pxy[tbins,:] + np.multiply((pxy[tbins,:] - pxy[tbins-1,:]),(np.vstack([s]*2)).T)

	return pt 


foil_origin = (0.0, 0.0)

def rotate(angle, O, P):
	x, y = P[0] - O[0], P[1] - O[1]
	x2 = cos(angle) * x - sin(angle) * y + O[0]
	y2 = sin(angle) * x - cos(angle) * y + O[1]

	return [x2, y2]


def shiftscale(points, shift=(0,0), scale=1.):
	'''
	shifts & scales points (np.array) 
	shift: tuple(x,y)
	scale: float
	'''
	global foil_origin
	O = foil_origin
	tmp = list(points)
	for p in tmp:
		p[0] -= O[0]
		p[1] -= O[1]
		p[0] *= scale
		p[1] *= scale
		p[0] += shift[0]
		p[1] += shift[1]
	return tmp


def NACA(serial, rot_angle, for_plotting = False, dX = 0.0125):
	# NACA-(A1,A2,A34)
	#serial = input("NACA WHAT?")
	#rot_angle = input("AOA?")
	global foil_origin

	#dX = 0.0125

	rot_angle = -radians(rot_angle)
	
	A1 = float(serial / 1000)
	A2 = float((serial / 100) % 10)
	A34 = float(serial % 100)

	T = A34 / 100
	X = 1.0
	C = 1.0
	M = A1 / 100

	P = A2 / 10
	O = foil_origin

	
	final = []

	while (X <= C):  # DO UPPER PART

		S1 = 0.2969 * ((X / C)**(0.5))
		S2 = (-0.126) * (X / C)
		S3 = (-0.3516) * (X / C)**2
		S4 = (0.2843) * (X / C)**3
		S5 = ((-0.1015) * ((X / C)**4))
		Yt = 5 * T * C * (S1 + S2 + S3 + S4 + S5)

		if 0 <= X < P * C:
			Yc = M * X / (P**2) * (2 * P - X / C)
			Tan_theta = (2 * M * (C * P - X)) / ((P**2) * C)

			Cos_theta = ((((P**2) * C)**2) / (((2 * M * (C * P - X))
											   ** 2) + (((P**2) * C)**2)))**(1.0 / 2)

			Sin_theta = (((2 * M * (C * P - X))**2) /
						 (((2 * M * (C * P - X))**2) + (((P**2) * C)**2)))**(1.0 / 2)
		else:
			Yc = M * (C - X) / ((1 - P)**2) * (1 + (X / C) - 2 * P)
			Tan_theta = (2 * M * (C * P - X)) / (((1 - P)**2) * C)

			Cos_theta = (((((1 - P)**2) * C)**2) / (((2 * M * (C * P - X))
													 ** 2) + ((((1 - P)**2) * C)**2)))**(1.0 / 2)

			Sin_theta = (((2 * M * (C * P - X))**2) / (((2 * M *
														 (C * P - X))**2) + ((((1 - P)**2) * C)**2)))**(1.0 / 2)

		X_U = X - Yt * Sin_theta
		Y_U = Yc + Yt * Cos_theta

		X_L = X + Yt * Sin_theta
		Y_L = Yc - Yt * Cos_theta

		coor = [X_U, -Y_U]
		if 0.999 < X/C < 1.001: 
			coor = [1., 0.]
		coor = rotate(angle=rot_angle, O=O, P=coor)
		final.append([round(coor[0], 5), round(coor[1], 5)])

		if (X < dX):
			break
		X -= dX

	X = 0
	while (X < C - dX):  # DO LOWER PART

		X += dX

		S1 = 0.2969 * ((X / C)**(0.5))
		S2 = (-0.126) * (X / C)
		S3 = (-0.3516) * (X / C)**2
		S4 = (0.2843) * (X / C)**3
		S5 = ((-0.1015) * ((X / C)**4))
		Yt = 5 * T * C * (S1 + S2 + S3 + S4 + S5)

		if 0 <= X < P * C:

			Yc = M * X / (P**2) * (2 * P - X / C)
			Tan_theta = (2 * M * (C * P - X)) / ((P**2) * C)

			Cos_theta = ((((P**2) * C)**2) / (((2 * M * (C * P - X))
											   ** 2) + (((P**2) * C)**2)))**(1.0 / 2)

			Sin_theta = (((2 * M * (C * P - X))**2) /
						 (((2 * M * (C * P - X))**2) + (((P**2) * C)**2)))**(1.0 / 2)

		else:

			Yc = M * (C - X) / ((1 - P)**2) * (1 + (X / C) - 2 * P)

			Tan_theta = (2 * M * (C * P - X)) / (((1 - P)**2) * C)

			Cos_theta = (((((1 - P)**2) * C)**2) / (((2 * M * (C * P - X))
													 ** 2) + ((((1 - P)**2) * C)**2)))**(1.0 / 2)

			Sin_theta = (((2 * M * (C * P - X))**2) / (((2 * M *
														 (C * P - X))**2) + ((((1 - P)**2) * C)**2)))**(1.0 / 2)

		X_U = X - Yt * Sin_theta
		Y_U = Yc + Yt * Cos_theta
		X_L = X + Yt * Sin_theta
		Y_L = Yc - Yt * Cos_theta

		coor = [X_L, -Y_L]
		if 0.999 < X/C < 1.001: 
			coor = [1., 0.]

		coor = rotate(angle=rot_angle, O=O, P=coor)
		final.append([round(coor[0], 5), round(coor[1], 5)])
		# print(round(coor[0],3),round(coor[1],3))
		# print(round(X_U,3),round(Y_U,5))
		# print(round(X_L,3),round(Y_L,5))
		# print(round(X,3),round(Yt,5))
		nparr = np.array(final)

	# print(final)
	shiftscale(nparr,shift = (-0.5,0))
	if for_plotting:
		return nparr.T
	return nparr

def dist(a,b):
	return sqrt((b[0] - a[0])**2 + (b[1]-a[1]) **2)

def evenNACA(number, angle, target_dist):
	global tmp
	#print("generating NACA {:_>4}".format(number))
	data = NACA(number,angle,True)
	N = len(data.T)
	#print "N = %f"%N
	even = interpcurve(N, data[0], data[1])
	d1 = dist(even[0],even[1])
	#print d1

	new_dX = 0.0125 * target_dist/d1
	#print "new_dX = ", new_dX
	new_data = NACA(number,angle,True,dX = new_dX)
	new_N = len(new_data.T)
	new_even = interpcurve(new_N, new_data[0], new_data[1])
	new_dist = dist(new_even[0],new_even[1])
	#print "new N = ", new_N 
	#print "target dist = %f, new dist = %f"%(target_dist, new_dist) 
	error = (new_dist/target_dist-1)*100
	if error > 1:
		print('NACA Serial: {}'.format(number))
		print('error rate: {} %'.format(error) )
		tmp+=1
	return new_even


def make_body_file(data, name):
	print( 'creating file: {}'.format(name) )
	file = open(name, 'w')
	file.write(str(len(data)))
	file.write('\n')
	for p in data:
		file.write('{} {}\n'.format(p[0],p[1]))
	file.close()


if __name__ == '__main__':
	global tmp

	# ====================================
	# NACA Airfoil:
	# 1st number = maximum camber height
	# 2st number = maximum camber at x = ?
	# 3,4th number = maximum thickness
	# ====================================


	#'''
	range1 = [i for i in range(9)]
	range34= [j for j in range(10,25)]
	before = NACA(7312,0,True)
	#after1 = evenNACA(2520,0,0.015)
	#print after1
	#make_body_file(after1,'body.001.inp')
	#after1 = after1.T
	plt.axis('equal')
	#plt.plot(data1[0],data1[1],'r-')
	#plt.plot(data2[0],data2[1],'r-')
	#plt.plot(after1[0],after1[1],'b+')
	plt.plot(before[0],before[1],'r-')
	#plt.plot(after2[0],after2[1],'b+')
	plt.show()
	#'''
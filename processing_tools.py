import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from functions import *

def time_zero (data, dtdx, t0=0.0):
	"""Place the initial time of your radargrams at t0 """

	dt, dx = dtdx
	iterations, taces = data.shape
	t = np.linspace(0, 1, iterations) * (iterations * dt)
	
	index = 0	
	while t0*1e-9 > t[index]: 
		index += 1
	
	return data[index:]
		
def user_gain (data, dtdx, gain='', start=0.0, stop=0.0):
	"""Add a user defined gain chosen beetween {'constant', 'linear', 'enxponential'} on data.
	start and stop define the number of traces you want modify.
	"""
	dt, dx = dtdx
	iterations, traces = data.shape
	t = np.linspace(0, 1, iterations) * (iterations * dt)
	
	width = stop - start
	
	if gain == 'constant':
		#check = input('Please enter a constant value : ')
		#constant = float(check)
		constant = 4
		fgain = [constant]*width
		
	elif gain == 'linear':
		check = input('Please enter gradient value: ')
		gradient = float(check)
		fgain = [gradient*1e9*t for time in t[start:stop]]

	elif gain == 'exponential':
		check = input('fgain = A * exp(B*t). Please enter A, B values: ')
		A = float(check[0])
		B = float(check[2])
		fgain = [A*np.exp(B*t) for time in t[start:stop]]
	else:
		fgain = [1]*width
		
	for model in range(traces):
		data[start:stop, model] *= np.array(fgain, dtype=data.dtype)

def velocity_analysis (data, dtdx, param, start, stop):
	"""
	Plot the radargram along with the hyperbolic function initialized by x0,t0,c and r.
	"""

	traces = data.shape[1]	
	dx = dtdx[1]
	x0, t0, v, r = param
	z0 = (t0 * v + 2 * r) / 2
	x = np.linspace(0, 1, traces) * (traces * dx)
	v *= 1e9
	
	hyperbol =  (2 / v) * (np.sqrt((x0-x[start:stop])**2 + z0**2) - r) 
	
	fig = plt.figure(num='radargrams', figsize=(20, 10), facecolor='w', edgecolor='w')
	plt.plot(x[start:stop], hyperbol)
	plot_radargrams(data, dtdx)
	
def stolt_migration(data, param):
	
	
	dx = param[1]
	dt = param[0]*1e9
	fs = 1/dt
	c = param[2]


	nt0, nx0 = data.shape
	
	t = np.linspace(0,nt0*dt,nt0) 
	x = np.linspace(0,nx0*dx,nx0)
	
	nt = 2 * nextpower(nt0)
	nx = 2 * nx0

	ERMv = c / 2 # One Emiter-Receiver scenario

	fftRF = np.fft.fftshift(np.fft.fft2(data, s=(nt,nx)))

	f = np.linspace(-nt/2, nt/2-1, nt) * fs / nt 
	kx = np.linspace(-nx/2,nx/2-1, nx) / dx / nx
	
	KX, F = np.meshgrid(kx, f)
	fkz = ERMv*np.sign(F)*np.sqrt(KX**2 + F**2/ERMv**2)
	
	# interpolate onto the new grid
	fftRF = griddata((KX.ravel(), F.ravel()), fftRF.ravel(), (KX, fkz), method='nearest')
	
	
	migRF = np.fft.ifft2(np.fft.ifftshift(fftRF))
	migRF = migRF[1:nt0,1:nx0]
	
	fig = plt.figure(num='', figsize=(20, 10), facecolor='w', edgecolor='w')
	plt.imshow(abs(migRF), extent=[np.amin(x), np.amax(x), np.amax(t), np.amin(t)], interpolation='nearest', aspect='auto', cmap='seismic', vmin=-np.amax(abs(migRF)), vmax=np.amax(abs(migRF)))
	plt.show()
	

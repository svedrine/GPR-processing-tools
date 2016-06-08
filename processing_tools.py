import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from functions import *


def dewow_filtering(data, dx_dt, time_window):
	
	dx, dt = dx_dt
	iterations, traces = data.shape
	t = np.linspace(0, 1, iterations) * (dt * iterations)
	width = int(time_window / dt)
	
	for trace in traces:
		index = 0
		
	
def time_zero (data, dx_dt, t0=0.0):
	"""Reshape start time of your radargrams at t0 """

	dx, dt = dx_dt
	iterations, traces = data.shape
	t = np.linspace(0, 1, iterations) * (iterations * dt)
	index = int(t0 / dt)

	return data[index:]
		
def user_gain (data, dx_dt, gain='', start=0.0, stop=0.0):
	"""Add a user defined gain chosen beetween {'constant', 'linear', 'enxponential'} on data.
	start and stop define the number of traces you want to modify.
	"""
	dx, dt = dx_dt
	iterations, traces = data.shape
	t = np.linspace(0, 1, iterations) * (iterations * dt)
	
	start = int(start / dt)
	stop = int(stop / dt)
	
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

def velocity_analysis (data, dx_dt, param, start, stop):
	"""
	Plot the radargram along with the hyperbolic function initialized by x0,t0,c and r.
	"""

	traces = data.shape[1]	
	dx = dx_dt[1]
	x0, t0, v, r = param
	z0 = (t0 * v + 2 * r) / 2
	x = np.linspace(0, 1, traces) * (traces * dx)
	
	hyperbol =  (2 / v) * (np.sqrt((x0-x[start:stop])**2 + z0**2) - r) 
	
	plt.plot(x[start:stop], hyperbol)
	plot_radargrams(data, dx_dt, 'velocity_analysis')

def stolt_migration(data, dx_dt, c):
	
	# parameters
	dx, dt = dx_dt
	fs = 1/dt
	eps = 1e-4
	nt0, nx0 = data.shape
	t = np.linspace(0,nt0*dt,nt0) 
	x = np.linspace(0,nx0*dx,nx0)
	
	# Zero-padding
	nt = 2 * nextpower(nt0)
	nx = 2 * nx0
	
	# One Emiter-Receiver scenario
	ERMv = c / 2
	
	# FFT & shift 
	fftdata = np.fft.fftshift(np.fft.fft2(data, s=(nt,nx)))
	
    
	# Linear interpolation
	f = np.linspace(-nt/2, nt/2-1, nt) * fs / nt 
	kx = np.linspace(-nx/2,nx/2-1, nx) / dx / nx
	kx, f = np.meshgrid(kx, f)
	
	# Remove evanescent parts
	evanescent = (abs(f)  / abs(kx+eps) > c).astype(int)
	fftdata *= evanescent
	
	# Remapping
	fkz = ERMv*np.sign(f)*np.sqrt(kx**2 + f**2/ERMv**2)
	fftdata = griddata((kx.ravel(), f.ravel()), fftdata.ravel(), (kx, fkz), method='nearest')
	
	# Jacombien
	fftdata *= f / (np.sqrt(kx**2 + f**2/ERMv**2)+eps)
	
	# IFFT & Migrated RF
	mig = np.fft.ifft2(np.fft.ifftshift(fftdata))
	mig = mig[1:nt0,1:nx0]
	
	z = t*c / 2
	
	fig = plt.figure(num='migrated_data', figsize=(20, 10), facecolor='w', edgecolor='w')
	plt.imshow(abs(mig), extent=[np.amin(x), np.amax(x), np.amax(z), np.amin(z)], interpolation='nearest', aspect='auto', cmap='seismic', vmin=-np.amax(abs(mig)), vmax=np.amax(abs(mig)))
	plt.xlabel('Distance [m]')
	plt.ylabel('Depth [m]')
	
	plt.show()
	

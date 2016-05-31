import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d
from functions import nextpower, ascii_to_nparray, plot_radargrams

def time_zero (path, filename, t0=0.0):
	"""Place the initial time of your radargrams at t0
	create a new ascii file in path"""
	
	data = ascii_to_nparray(path, filename)
	data = data.T
	new_name = filename.replace('.ASC', '_zero.ASC')
	f = open('%s%s' %(path,new_name), 'wb')
	
	dt = 0.018654
	iterations, traces = data.shape
	t = np.linspace(0, 1, iterations) * (iterations * dt)
	index = 0	
	while t0 > t[index]: index += 1
	print(t0)
	np.savetxt(f,data[index:].astype(int),fmt='%d', delimiter = '  ', newline='\n')
		
	
def user_gain (path, filename, gain='', start=0.0, stop=0.0):
	"""Add a user defined gain chosen beetween {'constant', 'linear', 'enxponential'} on data.
	start and stop define the number of traces you want modify.
	create a new ascii file gained in path.
	"""
	data = ascii_to_nparray(path, filename)
	data = data.T
	new_name = filename.replace('.ASC', '_%s.ASC' % gain)
	f = open('%s%s' %(path,new_name), 'wb')
	
	# init essential items
	dt = 0.018654
	dx = 0.00667
	iterations, traces = data.shape
	t = np.linspace(0, 1, iterations) * (iterations * dt)
	
	width = stop - start
	
	if gain == 'constant':
		check = input('Please enter a constant value : ')
		constant = float(check)
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
	
	np.savetxt(f,data.astype(int),fmt='%d', delimiter = '  ', newline='\n')
		

def velocity_analysis (data, x0, t0, v, r, start, stop):
	"""
	Plot the radargram along with the hyperbolic function initialized by x0,t0,c and r.
	"""
	data = data.T
	traces = data.shape[1]
	
	dx = 0.00667 
	z0 = (t0 * v + 2 * r) / 2
	x = np.linspace(0, 1, traces) * (traces * dx)
	hyperbol =  (2 / v) * (np.sqrt((x0-x[start:stop])**2 + z0**2) - r) 
	
	
	fig = plt.figure(num='radargrams', figsize=(20, 10), facecolor='w', edgecolor='w')
	plt.plot(x[start:stop], hyperbol)
	plot_radargrams(data)
	
	
def stolt_migration (data):
	
	data = data.T
	dt = 0.018654
	dx = 0.00667
	fs = 1/dt
	c = 0.25

	t = np.linspace(0,1024*dt,1024) 
	x = np.linspace(0,610*dx,610)

	nt0, nx0 = data.shape

	nt = 2 * nextpower(nt0)
	nx = 2 * nx0

	ERMv = c 

	fftRF = np.fft.fftshift(np.fft.fft2(data, s=(nt,nx)))

	f = np.linspace(-nt/2, nt/2-1, nt) * fs / nt 
	kx = np.linspace(-nx/2,nx/2-1, nx) / dx / nx 
	
	real = interp2d(kx, f, fftRF.real, 'linear')
	imag = interp2d(kx, f, fftRF.imag, 'linear')	
	
	kx, f = np.meshgrid(kx, f)
	
	fkz = ERMv*np.sign(f)*np.sqrt(kx**2 + f**2/ERMv**2)
	print(fkz)
	"""for row in range(nx):
		fftRF[row,:].real = real(kx[0,:], fkz[row,:])[row,:]
		fftRF[row,:].imag = imag(kx[0,:], fkz[row,:])[row,:]"""
	
	
	
	migRF = np.fft.ifft2(np.fft.ifftshift(fftRF))
	migRF = migRF[1:nt0,1:nx0]
	
	fig = plt.figure(num='', figsize=(20, 10), facecolor='w', edgecolor='w')
	plt.imshow(abs(fkz), extent=[np.amin(x), np.amax(x), np.amax(t), 0], interpolation='nearest', aspect='auto', cmap='seismic', vmin=-np.amax(abs(fkz)), vmax=np.amax(abs(fkz)))
	plt.show()
	
	return 0.0


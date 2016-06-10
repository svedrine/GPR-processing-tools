import numpy as np
import matplotlib.pyplot as plt
import h5py
import binascii


def plot_radargram (data, dx_dt, title):
	
	dx, dt = dx_dt
	iterations, traces = data.shape
	t = np.linspace(0, 1, iterations) * (iterations * dt)
	x = np.linspace(0, 1, traces) * (traces * dx)	
	
	fig = plt.figure(num=title, figsize=(20, 10), facecolor='w', edgecolor='w')
	plt.imshow(data, extent=[np.amin(x), np.amax(x), np.amax(t), np.amin(t)], interpolation='nearest', aspect='auto', cmap='seismic', vmin=-np.amax(abs(data)), vmax=np.amax(abs(data)))
	plt.xlabel('Distance [m]')
	plt.ylabel('Two-way travel time [ns]')
	
	plt.show()
	
def nextpower (n, base=2.0):
    """Return the next integral power of two greater than the given number.
    Specifically, return m such that
        m >= n
        m == 2**x
    where x is an integer. Use base argument to specify a base other than 2.
    This is useful for ensuring fast FFT sizes.
    """
    x = base**np.ceil(np.log(n) / np.log(base))
    if type(n) == np.ndarray:
        return np.asarray (x, dtype=int)
    else:
        return int (x)

def read_ascii (path, filename):
	"""
	Comments
	"""
	f = open('%s%s' % (path, filename), 'r')
	mylist = f.readlines()
	while '\n' in mylist: mylist.remove('\n')
	
	# Format data
	modelruns = np.shape(mylist)[0]
	data = []
	for row in range(modelruns):
		data.append(list(map(int, mylist[row].split())))
	f.close()
	
	dx_dt = (0.0185, 0.00667)
	
	return np.array(data).T, dx_dt

def read_hdf5 (path, filename):
	""" Convert an h5py data file into ascii file. """
	
	f = h5py.File('%s%s' % (path, filename), 'r')
	path = '/rxs/rx1/'
	modelruns = f.attrs['Modelruns']
	iterations = f.attrs['Iterations'] 
	dt = f.attrs['dt']*1e9
	positions = f.attrs['Positions'][:,0,0]
	dx = np.diff(positions)[0]
	data = np.ones((iterations, modelruns))
	for model in range(modelruns):
		data[:,model] = f['%s%s' % (path, 'Ez')][:,model]
	traces = modelruns	
	dx_dt = (iterations, traces, dx, dt)
	
	return data, dx_dt

def read_rd3 (path, filename):
	
	""" a modifier """
	f = open('%s%s' % (path, filename), 'rb').read()
	f.decode("latin-1")
	
	data = np.array(np.fromstring(f, dtype= 'int16'))
	iterations, traces = data.shape
	traces = int(np.shape(data)[0] / iterations)
	data = data.reshape((traces, iterations)
	
	return data.T
	
def read_header(path, filename):
	
	""" return dt (ns) and dx (m) """
	
	f = open('%s%s' %(path, filename), 'r')
	mylist = f.read().split('\n')
	
	iterations = int(mylist[0].split(':')[1])
	frequency = float(mylist[1].split(':')[1])*1e-3
	traces = int(mylist[22].split(':')[1])
	position = float(mylist[23].split(':')[1])
	
	dt = 1 / frequency
	dx = position / traces
	
	return (dx, dt)
	







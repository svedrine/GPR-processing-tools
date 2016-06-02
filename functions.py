import numpy as np
import matplotlib.pyplot as plt
import h5py

def plot_radargrams (data, params):
	
	dx, dt, traces, iterations = params
	t = np.linspace(0, 1, iterations) * (iterations * dt)
	x = np.linspace(0, 1, traces) * (traces * dx)	
	
	plt.imshow(data, extent=[np.amin(x), np.amax(x), np.amax(t), np.amin(t)], interpolation='nearest', aspect='auto', cmap='seismic', vmin=-np.amax(abs(data)), vmax=np.amax(abs(data)))
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

def ascii_to_nparray (path, filename):
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
	
	return np.array(data)


def convert_hdf5 (path, filename):
	""" Convert an h5py data file into ascii file. """
	
	f = h5py.File('%s%s' % (path, filename), 'r')
	path = '/rxs/rx1/'
	modelruns = f.attrs['Modelruns']
	iterations = f.attrs['Iterations'] 
	dt = f.attrs['dt']
	positions = f.attrs['Positions'][:,0,0]
	dx = np.diff(positions)[0]
	data = np.ones((iterations, modelruns))
	for model in range(modelruns):
		data[:,model] = f['%s%s' % (path, 'Ez')][:,model]
	params = (dx, dt, modelruns, iterations) 
	
	return data, params


import numpy as np

def plot_radargrams (data):
	
	iterations, traces = data.shape
	dt = 0.018654
	dx = 0.00667 
	t = np.linspace(0, 1, iterations) * (iterations * dt)
	positions = np.linspace(0, 1, traces) * (traces * dx)
	
	plt.imshow(data, extent=[np.amin(positions), np.amax(positions), np.amax(t), 0], interpolation='nearest', aspect='auto', cmap='seismic', vmin=-np.amax(abs(data)), vmax=np.amax(abs(data)))
	
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


def convert_hdf5 (filename):
	""" Convert an hf5py data file into ascii file.
	"""
	return 0.0


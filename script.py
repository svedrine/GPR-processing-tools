from processing_tools import *
from functions import *

# import and read 
path ='/home/svedrine/simon_models/'
filename = 'Bscan_merged.out' #'F119_JP1X.a.ASC'
data, dtdx =  read_hdf5(path, filename) #read_ascii(path, filename)
dt, dx = dtdx


data_zero = time_zero(data, dtdx, t0 = 0.5)
#user_gain(data_zero, dtdx, gain='constant', start = 150, stop = 500)

#param = (1.45, 5.5, 0.16, 0.1)
#velocity_analysis(data_gain, dtdx, param, start=5, stop=25)

param  = (dt, dx, 0.15)
stolt_migration(data_zero, param)










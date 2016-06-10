from processing_tools import *
from functions import *


# import and read 
path ='/home/svedrine/Desktop/160607-Merantaise/'
headername = 'Prof2.rad' #'Bscan_merged.out'
dataname = 'Prof2.rd3'

"""data, dtdx = read_ascii(path, filename) #read_hdf5(path, filename)
dt, dx = dtdx


data_zero = time_zero(data, dtdx, t0 = 0.5)
user_gain(data_zero, dtdx, gain='constant', start = 150, stop = 400)

param = (1.45, 5.5, 0.16, 0.1)
velocity_analysis(data_zero, dtdx, param, start=5, stop=25)

param  = (dt, dx, 0.195)
stolt_migration(data_zero, param)
plot_radargrams(data_zero, dtdx)"""

dx_dt = read_header(path, headername)
data = read_rd3(path, dataname)
dc_data = dewow_filtering(data, dx_dt, 25)




plot_radargram(dc_data, dx_dt, 'Allo')

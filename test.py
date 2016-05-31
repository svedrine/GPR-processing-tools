from processing_tools import *

path ='/home/svedrine/simon_data/'
filename = 'F119_JP1X.a.ASC'
data = ascii_to_nparray(path, filename)

stolt_migration(data)


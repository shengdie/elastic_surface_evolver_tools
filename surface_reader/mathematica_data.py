import numpy as np
import h5py

def import_hdf5(filepath):
    """import data from hdf5 file, convert to numpy array.
    
    Arguments:
        filepath {str} -- hdf5 file path
    
    Returns:
        array -- data read from hdf5 file
    """

    f = h5py.File(filepath, 'r')
    a_group_key = list(f.keys())[0]
    data = list(f[a_group_key])
    return np.array(data)

def data3d_for_plotly(data):
    x = data[:,:,0][:,0]
    y = data[:,:,1][0]
    z = np.transpose(data[:,:,2])
    return (x,y,z)


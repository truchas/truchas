import matplotlib.pyplot as plt
import numpy as np
import os
import h5py

## results from ghia et al for Re=100 (position, velocity)
ghia_u = (
    np.array([0.0, 0.0547, 0.0625, 0.0703, 0.1016, 0.1719, 0.2813,
              0.4531, 0.5, 0.6172, 0.7344, 0.8516, 0.9531,
              0.9609, 0.9688, 0.9766, 1.0]),
    np.array([0.0, -0.03717, -0.04192, -0.04775, -0.06434,
             -0.10150, -0.15662, -0.21090, -0.20581, -0.13641,
              0.00332, 0.23151, 0.68717, 0.73722, 0.78871,
              0.84123, 1.0]))
ghia_v = (
    np.array([0.0, 0.0625, 0.0703, 0.0781, 0.0938, 0.1563, 0.2266,
              0.2344, 0.5, 0.8047, 0.8594, 0.9063, 0.9453,
              0.9531, 0.9609, 0.9688, 1.0]),
    np.array([0.0, 0.09233, 0.10091, 0.10890, 0.12317, 0.16077,
              0.17507, 0.17527, 0.05454, -0.24533, -0.22445,
              -0.16914, -0.10313, -0.08864, -0.07391, -0.05906, 0.0]))

def truchas_output_from_inp(inp):
    """ generate name of truchas output for a given inp file """
    return os.path.join("_".join((inp, "output")),
                        ".".join((inp, "h5")))

def pbf_file(visit_dir):
    return os.path.join(visit_dir, "Level_0", "Cell_D_00000")

def read_pbf_scalar(visit_dir, nx=33):
    fname = pbf_file(visit_dir)
    with open(fname, 'rb') as f:
        double_bytes = nx*nx*8
        file_bytes = os.stat(fname).st_size
        f.seek(file_bytes-double_bytes)
        data = np.fromfile(f, count=double_bytes//8, dtype=np.float64)
        return data.reshape((nx, nx)).T


def read_pbf_vector(visit_dir, nx=33):
    fname = pbf_file(visit_dir)
    with open(fname, 'rb') as f:
        double_bytes = nx*nx*8*3
        file_bytes = os.stat(fname).st_size
        f.seek(file_bytes-double_bytes)
        ndoubles = nx*nx
        v1 = np.fromfile(f, count=ndoubles, dtype=np.float64)
        v2 = np.fromfile(f, count=ndoubles, dtype=np.float64)
        v3 = np.fromfile(f, count=ndoubles, dtype=np.float64)
        return (v1.reshape((nx, nx)).T,
                v2.reshape((nx, nx)).T,
                v3.reshape((nx, nx)).T)


def sorted_index_from_centroid(cent, f=1000):
    c = [ tuple(x)+(i,) for i, x in enumerate(cent) ]
    s = sorted(c, key=lambda x : tuple(int(f*x) for x in x))
    return [x[-1] for x in s]


def truchas_nseries(inp):
    truchas_h5_file = truchas_output_from_inp(inp)
    h5_series = h5py.File(truchas_h5_file, 'r')['Simulations/MAIN/Series Data']
    return len(list(h5_series.keys()))


def truchas_data_group(truchas_h5_file, series_index):
    h5_series = h5py.File(truchas_h5_file, 'r')['Simulations/MAIN/Series Data']
    group_name = sorted(list(h5_series.keys()),
                        key=lambda x: int(x.split()[-1]))[series_index]
    return h5_series[group_name]


def truchas_sorted_centroid(truchas_group, nx):
    raw_centroid = truchas_group['CENTROID'][:,:2] # only need x, y
    idx = sorted_index_from_centroid(raw_centroid)
    xyz = raw_centroid[idx].reshape((nx, nx, -1))
    return idx, xyz


def read_truchas_scalar(truchas_h5_file, var_name, series_index, nx=33):
    d = truchas_data_group(truchas_h5_file, series_index)

    idx, xyz = truchas_sorted_centroid(d, nx)

    raw_scalar = d[var_name][:]
    scalar = raw_scalar[idx].reshape((nx, nx))
    truchas_index = np.array(idx)+1
    return scalar, xyz, truchas_index.reshape((nx, nx))


def read_truchas_vector(truchas_h5_file, var_name, series_index, nx=33):
    d = truchas_data_group(truchas_h5_file, series_index)

    idx, xyz = truchas_sorted_centroid(d, nx)

    raw_vec = d[var_name][:,:2] # only 2d comparisons for now
    vec = raw_vec[idx].reshape((nx, nx, -1))
    v1 = vec[:,:,0]
    v2 = vec[:,:,1]

    truchas_index = np.array(idx)+1
    return v1, v2, xyz, truchas_index.reshape((nx, nx))


def compare_velocity_(inp, nx, series_index):
    """ compare velocities at single index """
    h5_new = truchas_output_from_inp(inp)
    t_new = dict(zip(('u', 'v', 'xyz', 'idx'),
                     read_truchas_vector(h5_new, 'Z_VC', series_index, nx)))
    h5_old = truchas_output_from_inp(inp + "-old")
    t_old = dict(zip(('u', 'v', 'xyz', 'idx'),
                     read_truchas_vector(h5_old, 'Z_VC', series_index, nx)))
    return t_new, t_old


def compare_pressure_(inp, nx, series_index):
    """ compare velocities at single index """
    h5_new = truchas_output_from_inp(inp)
    t_new = dict(zip(('p', 'xyz', 'idx'),
                     read_truchas_scalar(h5_new, 'Z_P', series_index, nx)))
    h5_old = truchas_output_from_inp(inp + "-old")
    t_old = dict(zip(('p', 'xyz', 'idx'),
                     read_truchas_scalar(h5_old, 'Z_P', series_index, nx)))
    return t_new, t_old

def compare_func_(f, inp, nx, series_index):
    if series_index is not None:
        return f(inp, nx, series_index)
    else:
        nseries = truchas_nseries(inp)
        res = []
        for i in range(nseries):
            res.append(f(inp, nx, i))
        return res


def compare_velocity(inp, nx, series_index=None):
    """ returns dict of velocity comparisons """
    return compare_func_(compare_velocity_, inp, nx, series_index)


def compare_pressure(inp, nx, series_index=None):
    """ returns dict of velocity comparisons """
    return compare_func_(compare_pressure_, inp, nx, series_index)

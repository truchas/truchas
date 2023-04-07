#!/usr/bin/env python3

import numpy as np

import truchas


def step(x, w=0.01): return 0.5 * (1 + np.tanh(x / w))
#def step(x, w=0): return np.where(x < 0, 0, 1)
def centered_box(sizes): return np.array([[-s / 2, s / 2] for s in sizes])


def E_analytic(x, z):
    # problem parameters
    f = 2.45e9 # [Hz] driving frequency
    d = 0.2 # [m] waveguide half-length
    a = 3.4 * 0.0254 # [m] WR340 waveguide dimension, long side
    b = 1.7 * 0.0254 # waveguide dimension, short side
    epsr = 4 + 0.1j # relative permittivity
    p = 100.0 # waveguide power

    # physical constants
    c = 299792458 # [m/s]
    mu0 = 4 * np.pi * 1e-7 # [H/m]

    # derived quantities
    omega = 2 * np.pi * f
    k0 = omega / c
    h0 = np.sqrt(k0**2 - np.pi**2 / a**2) # TE01 propagation constant in vacuum
    h1 = np.sqrt(epsr * k0**2 - np.pi**2 / a**2) # TE01 propagation constant in dielectric
    r01 = (h0 - h1) / (h0 + h1) # Reflection/transmission coefficients at each interface, here from reflection from 1 into 0
    r10 = -r01
    t01 = 2 * h0 / (h0 + h1)
    t10 = 2 * h1 / (h0 + h1)
    r = r01 - t01 * t10 * np.exp(2j*h1*d) / (1 + r10 * np.exp(2j*h1*d)) # total reflection coefficient
    A = (1 + r) / (1 - np.exp(2j*h1*d))  # amplitude of mode propagating into +z-dir in dielectric
    B = (1 + r) / (1 - np.exp(-2j*h1*d)) # amplitude of mode propagating into -z-dir in dielectric
    E0 = np.sqrt((4*p/(a*b))*(k0*mu0*c/h0)) # incoming amplitude

    print(f"The absolute value of the reflection coefficient is {np.abs(r)}.")

    Eyl = E0 * np.cos(np.pi * x / a) * (np.exp(1j * h0 * z) + r * np.exp(-1j * h0 * z))
    Eyr = E0 * np.cos(np.pi * x / a) * (A * np.exp(1j * h1 * z) + B * np.exp(-1j * h1 * z))
    Ey = step(-z) * Eyl + step(z) * Eyr

    return np.abs(Ey)


def project_data(xin, fin):
    """Project data to a uniform Cartesian mesh for plotting with matplotlib."""
    import scipy as sp

    a = 3.4 * 0.0254
    xlim = centered_box([a, a/2, 0.4])
    dx = np.array([a / 9, a / 2 / 4, 0.4 / 40]) / 10
    x = np.mgrid[xlim[0,0]:xlim[0,1]+dx[0]:dx[0],
                 xlim[1,0]:xlim[1,1]+dx[1]:dx[1],
                 xlim[2,0]:xlim[2,1]+dx[2]:dx[2]]

    xr = np.array([q.ravel() for q in x]).T
    f = sp.interpolate.griddata(xin, fin, xr, method='nearest').reshape(x[0].shape)

    return x, f


def reduce_data(x3, f3):
    """Reduce from 3D to 2D using an average in the y direction."""
    x2 = np.array([x3[0][:,0,:], x3[2][:,0,:]])
    f2 = np.average(f3, axis=1)
    #f2 = f3[:,0,:]
    return x2, f2


def plot_results(x0, Et0):
    import matplotlib as mpl
    import matplotlib.pyplot as plt

    # get data
    Eg0 = E_analytic(x0[:,0], x0[:,2])
    x3D, Et3D = project_data(x0, Et0)
    _, Eg3D = project_data(x0, Eg0)
    x, Et = reduce_data(x3D, Et3D)
    _, Eg = reduce_data(x3D, Eg3D)
    #Eg = E_analytic(x[0], x[1])

    # plot
    nrows = 3; ncols = 1
    aspect_ratio = 4/3 * ncols / nrows
    height = 6 * nrows
    fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(aspect_ratio*height, height))

    fig.suptitle("Waveguide")

    for axi, f, label in ((ax[0], Et, "|E| Truchas"),
                          (ax[1], Eg, "|E| Analytic"),
                          (ax[2], Et - Eg, "|E| difference")):
        axi.set_xlabel("x")
        axi.set_xlabel("x")

        p = axi.pcolormesh(x[0], x[1], f)
        #p = axi.contourf(x[0], x[1], f, levels=20)

        cbar = fig.colorbar(p)
        cbar.ax.set_ylabel(label)

    fig.tight_layout()
    fig.savefig("wg.png", dpi=200)


def run_test(tenv):
    nfail = 0
    stdout, output_base = tenv.truchas(4, "wg1.inp")
    golden_base = tenv.output("wg1_golden/wg1.h5")
    output = output_base.em_data()
    golden = golden_base.em_data()

    x = output.centroids()
    Ey_gold = E_analytic(x[:,0], x[:,2])
    Emax = np.max(np.abs(Ey_gold)) # reference for error measurement

    test = output.field("|E|")
    gold = golden.field("|E|")
    nfail += truchas.compare_max(test, gold, 1e-9*Emax, "E-golden", 0.0)

    Ey_test = test[:,1]
    nfail += truchas.compare_max(Ey_test, Ey_gold, 0.15*Emax, "Ey-analytic", 0.0)
    nfail += truchas.compare_max(test[:,0], 0.0,  0.15*Emax, "Ex-analytic", 0.0)
    nfail += truchas.compare_max(test[:,2], 0.0, 0.15*Emax, "Ez-analytic", 0.0)

    #plot_results(x, Ey_test)

    truchas.report_summary(nfail)
    return nfail


if __name__=="__main__":
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0

#!/usr/bin/env python3

import time
import os
import multiprocessing

import numpy as np
import numpy.linalg as npla
import scipy as sp
import scipy.special as spsp
import scipy.constants
import matplotlib as mpl
import matplotlib.pyplot as plt

import truchas


def run_truchas(input_parameters, fdme=False):
    input_parameters["fdme"] = "T" if fdme else "F"
    input_file = "fdme.inp" if fdme else "tdme.inp"
    tenv.generate_input_deck(input_parameters, "template-em.inp", input_file)
    t_elapsed = -time.time()
    stdout, tdata = tenv.truchas(1, input_file)
    t_elapsed += time.time()
    solver = "FDME" if fdme else "TDME"
    print(f"{solver} elapsed {t_elapsed:.2f} seconds.")
    return tdata


def get_heat_source(tdata):
    sid = tdata.num_series()
    xc = tdata.centroids()
    r_truchas = np.sqrt(xc[:,0]**2 + xc[:,1]**2)
    Q_truchas = tdata.field(sid, "Joule_P")
    return r_truchas, Q_truchas


def get_Emag(tdata):
    temdata = tdata.em_data()
    xc = temdata.centroids()
    r = np.sqrt(xc[:,0]**2 + xc[:,1]**2)
    Emagv = temdata.field("|E|")
    Emag = npla.norm(Emagv, axis=1)
    return r, Emag


def heat_analytic(H_bc, source_frequency, sigma, R,
                  electric_susceptibility, magnetic_susceptibility):
    r = np.linspace(0, R, 100)
    Q, Emag = analytic_soln(r, H_bc, source_frequency, sigma, R,
                            electric_susceptibility, magnetic_susceptibility,
                            verbose=True)
    return r, Q, Emag


def analytic_soln(r, H_bc, source_frequency, sigma, R,
                  electric_susceptibility, magnetic_susceptibility,
                  verbose=False):
    mu = sp.constants.mu_0 * (1 + magnetic_susceptibility)
    epsilon = sp.constants.epsilon_0 * (1 + electric_susceptibility)
    omega = 2 * np.pi * source_frequency
    B_bc = mu * H_bc

    c = mu * (1j * omega * sigma - omega**2 * epsilon)
    sqrt_c = np.sqrt(c)
    A = B_bc / spsp.jv(0, sqrt_c * R)
    B = A * spsp.jv(0, sqrt_c * r)
    E = -1j * omega * A / sqrt_c * spsp.jvp(0, sqrt_c * r)
    E2 = np.abs(E)**2
    Q = (sigma + omega * np.imag(epsilon)) * np.abs(E)**2 / 2
    Emag = np.sqrt(E2)

    if verbose:
        print(f"mu: {mu}")
        print(f"epsilon: {epsilon}")
        print(f"sigma: {sigma}")
        print(c)
        print(sqrt_c)
        print(B_bc, B[-1])

    return Q, Emag


def plot_results(r, Q, Emag, tdata, title):
    r_fd, Q_fd = get_heat_source(tdata[0])
    if len(tdata) > 1: r_td, Q_td = get_heat_source(tdata[1])
    r_fd_em, Emag_fd = get_Emag(tdata[0])

    nrows = 1; ncols = 2
    aspect_ratio = 4/3 * ncols / nrows
    height = 6 * nrows
    fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(aspect_ratio*height, height))
    fig.suptitle(title)

    ax[0].set_xlabel("r")
    ax[0].set_ylabel("Q")
    if len(tdata) > 1: ax[0].semilogy(r_td, Q_td, 'x', label="TDME")
    ax[0].semilogy(r_fd, Q_fd, '+', label="FDME")
    ax[0].semilogy(r, Q, '-', label="analytic")
    ax[0].legend(loc="upper left")

    ax[1].set_xlabel("r")
    ax[1].set_ylabel("|E|")
    #if len(tdata) > 1: ax[1].semilogy(r_td, Q_td, 'x', label="TDME")
    ax[1].semilogy(r_fd_em, Emag_fd, '+', label="FDME")
    ax[1].semilogy(r, Emag, '-', label="analytic")

    fig.tight_layout()
    fig.savefig("em.png", dpi=200)


def test(tenv, H_bc, source_frequency, sigma, R,
         electric_susceptibility, magnetic_susceptibility, tol):
    input_parameters = {"H_bc": H_bc,
                        "source_frequency": source_frequency,
                        "sigma": sigma,
                        "electric_susceptibility": electric_susceptibility,
                        "magnetic_susceptibility": magnetic_susceptibility,
                        }
    tdata = [run_truchas(input_parameters, fdme=True),
             run_truchas(input_parameters, fdme=False),
             ]

    r_fd, Q_fd = get_heat_source(tdata[0])
    r_td, Q_td = get_heat_source(tdata[1])
    r_fd_em, Emag = get_Emag(tdata[0])

    # filter out the interior where Q is very small
    cutoff = 0.005
    Q_fd = np.extract(r_fd > cutoff, Q_fd)
    Q_td = np.extract(r_fd > cutoff, Q_td)
    r_fd = np.extract(r_fd > cutoff, r_fd)

    Qex, _ = analytic_soln(r_fd, H_bc, source_frequency, sigma, R,
                           electric_susceptibility, magnetic_susceptibility)
    _, Emagex = analytic_soln(r_fd_em, H_bc, source_frequency, sigma, R,
                              electric_susceptibility, magnetic_susceptibility)

    nfail = 0
    nfail += truchas.compare_max_rel(Q_fd, Qex, tol, "Q-FD", 0.0)
    nfail += truchas.compare_max_rel(Q_td, Qex, tol, "Q-TD", 0.0)
    nfail += truchas.compare_max_rel(Q_fd, Q_td, 1e-2, "Q-FD-TD", 0.0)
    nfail += truchas.compare_max_rel(Emag, Emagex, tol, "E-FD", 0.0)

    return nfail


def run_test(tenv):
    R = 2e-2 # currently hard-coded via cubit
    electric_susceptibility = 0
    magnetic_susceptibility = 0
    H_bc = 1
    frequencies = [500.0, 2.45e9]
    tols = [5e-2, 1.0]
    sigma = 1.0

    nfail = 0
    for f, tol in zip(frequencies, tols):
        nfail += test(tenv, H_bc, f, sigma, R,
                      electric_susceptibility, magnetic_susceptibility, tol)

    return nfail


# for a direct run, not ctest
def plot_comparison(tenv):
    H_bc = 1
    #source_frequency = 1e10
    #source_frequency = 47713451.59236942
    #source_frequency = 47713451.59236942e2
    source_frequency = 2.45e9
    #source_frequency = 500.0
    #source_frequency = 1e1
    #sigma = 1e3 # graphite is 1e3 - 1e5
    #sigma = 1e5
    #sigma = 2.654418727984993e2
    sigma = 2.654418727984993
    #sigma = 0.002654418727984993
    #sigma = 1e-2
    #sigma = 0

    R = 2e-2 # currently hard-coded via cubit
    electric_susceptibility = 0
    magnetic_susceptibility = 0

    input_parameters = {"H_bc": H_bc,
                        "source_frequency": source_frequency,
                        "sigma": sigma,
                        "electric_susceptibility": electric_susceptibility,
                        "magnetic_susceptibility": magnetic_susceptibility,
                        }
    tdata = [run_truchas(input_parameters, fdme=True),
             run_truchas(input_parameters, fdme=False),
             ]

    r, Q, Emag = heat_analytic(H_bc, source_frequency, sigma, R,
                               electric_susceptibility, magnetic_susceptibility)
    plot_results(r, Q, Emag, tdata,
                 f"uniform_source = {H_bc:.1e}, source_frequency = {source_frequency:.1e}, sigma = {sigma:.1e}")


if __name__ == "__main__":
    # tenv = truchas.TruchasEnvironment.default(".", overwrite_output=True)
    # plot_comparison(tenv)
    tenv = truchas.TruchasEnvironment.default()
    nfail = run_test(tenv)
    assert nfail == 0

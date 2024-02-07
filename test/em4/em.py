#!/usr/bin/env python3

import time
import os
import multiprocessing

import numpy as np
import scipy as sp
import scipy.special as spsp
import scipy.constants
import matplotlib as mpl
import matplotlib.pyplot as plt

import truchas


def run_truchas(input_parameters, emfd=False):
    input_parameters["emfd"] = "T" if emfd else "F"
    input_file = "emfd.inp" if emfd else "emtd.inp"
    tenv.generate_input_deck(input_parameters, "template-em.inp", input_file)
    t_elapsed = -time.time()
    stdout, tdata = tenv.truchas(1, input_file)
    t_elapsed += time.time()
    solver = "EMFD" if emfd else "EMTD"
    print(f"{solver} elapsed {t_elapsed:.2f} seconds.")
    return tdata


def get_heat_source(tdata):
    sid = tdata.num_series()
    xc = tdata.centroids()
    r_truchas = np.sqrt(xc[:,0]**2 + xc[:,1]**2)
    Q_truchas = tdata.field(sid, "Joule_P")
    return r_truchas, Q_truchas


def heat_analytic(H_bc, source_frequency, sigma, R,
                  electric_susceptibility, magnetic_susceptibility):
    mu = sp.constants.mu_0 * (1 + magnetic_susceptibility)
    epsilon = sp.constants.epsilon_0 * (1 + electric_susceptibility)
    omega = 2 * np.pi * source_frequency
    B_bc = mu * H_bc

    c = mu * (1j * omega * sigma - omega**2 * epsilon)
    sqrt_c = np.sqrt(c)
    A = B_bc / spsp.jv(0, sqrt_c * R)
    r = np.linspace(0, R, 100)
    B = A * spsp.jv(0, sqrt_c * r)
    E = -1j * omega * A / sqrt_c * spsp.jvp(0, sqrt_c * r)
    E2 = np.abs(E)**2
    Q = (sigma + omega * np.imag(epsilon)) * np.abs(E)**2 / 2

    print(f"mu: {mu}")
    print(f"epsilon: {epsilon}")
    print(f"sigma: {sigma}")
    print(c)
    print(sqrt_c)
    print(B_bc, B[-1])

    return r, Q, E2


def plot_results(r, Q, tdata, title):
    r_fd, Q_fd = get_heat_source(tdata[0])
    if len(tdata) > 1: r_td, Q_td = get_heat_source(tdata[1])

    nrows = 1; ncols = 1
    aspect_ratio = 4/3 * ncols / nrows
    height = 6 * nrows
    fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(aspect_ratio*height, height))
    ax.set_title(title)
    ax.set_xlabel("r")
    ax.set_ylabel("Q")
    ax.plot(r_fd, Q_fd, '.', label="EMFD")
    if len(tdata) > 1: ax.plot(r_td, Q_td, '.', label="EMTD")
    ax.plot(r, Q, '-', label="analytic")
    ax.legend(loc="upper left")
    fig.tight_layout();
    fig.savefig("em.png", dpi=200)


if __name__ == "__main__":
    tenv = truchas.TruchasEnvironment.default(".", overwrite_output=True)

    H_bc = 1
    #source_frequency = 1e10
    #source_frequency = 47713451.59236942
    source_frequency = 47713451.59236942e2
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
    tdata = [run_truchas(input_parameters, emfd=True),
             #truchas.TruchasData("emfd_output/emfd.h5"),
             run_truchas(input_parameters, emfd=False),
             ]

    r, Q, E2 = heat_analytic(H_bc, source_frequency, sigma, R,
                             electric_susceptibility, magnetic_susceptibility)
    plot_results(r, Q, tdata,
                 f"uniform_source = {H_bc:.1e}, source_frequency = {source_frequency:.1e}, sigma = {sigma:.1e}")

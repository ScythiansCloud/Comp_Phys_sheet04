#settings


import numpy as np

def init():
    
    global nsteps            # number of time step to analyze
    nsteps = 20000
    global eqsteps          # eqsteps
    eqsteps = 10000
    global mass              # mass of the LJ particles (gram/mole)
    mass = 39.95
    global kb                # boltzmann's constant (g*nm^2/fs^2/mole/K) 
    kb =0.0019858775 * 4.184e-6
    global Tdesired          # temperature of the experiment in K
    Tdesired = 300.
    global eps               # eps in LJ (g*nm^2/fs^2/mole)
    eps = 0.29788162 * 4.184e-6         # CHANGE TO: 0.1488 kcal/mol
    global sig                # r0 in LJ (nm)
    sig = 0.188

    # Wall values here
    global epsWall
    epsWall = 1.4887 * 4.184e-6     # CAN BE CHANGED LATER
    global sigWall
    sigWall = 0.0376e-9             # CAN BE CHANGED LATER

    global cutoff            # cutoff arbitrary at 2.5 r0
    cutoff = 2.5*sig
    global dt          # time step (fs)
    dt = 1
    global dr 
    dr = sig/30
    global Nanalyze
    Nanalyze = 10

    
    # number of particle = n1*n2 distributed on s square lattice
    # global n
    n12 = 8
    n3 = 12
    global N #particle number
    N = n12 * n12 * n3
    global rho
    rho = 0.5*sig**(-3) 
    global L # box length
    L = (N/rho)**(1/3)

    
    global deltaxyz # lattice parameter to setup the initial configuration on a lattice
    # deltaxyz = L/n
    
    #rescaling of temperature
    global Trescale
    Trescale = 1 #1 = rescale temperature; 0 = no rescaling
    
    

    

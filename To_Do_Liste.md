To Do Liste

Inititialization adaption:
- remove pbc of z component in force
- adapt all constants accordingly in settings


Write ForceWall():
- 





############### force ###############

def ForceWall(z, N, L, epsWall, sigWall):
    cutoffWall = 2.5 * sigWall
    sigW3 = 3*sigWall**3
    sigW9 = 9*sigWall**9

    fzw = np.zeros(N)

    prefac = 3*np.sqrt(3) *0.5 * epsWall

    for i in prange(N):
        # ziw: z distance of particle i with wall w (left zw1=0 & right zw2=2L)
        ziw = z % (2*L)

        if ziw < cutoffWall:
            fzw[i] = prefac * sigW3/ziw**4 - sigW9/ziw**10



############# settings ################

def init():
    
    global nsteps            # number of time step to analyze
    nsteps = 10000
    global eqsteps          # eqsteps
    eqsteps = 10000
    global mass              # mass of the LJ particles (gram/mole)
    mass = 39.95
    global kb                # boltzmann's constant (g*nm^2/fs^2/mole/K) 
    kb =0.0019858775 * 4.184e-6
    global Tdesired          # temperature of the experiment in K
    Tdesired = 300.
    global eps               # eps in LJ (g*nm^2/fs^2/mole)
    eps = 0.29788162 * 4.184e-6     # CHANGE TO: 0.1488 kcal/mol
    global sig                # r0 in LJ (nm)
    sig = 0.188
    global cutoff            # cutoff arbitrary at 2.5 r0
    cutoff = 2.5*sig
    global dt          # time step (fs)
    dt = 1
    global dr 
    dr = sig/15
    global Nanalyze
    Nanalyze = 5

    epsWall = 1.4887 * 4.184e-6  # CAN BE CHANGED LATER
    sigWall = 0.0376e-9      # CAN BE CHANGED LATER

    
    # number of particle = n1*n2 distributed on s square lattice
    global n
    n12 = 8
    n3 = 12
    global N #particle number
    N = n12 * n12 * n3
    global rho
    rho = 0.25*sig**(-3) 
    global L # box length
    L = (N/rho)**(1/3)

    
    global deltaxyz # lattice parameter to setup the initial configuration on a lattice
    # needs discussion herer?
    # deltaxyz = L/n
    
    #rescaling of temperature
    global Trescale
    Trescale = 1 #1 = rescale temperature; 0 = no rescaling
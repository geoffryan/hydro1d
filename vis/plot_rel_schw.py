import sys
import math
import numpy as np
import matplotlib.pyplot as plt
import readGrid as rg
    
M = 1.0
gam = 5.0/3.0

def metric4_schw_ks(r, M):
    a = 1.0/np.sqrt(1.0+2*M/r)
    g00 = -1.0 - 2*M/r
    g0r = 2*M/r
    g0p = np.zeros(r.shape)
    grr = 1.0 - 2*M/r
    grp = np.zeros(r.shape)
    gpp = 1.0/(r*r)

    return a, g00, g0r, g0p, grr, grp, gpp

def metric3_schw_ks(r, M):
    grr = r/(r+2*M)
    grp = np.zeros(r.shape)
    gpp = 1.0/(r*r)

    return grr, grp, gpp

def plot_dat(ax, x, dat, fmt, xlabel=None, ylabel=None, xlim=None, ylim=None, 
            **kwargs):

    ax.plot(x, dat, fmt, **kwargs)

    if xlabel is not None:
        ax.set_xlabel(xlabel)
    if ylabel is not None:
        ax.set_ylabel(ylabel)
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)

def plot_grid(filename, show=False):

    ver, t, xf, r, rho, P, ur, up = rg.readGrid(filename)
    
    a, g00, g0r, g0p, grr, grp, gpp = metric4_schw_ks(r, M)
    g3rr, g3rp, g3pp = metric3_schw_ks(r, M)

    w = np.sqrt(1.0 + g3rr*ur*ur + g3pp*up*up + 2*g3rp*ur*up)
    u0 = (w/a - g0r*ur - g0p*up) / g00
    uR = g0r*u0 + grr*ur + grp*up
    uP = g0p*u0 + grp*ur + gpp*up

    rhoh = rho + gam/(gam-1) * P
    cs = np.sqrt(gam*P/rhoh)
    cs4 = cs/np.sqrt(1-cs*cs)
    u2 = g3rr*ur*ur + 2*g3rp*ur*up + g3pp*up*up

    mach = np.sqrt(u2) / cs4

    rin = r[r<6.0*M]
    rou = r[r>=6.0*M]
    UR = np.zeros(r.shape)
    UP = np.zeros(r.shape)
    Ur = np.zeros(r.shape)
    Up = np.zeros(r.shape)
    U0 = np.zeros(r.shape)

    UR[r<6.0*M] = -np.power(6.0*M/rin - 1.0, 1.5) / 3.0
    UP[r>=6.0*M] = np.sqrt(M/(rou*rou*(rou-3.0*M)))
    UP[r<6.0*M] = 2.0*math.sqrt(3.0) * M/(rin*rin)

    Ur[r>=6.0*M] = 2.0*M/(rou*np.sqrt(1.0-3*M/rou))
    Ur[r<6.0*M] = (4.0*math.sqrt(2)*M/rin - np.power(6*M/rin-1,1.5)) / (3.0*(1-2*M/rin))
    Up[r>=6.0*M] = np.sqrt(M*rou / (1.0-3*M/rou))
    Up[r<6.0*M] = 2.0*math.sqrt(3.0) * M
    U0[r>=6.0] = -(1.0-2*M/rou) / np.sqrt(1.0-3*M/rou)
    U0[r<6.0] = -math.sqrt(8.0)/3.0

    URR = g0r*U0 + grr*Ur + grp*Up


    fig = plt.figure(figsize=(12,10))

    xlim = (xf[0], xf[-1])

    ax = fig.add_subplot(3,3,1)
    plot_dat(ax, r, rho, 'k+', xlabel=r"$r$", ylabel=r"$\rho$", xlim=xlim)
    ax.set_yscale("log")
    ax = fig.add_subplot(3,3,2)
    plot_dat(ax, r, P, 'k+', xlabel=r"$r$", ylabel=r"$P$", xlim=xlim)
    ax.set_yscale("log")
    ax = fig.add_subplot(3,3,3)
    plot_dat(ax, r, mach, 'k+', xlabel=r"$r$", ylabel=r"$\mathcal{M}$",
                xlim=xlim)
    ax.set_yscale("log")

    ax = fig.add_subplot(3,3,4)
    plot_dat(ax, r, uR, 'k+', xlabel=r"$r$", ylabel=r"$u^r$", xlim=xlim)
    ax.plot(r, UR, 'r')
    ax.plot(r, URR, 'b')
    ax = fig.add_subplot(3,3,5)
    plot_dat(ax, r, uP, 'k+', xlabel=r"$r$", ylabel=r"$u^\phi$", xlim=xlim)
    ax.plot(r, UP, 'r')
    ax = fig.add_subplot(3,3,6)
    plot_dat(ax, r, w, 'k+', xlabel=r"$r$", ylabel=r"$\mathcal{W}$",
                xlim=xlim)

    ax = fig.add_subplot(3,3,7)
    plot_dat(ax, r, ur, 'k+', xlabel=r"$r$", ylabel=r"$u_r$", xlim=xlim)
    ax.plot(r, Ur, 'r')
    ax = fig.add_subplot(3,3,8)
    plot_dat(ax, r, up, 'k+', xlabel=r"$r$", ylabel=r"$u_\phi$", xlim=xlim)
    ax.plot(r, Up, 'r')

    ax = fig.add_subplot(3,3,9)
    plot_dat(ax, rho, P, 'k+', xlabel=r"$\rho$", ylabel=r"$P$") 
    ax.set_xscale("log")
    ax.set_yscale("log")
    
    plt.suptitle(r"$t = {0:.3g}$".format(t))

    plt.tight_layout()

    name = ".".join(filename.split(".")[:-1])
    fig.savefig(name+".png")
    if show:
        plt.show()
    plt.close()

if __name__== "__main__":

    if len(sys.argv) < 2:
        print("I need some files to nom.\n")
        sys.exit()

    if len(sys.argv) == 2:
        print("Plotting {0:s}".format(sys.argv[1]))
        plot_grid(sys.argv[1], show=True)

    else:
        for filename in sys.argv[1:]:
            print("Plotting {0:s}".format(filename))
            plot_grid(filename, show=False)

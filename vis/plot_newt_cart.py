import sys
import numpy as np
import matplotlib.pyplot as plt
import readGrid as rg

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

def plot_grid(filename):

    ver, t, xf, x, rho, P, v = rg.readGrid(filename)

    fig = plt.figure()

    ax = fig.add_subplot(2,2,1)
    plot_dat(ax, x, rho, 'k+', xlabel=r"$x$", ylabel=r"$\rho$")
    ax = fig.add_subplot(2,2,2)
    plot_dat(ax, x, v, 'k+', xlabel=r"$x$", ylabel=r"$v_x$")
    ax = fig.add_subplot(2,2,3)
    plot_dat(ax, x, P, 'k+', xlabel=r"$x$", ylabel=r"$P$")

    plt.suptitle(r"$t = {0}$".format(t))

    name = ".".join(filename.split(".")[:-1])
    fig.savefig(name+".png")
    plt.close()

if __name__== "__main__":

    if len(sys.argv) < 2:
        print("I need some files to nom.\n")
        sys.exit()

    for filename in sys.argv[1:]:
        print("Plotting {0:s}".format(filename))
        plot_grid(filename)

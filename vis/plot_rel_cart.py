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

def plot_grid(filename, show=False):

    ver, t, xf, x, rho, P, u1, u2 = rg.readGrid(filename)

    w = np.sqrt(1 + u1*u1 + u2*u2)
    v1 = u1/w
    v2 = u2/w

    fig = plt.figure()

    xlim = (xf[0], xf[-1])

    ax = fig.add_subplot(2,2,1)
    plot_dat(ax, x, rho, 'k+', xlabel=r"$x$", ylabel=r"$\rho$", xlim=xlim)
    ax.set_yscale("log")
    ax = fig.add_subplot(2,2,2)
    plot_dat(ax, x, P, 'k+', xlabel=r"$x$", ylabel=r"$P$", xlim=xlim)
    ax.set_yscale("log")
    ax = fig.add_subplot(2,2,3)
    plot_dat(ax, x, v1, 'k+', xlabel=r"$x$", ylabel=r"$v_x$", xlim=xlim)
    ax = fig.add_subplot(2,2,4)
    plot_dat(ax, x, v2, 'k+', xlabel=r"$x$", ylabel=r"$v_y$", xlim=xlim)

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

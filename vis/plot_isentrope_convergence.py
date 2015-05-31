import sys
import math
import numpy as np
import matplotlib.pyplot as plt
import readGrid as rg

rhoRef = 1.0
PRef = 1.0
xRef = 1.5
L = 0.4
a = 1.0
GAM = 5.0/3.0

K = PRef / math.pow(rhoRef, GAM)

def profile(x):

    rho = rhoRef * np.ones(x.shape)

    dx = np.fabs((x-xRef)/L)
    csRef = math.sqrt(GAM * PRef / rhoRef)
    
    J = -2.0 * csRef / (GAM-1.0)

    inds = (dx <= 1.0)

    rho[inds] += a * rhoRef * (dx[inds]*dx[inds]-1.0)**4

    P = K * np.power(rho, GAM)
    cs = np.sqrt(GAM * P / rho)

    v = J + 2.0 * cs / (GAM-1.0)

    return rho, P, v

def isentropeDirect(t, x0):

    rho, P, v = profile(x0)
    la = v + np.sqrt(GAM*P/rho)
    x = x0 + la*t

    return x, rho, P, v

def isentrope(t, x, err=1.0e-10):
    
    csMax = math.sqrt(GAM * K * math.pow(rhoRef*(1.0+a), GAM-1.0))
    xb = x.copy()
    xa = x - 10*csMax * t

    x0 = 0.5*(xa+xb)
    dx = xb-xa

    i = 0

    while np.abs(dx/x0).mean() > err:
        rho0, P0, v0 = profile(x0)
        
        la = v0 + np.sqrt(GAM*P0/rho0)

        errs = la*t - x + x0
        xa[errs<0.0] = x0[errs<0.0]
        xb[errs>0.0] = x0[errs>0.0]
        dx = xb-xa
        x0 = 0.5*(xa+xb)

        i += 1

    rho, P, v = profile(x0)

    return x0, rho, P, v

def calc_err(filename):

    ver, t, xf, x, rho, P, v1, v2 = rg.readGrid(filename)
    dx = xf[1:] - xf[:-1]

    x0, rho_exact, P_exact, v_exact = isentrope(t, x)

    err = (np.abs(rho[2:-2]-rho_exact[2:-2])*dx[2:-2]).sum() \
            / (np.abs(rho_exact[2:-2])*dx[2:-2]).sum()

    N = x.shape[0] - 4

    return N, err, t, x, (rho,P,v1), (rho_exact,P_exact,v_exact)

def plot_convergence(filenames):

    fig = plt.figure()

    ax1 = fig.add_subplot(2,2,1)
    ax2 = fig.add_subplot(2,2,2)
    ax3 = fig.add_subplot(2,2,3)
    ax4 = fig.add_subplot(2,2,4)

    Ns = []
    errs = []
    ts = []

    cols = [(31.0/255, 119.0/255, 180.0/255), (255.0/255, 127.0/255, 14.0/255),
        (44.0/255, 160.0/255, 44.0/255), (214.0/255, 39.0/255, 40.0/255),
        (148.0/255, 103.0/255, 189.0/255), (140.0/255, 86.0/255, 75.0/255),
        (227.0/255, 119.0/255, 194.0/255), (127.0/255, 127.0/255, 127.0/255),
        (188.0/255, 189.0/255, 34.0/255), (23.0/255, 190.0/255, 207.0/255)]

    lcols = [(174.0/255,199.0/255,232.0/255), (255.0/255,187.0/255,120.0/255),
        (152.0/255, 223.0/255, 138.0/255), (255.0/255, 152.0/255, 150.0/255),
        (197.0/255, 176.0/255, 213.0/255), (196.0/255, 156.0/255, 148.0/255),
        (247.0/255, 182.0/255, 210.0/255), (199.0/255, 199.0/255, 199.0/255),
        (219.0/255, 219.0/255, 141.0/255), (158.0/255, 218.0/255, 229.0/255)]


    for i,f in enumerate(filenames):
        N, err, t, x, dat, exact = calc_err(f)
        x0 = np.linspace(x.min(), x.max(), 1000)
        xe, rhoe, Pe, ve = isentropeDirect(t, x0)
        Ns.append(N)
        errs.append(err)
        ts.append(t)

        ax1.plot(x, dat[0], marker='x', linestyle=' ', color=cols[i%10])
        #ax1.plot(x, exact[0], marker='o', linestyle=' ', color=cols[i])
        ax1.plot(xe, rhoe, marker='', linestyle='-', color=lcols[i%10])
        ax2.plot(x, dat[1], marker='x', linestyle=' ', color=cols[i%10])
        #ax2.plot(x, exact[1], marker='o', linestyle=' ', color=cols[i])
        ax2.plot(xe, Pe, marker='', linestyle='-', color=lcols[i%10])
        ax3.plot(x, dat[2], marker='x', linestyle=' ', color=cols[i%10])
        #ax3.plot(x, exact[2], marker='o', linestyle=' ', color=cols[i])
        ax3.plot(xe, ve, marker='', linestyle='-', color=lcols[i%10])

    p = np.polyfit(np.log(Ns), np.log(errs), 1)
    print("Convergence Order: {0:.2f}".format(-p[0]))

    NNs = np.logspace(np.log10(Ns).min(),np.log10(Ns).max(),base=10.0,num=100)

    ax4.plot(NNs, math.exp(p[1]) * np.power(NNs, p[0]), ls='-', lw=2, 
                color=lcols[0])
    ax4.plot(Ns, errs, 'kx')

    ax1.set_xlabel(r'$x$')
    ax1.set_ylabel(r'$\rho$')
    ax2.set_xlabel(r'$x$')
    ax2.set_ylabel(r'$P$')
    ax3.set_xlabel(r'$x$')
    ax3.set_ylabel(r'$v$')

    ax4.set_xlabel(r'$N$')
    ax4.set_ylabel(r'$\epsilon$')
    ax4.set_xscale('log')
    ax4.set_yscale('log')

    plt.tight_layout()

    fig.savefig("isentrope_err.png")

if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("I need some files to nom.\n")
        sys.exit()

    plot_convergence(sys.argv[1:])

    plt.show()


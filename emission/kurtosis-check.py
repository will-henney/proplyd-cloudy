import numpy as np

def gauss(x, x0, sig):
    """
    Gaussian profile
    """
    return np.exp( -(x - x0)**2 / (2*sig**2) )


def twogauss(x, x0, sig):
    """
    Double Gaussian, separated by 2 x0
    """
    return gauss(x, x0, sig) + gauss(x, -x0, sig)


VMAX = -10.0
N = 1000

V = np.linspace(-VMAX, VMAX, N)

sigma1 = 1.0                    # width of individual Gaussians

Sep = np.linspace(0.0, VMAX/10, 100) # (half) separation between Gaussians

Variance = list()
for sep in Sep:
    Y = twogauss(V, sep, sigma1)
    sum_ = Y.sum()
    mean = np.sum(V*Y)/sum_
    variance = np.sum(Y*(V-mean)**2)/sum_
    Variance.append(variance)
Variance = np.array(Variance)
    
import matplotlib.pyplot as plt

plt.plot(Sep**2, Variance)

plt.savefig("kurtosis-check.png")


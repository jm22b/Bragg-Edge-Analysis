import matplotlib.pyplot as plt
import numpy as np
import scipy.special
import scipy.signal
import warnings
from scipy.optimize import curve_fit, OptimizeWarning

x = []
y = []
f = open('edgedata.txt')
for line in f:
    s = line.split()
    x.append(float(s[0]))
    y.append(float(s[1]))
f.close()

x = np.array(x)
y = np.array(y)

"""
Parameters:

p1 = a_left - ordinate of line for data < lambda_0
p2 = m_left - slope of line for data < lambda_0
p3 = a_right - ordinate of line for data > lambda_0
p4 = m_right - slope of line for data > lambda_0
p5 = lambda_0 - Bragg Edge Position
p6 = sigma - symmetric (gaussian) broadening of edge
p7 = tau - asymmetric (decaying exponential) broadening of edge
"""

# filtering noise to extract parameters
b, a = scipy.signal.butter(3, 0.05)
zi = scipy.signal.lfilter_zi(b, a)
z, _ = scipy.signal.lfilter(b, a, y, zi=zi*y[0])
z2, _ = scipy.signal.lfilter(b, a, y, zi=zi*z[0])
ynew = scipy.signal.filtfilt(b, a, y)

dy = np.diff(ynew)
dx = np.diff(x)
dydx = dy/dx
edgeIndex = np.argmax(dydx)

lambda_0 = x[edgeIndex]
print "lambda_0: %f" % lambda_0

m_left = (ynew[0:edgeIndex-20][-1] - ynew[0:edgeIndex-20][0])/(x[0:edgeIndex-20][-1] - x[0:edgeIndex-20][0])
a_left = np.median(y[0:edgeIndex])
print "a_left: %f" % a_left
print "m_left: %f" % m_left


m_right = (ynew[edgeIndex+20:][-1] - ynew[edgeIndex+20:][0])/(x[edgeIndex+20:][-1] - x[edgeIndex+20:][0])
a_right = np.median(y[edgeIndex:])
print "a_right: %f" % a_right
print "m_right: %f" % m_right


sigma = 0.004
tau = 0.04

# far right fitting
def rightFunc(x, m, c):
    return np.exp(-(m*x + c))

popt_r, pcov = curve_fit(rightFunc, x[edgeIndex+30:], y[edgeIndex+30:], p0=[m_right, a_right])
m_right = popt_r[0]
a_right = popt_r[1]

# far left fitting
def leftFunc(x, b, a):
    return np.exp(-(m_right*x + a_right))*np.exp(-(a + b*x))

popt_l, pcov = curve_fit(leftFunc, x[:edgeIndex-20], y[:edgeIndex-20], p0=[m_left, a_left])
m_left = popt_l[0]
a_left = popt_l[1]

def centralFunc(x, lambda_0, sigma, tau):

    x_sig = -(x - lambda_0) / (np.sqrt(2)*sigma)
    x_tau = -(x - lambda_0) / tau
    sig_tau = sigma / tau
    b = 0.5*(scipy.special.erfc(x_sig) - np.exp(x_tau + (0.5*sig_tau**2))*scipy.special.erfc(x_sig + sig_tau))
    return np.exp(-(a_right + m_right*x))*(np.exp(-(a_left + m_left*x)) + (1-np.exp(-(a_left + m_left*x)))*b)

popt_c, pcov = curve_fit(centralFunc, x, y, p0=[lambda_0, sigma, tau])

plt.plot(x, y, 'x')
plt.plot(x, centralFunc(x, popt_c[0], popt_c[1], popt_c[2]))
plt.show()
plt.close()

"""
warnings.simplefilter("error", OptimizeWarning)
try:
    
    popt, pcov = curve_fit(analyticalFunc, x, y, p0=[a_left, m_left, a_right, m_right, lambda_0, sigma, tau])
except (RuntimeError, OptimizeWarning):
    
    print "exception"
    plt.plot(x, y, 'x')
    plt.plot(x, ynew)
    plt.plot(x, analyticalFunc(x, a_left, m_left, a_right, m_right, lambda_0, sigma, tau))
    plt.show()
    plt.close()
"""

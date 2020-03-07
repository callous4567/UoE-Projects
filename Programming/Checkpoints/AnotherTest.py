import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# function that returns dy/dt
def model(y,t):
    k = 0.3
    dydt = -k * y
    return dydt

# initial condition
y0 = 5

# time points
t = np.linspace(0,20)

# solve ODE
y = odeint(model,y0,t)

# plot results
plt.plot(t,y)
plt.xlabel('time')
plt.ylabel('y(t)')
plt.show()# -*- coding: utf-8 -*-
from scipy import integrate
from pylab import * # for plotting commands

def rlc(A,t): 
    Vc,x=A 
    V = 1.0 #voltageSource
    R = 5.0
    L=100.0e-9 #100nH
    C = 1.0e-9 #1nF
    res=array([x,(V-Vc-(x*R*C))/(L*C)])
    return res

time = linspace(0.0,0.6e-6,1001)
vc,x = integrate.odeint(rlc,[0.0,0.0],time).T
i=1.0e-9*x
figure()
plot(time,vc)
xlabel('t')
ylabel('Vc')
show()

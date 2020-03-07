import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as pltick
import sys
import time
import random

class utility:
    def time_delay(text):
        for c in text:
            sys.stdout.write(c)
            sys.stdout.flush()
            resttime = random.uniform(0.001, 0.005)
            time.sleep(resttime)
        print()

    def value_getter(value_name):
        utility.time_delay(("Please provide the value for the quantity: {}").format(value_name))
        try:
            val_val = float(input())
            return val_val
        except:
            utility.time_delay("Hey. Didn't work.")

    def txt_reader():
        utility.time_delay("Please provide the file name. Extension is assumed .txt!)")
        filename = str(input())
        utility.time_delay("How many lines to clip?")
        clip = int(input())
        filefile = open(("{0}.{1}").format(filename, "txt"))
        filefile_lines = filefile.readlines()
        filefile_lines = filefile_lines[clip:len(filefile_lines)]
        v_list = []
        c_list = []
        for line in filefile_lines:
            newline = line.split(",")
            try:
                v_list.append(float(newline[0]))
                c_list.append(float(newline[1]))
            except:
                utility.time_delay("Went wrong here. ##1")
        filefile.close()
        return v_list, c_list
            
class v_to_r:
    def __init__(self, volume, radius, area):
        self.volume = volume
        self.radius = radius
        self.area = area
    def input(self):
        vol = utility.value_getter("Volume")
        self.volume = vol
        utility.time_delay("What's the unit exponent of the volume? mm = -3")
        exponent = utility.value_getter("Exponent")
        self.volume = vol*((10**exponent)**3)
    def converter(self):
        vol = self.volume
        # Spherical 4/3*pi*r^3 = v *** 3V/4pi to the 1/3 is r.
        self.rad = (3*vol/(4*np.pi))**(1/3)
        self.area = 4*np.pi*(self.rad**2)
    def printer(self):
        utility.time_delay(("Your radius is {0:.2e} metres, surface area {1:.2e} metres squared").format(self.rad, self.area))

class quadratic_solve:
    def __init__(self, a, b, c):
        self.a = a
        self.b = b
        self.c = c
    def value_retriever(self):
        utility.time_delay("ax^2 + bx + (c); provide the values.")
        self.a = utility.value_getter("a")
        self.b = utility.value_getter("b")
        self.c = utility.value_getter("c")
    def solver(self):
        discriminant = ((self.b)**2 - 4*(self.a)*(self.c))
        if discriminant > 0:
            value_1 = ((-self.b + np.sqrt(discriminant))/(2*self.a))
            value_2 = ((-self.b - np.sqrt(discriminant))/(2*self.a))
            utility.time_delay(("Two solutions. {0:.2e} and {1:.2e}").format(value_1, value_2))
        elif discriminant == 0:
            value_1 = -self.b/(2*self.a)
            utility.time_delay(("Simultaneous solution; {0:.2e}.").format(value_1))
        else:
            utility.time_delay("No solutions.")


class oscillator:
    def __init__(self, x, omega, gamma, a, b, p, t, k, m, n):
        self.x = []
        self.omega = []
        self.gamma = 0
        self.a = 0
        self.b = 0
        self.p = 0
        self.t = []
        self.k = 0
        self.m = 0
        self.n = 0
    def constraints(self):
        hej = utility.value_getter("Omega Initial")
        self.omega.append(hej)
        self.gamma = utility.value_getter("Gamma")
        self.n = utility.value_getter("Points")
        self.a = 1
        self.t = np.linspace(0, 5*np.pi*(1/self.omega[0]), self.n)
        if self.gamma > 2*self.omega[0]:
            self.p = (np.sqrt((self.gamma/2)**2 - self.omega[0]**2))
            self.a = 1
            self.b = self.gamma/(2*self.p)
        elif self.gamma == 2*self.omega[0]:
            self.a = 1
            self.b = self.gamma/2
        else:
            hij = np.sqrt(self.omega[0]**2 - (self.gamma/2)**2)
            self.omega.append(hij)
            self.a = 1
            self.b = self.gamma/(2*self.omega[1])
    def values(self):
        if self.gamma > 2*self.omega[0]:
            x_t = lambda t: (np.exp(-self.gamma*t*1/2))*(self.a*np.cosh(self.p*t) + self.b*np.sinh(self.p*t))
            self.x = [x_t(d) for d in self.t]
        elif self.gamma == 2*self.omega[0]:
            x_t = lambda t: (np.exp(-self.gamma*t*1/2))*(self.a + self.b*t)
            self.x = [x_t(d) for d in self.t]
        else:
            x_t = lambda t: (np.exp(-self.gamma*t*1/2))*(self.a*np.cos(self.omega[1]*t) + self.b*np.sin(self.omega[1]*t))
            self.x = [x_t(d) for d in self.t]
    def grapher(self):
        fig = plt.figure(111)
        ax1 = plt.subplot(111)
        ax1.plot(self.t, self.x, label = "x vs. t for oscillator", color = "black")
        plt.legend()
        plt.show()

class v_i:
    def __init__(self, filename, V, I, t, p):
        self.filename = filename
        self.V = V
        self.I = I
        self.t = []
        self.p = []
    def constraints(self):
        files = utility.txt_reader()
        self.V,self.I = files[0],files[1]
        for d in range(0, len(self.V)):
            p_t = np.log(self.V[d]*self.I[d])
            self.p.append(p_t)
            self.t.append(d*(1/25000))
            
        print(self.p)
        print(self.t)
    def grapher(self):
        fig = plt.figure(111)
        ax1 = plt.subplot(111)
        ax1.plot(self.t, self.p, color="black")
        plt.show()

class particle_iterator:
    def __init__(self, g, x, y, t, vel, velx, vely, theta, beta, step):
        self.vel = []
        self.velx = []
        self.vely = []
        self.theta = 0
        self.beta = 0
        self.step = 0
        self.x = [0]
        self.y = [0]
        self.t = [0]
        self.g = 9.806
    def constraints(self):
        self.beta = utility.value_getter("Beta constant")
        self.step = utility.value_getter("Time Step in seconds")
        vel0 = utility.value_getter("Initial Velocity")
        self.vel.append(vel0)
        theta0 = utility.value_getter("Initial Phase Angle")
        self.theta = (theta0)
        self.velx.append(vel0*np.cos(theta0))
        self.vely.append(vel0*np.sin(theta0))
    def iterator(self):
        ax = lambda B,v,vx: -B*v*vx
        ay = lambda B,v,vy,g: -B*v*vy - g
        vx = lambda vx,ax,dt: vx + ax*dt
        vy = lambda vy,ay,dt: vy + ay*dt
        velmag = lambda vx,vy: (vx**2 + vy**2)**0.5
        while True:
            self.x.append(self.x[-1] + self.velx[-1]*self.step)
            self.y.append(self.y[-1] + self.vely[-1]*self.step)
            self.t.append(self.step)
            self.velx.append(vx(self.velx[-1], ax(self.beta,self.vel[-1],self.velx[-1]), self.step))
            self.vely.append(vy(self.vely[-1], ay(self.beta,self.vel[-1],self.vely[-1], self.g), self.step))
            self.vel.append(velmag(self.velx[-1], self.vely[-1]))
            if self.y[-1] < 0:
                break
    def grapher(self):
        fig = plt.figure(111)
        ax1 = plt.subplot(111)
        ax1.plot(self.x, self.y, color="black")
        plt.show()

hej = particle_iterator(0,0,0,0,0,0,0,0,0,0)
hej.constraints()
hej.iterator()
hej.grapher()

import numpy as np
from numpy import linalg as LA
from Particle3D import Particle3D

# Gravforce Function, calculates the force and direction of the gravitational potential
def gravforce(p1, p2, G):
    sep = Particle3D.vect_sep(p1,p2)
    Norm = LA.norm(sep)
    force = G*p1.mass*p2.mass/(Norm**2)
    direct = sep/Norm
    return force*direct

#gravpotential function, calculates the potential energy of gravity
def gravpotential(p1,p2,G):
    sep = Particle3D.vect_sep(p1,p2)
    Norm = LA.norm(sep)
    potential = -G*p1.mass*p2.mass/Norm
    return potential

def main():

    # Open input and output files
    parameters = open("parameters.dat", "r")
    particle_file = open("particles.dat", "r")
    vmd_file = open("vmdoutput.xyz", "w")
    Energy = open("Energy.dat", "w")
    
    #create a list of particles from the file
    particle_list = Particle3D.new_particle(particle_file)
    for line in parameters.readlines():
            partline = line.split(", ")
            
    # Set up simulation parameters from "parameters.dat" file
    dt = float(partline[0])
    numstep = int(partline[1])
    time = 0.0
    G = float(partline[2])
    

    
    Bodylist = []
    # creates some Particle3D objects from the list of particles
    for i in range (len(particle_list)):
        Bodylist.append(Particle3D(particle_list[i]))
        print(Bodylist[i].label)
     
    #calculates the COM velocity
    P = 0.0
    M = 0.0
    for i in range (len(particle_list)):
        P += Bodylist[i].velocity*Bodylist[i].mass
        M += Bodylist[i].mass
    Vcom = P/M
        
    #accounts for the COM correction for the system
    for i in range (len(particle_list)):
        Bodylist[i].velocity = Bodylist[i].velocity -Vcom
    
    #initialises lists and data storage
    forcelist = []
    newforcelist = []
    filehandlelist = []
    potential = 0
    kinetic = 0
    
    #Calculate initial forces on each body before the start of the timeloop
    #For each body i, the gravational force between body i and every other body j is calculated
    #using the gravforce function, the forces are summed to resolve the resultant force
    #on body i, then the np array representing the force is appended as the ith element of the list
    for i in range(len(Bodylist)):
        force_initial = np.array([0.0,0.0,0.0])
        filehandlelist.append(Bodylist[i].label)
        for j in range(len(Bodylist)):
            if i!=j:
                force_initial += gravforce(Bodylist[i],Bodylist[j], G)
        forcelist.append(force_initial)
        
    #opens a list of outfiles to write in for the orbit calculation
    outfileparts = []
    for i in range (len(particle_list)):
        outfileparts.append(open(filehandlelist[i] + ".dat", "w"))
        
     #begins simulation   
    for t in range(numstep):
        #writes the VMD and planetary data files in the proper format
        vmd_file.write(str(len(Bodylist))+ "\n")
        vmd_file.write("Point = " + str(t) + "\n")
        for i in range(len(Bodylist)):
            write = Bodylist[i].stringy()
            outfileparts[i].write(str(write)+ "\n")
            vmd_file.write(Bodylist[i].label + " " + str(Bodylist[i].position[0]) + " " +
                                                   str(Bodylist[i].position[1]) + " "
                                                        + str(Bodylist[i].position[2]) + "\n")
            
        #uses the force and velocity calculated from the previous timestep to move each body               
        for i in range(len(Bodylist)):
            Bodylist[i].leap_pos2nd(dt, forcelist[i])
        
        #Caculated the resultant force on each body i, in a similar fashion to the initial force
        #calcuated before the start of the timeloop
        #For each body i a blank np array is created with 3 dimensions representing the 
        #the xyz components of the force
        #Once the resultant force on a body is calculated, the array is appended to a new force list
        #For use in the velocity verlet velocity calculation
        #The np array is then reset for the next body
        for i in range(len(Bodylist)):
            force_new_i = np.array([0.0,0.0,0.0])
            for j in range(len(Bodylist)):
                if i!=j:
                    force_new_i += gravforce(Bodylist[i],Bodylist[j],G)
            newforcelist.append(force_new_i)
        
        #Calculates the potential energy
        potential = 0   
        for i in range(len(Bodylist)):
            for j in range(len(Bodylist)):
                if j>i:
                    potential += gravpotential(Bodylist[i],Bodylist[j],G)
        
        #Calculates the kinetic energy
        kinetic = 0
        for i in range(len(Bodylist)):
            kinetic += Bodylist[i].kinetic_energy()

        #writes the total energy to an Energy File
        Energy.write(str(kinetic + potential) + "\n")
            
        #Velocity calculation, uses the previous force and current force calculated
        for i in range(len(Bodylist)):
            Bodylist[i].leap_velocity(dt,(forcelist[i]+newforcelist[i])/2)
            
            
        #Replaces the force for each body from the previous timestep to their current forces,
        #for use in the velocity calculation in the next timestep
        for i in range(len(Bodylist)):  
            forcelist[i] = newforcelist[i]
        #clears the previous force list
        newforcelist = []
        
        
        #Adds dt to the total time
        time += dt
    
            
    
if __name__ == "__main__":
    main()

print("Done")
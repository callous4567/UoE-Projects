"""
 CMod Ex2: Particle1D, a class to describe 1D particles
"""

class Particle1D(object):
    """
    Class to describe 1D particles.

    Properties:
    position(float) - position along the x axis
    velocity(float) - velocity along the x axis
    mass(float) - particle mass

    Methods:
    * formatted output
    * kinetic energy
    * first-order velocity update
    * first- and second order position updates
    """

    def __init__(self, pos, vel, mass):
        """
        Initialise a Particle1D instance
        
        :param pos: position as float
        :param vel: velocity as float
        :param mass: mass as float
        """
        self.position = pos
        self.velocity = vel
        self.mass = mass
    

    def __str__(self):
        """
        Define output format.
        For particle p=(2.0, 0.5, 1.0) this will print as
        "x = 2.0, v = 0.5, m = 1.0"
        """
        return "x = " + str(self.position) + ", v = " + str(self.velocity) + ", m = " + str(self.mass)

    
    def kinetic_energy(self):
        """
        Return kinetic energy as
        1/2*mass*vel^2
        """
        return 0.5*self.mass*self.velocity**2
        

    # Time integration methods
    def leap_velocity(self, dt, force):
        """
        First-order velocity update,
        v(t+dt) = v(t) + dt*F(t)

        :param dt: timestep as float
        :param force: force on particle as float
        """
        self.velocity += dt*force/self.mass


    def leap_pos1st(self, dt):
        """
        First-order position update,
        x(t+dt) = x(t) + dt*v(t)

        :param dt: timestep as float
        """
        self.position += dt*self.velocity


    def leap_pos2nd(self, dt, force):
        """
        Second-order position update,
        x(t+dt) = x(t) + dt*v(t) + 1/2*dt^2*F(t)

        :param dt: timestep as float
        :param force: current force as float
        """
        self.position += dt*self.velocity + 0.5*dt**2*force/self.mass


hey = Particle1D()
print(hey.position)

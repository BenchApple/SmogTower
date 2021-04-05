# Ben Chappell

# This file is meant to model newton's 2nd law for the electrostatic field that pulls point charges (smog particles)
# out of the air and into a filtration system.

import math
from forces import *
import matplotlib.pyplot as plt

# Sigma is charge density of plate
# q is charge of individual particle
# default plate voltage comes from https://en.wikipedia.org/wiki/Electrostatic_precipitator
# plate charge density as https://physics.stackexchange.com/questions/304828/surface-charge-density-of-parallel-plate-capacitor
# concentration default is given in micrograms per m^3 as given by the WHO
# TODO make sure i didn't mess up the reynold's number equation, bc that would be bad
class SmogCollector:
    def __init__(self, _rho_part=500, _d_part=0.0000025, _rho_fluid=1.225, _mu_fluid=0.000001803, \
                 _plate_voltage=20000, _plate_dist=0.05, _H=400, _x_init=100, _v_init=-0.05, _a_init=-0.01, \
                 _concentration=150, _t_step=0.002):
        e_0 = 8.854187 * math.pow(10, -12)

        # I have literally no clue how to find the charge of a particle so im gonna use this
        self.q = 1.906 * math.pow(10, -4)

        self.rho_part = _rho_part
        self.d_part = _d_part
        self.rho_fluid = _rho_fluid
        self.mu_fluid = _mu_fluid
        self.plate_voltage = _plate_voltage
        self.plate_dist = _plate_dist
        self.H = _H
        self.x_init = _x_init
        self.v_init = _v_init
        self.a_init = _a_init
        self.c_init = _concentration
        self.t_step = _t_step

        self.sigma = (e_0 * self.plate_voltage) / self.plate_dist

        self.v_apparent = -self.v_init

        # The following keep track of the current values of the force, accel, velo, and dist
        self.cur_a = self.a_init
        self.cur_v = self.v_init
        self.cur_x = self.x_init
        self.cur_c = self.c_init

        self.part_mass = mass(self.rho_part, self.d_part)
        self.part_mass = _rho_part * _concentration * math.pow(10, -6) * _x_init

        self.cur_force = self.part_mass * self.cur_a

        # Each of the following lists track the accel, velocity, position, and concentration over time
        self.accel = [self.cur_a]
        self.velo = [self.cur_v]
        self.pos = [self.cur_x]
        self.concen = [self.cur_c]
        self.force = [self.cur_force]

        # t stores what the current time is.
        self.t = 0

    # Gets the value of F where F is defined as the bouyancy force minus the gravitational 
    # minus the force from the electric field
    def get_F(self):
        return (bouyant_force(self.rho_part, self.d_part) - grav_force(self.rho_part, self.d_part) - electric_collector(self.q, self.sigma))
    
    # Prints all of the forces out.
    def print_forces(self):
        # This basically just tests the forces file.
        print("Grav")
        print(-grav_force(self.rho_part, self.d_part))
        print("Bouyant")
        print(bouyant_force(self.rho_part, self.d_part))
        print("Reynolds")
        print(reynolds(self.rho_fluid, self.d_part, self.v_apparent, self.mu_fluid))
        print("Drag Coeff")
        print(drag_coeff(self.rho_fluid, self.d_part, self.v_apparent, self.mu_fluid))
        print("Drag Force")
        print(drag_force(self.rho_fluid, self.d_part, self.v_apparent, self.mu_fluid, self.cur_v))
        print("Electric Collector")
        print(-electric_collector(self.q, self.sigma))
        print("Electric Other")
        print(electric_other(self.q, self.cur_c, self.cur_x, self.H))
        print("Mass")
        print(self.part_mass)


    # Iterates time by one step, updating all of the elements in the process.
    def step_time(self):
        new_x = self.cur_x + (self.cur_v * self.t_step)
        self.pos.append(new_x)

        #new_v = ((self.t_step / self.part_mass) * (self.get_F() + \
        #                                          drag_force(self.rho_fluid, self.d_part, self.v_apparent, self.mu_fluid, self.cur_v) + \
        #                                          electric_other(self.q, self.cur_c, self.cur_x, self.H))) + \
        #                                        self.cur_v
        new_v = ((self.t_step / self.part_mass)  * (self.get_F() + \
                                                  drag_force(self.rho_fluid, self.d_part, self.v_apparent, self.mu_fluid, self.cur_v) + \
                                                  electric_other(self.q, self.cur_c, self.cur_x, self.H))) + \
                                                self.cur_v
        self.velo.append(new_v)

        new_a = (1 / self.part_mass) * (self.get_F() + \
                                                  drag_force(self.rho_fluid, self.d_part, -new_v, self.mu_fluid, new_v) + \
                                                  electric_other(self.q, self.cur_c, new_x, self.H))
        self.accel.append(new_a)

        # For now, just model the concentration as going down a little bit every time time step.
        # We will find a more accurate version of this later.
        new_c = self.cur_c - 0.000005
        self.concen.append(new_c)

        new_F = (self.get_F() + drag_force(self.rho_fluid, self.d_part, -new_v, self.mu_fluid, new_v) + \
                                electric_other(self.q, self.cur_c, new_x, self.H))
        self.force.append(new_F)

        # Now update each of the current values
        self.cur_a = new_a
        self.cur_v = new_v
        self.cur_x = new_x
        self.cur_c = new_c
        self.cur_force = new_F
        self.v_apparent = - self.cur_v

        # Increment the time tracker
        self.t += self.t_step

    def print_stats(self):
        print("Current position is: " + str(self.cur_x))
        print("Current velocity is: " + str(self.cur_v))
        print("Current acceleration is: " + str(self.cur_a))
        print("Current concentration is: " + str(self.cur_c))
        print("Current force is: " + str(self.cur_force))
        print("Current time is: " + str(self.t))

def test1():
    test_system = SmogCollector()
    test_system.print_forces()
    test_system.print_stats()
    print()
    n=100000
    i = 1

    while test_system.cur_x > 0:
        test_system.step_time()

        if i % 1000 == 0:
            test_system.print_forces()
            test_system.print_stats()
            print()
        
        #input("Hit enter to step")
        i += 1

    test_system.print_stats()
    plt.plot(test_system.pos)
    plt.show()

def test2():
    t_tracker = []
    forces = []

    for i in range(0, 200):
        cur_system = SmogCollector(_x_init=i)

        while cur_system.cur_x > 0:
            cur_system.step_time()

        t = cur_system.t
        t_tracker.append(t)
        f = cur_system.cur_force
        forces.append(f)
        print(t)
        print(f)

    plt.plot(forces)
    plt.show()

def main():
    #test1()
    test2()

if __name__ == "__main__":
    main()
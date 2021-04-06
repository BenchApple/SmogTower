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
# Density found here: https://pubmed.ncbi.nlm.nih.gov/26151089/
# TODO check to make sure my accel and velo readings are right because they seem a little too high rn. 
class SmogCollector:
    def __init__(self, _rho_part=1650, _d_part=0.0000025, _rho_fluid=1.225, _mu_fluid=0.000001803, \
                 _plate_voltage=20000, _plate_dist=0.05, _H=500, _x_init=100, _v_init=-0.01, _a_init=-0.01, \
                 _concentration=150, _t_step=0.002, _q=1.906*math.pow(10,-6)):
        e_0 = 8.854187 * math.pow(10, -12)

        # I have literally no clue how to find the charge of a particle so im gonna use this
        self.q = _q

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
        print(self.part_mass)

        # TODO calculate the total mass of the particles we;re dealing with.
        self.part_per_vol = (_concentration * math.pow(10, -6)) / self.part_mass
        print(self.part_per_vol)

        self.total_particles = self.part_per_vol * self.x_init
        print(self.total_particles)
        self.cur_particles = self.total_particles

        self.total_mass = self.total_particles * self.part_mass
        print(self.total_mass)
        self.cur_mass = self.total_mass

        #self.part_mass = _rho_part * _concentration * math.pow(10, -6) * _x_init

        self.cur_force = self.part_mass * self.cur_a

        # Each of the following lists track the accel, velocity, position, and concentration over time
        self.accel = [self.cur_a]
        self.velo = [self.cur_v]
        self.pos = [self.cur_x]
        self.concen = [self.cur_c]
        self.force = [self.cur_force]
        self.masses = [self.cur_mass]

        # t stores what the current time is.
        self.t = 0

    # Gets the value of F where F is defined as the bouyancy force minus the gravitational 
    # minus the force from the electric field
    def get_F(self):
        return (bouyant_force(self.rho_part, self.d_part, self.cur_particles) - grav_force(self.cur_mass) - electric_collector(self.q, self.sigma))
    
    # Prints all of the forces out.
    def print_forces(self):
        # This basically just tests the forces file.
        print("Grav")
        print(-grav_force(self.cur_mass))
        print("Bouyant")
        print(bouyant_force(self.rho_part, self.d_part, self.cur_particles))
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
        print(self.cur_mass)


    # Iterates time by one step, updating all of the elements in the process.
    def step_time(self):
        new_x = self.cur_x + (self.cur_v * self.t_step)
        self.pos.append(new_x)

        #new_v = ((self.t_step / self.part_mass) * (self.get_F() + \
        #                                          drag_force(self.rho_fluid, self.d_part, self.v_apparent, self.mu_fluid, self.cur_v) + \
        #                                          electric_other(self.q, self.cur_c, self.cur_x, self.H))) + \
        #                                        self.cur_v
        new_v = ((self.t_step / self.total_mass)  * (self.get_F() + \
                                                  drag_force(self.rho_fluid, self.d_part, self.v_apparent, self.mu_fluid, self.cur_v) + \
                                                  electric_other(self.q, self.cur_c, self.cur_x, self.H))) + \
                                                self.cur_v
        self.velo.append(new_v)

        new_a = (1 / self.total_mass) * (self.get_F() + \
                                                  drag_force(self.rho_fluid, self.d_part, -new_v, self.mu_fluid, new_v) + \
                                                  electric_other(self.q, self.cur_c, new_x, self.H))
        self.accel.append(new_a)

        # For now, just model the concentration as going down a little bit every time time step.
        # We will find a more accurate version of this later.
        new_c = self.c_init * (new_x / self.x_init)
        self.concen.append(new_c)

        new_F = (self.get_F() + drag_force(self.rho_fluid, self.d_part, -new_v, self.mu_fluid, new_v) + \
                                electric_other(self.q, self.cur_c, new_x, self.H))
        self.force.append(new_F)

        total_particles = self.part_per_vol * new_x
        self.cur_particles = total_particles

        total_mass = self.cur_particles * self.part_mass
        self.cur_mass = total_mass

        # Now update each of the current values
        self.cur_a = new_a
        self.cur_v = new_v
        self.cur_x = new_x
        self.cur_c = new_c
        self.cur_force = new_F
        self.v_apparent = - self.cur_v

        # Increment the time tracker
        self.t += self.t_step

    def run_simulation(self):
        while self.cur_x > 0:
            self.step_time()

        return self.t

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

        if i % 100 == 0:
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

    plt.plot(t_tracker)
    plt.show()

def testing():
    test = SmogCollector()

def cotters():
    dia_lo = 0.0000025
    dia_hi = 0.000001

    volt_lo = 100000
    volt_hi = 1000000

    dist_lo = 0.01
    dist_hi = 0.18

    concen_lo = 20
    concen_hi = 200

    q_lo = math.pow(10, -10)
    q_hi = math.pow(10, -5)

    trial1 = SmogCollector(_d_part=dia_lo, _plate_voltage=volt_lo, _plate_dist=dist_lo, _concentration=concen_lo, _q=q_lo)
    trial2 = SmogCollector(_d_part=dia_hi, _plate_voltage=volt_lo, _plate_dist=dist_lo, _concentration=concen_lo, _q=q_lo)
    trial3 = SmogCollector(_d_part=dia_lo, _plate_voltage=volt_hi, _plate_dist=dist_lo, _concentration=concen_lo, _q=q_lo)
    trial4 = SmogCollector(_d_part=dia_lo, _plate_voltage=volt_lo, _plate_dist=dist_hi, _concentration=concen_lo, _q=q_lo)
    trial5 = SmogCollector(_d_part=dia_lo, _plate_voltage=volt_lo, _plate_dist=dist_lo, _concentration=concen_hi, _q=q_lo)
    trial6 = SmogCollector(_d_part=dia_lo, _plate_voltage=volt_lo, _plate_dist=dist_lo, _concentration=concen_lo, _q=q_hi)
    trial7 = SmogCollector(_d_part=dia_lo, _plate_voltage=volt_hi, _plate_dist=dist_hi, _concentration=concen_hi, _q=q_hi)
    trial8 = SmogCollector(_d_part=dia_hi, _plate_voltage=volt_lo, _plate_dist=dist_hi, _concentration=concen_hi, _q=q_hi)
    trial9 = SmogCollector(_d_part=dia_hi, _plate_voltage=volt_hi, _plate_dist=dist_lo, _concentration=concen_hi, _q=q_hi)
    trial10 = SmogCollector(_d_part=dia_hi, _plate_voltage=volt_hi, _plate_dist=dist_hi, _concentration=concen_lo, _q=q_hi)
    trial11 = SmogCollector(_d_part=dia_hi, _plate_voltage=volt_hi, _plate_dist=dist_hi, _concentration=concen_hi, _q=q_lo)
    trial12 = SmogCollector(_d_part=dia_hi, _plate_voltage=volt_hi, _plate_dist=dist_hi, _concentration=concen_hi, _q=q_hi)

    t1 = trial1.run_simulation()
    t2 = trial2.run_simulation()
    t3 = trial3.run_simulation()
    t4 = trial4.run_simulation()
    t5 = trial5.run_simulation()
    t6 = trial6.run_simulation()
    t7 = trial7.run_simulation()
    t8 = trial8.run_simulation()
    t9 = trial9.run_simulation()
    t10 = trial10.run_simulation()
    t11 = trial11.run_simulation()
    t12 = trial12.run_simulation()

    print("Trial 1: " + str(t1))
    print("Trial 2: " + str(t2))
    print("Trial 3: " + str(t3))
    print("Trial 4: " + str(t4))
    print("Trial 5: " + str(t5))
    print("Trial 6: " + str(t6))
    print("Trial 7: " + str(t7))
    print("Trial 8: " + str(t8))
    print("Trial 9: " + str(t9))
    print("Trial 10: " + str(t10))
    print("Trial 11: " + str(t11))
    print("Trial 12: " + str(t12))

def graphing():
    # Graph the plate voltage across the different x values and the resultant time.
    x_vals = [25, 50, 100, 150, 200, 300]
    plate_voltages = [i for i in range(10000, 100001, 1000)]
    time_vals = []
    q_default = math.pow(10, -7)

    for x in x_vals:
        for v in plate_voltages:
            cur_trial = SmogCollector(_plate_voltage=v, _x_init=x, _q=q_default)
            cur_t = cur_trial.run_simulation()
            time_vals.append(cur_t)

        plt.plot(plate_voltages, time_vals)
        plt.ylabel("Time to Collect all Particles (seconds)")
        plt.xlabel("Voltage Running Through Base Plate (Volts)")
        plt.title("Voltage Running through Base Plate vs Time to Collect Particles at " + str(x) + " meters above Collector Plate (Volts vs seconds)")
        plt.show()

        # Wipe the time values
        time_vals = []

def graphing_q():
    # Graph the plate voltage across the different x values and the resultant time.
    x_vals = [25, 50, 100, 150, 200, 300]
    time_vals = []
    v_default = 30000

    # Create the q values
    q_vals = []
    for n in range(15, 5 - 1, -1):
        q_vals.append(1 * math.pow(10, -n))
        q_vals.append(2.5 * math.pow(10, -n))
        q_vals.append(5 * math.pow(10, -n))
        q_vals.append(7.5 * math.pow(10, -n))

    for x in x_vals:
        for q in q_vals:
            cur_trial = SmogCollector(_plate_voltage=v_default, _x_init=x, _q=q)
            cur_t = cur_trial.run_simulation()
            time_vals.append(cur_t)

        plt.semilogx(q_vals, time_vals)
        plt.ylabel("Time to Collect all Particles (seconds)")
        plt.xlabel("Ionization Charge of Smog Particles (C)")
        plt.title("Ionization Charge of Smog Particles vs Time to Collect Particles at " + str(x) + " meters above Collector Plate (Coulombs vs seconds)")
        plt.show()

        # Wipe the time values
        time_vals = []
    

def main():
    #testing()
    #test1()
    #test2()
    #cotters()
    #graphing()
    graphing_q()

if __name__ == "__main__":
    main()

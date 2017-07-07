# -*- coding: utf-8 -*-
"""
Created on Fri Dec  2 16:33:35 2016

@author: Anton
"""

#!usr/bin/python

import numpy
import matplotlib.pyplot as plt
from time import time

t0 = time()

def initial_parameters(total_seconds, time_per_iteration):
    global total_time, time_increment, iterations_total
    total_time = total_seconds
    time_increment = time_per_iteration
    iterations_total = total_time / time_increment

Grav_const = 6.6740831e-11

def angle_between(x1,y1,x2,y2):
    global angle1_2
    float(x1),float(x2),float(y1),float(y2)
    if x2>x1 and y2>y1:#comparatively Q1
        angle1_2 = (numpy.arctan((y2-y1)/(x2-x1)))
    elif x2<x1 and y2>y1:#comparatively Q2
        angle1_2 = numpy.pi - (numpy.arctan((y2-y1)/(x1-x2)))
    elif x2<x1 and y2<y1:#comparatively Q3
        angle1_2 = -(numpy.pi - (numpy.arctan((y1-y2)/(x1-x2))))
    elif x2>x1 and y2<y1:#comparatively Q4
        angle1_2 = -(numpy.arctan((y1-y2)/(x2-x1)))
    elif y1==y2 and x1>x2:
        angle1_2 = numpy.pi
    elif y1==y2 and x2>x1:
        angle1_2 = 0
    elif x1==x2 and y1>y2:
        angle1_2 = -numpy.pi/2
    elif x1==x2 and y2>y1:
        angle1_2 = numpy.pi/2
    else:
        print("theres been a problem in 'angle_between'")
    return angle1_2

def sq_distance_between(x1,y1,x2,y2):
    global sdb
    sdb = (x1-x2)**2 + (y1-y2)**2
    return sdb

class Particle(object):
    def __init__(self, start_x, start_y, mass, initial_vel, vel_theta, model_gravity = 0, static = 0):
        self.mass = mass
        self.iv_x = initial_vel * numpy.cos((vel_theta)*(numpy.pi/180))
        self.iv_y = initial_vel * numpy.sin((vel_theta)*(numpy.pi/180))
        self.ipos_x = start_x
        self.ipos_y = start_y
        self.gravity = model_gravity
        self.static = static
        self.acc_list_x = [0]*int(round(iterations_total,1))
        self.acc_list_y = [0]*int(round(iterations_total,1))
    def add_forces(self, *args):
        self.args2 = [list(i) for i in args]
        self.args2 = list([item for sublist in args for item in sublist])
        n = 0
        self.F_dict = {}
        for k in range(int(len(self.args2)/4)):
            for self.F_vars in args:
                self.F_time_start=self.F_vars[n]#timestart
                self.F_time_end=self.F_vars[n+1]#timeend
                self.F_mag=self.F_vars[n+2]#magnitudeofforcevector
                self.F_theta=self.F_vars[n+3]#angleofforcevector
                self.F_mag_x=self.F_mag * numpy.cos((self.F_theta)*(numpy.pi/180))#resolving
                self.F_mag_y=self.F_mag * numpy.sin((self.F_theta)*(numpy.pi/180))#()
                self.F_dict.update({n:self.F_time_start, n+1:self.F_time_end, n+2:self.F_mag_x, n+3:self.F_mag_y})#package dict
                if len(self.F_dict)%4 != 0:
                    print("Warning: check number of arguments for ADDEDFORCES")
                if self.F_time_start>=self.F_time_end:
                    print("Warning: check ADDEDFORCES for inconsistencies")
                n = n + 4
        return(self.F_dict)
    def instantiate_static_forces(self):
        o = 0
        n = 0
        for k in range(int(len(self.F_dict)/4)):#lets all added forces get 'added'(accelerationloops)
            self.F_start =self.F_dict[n]
            for m in range(int(round((self.F_dict[n+1]-self.F_dict[n])/time_increment,1))):#fills up a list with 1 'added force' at a time
                self.acc_list_x[(o+int(round(self.F_start/time_increment)))] = self.acc_list_x[(o+int(round(self.F_start/time_increment)))] + self.F_dict[n+2]/self.mass
                self.acc_list_y[(o+int(round(self.F_start/time_increment)))] = self.acc_list_y[(o+int(round(self.F_start/time_increment)))] + self.F_dict[n+3]/self.mass
                o = o+1
            o = 0
            n = n+4
        return self.acc_list_x, self.acc_list_y
        
def simulate_particles(*args):
    const = Grav_const
    n = 0
    D_O_L = {}
    for arg in args:
        pos_list_x = []
        pos_list_y = []
        vel_list_x = []
        vel_list_y = []
        a_list_x = []
        a_list_y = []
        pos_list_x.append(arg.ipos_x)
        pos_list_y.append(arg.ipos_y)
        vel_list_x.append(arg.iv_x)
        vel_list_y.append(arg.iv_y)
        a_list_x = arg.acc_list_x
        a_list_y = arg.acc_list_y
        D_O_L.update({
                 n:pos_list_x, n+1:pos_list_y,
                 n+2:vel_list_x, n+3:vel_list_y,
                 n+4:a_list_x, n+5:a_list_y,
                 n+6:arg.mass
                    })
        n = n + 7
    n = 0
    k = 0
    for i in range(int(round(iterations_total,1))):#for amount of iterations i, where i is the iteration
        k = 0
        for n in range(len(args)):
            for j in range(len(args)-1):
                if 7*k == 7*n:
                    k = k+1
                sq_distance_between(D_O_L[7*n][i], D_O_L[7*n+1][i], D_O_L[7*k][i], D_O_L[7*k+1][i])
                angle_between(D_O_L[7*n][i], D_O_L[7*n+1][i], D_O_L[7*k][i], D_O_L[7*k+1][i])
                D_O_L[7*n+4][i] = (D_O_L[7*n+4][i] + numpy.cos(angle1_2) *(const * D_O_L[7*n+6] * D_O_L[7*k+6])/sdb)/D_O_L[7*n+6]#x-accel#
                D_O_L[7*n+5][i] = (D_O_L[7*n+5][i] + numpy.sin(angle1_2) *(const * D_O_L[7*n+6] * D_O_L[7*k+6])/sdb)/D_O_L[7*n+6]#y-accel#
                k = k+1#for each particle to be in relation to every other particle
            k = 0
        for t in range(len(args)):
            D_O_L[7*t].append(D_O_L[7*t][i] + (D_O_L[7*t+2][i]+D_O_L[7*t+4][i]*time_increment)*time_increment#x-position###[x(n+1) = xn + (u+at)t+0.5at^2]
                  +0.5*D_O_L[7*t+4][i]*time_increment**2)
            D_O_L[7*t+1].append(D_O_L[7*t+1][i] + (D_O_L[7*t+3][i]+D_O_L[7*t+5][i]*time_increment)*time_increment#y-position
                  +0.5*D_O_L[7*t+5][i]*time_increment**2)
            D_O_L[7*t+2].append(D_O_L[7*t+2][i]+D_O_L[7*t+4][i]*time_increment)#x-vel
            D_O_L[7*t+3].append(D_O_L[7*t+3][i]+D_O_L[7*t+5][i]*time_increment)#x-vel
        for b in range(len(args)):
            c = numpy.random.rand(3,1)
            plt.figure(1)
            plt.plot(D_O_L[7*b], D_O_L[7*b+1], color = c, linestyle = '-')
            plt.figure(2)
            plt.scatter(D_O_L[7*b], D_O_L[7*b+1], color = c, edgecolors = 'k', alpha = 0.5)
    plt.show()
    t1 = time()
    print("#####################Complete!#####################")
    print("Everything took " + str(t1-t0) + "s")

###########################################################
#initial_parameters(total_time, time_between_iterations, filename)
###########################################################
#Particle-n = Particle(start_x, start_y, mass, initial_vel, vel_theta)
###########################################################
#Particle-n.add_forces([time_start_force1,time_end_force1,magnitude_force1,angle_force1,time_start_force2,time_end_force2,magnitude_force2,angle_force2,...n])
#Particle-n.instantiate_static_forces()
###########################################################
#!!!Computations can take time, depending on number of particles, and time between iterations(accuracy of simulation)!!!#
#!!!Accuracy of simulation greatly increased by adjusting 'time_between_iterations'!!!#
###########################################################
###################Example#################################
'''
initial_parameters(100000,1000)
obj1 = Particle(0,0,2900000000000000000000000,0,0)
obj3 = Particle(10000000,0,5000,5,90)
obj2 = Particle(5000000,0,500,50,-45)
obj2.add_forces([0,10000,10,-90,5000,10000,10000000,0,5000,10000,1000000,90])
obj2.instantiate_static_forces()
simulate_particles(obj1, obj2, obj3)
'''
###########################################################
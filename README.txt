Overview:
2D Particle motion simulator based on Euler's method with options for gravity between particles and extra forces acting on particles.
Supports any number of particles with any number of extra forces on them.
Output styles include scatter graph of positions of particles and/or file of lists of position, velocities, accelerations and masses for each particle.

Functions:
#####################INITIAL-PARAMETERS#####################(required)
#initial_parameters(total_time, time_between_iterations, filename)
#####################OUTPUT-STYLE###########################(required)
#output(graph = True/False, outputfile = 'outputfile.txt'(if needed))
#####################PARTICLE-CREATION######################(required)
#Particle-n = Particle(start_x, start_y, mass, initial_vel, vel_theta)
#####################EXTRA-ADDED-FORCES######################(not required)
#Particle-n.add_forces([time_start_force1,time_end_force1,magnitude_force1,angle_force1,time_start_force2,time_end_force2,magnitude_force2,angle_force2,...n])
#####################RUN-THE-SIMULATION######################(required)
#simulate_particles(particle-1,particle-2,particle-3,...,particle-n)

Accuracy:
Accuracy of simulation can be adjusted by changing 'time_between_iterations'.

For many iterations:
#######################!!!WARNINGS!!!########################
#!!!Computations can take time, depending on number of particles, and time between iterations(accuracy of simulation)!!!
#!!!Accuracy of simulation greatly increased by adjusting 'time_between_iterations'!!!
#!!!output files can be huge, >0.5Gb depending on accuracy of simulation!!!
#!!!DO NOT OUTPUT AS GRAPH IF MANY ITERATIONS ARE NEEDED, AS THIS TAKES HUGE AMOUNTS OF MEMORY!!!
#For an n number of particles, algorithm time complexity = O(n^2).
#############################################################

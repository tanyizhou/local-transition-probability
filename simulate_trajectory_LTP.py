# generate Brownian motion trajectories in int for LTP analysis.
import random
import math
import csv
import os


def bias_force(x):
    return -1/10


def stable_doublewell_force(x):
    return -1/40 - x/450 + 0.15*math.sin(x/10)


trajectory_file = "C:/Users/tanyi/PycharmProjects/test1/bias_force.csv"  # generated new trajectories

n = 500  # particle numbers
D = 0.25*20*20  # diffusion coefficient
channel_length = [-60, 60]
dt = 0.001  # frame time
dt_per_frame = round(1.0/dt)
ita = math.sqrt(2*D*dt)

if os.path.exists(trajectory_file):
    os.remove(trajectory_file)
    print("File deleted")
else:
    print("The file does not exist")


for j in range(n):
    trajectory = []
    i = 0
    x = random.uniform(channel_length[0], channel_length[1])
    while channel_length[0] <= x <= channel_length[1]:
        delta_x = bias_force(x) * D * dt + ita * random.normalvariate(0, 1)  # Langevin equation
        x += delta_x
        i += 1
        trajectory.append([j, i, round(x)])

    with open(trajectory_file, "a", newline='') as my_csv:
        csvWriter = csv.writer(my_csv)
        csvWriter.writerows(trajectory)


import matplotlib.pyplot as plt
import numpy as np


mp_strong = "strong"
mpi_strong = "mpi_strong"
mpi_basic_strong = "mpi_basic_strong"

mp_weak = "weak"
mpi_weak = "mpi_weak"
mpi_basic_weak = "mpi_basic_weak"


strong = [28.5422, 23.7999, 20.5307, 17.8587, 15.9314, 14.3476, 13.1119, 11.8945]
mpi_strong = [
    2.8083,
    1.24169,
    0.480095,
    0.310198,
    0.231433,
    0.196823,
    0.160241,
    0.140224,
]

mpi_basic_strong = [
    30.6274,
    15.6264,
    8.59633,
    6.27616,
    4.3309,
    3.16183,
    2.20938,
    1.76609,
]

weak = [1.22172, 2.23398, 13.0045, 20.2902, 28.4865, 38.9071, 58.5644, 71.7907]

mpi_weak = [
    0.062465,
    0.098162,
    0.138647,
    0.182928,
    0.230239,
    0.287746,
    0.347935,
    0.41029,
]

mpi_basic_weak = [
    0.119333,
    0.276746,
    2.05449,
    3.50382,
    4.49671,
    5.10296,
    5.57998,
    6.73858,
]

# Note that even in the OO-style, we use `.pyplot.figure` to create the figure.


def plot_strong():
    fig, ax = plt.subplots()  # Create a figure and an axes.
    ax.plot(
        [i for i in range(1, 9)], strong, label="OpenMP"
    )  # Plot some data on the axes.
    ax.plot(
        [i for i in range(1, 9)], mpi_basic_strong, label="MPI"
    )  # ... and some more.
    ax.plot([i for i in range(1, 9)], mpi_strong, label="MPI-Cannon")
    ax.set_xlabel("Compute Size")  # Add an x-label to the axes.
    ax.set_ylabel("Time")  # Add a y-label to the axes.
    ax.set_title("Strong Scaling on n=1800")  # Add a title to the axes.
    ax.legend()  # Add a legend.
    plt.show()


def plot_weak():
    fig, ax = plt.subplots()  # Create a figure and an axes.
    ax.plot(
        [i for i in range(600, 300 * 10, 300)], weak, label="OpenMP"
    )  # Plot some data on the axes.
    ax.plot(
        [i for i in range(600, 300 * 10, 300)], mpi_basic_weak, label="MPI"
    )  # ... and some more.
    ax.plot([i for i in range(600, 300 * 10, 300)], mpi_weak, label="MPI-Cannon")
    ax.set_xlabel("Grid Size")  # Add an x-label to the axes.
    ax.set_ylabel("Time")  # Add a y-label to the axes.
    ax.set_title("Weak Scaling")  # Add a title to the axes.
    ax.legend()  # Add a legend.
    plt.show()


plot_strong()
plot_weak()

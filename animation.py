"""
Script for simulating vortex interactions in 2D.

@author Yuliya Shpunarska
@collab Maya Tatarelli

10 Feb 2022
"""

import numpy as np
import matplotlib.pyplot as plt

dt = 10 # Adjust this
Nsteps = 100 # Adjust this

## Setting up initial conditions
# Vortex rings

y_v = np.array([-50,50,-50,50], dtype="f") # y positions of vortices
x_v = np.array([-50,-50,50,50], dtype="f") # x positions of vortices

k_v = np.array(50) # Line vortex constant

# Set orientations of vortices: +1 is out of the page, -1 is into the page
orientations = np.array([1, 1, -1, -1])

# Set up the plot
plt.ion()
fig, ax = plt.subplots( 1 , 1 )

fig.set_size_inches(10, 10)

# Mark position of vortices
p, = ax.plot(x_v, y_v, "k.", markersize=10)

ngrid = 200
res = 360j
increment = (2*ngrid)/res.imag # step size of grid
Y, X = np.mgrid[-ngrid:ngrid:res, -ngrid:ngrid:res]

# Velocities
vel_y = np.zeros(np.shape(Y))
vel_x = np.zeros(np.shape(X))

# masking radius
r_mask = 5

## Helper function to calculate velocity field

def calculate_velocity(vortex_x, vortex_y, x, y, r_mask=0, orientation=1):
    """
    This function calculates the x and y components of the velocity everywhere
    due to one line vortex.

    Parameters
    ----------
    vortex_x : float or int
        x-position of the vortex
    vortex_y : float or int
        x-position of the vortex
    x : array
        x-positions on the simulation grid as outputted by mgrid
    y : array
        y-positions on the simulation grid as outputted by mgrid
    r_mask : float or int
        radius around the vortex to mask (values of velocity within that radius
        will be nan)
    orientation : -1 or 1
        the vortex is oriented out of the page (1) or into the page (-1)

    """
    # Get position of this vortex
    r_vortex = np.sqrt(vortex_x**2+vortex_y**2) # distance between vortex and origin

    # Get distance between all points and vortex
    r = np.sqrt((x-vortex_x)**2 + (y-vortex_y)**2) # unmasked version
    r_masked = r.copy()
    r_masked[r<r_mask] = np.nan # Mask values within radius r_mask

    ## Compute speeds
    """
    Azimuthal velocity k_v/r can be decomposed into x and y components using
    the cos and sin of the angle formed (delta y/r and delta x/r respectively)
    """
    u_x_masked = -k_v * (y-vortex_y)/r_masked**2 * orientation
    u_y_masked = k_v * (x-vortex_x)/r_masked**2 * orientation

    u_x = -k_v/r * (y-vortex_y)/r * orientation
    u_y = k_v/r * (x-vortex_x)/r * orientation

    return {"vel_x" : u_x_masked,
           "vel_y" : u_y_masked,
           "vel_x_no_mask": u_x,
           "vel_y_no_mask": u_y}

## Evolution

count = 0
while count < Nsteps:

    # initialize velocities
    vel_x = np.zeros(np.shape(X))
    vel_y = np.zeros(np.shape(Y))

    # Store velocities at position of the vortices in this
    u_x_v = np.zeros(np.shape(x_v))
    u_y_v = np.zeros(np.shape(y_v))

    ## compute and update total velocity
    for i in range(len(x_v)):
        v = calculate_velocity(x_v[i], y_v[i], X, Y, orientation=orientations[i], r_mask=r_mask)
        vel_x += v["vel_x"]
        vel_y += v["vel_y"]

        for k in range(len(x_v)): ## find velocities affecting every other vortex
            if k != i:
                # get indices corresponding to kth vortex in the grid
                x_index = int(x_v[k]/increment)
                y_index = int(y_v[k]/increment)

                ## compute and update velocity components affecting each vortex
                if y_v[k] >= y_v[i]: # if vortex is higher than the main one
                    u_x_v[k] -= v["vel_x"][x_index, y_index]
                else: u_x_v[k] += v["vel_x"][x_index, y_index]

                if x_v[k] >= x_v[i]: # if vortex is to the right of main one
                    u_y_v[k] += v["vel_y"][x_index, y_index]
                else: u_y_v[k] -= v["vel_y"][x_index, y_index]

    # Plot streamlines
    ax.streamplot(X, Y, vel_x, vel_y, density=[1,1], color="b")

    # update positions of vortices
    y_v += u_y_v*dt
    x_v += u_x_v*dt

    p.set_xdata(x_v)
    p.set_ydata(y_v)

    # plot
    fig.canvas.draw()
    plt.pause(0.0001)

    # clear screen
    ax.collections = []
    ax.patches = []

    count += 1

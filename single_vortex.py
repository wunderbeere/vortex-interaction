"""
This is a script for testing the calculate_velocity function. Run this to see
the velocity field of a single vortex.
"""

import numpy as np
import matplotlib.pyplot as plt

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
    #r_vortex = np.sqrt(vortex_x**2+vortex_y**2) # distance between vortex and origin

    # Get distance between all points and vortex
    r = np.sqrt((x-vortex_x)**2 + (y-vortex_y)**2) # unmasked version
    r_masked = r.copy()
    r_masked[r<r_mask] = np.nan # Mask values within radius r_mask

    ## Compute speeds
    #u_phi = k_v/r * orientation # azimuthal velocity
    #u_phi_masked = k_v/r_masked * orientation # azimuthal velocity

    # Azimuthal velocity k_v/r can be decomposed into x and y components using
    # the cos and sin of the angle formed (delta y/r and delta x/r respectively)
    u_x_masked = -k_v * (y-vortex_y)/r_masked**2 * orientation
    u_y_masked = k_v * (x-vortex_x)/r_masked**2 * orientation

    u_x = k_v/r * (y-vortex_y)/r * orientation
    u_y = k_v/r * (x-vortex_x)/r * orientation

    ## decompose velocity into x and y components
    # find angles of velocities w.r.t. the x-axis
    #delta_x = x-vortex_x
    #alpha = np.arccos(delta_x/r)
    #gamma = np.pi/2 + alpha # this is the angle between u_phi and the x-axis

    # same but masked version
    #alpha_masked = np.arccos(delta_x/r_masked)
    #gamma_masked = np.pi/2 + alpha_masked # this is the angle between u_phi and the x-axis

    # get x and y components of velocity for this vortex
    #u_x = u_phi*np.cos(gamma)
    #u_y = u_phi*np.sin(gamma)

    #u_x_masked = u_phi*np.cos(gamma_masked)
    #u_y_masked = u_phi*np.sin(gamma_masked)

    return {"vel_x" : u_x_masked,
           "vel_y" : u_y_masked,
           "vel_x_no_mask": u_x,
           "vel_y_no_mask": u_y}


# Set up the plot
plt.ion()
fig, ax = plt.subplots( 1 , 1 )

fig.set_size_inches(10, 10)

ngrid = 200
res = 360j
increment = (2*ngrid)/res.imag # step size of grid
Y, X = np.mgrid[-ngrid:ngrid:res, -ngrid:ngrid:res]

# masking radius
r_mask = 5

k_v = 10
x_v = 50
y_v = 50
# Velocities
v = calculate_velocity(x_v, y_v, X, Y, r_mask=r_mask, orientation=-1)
vel_y = v["vel_y"]
vel_x = v["vel_x"]

# Mark position of vortices
p, = ax.plot(x_v, y_v, "k.", markersize=10)

ax.streamplot(X, Y, vel_x, vel_y, density=[1,1], color="b")

fig.canvas.draw()
plt.pause(20)
print("here")

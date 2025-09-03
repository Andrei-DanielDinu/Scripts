import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Parameters
k = 2 * np.pi / 5      # wave number
omega = 2 * np.pi / 20 # angular frequency
E0 = 1.0               # electric field amplitude
B0 = 1.0               # magnetic field amplitude
z_max = 20
num_points = 20

z = np.linspace(0, z_max, num_points)

def fields(z, t):
    """Return electric and magnetic field values at position z and time t."""
    phase = k * z - omega * t
    Ex = E0 * np.cos(phase)   # Electric field along x
    By = B0 * np.cos(phase)   # Magnetic field along y
    return Ex, By

# --- Plot setup ---
fig = plt.figure()
ax = fig.add_subplot(111, projection="3d")

def init():
    """Initialize the plot."""
    ax.set_xlim(-1.5, 1.5)
    ax.set_ylim(-1.5, 1.5)
    ax.set_zlim(0, z_max)
    ax.set_xlabel("X (Electric Field)")
    ax.set_ylabel("Y (Magnetic Field)")
    ax.set_zlabel("Z (Propagation)")
    ax.set_title("3D Electromagnetic Wave (Maxwell)")
    return []

def update(frame):
    """Update function for animation."""
    ax.cla()  # Clear and redraw everything

    # Axes & labels
    ax.set_xlim(-1.5, 1.5)
    ax.set_ylim(-1.5, 1.5)
    ax.set_zlim(0, z_max)
    ax.set_xlabel("X (Electric Field)")
    ax.set_ylabel("Y (Magnetic Field)")
    ax.set_zlabel("Z (Propagation)")
    ax.set_title("3D Electromagnetic Wave (Maxwell)")

    Ex, By = fields(z, frame)

    # Electric field (red arrows along x)
    ax.quiver(np.zeros(num_points), np.zeros(num_points), z,
              Ex, np.zeros(num_points), np.zeros(num_points),
              color="red", length=0.5, normalize=False)

    # Magnetic field (blue arrows along y)
    ax.quiver(np.zeros(num_points), np.zeros(num_points), z,
              np.zeros(num_points), By, np.zeros(num_points),
              color="blue", length=0.5, normalize=False)

    # Propagation direction (green arrow)
    ax.quiver(0, 0, 0, 0, 0, z_max, color="green", linewidth=2)

    return []

# --- Animate ---
ani = FuncAnimation(fig, update, frames=100, init_func=init,
                    interval=100, blit=False)
plt.show()

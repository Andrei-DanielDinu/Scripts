import numpy as np
import matplotlib.pyplot as plt

def magnetic_field_dipole(m, x0, y0, X, Y):
    """
    Compute the 2D magnetic field of a dipole located at (x0, y0) with moment m=(mx, my).
    Formula comes from magnetic dipole field equations in physics.
    """
    dx = X - x0
    dy = Y - y0
    r2 = dx**2 + dy**2
    r = np.sqrt(r2) + 1e-9  # avoid division by zero

    # unit vectors
    ex, ey = dx/r, dy/r

    # dot product m·r̂
    mdotr = m[0]*ex + m[1]*ey

    # Dipole field formula (2D simplification)
    Bx = (3*mdotr*ex - m[0]) / (r**3)
    By = (3*mdotr*ey - m[1]) / (r**3)
    return Bx, By

def plot_two_magnets(m1, pos1, m2, pos2, title):
    # Grid for field calculation
    x = np.linspace(-3, 3, 80)
    y = np.linspace(-3, 3, 80)
    X, Y = np.meshgrid(x, y)

    # Field from both dipoles
    Bx1, By1 = magnetic_field_dipole(m1, *pos1, X, Y)
    Bx2, By2 = magnetic_field_dipole(m2, *pos2, X, Y)
    Bx, By = Bx1 + Bx2, By1 + By2
    Bmag = np.sqrt(Bx**2 + By**2)

    # Plot
    plt.figure(figsize=(7,7))
    plt.streamplot(X, Y, Bx, By, color=Bmag, cmap="plasma", density=1.4, linewidth=1)
    plt.scatter(*pos1, color="red", s=200, marker="^", label="Magnet 1 (dipole)")
    plt.scatter(*pos2, color="blue", s=200, marker="v", label="Magnet 2 (dipole)")
    plt.title(title)
    plt.axis("equal")
    plt.legend()
    plt.show()

def main():
    print("🧲 Why do magnets attract or repel?")
    print("\nPhysics explanation:")
    print(" - A bar magnet can be thought of as a dipole (North-South).")
    print(" - Magnetic field lines always leave the North and enter the South.")
    print(" - If opposite poles face each other, the field lines connect -> Attraction.")
    print(" - If like poles face each other, the field lines clash -> Repulsion.\n")

    # Example 1: Opposite poles facing each other (attraction)
    m1 = (0, 1)   # dipole moment pointing up
    pos1 = (-1, 0)
    m2 = (0, -1)  # dipole moment pointing down
    pos2 = (1, 0)
    plot_two_magnets(m1, pos1, m2, pos2, "Opposite Poles → Attraction")

    # Example 2: Like poles facing each other (repulsion)
    m1 = (0, 1)   # both pointing up
    pos1 = (-1, 0)
    m2 = (0, 1)
    pos2 = (1, 0)
    plot_two_magnets(m1, pos1, m2, pos2, "Like Poles → Repulsion")

if __name__ == "__main__":
    main()

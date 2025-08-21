import numpy as np
import matplotlib.pyplot as plt

# Magnetic constant (vacuum permeability)
mu0 = 4 * np.pi * 1e-7  

def magnetic_field_dipole(m, x, y):
    """
    Magnetic dipole field in 2D plane (z=0).
    Approximates a bar magnet centered at origin along y-axis.
    """
    r2 = x**2 + y**2
    r = np.sqrt(r2)
    r5 = r**5 + 1e-20  # avoid division by zero

    # Magnetic dipole moment vector (pointing upwards)
    mx, my = 0, m  

    # Dipole field equations
    Bx = mu0/(4*np.pi) * (3*x*(mx*x + my*y)/r5 - mx/r**3)
    By = mu0/(4*np.pi) * (3*y*(mx*x + my*y)/r5 - my/r**3)
    return Bx, By


def main():
    print("=== Magnetic Dipole Field Simulation ===")
    m = float(input("Enter magnetic dipole strength (A·m², e.g. 1): "))

    # Grid
    x = np.linspace(-2, 2, 40)
    y = np.linspace(-2, 2, 40)
    X, Y = np.meshgrid(x, y)
    Bx, By = magnetic_field_dipole(m, X, Y)

    # Normalize for quiver plotting
    B_magnitude = np.sqrt(Bx**2 + By**2)
    Bx_norm = Bx / (B_magnitude + 1e-20)
    By_norm = By / (B_magnitude + 1e-20)

    # Plot
    plt.figure(figsize=(8, 8))
    plt.streamplot(X, Y, Bx, By, color=B_magnitude, cmap="plasma", density=2)
    plt.quiver(X, Y, Bx_norm, By_norm, color="black", alpha=0.5)
    plt.scatter(0, 0, c="red", s=200, label="Magnet center")
    plt.title("Magnetic Field of a Bar Magnet (Dipole Approximation)")
    plt.xlabel("x (m)")
    plt.ylabel("y (m)")
    plt.legend()
    plt.colorbar(label="Magnetic field strength (T, arbitrary scale)")
    plt.axis("equal")
    plt.show()

    # Explanations
    print("\n=== Explanation ===")
    print("This plot shows the magnetic field lines around a bar magnet.")
    print("- The magnet is placed at the origin (red dot).")
    print("- Field lines emerge from the North pole and curve back into the South pole.")
    print("- The arrows show direction, the color shows relative strength.")
    print("- Near the magnet, the field is strongest (yellow/white).")
    print("- Farther away, the field gets weaker (dark blue).")
    print("\nThis is why iron filings form curved patterns around a magnet!")

if __name__ == "__main__":
    main()

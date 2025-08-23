import numpy as np
import matplotlib.pyplot as plt

# Constants
mu0 = 4 * np.pi * 1e-7  # T·m/A

def magnetic_field_wire(I, X, Y):
    """Compute Bx, By for a straight wire along z-axis."""
    R = np.sqrt(X**2 + Y**2)
    B_magnitude = np.where(R != 0, mu0 * I / (2 * np.pi * R), 0)

    # Tangential direction (azimuthal around wire)
    Bx = -Y / (R + 1e-12) * B_magnitude
    By = X / (R + 1e-12) * B_magnitude
    return Bx, By

def main():
    print("🔵 2D Magnetic Field Around a Straight Wire 🔵")

    I = float(input("Enter current I in Amperes (e.g. 5): "))
    extent = float(input("Enter spatial extent in meters (e.g. 1.0 for -1 to 1 m): "))
    grid_size = int(input("Enter grid resolution (e.g. 20 for coarse, 40 for fine): "))
    step = int(input("Enter arrow downsampling step (higher = fewer arrows, e.g. 2): "))

    # Create 2D grid
    x = np.linspace(-extent, extent, grid_size)
    y = np.linspace(-extent, extent, grid_size)
    X, Y = np.meshgrid(x, y)

    # Compute magnetic field
    Bx, By = magnetic_field_wire(I, X, Y)
    B_mag = np.sqrt(Bx**2 + By**2)

    print("\n📘 Explanation:")
    print(f" - Current I = {I} A flows along the z-axis (perpendicular to this plane).")
    print(" - Magnetic field forms circular loops around the wire.")
    print(" - Field strength decreases with distance (1/r).")
    print(" - Arrows show field direction; color shows field magnitude.")

    # Plot 2D vector field
    plt.figure(figsize=(7,7))
    plt.streamplot(X, Y, Bx, By, color=B_mag, cmap="plasma", density=1.5, linewidth=1)
    plt.quiver(X[::step, ::step], Y[::step, ::step],
               Bx[::step, ::step], By[::step, ::step],
               B_mag[::step, ::step], cmap="plasma", scale=30)
    plt.colorbar(label="Magnetic Field Strength (T)")
    plt.title(f"Magnetic Field Around a Straight Wire (I={I} A)")
    plt.xlabel("X (m)")
    plt.ylabel("Y (m)")
    plt.axis("equal")
    plt.show()

if __name__ == "__main__":
    main()

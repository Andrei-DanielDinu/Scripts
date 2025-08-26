import numpy as np
import matplotlib.pyplot as plt

def main():
    print("⚡ Thomson's Cathode Ray Experiment ⚡\n")
    print("Goal: Measure the charge-to-mass ratio (e/m) of the electron.\n")
    
    print("Principle:")
    print(" - Electrons (cathode rays) are accelerated and pass through crossed E and B fields.")
    print(" - When the electric force equals the magnetic force, the beam is undeflected:")
    print("     eE = e v B  →  v = E / B")
    print(" - In magnetic field only, the electron beam curves in a circle:")
    print("     r = m v / (e B)")
    print(" - Combining gives:")
    print("     e/m = E / (B^2 r)\n")
    
    # Example setup
    V = 2000       # accelerating voltage (V)
    d = 0.01       # plate separation for electric field (m)
    E = V / d      # electric field (V/m)
    B = 5e-3       # magnetic field (Tesla)
    r = 0.05       # measured curvature radius (m)

    # Velocity of electrons when beam is balanced
    v = E / B
    # Charge-to-mass ratio
    e_over_m = E / (B**2 * r)

    print("Example Calculation:")
    print(f" - Accelerating Voltage V = {V} V")
    print(f" - Plate separation d = {d} m")
    print(f" - Electric field E = {E:.2e} V/m")
    print(f" - Magnetic field B = {B:.2e} T")
    print(f" - Radius of curvature r = {r:.2f} m")
    print(f" - Electron velocity v = {v:.2e} m/s")
    print(f" - Measured e/m ratio = {e_over_m:.2e} C/kg")
    
    print("\nKnown value: e/m ≈ 1.76×10^11 C/kg")
    
    # Visualization of electron path
    fig, ax = plt.subplots(figsize=(7,5))
    
    # Draw parallel plates
    ax.plot([0,0.2],[0.05,0.05],"k",linewidth=3)
    ax.plot([0,0.2],[-0.05,-0.05],"k",linewidth=3)
    ax.text(0.21,0.05,"+ Plate",va="center")
    ax.text(0.21,-0.05,"- Plate",va="center")

    # Electron initial path
    x = np.linspace(0,0.05,100)
    y = np.zeros_like(x)
    ax.plot(x,y,"b--",label="Without fields")

    # Curved trajectory (magnetic only)
    theta = np.linspace(0, np.pi/4, 100)
    x_curve = r*np.sin(theta)
    y_curve = r*(1-np.cos(theta))
    ax.plot(0.05+x_curve, y_curve, "r", label="Deflected path")

    # Mark electron
    ax.scatter([0.05],[0],s=80,color="blue",label="Electron")

    ax.set_aspect("equal")
    ax.set_xlim(0,0.2)
    ax.set_ylim(-0.02,0.08)
    ax.set_xlabel("Distance (m)")
    ax.set_ylabel("Height (m)")
    ax.set_title("Thomson's Cathode Ray Tube Experiment")
    ax.legend()
    plt.show()

if __name__ == "__main__":
    main()

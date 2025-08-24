import numpy as np
import matplotlib.pyplot as plt

def main():
    print("⚡ Millikan's Oil Drop Experiment ⚡")
    print("\nGoal: Measure the charge of a single electron (e).")
    print("\nPrinciple:")
    print(" - Tiny oil droplets are sprayed into a chamber.")
    print(" - Gravity pulls them down (Fg = mg).")
    print(" - Between two capacitor plates, an electric field E applies a force Fe = qE.")
    print(" - By adjusting the voltage, droplets can be suspended (Fg = Fe).")
    print(" - From this balance, q = mg / E.")
    print(" - Repeating the experiment shows that charges are always multiples of a smallest value: e ≈ 1.6×10^-19 C.\n")

    # Constants
    g = 9.81              # gravity (m/s^2)
    rho_oil = 860         # density of oil (kg/m^3)
    r = 1e-6              # droplet radius ~ 1 micrometer
    V = 200               # voltage applied across plates (V)
    d = 0.005             # distance between plates (5 mm)
    E = V/d               # electric field (V/m)

    # Mass of the droplet (sphere)
    volume = (4/3)*np.pi*r**3
    m = rho_oil * volume
    Fg = m*g

    # Suppose droplet is balanced (hovering)
    q = Fg / E

    print(f"Example Calculation:")
    print(f" - Droplet radius r = {r:.1e} m")
    print(f" - Droplet mass m = {m:.2e} kg")
    print(f" - Gravitational force Fg = {Fg:.2e} N")
    print(f" - Electric field E = {E:.2e} V/m")
    print(f" - Measured droplet charge q = {q:.2e} C")

    # Compare to electron charge
    e = 1.602e-19
    n = q / e

    print(f" - Charge multiples of e: q/e ≈ {n:.1f}")
    print("→ Meaning: the droplet’s charge is an integer multiple of the fundamental electron charge e.")

    # Simple visualization
    fig, ax = plt.subplots(figsize=(6,6))
    ax.set_xlim(-1, 1)
    ax.set_ylim(-2, 2)

    # Draw plates
    ax.plot([-1,1], [1,1], color="black", linewidth=3)
    ax.plot([-1,1], [-1,-1], color="black", linewidth=3)
    ax.text(1.1,1,"+ Plate", va="center", fontsize=10)
    ax.text(1.1,-1,"- Plate", va="center", fontsize=10)

    # Draw droplet
    ax.scatter(0, 0, s=200, color="orange", edgecolor="black", label="Oil droplet")

    # Force arrows
    ax.arrow(0, 0, 0, -0.8, head_width=0.05, head_length=0.1, color="red")
    ax.text(0.1, -0.4, "Fg", color="red")

    ax.arrow(0, 0, 0, 0.8, head_width=0.05, head_length=0.1, color="blue")
    ax.text(0.1, 0.4, "Fe", color="blue")

    ax.set_title("Millikan Oil Drop: Balance of Forces")
    ax.axis("off")
    plt.legend()
    plt.show()

if __name__ == "__main__":
    main()

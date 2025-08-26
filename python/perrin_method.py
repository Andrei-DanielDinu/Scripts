import numpy as np
import matplotlib.pyplot as plt

def main():
    print("🔬 Perrin’s Experiment and the Discovery of Avogadro’s Number 🔬\n")
    print("Principle:")
    print(" - Colloidal particles suspended in liquid settle under gravity.")
    print(" - Their concentration follows Boltzmann’s distribution law:")
    print("     n(z) = n(0) * exp(-m g z / (k_B T))")
    print(" - From the slope, we find Boltzmann’s constant k_B.")
    print(" - Then, using R = N_A * k_B, Avogadro’s number N_A is determined.\n")
    
    # Constants
    R = 8.314      # J/mol·K (Gas constant)
    T = 298        # K (room temperature)
    g = 9.81       # m/s²
    rho_particle = 1.05e3   # density of particle (kg/m³) ~ latex particle
    rho_fluid = 1.00e3      # density of water (kg/m³)
    r = 0.5e-6     # particle radius (m)

    # Effective mass of one particle
    volume = (4/3) * np.pi * r**3
    m = (rho_particle - rho_fluid) * volume

    print("Experimental Setup:")
    print(f" - Particle radius r = {r:.2e} m")
    print(f" - Density difference Δρ = {(rho_particle - rho_fluid):.1f} kg/m³")
    print(f" - Effective mass m = {m:.2e} kg\n")

    # Generate concentration profile
    z = np.linspace(0, 20e-6, 200)  # height (m)
    kB_real = 1.38e-23              # true Boltzmann constant
    n0 = 1.0                        # arbitrary concentration
    n = n0 * np.exp(-m*g*z/(kB_real*T))

    # Fit slope to extract k_B
    slope, _ = np.polyfit(z, np.log(n), 1)
    kB_measured = -m*g/slope
    N_A_measured = R / kB_measured

    print("Results from Simulation:")
    print(f" - True Boltzmann constant k_B = {kB_real:.3e} J/K")
    print(f" - Measured Boltzmann constant k_B = {kB_measured:.3e} J/K")
    print(f" - Avogadro’s number N_A = {N_A_measured:.3e} mol⁻¹")
    print("\nKnown value: N_A ≈ 6.022×10^23 mol⁻¹")

    # Visualization
    plt.figure(figsize=(7,5))
    plt.plot(z*1e6, n, "b", label="Particle concentration")
    plt.xlabel("Height above bottom (µm)")
    plt.ylabel("Relative concentration n(z)")
    plt.title("Perrin’s Method to Measure Avogadro’s Number")
    plt.legend()
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    main()

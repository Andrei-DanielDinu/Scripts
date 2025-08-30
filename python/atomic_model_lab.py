"""
Interactive educational simulator inspired by Rutherford-era atomic experiments
and later models (Bohr, quantum 1s orbital). It includes:
  • Rutherford scattering Monte Carlo + analytic Rutherford cross-section
  • Bohr model radius & energy for hydrogen-like atoms
  • 2D slice visualization of 1s electron probability density
  • Visual comparison of nucleus size vs atomic size
  • Explanatory text about electron spin and limitations of classical pictures

"""

import math
import numpy as np
import matplotlib.pyplot as plt

# ----------------------------
# Physical constants (SI)
# ----------------------------
e = 1.602176634e-19         # elementary charge C
epsilon0 = 8.8541878128e-12
k_e = 1 / (4 * math.pi * epsilon0)  # Coulomb constant ~ 8.9875e9 N·m^2/C^2
m_e = 9.1093837015e-31      # electron mass kg
h = 6.62607015e-34          # Planck constant J·s
hbar = h / (2*math.pi)
a0 = 5.29177210903e-11      # Bohr radius (m)
Rydberg = 13.605693122994   # eV (ionization energy of hydrogen)
amu = 1.66053906660e-27     # atomic mass unit kg
m_alpha = 4.001506 * amu    # alpha particle mass approx

# ----------------------------
# Utilities
# ----------------------------
def ev_to_joules(E_eV):
    return E_eV * e

def joules_to_ev(E_J):
    return E_J / e

# Rutherford scattering formula:
# for a particle with charge q1 (alpha = +2e) and target nucleus charge q2 = Z*e:
# scattering angle theta from an impact parameter b and energy E:
# theta(b) = 2 * arctan( (k * q1 * q2) / (2 * E * b) )
def scattering_angle_from_b(b, q1, q2, E_J):
    num = k_e * q1 * q2
    denom = 2.0 * E_J * b
    # avoid division by zero
    val = num / (denom + 1e-300)
    theta = 2.0 * np.arctan(val)
    return theta

# Rutherford differential cross-section:
# dσ/dΩ = ( (k*q1*q2) / (4E) )^2 * 1 / sin^4(theta/2)
def rutherford_differential(theta, q1, q2, E_J):
    pref = (k_e * q1 * q2) / (4.0 * E_J)
    denom = (math.sin(theta/2.0)**4 + 1e-300)
    return (pref**2) / denom

# Bohr model radius and energy for hydrogen-like atom (n, Z)
def bohr_radius(n, Z=1):
    return (n**2) * a0 / Z

def bohr_energy_ev(n, Z=1):
    # energy (negative) in eV
    return - Rydberg * (Z**2) / (n**2)

# 1s hydrogenic radial probability density (3D) -> we'll compute a 2D slice in xy plane:
# |psi_1s(r)|^2 = (1/pi) * (1/a0)^3 * exp(-2 r / a0)
def psi_1s_probability_density(x, y, z=0, Z=1):
    # For hydrogenic with nuclear charge Z effective a0/Z
    a = a0 / Z
    r = np.sqrt(x**2 + y**2 + z**2)
    pref = 1.0 / (math.pi * a**3)
    return pref * np.exp(-2.0 * r / a)

# ----------------------------
# Main interactive function
# ----------------------------
def main():
    print("\n=== Atomic Model Lab: Rutherford, Bohr, Quantum 1s & 'Spin' ===\n")
    print("This script will let you:")
    print(" 1) simulate Rutherford scattering (alpha particles vs nucleus)")
    print(" 2) compute Bohr radii & energies for H-like atoms")
    print(" 3) visualize a 2D slice of the 1s electron probability density")
    print(" 4) compare nucleus size vs atomic size and discuss spin\n")

    # ---- Rutherford scattering ----
    print("---- Rutherford scattering simulation ----")
    try:
        Z = int(input("Enter target nucleus atomic number Z (e.g. 79 for Au): ") or "79")
        E_keV = float(input("Alpha particle energy (keV) (e.g. 5000 for 5 MeV): ") or "5000")
        N_particles = int(input("Number of simulated alpha particles (e.g. 20000): ") or "20000")
        b_max = float(input("Max impact parameter b_max in fm (e.g. 1000): ") or "1000")  # fm
    except Exception:
        print("Invalid input — using defaults: Z=79, E=5000 keV, N=20000, b_max=1000 fm")
        Z = 79
        E_keV = 5000.0
        N_particles = 20000
        b_max = 1000.0

    E_J = ev_to_joules(E_keV * 1e3)  # convert keV to eV then to J
    q_alpha = 2.0 * e
    q_nucleus = Z * e

    # Sample impact parameters with uniform probability in area: p(b) db ∝ b db => b = sqrt(u) * b_max
    rng = np.random.default_rng(1)
    u = rng.random(N_particles)
    b_meters = (np.sqrt(u) * b_max) * 1e-15  # convert fm -> m

    theta = scattering_angle_from_b(b_meters, q_alpha, q_nucleus, E_J)  # in radians
    theta_deg = np.degrees(theta)

    # Build histogram of scattering angles
    bins = np.linspace(0, 180, 181)
    hist, _ = np.histogram(theta_deg, bins=bins)

    # Analytic Rutherford differential cross-section for a set of theta bins (midpoints)
    bin_centers = 0.5 * (bins[:-1] + bins[1:])
    theta_rad_centers = np.radians(bin_centers)
    dsdo = np.array([rutherford_differential(t, q_alpha, q_nucleus, E_J) for t in theta_rad_centers])

    # Normalize analytic curve for plotting overlay
    # convert dσ/dΩ to expected counts per bin ~ (dσ/dΩ) * ΔΩ ; ΔΩ ~ 2π * sinθ * Δθ for azimuthally symmetric measurement
    dtheta = np.radians(bins[1] - bins[0])
    delta_omega = 2 * math.pi * np.sin(theta_rad_centers) * dtheta
    expected_counts = dsdo * delta_omega
    # normalize to the number of simulated particles for plotting
    expected_counts_scaled = expected_counts / np.sum(expected_counts) * N_particles

    # Print informative numbers
    print("\nRutherford simulation summary:")
    print(f" - N particles simulated: {N_particles}")
    print(f" - Target Z = {Z}, alpha charge = +2e")
    print(f" - Alpha energy = {E_keV:.1f} keV ({E_J:.3e} J)")
    print(f" - Impact parameter range: 0 .. {b_max:.1f} fm")
    print(" - Observations: Rutherford predicted many small-angle scatterings but rare large-angle events.\n")

    # Plots for Rutherford
    plt.figure(figsize=(10,4))
    plt.subplot(1,2,1)
    plt.title("Rutherford scattering: simulated angle histogram")
    plt.bar(bin_centers, hist, width=bins[1]-bins[0], color="C0", alpha=0.7, label="Monte Carlo")
    plt.plot(bin_centers, expected_counts_scaled, 'r-', lw=2, label="Rutherford analytic (scaled)")
    plt.xlabel("Scattering angle (deg)")
    plt.ylabel("Counts")
    plt.legend()
    plt.grid(True)

    # Show polar scatter (angles)
    plt.subplot(1,2,2, projection='polar')
    plt.title("Polar scatter (angles) - sample of particles")
    sample_idx = rng.choice(N_particles, size=min(800, N_particles), replace=False)
    # To plot, need radius positive; use min(1, scaled) to keep visualization
    r_plot = np.ones_like(sample_idx)  # unit radius for each point
    plt.scatter(np.radians(theta_deg[sample_idx]), r_plot, s=4, alpha=0.6)
    plt.yticks([])
    plt.show()

    # ---- Bohr model and atom size ----
    print("---- Bohr model & atom size estimates ----")
    try:
        Z_bohr = int(input("Enter Z for Bohr-like hydrogenic atom (e.g. 1 = H, 2 = He+): ") or "1")
        n_level = int(input("Enter principal quantum number n (e.g. 1): ") or "1")
    except Exception:
        Z_bohr = 1
        n_level = 1

    r_bohr = bohr_radius(n_level, Z=Z_bohr)
    E_bohr = bohr_energy_ev(n_level, Z=Z_bohr)

    print(f"Bohr radius for n={n_level}, Z={Z_bohr}: r = {r_bohr:.3e} m ({r_bohr*1e10:.3f} Å)")
    print(f"Bohr energy level: E = {E_bohr:.3f} eV (negative => bound)\n")

    # Typical sizes
    atom_size = 1e-10   # ~ 1 Å
    nucleus_size = 1e-15  # ~ 1 fm
    print("Scale comparison (typical):")
    print(f" - Typical atomic radius ~ {atom_size:.1e} m (Å scale)")
    print(f" - Typical nucleus radius ~ {nucleus_size:.1e} m (fm scale)")
    print(" -> Nucleus is ~10^5 times smaller in linear scale than atom; most of atom is empty space.\n")

    # ---- Quantum 1s probability density visualization (2D slice) ----
    print("---- 1s electron probability density (2D slice) ----")
    try:
        plot_extent_A = float(input("Enter half-width for 1s density plot in Å (e.g. 2.0): ") or "2.0")
        resolution = int(input("Grid resolution (e.g. 300): ") or "300")
    except Exception:
        plot_extent_A = 2.0
        resolution = 300

    extent_m = plot_extent_A * 1e-10
    x = np.linspace(-extent_m, extent_m, resolution)
    y = np.linspace(-extent_m, extent_m, resolution)
    X, Y = np.meshgrid(x, y, indexing='xy')
    Z0 = 0.0
    prob = psi_1s_probability_density(X, Y, Z0, Z=Z_bohr)  # Z param here is nuclear charge number for hydrogenic

    # Normalize probability for plotting (not necessary, but better colors)
    prob_norm = prob / prob.max()

    plt.figure(figsize=(6,5))
    im = plt.imshow(prob_norm.T, extent=[-plot_extent_A, plot_extent_A, -plot_extent_A, plot_extent_A],
                    origin='lower', cmap='inferno', norm=None)
    plt.colorbar(im, label='Normalized |ψ|^2 (arbitrary units)')
    plt.title(f"1s probability density (slice z=0), Z={Z_bohr}")
    plt.xlabel("x (Å)")
    plt.ylabel("y (Å)")
    # Overplot Bohr radius for reference
    circle_r_A = (a0 / Z_bohr) * 1e10
    theta = np.linspace(0, 2*np.pi, 400)
    plt.contour((X*1e10), (Y*1e10), prob_norm.T, levels=[0.1, 0.3, 0.5], colors=['white','cyan','lime'])
    plt.gca().add_patch(plt.Circle((0,0), circle_r_A, color='white', fill=False, lw=1.2, alpha=0.7, linestyle='--'))
    plt.text(0.5*plot_extent_A, -plot_extent_A*0.9, f"Bohr radius a0/Z ≈ {circle_r_A:.2f} Å", color='white')
    plt.show()

    # ---- Simple "electron motion" classical toy and spin discussion ----
    print("---- Electrons, motion, and spin (conceptual) ----\n")
    print("Classical 'orbit' picture (Bohr) vs quantum reality:")
    print(" - Bohr model: electrons in circular orbits with quantized angular momentum; gives correct H energies.")
    print(" - Quantum model: electrons are wavefunctions (probability clouds). No definite orbit; only probability densities.")
    print("\nSpin & magnetic moment (qualitative):")
    print(" - 'Spin' is an intrinsic quantum property; it is NOT a tiny classical spinning ball.")
    print(" - Spin quantum number s=1/2 for electrons; spin projection m_s = ±1/2.")
    print(" - Spin produces a magnetic moment μ = -g (e/2m) S (quantum operator). In practice, this gives fine structure, Zeeman splitting.")
    print(" - In experiments (Stern–Gerlach), spin causes beam splitting in magnetic field gradients.\n")

    # Provide a tiny visualization of classical vs quantum size scale
    fig, ax = plt.subplots(figsize=(6,4))
    ax.set_title("Size scales: nucleus vs atom")
    ax.scatter(0, 0, s=200, color='red', label='nucleus (fm scale)')
    ax.scatter(0.5, 0, s=20, color='blue', label='electron cloud (Å scale)')
    ax.set_xlim(-1, 2)
    ax.set_ylim(-1, 1)
    ax.axis('off')
    ax.legend(loc='upper right')
    plt.show()

    print("---- Final remarks ----")
    print(" • Rutherford's scattering showed a tiny, dense nucleus; alpha particles sometimes scatter at large angles.")
    print(" • Bohr explained hydrogen spectral lines with quantized orbits; quantum mechanics replaced orbits with wavefunctions.")
    print(" • Spin is intrinsically quantum and gives rise to magnetic properties and statistics (Pauli exclusion).")
    print("\nYou can re-run sections by changing inputs. Try varying Z, energy, and see how rare large-angle scattering becomes!")

if __name__ == "__main__":
    main()

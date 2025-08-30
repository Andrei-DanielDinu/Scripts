"""
atomic_extensions.py

Menu-driven script implementing:
 1) Monte Carlo Rutherford scattering + fit for Z
 2) Animate 1s -> n hydrogenic orbital changes (probability density)
 3) Stern-Gerlach spin-1/2 beam demo (animation + histogram)
 4) Spin precession (Bloch sphere Larmor precession)

Dependencies: numpy, matplotlib
Run: python atomic_extensions.py
"""

import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

# Physical constants
e = 1.602176634e-19
epsilon0 = 8.8541878128e-12
k_e = 1 / (4 * math.pi * epsilon0)
hbar = 1.054571817e-34
m_e = 9.1093837015e-31
a0 = 5.29177210903e-11
Rydberg = 13.605693122994  # eV

# ----------------------------
# 1) Rutherford Monte Carlo + fit for Z
# ----------------------------
def rutherford_angle_from_b(b, q1, q2, E_J):
    val = (k_e * q1 * q2) / (2.0 * E_J * b + 1e-300)
    return 2.0 * np.arctan(val)

def simulate_rutherford(Z=79, E_keV=5000.0, N=15000, b_max_fm=1000.0, noise_frac=0.05, seed=1):
    """
    Simulate N alpha particles with impact parameter distribution uniform in area up to b_max (fm).
    Returns histogram counts (per degree bin) and analytic expectation (scaled).
    """
    rng = np.random.default_rng(seed)
    q_alpha = 2*e
    q_nucleus = Z * e
    E_J = (E_keV * 1e3) * e

    # sample b ~ sqrt(u) * b_max
    u = rng.random(N)
    b_m = np.sqrt(u) * b_max_fm * 1e-15
    theta = rutherford_angle_from_b(b_m, q_alpha, q_nucleus, E_J)  # radians
    theta_deg = np.degrees(theta)
    # histogram
    bins = np.linspace(0.0, 180.0, 181)
    hist, _ = np.histogram(theta_deg, bins=bins)

    # add multiplicative noise to histogram counts (simulate measurement noise)
    hist_noisy = hist * (1.0 + rng.normal(0, noise_frac, size=hist.size))
    hist_noisy = np.maximum(hist_noisy, 0.0)

    # analytic differential function (unnormalized): g(theta) = 1/sin^4(theta/2)
    bin_centers = 0.5 * (bins[:-1] + bins[1:])
    theta_r = np.radians(bin_centers)
    g = 1.0 / (np.sin(theta_r / 2.0)**4 + 1e-300)

    # return data
    return bins, bin_centers, hist_noisy, g, (q_alpha, E_J)

def fit_Z_from_histogram(bin_centers, counts, g_func, q_alpha, E_J, Z_search=np.arange(10, 120, 0.5)):
    """
    Fit Z by minimizing residuals where model = S * (Z^2 * g(theta)).
    For fixed Z, optimal S = sum(counts * gZ) / sum(gZ^2).
    Choose Z that minimizes weighted residual sum-of-squares.
    """
    best = None
    rss_list = []
    for Z in Z_search:
        model_basis = (Z**2) * g_func  # vector
        S_opt = np.sum(counts * model_basis) / (np.sum(model_basis**2) + 1e-300)
        model = S_opt * model_basis
        rss = np.sum((counts - model)**2)
        rss_list.append(rss)
        if best is None or rss < best[0]:
            best = (rss, Z, S_opt, model)
    rss_best, Z_best, S_best, model_best = best
    return Z_best, S_best, model_best, np.array(rss_list)

def run_rutherford_demo():
    print("\n--- Rutherford Monte Carlo + Fit Demo ---")
    print("This simulates alpha particles scattering off a nucleus and attempts to recover Z.")
    try:
        Z = int(input("Enter true nucleus Z (e.g. 79): ") or "79")
    except:
        Z = 79
    E_keV = float(input("Alpha energy (keV, e.g. 5000): ") or "5000")
    N = int(input("Number of particles (e.g. 15000): ") or "15000")
    b_max = float(input("Max impact parameter (fm, e.g. 1000): ") or "1000")
    noise = float(input("Histogram noise fraction (e.g. 0.05): ") or "0.05")

    bins, centers, hist_noisy, g, (q_alpha, E_J) = simulate_rutherford(Z=Z, E_keV=E_keV, N=N, b_max_fm=b_max, noise_frac=noise)
    # Fit Z over a search range
    Z_grid = np.arange(max(1, Z-40), Z+40, 0.25)
    Z_fit, S_fit, model_fit, rss = fit_Z_from_histogram(centers, hist_noisy, g, q_alpha, E_J, Z_search=Z_grid)

    print(f"\nTrue Z = {Z}")
    print(f"Fitted Z (grid search) = {Z_fit:.2f}")
    print("Interpretation: Rutherford's analysis could extract nucleus charge from angular scattering patterns.\n")

    # Plots
    plt.figure(figsize=(10,4))
    plt.subplot(1,2,1)
    plt.title("Angle histogram (simulated data)")
    plt.bar(centers, hist_noisy, width=bins[1]-bins[0], alpha=0.6, label="Simulated (noisy)")
    plt.plot(centers, model_fit, 'r-', lw=2, label=f"Best-fit model (Z={Z_fit:.2f})")
    plt.xlabel("Scattering angle (deg)")
    plt.ylabel("Counts")
    plt.legend()
    plt.grid(True)

    plt.subplot(1,2,2)
    plt.title("Residuals and RSS vs Z (grid search)")
    plt.plot(Z_grid, rss / np.max(rss), '-k')
    plt.axvline(Z_fit, color='r', linestyle='--', label=f"Z_fit={Z_fit:.2f}")
    plt.xlabel("Z (tested)")
    plt.ylabel("Normalized RSS")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

# ----------------------------
# 2) Animate 1s -> n hydrogenic orbital probability density
# ----------------------------
def psi_1s_prob_2d(X, Y, Z_nuclear=1):
    # hydrogenic 1s probability density (normalized): pref * exp(-2r/a)
    a_eff = a0 / Z_nuclear
    r = np.sqrt(X**2 + Y**2)
    pref = 1.0 / (math.pi * a_eff**3)
    return pref * np.exp(-2.0 * r / a_eff)

def psi_n_radial_probability_radial(r, n, Z_nuclear=1):
    # Not used explicitly; we'll show 1s and scale for higher n roughly
    a_eff = a0 / Z_nuclear
    # approximate radial scaling: peak radius ~ n^2 a_eff
    return np.exp(-2.0 * r / (n**2 * a_eff))

def animate_orbitals():
    print("\n--- Animate hydrogenic probability density changes (1s -> n) ---")
    try:
        Z_nuc = int(input("Nuclear charge Z (1=H, 2=He+,...): ") or "1")
    except:
        Z_nuc = 1
    max_n = int(input("Max principal quantum number n to show (e.g. 4): ") or "4")
    extent_A = float(input("Half-width to show (Å) (e.g. 5): ") or "5")
    res = int(input("Grid resolution (e.g. 200): ") or "200")

    extent = extent_A * 1e-10
    x = np.linspace(-extent, extent, res)
    y = np.linspace(-extent, extent, res)
    X, Y = np.meshgrid(x, y, indexing='xy')

    fig, ax = plt.subplots(figsize=(6,5))
    ax.set_title("Hydrogenic probability density (slice z=0)")
    im = ax.imshow(np.zeros_like(X).T, origin='lower', extent=[-extent_A, extent_A, -extent_A, extent_A],
                   cmap='inferno', vmin=0, vmax=1)
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label('Normalized |ψ|^2')

    ns = list(range(1, max_n+1))
    def frame(i):
        n = ns[i % len(ns)]
        # For simplicity, show 1s shape scaled radially for n: we approximate by scaling a
        a_eff = a0 / Z_nuc
        r = np.sqrt(X**2 + Y**2)
        # Use approximate hydrogenic radial probability scaling for n levels (qualitative)
        prob = (1.0/(math.pi * (n**2 * a_eff)**3)) * np.exp(-2.0 * r/(n**2 * a_eff))
        prob_norm = prob / (prob.max() + 1e-300)
        im.set_data(prob_norm.T)
        ax.set_title(f"Hydrogenic |ψ|² slice z=0 — n={n}, Z={Z_nuc}")
        return (im,)

    anim = FuncAnimation(fig, frame, frames=len(ns), interval=900, blit=False, repeat=True)
    plt.show()
    print("Animation finished: each frame shows approximate 2D probability density for n increasing.")

# ----------------------------
# 3) Stern–Gerlach spin-1/2 beam demo
# ----------------------------
def run_stern_gerlach_demo():
    print("\n--- Stern–Gerlach Spin-1/2 Beam Demo ---")
    print("We simulate many spin-1/2 particles with random initial spin orientation (Bloch vectors),")
    print("they pass through a magnetic gradient that splits them by spin projection along z into two spots.")
    N = int(input("Number of atoms/particles (e.g. 2000): ") or "2000")
    spread = float(input("Initial transverse spread (mm, e.g. 1.0): ") or "1.0")
    grad_strength = float(input("Gradient strength factor (arb. units, e.g. 1.0): ") or "1.0")
    random_seed = int(input("Random seed (e.g. 1): ") or "1")

    rng = np.random.default_rng(random_seed)
    # initial positions (x) around beam center; we show vertical splitting (y-axis)
    x0 = rng.normal(0.0, spread*1e-3, N)
    # initial spin directions: sample uniformly on Bloch sphere
    u = rng.random(N)
    v = rng.random(N)
    theta = np.arccos(1 - 2*u)  # polar
    phi = 2*np.pi*v
    sz = np.cos(theta)  # projection along z = cos(theta)
    # Particles deflect proportional to sz * grad_strength; map to y displacement
    deflection = sz * grad_strength
    # Add some thermal/turbulent perpendicular spread
    noise = rng.normal(0.0, 0.5*1e-3, N)
    y_final = deflection + noise

    # Plot initial & final distributions and animated splitting
    plt.figure(figsize=(8,4))
    plt.subplot(1,2,1)
    plt.title("Initial beam (x vs y~0)")
    plt.scatter(x0*1e3, np.zeros_like(x0), s=4, alpha=0.6)
    plt.xlabel("x (mm)")

    plt.subplot(1,2,2)
    plt.title("After Stern-Gerlach (split by spin z)")
    plt.scatter(x0*1e3, y_final*1e3, s=4, alpha=0.6)
    plt.xlabel("x (mm)")
    plt.ylabel("y deflection (mm)")
    plt.tight_layout()
    plt.show()

    # Histogram of deflections -> two peaks around ±grad_strength
    plt.figure(figsize=(6,4))
    plt.hist(y_final*1e3, bins=60, color='C2', alpha=0.8)
    plt.title("Histogram of final vertical deflections (mm)")
    plt.xlabel("deflection (mm)")
    plt.grid(True)
    plt.show()

    print("\nExplanation:")
    print(" - Spin-1/2 particles have projections +ħ/2 or -ħ/2 along any quantization axis.")
    print(" - The magnetic-field gradient gives differing forces for those projections and splits the beam.")
    print(" - The simulation sampled random spin directions; projection onto z determines deflection sign.")
    print(" - Observed two peaks correspond to the two spin projection outcomes (quantization).")

# ----------------------------
# 4) Spin precession (Bloch sphere)
# ----------------------------
def run_bloch_precession():
    print("\n--- Spin Precession (Bloch Sphere) ---")
    print("Simulates Larmor precession of a magnetic moment around an applied B field.")
    try:
        B_mag = float(input("Magnetic field magnitude (arb units, e.g. 1.0): ") or "1.0")
        omega = float(input("Larmor angular frequency ω (rad/s) (e.g. 2.0): ") or "2.0")
        t_max = float(input("Total time (s) to animate (e.g. 6.28): ") or "6.28")
    except:
        B_mag, omega, t_max = 1.0, 2.0, 6.28

    # initial Bloch vector (unit)
    theta0 = np.radians(30.0)
    phi0 = 0.0
    S0 = np.array([np.sin(theta0)*np.cos(phi0), np.sin(theta0)*np.sin(phi0), np.cos(theta0)])

    times = np.linspace(0, t_max, 200)
    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111, projection='3d')
    # Draw unit sphere (wireframe)
    u = np.linspace(0, 2*np.pi, 40)
    v = np.linspace(0, np.pi, 20)
    xs = np.outer(np.cos(u), np.sin(v))
    ys = np.outer(np.sin(u), np.sin(v))
    zs = np.outer(np.ones_like(u), np.cos(v))
    ax.plot_wireframe(xs, ys, zs, color='lightgray', alpha=0.6)

    vec_line, = ax.plot([], [], [], 'o-', lw=2, color='C1')
    path_line, = ax.plot([], [], [], '-', lw=1, color='C0', alpha=0.6)

    ax.set_xlim([-1,1]); ax.set_ylim([-1,1]); ax.set_zlim([-1,1])
    ax.set_xlabel('Sx'); ax.set_ylabel('Sy'); ax.set_zlabel('Sz')
    ax.set_title('Bloch sphere: spin precession')

    # precompute trajectory: rotate around z-axis for simplicity (B along z)
    trajectory = np.zeros((len(times), 3))
    for i, t in enumerate(times):
        phi = omega * t
        # rotation around z by phi applied to initial transverse component
        Sx = S0[0]*np.cos(phi) - S0[1]*np.sin(phi)
        Sy = S0[0]*np.sin(phi) + S0[1]*np.cos(phi)
        Sz = S0[2]
        trajectory[i] = [Sx, Sy, Sz]

    def init():
        vec_line.set_data([], [])
        vec_line.set_3d_properties([])
        path_line.set_data([], [])
        path_line.set_3d_properties([])
        return vec_line, path_line

    def update(i):
        v = trajectory[i]
        vec_line.set_data([0, v[0]], [0, v[1]])
        vec_line.set_3d_properties([0, v[2]])
        path_line.set_data(trajectory[:i,0], trajectory[:i,1])
        path_line.set_3d_properties(trajectory[:i,2])
        return vec_line, path_line

    anim = FuncAnimation(fig, update, frames=len(times), init_func=init, interval=50, blit=False)
    plt.show()
    print("Larmor precession shown: spin vector precesses around B direction (here chosen as z-axis).")

# ----------------------------
# Main menu
# ----------------------------
def main_menu():
    menu = """
Atomic extensions menu:
 1) Rutherford Monte Carlo + fit for Z
 2) Animate hydrogenic orbitals (1s -> n)
 3) Stern–Gerlach spin-1/2 demo
 4) Spin precession (Bloch sphere)
 0) Quit
Choose option: """
    while True:
        choice = input(menu).strip()
        if choice == "1":
            run_rutherford_demo()
        elif choice == "2":
            animate_orbitals()
        elif choice == "3":
            run_stern_gerlach_demo()
        elif choice == "4":
            run_bloch_precession()
        elif choice == "0" or choice.lower() in ("q", "quit", "exit"):
            print("Goodbye — enjoy exploring atomic physics!")
            break
        else:
            print("Invalid choice. Enter 1..4 or 0 to quit.")

if __name__ == "__main__":
    main_menu()

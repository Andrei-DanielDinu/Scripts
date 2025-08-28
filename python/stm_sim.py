import numpy as np
import matplotlib.pyplot as plt

# -----------------------------
# STM CONSTANTS & UTILITIES
# -----------------------------
e = 1.602176634e-19      # Elementary charge [C]
hbar = 1.054_571_817e-34 # Reduced Planck constant [J·s]
m_e = 9.109_383_7015e-31 # Electron mass [kg]

def kappa(work_function_eV):
    """
    Tunneling decay constant κ = sqrt(2 m Φ) / ħ
    Φ given in eV → convert to Joules: Φ[J] = Φ[eV] * e
    """
    Phi_J = work_function_eV * e
    return np.sqrt(2 * m_e * Phi_J) / hbar  # [1/m]

def tunneling_current(V, z_gap, kappa_value, I0=1e-6):
    """
    Simple STM current model: I = I0 * V * exp(-2 κ z)
    - V: bias [V]
    - z_gap: tip-sample separation [m]
    - κ: decay constant [1/m]
    - I0: scale factor to keep numbers nice (arbitrary)
    Returns current [A]
    """
    return I0 * V * np.exp(-2 * kappa_value * z_gap)

# -----------------------------
# SURFACE MODELS
# -----------------------------
def make_surface(nx, ny, dx_nm, dy_nm, mode="lattice+step+defects",
                 amp_pm=50, lattice_nm=0.35, step_height_pm=200, defect_density=0.01):
    """
    Generate a synthetic sample surface height h(x,y).
    Units:
      - dx_nm, dy_nm: pixel size [nm] (display convenience)
      - Output height in meters (internally use SI)
    'amp_pm' = sinusoidal corrugation amplitude in picometers
    """
    x = (np.arange(nx) - nx/2) * dx_nm  # nm
    y = (np.arange(ny) - ny/2) * dy_nm  # nm
    X, Y = np.meshgrid(x, y, indexing='ij')

    h_pm = np.zeros((nx, ny))  # picometers

    if "lattice" in mode:
        a = lattice_nm
        h_pm += amp_pm * (np.sin(2*np.pi*X/a) + np.sin(2*np.pi*Y/a)) * 0.5

    if "step" in mode:
        # Single mono-atomic step across y=0
        h_pm += (Y > 0).astype(float) * step_height_pm

    if "defect" in mode or "defects" in mode:
        rng = np.random.default_rng(42)
        num_def = int(defect_density * nx * ny)
        for _ in range(num_def):
            i = rng.integers(0, nx)
            j = rng.integers(0, ny)
            radius = rng.integers(2, 6)
            bump = rng.choice([+1, -1]) * rng.uniform(20, 80)  # pm
            ii = slice(max(i-radius, 0), min(i+radius+1, nx))
            jj = slice(max(j-radius, 0), min(j+radius+1, ny))
            # Smooth gaussian-like bump
            ii_idx = np.arange(ii.start, ii.stop)[:, None]
            jj_idx = np.arange(jj.start, jj.stop)[None, :]
            rr2 = (ii_idx - i)**2 + (jj_idx - j)**2
            h_pm[ii, jj] += bump * np.exp(-rr2 / (2*(radius/2.0)**2))

    # Convert pm → meters
    h_m = h_pm * 1e-12
    return h_m, X, Y

# -----------------------------
# STM SCAN (CONSTANT-CURRENT)
# -----------------------------
def stm_scan_constant_current(h_surface_m,
                              V_bias=0.1,
                              work_function_eV=4.5,
                              I_set_nA=1.0,
                              z_init_A=6.0,
                              noise_frac=0.02,
                              lateral_kappa_variation=False):
    """
    Simulate a constant-current STM topography:
      - Given true surface height h(x,y) [m],
      - Find tip height z_tip(x,y) [m] so that I ≈ I_set at each (x,y).

    We mimic a simple feedback:
      I = I0 * V * exp(-2 κ (z_tip - h))
      Solve for z_tip:
        z_tip = h + (1/(2κ)) * ln(I0*V / I_set) + noise
    - Optionally add small lateral variations in κ to mimic local work function changes.
    """
    kappa0 = kappa(work_function_eV)
    I_set = I_set_nA * 1e-9
    I0 = 1e-6  # scale in current model

    nx, ny = h_surface_m.shape
    z_tip = np.zeros_like(h_surface_m)

    # Optional tiny spatial variation in barrier height/work function
    if lateral_kappa_variation:
        rng = np.random.default_rng(123)
        dPhi = 0.10  # ±0.10 eV local variation
        Phi_map = work_function_eV + dPhi * (rng.random((nx, ny)) - 0.5)
        kappa_map = kappa(Phi_map)
    else:
        kappa_map = np.full((nx, ny), kappa0)

    # Exact inversion for constant-current setpoint (plus small measurement noise)
    rng_noise = np.random.default_rng(7)
    meas_noise = rng_noise.normal(0.0, noise_frac, size=(nx, ny))  # fractional noise

    # Avoid log of zero/negative
    eps = 1e-24
    target_ratio = (I0 * V_bias) / np.maximum(I_set, eps)
    log_term = np.log(np.maximum(target_ratio, eps))

    # z_tip (m)
    z_tip = h_surface_m + (log_term / (2.0 * kappa_map))

    # Add small noise as random offset in measured current → equivalent height noise:
    # If I_meas = I_set * (1 + δ), then Δz ≈ (1/(2κ)) ln(1/(1+δ)) ≈ -(δ/(2κ)) for small δ.
    z_tip_noisy = z_tip - (meas_noise / (2.0 * kappa_map))

    # For reference: initial gap if surface were flat
    z_init_m = z_init_A * 1e-10

    info = {
        "kappa_mean": float(np.mean(kappa_map)),
        "kappa_std": float(np.std(kappa_map)),
        "z_offset_mean_A": float(np.mean((z_tip - h_surface_m)) * 1e10),
        "z_offset_std_A": float(np.std((z_tip - h_surface_m)) * 1e10),
        "z_init_A": z_init_A,
        "I_set_nA": I_set_nA,
        "V_bias_V": V_bias,
        "Phi_eV": work_function_eV,
        "noise_frac": noise_frac,
        "lateral_kappa_variation": lateral_kappa_variation,
    }

    return z_tip_noisy, info

# -----------------------------
# MAIN (USER INPUT + PLOTS)
# -----------------------------
def main():
    print("🧪 Scanning Tunneling Microscope (STM) — Constant-Current Simulation")
    print("\nConcept:")
    print(" - A sharp metallic tip is brought nanometers above a conducting surface.")
    print(" - Electrons quantum-tunnel across the gap; current I ∝ V · exp(-2 κ z).")
    print(" - A feedback loop moves the tip up/down to keep I = I_set constant.")
    print(" - The recorded tip height z_tip(x,y) is the STM 'topography'.\n")

    # ----- User inputs -----
    try:
        nx = int(input("Grid points along X (e.g., 128): ").strip() or "128")
        ny = int(input("Grid points along Y (e.g., 128): ").strip() or "128")
        dx_nm = float(input("Pixel size X [nm] (e.g., 0.05): ").strip() or "0.05")
        dy_nm = float(input("Pixel size Y [nm] (e.g., 0.05): ").strip() or "0.05")

        surface_mode = input("Surface type [lattice / lattice+step / lattice+step+defects] (default lattice+step+defects): ").strip() or "lattice+step+defects"
        amp_pm = float(input("Lattice corrugation amplitude [pm] (e.g., 40): ").strip() or "40")
        step_height_pm = float(input("Step height [pm] (e.g., 200): ").strip() or "200")
        defect_density = float(input("Defect density (0–0.05 recommended, e.g., 0.01): ").strip() or "0.01")
        lattice_nm = float(input("Lattice period [nm] (e.g., 0.35): ").strip() or "0.35")

        V_bias = float(input("Bias voltage V [V] (e.g., 0.1): ").strip() or "0.1")
        Phi_eV = float(input("Work function Φ [eV] (e.g., 4.5): ").strip() or "4.5")
        I_set_nA = float(input("Current setpoint [nA] (e.g., 1.0): ").strip() or "1.0")
        z_init_A = float(input("Initial gap guess [Å] (for info only, e.g., 6.0): ").strip() or "6.0")
        noise_frac = float(input("Measurement noise (fraction, e.g., 0.02 = 2%): ").strip() or "0.02")
        vary_kappa = input("Simulate lateral work-function variation? [y/N]: ").strip().lower() == "y"
    except Exception:
        print("Using defaults due to input error.")
        nx, ny = 128, 128
        dx_nm = dy_nm = 0.05
        surface_mode = "lattice+step+defects"
        amp_pm, step_height_pm = 40, 200
        defect_density, lattice_nm = 0.01, 0.35
        V_bias, Phi_eV, I_set_nA, z_init_A, noise_frac = 0.1, 4.5, 1.0, 6.0, 0.02
        vary_kappa = False

    # ----- Build surface -----
    mode_flag = surface_mode.lower()
    use_defects = "defect" in mode_flag
    use_step = "step" in mode_flag
    base_mode = "lattice"
    if use_step and use_defects:
        mode_used = "lattice+step+defects"
    elif use_step:
        mode_used = "lattice+step"
    else:
        mode_used = "lattice"

    h_surface_m, X_nm, Y_nm = make_surface(
        nx, ny, dx_nm, dy_nm,
        mode=mode_used,
        amp_pm=amp_pm,
        lattice_nm=lattice_nm,
        step_height_pm=step_height_pm,
        defect_density=defect_density
    )

    # ----- Simulate STM constant-current topography -----
    z_topo_m, info = stm_scan_constant_current(
        h_surface_m,
        V_bias=V_bias,
        work_function_eV=Phi_eV,
        I_set_nA=I_set_nA,
        z_init_A=z_init_A,
        noise_frac=noise_frac,
        lateral_kappa_variation=vary_kappa
    )

    # ----- Explanations / numbers -----
    print("\n📘 Explanations & Key Numbers")
    print(f" - Work function Φ = {Phi_eV:.2f} eV → κ = {info['kappa_mean']:.3e} 1/m (mean)")
    if info["lateral_kappa_variation"]:
        print(f"   (with lateral variation; κ std = {info['kappa_std']:.3e} 1/m)")
    print(f" - Current setpoint I_set = {info['I_set_nA']:.3f} nA at bias V = {info['V_bias_V']:.3f} V")
    print(f" - Constant-current offset <z_tip - h> ≈ {info['z_offset_mean_A']:.2f} Å (std {info['z_offset_std_A']:.2f} Å)")
    print(" - In constant-current mode, the apparent height includes true topography")
    print("   + electronic effects (work-function/barrier changes) + noise.")
    print(" - Because I ~ exp(-2κz), even picometer height changes cause sizable current changes.")

    # ----- Plots -----
    extent_nm = [X_nm.min(), X_nm.max(), Y_nm.min(), Y_nm.max()]

    fig, axs = plt.subplots(1, 3, figsize=(15, 4.8), constrained_layout=True)

    im0 = axs[0].imshow(h_surface_m.T * 1e12, origin='lower', extent=extent_nm, cmap='viridis', aspect='equal')
    axs[0].set_title("True Surface Height h(x,y) [pm]")
    axs[0].set_xlabel("x [nm]"); axs[0].set_ylabel("y [nm]")
    fig.colorbar(im0, ax=axs[0], fraction=0.046, pad=0.04)

    im1 = axs[1].imshow(z_topo_m.T * 1e12, origin='lower', extent=extent_nm, cmap='magma', aspect='equal')
    axs[1].set_title("STM Topography z_tip(x,y) [pm]\n(Constant-current)")
    axs[1].set_xlabel("x [nm]"); axs[1].set_ylabel("y [nm]")
    fig.colorbar(im1, ax=axs[1], fraction=0.046, pad=0.04)

    # Line profile through the center
    mid = ny // 2
    x_axis_nm = X_nm[:, 0]
    axs[2].plot(x_axis_nm, (h_surface_m[:, mid] * 1e12), label="True surface h", lw=2)
    axs[2].plot(x_axis_nm, (z_topo_m[:, mid] * 1e12), label="STM topography z_tip", lw=2, alpha=0.8)
    axs[2].set_title("Center Line Profile [pm]")
    axs[2].set_xlabel("x [nm]"); axs[2].set_ylabel("height [pm]")
    axs[2].legend()
    axs[2].grid(True, alpha=0.3)

    plt.suptitle("STM Simulation — Constant-Current Imaging", y=1.03, fontsize=14)
    plt.show()

if __name__ == "__main__":
    main()

# coulomb_k_fit.py
# Numerically estimate Coulomb's constant k from force/charge/distance data.
# Includes: manual data entry OR synthetic data generator, weighted least squares,
# inverse-square (exponent) check via log fit, and explanatory plots.

import numpy as np
import matplotlib.pyplot as plt

K_TRUE = 8.9875517923e9  # CODATA 2018 Coulomb constant (N·m^2/C^2), for reference only

def wls_through_origin(X, Y, sigma_Y=None):
    """
    Weighted least-squares fit of Y = k * X, forced through origin.
    Returns k_hat, sigma_k, residuals, R^2.
    If sigma_Y is None, uses equal weights.
    """
    X = np.asarray(X, dtype=float)
    Y = np.asarray(Y, dtype=float)
    n = len(X)

    if sigma_Y is None:
        w = np.ones_like(Y)
    else:
        w = 1.0 / np.asarray(sigma_Y, dtype=float)**2

    Sxx = np.sum(w * X * X)
    Sxy = np.sum(w * X * Y)

    k_hat = Sxy / Sxx

    # Residuals and variance estimate
    Y_hat = k_hat * X
    residuals = Y - Y_hat
    # Weighted residual sum of squares
    WRSS = np.sum(w * residuals**2)
    # Degrees of freedom: n - 1 (one parameter)
    dof = max(n - 1, 1)
    sigma2 = WRSS / dof
    # Variance of slope (through origin, weighted)
    var_k = sigma2 / Sxx
    sigma_k = np.sqrt(var_k)

    # R^2 using unweighted total variance around mean
    Y_mean = np.mean(Y)
    SS_tot = np.sum((Y - Y_mean)**2)
    SS_res = np.sum((Y - Y_hat)**2)
    R2 = 1.0 - SS_res / SS_tot if SS_tot > 0 else np.nan

    return k_hat, sigma_k, residuals, R2

def log_fit_inverse_square(q1, q2, r, F, sigma_F=None):
    """
    Estimate the inverse power n and k from:
      F = k * |q1 q2| / r^n
    Take logs and rearrange:
      ln F - ln|q1 q2| = ln k - n ln r
    Do an ordinary least squares fit of Y' vs X = ln r:
      Y' = a + b * X, where a = ln k and b = -n
    Returns: n_hat, k_hat_log, slope_std, intercept_std
    """
    qprod_abs = np.abs(q1 * q2)
    mask = (qprod_abs > 0) & (r > 0) & (F > 0)
    if not np.all(mask):
        qprod_abs = qprod_abs[mask]
        r = r[mask]
        F = F[mask]
        if sigma_F is not None:
            sigma_F = sigma_F[mask]

    X = np.log(r)
    Yp = np.log(F) - np.log(qprod_abs)

    # Simple OLS (unweighted) linear regression
    xbar = np.mean(X)
    ybar = np.mean(Yp)
    Sxx = np.sum((X - xbar)**2)
    Sxy = np.sum((X - xbar) * (Yp - ybar))
    slope = Sxy / Sxx
    intercept = ybar - slope * xbar

    # Residuals and standard errors
    Yp_hat = intercept + slope * X
    res = Yp - Yp_hat
    dof = max(len(X) - 2, 1)
    s2 = np.sum(res**2) / dof
    slope_std = np.sqrt(s2 / Sxx)
    intercept_std = np.sqrt(s2 * (1.0/len(X) + xbar**2 / Sxx))

    n_hat = -slope
    k_hat_log = np.exp(intercept)

    return n_hat, k_hat_log, slope_std, intercept_std, (X, Yp, Yp_hat)

def get_user_data():
    print("\n=== Coulomb’s Constant Experiment (Numerical Fit) ===")
    print("Choose data mode:")
    print("  1) Manual entry (you type charge1, charge2, distance, measured force)")
    print("  2) Synthetic data (script generates data with realistic noise)")
    mode = input("Enter 1 or 2: ").strip()

    if mode == "1":
        m = int(input("\nHow many measurements will you enter? "))
        q1, q2, r, F, sF = [], [], [], [], []
        print("\nEnter each measurement:")
        print("  q1, q2 in Coulombs; r in meters; F in Newtons; optional σ_F (press Enter to skip).")
        for i in range(m):
            q1_i = float(input(f"[{i+1}] q1 (C): "))
            q2_i = float(input(f"[{i+1}] q2 (C): "))
            r_i  = float(input(f"[{i+1}] r  (m): "))
            F_i  = float(input(f"[{i+1}] F_measured (N): "))
            sF_in = input(f"[{i+1}] σ_F (N, optional): ").strip()
            q1.append(q1_i); q2.append(q2_i); r.append(r_i); F.append(F_i)
            sF.append(float(sF_in) if sF_in else np.nan)

        q1 = np.array(q1); q2 = np.array(q2); r = np.array(r); F = np.array(F)
        sigma_F = np.array(sF)
        if np.all(np.isnan(sigma_F)):
            sigma_F = None  # use equal weights
        return q1, q2, r, F, sigma_F, "manual"

    # Synthetic data
    print("\nSynthetic data generator:")
    print("You can fix charges and vary the distance, or vary charges as well.")
    n = int(input("Number of measurements to generate (e.g., 12): ") or "12")
    q1_val = float(input("Baseline q1 (C) (e.g., 1e-6): ") or "1e-6")
    q2_val = float(input("Baseline q2 (C) (e.g., 1e-6): ") or "1e-6")
    r_min = float(input("Min distance r_min (m) (e.g., 0.02): ") or "0.02")
    r_max = float(input("Max distance r_max (m) (e.g., 0.20): ") or "0.20")
    vary_q = input("Vary charges a bit? (y/n) [y]: ").strip().lower() or "y"
    noise_pct = float(input("Multiplicative noise on F (percent, e.g., 5): ") or "5")
    add_bias = float(input("Add constant force bias (N, e.g., 0): ") or "0")

    rng = np.random.default_rng(42)
    r = np.linspace(r_min, r_max, n)
    if vary_q == "y":
        q1 = q1_val * (1 + rng.normal(0, 0.1, size=n))
        q2 = q2_val * (1 + rng.normal(0, 0.1, size=n))
    else:
        q1 = np.full(n, q1_val)
        q2 = np.full(n, q2_val)

    F_clean = K_TRUE * q1 * q2 / (r**2)
    F = F_clean * (1 + rng.normal(0, noise_pct/100.0, size=n)) + add_bias
    # Optional measurement uncertainty: assume known fractional uncertainty
    sigma_F = np.abs(F_clean) * (noise_pct/100.0)

    return q1, q2, r, F, sigma_F, "synthetic"

def main():
    q1, q2, r, F, sigma_F, mode = get_user_data()

    # Build the "predictor" X = q1*q2 / r^2 (so model is F = k * X)
    X = q1 * q2 / (r**2)

    print("\n--- Running fits ---")
    print("1) Weighted fit of F = k * (q1*q2/r^2) (through origin).")
    k_hat, sigma_k, residuals, R2 = wls_through_origin(X, F, sigma_Y=sigma_F)

    print("2) Inverse-square check: fit exponent n and k via log form.")
    n_hat, k_hat_log, slope_std, intercept_std, logpack = log_fit_inverse_square(q1, q2, r, F, sigma_F)

    # Compare to CODATA for context (not used in fitting)
    err1_pct = (k_hat - K_TRUE) / K_TRUE * 100.0
    err2_pct = (k_hat_log - K_TRUE) / K_TRUE * 100.0

    # ---------- Explanations ----------
    print("\n================ RESULTS & EXPLANATIONS ================\n")

    print("A) Direct Fit (through origin):  F = k · (q1 q2 / r^2)")
    print(f"   Estimated k = {k_hat:.6e} N·m^2/C^2")
    print(f"   1σ uncertainty (fit) ≈ {sigma_k:.3e}  →  ~{(sigma_k/k_hat*100):.2f}% relative")
    print(f"   Coefficient of determination R^2 = {R2:.4f} (closer to 1 means better fit)")
    print(f"   Context (CODATA reference): k_true = {K_TRUE:.6e} N·m^2/C^2  →  error ≈ {err1_pct:+.2f}%\n")
    print("   Meaning: We regress measured forces against (q1 q2 / r^2).")
    print("            The slope of that line is an estimate of Coulomb’s constant k.")
    print("            We force the line through the origin because F should be 0 when X=0.\n")

    print("B) Inverse-Square Law Check (log fit):  F = k · |q1 q2| / r^n")
    print(f"   Estimated exponent n = {n_hat:.4f}  (ideal inverse-square → n = 2)")
    print(f"   Estimated k (from log fit) = {k_hat_log:.6e} N·m^2/C^2  →  error vs CODATA ≈ {err2_pct:+.2f}%\n")
    print("   Meaning: We subtract ln|q1 q2| from ln F and fit a straight line vs ln r.")
    print("            The slope gives −n and the intercept gives ln k.")
    print("            If n is close to 2, your data supports the inverse-square law.\n")

    if mode == "synthetic":
        print("Note: You generated synthetic data around k_true with noise;")
        print("      small deviations in k and n are expected due to the added noise/bias.\n")
    else:
        print("Note: For manual data, large deviations might indicate:")
        print("      • mis-measured distances/charges,")
        print("      • systematic bias (e.g., sensor offset),")
        print("      • or that forces were too small compared to noise.\n")

    # ---------- Plots ----------
    # 1) F vs X with fitted line
    plt.figure(figsize=(10, 4))
    plt.subplot(1, 2, 1)
    plt.title("Force vs (q1 q2 / r²) with Fitted k")
    if sigma_F is not None:
        plt.errorbar(X, F, yerr=sigma_F, fmt="o", ms=4, label="Data", alpha=0.8)
    else:
        plt.plot(X, F, "o", ms=4, label="Data", alpha=0.8)
    X_line = np.linspace(min(X), max(X), 200)
    plt.plot(X_line, k_hat * X_line, "-", label=f"Fit: F = k·X\nk ≈ {k_hat:.3e}")
    plt.xlabel("X = q1·q2 / r²  [C²/m²]")
    plt.ylabel("Force F  [N]")
    plt.legend()
    plt.grid(True)

    # Residuals
    plt.subplot(1, 2, 2)
    plt.title("Residuals (F_meas − F_fit)")
    F_fit = k_hat * X
    res = F - F_fit
    plt.axhline(0, lw=1, color="k")
    plt.plot(X, res, "o", ms=4, alpha=0.8)
    plt.xlabel("X = q1·q2 / r²")
    plt.ylabel("Residual [N]")
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    # 2) Log–log inverse-square check
    Xlog, Yp, Yp_hat = logpack  # Xlog = ln r; Yp = ln F − ln|q1 q2|
    plt.figure(figsize=(8, 4))
    plt.title("Inverse-Square Check: ln F − ln|q1 q2| vs ln r")
    plt.plot(Xlog, Yp, "o", ms=4, label="Transformed data", alpha=0.8)
    xx = np.linspace(min(Xlog), max(Xlog), 200)
    # Straight line with slope = -(n_hat) and intercept = ln(k_hat_log)
    plt.plot(xx, np.log(k_hat_log) - n_hat * xx, "-", label=f"Fit: slope = −{n_hat:.3f}")
    plt.xlabel("ln r")
    plt.ylabel("ln F − ln|q1 q2|")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()

import numpy as np

def electrolysis_experiment():
    print("⚡ Faraday’s Electrolysis and Avogadro’s Number ⚡\n")
    print("Principle:")
    print(" - Passing current through an electrolyte deposits/consumes matter at the electrodes.")
    print(" - The deposited mass depends on the total charge Q = I * t.")
    print(" - From Faraday’s constant (F), we can deduce Avogadro’s number: N_A = F / e.\n")

    # Known constants
    e = 1.602e-19      # C (elementary charge)
    M = 63.546e-3      # kg/mol (Molar mass of Copper, example)
    z = 2              # Cu²⁺ needs 2 electrons
    I = 0.5            # A (current)
    t = 1800           # s (time = 30 min)

    # Experimentally measured deposited mass (simulation)
    Q = I * t                      # Total charge
    F_exp = (Q * M) / (z * 1e-3)   # From formula rearranged (m=1g assumption)

    # If we assume 1 g of copper deposited:
    m = 1e-3   # kg
    F_measured = (Q * M) / (z * m)

    # Calculate Avogadro number
    N_A_measured = F_measured / e

    print("Experimental Setup:")
    print(f" - Electrolyte: Copper(II) sulfate (CuSO₄)")
    print(f" - Current I = {I} A, Time t = {t} s → Q = {Q:.1f} C")
    print(f" - Molar mass M = {M*1000:.3f} g/mol, z = {z}")
    print(f" - Deposited mass (assumed measurement) = {m*1000:.3f} g\n")

    print("Results:")
    print(f" - Measured Faraday constant F = {F_measured:.3e} C/mol")
    print(f" - Avogadro’s number N_A = {N_A_measured:.3e} mol⁻¹")
    print("\nKnown values:")
    print(f" - True F ≈ 96485 C/mol")
    print(f" - True N_A ≈ 6.022×10^23 mol⁻¹")

if __name__ == "__main__":
    electrolysis_experiment()

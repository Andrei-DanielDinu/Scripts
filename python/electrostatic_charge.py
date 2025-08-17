import math
import matplotlib.pyplot as plt
import random

# --- CONSTANTS ---
k_e = 8.9875517923e9  # Coulomb's constant (N·m²/C²)
e_charge = 1.602e-19  # elementary charge (Coulombs)

# --- USER INPUT ---
material1 = input("Enter first material (e.g. wool, glass, rubber): ").lower()
material2 = input("Enter second material (e.g. plastic, silk, amber): ").lower()
friction_time = float(input("Friction duration (seconds): "))
distance = float(input("Distance between objects (meters): "))

# --- TRIBOELECTRIC COEFFICIENTS (very simplified model) ---
tribo_strength = {
    ("wool", "plastic"): 5e12,
    ("glass", "silk"): 4e12,
    ("rubber", "fur"): 3e12,
    ("amber", "wool"): 2e12
}

def estimate_charge(m1, m2, time):
    """Return charge generated (Coulombs) after friction."""
    coeff = tribo_strength.get((m1, m2), 1e12)  # default if not listed
    # Electrons transferred proportional to time
    transferred_electrons = coeff * time
    return transferred_electrons * e_charge

# --- CALCULATIONS ---
q1 = estimate_charge(material1, material2, friction_time)
q2 = -q1  # equal and opposite charges

# Add noise to simulate measurement error
noise_factor = 1 + random.uniform(-0.05, 0.05)  
q1 *= noise_factor
q2 *= noise_factor

force = k_e * abs(q1 * q2) / (distance**2)

# --- OUTPUT ---
print("\n--- Electrostatic Simulation ---")
print(f"Charge on {material1}: {q1:.3e} C")
print(f"Charge on {material2}: {q2:.3e} C")
print(f"Electrostatic force at {distance} m: {force:.3e} N")

# --- PLOT: Charge growth with friction time ---
times = [t for t in range(1, int(friction_time)+1)]
charges = [estimate_charge(material1, material2, t) for t in times]

plt.figure(figsize=(8,5))
plt.plot(times, charges, marker="o")
plt.xlabel("Friction time (s)")
plt.ylabel("Charge (Coulombs)")
plt.title("Charge buildup due to friction")
plt.grid(True)
plt.show()

import math
import matplotlib.pyplot as plt

# --- CONSTANTS ---
G_true = 6.67430e-11  # gravitational constant (m^3·kg^-1·s^-2)
pi = math.pi

# --- USER INPUT ---
m_small = float(input("Enter small sphere mass (kg): "))
m_large = float(input("Enter large sphere mass (kg): "))
distance_between_centers = float(input("Enter distance between centers (m): "))

torsion_constant = 2e-7   # torsion wire constant (N·m/rad)
damping = 1e-7            # damping factor
I = 1e-4                  # moment of inertia (kg·m^2)

# --- FORCE AND TORQUE ---
def gravitational_force(m1, m2, r):
    return G_true * m1 * m2 / (r**2)

def torque_gravity(m1, m2, r, lever_arm):
    F = gravitational_force(m1, m2, r)
    return F * lever_arm

# --- SIMULATION PARAMETERS ---
dt = 0.1   # time step (s)
T = 2000   # total time (s)

theta = 0.0   # angular displacement (rad)
omega = 0.0   # angular velocity (rad/s)

lever_arm = distance_between_centers / 2
torque_ext = torque_gravity(m_small, m_large, distance_between_centers, lever_arm)

times = []
angles = []

# --- TIME EVOLUTION ---
for step in range(int(T/dt)):
    t = step * dt
    
    # Equation of motion: I*theta'' + b*theta' + k*theta = torque
    alpha = (torque_ext - torsion_constant*theta - damping*omega) / I
    
    omega += alpha * dt
    theta += omega * dt
    
    times.append(t)
    angles.append(theta)

# --- FINAL EQUILIBRIUM METHOD ---
theta_equilibrium = angles[-1]
F_measured = torsion_constant * theta_equilibrium / lever_arm
G_est_equilibrium = F_measured * (distance_between_centers**2) / (m_small * m_large)

# --- OSCILLATION PERIOD METHOD ---
def measure_period(I, k):
    """Return natural oscillation period (s)."""
    return 2 * math.pi * math.sqrt(I / k)

T_osc = measure_period(I, torsion_constant)

# Use oscillation period in Cavendish method
# Formula (simplified model):
# G = (4 * pi^2 * lever_arm * theta) / (T^2 * m1 * m2 / r^2)
G_est_period = (4 * math.pi**2 * lever_arm * theta_equilibrium * distance_between_centers**2) / (T_osc**2 * m_small * m_large)

# --- OUTPUT ---
print("\n--- Cavendish Simulation Results ---")
print(f"Final equilibrium angle: {theta_equilibrium:.3e} rad ({math.degrees(theta_equilibrium):.6f}°)")

print("\n[Method 1: Equilibrium displacement]")
print(f"Estimated G: {G_est_equilibrium:.5e}")
print(f"Error: {abs(G_est_equilibrium - G_true)/G_true*100:.3f}%")

print("\n[Method 2: Oscillation period]")
print(f"Oscillation period T: {T_osc:.2f} s")
print(f"Estimated G: {G_est_period:.5e}")
print(f"Error: {abs(G_est_period - G_true)/G_true*100:.3f}%")

print(f"\nTrue G: {G_true:.5e}")

# --- PLOT ---
plt.figure(figsize=(10,5))
plt.plot(times, angles, label="Angular displacement (rad)")
plt.axhline(theta_equilibrium, color="red", linestyle="--", label="Equilibrium angle")
plt.xlabel("Time (s)")
plt.ylabel("Angle (rad)")
plt.title("Cavendish Experiment Simulation (with Oscillation)")
plt.legend()
plt.grid(True)
plt.show()

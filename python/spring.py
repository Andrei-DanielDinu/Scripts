import numpy as np
import matplotlib.pyplot as plt

# ------------------------
# User input
# ------------------------
m = float(input("Mass (kg): "))
k = float(input("Spring constant k (N/m): "))
c = float(input("Damping coefficient c (kg/s): "))
x0 = float(input("Initial displacement from equilibrium (m, +down): "))
v0 = float(input("Initial velocity (m/s): "))
g = 9.81  # gravity

# ------------------------
# Parameters
# ------------------------
omega0 = np.sqrt(k / m)
zeta = c / (2 * np.sqrt(m * k))
x_eq = m * g / k

if zeta < 1:
    stability = "Underdamped"
elif np.isclose(zeta, 1):
    stability = "Critically damped"
else:
    stability = "Overdamped"

print("\n--- Analysis ---")
print(f"Natural frequency: {omega0:.3f} rad/s")
print(f"Damping ratio: {zeta:.3f} → {stability}")
print(f"Equilibrium position: {x_eq:.3f} m below the ceiling\n")

# ------------------------
# Simulation
# ------------------------
t_max, dt = 10, 0.001
t = np.arange(0, t_max, dt)
x = np.zeros_like(t)
v = np.zeros_like(t)
a = np.zeros_like(t)
F_spring = np.zeros_like(t)

x[0] = x0
v[0] = v0

for i in range(1, len(t)):
    a[i-1] = (-c*v[i-1] - k*x[i-1]) / m
    v[i] = v[i-1] + a[i-1]*dt
    x[i] = x[i-1] + v[i]*dt
    F_spring[i] = -k*x[i]

a[-1] = (-c*v[-1] - k*x[-1]) / m

# ------------------------
# Metrics
# ------------------------
amplitude = np.max(np.abs(x))
max_speed = np.max(np.abs(v))
variance_x = np.var(x)
print(f"Initial amplitude: {amplitude:.4f} m")
print(f"Max speed: {max_speed:.4f} m/s")
print(f"Variance in displacement: {variance_x:.6f} m²")

# ------------------------
# Energies
# ------------------------
KE = 0.5 * m * v**2
PE_spring = 0.5 * k * x**2
E_total = KE + PE_spring

# ------------------------
# Plots
# ------------------------
fig, axs = plt.subplots(2, 2, figsize=(12, 8))
fig.suptitle("Vertical Spring-Mass System", fontsize=16)

# Displacement
axs[0, 0].plot(t, x + x_eq)
axs[0, 0].set_title("Displacement (absolute position)")
axs[0, 0].set_xlabel("Time (s)")
axs[0, 0].set_ylabel("Position (m)")
axs[0, 0].grid()

# Velocity
axs[0, 1].plot(t, v, color="orange")
axs[0, 1].set_title("Velocity")
axs[0, 1].set_xlabel("Time (s)")
axs[0, 1].set_ylabel("Velocity (m/s)")
axs[0, 1].grid()

# Force
axs[1, 0].plot(t, F_spring, color="green")
axs[1, 0].set_title("Spring Force")
axs[1, 0].set_xlabel("Time (s)")
axs[1, 0].set_ylabel("Force (N)")
axs[1, 0].grid()

# Energy
axs[1, 1].plot(t, KE, label="Kinetic")
axs[1, 1].plot(t, PE_spring, label="Potential")
axs[1, 1].plot(t, E_total, label="Total", linestyle="--")
axs[1, 1].set_title("Energy")
axs[1, 1].set_xlabel("Time (s)")
axs[1, 1].set_ylabel("Energy (J)")
axs[1, 1].legend()
axs[1, 1].grid()

plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.show()

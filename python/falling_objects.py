import numpy as np
import matplotlib.pyplot as plt

# === USER INPUTS ===
mode = input("Choose input mode ('speed' or 'force'): ").strip().lower()

mass1 = float(input("Enter mass of Object 1 (kg): "))
mass2 = float(input("Enter mass of Object 2 (kg): "))
angle_deg = float(input("Enter launch angle (degrees): "))
height = float(input("Enter initial height from ground (m): "))

# Constants
g = 9.81  # m/s²
theta = np.radians(angle_deg)  # convert to radians

if mode == "speed":
    v0 = float(input("Enter initial speed (m/s): "))
    v1 = v2 = v0
elif mode == "force":
    force = float(input("Enter initial applied force (N): "))
    # F = m * a ⇒ a = F/m ⇒ v = a * Δt ⇒ we'll assume Δt = 1s impulse for simplicity
    a1 = force / mass1
    a2 = force / mass2
    v1 = a1 * 1
    v2 = a2 * 1
    print(f"Derived speed for Object 1: {v1:.2f} m/s")
    print(f"Derived speed for Object 2: {v2:.2f} m/s")
else:
    raise ValueError("Invalid mode. Use 'speed' or 'force'.")

# === TIME CALCULATIONS ===
def flight_time(v, h, theta):
    vy = v * np.sin(theta)
    discriminant = vy**2 + 2 * g * h
    return (vy + np.sqrt(discriminant)) / g

t1_land = flight_time(v1, height, theta)
t2_land = flight_time(v2, height, theta)
t_max = max(t1_land, t2_land)
t = np.linspace(0, t_max, 500)

# === TRAJECTORIES ===
def trajectory(v, theta, h):
    x = v * np.cos(theta) * t
    y = h + v * np.sin(theta) * t - 0.5 * g * t**2
    y = np.maximum(y, 0)
    return x, y

x1, y1 = trajectory(v1, theta, height)
x2, y2 = trajectory(v2, theta, height)

# === SPEED & ACCELERATION ===
def velocity_and_acceleration(v, theta):
    vx = v * np.cos(theta)
    vy = v * np.sin(theta) - g * t
    v_total = np.sqrt(vx**2 + vy**2)
    a_y = -g * np.ones_like(t)
    a_total = np.sqrt(0**2 + a_y**2)
    return v_total, a_total, vx, vy

v1_t, a1_t, vx1, vy1 = velocity_and_acceleration(v1, theta)
v2_t, a2_t, vx2, vy2 = velocity_and_acceleration(v2, theta)

# === ENERGIES ===
def energy(m, v, y):
    ke = 0.5 * m * v**2
    pe = m * g * y
    return ke, pe

ke1, pe1 = energy(mass1, v1_t, y1)
ke2, pe2 = energy(mass2, v2_t, y2)

# === IMPACT SPEEDS ===
# Define a small threshold instead of exact 0 to avoid numerical issues
epsilon = 1e-3

idx1 = np.where(y1 <= epsilon)[0]
idx2 = np.where(y2 <= epsilon)[0]

impact_speed1 = v1_t[idx1[0]] if len(idx1) > 0 else v1_t[-1]
impact_speed2 = v2_t[idx2[0]] if len(idx2) > 0 else v2_t[-1]


# === PLOTTING ===

plt.figure(figsize=(10, 6))
plt.plot(x1, y1, label="Object 1")
plt.plot(x2, y2, label="Object 2")
plt.title("Trajectory")
plt.xlabel("Horizontal Distance (m)")
plt.ylabel("Height (m)")
plt.legend()
plt.grid(True)
plt.axis('equal')
plt.show()

plt.figure(figsize=(10, 6))
plt.plot(t, v1_t, label="Speed Object 1")
plt.plot(t, v2_t, label="Speed Object 2")
plt.title("Total Speed vs Time")
plt.xlabel("Time (s)")
plt.ylabel("Speed (m/s)")
plt.legend()
plt.grid(True)
plt.show()

plt.figure(figsize=(10, 6))
plt.plot(t, a1_t, label="Acceleration Object 1")
plt.plot(t, a2_t, label="Acceleration Object 2")
plt.title("Acceleration vs Time")
plt.xlabel("Time (s)")
plt.ylabel("Acceleration (m/s²)")
plt.legend()
plt.grid(True)
plt.show()

plt.figure(figsize=(10, 6))
plt.plot(t, ke1, label="Kinetic Energy Object 1")
plt.plot(t, ke2, label="Kinetic Energy Object 2")
plt.plot(t, pe1, '--', label="Potential Energy Object 1")
plt.plot(t, pe2, '--', label="Potential Energy Object 2")
plt.title("Energy vs Time")
plt.xlabel("Time (s)")
plt.ylabel("Energy (J)")
plt.legend()
plt.grid(True)
plt.show()

# === RESULTS ===
print("\n--- RESULTS ---")
print(f"Time to land (Object 1): {t1_land:.2f} s")
print(f"Time to land (Object 2): {t2_land:.2f} s")
print(f"Impact speed (Object 1): {impact_speed1:.2f} m/s")
print(f"Impact speed (Object 2): {impact_speed2:.2f} m/s")

# === EXPLANATION ===
print("""
Why doesn't speed reach 0?
- If the object is launched upward, the vertical speed (vy) decreases due to gravity.
- But total speed (vector sum of vx and vy) doesn't necessarily reach 0 unless the object stops at its peak with vx = 0 too — which it doesn't in projectile motion.

Acceleration:
- Gravity is always downward and constant at -9.81 m/s².
- If upward motion: gravity resists (reduces vy)
- If downward motion: gravity increases vy
- Acceleration vector does not change unless forces change (we ignore air friction here)

Galileo's Principle:
- Mass doesn’t affect the trajectory (only initial speed, angle, and gravity do).
""")

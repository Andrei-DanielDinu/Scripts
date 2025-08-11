#!/usr/bin/env python3
"""
Planet escape/orbit table + optional simulation.

- Prints orbital & escape velocities for Solar System bodies (Pluto included).
- Optionally simulates a user-chosen launch (2D RK4, gravity-only).
- Simulation uses safe dt chosen from estimated orbital period to avoid huge arrays.
"""
import math
import numpy as np
import matplotlib.pyplot as plt

# Planetary data (mass in kg, radius in m)
planets = {
    "Mercury": {"mass": 3.3011e23, "radius": 2.4397e6},
    "Venus":   {"mass": 4.8675e24, "radius": 6.0518e6},
    "Earth":   {"mass": 5.97237e24, "radius": 6.371e6},
    "Mars":    {"mass": 6.4171e23, "radius": 3.3895e6},
    "Jupiter": {"mass": 1.8982e27, "radius": 6.9911e7},
    "Saturn":  {"mass": 5.6834e26, "radius": 5.8232e7},
    "Uranus":  {"mass": 8.6810e25, "radius": 2.5362e7},
    "Neptune": {"mass": 1.02413e26, "radius": 2.4622e7},
    "Pluto":   {"mass": 1.303e22, "radius": 1.1883e6}
}

G = 6.67430e-11

def v_escape(M, R):
    return math.sqrt(2 * G * M / R)

def v_orbital(M, R):
    return math.sqrt(G * M / R)

# Print table
print("Planetary orbital & escape velocity table (at surface level):")
print("{:<10s} {:>15s} {:>15s} {:>12s}".format("Planet", "v_orb (m/s)", "v_esc (m/s)", "Radius (km)"))
print("-"*56)
for name, data in planets.items():
    M = data["mass"]
    R = data["radius"]
    vo = v_orbital(M, R)
    ve = v_escape(M, R)
    print(f"{name:<10s} {vo:15.2f} {ve:15.2f} {R/1000:12.1f}")
print()

# Ask user whether to run a simulation
do_sim = input("Do you want to simulate a launch for one planet? (y/n): ").strip().lower()
if do_sim not in ("y", "yes"):
    print("Done. Table printed. Re-run with simulation if needed.")
    raise SystemExit(0)

# --- Simulation parameters from user ---
planet_choice = input(f"Choose planet {list(planets.keys())}: ").strip()
if planet_choice not in planets:
    print("Invalid planet. Exiting.")
    raise SystemExit(1)

M = planets[planet_choice]["mass"]
R = planets[planet_choice]["radius"]
h = float(input("Launch height above surface (m): "))
angle_deg = float(input("Launch angle above local horizontal (deg): "))
mode = input("Initial input mode ('speed' or 'force'): ").strip().lower()
if mode == "speed":
    v0 = float(input("Initial speed (m/s): "))
elif mode == "force":
    F = float(input("Initial applied force (N): "))
    m_obj = float(input("Object mass (kg): "))
    # simple impulse model: delta_v = F * dt_impulse / m_obj; assume dt_impulse = 1s (user can change)
    dt_impulse = 1.0
    v0 = F * dt_impulse / m_obj
    print(f"Approximated initial speed from impulse: {v0:.6f} m/s")
else:
    print("Invalid mode. Exiting.")
    raise SystemExit(1)

angle_rad = math.radians(angle_deg)
r0 = R + h
mu = G * M

# Analytical velocities
v_circ = v_orbital(M, r0)
v_esc = v_escape(M, r0)

print(f"\nPlanet: {planet_choice} (M={M:.3e} kg, R={R/1000:.1f} km)")
print(f"Launch radius r0 = R + h = {r0:.3e} m")
print(f"Analytic circular speed at r0: {v_circ:.3f} m/s")
print(f"Analytic escape speed at r0:   {v_esc:.3f} m/s\n")

# --- RK4 integrator functions ---
def accel(px, py, mu_local):
    r2 = px*px + py*py
    r = math.sqrt(r2)
    if r == 0:
        return 0.0, 0.0
    a = -mu_local / (r2 * r)
    return a * px, a * py

def rk4_step(px, py, vx, vy, dt, mu_local):
    ax1, ay1 = accel(px, py, mu_local)
    k1vx, k1vy = ax1 * dt, ay1 * dt
    k1x, k1y = vx * dt, vy * dt

    ax2, ay2 = accel(px + 0.5 * k1x, py + 0.5 * k1y, mu_local)
    k2vx, k2vy = ax2 * dt, ay2 * dt
    k2x, k2y = (vx + 0.5 * k1vx) * dt, (vy + 0.5 * k1vy) * dt

    ax3, ay3 = accel(px + 0.5 * k2x, py + 0.5 * k2y, mu_local)
    k3vx, k3vy = ax3 * dt, ay3 * dt
    k3x, k3y = (vx + 0.5 * k2vx) * dt, (vy + 0.5 * k2vy) * dt

    ax4, ay4 = accel(px + k3x, py + k3y, mu_local)
    k4vx, k4vy = ax4 * dt, ay4 * dt
    k4x, k4y = (vx + k3vx) * dt, (vy + k3vy) * dt

    vx_new = vx + (k1vx + 2*k2vx + 2*k3vx + k4vx) / 6.0
    vy_new = vy + (k1vy + 2*k2vy + 2*k3vy + k4vy) / 6.0
    px_new = px + (k1x + 2*k2x + 2*k3x + k4x) / 6.0
    py_new = py + (k1y + 2*k2y + 2*k3y + k4y) / 6.0

    return px_new, py_new, vx_new, vy_new

# --- Simulation routine (uses lists to avoid huge prealloc) ---
def simulate(px0, py0, vx0, vy0, mu_local, R_local, dt, max_steps=200000):
    px = px0; py = py0; vx = vx0; vy = vy0
    positions = []
    speeds = []
    energies = []
    times = []

    for step in range(max_steps):
        t = step * dt
        r = math.hypot(px, py)
        v = math.hypot(vx, vy)
        # specific energy (per unit mass)
        eps = 0.5 * v*v - mu_local / r

        positions.append((px, py))
        speeds.append(v)
        energies.append(eps)
        times.append(t)

        # Crash detection
        if r <= R_local:
            return {"status": "CRASH", "positions": positions, "speeds": speeds, "energies": energies, "times": times}

        # Escape heuristic (moving away and positive energy)
        if (r > 20 * r0) and (eps > 0):
            return {"status": "ESCAPE", "positions": positions, "speeds": speeds, "energies": energies, "times": times}

        # Orbit completion detection: track angle sweep
        # We'll check every few steps to avoid cost; unwrap angle array and check delta
        if step % 500 == 0 and len(positions) > 1000:
            arr = np.array(positions)
            ang = np.unwrap(np.arctan2(arr[:,1], arr[:,0]))
            delta_ang = ang[-1] - ang[0]
            if abs(delta_ang) >= 2*math.pi:
                return {"status": "ORBIT", "positions": positions, "speeds": speeds, "energies": energies, "times": times, "delta_ang": delta_ang}

        # advance one RK4 step
        px, py, vx, vy = rk4_step(px, py, vx, vy, dt, mu_local)

    return {"status": "TIMEOUT", "positions": positions, "speeds": speeds, "energies": energies, "times": times}

# choose dt based on estimated orbital period at r0
T_est = 2 * math.pi * math.sqrt(r0**3 / (G * M))
# choose dt so that ~2000 steps per period, but not smaller than 0.01s
dt = max(0.01, T_est / 2000.0)
max_steps = int(min(2000000, max(200000, 5 * (T_est / dt))))  # cap total steps

print(f"Using dt = {dt:.4f} s, max_steps = {max_steps}, estimated orbital period T ≈ {T_est:.1f} s")

# initial conditions: place launch point at (r0, 0), local horizontal = +y direction
px0 = r0
py0 = 0.0
vx0 = v0 * math.cos(angle_rad)
vy0 = v0 * math.sin(angle_rad)

print("\nSimulating...")
result = simulate(px0, py0, vx0, vy0, mu, R, dt, max_steps=max_steps)
status = result["status"]
print("Simulation status:", status)
if status == "ORBIT":
    print(f"Completed revolution; angle swept = {result.get('delta_ang', 0):.3f} rad")
elif status == "TIMEOUT":
    print("Simulation timed out (max steps reached).")
elif status == "ESCAPE":
    print("Object escaped (radial distance grew large and energy > 0).")
elif status == "CRASH":
    print("Object crashed on the surface.")

# print simple summary of initial energies
v_initial = v0
eps_initial = 0.5 * v_initial**2 - mu / r0
print(f"\nInitial speed = {v_initial:.3f} m/s")
print(f"Initial specific energy ε = {eps_initial:.3e} J/kg (ε >= 0 => unbound)")

# Plot results
positions = np.array(result["positions"])
times = np.array(result["times"])
speeds = np.array(result["speeds"])
energies = np.array(result["energies"])

# Plot trajectory and planet
fig, axs = plt.subplots(2, 2, figsize=(12,10))

ax = axs[0,0]
ax.plot(positions[:,0]/1000.0, positions[:,1]/1000.0, '-b', label='Trajectory')
# planet circle
theta = np.linspace(0, 2*math.pi, 400)
ax.fill((R*np.cos(theta))/1000.0, (R*np.sin(theta))/1000.0, color='saddlebrown', alpha=0.6, label=planet_choice)
ax.set_title(f"Trajectory around {planet_choice}")
ax.set_xlabel("x (km)"); ax.set_ylabel("y (km)")
ax.set_aspect('equal')
ax.legend()
ax.grid(True)

# altitude vs time
ax2 = axs[0,1]
r_arr = np.hypot(positions[:,0], positions[:,1])
ax2.plot(times, r_arr/1000.0 - R/1000.0)
ax2.set_title("Altitude above surface vs time")
ax2.set_xlabel("Time (s)"); ax2.set_ylabel("Altitude (km)")
ax2.grid(True)

# speed vs time
ax3 = axs[1,0]
ax3.plot(times, speeds)
ax3.axhline(v_circ, color='g', linestyle='--', label=f"v_circ={v_circ:.1f} m/s")
ax3.axhline(v_esc, color='r', linestyle='--', label=f"v_esc={v_esc:.1f} m/s")
ax3.set_title("Speed vs time")
ax3.set_xlabel("Time (s)"); ax3.set_ylabel("Speed (m/s)")
ax3.legend()
ax3.grid(True)

# specific energy vs time
ax4 = axs[1,1]
ax4.plot(times, energies)
ax4.axhline(0.0, color='r', linestyle='--', label='ε=0 (escape threshold)')
ax4.set_title("Specific mechanical energy (ε) vs time")
ax4.set_xlabel("Time (s)"); ax4.set_ylabel("ε (J/kg)")
ax4.legend()
ax4.grid(True)

plt.tight_layout()
plt.show()

print("\nDone. You can re-run the script to test other planets or input modes.")

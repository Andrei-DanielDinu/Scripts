"""
projectile_and_recoil_simulator.py

Comprehensive educational simulator (physics-only) for a 2D projectile with quadratic air drag
and steady wind, plus a detailed recoil/impulse calculator for the shooter.

Features:
- Numerically integrate projectile motion with quadratic drag using scipy.solve_ivp.
- Compute momentum of projectile and propellant gases, recoil velocity, recoil energy,
  impulse, average and peak forces on the shooter.
- Model a realistic force-vs-time profile for recoil (half-sine or triangular) and plot it.
- Print step-by-step results and draw publication-quality plots.

NOT INTENDED FOR WEAPON USE. This is an educational physics tool only.

Formulas (implemented below):
- Relative velocity to air: u = v - v_wind
- Quadratic drag: F_drag = -0.5 * rho * C_d * A * u * u_vector
- Total forward momentum: p_total = m_b * v_b + m_g * v_g
  (we allow user to set gas mass m_g and gas velocity factor f_g such that v_g = f_g * v_b)
- Recoil velocity: v_recoil = p_total / M_system
- Recoil energy: E_recoil = 0.5 * M_system * v_recoil**2
- Impulse: J = p_total
- Average force: F_avg = J / delta_t
- Peak force (for half-sine profile): F_peak = J * (pi / (2 * delta_t)) = (pi/2) * F_avg

Requirements: numpy, scipy, matplotlib

Usage: run this script with Python 3.8+.

"""

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import math

# ----------------------------- Helper: pretty print -----------------------------
def fmt(x, unit=""):
    return f"{x:.4g} {unit}" if isinstance(x, (float, int)) else str(x)

# ----------------------------- Physical & projectile defaults -----------------------------
# Environment
g = 9.80665               # gravity (m/s^2)
rho = 1.225               # air density at sea level (kg/m^3)

# Projectile (editable)
m_b = 0.010               # bullet mass (kg) -- 10 g (example)
Cd = 0.295                # drag coefficient (typical bullet/streamlined ~0.2-0.3)
diameter = 0.007         # m (7 mm)
A = math.pi * (diameter/2)**2

# Launch conditions (editable)
v0 = 1000.0               # muzzle velocity (m/s) -- choose responsibly
angle_deg = 30.0          # launch angle in degrees
theta = math.radians(angle_deg)
vx0 = v0 * math.cos(theta)
vy0 = v0 * math.sin(theta)

# Wind (air velocity)
wind_speed = 0.0          # m/s (positive = tailwind in +x)
wind_dir_deg = 0.0
wx = wind_speed * math.cos(math.radians(wind_dir_deg))
wy = wind_speed * math.sin(math.radians(wind_dir_deg))
v_wind = np.array([wx, wy])

# Shooter + gun system
M_system = 80.0           # kg (combined mass of shooter+gun that absorbs recoil momentum)

# Propellant gas model (editable)
m_g = 0.003               # kg of gas exiting the barrel (example: 3 g)
f_g = 1.5                 # gas exit speed factor relative to muzzle velocity (v_g = f_g * v0)

# Impulse time model (editable)
# Typical impulse times for shoulder-fired weapons are on the order of 0.005-0.02 s.
# Choose delta_t consistent with the system; shorter = sharper, larger peak forces.
delta_t = 0.010          # seconds (duration over which shoulder receives the impulse)

# Force profile shape: 'half-sine' or 'triangular'
force_profile = 'half-sine'  # options: 'half-sine', 'triangular'

# ----------------------------- ODE system (projectile motion) -----------------------------

def projectile_deriv(t, state):
    # state = [x, y, vx, vy]
    x, y, vx, vy = state
    # relative velocity w.r.t. the air
    ux = vx - v_wind[0]
    uy = vy - v_wind[1]
    u = math.hypot(ux, uy)
    if u == 0.0:
        drag_x = drag_y = 0.0
    else:
        coeff = 0.5 * rho * Cd * A
        drag_x = -coeff * u * ux
        drag_y = -coeff * u * uy

    ax = drag_x / m_b
    ay = -g + (drag_y / m_b)
    return [vx, vy, ax, ay]

# Stop event when projectile returns to ground (y=0) going downwards

def hit_ground_event(t, state):
    return state[1]

hit_ground_event.terminal = True
hit_ground_event.direction = -1

# Integrate
initial_state = [0.0, 0.0, vx0, vy0]
max_t = 300.0  # safety cap on integration time
sol = solve_ivp(projectile_deriv, (0.0, max_t), initial_state, events=hit_ground_event, max_step=0.02, rtol=1e-8, atol=1e-9)

t = sol.t
x = sol.y[0]
y = sol.y[1]
vx = sol.y[2]
vy = sol.y[3]
speed = np.hypot(vx, vy)

# ----------------------------- Recoil calculations -----------------------------
# Momentum of bullet
p_bullet = m_b * v0  # note: we use muzzle speed v0 for approximation of forward momentum
# Momentum from gases
v_g = f_g * v0
p_gases = m_g * v_g
# Total forward momentum (conserved, opposite sign for shooter)
p_total = p_bullet + p_gases
# Recoil velocity of the shooter+gun system
v_recoil = p_total / M_system
# Recoil kinetic energy
E_recoil = 0.5 * M_system * v_recoil**2
# Impulse and average force
impulse = p_total
F_avg = impulse / delta_t
# Peak force depending on profile
if force_profile == 'half-sine':
    # For a half-sine F(t)=F_peak*sin(pi*t/delta_t) over 0..delta_t, impulse J = (2/ pi) * F_peak * delta_t
    # So F_peak = J * (pi / (2*delta_t)) = (pi/2) * F_avg
    F_peak = impulse * (math.pi / (2.0 * delta_t))
elif force_profile == 'triangular':
    # Triangular pulse (rise then fall) with peak F_peak has impulse J = 0.5 * F_peak * delta_t
    F_peak = 2.0 * F_avg
else:
    F_peak = F_avg

# ----------------------------- Outputs (printed neatly) -----------------------------
print("\n=== Projectile simulation summary ===")
print(f"Initial speed v0 = {v0:.3f} m/s, angle = {angle_deg:.2f} deg")
print(f"Bullet mass m_b = {m_b:.4f} kg, Cd = {Cd:.3f}, diameter = {diameter:.4f} m, area A = {A:.6e} m^2")
print(f"Air density rho = {rho} kg/m^3, wind = ({wx:.2f}, {wy:.2f}) m/s")
print(f"Range (downrange when y hits 0) = {x[-1]:.3f} m")
print(f"Time of flight = {t[-1]:.4f} s, max height = {np.max(y):.4f} m")

print("\n=== Recoil / impulse summary ===")
print(f"Bullet momentum p_bullet = m_b * v0 = {p_bullet:.6f} kg·m/s")
print(f"Gases momentum p_gases = m_g * v_g = {p_gases:.6f} kg·m/s (v_g = f_g * v0, f_g={f_g})")
print(f"Total forward momentum p_total = {p_total:.6f} kg·m/s")
print(f"Recoil velocity v_recoil = p_total / M_system = {v_recoil:.6f} m/s")
print(f"Recoil kinetic energy E_recoil = {E_recoil:.6f} J")
print(f"Impulse J = {impulse:.6f} N·s")
print(f"Assumed impulse duration delta_t = {delta_t:.6f} s")
print(f"Average force F_avg = J / delta_t = {F_avg:.3f} N (~{F_avg/9.80665:.3f} kgf)")
print(f"Modeled peak force F_peak = {F_peak:.3f} N (~{F_peak/9.80665:.3f} kgf) using profile '{force_profile}'")

# ----------------------------- Construct force-vs-time profile for plotting -----------------------------
# We'll model F(t) starting at t=0 (shot) to t=delta_t; after that force is zero (for recoil impulse)
num_pts = 500
t_force = np.linspace(0.0, delta_t, num_pts)
if force_profile == 'half-sine':
    F_t = F_peak * np.sin(math.pi * t_force / delta_t)
elif force_profile == 'triangular':
    # rise to peak at delta_t/2 then fall
    F_t = np.where(t_force <= delta_t/2.0,
                   (2.0 * F_peak / delta_t) * t_force,
                   (2.0 * F_peak / delta_t) * (delta_t - t_force))
else:
    F_t = np.ones_like(t_force) * F_avg

# ----------------------------- Plots -----------------------------
plt.rcParams.update({'figure.dpi': 120})

# Trajectory
plt.figure(figsize=(10,6))
plt.plot(x, y, linewidth=2)
plt.xlabel('Downrange distance (m)')
plt.ylabel('Height (m)')
plt.title(f'Trajectory — v0={v0} m/s, angle={angle_deg}°, Cd={Cd}, m_b={m_b} kg')
plt.grid(True)

# Speed & vertical velocity vs time
plt.figure(figsize=(10,4))
plt.plot(t, speed, label='speed (m/s)')
plt.plot(t, vy, label='vertical velocity (m/s)')
plt.xlabel('Time (s)')
plt.legend()
plt.grid(True)

# Downrange vs time
plt.figure(figsize=(10,4))
plt.plot(t, x, label='downrange x(t)')
plt.xlabel('Time (s)')
plt.ylabel('Downrange distance (m)')
plt.grid(True)

# Force vs time (recoil)
plt.figure(figsize=(10,4))
plt.plot(t_force * 1000.0, F_t, linewidth=2)  # convert time to milliseconds for readability
plt.xlabel('Time since shot (ms)')
plt.ylabel('Force on shoulder (N)')
plt.title(f'Recoil force profile — impulse J={impulse:.3f} N·s, duration={delta_t*1000:.2f} ms')
plt.grid(True)

plt.tight_layout()
plt.show()

# ----------------------------- Example additional print for a 2-mile question -----------------------------
# If user is curious about a hypothetical projectile "that goes 2 miles", we can show how to plug numbers.
# Note: range depends strongly on drag and angle; achieving 2 miles (3218.69 m) in real air may require very high v0.

example_range_target = 3218.69  # meters (2 miles)
print('\nNote: achieving a 2-mile downrange distance depends on drag, angle, and v0.\n'
      'This script prints the actual range for the provided v0 and parameters above.\n'
      'If you want, I can add a solver that searches for the necessary v0 to reach a target range '\
      '(it will be an iterative process that accounts for drag).')

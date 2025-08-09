import math
import numpy as np
import matplotlib.pyplot as plt

# ------------------ USER INPUTS ------------------
# Pendulum inputs
L_p = float(input("Pendulum arm length (m): "))
theta0_deg = float(input("Initial pendulum angle (degrees): "))
efficiency_p = float(input("Pendulum efficiency per swing (0-1): "))

# Ball inputs
m_b = float(input("Ball mass (kg): "))
r_b = float(input("Ball radius (m): "))
ramp_mode = input("Ramp input type ('height+base' or 'height+angle'): ").strip().lower()
h_r = float(input("Ramp height (m): "))

if ramp_mode == "height+base":
    base_r = float(input("Ramp base length (m): "))
    ramp_length = math.sqrt(h_r**2 + base_r**2)
    angle_r = math.atan2(h_r, base_r)
elif ramp_mode == "height+angle":
    angle_r = math.radians(float(input("Ramp angle (degrees): ")))
    ramp_length = h_r / math.sin(angle_r)
else:
    raise ValueError("Invalid ramp mode.")

efficiency_ramp = float(input("Ramp-to-ground kinetic energy efficiency (0-1): "))
rolling_resist = float(input("Rolling resistance coefficient (ground): "))

g = 9.81

# ------------------ PENDULUM SIMULATION ------------------
theta0 = math.radians(theta0_deg)
T0 = 2 * math.pi * math.sqrt(L_p / g)  # small angle approximation

pendulum_times = []
pendulum_thetas = []
pendulum_energies_pot = []
pendulum_energies_kin = []
time_elapsed = 0
theta_amp = theta0
swing_count = 0

while theta_amp > math.radians(0.5):  # stop when amplitude small
    swing_time = T0 * (math.sin(theta_amp / 2) / math.sin(theta0 / 2))
    t = np.linspace(0, swing_time, 100)
    theta_t = theta_amp * np.cos((2 * math.pi / swing_time) * t)
    v_t = -theta_amp * (2 * math.pi / swing_time) * np.sin((2 * math.pi / swing_time) * t)

    for i in range(len(t)):
        pendulum_times.append(time_elapsed + t[i])
        pendulum_thetas.append(theta_t[i])
        h = L_p * (1 - math.cos(theta_t[i]))
        v = L_p * v_t[i]
        pot = m_b * g * h
        kin = 0.5 * m_b * v**2
        pendulum_energies_pot.append(pot)
        pendulum_energies_kin.append(kin)

    time_elapsed += swing_time
    theta_amp *= efficiency_p
    swing_count += 1

print(f"Pendulum swings until stop: {swing_count}")
print(f"Total pendulum time: {time_elapsed:.2f} s")

# ------------------ BALL RAMP + GROUND SIMULATION ------------------
# Ramp phase
a_ramp = g * math.sin(angle_r) * (1 / (1 + (2/5)))  # rolling sphere
t_ramp = math.sqrt(2 * ramp_length / a_ramp)
v_ramp_end = a_ramp * t_ramp

# Ground phase (rolling resistance decelerates)
v_ground_start = v_ramp_end * math.sqrt(efficiency_ramp)
a_ground = -rolling_resist * g
t_ground = abs(v_ground_start / a_ground)
dist_ground = v_ground_start * t_ground + 0.5 * a_ground * t_ground**2

# Ball time array
ball_times = np.concatenate((np.linspace(0, t_ramp, 100), 
                              np.linspace(t_ramp, t_ramp + t_ground, 100)))
ball_speeds = []
ball_pot_energy = []
ball_kin_energy = []

for t in ball_times:
    if t <= t_ramp:
        s = 0.5 * a_ramp * t**2
        v = a_ramp * t
        h = h_r - s * math.sin(angle_r)
    else:
        t_rel = t - t_ramp
        v = max(v_ground_start + a_ground * t_rel, 0)
        h = 0
    ball_speeds.append(v)
    pot = m_b * g * h
    kin = 0.5 * m_b * v**2
    ball_pot_energy.append(pot)
    ball_kin_energy.append(kin)

print(f"Ramp length: {ramp_length:.2f} m")
print(f"Ramp time: {t_ramp:.3f} s")
print(f"Speed at ramp end: {v_ramp_end:.3f} m/s")
print(f"Ground distance: {dist_ground:.2f} m")
print(f"Total ball travel time: {t_ramp + t_ground:.2f} s")

# ------------------ ENERGY EQUALITY CHECK ------------------
energy_equal_points = []
min_len = min(len(pendulum_times), len(ball_times))
for i in range(min_len):
    if abs(pendulum_energies_pot[i] - pendulum_energies_kin[i]) < 1e-2:
        energy_equal_points.append(("Pendulum P=K", pendulum_times[i]))
    if abs(ball_pot_energy[i] - ball_kin_energy[i]) < 1e-2:
        energy_equal_points.append(("Ball P=K", ball_times[i]))
    total_pend = pendulum_energies_pot[i] + pendulum_energies_kin[i]
    total_ball = ball_pot_energy[i] + ball_kin_energy[i]
    if abs(total_pend - total_ball) < 1e-2:
        energy_equal_points.append(("Pendulum=Ball", pendulum_times[i]))

print("\nEnergy equalities found:")
for label, t in energy_equal_points:
    print(f"{label} at t={t:.3f}s")

# ------------------ PLOTS ------------------
fig, axs = plt.subplots(3, 2, figsize=(12, 10))

# Pendulum theta
axs[0, 0].plot(pendulum_times, np.degrees(pendulum_thetas))
axs[0, 0].set_title("Pendulum Angle vs Time")
axs[0, 0].set_xlabel("Time (s)")
axs[0, 0].set_ylabel("Angle (deg)")

# Ball speed
axs[0, 1].plot(ball_times, ball_speeds)
axs[0, 1].set_title("Ball Speed vs Time")
axs[0, 1].set_xlabel("Time (s)")
axs[0, 1].set_ylabel("Speed (m/s)")

# Pendulum energies
axs[1, 0].plot(pendulum_times, pendulum_energies_pot, label="Potential")
axs[1, 0].plot(pendulum_times, pendulum_energies_kin, label="Kinetic")
axs[1, 0].set_title("Pendulum Energies")
axs[1, 0].legend()

# Ball energies
axs[1, 1].plot(ball_times, ball_pot_energy, label="Potential")
axs[1, 1].plot(ball_times, ball_kin_energy, label="Kinetic")
axs[1, 1].set_title("Ball Energies")
axs[1, 1].legend()

# Pendulum arc length
arc_lengths = [L_p * abs(theta) for theta in pendulum_thetas]
axs[2, 0].plot(pendulum_times, arc_lengths)
axs[2, 0].set_title("Pendulum Arc Length vs Time")

# Ramp diagram
axs[2, 1].set_aspect('equal', adjustable='box')
axs[2, 1].plot([0, ramp_length * math.cos(angle_r)], [h_r, 0], 'r-')
axs[2, 1].plot(ramp_length*math.cos(angle_r), 0, 'bo')
axs[2, 1].annotate("Impact", (ramp_length*math.cos(angle_r), 0))
axs[2, 1].set_title("Ramp Diagram")

plt.tight_layout()
plt.show()

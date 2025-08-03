import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colormaps

# === Planetary orbital data ===
planet_data = {
    "Mercury": {"a": 0.39, "period": 0.24},
    "Venus": {"a": 0.72, "period": 0.615},
    "Earth": {"a": 1.0, "period": 1.0},
    "Mars": {"a": 1.524, "period": 1.881},
    "Jupiter": {"a": 5.203, "period": 11.86},
    "Saturn": {"a": 9.537, "period": 29.46},
}

# === Ask the user for a reference planet ===
print("Available planets:", ", ".join(planet_data.keys()))
reference_planet = input("Choose a reference planet to observe from (e.g. Earth): ").strip().capitalize()

if reference_planet not in planet_data:
    print(f"Error: '{reference_planet}' is not a valid planet.")
    exit()

# === Simulation parameters ===
longest_period = max(p["period"] for p in planet_data.values())
total_days = int(longest_period * 365.25 * 2)
days = np.arange(0, total_days, 1)

# === Compute positions for all planets ===
planet_positions = {}
for name, data in planet_data.items():
    theta = 2 * np.pi * (days / (data["period"] * 365.25))
    x = data["a"] * np.cos(theta)
    y = data["a"] * np.sin(theta)
    planet_positions[name] = np.column_stack((x, y))

# === Compute reference POV ===
reference_positions = planet_positions[reference_planet]

# === Plot settings ===
fig, ax = plt.subplots(figsize=(12, 6))
cmap = colormaps.get_cmap('tab10')
planet_names = list(planet_data.keys())

for i, planet in enumerate(planet_names):
    if planet == reference_planet:
        continue

    # Relative position
    rel_positions = planet_positions[planet] - reference_positions
    angles = np.arctan2(rel_positions[:, 1], rel_positions[:, 0])
    angles_unwrapped = np.unwrap(angles)
    angles_deg = np.degrees(angles_unwrapped)

    # Detect retrograde motion (angle decreasing)
    retro_mask = np.diff(angles_deg) < 0
    retro_days = days[1:][retro_mask]

    color = cmap(i % 10)  # safely cycle colors
    ax.plot(days, angles_deg, label=f"{planet} from {reference_planet}", color=color)

    for day in retro_days:
        ax.axvline(x=day, color=color, linestyle="--", alpha=0.1)

ax.set_title(f"Apparent Motion from {reference_planet}'s POV")
ax.set_xlabel("Days")
ax.set_ylabel("Ecliptic Longitude (degrees)")
ax.grid(True)
ax.legend()
plt.tight_layout()
plt.show()

# === Educational Explanation ===
print(f"\n--- What is Retrograde Motion? ---")
print(f"From the perspective of {reference_planet}, other planets may appear to move backwards in the sky.")
print("This illusion happens when Earth (or any observer planet) overtakes a slower-moving planet in orbit.")
print("The apparent reversal is due to the geometry of the orbits, not actual reversal of motion.")

print(f"\n--- Summary of Retrograde Motion (as seen from {reference_planet}) ---")
for planet in planet_names:
    if planet == reference_planet:
        continue
    angles = np.unwrap(np.arctan2(
        planet_positions[planet][:, 1] - reference_positions[:, 1],
        planet_positions[planet][:, 0] - reference_positions[:, 0]
    ))
    if np.any(np.diff(angles) < 0):
        print(f"✓ {planet} shows retrograde motion.")
    else:
        print(f"✗ {planet} does NOT show retrograde motion.")

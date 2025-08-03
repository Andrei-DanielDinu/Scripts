import math
import matplotlib.pyplot as plt

# Planetary radii (km)
planet_radii = {
    "Mercury": 2439.7,
    "Venus": 6051.8,
    "Earth": 6371.0,
    "Mars": 3389.5,
    "Jupiter": 69911,
    "Saturn": 58232,
    "Uranus": 25362,
    "Neptune": 24622,
}

# === User Inputs ===
stick_height = float(input("Stick height in meters: "))
error_margin = float(input("Target error margin in % (e.g., 0.5): "))
planet = input("Planet name: ").capitalize()
daily_walk_km = float(input("Distance walked per day (km): "))

if planet not in planet_radii:
    print(f"Unknown planet: {planet}")
    exit()

true_radius = planet_radii[planet]
print(f"\n🪐 Simulating on {planet} (true radius: {true_radius} km)")

# === Simulation ===
distance_walked = 0
estimates = []
errors = []
distances = []
angles_deg = []
shadow_lengths = []

latitude_deg = 0.0  # starting at equator
day = 0
measured = False
reference_angle = None
reference_lat = None

print("\nDay   Dist(km)   Shadow(cm)   Sun Ang(°)   Est.R(km)    Error(%)")

while True:
    day += 1
    distance_walked += daily_walk_km
    latitude_deg += (daily_walk_km / (2 * math.pi * true_radius)) * 360
    latitude_rad = math.radians(latitude_deg)

    # Measure shadow at noon
    shadow_length = stick_height * math.tan(latitude_rad)
    angle = math.degrees(math.atan(shadow_length / stick_height))

    if not measured:
        # Store first real measurement above equator
        reference_angle = angle
        reference_lat = latitude_deg
        measured = True
        continue

    delta_angle_rad = math.radians(angle - reference_angle)
    delta_lat_rad = math.radians(latitude_deg - reference_lat)

    if delta_angle_rad == 0:
        estimated_radius = float('inf')
    else:
        estimated_radius = daily_walk_km * day / delta_angle_rad

    error = abs((estimated_radius - true_radius) / true_radius) * 100

    print(f"{day:<5} {distance_walked:<10.2f} {shadow_length*100:<12.1f} {90-angle:<12.2f} {estimated_radius:<12.1f} {error:<8.2f}")

    # Track data for plotting
    estimates.append(estimated_radius)
    errors.append(error)
    distances.append(distance_walked)
    angles_deg.append(90 - angle)
    shadow_lengths.append(shadow_length * 100)  # to cm

    if error <= error_margin:
        print(f"\n🎯 Accuracy reached at Day {day}, Distance: {distance_walked:.2f} km")
        break

# === Plotting ===
plt.figure(figsize=(12, 5))

plt.subplot(1, 2, 1)
plt.plot(distances, estimates, label="Estimated Radius", color='blue')
plt.axhline(true_radius, color='green', linestyle='--', label="True Radius")
plt.xlabel("Distance Walked (km)")
plt.ylabel("Radius Estimate (km)")
plt.title("Estimated Planet Radius")
plt.legend()
plt.grid(True)

plt.subplot(1, 2, 2)
plt.plot(distances, errors, label="Error %", color='red')
plt.axhline(error_margin, color='gray', linestyle='--', label="Target Error Margin")
plt.xlabel("Distance Walked (km)")
plt.ylabel("Error (%)")
plt.title("Error Margin Over Time")
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.show()

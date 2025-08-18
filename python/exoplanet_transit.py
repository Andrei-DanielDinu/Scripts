import numpy as np
import matplotlib.pyplot as plt

# ----------------------------
# Transit Simulation Function
# ----------------------------
def simulate_light_curve(star_radius, planet_radius, orbital_period, duration_days, noise_level=0.001):
    """
    Simulates a star's light curve when a planet transits in front of it.
    
    Parameters:
    - star_radius: radius of the star (arbitrary units, e.g. solar radii)
    - planet_radius: radius of the planet (same units as star)
    - orbital_period: orbital period of the planet in days
    - duration_days: how many days to simulate
    - noise_level: random noise to mimic observations
    """
    time = np.linspace(0, duration_days, 2000)  # time in days
    flux = np.ones_like(time)

    # Transit depth = fraction of star’s light blocked
    transit_depth = (planet_radius / star_radius) ** 2  

    # Phase of orbit
    phase = (time % orbital_period) / orbital_period

    # Simple box-shaped transit model
    in_transit = (phase > 0.49) & (phase < 0.51)
    flux[in_transit] -= transit_depth

    # Add observational noise
    flux_noisy = flux + np.random.normal(0, noise_level, size=len(time))

    return time, flux_noisy, transit_depth


# ----------------------------
# User Inputs
# ----------------------------
print("\n🔭 Exoplanet Transit Method Simulator 🔭\n")

star_radius = float(input("Enter star radius (e.g. 1.0 for Sun): "))
planet_radius = float(input("Enter planet radius (e.g. 0.1 for Jupiter-size relative to Sun): "))
orbital_period = float(input("Enter orbital period (in days, e.g. 10): "))
duration_days = float(input("Enter duration of observation (days, e.g. 30): "))

# ----------------------------
# Run Simulation
# ----------------------------
time, flux, transit_depth = simulate_light_curve(star_radius, planet_radius, orbital_period, duration_days)

# ----------------------------
# Explanations
# ----------------------------
print("\n📊 Results and Explanations 📊\n")

print(f"➡ Orbital Period: {orbital_period:.2f} days")
print("   Meaning: This is how long the planet takes to complete one orbit around its star.\n")

print(f"➡ Transit Depth: {transit_depth*100:.3f}% dimming of starlight")
print("   Meaning: This is how much the star's brightness decreases when the planet passes in front of it.\n")

# Estimate planet-to-star size ratio
size_ratio = planet_radius / star_radius
print(f"➡ Planet-to-Star Radius Ratio: {size_ratio:.3f}")
print("   Meaning: A larger ratio means the planet blocks more light. "
      "For example, Earth blocks ~0.01% of the Sun’s light, Jupiter ~1%.\n")

print("✅ By observing the periodic dips in brightness, astronomers can estimate:")
print("   - The planet’s orbital period (from spacing between dips).")
print("   - The planet’s relative size (from depth of dips).")
print("   - Even possible hints of multiple planets (if multiple dips appear with different patterns).\n")

# ----------------------------
# Plotting
# ----------------------------
plt.figure(figsize=(10, 5))
plt.scatter(time, flux, s=5, color="blue", alpha=0.6, label="Simulated Observations")
plt.xlabel("Time (days)")
plt.ylabel("Relative Brightness")
plt.title("Simulated Exoplanet Transit Light Curve")
plt.legend()
plt.show()

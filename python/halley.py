import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import newton

# Constants
AU = 1.496e11        # meters (astronomical unit)
DAY = 86400          # seconds in a day
YEAR = 365.25 * DAY   # seconds in a year
G = 6.67430e-11      # gravitational constant m^3/kg/s^2
M_sun = 1.989e30     # mass of the Sun in kg

# Orbital Elements class for 2D orbits
class Orbit2D:
    def __init__(self, a_AU, e, T_years, perihelion_time=0):
        """
        a_AU : semi-major axis in AU
        e : eccentricity
        T_years : orbital period in years
        perihelion_time : time of perihelion passage in seconds (epoch)
        """
        self.a = a_AU * AU
        self.e = e
        self.T = T_years * YEAR
        self.n = 2 * np.pi / self.T    # mean motion (rad/s)
        self.tau = perihelion_time    # perihelion passage time

    def kepler_E(self, M, tol=1e-10):
        """
        Solve Kepler's equation M = E - e sin(E) for eccentric anomaly E
        using Newton-Raphson method.
        """
        func = lambda E: E - self.e * np.sin(E) - M
        func_prime = lambda E: 1 - self.e * np.cos(E)
        # Initial guess
        E0 = M if self.e < 0.8 else np.pi
        E = newton(func, E0, fprime=func_prime, tol=tol)
        return E

    def position(self, t):
        """
        Compute (x,y) position at time t (seconds)
        relative to perihelion epoch.
        """
        M = self.n * (t - self.tau)  # mean anomaly
        M = M % (2 * np.pi)          # wrap between 0 and 2pi
        E = self.kepler_E(M)
        # True anomaly
        nu = 2 * np.arctan2(np.sqrt(1+self.e)*np.sin(E/2),
                            np.sqrt(1-self.e)*np.cos(E/2))
        # Radius
        r = self.a * (1 - self.e*np.cos(E))
        # Cartesian coordinates (Sun at origin)
        x = r * np.cos(nu)
        y = r * np.sin(nu)
        return np.array([x, y])

    def orbital_period_years(self):
        return self.T / YEAR


# Define Earth's orbit (approximate)
earth_orbit = Orbit2D(a_AU=1.000, e=0.0167, T_years=1.000)

# Halley's comet orbital parameters (2D approximation)
halley_orbit = Orbit2D(a_AU=17.8, e=0.967, T_years=75.3, perihelion_time=0)

# Simulation parameters
total_years = 100
dt = DAY  # time step in seconds (1 day)
steps = int(total_years * YEAR / dt)
times = np.linspace(0, total_years * YEAR, steps)

# Arrays to store positions and distances
earth_positions = np.zeros((steps, 2))
halley_positions = np.zeros((steps, 2))
distances = np.zeros(steps)

# Calculate positions and distances over time
for i, t in enumerate(times):
    earth_pos = earth_orbit.position(t)
    halley_pos = halley_orbit.position(t)
    earth_positions[i] = earth_pos
    halley_positions[i] = halley_pos
    distances[i] = np.linalg.norm(earth_pos - halley_pos)

# Find minimum distance and when it happens
min_index = np.argmin(distances)
min_distance_m = distances[min_index]
min_distance_AU = min_distance_m / AU
min_time_years = times[min_index] / YEAR

# Print results
print(f"Minimum distance between Earth and Halley's comet in {total_years} years:")
print(f"  Distance = {min_distance_AU:.4f} AU ({min_distance_m/1e9:.2f} million km)")
print(f"  Time = {min_time_years:.2f} years from start")

# Plot orbits in XY plane
plt.figure(figsize=(10,10))
plt.plot(earth_positions[:,0]/AU, earth_positions[:,1]/AU, label="Earth Orbit")
plt.plot(halley_positions[:,0]/AU, halley_positions[:,1]/AU, label="Halley's Comet Orbit")
plt.scatter([0], [0], color='yellow', s=200, label="Sun (focus)")
plt.scatter(earth_positions[min_index,0]/AU, earth_positions[min_index,1]/AU, color='blue', s=100, label="Earth at closest approach")
plt.scatter(halley_positions[min_index,0]/AU, halley_positions[min_index,1]/AU, color='red', s=100, label="Halley at closest approach")
plt.axis('equal')
plt.xlabel("X (AU)")
plt.ylabel("Y (AU)")
plt.title(f"2D Orbits of Earth and Halley's Comet (100 years)\nMin distance = {min_distance_AU:.4f} AU at {min_time_years:.2f} years")
plt.legend()
plt.grid(True)
plt.show()

# Plot distance over time
plt.figure(figsize=(12,5))
plt.plot(times / YEAR, distances / AU)
plt.xlabel("Time (years)")
plt.ylabel("Distance between Earth and Halley (AU)")
plt.title("Distance between Earth and Halley's comet over 100 years")
plt.grid(True)
plt.show()

# Show orbital period of Halley's comet
print(f"Halley's comet orbital period (Kepler's law) = {halley_orbit.orbital_period_years():.2f} years")

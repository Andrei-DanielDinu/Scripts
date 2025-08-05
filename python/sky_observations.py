from skyfield.api import load, Topos
from datetime import datetime
import matplotlib.pyplot as plt
import os

# Set observer location (Bucharest)
observer = Topos(latitude_degrees=44.4268, longitude_degrees=26.1025)

# Load updated ephemeris with planetary data (de440s has full coverage)
eph = load('de440s.bsp')

# Define celestial bodies
sun = eph['SUN']
moon = eph['MOON']
planets = {
    'Mercury': eph['MERCURY BARYCENTER'],
    'Venus': eph['VENUS BARYCENTER'],
    'Mars': eph['MARS BARYCENTER'],
    'Jupiter': eph['JUPITER BARYCENTER'],
    'Saturn': eph['SATURN BARYCENTER'],
    'Uranus': eph['URANUS BARYCENTER'],
    'Neptune': eph['NEPTUNE BARYCENTER'],
}


# Time setup
ts = load.timescale()
today = datetime.utcnow()
t = ts.utc(today.year, today.month, today.day)

# Observer setup
earth = eph['EARTH']
observer_position = earth + observer

# Observation data
astro_data = {}
moon_astrometric = observer_position.at(t).observe(moon).apparent()
moon_alt, moon_az, _ = moon_astrometric.altaz()
astro_data['Moon'] = (moon_alt.degrees, moon_az.degrees)

for name, body in planets.items():
    astrometric = observer_position.at(t).observe(body).apparent()
    alt, az, _ = astrometric.altaz()
    if alt.degrees > 0:
        astro_data[name] = (alt.degrees, az.degrees)

# Save to TXT file
log_dir = 'astro_logs'
os.makedirs(log_dir, exist_ok=True)
log_filename = os.path.join(log_dir, f"{today.strftime('%Y-%m-%d')}.txt")

with open(log_filename, 'w') as file:
    file.write(f"Date: {today.strftime('%Y-%m-%d')}\n")
    file.write(f"Location: Bucharest (Lat: 44.4268, Lon: 26.1025)\n")
    file.write("Visible Celestial Bodies at Sunrise:\n")
    for body, (alt, az) in astro_data.items():
        file.write(f"  {body:<8} - Altitude: {alt:.2f}°, Azimuth: {az:.2f}°\n")

# Polar plot
fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, polar=True)
for body, (alt, az) in astro_data.items():
    r = 90 - alt  # Inverted for plotting
    theta = az * (3.1415926535 / 180)
    ax.plot(theta, r, 'o', label=body)

ax.set_theta_zero_location('N')
ax.set_theta_direction(-1)
ax.set_rlim(0, 90)
plt.title(f"Visible Celestial Bodies at Sunrise ({today.strftime('%Y-%m-%d')}) - Bucharest")
plt.legend(loc='upper right', bbox_to_anchor=(1.3, 1.1))
plt.tight_layout()
plt.show()

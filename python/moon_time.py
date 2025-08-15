from skyfield.api import load, Topos, utc
from datetime import datetime, timedelta, timezone
import math

# --- CONFIGURATION ---
latitude = 40.7128   # Your latitude
longitude = -74.0060 # Your longitude
date = datetime.now(timezone.utc).date()  # Today's date in UTC

# --- LOAD ASTRONOMY DATA ---
ts = load.timescale()
eph = load('de421.bsp')  # Ephemeris data
moon = eph['moon']
earth = eph['earth']

# Observer's location
observer = earth + Topos(latitude_degrees=latitude, longitude_degrees=longitude)

# --- FUNCTION TO GET MOON PHASE ---
def moon_phase_percentage(t):
    sun = eph['sun']
    e = earth.at(t)
    _, slon, _ = e.observe(sun).apparent().ecliptic_latlon()
    _, mlon, _ = e.observe(moon).apparent().ecliptic_latlon()
    phase_angle = (mlon.degrees - slon.degrees) % 360
    return (1 - math.cos(math.radians(phase_angle))) / 2 * 100

# --- FIND BEST MATCH FOR GIVEN ALTITUDE + PHASE ---
def estimate_time(target_altitude, target_phase, tolerance=1):
    best_match = None
    best_diff = float('inf')
    start_time = datetime.combine(date, datetime.min.time(), tzinfo=timezone.utc)
    
    for minutes in range(0, 24*60, 5):  # Check every 5 minutes
        current_time = start_time + timedelta(minutes=minutes)
        t = ts.utc(current_time)
        
        moon_pos = observer.at(t).observe(moon).apparent().altaz()
        altitude = moon_pos[0].degrees
        phase = moon_phase_percentage(t)
        
        diff = abs(altitude - target_altitude) + abs(phase - target_phase)
        if diff < best_diff:
            best_diff = diff
            best_match = current_time
    
    return best_match

# Example: Moon is 45° high and 60% illuminated
target_alt = 45
target_phase = 60
estimated_time = estimate_time(target_alt, target_phase)

print(f"Estimated time for altitude {target_alt}° and phase {target_phase}%: {estimated_time}")

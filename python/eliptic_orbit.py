import numpy as np
import matplotlib.pyplot as plt

# Constants
AU_km = 149_597_870.7  # 1 AU in km
days_per_month = 30.44

planets = {
    1: {"name": "Mercury", "a": 0.387, "e": 0.2056, "period_days": 87.97},
    2: {"name": "Venus", "a": 0.723, "e": 0.0068, "period_days": 224.7},
    3: {"name": "Earth", "a": 1.000, "e": 0.0167, "period_days": 365.25},
    4: {"name": "Mars", "a": 1.524, "e": 0.0934, "period_days": 687},
    5: {"name": "Jupiter", "a": 5.203, "e": 0.0489, "period_days": 4331},
    6: {"name": "Saturn", "a": 9.537, "e": 0.0565, "period_days": 10747},
    7: {"name": "Uranus", "a": 19.191, "e": 0.0463, "period_days": 30589},
    8: {"name": "Neptune", "a": 30.07, "e": 0.0086, "period_days": 59800},
}

def choose_planet():
    print("Pick a planet by number or name:")
    for num, data in planets.items():
        print(f"{num}: {data['name']}")
    choice = input("Your choice: ").strip()

    if choice.isdigit():
        num = int(choice)
        if num in planets:
            return planets[num]
        else:
            print("Number out of range!")
            return None
    else:
        for p in planets.values():
            if p["name"].lower() == choice.lower():
                return p
        print("Planet name not found!")
        return None

def plot_orbit_and_distance_time_centered(planet):
    a = planet["a"]  # in AU
    e = planet["e"]
    period_days = planet["period_days"]
    period_months = period_days / days_per_month

    b = a * np.sqrt(1 - e**2)
    c = a * e

    # t from -pi to pi, centered so t=0 is perihelion (closest)
    t = np.linspace(-np.pi, np.pi, 500)

    # Ellipse param, position in AU
    x = a * np.cos(t)
    y = b * np.sin(t)

    # Distance from Sun at each t (focus at (-c,0))
    dist_AU = np.sqrt((x + c)**2 + y**2)
    dist_km = dist_AU * AU_km

    # Mean distance (semi-major axis)
    mean_dist_km = a * AU_km

    # Distance variance from mean (positive and negative)
    dist_variance = dist_km - mean_dist_km

    # Map t (-pi to pi) to time (-period/2 to period/2)
    time_months = (t / (2 * np.pi)) * period_months

    import matplotlib.ticker as ticker

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

    # Orbit plot
    ax1.plot(x, y, label=f"{planet['name']} Orbit", color='dodgerblue')
    ax1.scatter(-c, 0, color='orange', s=250, label="Sun (Focus)")
    ax1.scatter(x[0], y[0], color='green', s=100, label=f"{planet['name']} at Perihelion")
    ax1.set_title(f"Orbit of {planet['name']}")
    ax1.set_xlabel("Distance (AU)")
    ax1.set_ylabel("Distance (AU)")
    ax1.legend()
    ax1.grid(True)
    ax1.axis('equal')

    # Distance variance plot vs time
    ax2.plot(time_months, dist_variance / 1e6, color='crimson')  # million km
    ax2.axhline(0, color='black', linestyle='--', linewidth=1)
    ax2.set_title(f"Distance from Sun Variance vs Time ({planet['name']})")
    ax2.set_xlabel("Time (months from Perihelion)")
    ax2.set_ylabel("Distance Variance (million km)")
    ax2.grid(True)

    # Beautify x-axis: show negative and positive months clearly
    ax2.xaxis.set_major_locator(ticker.MultipleLocator(period_months / 6))
    ax2.xaxis.set_minor_locator(ticker.MultipleLocator(period_months / 12))

    plt.tight_layout()
    plt.show()

planet = None
while not planet:
    planet = choose_planet()

plot_orbit_and_distance_time_centered(planet)

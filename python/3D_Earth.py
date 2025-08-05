import numpy as np
import plotly.graph_objects as go
from skyfield.api import load
from datetime import datetime

# Load planetary data
eph = load('de440.bsp')
ts = load.timescale()
t = ts.now()
earth = eph['Earth']

# Planets to track
planets = {
    'Mercury': eph['MERCURY BARYCENTER'],
    'Venus': eph['VENUS BARYCENTER'],
    'Mars': eph['MARS BARYCENTER'],
    'Jupiter': eph['JUPITER BARYCENTER'],
    'Saturn': eph['SATURN BARYCENTER'],
    'Uranus': eph['URANUS BARYCENTER'],
    'Neptune': eph['NEPTUNE BARYCENTER'],
}

# Earth radius for scaling
R = 1  # Unit sphere

# Get direction vectors from Earth to planets
planet_vectors = {}
for name, planet in planets.items():
    vec = earth.at(t).observe(planet).apparent().position.km
    vec = vec / np.linalg.norm(vec)  # Normalize
    planet_vectors[name] = vec

# Sphere coordinates for Earth
theta, phi = np.mgrid[0:np.pi:100j, 0:2*np.pi:100j]
x = R * np.sin(theta) * np.cos(phi)
y = R * np.sin(theta) * np.sin(phi)
z = R * np.cos(theta)

earth_surface = go.Surface(
    x=x, y=y, z=z,
    surfacecolor=np.zeros_like(z),
    colorscale='Earth',
    showscale=False,
    opacity=0.8
)

# Arrow spikes for each planet
spikes = []
for name, vec in planet_vectors.items():
    x_arrow = [0, vec[0] * 1.15]
    y_arrow = [0, vec[1] * 1.15]
    z_arrow = [0, vec[2] * 1.15]
    spikes.append(go.Scatter3d(
        x=x_arrow,
        y=y_arrow,
        z=z_arrow,
        mode='lines+text',
        line=dict(width=5),
        name=name,
        text=[None, name],
        textposition="top center"
    ))

# North–South Axis
ns_axis = go.Scatter3d(
    x=[0, 0],
    y=[0, 0],
    z=[-1.3, 1.3],
    mode='lines+text',
    line=dict(color='red', width=4, dash='dot'),
    name='N-S Axis',
    text=["South", "North"],
    textposition="top center"
)

# Prime Meridian (longitude = 0)
phi_meridian = 0
x_mer = R * np.sin(theta[:,0]) * np.cos(phi_meridian)
y_mer = R * np.sin(theta[:,0]) * np.sin(phi_meridian)
z_mer = R * np.cos(theta[:,0])
prime_meridian = go.Scatter3d(
    x=x_mer, y=y_mer, z=z_mer,
    mode='lines',
    line=dict(color='green', width=3),
    name='Prime Meridian'
)

# Equator
phi_eq = np.linspace(0, 2*np.pi, 200)
x_eq = R * np.cos(phi_eq)
y_eq = R * np.sin(phi_eq)
z_eq = np.zeros_like(phi_eq)
equator = go.Scatter3d(
    x=x_eq, y=y_eq, z=z_eq,
    mode='lines',
    line=dict(color='blue', width=3),
    name='Equator'
)

# Combine all plots
fig = go.Figure(data=[earth_surface, ns_axis, prime_meridian, equator] + spikes)

# Layout
fig.update_layout(
    title=f"Planet Directions from Earth's Core ({t.utc_datetime().strftime('%Y-%m-%d %H:%M UTC')})",
    scene=dict(
        xaxis=dict(visible=False),
        yaxis=dict(visible=False),
        zaxis=dict(visible=False),
        aspectmode='data'
    ),
    showlegend=True
)

fig.show()

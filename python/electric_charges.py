import numpy as np
import matplotlib.pyplot as plt

# --- User Input ---
print("=== Electric Charges & Coulomb Interaction Visualizer ===")
n = int(input("Enter number of charges: "))

charges = []
for i in range(n):
    q = float(input(f"Enter charge magnitude q{i+1} (Coulombs, use negative for -q): "))
    x = float(input(f"Enter x position of charge q{i+1} (meters): "))
    y = float(input(f"Enter y position of charge q{i+1} (meters): "))
    charges.append((q, np.array([x, y])))

print("\n--- Setup Complete ---")
print("You defined the following charges:")
for i, (q, pos) in enumerate(charges):
    sign = "positive" if q > 0 else "negative"
    print(f"Charge q{i+1}: {q} C ({sign}), at position {pos}")

# --- Constants ---
k = 8.988e9  # Coulomb's constant N·m²/C²

# --- Grid Setup ---
grid_size = 100
x = np.linspace(-5, 5, grid_size)
y = np.linspace(-5, 5, grid_size)
X, Y = np.meshgrid(x, y)
Ex, Ey, V = np.zeros_like(X), np.zeros_like(Y), np.zeros_like(X)

# --- Compute Field & Potential ---
for q, pos in charges:
    dx = X - pos[0]
    dy = Y - pos[1]
    r2 = dx**2 + dy**2
    r = np.sqrt(r2) + 1e-9  # avoid division by zero

    Ex += k * q * dx / r**3
    Ey += k * q * dy / r**3
    V += k * q / r

# --- Plot 1: Force Vectors ---
plt.figure(figsize=(12, 4))
plt.subplot(1, 3, 1)
plt.title("Force Vectors (Electric Field)")
plt.streamplot(X, Y, Ex, Ey, color="red", density=1, linewidth=0.7, arrowsize=1)
for q, pos in charges:
    plt.scatter(*pos, c="blue" if q > 0 else "black", s=100, marker="o")
plt.xlabel("x [m]")
plt.ylabel("y [m]")

# --- Plot 2: Field Lines ---
plt.subplot(1, 3, 2)
plt.title("Electric Field Lines")
plt.streamplot(X, Y, Ex, Ey, color=np.log(np.hypot(Ex, Ey)), cmap="inferno", density=2)
for q, pos in charges:
    plt.scatter(*pos, c="blue" if q > 0 else "black", s=100, marker="o")
plt.xlabel("x [m]")
plt.ylabel("y [m]")

# --- Plot 3: Potential Heatmap ---
plt.subplot(1, 3, 3)
plt.title("Electric Potential (Equipotential Lines)")
plt.contourf(X, Y, V, levels=50, cmap="coolwarm")
plt.colorbar(label="Potential [V]")
for q, pos in charges:
    plt.scatter(*pos, c="blue" if q > 0 else "black", s=100, marker="o")
plt.xlabel("x [m]")
plt.ylabel("y [m]")

plt.tight_layout()
plt.show()

# --- Explanations ---
print("\n--- Explanations ---")
print("1. Force Vectors: Arrows show the direction and relative strength of the electric field.")
print("   - Positive charges push field lines away.")
print("   - Negative charges pull field lines inward.")
print("2. Field Lines: Continuous curves showing where a small positive test charge would move.")
print("   - They always start on positive charges and end on negative charges.")
print("3. Potential Heatmap: Colors show potential energy per unit charge.")
print("   - Red = high potential (near positive charges).")
print("   - Blue = low potential (near negative charges).")
print("   - Equipotential lines (contours) are like 'heights' on a topographic map.")
print("\nThis combines Coulomb’s law with vector field visualization!")

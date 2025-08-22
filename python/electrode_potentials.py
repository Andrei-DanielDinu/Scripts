import matplotlib.pyplot as plt

# --- Standard Electrode Potentials (V) at 25°C ---
# Source: standard reduction potentials
electrode_potentials = {
    "Li": -3.04,
    "K": -2.93,
    "Ca": -2.87,
    "Na": -2.71,
    "Mg": -2.37,
    "Al": -1.66,
    "Zn": -0.76,
    "Fe": -0.44,
    "Ni": -0.25,
    "Sn": -0.14,
    "Pb": -0.13,
    "H": 0.00,      # Standard hydrogen electrode (SHE)
    "Cu": +0.34,
    "Ag": +0.80,
    "Hg": +0.85,
    "Au": +1.50,
    "F": +2.87
}

def calculate_cell_voltage(anode, cathode):
    """Compute cell voltage given anode and cathode metals."""
    if anode not in electrode_potentials or cathode not in electrode_potentials:
        raise ValueError("Unknown element. Please choose from the list.")

    E_anode = electrode_potentials[anode]   # oxidation happens here
    E_cathode = electrode_potentials[cathode]  # reduction happens here

    E_cell = E_cathode - E_anode
    return E_cell, E_anode, E_cathode

def main():
    print("🔋 Galvanic Cell Voltage Calculator 🔋")
    print("Available elements:", ", ".join(electrode_potentials.keys()))
    print("Note: Anode = where oxidation happens (negative electrode).")
    print("      Cathode = where reduction happens (positive electrode).\n")

    anode = input("Enter the ANODE metal (oxidation): ").capitalize()
    cathode = input("Enter the CATHODE metal (reduction): ").capitalize()

    try:
        E_cell, E_anode, E_cathode = calculate_cell_voltage(anode, cathode)
        print("\n📘 Results:")
        print(f"  Anode ({anode}):  E° = {E_anode:.2f} V  (oxidation)")
        print(f"  Cathode ({cathode}):  E° = {E_cathode:.2f} V  (reduction)")
        print(f"  → Cell Voltage: E°cell = {E_cathode:.2f} - ({E_anode:.2f}) = {E_cell:.2f} V")

        if E_cell > 0:
            print("✅ This reaction is SPONTANEOUS (battery can work).")
        else:
            print("❌ This reaction is NON-SPONTANEOUS (will not produce voltage).")

        # --- Simple visual representation ---
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.bar(["Anode", "Cathode"], [E_anode, E_cathode], color=["red", "green"])
        ax.set_ylabel("Electrode Potential (V)")
        ax.set_title(f"Galvanic Cell: {anode} | {cathode}\nCell Voltage = {E_cell:.2f} V")
        plt.show()

    except ValueError as e:
        print("Error:", e)

if __name__ == "__main__":
    main()

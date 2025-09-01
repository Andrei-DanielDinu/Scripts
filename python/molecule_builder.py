# Molecule Builder: Estimate State of Matter
# Simple Chemistry Simulator

periodic_table = {
    "H": {"name": "Hydrogen", "valence": 1, "mass": 1.008},
    "O": {"name": "Oxygen", "valence": 2, "mass": 16.00},
    "C": {"name": "Carbon", "valence": 4, "mass": 12.01},
    "N": {"name": "Nitrogen", "valence": 3, "mass": 14.01},
    "Cl": {"name": "Chlorine", "valence": 1, "mass": 35.45},
    "Na": {"name": "Sodium", "valence": 1, "mass": 22.99},
    "K": {"name": "Potassium", "valence": 1, "mass": 39.10}
    # Add more if you want
}

def molecule_builder():
    print("=== Molecule Builder ===")
    print("Pick elements and their counts. Type 'done' when finished.\n")

    molecule = {}
    total_valence = 0

    while True:
        element = input("Enter element symbol (H, O, C, N, Cl, Na, K) or 'done': ")
        if element.lower() == "done":
            break
        if element not in periodic_table:
            print("Unknown element. Try again.")
            continue

        count = int(input(f"How many atoms of {periodic_table[element]['name']}? "))
        molecule[element] = molecule.get(element, 0) + count
        total_valence += periodic_table[element]["valence"] * count

    # Saturate with hydrogen if bonds left
    if "H" not in molecule:
        molecule["H"] = 0
    if total_valence % 2 != 0:  # odd valence → add hydrogen
        molecule["H"] += 1
        total_valence += 1

    print("\nFinal Molecule:")
    formula = "".join([f"{el}{molecule[el] if molecule[el]>1 else ''}" for el in molecule])
    print(f"Formula: {formula}")

    # Estimate molar mass
    molar_mass = sum(periodic_table[el]["mass"] * count for el, count in molecule.items())
    print(f"Molar Mass ≈ {molar_mass:.2f} g/mol")

    # Rough rule for state at room temperature
    if molar_mass < 40 and "H" in molecule and len(molecule) <= 2:
        state = "Gas"
    elif molar_mass < 150:
        state = "Liquid"
    else:
        state = "Solid"

    print(f"Estimated State at Room Temperature: {state}")

molecule_builder()

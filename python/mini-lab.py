import plotly.graph_objects as go
import math
from itertools import combinations

# --- Periodic Table ---
periodic_table = {
    "H": {"name": "Hydrogen", "valence": 1, "mass": 1.008, "color": "white", "radius": 0.2},
    "O": {"name": "Oxygen", "valence": 2, "mass": 16.00, "color": "red", "radius": 0.3, "geometry":"bent", "angle":104.5},
    "C": {"name": "Carbon", "valence": 4, "mass": 12.01, "color": "black", "radius": 0.3, "geometry":"tetrahedral", "angle":109.5},
    "N": {"name": "Nitrogen", "valence": 3, "mass": 14.01, "color": "blue", "radius": 0.3, "geometry":"trigonal", "angle":107},
    "Cl": {"name": "Chlorine", "valence": 1, "mass": 35.45, "color": "green", "radius": 0.35, "geometry":"linear"}
}

# --- Build molecule ---
def build_molecule():
    molecule = {}
    print("\n--- Molecule Builder ---")
    while True:
        element = input("Element symbol (H, O, C, N, Cl) or 'done': ").capitalize()
        if element.lower() == "done":
            break
        if element not in periodic_table:
            print("❌ Unknown element.")
            continue
        try:
            count = int(input(f"Number of {element} atoms: "))
        except:
            print("❌ Invalid number.")
            continue
        molecule[element] = molecule.get(element,0) + count
    return molecule

# --- Analyze molecule ---
def analyze_molecule(molecule):
    formula = "".join([f"{el}{molecule[el] if molecule[el]>1 else ''}" for el in molecule])
    molar_mass = sum(periodic_table[el]["mass"]*count for el,count in molecule.items())
    return formula, molar_mass

# --- Generate 3D coordinates respecting valence ---
def generate_coords_valence(molecule):
    coords = []
    atoms = []
    bonds = []
    # Start with central atoms: C, N, O (highest valence)
    central_atoms=[]
    for el in sorted(molecule, key=lambda x: -periodic_table[x]["valence"]):
        for _ in range(molecule[el]):
            central_atoms.append(el)
    # Place first central atom at origin
    coords.append([0,0,0])
    atoms.append(central_atoms[0])
    valence_left = [periodic_table[central_atoms[0]]["valence"]]
    # Place remaining atoms
    for el in central_atoms[1:]:
        coords.append([len(coords)*1.5,0,0])
        atoms.append(el)
        valence_left.append(periodic_table[el]["valence"])
    # Add hydrogens or remaining atoms around central atoms
    for i,atom in enumerate(atoms):
        while valence_left[i]>0:
            # find next atom in molecule with remaining valence
            added=False
            for el2 in molecule:
                count2 = molecule[el2] - atoms.count(el2)
                if count2>0 and el2!="H":
                    atoms.append(el2)
                    x,y,z = coords[i]
                    coords.append([x+1.0, y+1.0, z])
                    bonds.append((i,len(atoms)-1))
                    valence_left[i]-=1
                    valence_left.append(periodic_table[el2]["valence"]-1)
                    added=True
                    break
            if not added:
                break
    # Add hydrogens to saturate remaining valences
    for i,vl in enumerate(valence_left):
        x,y,z = coords[i]
        for _ in range(vl):
            atoms.append("H")
            coords.append([x+0.7, y+0.7, z+0.7])
            bonds.append((i,len(atoms)-1))
    return coords, atoms, bonds

# --- Visualize molecule ---
def visualize_3d(coords, atoms, bonds):
    fig = go.Figure()
    for i,atom in enumerate(atoms):
        x,y,z = coords[i]
        fig.add_trace(go.Scatter3d(x=[x],y=[y],z=[z],mode='markers',
                                   marker=dict(size=20*periodic_table[atom]["radius"], color=periodic_table[atom]["color"]),
                                   name=atom, text=atom))
    for i,j in bonds:
        x0,y0,z0 = coords[i]
        x1,y1,z1 = coords[j]
        fig.add_trace(go.Scatter3d(x=[x0,x1],y=[y0,y1],z=[z0,z1],
                                   mode='lines', line=dict(color='gray', width=5), showlegend=False))
    fig.update_layout(scene=dict(xaxis=dict(title='X'),yaxis=dict(title='Y'),zaxis=dict(title='Z')),
                      title="3D Molecule (Valence-aware)", showlegend=True)
    fig.show()

# --- Mini Lab 7.0 ---
def mini_lab():
    molecules=[]
    print("🔬 Mini Lab 7.0 - Valence-aware 3D Molecules")
    while True:
        cmd=input("\nCommand ('new','list','exit'): ").lower()
        if cmd=="exit":
            print("Goodbye!")
            break
        elif cmd=="new":
            mol = build_molecule()
            formula,mass = analyze_molecule(mol)
            molecules.append({"formula":formula,"mass":mass})
            print(f"✅ Molecule: {formula} | Molar Mass: {mass:.2f} g/mol")
            coords, atoms, bonds = generate_coords_valence(mol)
            visualize_3d(coords, atoms, bonds)
        elif cmd=="list":
            for i,m in enumerate(molecules,1):
                print(f"{i}. {m['formula']} | {m['mass']:.2f} g/mol")
        else:
            print("Unknown command.")

if __name__=="__main__":
    mini_lab()

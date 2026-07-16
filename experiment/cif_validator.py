import os
import glob
from pymatgen.core import Structure
from chgnet.model import CHGNet
import pandas as pd

def validate_cifs(cif_directory):
    print("Loading Universal ML Potential (CHGNet)...")
    model = CHGNet.load()

    cif_files = glob.glob(os.path.join(cif_directory, "*.cif"))
    results = []

    print(f"Found {len(cif_files)} CIF files. Starting thermodynamic evaluation...")

    for cif_path in cif_files:
        filename = os.path.basename(cif_path)
        try:

            struct = Structure.from_file(cif_path)

            prediction = model.predict_structure(struct)

            energy_per_atom = prediction['e']

            max_force = max([abs(f) for force in prediction['f'] for f in force])

            results.append({
                "Filename": filename,
                "Energy_per_atom_eV": energy_per_atom,
                "Max_Force_eV_A": max_force,
                "Status": "Success"
            })
            print(f"Evaluated {filename}: Energy = {energy_per_atom:.4f} eV/atom")

        except Exception as e:
            print(f"Error evaluating {filename}: {e}")
            results.append({
                "Filename": filename,
                "Energy_per_atom_eV": None,
                "Max_Force_eV_A": None,
                "Status": f"Failed: {str(e)}"
            })

    df = pd.DataFrame(results)
    output_csv = "thermodynamic_validation_results.csv"
    df.to_csv(output_csv, index=False)
    print(f"\nValidation complete. Real physical energies saved to {output_csv}")

if __name__ == "__main__":

    CIF_DIR = "data/mof"
    validate_cifs(CIF_DIR)
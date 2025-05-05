# -*- coding: utf-8 -*-
"""
example_workflow.py: Demonstrates using the pdb_utils module
to prepare a PDB file for molecular dynamics simulation setup.

This script performs a sequence of cleaning and modification steps:
1. Creates a dummy input PDB file for demonstration.
2. Standardizes residue names (e.g., HIP->HIE).
3. Renames the assumed ligand residue to 'LIG'.
4. Makes atom names within the ligand residue unique.
5. Adds TER cards for AMBER tleap compatibility.
6. Sanitizes the file by removing CONECT records and fixing element symbols.
7. (Optional) Converts the final PDB to MOL2 using Open Babel.
"""

import os
import pdb_utils # Import the module we created

# --- Configuration ---
INPUT_PDB_NAME = "dummy_protein_ligand.pdb"
OUTPUT_DIR = "prepared_structure_example"
LIGAND_RESNAME = "LIG" # The standard 3-letter code we want for our ligand

# --- Create a Dummy Input PDB File ---
# In a real scenario, INPUT_PDB_NAME would be your actual PDB file.
dummy_pdb_content = """\
HEADER    DUMMY STRUCTURE FOR EXAMPLE
ATOM      1  N   MET A   1      10.000  10.000  10.000  1.00  0.00           N
ATOM      2  CA  MET A   1      11.000  10.500  10.000  1.00  0.00           C
ATOM      3  C   MET A   1      11.500  10.500  11.500  1.00  0.00           C
ATOM      4  O   MET A   1      11.000  10.000  12.500  1.00  0.00           O
ATOM      5  CB  MET A   1       9.000  10.000  10.000  1.00  0.00           C
ATOM      6  N   HIP B   2      12.000  11.000  11.500  1.00  0.00           N
ATOM      7  CA  HIP B   2      12.500  11.500  12.500  1.00  0.00           C
ATOM      8  HD1 HIP B   2      13.000  12.000  13.000  1.00  0.00           H
ATOM      9  C   HIP B   2      13.000  10.500  13.500  1.00  0.00           C
ATOM     10  O   HIP B   2      13.500   9.500  13.000  1.00  0.00           O
TER
HETATM   11 C1   XYZ C   1      15.000  15.000  15.000  1.00 20.00           C
HETATM   12 CL   XYZ C   1      16.000  16.000  16.000  1.00 20.00          Cl
HETATM   13 MG    MG D   2      20.000  20.000  20.000  1.00 30.00          MG
CONECT   11   12
END
"""
try:
    with open(INPUT_PDB_NAME, "w") as f:
        f.write(dummy_pdb_content)
    print(f"Created dummy input file: {INPUT_PDB_NAME}")
except IOError as e:
    print(f"Error creating dummy input file: {e}", file=sys.stderr)
    exit(1) # Exit if we can't even create the input

# --- Workflow Steps ---

print(f"\nStarting PDB preparation workflow...")
print(f"Output directory: {OUTPUT_DIR}")
os.makedirs(OUTPUT_DIR, exist_ok=True) # Create output dir if it doesn't exist

current_pdb = INPUT_PDB_NAME
step_counter = 1

# Step 1: Standardize residue names
print(f"\nStep {step_counter}: Standardizing residue names for AMBER...")
output_step1 = os.path.join(OUTPUT_DIR, f"step{step_counter}_amber_res.pdb")
current_pdb = pdb_utils.standardize_residues_for_amber(
    pdb_file=current_pdb,
    output_path=output_step1
)
if not current_pdb: # Check if the function returned None (error)
    print("Workflow aborted due to error in step 1.", file=sys.stderr)
    exit(1)
step_counter += 1

# Step 2: Rename the ligand residue (assuming 'XYZ' is the ligand here)
print(f"\nStep {step_counter}: Renaming ligand residue to '{LIGAND_RESNAME}'...")
output_step2 = os.path.join(OUTPUT_DIR, f"step{step_counter}_ligname.pdb")
current_pdb = pdb_utils.rename_hetatm_ligand(
    pdb_file=current_pdb,
    ligand_residue_name=LIGAND_RESNAME,
    # Keep the MG ion, exclude others based on default COMMON_HETATM_EXCLUDE
    exclude_list=pdb_utils.COMMON_HETATM_EXCLUDE, # Use the default list from the module
    output_path=output_step2
)
if not current_pdb:
    print("Workflow aborted due to error in step 2.", file=sys.stderr)
    exit(1)
step_counter += 1

# Step 3: Make atom names unique within the specified ligand residue
print(f"\nStep {step_counter}: Making atom names unique for ligand '{LIGAND_RESNAME}'...")
output_step3 = os.path.join(OUTPUT_DIR, f"step{step_counter}_ligunique.pdb")
current_pdb = pdb_utils.make_unique_atom_names_ligand(
    pdb_file=current_pdb,
    ligand_residue_name=LIGAND_RESNAME,
    output_path=output_step3
)
if not current_pdb:
    print("Workflow aborted due to error in step 3.", file=sys.stderr)
    exit(1)
step_counter += 1

# Step 4: Add TER cards between chains/molecule types
print(f"\nStep {step_counter}: Adding TER cards...")
output_step4 = os.path.join(OUTPUT_DIR, f"step{step_counter}_ter.pdb")
current_pdb = pdb_utils.add_ter_records_simple(
    pdb_file=current_pdb,
    output_path=output_step4
)
if not current_pdb:
    print("Workflow aborted due to error in step 4.", file=sys.stderr)
    exit(1)
step_counter += 1

# Step 5: Sanitize for tleap (remove CONECT, fix elements like CL->Cl)
print(f"\nStep {step_counter}: Sanitizing for tleap (removing CONECT, fixing elements)...")
output_step5 = os.path.join(OUTPUT_DIR, f"step{step_counter}_final_clean.pdb")
current_pdb = pdb_utils.sanitize_pdb_for_tleap(
    pdb_file=current_pdb,
    output_path=output_step5
)
if not current_pdb:
    print("Workflow aborted due to error in step 5.", file=sys.stderr)
    exit(1)
step_counter += 1


# --- Optional Step: Convert to MOL2 ---
# This requires Open Babel ('obabel') to be installed and in the system PATH.
CONVERT_TO_MOL2 = False # Set to True to enable MOL2 conversion attempt

if CONVERT_TO_MOL2:
    print(f"\nStep {step_counter}: Attempting conversion to MOL2 (requires Open Babel)...")
    output_mol2 = os.path.splitext(current_pdb)[0] + ".mol2" # Place in same dir as last PDB
    mol2_file = pdb_utils.pdb_to_mol2(
        pdb_file=current_pdb,
        output_path=output_mol2,
        add_hydrogens=False # Typically add H later during parametrization
    )
    if mol2_file:
        print(f"MOL2 conversion successful: {mol2_file}")
    else:
        print("MOL2 conversion failed or skipped.", file=sys.stderr)
    step_counter += 1


# --- Cleanup ---
# You might want to remove the intermediate files, or keep them for debugging.
# For this example, let's remove the initial dummy input file.
try:
    os.remove(INPUT_PDB_NAME)
    print(f"\nRemoved dummy input file: {INPUT_PDB_NAME}")
except OSError as e:
    print(f"\nWarning: Could not remove dummy input file {INPUT_PDB_NAME}. {e}", file=sys.stderr)

print("\nWorkflow completed.")
print(f"Final prepared PDB file: {current_pdb}")

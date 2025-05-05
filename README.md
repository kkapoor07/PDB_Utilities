# PDB Utilities for Simulation Setup

## Purpose

This repository contains a collection of Python functions (`pdb_utils.py`) designed to perform common cleaning, modification, and preparation tasks on Protein Data Bank (PDB) files. These utilities are particularly helpful for preparing structures for molecular dynamics (MD) simulations, especially using packages like AMBER, but can also be useful for general PDB standardization.

The goal is to provide simple, scriptable solutions to automate frequent preprocessing hurdles encountered with raw PDB files obtained from structural databases or initial modeling steps.

## Disclaimer & Limitations

**Important:** These utilities primarily rely on **parsing fixed-width PDB column formats** using basic string manipulation. They are **not guaranteed to work correctly on non-standard or malformed PDB files**.

*   For robust parsing and manipulation of complex or potentially non-standard PDB files, consider using dedicated cheminformatics libraries such as:
    *   [BioPython](https://biopython.org/)
    *   [MDAnalysis](https://www.mdanalysis.org/)
    *   [MDTraj](https://www.mdtraj.org/)

*   The functions modify PDB files based on common conventions (e.g., AMBER residue naming). Ensure these conventions match your specific force field and simulation requirements.

## Installation

No formal installation is required. Simply clone or download this repository. The core functionality is provided in the `pdb_utils.py` module file.

```bash
git clone https://github.com/kkapoor07/PDB_Utilities.git
cd PDB_Utilities
```

## Dependencies

*   **Python:** Version 3.6 or higher (due to f-strings).
*   **Open Babel (Optional):** The `pdb_to_mol2` function requires the Open Babel command-line tool (`obabel`) to be installed and accessible in your system's PATH. If Open Babel is not installed or found, this specific function will raise an error, but other functions will remain usable.

## Usage

Import the functions from the `pdb_utils` module into your own Python scripts.

```python
import pdb_utils
import os

# --- Example Workflow ---

# Define input/output files
input_pdb = "path/to/your/raw_structure.pdb"
base_name = os.path.splitext(input_pdb)[0]
output_dir = "prepared_structure"
os.makedirs(output_dir, exist_ok=True)

# 1. Standardize common residue names for AMBER (HIP->HIE, etc.)
amber_res_pdb = pdb_utils.standardize_residues_for_amber(
    pdb_file=input_pdb,
    output_path=os.path.join(output_dir, f"{os.path.basename(base_name)}_amber_res.pdb")
)

if amber_res_pdb:
    # 2. Rename the primary ligand HETATM to 'LIG' (excluding common ions/water)
    lig_renamed_pdb = pdb_utils.rename_hetatm_ligand(
        pdb_file=amber_res_pdb,
        ligand_residue_name="LIG",
        output_path=os.path.join(output_dir, f"{os.path.basename(base_name)}_ligname.pdb")
    )

    if lig_renamed_pdb:
        # 3. Make atom names within the 'LIG' residue unique
        lig_unique_atoms_pdb = pdb_utils.make_unique_atom_names_ligand(
            pdb_file=lig_renamed_pdb,
            ligand_residue_name="LIG",
            output_path=os.path.join(output_dir, f"{os.path.basename(base_name)}_ligunique.pdb")
        )

        if lig_unique_atoms_pdb:
            # 4. Add TER cards between chains / molecule types for tleap
            ter_added_pdb = pdb_utils.add_ter_records_simple(
                 pdb_file=lig_unique_atoms_pdb,
                 output_path=os.path.join(output_dir, f"{os.path.basename(base_name)}_ter.pdb")
            )

            if ter_added_pdb:
                # 5. Final cleaning (remove CONECT, fix element symbols) for tleap
                final_pdb = pdb_utils.sanitize_pdb_for_tleap(
                    pdb_file=ter_added_pdb,
                    output_path=os.path.join(output_dir, f"{os.path.basename(base_name)}_final.pdb")
                )

                if final_pdb:
                    print(f"Successfully prepared PDB: {final_pdb}")

                    # Optionally convert the final PDB to MOL2 (requires obabel)
                    # mol2_file = pdb_utils.pdb_to_mol2(final_pdb, add_hydrogens=False)
                    # if mol2_file:
                    #     print(f"Converted to MOL2: {mol2_file}")

```

See the `example_workflow.py` script in this repository for a runnable demonstration.

## Provided Functions (`pdb_utils.py`)

### General PDB Modifications

*   `make_unique_atom_names()`: Appends sequential numbers to all ATOM/HETATM names.
*   `rename_residue_name()`: Changes residue name for all ATOM/HETATM records.
*   `set_chain_name()`: Sets chain ID for all ATOM/HETATM records.
*   `set_residue_number()`: Sets residue number for all ATOM/HETATM records to a single value.
*   `remove_ter_records()`: Removes all TER cards.
*   `atom_to_hetatm()`: Converts all ATOM records to HETATM.

### AMBER / MD Preparation Specific

*   `standardize_residues_for_amber()`: Renames common non-standard residues (HIP->HIE, LYN->LYS, etc.) and removes problematic atoms (e.g., HD1 from HIP->HIE).
*   `delete_ligand()`: Removes HETATM records assumed to be ligands, keeping specified common ions/water/cofactors. Removes CONECT records.
*   `add_ter_records_simple()`: Adds TER cards between chains and between ATOM/HETATM blocks (useful for AMBER `tleap`).
*   `sanitize_pdb_for_tleap()`: Removes records problematic for `tleap` (CONECT, SEQRES, etc.) and corrects common element symbol formatting (e.g., CL -> Cl).
*   `rename_hetatm_ligand()`: Renames HETATMs *not* in an exclude list (ions, water) to a specified ligand name (e.g., 'LIG').
*   `make_unique_atom_names_ligand()`: Makes atom names unique *only* for HETATMs matching a specific ligand residue name.

### Specialized / Advanced

*   `prepare_for_zaff()`: Renames specific CYS, HIS, and ZN residues according to Zinc AMBER Force Field (ZAFF) conventions (Type 1 or 2). Requires careful input of binding residue numbers.
*   `replace_cysteine_with_ligand()`: Replaces atoms of a specified CYS residue with HETATM records from a separate ligand PDB file (for covalent modeling setup).

### External Tool Wrappers

*   `pdb_to_mol2()`: Converts a PDB file to MOL2 format using Open Babel (`obabel`). Requires Open Babel to be installed.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

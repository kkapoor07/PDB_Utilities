# -*- coding: utf-8 -*-
"""
pdb_utils.py: A collection of utility functions for cleaning, modifying,
and preparing PDB files, particularly for use in molecular dynamics
simulations (e.g., with AMBER).

This script relies on standard Python libraries and simple string manipulation
based on fixed PDB column formats. It may not handle non-standard PDB files
robustly. For more complex PDB parsing and manipulation, consider libraries
like BioPython, MDTraj, or MDAnalysis.

One function, pdb_to_mol2, requires the external program Open Babel ('obabel').
"""

import os
import subprocess
import sys

# --- Constants ---

# Map for standardizing common residue names for AMBER force fields
# Adjust or expand this map based on your specific needs/force field conventions
AMBER_RESIDUE_MAP = {
    'HIP': 'HIE',  # Histidine delta protonated
    'HID': 'HID',  # Histidine epsilon protonated (if needed)
    'HIE': 'HIE',  # Histidine epsilon protonated
    'LYN': 'LYS',  # N-protonated Lysine
    'ASH': 'ASP',  # Protonated Aspartic Acid
    'GLH': 'GLU',  # Protonated Glutamic Acid
    'CYX': 'CYS',  # Disulfide-bonded Cysteine (treat as CYS for setup)
    'CYM': 'CYS',  # Negatively charged Cysteine (treat as CYS for setup)
}

# Common HETATM residues to typically exclude when identifying 'the' ligand
# This list can be customized.
COMMON_HETATM_EXCLUDE = [
    # Ions
    'NA', 'Na', 'K', 'CL', 'Cl', 'CA', 'Ca', 'MG', 'Mg', 'ZN', 'Zn',
    'FE', 'Fe', 'MN', 'Mn', 'CU', 'Cu',
    # Water
    'HOH', 'WAT',
    # Common modified residues/caps
    'ACE', 'NME',
    # Common cofactors/buffers (add more as needed)
    'ADP', 'ATP', 'NAD', 'NDP', 'GOL', 'PO4', 'SO4', 'EDO',
    # DNA/RNA (can sometimes appear as HETATM)
    'DA', 'DC', 'DG', 'DT', 'A', 'C', 'G', 'U', 'RA', 'RC', 'RG', 'RU'
]

# --- Helper Function ---

def _check_file_exists(filepath, func_name):
    """Internal helper to check file existence."""
    if not os.path.exists(filepath):
        print(f"Error in {func_name}: Input file not found: {filepath}", file=sys.stderr)
        return False
    return True

def _get_output_path(input_path, output_path, suffix):
    """Internal helper to determine the output file path."""
    if output_path:
        return output_path
    else:
        base, _ = os.path.splitext(input_path)
        return f"{base}{suffix}.pdb"

# --- Core PDB Manipulation Functions ---

def make_unique_atom_names(pdb_file, output_path=None, suffix="_atomname"):
    """
    Makes atom names unique within a PDB file by appending sequential numbers.

    Strips existing digits/underscores from the base atom name first. Applies
    to both ATOM and HETATM records. Preserves PDB column formatting.

    Args:
        pdb_file (str): Path to the input PDB file.
        output_path (str, optional): Explicit path for the output file.
                                     If None, appends suffix to input name. Defaults to None.
        suffix (str, optional): Suffix used if output_path is None. Defaults to "_atomname".

    Returns:
        str: Path to the modified PDB file, or None if an error occurred.
    """
    func_name = "make_unique_atom_names"
    if not _check_file_exists(pdb_file, func_name):
        return None

    final_output_path = _get_output_path(pdb_file, output_path, suffix)
    atom_count = 1
    output_lines = []

    try:
        with open(pdb_file, 'r') as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")):
                    # Atom name is columns 13-16 (index 12-15)
                    # Get base name, removing digits/underscores
                    base_atom_name = ''.join(c for c in line[12:16].strip() if not (c.isdigit() or c == '_'))
                    # Create unique name like 'CA1', 'O10', 'FE100'
                    unique_name = f"{base_atom_name}{atom_count}"

                    # Format correctly for PDB width (4 chars, left-justified)
                    # Standard: Element right justified (cols 13-14), rest left (15-16)
                    # Simple approach: Left justify the whole unique name within 4 chars
                    formatted_name = unique_name.ljust(4)
                    if len(formatted_name) > 4: # Truncate if name + number > 4 chars
                        formatted_name = formatted_name[:4]

                    new_line = line[:12] + formatted_name + line[16:]
                    atom_count += 1
                    output_lines.append(new_line)
                else:
                    output_lines.append(line) # Keep non-atomic lines

        with open(final_output_path, "w") as f:
            f.writelines(output_lines)

        print(f"Success ({func_name}): Created file {final_output_path}")
        return final_output_path
    except IOError as e:
        print(f"Error in {func_name}: Could not read/write file. {e}", file=sys.stderr)
        return None
    except Exception as e:
        print(f"An unexpected error occurred in {func_name}: {e}", file=sys.stderr)
        return None


def rename_residue_name(pdb_file, new_residue_name, output_path=None, suffix="_resname"):
    """
    Changes the residue name for all ATOM and HETATM records in a PDB file.

    Args:
        pdb_file (str): Path to the input PDB file.
        new_residue_name (str): The new 3-letter residue name.
        output_path (str, optional): Explicit path for the output file.
                                     If None, appends suffix to input name. Defaults to None.
        suffix (str, optional): Suffix used if output_path is None. Defaults to "_resname".

    Returns:
        str: Path to the modified PDB file, or None if an error occurred.
    """
    func_name = "rename_residue_name"
    if not _check_file_exists(pdb_file, func_name):
        return None
    if len(new_residue_name) > 3:
        print(f"Warning in {func_name}: New residue name '{new_residue_name}' > 3 characters.", file=sys.stderr)
        # It might still work depending on the downstream tool, but PDB standard is 3.

    final_output_path = _get_output_path(pdb_file, output_path, suffix)
    output_lines = []
    # Residue name cols 18-20 (index 17-19), ensure 3 chars width
    formatted_res_name = new_residue_name.ljust(3)[:3] # Pad or truncate to 3

    try:
        with open(pdb_file, 'r') as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")):
                    new_line = line[:17] + formatted_res_name + line[20:]
                    output_lines.append(new_line)
                else:
                    output_lines.append(line)

        with open(final_output_path, "w") as f:
            f.writelines(output_lines)

        print(f"Success ({func_name}): Renamed residues to '{new_residue_name}', created {final_output_path}")
        return final_output_path
    except IOError as e:
        print(f"Error in {func_name}: Could not read/write file. {e}", file=sys.stderr)
        return None
    except Exception as e:
        print(f"An unexpected error occurred in {func_name}: {e}", file=sys.stderr)
        return None


def set_chain_name(pdb_file, new_chain_name, output_path=None, suffix="_chain"):
    """
    Sets the chain identifier for all ATOM and HETATM records.

    Args:
        pdb_file (str): Path to the input PDB file.
        new_chain_name (str): The new chain identifier (typically 1 character).
        output_path (str, optional): Explicit path for the output file.
                                     If None, appends suffix to input name. Defaults to None.
        suffix (str, optional): Suffix used if output_path is None. Defaults to "_chain".

    Returns:
        str: Path to the modified PDB file, or None if an error occurred.
    """
    func_name = "set_chain_name"
    if not _check_file_exists(pdb_file, func_name):
        return None
    if len(new_chain_name) != 1:
        print(f"Warning in {func_name}: Chain ID '{new_chain_name}' is not 1 character.", file=sys.stderr)
        # PDB standard is 1 char, but we allow setting it anyway.

    final_output_path = _get_output_path(pdb_file, output_path, suffix)
    output_lines = []
    # Chain ID col 22 (index 21)
    chain_char = new_chain_name[0] if new_chain_name else ' ' # Use first char or space

    try:
        with open(pdb_file, 'r') as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")):
                    new_line = line[:21] + chain_char + line[22:]
                    output_lines.append(new_line)
                else:
                    output_lines.append(line)

        with open(final_output_path, "w") as f:
            f.writelines(output_lines)

        print(f"Success ({func_name}): Set chain ID to '{chain_char}', created {final_output_path}")
        return final_output_path
    except IOError as e:
        print(f"Error in {func_name}: Could not read/write file. {e}", file=sys.stderr)
        return None
    except Exception as e:
        print(f"An unexpected error occurred in {func_name}: {e}", file=sys.stderr)
        return None


def set_residue_number(pdb_file, new_residue_number, output_path=None, suffix="_resnum"):
    """
    Sets the residue number for all ATOM and HETATM records to a single value.

    Args:
        pdb_file (str): Path to the input PDB file.
        new_residue_number (int): The new residue number.
        output_path (str, optional): Explicit path for the output file.
                                     If None, appends suffix to input name. Defaults to None.
        suffix (str, optional): Suffix used if output_path is None. Defaults to "_resnum".

    Returns:
        str: Path to the modified PDB file, or None if an error occurred.
    """
    func_name = "set_residue_number"
    if not _check_file_exists(pdb_file, func_name):
        return None

    final_output_path = _get_output_path(pdb_file, output_path, suffix)
    output_lines = []
    # Residue number cols 23-26 (index 22-25), right justified, 4 chars
    formatted_res_num = str(new_residue_number).rjust(4)
    if len(formatted_res_num) > 4:
        print(f"Warning in {func_name}: Residue number {new_residue_number} exceeds 4 digits.", file=sys.stderr)
        formatted_res_num = formatted_res_num[-4:] # Truncate to last 4

    try:
        with open(pdb_file, 'r') as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")):
                    new_line = line[:22] + formatted_res_num + line[26:]
                    output_lines.append(new_line)
                else:
                    output_lines.append(line)

        with open(final_output_path, "w") as f:
            f.writelines(output_lines)

        print(f"Success ({func_name}): Set residue number to {new_residue_number}, created {final_output_path}")
        return final_output_path
    except IOError as e:
        print(f"Error in {func_name}: Could not read/write file. {e}", file=sys.stderr)
        return None
    except Exception as e:
        print(f"An unexpected error occurred in {func_name}: {e}", file=sys.stderr)
        return None


def remove_ter_records(pdb_file, output_path=None, suffix="_noTER"):
    """
    Removes all lines starting with "TER" from a PDB file.

    Args:
        pdb_file (str): Path to the input PDB file.
        output_path (str, optional): Explicit path for the output file.
                                     If None, appends suffix to input name. Defaults to None.
        suffix (str, optional): Suffix used if output_path is None. Defaults to "_noTER".

    Returns:
        str: Path to the modified PDB file, or None if an error occurred.
    """
    func_name = "remove_ter_records"
    if not _check_file_exists(pdb_file, func_name):
        return None

    final_output_path = _get_output_path(pdb_file, output_path, suffix)
    output_lines = []

    try:
        with open(pdb_file, 'r') as f:
            for line in f:
                if not line.strip().startswith("TER"): # Use strip() for flexibility
                    output_lines.append(line)

        with open(final_output_path, "w") as f:
            f.writelines(output_lines)

        print(f"Success ({func_name}): Removed TER records, created {final_output_path}")
        return final_output_path
    except IOError as e:
        print(f"Error in {func_name}: Could not read/write file. {e}", file=sys.stderr)
        return None
    except Exception as e:
        print(f"An unexpected error occurred in {func_name}: {e}", file=sys.stderr)
        return None


def atom_to_hetatm(pdb_file, output_path=None, suffix="_HETATM"):
    """
    Changes all ATOM records to HETATM records.

    Args:
        pdb_file (str): Path to the input PDB file.
        output_path (str, optional): Explicit path for the output file.
                                     If None, appends suffix to input name. Defaults to None.
        suffix (str, optional): Suffix used if output_path is None. Defaults to "_HETATM".

    Returns:
        str: Path to the modified PDB file, or None if an error occurred.
    """
    func_name = "atom_to_hetatm"
    if not _check_file_exists(pdb_file, func_name):
        return None

    final_output_path = _get_output_path(pdb_file, output_path, suffix)
    output_lines = []

    try:
        with open(pdb_file, 'r') as f:
            for line in f:
                if line.startswith("ATOM"):
                    new_line = "HETATM" + line[6:] # Replace first 6 chars
                    output_lines.append(new_line)
                else:
                    output_lines.append(line)

        with open(final_output_path, "w") as f:
            f.writelines(output_lines)

        print(f"Success ({func_name}): Converted ATOM to HETATM, created {final_output_path}")
        return final_output_path
    except IOError as e:
        print(f"Error in {func_name}: Could not read/write file. {e}", file=sys.stderr)
        return None
    except Exception as e:
        print(f"An unexpected error occurred in {func_name}: {e}", file=sys.stderr)
        return None

# --- Functions specific to MD Setup (especially AMBER) ---

def standardize_residues_for_amber(pdb_file, output_path=None, suffix="_amber_res"):
    """
    Standardizes common residue names and removes specific atoms for AMBER compatibility.

    Uses the `AMBER_RESIDUE_MAP` constant to rename residues (e.g., HIP->HIE).
    Also removes atoms often problematic for standard AMBER force fields, such as:
      - HD1 from HIP when converting to HIE.
      - N-terminal capping group atoms like 'HC' (adjust if needed).

    Args:
        pdb_file (str): Path to the input PDB file.
        output_path (str, optional): Explicit path for the output file.
                                     If None, appends suffix to input name. Defaults to None.
        suffix (str, optional): Suffix used if output_path is None. Defaults to "_amber_res".

    Returns:
        str: Path to the modified PDB file, or None if an error occurred.
    """
    func_name = "standardize_residues_for_amber"
    if not _check_file_exists(pdb_file, func_name):
        return None

    final_output_path = _get_output_path(pdb_file, output_path, suffix)
    output_lines = []

    try:
        with open(pdb_file, 'r') as f:
            lines = f.readlines()

        for line in lines:
            if line.startswith('ATOM'):
                residue_name = line[17:20].strip()
                atom_name = line[12:16].strip()
                keep_line = True
                modified_line = line # Start with original line

                # Check if residue needs renaming
                if residue_name in AMBER_RESIDUE_MAP:
                    new_res_name = AMBER_RESIDUE_MAP[residue_name]
                    formatted_new_res = new_res_name.ljust(3)

                    # Specific atom removal rules based on common conversions
                    # Example: If converting HIP->HIE, the HD1 atom is usually removed
                    if residue_name == 'HIP' and new_res_name == 'HIE' and atom_name == 'HD1':
                        keep_line = False
                    # Add more rules here if needed, e.g., for ASH->ASP, GLH->GLU conversions

                    # Apply the residue name change if we are keeping the line
                    if keep_line:
                        modified_line = line[:17] + formatted_new_res + line[20:]

                # Check for general atom removals (independent of residue rename)
                # Example: Removing N-terminal acetyl cap hydrogen 'HC' if not parametrized
                if atom_name == 'HC':
                     print(f"Info in {func_name}: Removing atom '{atom_name}' from residue '{residue_name}' (assumed problematic cap atom).")
                     keep_line = False

                # Append the line if it wasn't marked for removal
                if keep_line:
                    output_lines.append(modified_line)

            else: # Keep non-ATOM lines (HETATM, TER, etc.) unmodified
                output_lines.append(line)

        with open(final_output_path, "w") as f:
            f.writelines(output_lines)

        print(f"Success ({func_name}): Standardized residues for AMBER, created {final_output_path}")
        return final_output_path
    except IOError as e:
        print(f"Error in {func_name}: Could not read/write file. {e}", file=sys.stderr)
        return None
    except Exception as e:
        print(f"An unexpected error occurred in {func_name}: {e}", file=sys.stderr)
        return None


def delete_ligand(pdb_file, keep_residues=COMMON_HETATM_EXCLUDE, output_path=None, suffix="_delLIG"):
    """
    Deletes HETATM records assumed to be ligands.

    Keeps HETATM residues specified in `keep_residues` (e.g., ions, water,
    common cofactors). Also removes all CONECT records, as they often refer
    to deleted ligands.

    Args:
        pdb_file (str): Path to the input PDB file.
        keep_residues (list): List of 3-letter residue names (HETATM) to preserve.
                              Defaults to `COMMON_HETATM_EXCLUDE`.
        output_path (str, optional): Explicit path for the output file.
                                     If None, appends suffix to input name. Defaults to None.
        suffix (str, optional): Suffix used if output_path is None. Defaults to "_delLIG".

    Returns:
        str: Path to the modified PDB file (protein + kept HETATMs), or None on error.
    """
    func_name = "delete_ligand"
    if not _check_file_exists(pdb_file, func_name):
        return None

    final_output_path = _get_output_path(pdb_file, output_path, suffix)
    output_lines = []
    keep_set = set(keep_residues) # Use a set for faster lookups
    ligand_deleted_count = 0

    try:
        with open(pdb_file, 'r') as f:
            lines = f.readlines()

        for line in lines:
            record_type = line[0:6].strip()

            if record_type == 'HETATM':
                residue_name = line[17:20].strip()
                # Keep if the residue name is in our explicit keep list
                if residue_name in keep_set:
                    output_lines.append(line)
                else:
                    # This HETATM is assumed to be the ligand, so skip it
                    ligand_deleted_count += 1
                    continue
            elif record_type == 'CONECT':
                # Skip CONECT records entirely when deleting ligands
                continue
            else:
                # Keep all other record types (ATOM, TER, HEADER, etc.)
                output_lines.append(line)

        with open(final_output_path, "w") as f:
            f.writelines(output_lines)

        if ligand_deleted_count > 0:
             print(f"Success ({func_name}): Deleted {ligand_deleted_count} assumed ligand HETATM atoms, created {final_output_path}")
        else:
             print(f"Warning in {func_name}: No HETATM records deleted (maybe none present or all were in keep_list?). Created {final_output_path}")
        return final_output_path

    except IOError as e:
        print(f"Error in {func_name}: Could not read/write file. {e}", file=sys.stderr)
        return None
    except Exception as e:
        print(f"An unexpected error occurred in {func_name}: {e}", file=sys.stderr)
        return None


def add_ter_records_simple(pdb_file, output_path=None, suffix="_TER"):
    """
    Adds TER records between chains and between the last ATOM and first HETATM.

    This implements a simple logic often sufficient for AMBER's tleap:
    1. Adds TER if the chain ID changes between consecutive ATOM/HETATM records.
    2. Adds TER between the last ATOM record and the first subsequent HETATM record.
    3. Adds a final TER record if the file ends with ATOM or HETATM.

    Note: Does not strictly follow all PDB conventions for TER placement but
    works well for separating molecules/chains for tleap.

    Args:
        pdb_file (str): Path to the input PDB file.
        output_path (str, optional): Explicit path for the output file.
                                     If None, appends suffix to input name. Defaults to None.
        suffix (str, optional): Suffix used if output_path is None. Defaults to "_TER".

    Returns:
        str: Path to the modified PDB file, or None if an error occurred.
    """
    func_name = "add_ter_records_simple"
    if not _check_file_exists(pdb_file, func_name):
        return None

    final_output_path = _get_output_path(pdb_file, output_path, suffix)
    output_lines = []
    last_chain_id = None
    last_record_type = None # 'ATOM', 'HETATM', or None
    ter_added_count = 0

    try:
        with open(pdb_file, 'r') as f:
            lines = f.readlines()

        for i, line in enumerate(lines):
            # Skip blank lines
            if not line.strip():
                continue

            current_record_type = None
            current_chain_id = None
            add_ter_before_this_line = False

            if line.startswith('ATOM'):
                current_record_type = 'ATOM'
                current_chain_id = line[21] # Chain ID
            elif line.startswith('HETATM'):
                current_record_type = 'HETATM'
                current_chain_id = line[21] # Chain ID
            elif line.strip().startswith('TER'):
                # If we encounter an existing TER, reset state and keep it
                output_lines.append(line)
                last_record_type = None
                last_chain_id = None
                continue # Skip TER checking logic for this line

            # --- Logic to add TER before the current line ---
            # 1. Chain break: If the last line was ATOM/HETATM and chain ID changes
            if last_record_type and current_record_type and current_chain_id != last_chain_id:
                 add_ter_before_this_line = True

            # 2. Protein to non-protein transition: If last was ATOM and current is HETATM
            if last_record_type == 'ATOM' and current_record_type == 'HETATM':
                 add_ter_before_this_line = True

            # Add the TER record if needed, ensuring no double TER
            if add_ter_before_this_line:
                # Check if the *previous actual line written* was already a TER
                if output_lines and not output_lines[-1].strip().startswith('TER'):
                    output_lines.append('TER\n')
                    ter_added_count += 1

            # Append the current line itself
            output_lines.append(line)

            # Update status for the *next* iteration
            if current_record_type: # Only update if it was ATOM/HETATM
                last_record_type = current_record_type
                last_chain_id = current_chain_id
            # If it was something else (HEADER, etc.), don't change last_* state

        # Add a final TER if the *very last* record processed was ATOM or HETATM
        if last_record_type in ['ATOM', 'HETATM']:
             # Check if the last line added wasn't already a TER
             if not output_lines[-1].strip().startswith('TER'):
                  output_lines.append('TER\n')
                  ter_added_count += 1

        with open(final_output_path, "w") as f:
            f.writelines(output_lines)

        print(f"Success ({func_name}): Added {ter_added_count} TER records, created {final_output_path}")
        return final_output_path
    except IOError as e:
        print(f"Error in {func_name}: Could not read/write file. {e}", file=sys.stderr)
        return None
    except Exception as e:
        print(f"An unexpected error occurred in {func_name}: {e}", file=sys.stderr)
        return None


def prepare_for_zaff(pdb_file, cysteine_residues, histidine_residues=None,
                     zinc_resname='ZN', zaff_type=1, output_path=None, suffix="_zaff"):
    """
    Prepares specific CYS, HIS, and ZN residues for Zinc AMBER Force Field (ZAFF).

    Renames residues based on ZAFF type (e.g., CYS->CY1, ZN->ZN1 for type 1;
    CYS->CY2, HIS->HE1, ZN->ZN2 for type 2). Also renames the backbone 'H' atom
    to 'HN' in the modified residues for AMBER compatibility.

    Requires careful identification of residue numbers involved in zinc binding.

    Args:
        pdb_file (str): Path to the input PDB file.
        cysteine_residues (list): List of integer residue numbers for CYS involved
                                  in ZN binding for this ZAFF site.
        histidine_residues (list, optional): List of integer residue numbers for HIS
                                             involved in ZN binding. Required only
                                             if `zaff_type=2`. Defaults to None.
        zinc_resname (str): Original residue name of the Zinc ion(s) to modify.
                            Defaults to 'ZN'.
        zaff_type (int): ZAFF type (1 or 2). Determines target residue names.
                         Defaults to 1.
        output_path (str, optional): Explicit path for the output file.
                                     If None, appends suffix to input name. Defaults to None.
        suffix (str, optional): Suffix used if output_path is None. Defaults to "_zaff".

    Returns:
        str: Path to the modified PDB file, or None if an error occurred.
    """
    func_name = "prepare_for_zaff"
    if not _check_file_exists(pdb_file, func_name):
        return None
    if zaff_type not in [1, 2]:
        print(f"Error in {func_name}: zaff_type must be 1 or 2.", file=sys.stderr)
        return None
    if zaff_type == 2 and not histidine_residues:
        print(f"Error in {func_name}: histidine_residues list must be provided for zaff_type=2.", file=sys.stderr)
        return None

    final_output_path = _get_output_path(pdb_file, output_path, suffix)
    output_lines = []
    cys_set = set(cysteine_residues)
    his_set = set(histidine_residues) if histidine_residues else set()

    # Determine target names based on ZAFF type
    target_cys_name = f'CY{zaff_type}'
    target_his_name = 'HE1' # Assuming HE1 (epsilon-protonated like HIE) for ZAFF type 2 binding
    target_zn_name = f'ZN{zaff_type}'

    modified_res_count = 0

    try:
        with open(pdb_file, 'r') as f:
            lines = f.readlines()

        for line in lines:
            original_line = line # Keep original for appending if no changes
            line_to_append = line # Assume we append this line unless modified or skipped

            if line.startswith('ATOM'):
                residue_name = line[17:20].strip()
                atom_name = line[12:16].strip()
                try:
                    residue_number = int(line[22:26].strip())
                except ValueError:
                    # Keep lines where residue number cannot be parsed
                    output_lines.append(original_line)
                    continue

                # --- Cysteine Modification ---
                if residue_name in ('CYS', 'CYX', 'CYM') and residue_number in cys_set:
                    line_to_append = line[:17] + target_cys_name.ljust(3) + line[20:]
                    # Rename backbone H -> HN for AMBER compatibility
                    if atom_name == 'H':
                        line_to_append = line_to_append[:12] + ' HN '.ljust(4) + line_to_append[16:]
                    modified_res_count +=1 # Count atoms modified

                # --- Histidine Modification (only for ZAFF type 2) ---
                elif zaff_type == 2 and residue_name in ('HIS', 'HIE', 'HID', 'HIP') and residue_number in his_set:
                     # Rename residue to target_his_name (HE1)
                     line_to_append = line[:17] + target_his_name.ljust(3) + line[20:]
                     # Rename backbone H -> HN
                     if atom_name == 'H':
                          line_to_append = line_to_append[:12] + ' HN '.ljust(4) + line_to_append[16:]
                     # Remove the unnecessary HD1 proton if original was HID or HIP
                     # (HE1 implies HIE-like protonation state)
                     if atom_name == 'HD1' and residue_name in ('HID', 'HIP'):
                         line_to_append = None # Mark line for skipping
                     modified_res_count +=1 # Count atoms modified

            elif line.startswith('HETATM'):
                residue_name = line[17:20].strip()
                # --- Zinc Renaming ---
                if residue_name == zinc_resname:
                    line_to_append = line[:17] + target_zn_name.ljust(3) + line[20:]
                    modified_res_count +=1

            # Append the line if it wasn't marked for skipping (e.g., HD1 removal)
            if line_to_append is not None:
                output_lines.append(line_to_append)

        with open(final_output_path, "w") as f:
            f.writelines(output_lines)

        print(f"Success ({func_name}): Prepared {modified_res_count} atoms in specified residues for ZAFF Type {zaff_type}, created {final_output_path}")
        return final_output_path

    except IOError as e:
        print(f"Error in {func_name}: Could not read/write file. {e}", file=sys.stderr)
        return None
    except Exception as e:
        print(f"An unexpected error occurred in {func_name}: {e}", file=sys.stderr)
        return None


def replace_cysteine_with_ligand(protein_file, ligand_file, cysteine_residue_number,
                                 output_path=None, suffix="_replaceCYS"):
    """
    Replaces a specified Cysteine residue in a protein PDB with HETATM records from a ligand PDB.

    This is a specialized function useful for setting up covalent docking models
    or simulations where a ligand replaces a residue side chain.

    Args:
        protein_file (str): Path to the protein PDB file.
        ligand_file (str): Path to the ligand PDB file (must contain HETATM records).
        cysteine_residue_number (int): Residue number of the CYS in the protein file
                                       to be replaced by the ligand.
        output_path (str, optional): Explicit path for the output file.
                                     If None, appends suffix to input name. Defaults to None.
        suffix (str, optional): Suffix used if output_path is None. Defaults to "_replaceCYS".

    Returns:
        str: Path to the modified PDB file, or None if an error occurred.
    """
    func_name = "replace_cysteine_with_ligand"
    if not _check_file_exists(protein_file, func_name) or not _check_file_exists(ligand_file, func_name):
        return None

    final_output_path = _get_output_path(protein_file, output_path, suffix)
    output_lines = []
    ligand_hetatm_lines = []
    ligand_inserted = False
    cys_found = False

    try:
        # 1. Read all HETATM lines from the ligand file
        with open(ligand_file, 'r') as lf:
            for line in lf:
                if line.startswith('HETATM'):
                    ligand_hetatm_lines.append(line)
        if not ligand_hetatm_lines:
             print(f"Warning in {func_name}: No HETATM records found in ligand file: {ligand_file}", file=sys.stderr)
             # Function will proceed, effectively just deleting the CYS

        # 2. Process the protein file
        with open(protein_file, 'r') as pf:
            for line in pf:
                if line.startswith('ATOM'):
                    residue_name = line[17:20].strip()
                    try:
                        residue_number = int(line[22:26].strip())
                    except ValueError:
                         output_lines.append(line) # Keep lines with unparseable resnum
                         continue

                    # Check if this is the Cysteine residue to be replaced
                    if residue_name in ('CYS', 'CYX', 'CYM') and residue_number == cysteine_residue_number:
                        cys_found = True
                        # If this is the first atom of the target CYS encountered,
                        # insert the ligand lines *before* skipping the CYS atom.
                        if not ligand_inserted:
                            output_lines.extend(ligand_hetatm_lines)
                            ligand_inserted = True
                        # Skip all ATOM lines belonging to this Cysteine residue
                        continue
                    else:
                        # Keep all other ATOM lines
                        output_lines.append(line)
                else:
                    # Keep all non-ATOM lines (TER, other HETATMs, HEADER, etc.)
                    output_lines.append(line)

        # Check if the target Cysteine was actually found
        if not cys_found:
             print(f"Warning in {func_name}: Cysteine residue {cysteine_residue_number} not found in {protein_file}.", file=sys.stderr)
        elif not ligand_inserted and ligand_hetatm_lines:
             # This case should ideally not happen if cys_found is True, but as a safeguard:
             print(f"Warning in {func_name}: CYS {cysteine_residue_number} found, but ligand insertion failed unexpectedly.", file=sys.stderr)

        # 3. Write the output file
        with open(final_output_path, "w") as f:
            f.writelines(output_lines)

        print(f"Success ({func_name}): Processed CYS {cysteine_residue_number} replacement, created {final_output_path}")
        return final_output_path

    except IOError as e:
        print(f"Error in {func_name}: Could not read/write file. {e}", file=sys.stderr)
        return None
    except Exception as e:
        print(f"An unexpected error occurred in {func_name}: {e}", file=sys.stderr)
        return None


def sanitize_pdb_for_tleap(pdb_file, output_path=None, suffix="_clean"):
    """
    Performs common cleaning steps often needed before using AMBER's tleap.

    Specifically:
    - Corrects common element symbol issues in HETATM atom names (e.g., CL -> Cl).
    - Removes PDB record types known to cause issues with tleap
      (CONECT, SEQRES, HELIX, SHEET, LINK, SSBOND).

    Args:
        pdb_file (str): Path to the input PDB file.
        output_path (str, optional): Explicit path for the output file.
                                     If None, appends suffix to input name. Defaults to None.
        suffix (str, optional): Suffix used if output_path is None. Defaults to "_clean".

    Returns:
        str: Path to the modified PDB file, or None if an error occurred.
    """
    func_name = "sanitize_pdb_for_tleap"
    if not _check_file_exists(pdb_file, func_name):
        return None

    final_output_path = _get_output_path(pdb_file, output_path, suffix)
    output_lines = []
    # Record types to remove completely
    records_to_remove = {'CONECT', 'SEQRES', 'HELIX', 'SHEET', 'LINK', 'SSBOND', 'CISPEP'}
    # Element corrections map (uppercase start -> proper case)
    element_fixes = {'CL': 'Cl', 'BR': 'Br', 'ZN': 'Zn', 'MG': 'Mg',
                     'CA': 'Ca', 'FE': 'Fe', 'MN': 'Mn', 'CU': 'Cu',
                     'NA': 'Na', 'K ': 'K '} # Handle potential space in K

    try:
        with open(pdb_file, 'r') as f:
            lines = f.readlines()

        for line in lines:
            record_type = line[0:6].strip()

            # 1. Skip records known to cause issues with tleap
            if record_type in records_to_remove:
                continue

            # 2. Fix element symbols in HETATM atom names
            #    Common issue: CL1 instead of Cl1
            elif line.startswith('HETATM'):
                atom_name_raw = line[12:16] # Keep spaces for now
                element_symbol_part = atom_name_raw[:2].strip().upper() # First two chars, uppercase

                corrected = False
                for bad_symbol, good_symbol in element_fixes.items():
                    # Check if the first 1 or 2 chars match a key in our map
                    # Handle both 'CL ' and 'CL1 ' potential matches
                    if atom_name_raw.strip().startswith(bad_symbol):
                         # Reconstruct: good symbol + rest of original name
                         rest_of_name = atom_name_raw.strip()[len(bad_symbol):]
                         new_atom_name = (good_symbol + rest_of_name).ljust(4)
                         # Check length just in case
                         if len(new_atom_name)>4: new_atom_name = new_atom_name[:4]
                         line = line[:12] + new_atom_name + line[16:]
                         corrected = True
                         break # Apply only the first matching fix

            # 3. Append the (potentially modified) line
            output_lines.append(line)

        # 4. Write the output file
        with open(final_output_path, "w") as f:
            f.writelines(output_lines)

        print(f"Success ({func_name}): Sanitized PDB for tleap, created {final_output_path}")
        return final_output_path

    except IOError as e:
        print(f"Error in {func_name}: Could not read/write file. {e}", file=sys.stderr)
        return None
    except Exception as e:
        print(f"An unexpected error occurred in {func_name}: {e}", file=sys.stderr)
        return None


def rename_hetatm_ligand(pdb_file, ligand_residue_name,
                         exclude_list=COMMON_HETATM_EXCLUDE, output_path=None, suffix="_ligname"):
    """
    Renames HETATM entries assumed to be the primary ligand(s) to a specified name.

    It identifies potential ligand atoms by checking HETATM records whose residue
    name is *not* in the `exclude_list` (which contains common ions, water, etc.).

    Args:
        pdb_file (str): Path to the input PDB file.
        ligand_residue_name (str): The desired 3-letter residue name for the ligand(s).
        exclude_list (list): List of 3-letter residue names (HETATM) to ignore.
                               Defaults to `COMMON_HETATM_EXCLUDE`.
        output_path (str, optional): Explicit path for the output file.
                                     If None, appends suffix to input name. Defaults to None.
        suffix (str, optional): Suffix used if output_path is None. Defaults to "_ligname".

    Returns:
        str: Path to the modified PDB file, or None if an error occurred.
    """
    func_name = "rename_hetatm_ligand"
    if not _check_file_exists(pdb_file, func_name):
        return None

    final_output_path = _get_output_path(pdb_file, output_path, suffix)
    output_lines = []
    exclude_set = set(exclude_list) # Use set for faster lookups
    # Ensure ligand name is 3 chars, left-justified
    formatted_lig_name = ligand_residue_name.ljust(3)[:3]
    rename_count = 0

    try:
        with open(pdb_file, 'r') as f:
            lines = f.readlines()

        for line in lines:
            if line.startswith('HETATM'):
                residue_name = line[17:20].strip()
                # Check if this HETATM residue is NOT in the exclude list
                if residue_name not in exclude_set:
                    # If not excluded, assume it's the ligand and rename it
                    line = line[:17] + formatted_lig_name + line[20:]
                    rename_count += 1

            output_lines.append(line)

        with open(final_output_path, "w") as f:
            f.writelines(output_lines)

        if rename_count > 0:
             print(f"Success ({func_name}): Renamed {rename_count} assumed ligand HETATM atoms to '{ligand_residue_name}', created {final_output_path}")
        else:
             print(f"Warning in {func_name}: No HETATM records renamed (maybe none present or all were in exclude_list?). Created {final_output_path}")
        return final_output_path

    except IOError as e:
        print(f"Error in {func_name}: Could not read/write file. {e}", file=sys.stderr)
        return None
    except Exception as e:
        print(f"An unexpected error occurred in {func_name}: {e}", file=sys.stderr)
        return None


def make_unique_atom_names_ligand(pdb_file, ligand_residue_name,
                                  output_path=None, suffix="_ligatomname"):
    """
    Makes atom names unique *only* for HETATM records matching a specific ligand residue name.

    Appends sequential numbers (1, 2, 3...) to the base atom name after removing
    any existing digits. Useful for preparing ligands for parameterization tools
    (like antechamber) that require unique atom names within the residue.

    Args:
        pdb_file (str): Path to the input PDB file.
        ligand_residue_name (str): The 3-letter residue name of the ligand whose atoms
                                   need unique names.
        output_path (str, optional): Explicit path for the output file.
                                     If None, appends suffix to input name. Defaults to None.
        suffix (str, optional): Suffix used if output_path is None. Defaults to "_ligatomname".

    Returns:
        str: Path to the modified PDB file, or None if an error occurred.
    """
    func_name = "make_unique_atom_names_ligand"
    if not _check_file_exists(pdb_file, func_name):
        return None

    final_output_path = _get_output_path(pdb_file, output_path, suffix)
    atom_count = 1 # Start numbering from 1 for this specific ligand
    output_lines = []

    try:
        with open(pdb_file, 'r') as f:
            for line in f:
                # Check if it's a HETATM *and* matches the target ligand residue name
                if line.startswith("HETATM") and line[17:20].strip() == ligand_residue_name:
                    # Extract base atom name (remove digits/underscores)
                    base_atom_name = ''.join(c for c in line[12:16].strip() if not (c.isdigit() or c == '_'))
                    # Create unique name like 'C1', 'N2', 'O10'
                    unique_name = f"{base_atom_name}{atom_count}"

                    # Format correctly for PDB width (4 chars, left-justified)
                    formatted_name = unique_name.ljust(4)
                    if len(formatted_name) > 4: # Truncate if needed
                        formatted_name = formatted_name[:4]

                    new_line = line[:12] + formatted_name + line[16:]
                    atom_count += 1
                    output_lines.append(new_line)
                else:
                    # Keep all other lines (ATOMs, other HETATMs, TER, etc.) unchanged
                    output_lines.append(line)

        # Check if any atoms were actually renamed for the specified ligand
        if atom_count == 1: # If counter is still 1, no matching atoms were found
            print(f"Warning in {func_name}: No HETATM atoms found for ligand residue '{ligand_residue_name}' in {pdb_file}.", file=sys.stderr)

        with open(final_output_path, "w") as f:
            f.writelines(output_lines)

        print(f"Success ({func_name}): Processed unique atom names for ligand '{ligand_residue_name}', created {final_output_path}")
        return final_output_path
    except IOError as e:
        print(f"Error in {func_name}: Could not read/write file. {e}", file=sys.stderr)
        return None
    except Exception as e:
        print(f"An unexpected error occurred in {func_name}: {e}", file=sys.stderr)
        return None

# --- External Tool Wrappers ---

def pdb_to_mol2(pdb_file, output_path=None, obabel_path="obabel", add_hydrogens=False):
    """
    Converts a PDB file to MOL2 format using the Open Babel command-line tool.

    Requires 'obabel' to be installed and accessible in the system PATH or
    specified via `obabel_path`.

    Args:
        pdb_file (str): Path to the input PDB file.
        output_path (str, optional): Explicit path for the output MOL2 file.
                                     If None, changes input extension to .mol2. Defaults to None.
        obabel_path (str): Path to the obabel executable (if not in PATH). Defaults to "obabel".
        add_hydrogens (bool): If True, adds the '-h' flag to the obabel command
                              to add hydrogens according to pH 7.4 model. Defaults to False.

    Returns:
        str: Path to the output MOL2 file, or None if an error occurred.
    """
    func_name = "pdb_to_mol2"
    if not _check_file_exists(pdb_file, func_name):
        return None

    # Determine output path
    if output_path:
        final_output_path = output_path
    else:
        base, _ = os.path.splitext(pdb_file)
        final_output_path = f"{base}.mol2"

    # Construct obabel command
    command = [obabel_path, pdb_file, "-O", final_output_path]
    if add_hydrogens:
        command.append("-h") # Add hydrogens appropriate for pH ~7.4

    try:
        print(f"Info ({func_name}): Running command: {' '.join(command)}")
        # Execute the command
        result = subprocess.run(command, capture_output=True, text=True, check=True, encoding='utf-8')

        print(f"Success ({func_name}): Converted PDB to MOL2, created {final_output_path}")
        # Print warnings from obabel if any
        if result.stderr:
             print(f"Open Babel stderr output:\n{result.stderr}", file=sys.stderr)
        return final_output_path

    except FileNotFoundError:
         # Specific error if obabel executable is not found
         print(f"Error in {func_name}: '{obabel_path}' command not found.", file=sys.stderr)
         print("Please ensure Open Babel is installed and its path is correct.", file=sys.stderr)
         return None
    except subprocess.CalledProcessError as e:
        # Error during obabel execution (non-zero exit code)
        print(f"Error in {func_name}: Open Babel execution failed.", file=sys.stderr)
        print(f"Command: {' '.join(e.cmd)}", file=sys.stderr)
        print(f"Return Code: {e.returncode}", file=sys.stderr)
        # Print the error message from obabel (usually in stderr)
        print(f"Stderr:\n{e.stderr}", file=sys.stderr)
        return None
    except IOError as e:
        # Error reading/writing files before/after command execution
        print(f"Error in {func_name}: File I/O error during conversion. {e}", file=sys.stderr)
        return None
    except Exception as e:
        # Catch any other unexpected errors
        print(f"An unexpected error occurred in {func_name}: {e}", file=sys.stderr)
        return None

# End of pdb_utils.py module
# No __main__ block; intended for import.

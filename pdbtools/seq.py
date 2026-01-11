#!/usr/bin/env python

# Copyright 2007, Michael J. Harms
# This program is distributed under General Public License v. 3.  See the file
# COPYING for a copy of the license.

__description__ = \
"""
pdb_seq.py

Reads sequence from pdb and returns in single amino acid or nucleotide code FASTA format.
First tries to read SEQRES entries; if this fails, uses pdb coordinates.
Supports both protein sequences (amino acids) and nucleic acid sequences (DNA/RNA).
"""
__author__ = "Michael J. Harms"
__date__ = "080123"

import os
from .data.common import *

class PdbSeqError(Exception):
    """
    Error class for this module.
    """

    pass


def pdbSeq(pdb,use_atoms=False,molecule_type="auto"):
    """
    Parse the SEQRES entries in a pdb file.  If this fails, use the ATOM
    entries.  Return dictionary of sequences keyed to chain and type of
    sequence used.
    
    molecule_type: "auto" (detect), "protein", or "nucleic"
    """

    # Try using SEQRES
    seq = [l for l in pdb if l[0:6] == "SEQRES"]
    if len(seq) != 0 and not use_atoms:
        seq_type = "SEQRES"
        chain_dict = dict([(l[11],[]) for l in seq])
        for c in list(chain_dict.keys()):
            chain_seq = [l[19:70].split() for l in seq if l[11] == c]
            for x in chain_seq:
                chain_dict[c].extend(x)
        
        # Auto-detect molecule type from SEQRES if needed
        if molecule_type == "auto" and chain_dict:
            # Sample some residues to determine type
            sample_residues = []
            for c in list(chain_dict.keys())[:3]:  # Check up to 3 chains
                sample_residues.extend(chain_dict[c][:10])  # First 10 residues
            
            # Count matches
            from .data.common import AA3_TO_AA1, NUC3_TO_NUC1
            protein_matches = sum(1 for res in sample_residues if res in AA3_TO_AA1)
            nucleic_matches = sum(1 for res in sample_residues if res in NUC3_TO_NUC1)
            
            if nucleic_matches > protein_matches:
                molecule_type = "nucleic"
            else:
                molecule_type = "protein"

    # Otherwise, use ATOM
    else:

        seq_type = "ATOM  "

        # Check to see if there are multiple models.  If there are, only look
        # at the first model.
        models = [i for i, l in enumerate(pdb) if l.startswith("MODEL")]
        if len(models) > 1:
            pdb = pdb[models[0]:models[1]]

        # Determine what type of atoms to look for based on molecule_type
        if molecule_type == "auto":
            # Auto-detect: check for protein CA atoms and nucleic C1' atoms
            has_ca = any(l[0:6] == "ATOM  " and l[13:16] == "CA " for l in pdb)
            has_c1p = any(l[0:6] == "ATOM  " and l[13:16] == "C1'" for l in pdb)
            
            if has_ca and not has_c1p:
                molecule_type = "protein"
            elif has_c1p and not has_ca:
                molecule_type = "nucleic"
            elif has_ca and has_c1p:
                # Both present - default to protein but will detect per-chain later
                molecule_type = "protein"
            else:
                # Neither found - default to protein
                molecule_type = "protein"
        
        # Grab atoms based on molecule type
        atoms = []
        if molecule_type == "nucleic":
            # For nucleic acids, use C1' (sugar carbon present in all nucleotides)
            for l in pdb:
                if l[0:6] == "ATOM  " and l[13:16] == "C1'":
                    # Check to see if this is a second conformation of the previous atom
                    if len(atoms) != 0:
                        if atoms[-1][17:26] == l[17:26]:
                            continue
                    atoms.append(l)
                elif l[0:6] == "HETATM" and l[13:16] == "C1'":
                    # Check to see if this is a second conformation of the previous atom
                    if len(atoms) != 0:
                        if atoms[-1][17:26] == l[17:26]:
                            continue
                    atoms.append(l)
        else:
            # For proteins, use CA atoms
            for l in pdb:
                if l[0:6] == "ATOM  " and l[13:16] == "CA ":
                    # Check to see if this is a second conformation of the previous atom
                    if len(atoms) != 0:
                        if atoms[-1][17:26] == l[17:26]:
                            continue
                    atoms.append(l)
                elif l[0:6] == "HETATM" and l[13:16] == "CA " and l[17:20] == "MSE":
                    # Check to see if this is a second conformation of the previous atom
                    if len(atoms) != 0:
                        if atoms[-1][17:26] == l[17:26]:
                            continue
                    atoms.append(l)

        chain_dict = dict([(l[21],[]) for l in atoms])
        for c in list(chain_dict.keys()):
            chain_dict[c] = [l[17:20].strip() for l in atoms if l[21] == c]

    return chain_dict, seq_type, molecule_type


def convertModifiedAA(chain_dict,pdb):
    """
    Convert modified amino acids to their normal counterparts.
    """

    # See if there are any non-standard amino acids in the pdb file.  If there
    # are not, return
    modres = [l for l in pdb if l[0:6] == "MODRES"]
    if len(modres) == 0:
        return chain_dict

    # Create list of modified residues
    mod_dict = dict([(l[12:15].strip(),l[24:27].strip()) for l in modres])

    # Replace all entries in chain_dict with their unmodified counterparts.
    for c in list(chain_dict.keys()):
        for i, a in enumerate(chain_dict[c]):
            if a in mod_dict:
                chain_dict[c][i] = mod_dict[a]

    return chain_dict


def pdbSeq2Fasta(pdb,pdb_id="",chain="all",use_atoms=False):
    """
    Extract sequence from pdb file and write out in FASTA format.
    Auto-detects whether chains contain protein or nucleic acid sequences.
    """

    # Grab sequences with auto-detection
    chain_dict, seq_type, detected_type = pdbSeq(pdb,use_atoms,"auto")

    # Convert modified amino acids to their natural counterparts
    chain_dict = convertModifiedAA(chain_dict,pdb)

    # Determine which chains are being written out
    if chain == "all":
        chains_to_write = list(chain_dict.keys())
        chains_to_write.sort()
    else:
        if chain in list(chain_dict.keys()):
            chains_to_write = [chain]
        else:
            err = "Chain \"%s\" not in pdb!" % chain
            raise PdbSeqError(err)

    # Detect molecule type per-chain and store original residue names for DNA/RNA detection
    # This handles mixed complexes (e.g., protein-DNA complexes)
    chain_types = {}
    chain_original_residues = {}
    for c in chains_to_write:
        # Store original residues for DNA/RNA detection
        chain_original_residues[c] = chain_dict[c][:]
        
        # Count matches in protein vs nucleic dictionaries
        protein_matches = sum(1 for res in chain_dict[c] if res in AA3_TO_AA1)
        nucleic_matches = sum(1 for res in chain_dict[c] if res in NUC3_TO_NUC1)
        
        # Detect mixed chains (both protein and nucleic in same chain - rare!)
        if protein_matches > 0 and nucleic_matches > 0:
            import sys
            print("WARNING: Chain %s appears to contain both protein and nucleic acid residues!" % c, file=sys.stderr)
            print("         Writing nucleic acids in lowercase for this chain.", file=sys.stderr)
            chain_types[c] = "mixed"
        # Use per-chain counts to determine type
        elif nucleic_matches > protein_matches:
            chain_types[c] = "nucleic"
        else:
            chain_types[c] = "protein"
    
    # Convert sequences to 1-letter format and join strings
    for c in chains_to_write:
        for aa_index, aa in enumerate(chain_dict[c]):
            if chain_types[c] == "nucleic":
                try:
                    chain_dict[c][aa_index] = NUC3_TO_NUC1[aa]
                except KeyError:
                    chain_dict[c][aa_index] = "X"
            elif chain_types[c] == "mixed":
                # Try nucleic first, then protein, use case to distinguish
                if aa in NUC3_TO_NUC1:
                    chain_dict[c][aa_index] = NUC3_TO_NUC1[aa].lower()
                elif aa in AA3_TO_AA1:
                    chain_dict[c][aa_index] = AA3_TO_AA1[aa]
                else:
                    chain_dict[c][aa_index] = "X"
            else:  # protein
                try:
                    chain_dict[c][aa_index] = AA3_TO_AA1[aa]
                except KeyError:
                    chain_dict[c][aa_index] = "X"

    out = []
    for c in chains_to_write:
        # Add molecule type to header for DNA/RNA
        mol_label = ""
        if chain_types[c] == "nucleic":
            # Determine if DNA or RNA based on original residue names
            # Check for RNA-specific bases (U and derivatives)
            has_rna = any(res in ['U', 'ADE', 'CYT', 'GUA', 'URA', 'PSU', '5MU', 'H2U', 'M2G', '2MG', '7MG', 'OMC', 'OMG', 'YYG', 'GTP', 'ATP', 'CTP', 'UTP', '1MA', 'A', 'C', 'G'] 
                          for res in chain_original_residues[c])
            # Check for DNA-specific bases (deoxy forms)
            has_dna = any(res in ['DA', 'DC', 'DG', 'DT', 'DU', '5MC', '5CM', '2MG', '6MA', 'M5M'] 
                          for res in chain_original_residues[c])
            # Look specifically for U (RNA marker) and DT (DNA marker)
            has_u = any(res in ['U', 'URA', 'PSU', '5MU', 'H2U', 'UTP'] for res in chain_original_residues[c])
            has_dt = any(res in ['DT'] for res in chain_original_residues[c])
            
            if has_u and not has_dt:
                mol_label = " RNA"
            elif has_dt and not has_u:
                mol_label = " DNA"
            elif has_u and has_dt:
                mol_label = " DNA/RNA"  # Mixed
            elif has_dna:
                mol_label = " DNA"  # DNA bases without explicit DT
            elif has_rna:
                mol_label = " RNA"  # RNA bases without explicit U
            else:
                mol_label = " DNA/RNA"  # Ambiguous
        elif chain_types[c] == "mixed":
            mol_label = " MIXED"
        
        out.append(">%s%s_%s%s" % (pdb_id,c,seq_type,mol_label))

        # Write output in lines 80 characters long
        seq_length = len(chain_dict[c])
        num_lines = seq_length // 80

        for i in range(num_lines+1):
            out.append("".join([aa for aa in chain_dict[c][80*i:80*(i+1)]]))
        out.append("".join([aa for aa in chain_dict[c][80*(i+1):]]))
        if out[-1] == "":
            out.pop(-1)


    return "\n".join(out)

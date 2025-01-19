#!/usr/bin/env python

# Copyright 2007, Michael J. Harms
# This program is distributed under General Public License v. 3.  See the file
# COPYING for a copy of the license.

"""
pdb_download.py

Download a pdb or set of pdb/mmCIF files and uncompress them using HTTPS.
"""

__author__ = "Michael J. Harms"
__date__ = "070521"
__usage__ = "pdb_download.py pdb_id or file w/ list of ids"

import os
import sys
import gzip
import requests

# Base URLs for HTTPS downloads
PDB_BASE_URL = "https://files.rcsb.org/download/"
CIF_BASE_URL = "https://files.rcsb.org/download/"

SUFFIXES = {"pdb": ".pdb", "cif": ".cif"}

def unZip(input_file, output_file):
    """
    Unzip input_file using the gzip library and write to output_file.
    """
    with gzip.open(input_file, 'rb') as f_in, open(output_file, 'wb') as f_out:
        f_out.write(f_in.read())
    os.remove(input_file)

def pdbDownload(file_list, file_format="pdb"):
    """
    Download all PDB or mmCIF files in file_list using HTTPS and unzip them.
    
    file_format can be "pdb" or "cif".
    """
    if file_format not in SUFFIXES:
        raise ValueError("Invalid file format. Use 'pdb' or 'cif'.")

    base_url = CIF_BASE_URL if file_format == "cif" else PDB_BASE_URL
    suffix = SUFFIXES[file_format]
    success = True

    for pdb_id in file_list:
        file_base_url=base_url
        file_suffix=suffix
        for ext in (".pdb", ".cif"):
            if pdb_id.endswith(ext):
                pdb_id = pdb_id[:pdb_id.index(ext)]
            file_format = ext[1:]
            file_base_url = CIF_BASE_URL if file_format == "cif" else PDB_BASE_URL
            file_suffix = SUFFIXES[file_format]
        
        pdb_id = pdb_id.strip().lower()
        if len(pdb_id) != 4:
            print(f"ERROR! {pdb_id} is not a valid PDB ID!")
            success = False
            continue

        url = f"{file_base_url}{pdb_id}{file_suffix}.gz"
        compressed_file = f"{pdb_id}{file_suffix}.gz"
        output_file = f"{pdb_id}{file_suffix}"

        try:
            print(f"Downloading {pdb_id} from {url}...")
            response = requests.get(url, stream=True)
            if response.status_code == 200:
                with open(compressed_file, 'wb') as f:
                    f.write(response.content)
                unZip(compressed_file, output_file)
                print(f"{output_file} retrieved successfully.")
            else:
                print(f"ERROR! {pdb_id} could not be retrieved! HTTP {response.status_code}")
                success = False
        except Exception as e:
            print(f"ERROR! {pdb_id} could not be retrieved! {e}")
            success = False

    return success

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: pdb_download.py pdb_id or file w/ list of ids")
        sys.exit(1)

    input_arg = sys.argv[1]

    if os.path.isfile(input_arg):
        with open(input_arg) as f:
            file_list = [line.strip() for line in f.readlines()]
    else:
        file_list = [input_arg.strip()]

    file_format = "cif" if len(sys.argv) > 2 and sys.argv[2].lower() == "cif" else "pdb"

    download_files(file_list, file_format=file_format)

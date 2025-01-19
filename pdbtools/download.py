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
BASE_URL = "https://files.rcsb.org/download"

FORMATS={"pdb", "cif"}

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
    if file_format not in FORMATS:
        raise ValueError("Invalid file format. Use 'pdb' or 'cif'.")

    suffix = file_format
    success = True

    for pdb_id in file_list:
        file_suffix=suffix
        for ext in (".pdb", ".cif"):
            if pdb_id.endswith(ext):
                pdb_id = pdb_id[:pdb_id.index(ext)]
            file_suffix = ext[1:]
        
        pdb_id = pdb_id.strip().lower()
        if len(pdb_id) != 4:
            print(f"ERROR! {pdb_id} is not a valid PDB ID!")
            success = False
            continue

        url = f"{BASE_URL}/{pdb_id}.{file_suffix}.gz"
        compressed_file = f"{pdb_id}.{file_suffix}.gz"
        output_file = f"{pdb_id}.{file_suffix}"

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

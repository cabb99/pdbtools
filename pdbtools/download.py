#!/usr/bin/env python

# Copyright 2007, Michael J. Harms
# This program is distributed under General Public License v. 3.  See the file
# COPYING for a copy of the license.

"""
pdb_download.py

Download a pdb or set of pdb files and uncompress them.
"""

__author__ = "Michael J. Harms"
__date__ = "070521"
__usage__ = "pdb_download.py pdb_id or file w/ list of ids"

import os, sys, ftplib, shutil, gzip

HOSTNAME="ftp.ebi.ac.uk"
DIRECTORY="/pub/databases/pdb/data/structures/all/pdb/"
PREFIX="pdb"
SUFFIX=".ent.gz"

# Added constants for mmCIF
CIF_DIRECTORY = "/pub/databases/pdb/data/structures/all/mmCIF/"
CIF_PREFIX = ""
CIF_SUFFIX = ".cif.gz"

def unZip(some_file, some_output):
    """
    Unzip some_file using the gzip library and write to some_output.
    """

    g = open(some_output,'wb')
    with gzip.open(some_file,'rb') as f:
        g.write(f.read()) #.decode("utf-8"))
    g.close()
    #g.writelines(f.readlines())
    #f.close()
    #g.close()

    os.remove(some_file)

def pdbDownload(file_list,hostname=HOSTNAME,directory=DIRECTORY,prefix=PREFIX,
                suffix=SUFFIX, file_format="pdb"):
    """
    Download all PDB or mmCIF files in file_list and unzip them.
    
    By default, downloads PDB files from EBI (file_format="pdb").
    To download mmCIF instead, call with file_format="cif".
    """

    # If user wants CIF, override defaults
    if file_format.lower() == "cif":
        directory = CIF_DIRECTORY
        prefix = CIF_PREFIX
        suffix = CIF_SUFFIX

    success = True

    # Log into server
    print("Connecting...")
    ftp = ftplib.FTP()
    ftp.connect(hostname)
    ftp.login()

    # Remove .pdb or .cif extensions from file_list
    for file_index, file in enumerate(file_list):
        for ext in (".pdb", ".cif"):
            if file.endswith(ext):
                file_list[file_index] = file_list[file_index][:file.index(ext)]

    # Decide final extension based on file_format
    if file_format.lower() == "cif":
        final_ext = ".cif"
    else:
        final_ext = ".pdb"

    # Download all files in file_list
    to_get = ["%s/%s%s%s" % (directory,prefix,f,suffix) for f in file_list]
    to_write = ["%s%s" % (f,suffix) for f in file_list]
    for i in range(len(to_get)):
        try:
            ftp.retrbinary("RETR %s" % to_get[i],open(to_write[i],"wb").write)
            final_name = "%s%s" % (to_write[i][:to_write[i].index(".")],final_ext)
            unZip(to_write[i],final_name)
            print("%s retrieved successfully." % final_name)
        except ftplib.error_perm:
            if os.path.exists(to_write[i]):
                os.remove(to_write[i])
            print("ERROR!  %s could not be retrieved!" % file_list[i])
            success = False

    # Log out
    ftp.quit()

    return success

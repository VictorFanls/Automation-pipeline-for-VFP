# VFP_Data_Extraction_Utils.py

import re
import os
from io import StringIO
import pandas as pd

# Extract a single numeric value by label from a line
def extract_values(line, label):
    pattern = rf"{re.escape(label)}\s*=\s*(-?\d+\.\d+)"
    match = re.search(pattern, line)
    if match:
        return float(match.group(1))
    return None

# Extract a string-formatted numeric value preserving decimal places
def extract_value_str_preserve_precision(line, label):
    pattern = rf"{re.escape(label)}\s*=\s*(-?\d+\.\d+)"
    match = re.search(pattern, line)
    if match:
        return match.group(1)  # Return as string to preserve formatting
    return None

# Extract angle of attack from a wavedrag filename
def extract_aoa_from_filename(filepath):
    filename = os.path.basename(filepath)
    match = re.search(r"M\d+Re\d+p\d+ma([-+]?\d+p\d+)", filename)
    if match:
        aoa_str = match.group(1).replace('p', '.')
        return round(float(aoa_str), 4)
    return None

# Extract angle of attack from a .forces filename
def extract_aoa_from_forces_filename(filepath):
    filename = os.path.basename(filepath)
    match = re.search(r"M\d+Re\d+p\d+ma([-+]?\d+p\d+)", filename)
    if match:
        aoa_str = match.group(1).replace('p', '.')
        return round(float(aoa_str), 4)
    return None

# Read a .forces file and a wavedrag file, extract fields
def extract_data(file_path_forces, file_path_wavedrag):
    data = {
        'Filename': os.path.basename(file_path_forces),
        'ALPHA': None,
        'MACH NO': None,
        'Cdv': None,
        'Cdi': None,
        'CL': None,
        'CDtotVFP': None,
        'CDtotIBE': None,
        'CMTOT(VFP)' : None,
        'CDW_Upper' : None,
        'CDW_Lower' : None,
        'CDW(tot)' : None,
    }

    with open(file_path_forces, 'r') as f:
        parsing = False
        for line in f:
            if not parsing:
                if "LEV=1" in line.replace(" ", ""):
                    parsing = True  # Enable parsing
                else:
                    continue  # Skip line unless LEV=1 is detected

            # Parse line regardless if it's the LEV=1 line or after
            if data['ALPHA'] is None and 'ALPHA' in line:
                data['ALPHA'] = extract_values(line, 'ALPHA')
            if data['MACH NO'] is None and 'MACH NO' in line:
                data['MACH NO'] = extract_values(line, 'MACH NO')
            if data['Cdv'] is None and 'Total viscous drag' in line:
                data['Cdv'] = extract_values(line, 'Total viscous drag')
            if data['Cdi'] is None and 'CD(vortd)' in line:
                data['Cdi'] = extract_values(line, 'CD(vortd)')
            if data['CL'] is None and 'CLTOT(IBE)' in line:
                data['CL'] = extract_values(line, 'CLTOT(IBE)')
            if data['CDtotVFP'] is None and 'CDTOT(VFP)' in line:
                data['CDtotVFP'] = extract_values(line, 'CDTOT(VFP)')
            if data['CDtotIBE'] is None and 'CDTOT(IBE)' in line:
                data['CDtotIBE'] = extract_values(line, 'CDTOT(IBE)')
            if data['CMTOT(VFP)'] is None and 'CMTOT(VFP)' in line:
                data['CMTOT(VFP)'] = extract_values(line, 'CMTOT(VFP)')

    with open(file_path_wavedrag, 'r') as f:
        count = 0
        cdw_upper_str = None
        cdw_lower_str = None
        for line in f:
            if 'Total wave drag for block is CDW(tot)' in line:
                val_str = extract_value_str_preserve_precision(line, 'Total wave drag for block is CDW(tot)')
                val = float(val_str)
                if count == 0:
                    data['CDW_Upper'] = val_str
                    cdw_upper_val = val
                elif count == 1:
                    data['CDW_Lower'] = val_str
                    cdw_lower_val = val
                count += 1
                if count == 2:
                    break

    if data['CDW_Upper'] is not None and data['CDW_Lower'] is not None:
        data['CDW(tot)'] = cdw_upper_val + cdw_lower_val

    return data
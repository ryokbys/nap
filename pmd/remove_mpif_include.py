#!/usr/bin/env python3
"""
Remove 'include mpif.h' from files that already have 'use pmdmpi'.
Since pmdmpi does 'include mpif.h', the MPI symbols are already in scope.
"""
import re
import os
import glob

PMD_DIR = os.path.dirname(os.path.abspath(__file__))

MPIF_PATTERN = re.compile(r"""^\s*include\s+['"]mpif\.h['"]\s*$""", re.IGNORECASE)

def process_file(filepath):
    with open(filepath, 'r') as f:
        lines = f.readlines()

    content = ''.join(lines)
    if 'use pmdmpi' not in content.lower():
        return False
    if not any(MPIF_PATTERN.match(l) for l in lines):
        return False

    new_lines = []
    removed = 0
    for line in lines:
        if MPIF_PATTERN.match(line):
            removed += 1
        else:
            new_lines.append(line)

    if removed > 0:
        with open(filepath, 'w') as f:
            f.writelines(new_lines)
        return True
    return False

if __name__ == '__main__':
    all_files = sorted(
        glob.glob(os.path.join(PMD_DIR, '*.F90')) +
        glob.glob(os.path.join(PMD_DIR, '*.F'))
    )
    for fp in all_files:
        if process_file(fp):
            print(f"Updated: {os.path.basename(fp)}")

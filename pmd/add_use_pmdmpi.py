#!/usr/bin/env python3
"""
Add 'use pmdmpi' to scoping units that use mpi_real_rp but lack the use statement.
"""
import re
import os
import glob

PMD_DIR = os.path.dirname(os.path.abspath(__file__))
USE_LINE = 'use pmdmpi'

def get_indent(line):
    m = re.match(r'^(\s*)', line)
    return m.group(1) if m else '  '

def already_has_use(lines, start, end):
    for i in range(start, end):
        if re.match(r'\s*use\s+pmdmpi\b', lines[i], re.IGNORECASE):
            return True
    return False

def find_insert_point(lines, start, end):
    for i in range(start, end):
        s = lines[i].strip().lower()
        if re.match(r'use\s+\w', s):
            return i
    for i in range(start, end):
        s = lines[i].strip().lower()
        if s.startswith('implicit'):
            return i
    return start

def find_scoping_units(lines):
    units = []
    i = 0
    n = len(lines)
    depth = 0

    while i < n:
        stripped = lines[i].strip().lower()
        stripped_nocomment = re.split(r'\s*!', stripped)[0].strip()

        if re.match(r'end\s*(module|program|subroutine|function|block\s*data)\b',
                    stripped_nocomment):
            if depth > 0:
                depth -= 1
            i += 1
            continue

        if stripped_nocomment == 'contains' and depth == 1:
            i += 1
            continue

        m = re.match(
            r'(module|program|subroutine|function|block\s*data)\s+(\w)',
            stripped_nocomment
        )
        if m:
            kind = m.group(1).replace(' ', '')
            if kind == 'module' and re.match(r'module\s+procedure\b', stripped_nocomment):
                i += 1
                continue
            depth += 1
            if depth == 1:
                body_start = i + 1
                while body_start < n and lines[body_start - 1].rstrip().endswith('&'):
                    body_start += 1
                units.append(body_start)
        i += 1

    return units

def process_file(filepath):
    basename = os.path.basename(filepath)
    # Skip the module that defines mpi_real_rp itself
    if basename == 'mod_pmdmpi.F90':
        return False

    with open(filepath, 'r') as f:
        lines = f.readlines()

    content = ''.join(lines)
    if 'mpi_real_rp' not in content:
        return False

    body_starts = find_scoping_units(lines)
    if not body_starts:
        return False

    insertions = []
    n = len(lines)

    for body_start in body_starts:
        end = n
        for j in range(body_start, n):
            s = lines[j].strip().lower().split('!')[0].strip()
            if re.match(r'end\s*(module|program|subroutine|function)\b', s):
                end = j
                break
            if s == 'contains':
                end = j
                break

        if already_has_use(lines, body_start, end):
            continue

        # For module files the MPI calls are inside contained subroutines
        # (after 'contains'), so check the entire file content
        if 'mpi_real_rp' not in content:
            continue

        insert_idx = find_insert_point(lines, body_start, end)
        indent = get_indent(lines[insert_idx]) if insert_idx < n else '  '
        if not indent:
            indent = '  '
        insertions.append((insert_idx, indent))

    if not insertions:
        return False

    for insert_idx, indent in sorted(insertions, reverse=True):
        lines.insert(insert_idx, f'{indent}{USE_LINE}\n')

    with open(filepath, 'w') as f:
        f.writelines(lines)
    return True

if __name__ == '__main__':
    f90_files = glob.glob(os.path.join(PMD_DIR, '*.F90'))
    f_files = glob.glob(os.path.join(PMD_DIR, '*.F'))
    all_files = sorted(f90_files + f_files)

    updated = []
    for fp in all_files:
        if process_file(fp):
            updated.append(os.path.basename(fp))

    print(f"Updated ({len(updated)}):")
    for f in updated:
        print(f"  {f}")

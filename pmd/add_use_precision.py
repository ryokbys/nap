#!/usr/bin/env python3
"""
Add 'use mod_precision' to each scoping unit in Fortran files.
"""
import re
import os
import glob

PMD_DIR = os.path.dirname(os.path.abspath(__file__))
USE_LINE = 'use mod_precision'

def get_indent(line):
    m = re.match(r'^(\s*)', line)
    return m.group(1) if m else '  '

def already_has_use(lines, start, end):
    for i in range(start, end):
        if re.match(r'\s*use\s+mod_precision\b', lines[i], re.IGNORECASE):
            return True
    return False

def find_insert_point(lines, start, end):
    """
    Return the index to insert 'use mod_precision' in lines[start:end].
    Priority: before first 'use', then before 'implicit', then after start.
    """
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
    """
    Find top-level scoping units (module, program, or standalone sub/function).
    Returns list of (start_body, end_exclusive) tuples.
    """
    units = []
    i = 0
    n = len(lines)
    depth = 0  # nesting depth of module/program/sub/function

    while i < n:
        stripped = lines[i].strip().lower()
        # Remove inline comments
        stripped_nocomment = re.split(r'\s*!', stripped)[0].strip()

        # Detect end of a scoping unit
        if re.match(r'end\s*(module|program|subroutine|function|block\s*data)\b',
                    stripped_nocomment):
            if depth > 0:
                depth -= 1
            i += 1
            continue

        # Detect 'contains' — subroutines after this are contained, skip them
        if stripped_nocomment == 'contains' and depth == 1:
            # Skip to end of the enclosing unit
            i += 1
            continue

        # Detect start of a scoping unit
        # Strip leading keywords: recursive, pure, elemental, impure
        bare = re.sub(
            r'^(recursive|pure|impure|elemental)\s+', '', stripped_nocomment
        )
        m = re.match(
            r'(module|program|subroutine|function|block\s*data)\s+(\w)',
            bare
        )
        if m:
            kind = m.group(1).replace(' ', '')
            # Skip 'module procedure'
            if kind == 'module' and re.match(r'module\s+procedure\b', stripped_nocomment):
                i += 1
                continue

            depth += 1
            if depth == 1:
                # Top-level scoping unit: body starts on next line
                body_start = i + 1
                # Skip continuation lines of the declaration
                while body_start < n and lines[body_start - 1].rstrip().endswith('&'):
                    body_start += 1
                units.append(body_start)
        i += 1

    return units

def process_file(filepath):
    if os.path.basename(filepath) == 'mod_precision.F90':
        return False

    with open(filepath, 'r') as f:
        lines = f.readlines()

    # Check real(rp) usage — if no real(rp) at all, skip
    content = ''.join(lines)
    if 'real(rp)' not in content and 'real(rp)' not in content.lower():
        return False

    body_starts = find_scoping_units(lines)
    if not body_starts:
        return False

    insertions = []  # (line_index, indent)
    n = len(lines)

    for body_start in body_starts:
        # Determine the end of this unit's declaration section
        # (before contains or end of unit)
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

        insert_idx = find_insert_point(lines, body_start, end)
        indent = get_indent(lines[insert_idx]) if insert_idx < n else '  '
        if not indent:
            indent = '  '
        insertions.append((insert_idx, indent))

    if not insertions:
        return False

    # Apply in reverse to preserve indices
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
    skipped = []
    for fp in all_files:
        if process_file(fp):
            updated.append(os.path.basename(fp))
        else:
            skipped.append(os.path.basename(fp))

    print(f"Updated ({len(updated)}):")
    for f in updated:
        print(f"  {f}")
    print(f"\nSkipped ({len(skipped)}):")
    for f in skipped:
        print(f"  {f}")

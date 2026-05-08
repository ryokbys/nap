#!/usr/bin/env python3
"""
Convert Fortran d-notation double-precision literals to real(value, rp)
in .h files. Skips comment lines.
"""
import re
import glob
import os

PMD_DIR = os.path.dirname(os.path.abspath(__file__))

# Match d-notation literals like 3.14d0, 1d-15, 2.5d+3
# Negative lookahead/behind to avoid matching partial identifiers
D_LITERAL = re.compile(
    r'(?<![a-zA-Z_])(\d+\.?\d*|\d*\.\d+)[dD]([+-]?\d+)(?![a-zA-Z_0-9])'
)

def convert_line(line):
    # Don't touch comment-only lines
    stripped = line.lstrip()
    if stripped.startswith('!') or stripped.startswith('c ') or stripped.startswith('C '):
        return line
    # Split at first inline comment to protect comment text
    if '!' in line:
        code_part, comment_part = line.split('!', 1)
        comment_part = '!' + comment_part
    else:
        code_part = line
        comment_part = ''
    new_code = D_LITERAL.sub(lambda m: f'real({m.group(0)}, rp)', code_part)
    return new_code + comment_part

def process_file(filepath):
    with open(filepath, 'r') as f:
        lines = f.readlines()
    new_lines = [convert_line(l) for l in lines]
    if new_lines != lines:
        with open(filepath, 'w') as f:
            f.writelines(new_lines)
        return True
    return False

if __name__ == '__main__':
    h_files = sorted(glob.glob(os.path.join(PMD_DIR, '*.h')))
    for fp in h_files:
        if process_file(fp):
            print(f"Updated: {os.path.basename(fp)}")
        else:
            print(f"No change: {os.path.basename(fp)}")

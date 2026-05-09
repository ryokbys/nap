#!/usr/bin/env python3
"""
Replace Fortran double-precision-specific intrinsics with generic ones.

  dsqrt -> sqrt,  dexp -> exp,  dlog -> log,  dabs -> abs,
  dsin/dcos/dtan -> sin/cos/tan,  dasin/dacos/datan/datan2,
  dsinh/dcosh/dtanh,  dnint -> nint,  dint -> aint,
  dmax1 -> max,  dmin1 -> min,  dsign -> sign,  dmod -> mod

  dble(expr) -> real(expr, rp)   (parenthesis-aware)
"""

import re
import os
import glob

PMD_DIR = os.path.dirname(os.path.abspath(__file__))

# Simple name-only replacements (longest first to avoid partial matches)
SIMPLE = [
    ('datan2',  'atan2'),
    ('dasin',   'asin'),
    ('dacos',   'acos'),
    ('datan',   'atan'),
    ('dsinh',   'sinh'),
    ('dcosh',   'cosh'),
    ('dtanh',   'tanh'),
    ('dsqrt',   'sqrt'),
    ('dexp',    'exp'),
    ('dlog10',  'log10'),
    ('dlog',    'log'),
    ('dabs',    'abs'),
    ('dsin',    'sin'),
    ('dcos',    'cos'),
    ('dtan',    'tan'),
    ('dnint',   'nint'),
    ('dint',    'aint'),
    ('dmax1',   'max'),
    ('dmin1',   'min'),
    ('dsign',   'sign'),
    ('dmod',    'mod'),
    ('dfloat',  'real'),
]

# Build one regex per simple replacement, with word-boundary on both sides
SIMPLE_PATS = [
    (re.compile(r'(?<![a-zA-Z_0-9])' + old + r'(?=[(\s])', re.IGNORECASE), new)
    for old, new in SIMPLE
]

DBLE_PAT = re.compile(r'(?<![a-zA-Z_0-9])dble\s*\(', re.IGNORECASE)


def replace_dble(code):
    """Replace dble(expr) -> real(expr, rp) throughout a code string."""
    result = []
    i = 0
    n = len(code)
    while i < n:
        m = DBLE_PAT.search(code, i)
        if m is None:
            result.append(code[i:])
            break
        result.append(code[i:m.start()])
        # Find matching closing paren
        depth = 1
        j = m.end()
        while j < n and depth > 0:
            if code[j] == '(':
                depth += 1
            elif code[j] == ')':
                depth -= 1
            j += 1
        # code[m.end():j-1] is the inner expression
        inner = code[m.end():j-1]
        result.append(f'real({inner}, rp)')
        i = j
    return ''.join(result)


def convert_line(line, is_fixed_format=False):
    if is_fixed_format and line and line[0] in 'cC*':
        return line

    code, comment = split_code_comment(line)

    # Apply simple replacements
    for pat, new in SIMPLE_PATS:
        code = pat.sub(new, code)

    # Apply dble() -> real(..., rp)
    code = replace_dble(code)

    return code + comment


def split_code_comment(line):
    in_single = False
    in_double = False
    for i, ch in enumerate(line):
        if ch == "'" and not in_double:
            in_single = not in_single
        elif ch == '"' and not in_single:
            in_double = not in_double
        elif ch == '!' and not in_single and not in_double:
            return line[:i], line[i:]
    return line, ''


def process_file(filepath):
    ext = os.path.splitext(filepath)[1].lower()
    is_fixed = ext == '.f'
    with open(filepath, 'r') as f:
        lines = f.readlines()
    new_lines = [convert_line(l, is_fixed_format=is_fixed) for l in lines]
    if new_lines != lines:
        with open(filepath, 'w') as f:
            f.writelines(new_lines)
        return True
    return False


if __name__ == '__main__':
    targets = sorted(
        glob.glob(os.path.join(PMD_DIR, '*.F90')) +
        glob.glob(os.path.join(PMD_DIR, '*.F'))
    )
    for fp in targets:
        if process_file(fp):
            print(f"Updated: {os.path.basename(fp)}")

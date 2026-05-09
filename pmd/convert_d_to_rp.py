#!/usr/bin/env python3
"""
Convert Fortran d-notation double-precision literals to _rp kind suffix.

  0.5d0        -> 0.5_rp
  1d0          -> 1.0_rp   (integer mantissa: add .0)
  1.d0         -> 1._rp
  1d-15        -> 1e-15_rp
  3.14d+2      -> 3.14e+2_rp

After conversion, also clean up redundant wrappers in .h files:
  real(0.5_rp, rp)  -> 0.5_rp
"""

import re
import os
import glob

PMD_DIR = os.path.dirname(os.path.abspath(__file__))

# Match d-notation literals; avoid matching inside identifiers
D_LIT = re.compile(
    r'(?<![a-zA-Z_])'          # not preceded by letter/underscore
    r'(\d+\.?\d*|\d*\.\d+)'    # mantissa: 1, 1.5, .5, 1.
    r'[dD]'                     # d or D
    r'([+-]?\d+)'               # exponent: 0, -15, +3
    r'(?![a-zA-Z_0-9])',        # not followed by alphanumeric
    re.IGNORECASE
)

# Match real(X_rp, rp) pattern to unwrap redundant conversions
REAL_WRAP = re.compile(
    r'\breal\(([^,()]+_rp)\s*,\s*rp\)',
    re.IGNORECASE
)


def d_to_rp(m):
    mantissa = m.group(1)
    exp = m.group(2)
    exp_val = int(exp)

    if exp_val == 0:
        if '.' in mantissa:
            return f'{mantissa}_rp'
        else:
            return f'{mantissa}.0_rp'
    else:
        return f'{mantissa}e{exp}_rp'


def convert_line(line, is_fixed_format=False):
    """Convert d-notation in code portion of one line."""

    # Fixed-format: column 1 is c/C/* => full-line comment
    if is_fixed_format and line and line[0] in 'cC*':
        return line

    # Split off inline comment (careful about ! inside strings)
    # Simple heuristic: find first ! not inside a string
    code, comment = split_code_comment(line)
    new_code = convert_code_part(code)
    return new_code + comment


def split_code_comment(line):
    """Return (code_part, comment_part) splitting on first ! outside strings."""
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


def convert_code_part(code):
    return D_LIT.sub(d_to_rp, code)


def process_file(filepath, cleanup_wrap=False):
    ext = os.path.splitext(filepath)[1].lower()
    is_fixed = ext == '.f'

    with open(filepath, 'r') as f:
        lines = f.readlines()

    new_lines = [convert_line(l, is_fixed_format=is_fixed) for l in lines]

    if cleanup_wrap:
        new_lines = [REAL_WRAP.sub(r'\1', l) for l in new_lines]

    if new_lines != lines:
        with open(filepath, 'w') as f:
            f.writelines(new_lines)
        return True
    return False


if __name__ == '__main__':
    targets = (
        glob.glob(os.path.join(PMD_DIR, '*.F90')) +
        glob.glob(os.path.join(PMD_DIR, '*.F')) +
        glob.glob(os.path.join(PMD_DIR, '*.h'))
    )

    for fp in sorted(targets):
        is_h = fp.endswith('.h')
        if process_file(fp, cleanup_wrap=is_h):
            print(f"Updated: {os.path.basename(fp)}")

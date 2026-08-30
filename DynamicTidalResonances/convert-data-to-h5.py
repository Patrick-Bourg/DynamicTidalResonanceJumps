"""
Routines to convert Patrick's tidal resonance data into h5py format
"""

import numpy as np
import h5py
import re
import os

def parse_mathematica_data(data_str: str) -> dict:
    """
    Parses a Mathematica-formatted dataset into a dictionary.

    Generated with Claude Sonnet 4.6

    Returns:
        dict mapping (a, p, e, x) -> (dE, dL, dQ) or "No resonance"
    """

    def convert_mathematica_number(s: str):
        """Convert Mathematica number strings like '1.*^-15' to Python floats."""
        s = s.strip()
        # Mathematica scientific notation: 1.*^-15 -> 1e-15
        s = re.sub(r'\*\^', 'e', s)
        try:
            return float(s)
        except ValueError:
            return s  # return as-is if not a number

    def parse_value(s: str):
        """Parse a single value: number, quoted string, or bare word."""
        s = s.strip()
        if s.startswith('"') and s.endswith('"'):
            return s[1:-1]  # strip quotes
        return convert_mathematica_number(s)

    def parse_list(s: str) -> list:
        """Parse a flat Mathematica list like {a, b, c} into a Python list."""
        s = s.strip()
        assert s.startswith('{') and s.endswith('}'), f"Expected list, got: {s}"
        inner = s[1:-1]
        items = split_top_level(inner)
        return [parse_value(item) for item in items]

    def split_top_level(s: str) -> list[str]:
        """Split a string by commas at the top nesting level only."""
        parts = []
        depth = 0
        current = []
        in_string = False
        for ch in s:
            if ch == '"':
                in_string = not in_string
            if not in_string:
                if ch == '{':
                    depth += 1
                elif ch == '}':
                    depth -= 1
                elif ch == ',' and depth == 0:
                    parts.append(''.join(current).strip())
                    current = []
                    continue
            current.append(ch)
        if current:
            parts.append(''.join(current).strip())
        return parts

    def parse_entry(s: str):
        """
        Parse one full entry: {{n,k,m}, {a,p,e,x}, {dE,dL,dQ}} or "No resonance" 
        as the third element.
        Returns ((a, p, e, x), (dE, dL, dQ) or "No resonance")
        """
        s = s.strip()
        assert s.startswith('{') and s.endswith('}')
        inner = s[1:-1]
        parts = split_top_level(inner)
        assert len(parts) == 3, f"Expected 3 parts in entry, got {len(parts)}: {parts}"

        # parts[1] is {a, p, e, x}
        apex = tuple(parse_list(parts[1]))

        # parts[2] is either {dE, dL, dQ} or "No resonance"
        third = parts[2].strip()
        if third.startswith('"'):
            value = third[1:-1]  # strip quotes -> "No resonance"
        else:
            value = tuple(parse_list(third))

        return apex, value

    # Strip outer braces and split into top-level entries
    data_str = data_str.strip()
    assert data_str.startswith('{') and data_str.endswith('}')
    inner = data_str[1:-1]
    entries = split_top_level(inner)

    result = {}
    for i, entry in enumerate(entries):
        entry = entry.strip()
        if not entry.startswith('{'):
            continue  # skip the header row (which starts with {{"n","k","m"},...})
        # Skip header: it contains string keys like "n", "k", "m"
        # Detect header by checking if first sub-list has quoted strings
        first_brace_content = entry[1:entry.index('}')]
        if '"' in first_brace_content:
            continue  # this is the header row

        apex, value = parse_entry(entry)
        result[apex] = value

    return result

def parse_mathematica_complex(s: str) -> complex:
    s = s.strip()

    # Replace Mathematica scientific notation: 1.*^-15 -> 1e-15
    s = re.sub(r'\*\^', 'e', s)

    # Replace Mathematica imaginary unit: *I -> j, bare I -> 1j
    s = re.sub(r'\*I\b', 'j', s)
    s = re.sub(r'\bI\b', '1j', s)

    # Normalise trailing decimal points: '0.' -> '0.0', '1.' -> '1.0'
    s = re.sub(r'(\d+\.)(?=[^0-9]|$)', r'\g<1>0', s)

    # Rewrite the operator between real and imaginary parts only,
    # using a lookbehind for a digit to avoid touching exponent signs.
    # Collapse surrounding whitespace but don't introduce '+' before '-'
    s = re.sub(r'(?<=\d)\s*\+\s*(?=[\d-])', '+', s)
    s = re.sub(r'(?<=\d)\s*-\s*(?=[\d])', '-', s)

    s = s.lstrip('+')

    return complex(s)

def convert_file_to_dict(filename: str) -> tuple[tuple[int, int, int], dict]:

    # Extract (n, k, m) from filename    
    parts = filename.split('/')[-1].split('_')
    n = int(parts[0][1:])
    k = int(parts[1][1:])
    m = int(parts[2][1:-2])

    with open(filename, 'r') as f:
        data = f.readlines()
        data = [line.rstrip() for line in data]
        out = ''
        for line in data:
            out+= line
    
    data_dict = parse_mathematica_data(out[57:]) # skip preamble

    avals = []
    pvals = []
    evals = []
    xvals = []
    for a, p, e, x in data_dict.keys():
        avals.append(a)
        pvals.append(p if p != "No resonance" else np.nan)
        evals.append(e if e != 1e-15 else 0.0)
        xvals.append(x)

    Ejump = []
    Ljump = []
    Qjump = []
    for jump in data_dict.values():
        if isinstance(jump, str) and jump == "No resonance":
            Ejump.append(np.nan)
            Ljump.append(np.nan)
            Qjump.append(np.nan)
        else:
            dE, dL, dQ = jump
            Ejump.append(parse_mathematica_complex(dE))
            Ljump.append(parse_mathematica_complex(dL))
            Qjump.append(parse_mathematica_complex(dQ))

    return (m, k, n), {
        'a': np.array(avals),
        'p': np.array(pvals),
        'e': np.array(evals),
        'x': np.array(xvals),
        'dE': np.array(Ejump),
        'dL': np.array(Ljump),
        'dQ': np.array(Qjump),
    }

if __name__ == "__main__":
    data_dir = "./Tidal_Resonance_Jumps/"
    resonance_data_dicts = {}
    
    for filename in os.listdir(data_dir):
        if ".m" not in filename:
            continue
        print(f"Processing {filename}...")
        mkn, data_dict = convert_file_to_dict(os.path.join(data_dir, filename))
        h5_filename = filename.replace('.txt', '.h5')
        resonance_data_dicts[mkn] = data_dict
    
    # construct grid from first one
    first_grid = resonance_data_dicts[list(resonance_data_dicts.keys())[0]]
    av = list(set(first_grid['a']))
    ev = list(set(first_grid['e']))
    xv = list(set(first_grid['x']))

    reshape_shape = (len(av), len(ev), len(xv))

    f = h5py.File("tidal_resonance_data.h5", 'w')
    for mkn, data_dict in resonance_data_dicts.items():
        group_name = f"m{mkn[0]}_k{mkn[1]}_n{mkn[2]}"
        group = f.create_group(group_name)
        for key, value in data_dict.items():
            group.create_dataset(key, data=value)
        
        grid = group.create_group("grid")
        grid.create_dataset("a", data=np.array(sorted(av)))
        grid.create_dataset("e", data=np.array(sorted(ev)))
        grid.create_dataset("x", data=np.array(sorted(xv)))
        grid.create_dataset("p", data=data_dict['p'].reshape(reshape_shape))
        grid.create_dataset("dE", data=data_dict['dE'].reshape(reshape_shape))
        grid.create_dataset("dL", data=data_dict['dL'].reshape(reshape_shape))
        grid.create_dataset("dQ", data=data_dict['dQ'].reshape(reshape_shape))

    f.close()
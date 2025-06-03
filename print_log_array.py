#!/usr/bin/env python3
import struct
import sys

def read_floats(filename):
    data = open(filename, 'rb').read()
    # ensure length is a multiple of 4
    if len(data) % 4 != 0:
        print(f"Warning: file length {len(data)} not a multiple of 4 bytes")
    floats = []
    for i in range(0, len(data), 4):
        word = data[i:i+4]
        # little-endian 32-bit IEEE‑754 float
        val = struct.unpack('<f', word)[0]
        floats.append(val)
    return floats

def main():
    # grab all args
    args = sys.argv[1:]
    # if first arg is '-1', drop it (means “skip bitrev output”)
    if args and args[0] == '-1':
        args = args[1:]
    # drop any V/NV flags
    args = [a for a in args if a.upper() not in ('V', 'NV')]

    # default files if none specified
    if not args:
        files = ['bitrev_output.hex', 'fft_output.hex']
    else:
        files = args

    for fn in files:
        try:
            vals = read_floats(fn)
        except FileNotFoundError:
            print(f"File not found: {fn}")
            continue

        print(f"\n{fn}:")
        # print as rows of complex pairs
        for i in range(0, len(vals), 2):
            r, im = vals[i], vals[i+1]
            print(f"  [{i//2:2d}] = {r:.6f} + {im:.6f}j")

if __name__ == "__main__":
    main()
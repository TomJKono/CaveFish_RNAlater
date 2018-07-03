#!/usr/bin/env python
"""Get the significantly differentially expressed genes from the CSV summary
table that compares the RNAlater to LN storage. Takes three arguments:
    1) Summary CSV
    2) P-value threshold for significance
    3) Up, down, or both
Writes *transcript IDs* to stdout."""

import sys


def main(pvalues, thresh, direction):
    """Main function."""
    # Set the appropriate arguments
    try:
        t = float(thresh)
        assert direction in ['Up', 'Down', 'Both']
    except ValueError:
        sys.stderr.write('Please supply a floating point number for the threshold.\n')
        exit(2)
    except AssertionError:
        sys.stderr.write('Please supply one of "Up," "Down," or "Both" for the direction. Case-sensitive.\n')
    # Iterate through the file and get the genes that pass the threshold
    with open(pvalues, 'r') as f:
        for index, line in enumerate(f):
            if index == 0:
                continue
            else:
                tmp = line.strip().split(',')
                try:
                    p = float(tmp[6])
                except ValueError:
                    continue
                if p < t:
                    txid = tmp[0].replace('G', 'T') + '.1'
                    # Check the direction!
                    if direction == 'Both':
                        print(txid)
                    else:
                        if tmp[7] == direction:
                            print(txid)
                        else:
                            continue
    return


if len(sys.argv) != 4:
    print("""Get the significantly differentially expressed genes from the CSV summary
table that compares the RNAlater to LN storage. Takes three arguments:
    1) Summary CSV
    2) P-value threshold for significance
    3) Up, down, or both
Writes *transcript IDs* to stdout.""")
    exit(1)
else:
    main(sys.argv[1], sys.argv[2], sys.argv[3])

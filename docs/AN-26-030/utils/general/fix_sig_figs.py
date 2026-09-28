#!/usr/bin/python3

# A simple script to round numbers to two sig figs (or another specified precision) in lists
# of numbers (for example, systematics tables). 
# We have to use the decimal package in order to get numbers to round the way we want them to.
# Flags:
#  -n <num>: use to change the number of sig figs to round to (default: 2)
#  -e: by default, this script will change numbers in exponential notation to LaTeX (e.g., 1.2\ten{6} instead
#      of 1.2E+6). If you want to keep the exponential notation, use this flag.

import sys
import re
import argparse
import decimal
import io

parser = argparse.ArgumentParser(description="Fix number of significant figures in input text.")
# make this optional so that we can use test mode without a file, but then we need to check below
parser.add_argument("input_file", help="Name of input file", nargs='?')
parser.add_argument("-n", "--sig-figs", help="Number of significant figures to round to (default=2)", type=int, default=2)
parser.add_argument("-e", "--exponential", help="Use exponential notation instead of LaTeX formatting for large or small numbers (i.e. 1.2E+6 instead of 1.2\\ten{6})", action="store_true")
parser.add_argument("-t", "--test-mode", help="Test mode: show some sample input and output (no input file is used)", action="store_true")
args = parser.parse_args()
if not args.test_mode and not args.input_file:
    parser.print_help()
    sys.exit(0)

# The decimal module uses banker's rounding by default, but for our purposes it probably makes more sense to
# use the traditional "always round 5 up" rule.
decimal.getcontext().rounding = decimal.ROUND_HALF_UP
decimal.getcontext().prec = args.sig_figs

# Data for testing
test_data = """1   12   123   1234   12345   123456   1234567   12345678
0.1   0.12   0.123   0.0123   0.00123   0.000123   0.0000123   0.00000123
1234   1234.567   123   123.456   12   12.3456   1   1.2   1.234
1250   1250.123   125   125.123   15   15.1234   1   1.5   1.512
0.1   0.12      0.12345   0.01   0.012   0.0123   0.001   0.0012   0.00123456
0.5   0.25      0.12512   0.05   0.015   0.0175   0.002   0.0025   0.00125123
.1      .12      .12345   .01   .012   .0123      .001   .0012   .00123456
1.23e-7   1.23e-5   1.23e-3   1.23E0   1.23E2   1.23E+7
+1234   +1234.567   +123   +123.456   +12   +12.3456   +1   +1.2   +1.234
-1234   -1234.567   -123   -123.456   -12   -12.3456   -1   -1.2   -1.234
"""

if args.test_mode:
    infile = io.StringIO(test_data)
else:
    infile = open(args.input_file)

for line in infile:
    if args.test_mode:
        print("Input: ", line+"Output: ", end='')

    # This RE probably misses some wacky corner cases but should get the vast majority of numbers.
    # Note that we have to consider the 123(.456) and .123 cases separately, since we have can't have a
    # number with digits neither before the decimal point nor after!
    for n in re.findall(r"""[-+]?                      # optionally begin with a sign
                            (?:[0-9]+(?:\.[0-9]+)?     # one or more digits, optionally followed by a decimal
                                                       # point and more digits
                            |\.[0-9]+)                 # or, just the decimal point and some digits
                            (?:[eE][-+]?[0-9]+)?       # optionally followed by E and an exponent
                         """,
                        line, flags=re.VERBOSE):
        # the unary plus here triggers the rounding of the value to the number of sig figs above
        new_n = +decimal.Decimal(n)

        # For numbers with more digits in the integer part than sig figs, the decimal module uses
        # scientific notation. This is more technically correct, but looks kind of awkward -- I think
        # most people would prefer 120 to 1.2E+02. So in this case, strip the exponent. For numbers
        # greater than a million or so, though, we can keep using scientific notation.
        if abs(new_n) >= 10**args.sig_figs and abs(new_n) < 1e6:
            new_n = new_n.quantize(decimal.Decimal(1), context=decimal.Context(prec=6))

        new_n = str(new_n)
        # Convert exponential notation to LaTeX, since that seems like the more likely use case,
        # unless we've been asked not to. Note that minus signs get saved, but plus signs do not.
        if not args.exponential:
            new_n = re.sub(r"E\+?(-?[0-9]+)", r"\\ten{\1}", new_n)
        line = line.replace(n, new_n, 1)

    print(line.rstrip())

infile.close()

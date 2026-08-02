#!/usr/bin/env python
"""The C++ header and the Python data file must not drift.

They are the same table for two languages. Nothing compiled or tested the
header before, so it could diverge from the list the oracle reads with nothing
failing -- and the index IS the feature position.
"""
import os
import re
import sys

root = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")

data = [ln.strip() for ln in open(f"{root}/data/peptdeep_mod_elements.txt")
        if ln.strip() and not ln.startswith("#")]
header = open(f"{root}/include/odia/PeptDeepElements.h").read()
body = header[header.index("PEPTDEEP_MOD_ELEMENTS{"):header.index("};")]
from_header = re.findall(r'"([^"]+)"', body)

failures = []
if len(data) != 109:
    failures.append(f"data file has {len(data)} elements, must be 109")
if from_header != data:
    for i, (a, b) in enumerate(zip(from_header, data)):
        if a != b:
            failures.append(f"header and data diverge at index {i}: {a!r} vs {b!r}")
            break
    else:
        failures.append(f"header has {len(from_header)} elements, data has {len(data)}")
if len(set(data)) != len(data):
    failures.append("duplicate symbols would make index lookup ambiguous")
# 'No' must survive as a string; YAML 1.1 reads a bare No as false.
if "No" not in data:
    failures.append("element 'No' is missing -- a YAML round-trip may have made it False")

for f in failures:
    print(f"  FAIL {f}")
print(f"element table: {len(failures)} failures ({len(data)} elements)")
sys.exit(1 if failures else 0)

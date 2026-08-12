# Convert an ODIA/DIA-NN library TSV to OpenSWATH's transition TSV convention.
#
# OpenSwathAssayGenerator reads by column NAME and ours are DIA-NN's, so it
# read Precursor.Id where it wanted a number and died with
#   Could not convert string 'AAAAAAAAAAAAAAAAGATC(Carbamidomethyl)LER2'
#     to a double value
#
# awk, not the csv module: the library is 59.8M transitions / 8.6 GB, and this
# runs in a couple of minutes rather than half an hour.
#
# The PeptideSequence column must be the BARE sequence. Deleting non-alpha
# characters is NOT how you get it -- `C(Carbamidomethyl)LER` would become
# `CCarbamidomethylLER`, silently turning the modification name into residues.
# Parenthesised and bracketed groups have to go first, as whole groups.

BEGIN { FS = OFS = "\t" }

NR == 1 {
  for (i = 1; i <= NF; ++i) { c[$i] = i }
  # Fail loudly on a layout we do not recognise rather than emitting a file
  # full of empty columns that OpenSWATH would read as zeros.
  split("Precursor.Id Modified.Sequence Precursor.Charge Decoy RT " \
        "Precursor.Mz Product.Mz Relative.Intensity Fragment.Type " \
        "Fragment.Charge Fragment.Series.Number Protein.Group", need, " ")
  for (i in need) {
    if (!(need[i] in c)) { print "missing column: " need[i] > "/dev/stderr"; exit 2 }
  }
  print "PrecursorMz", "ProductMz", "LibraryIntensity", "NormalizedRetentionTime",
        "ProteinId", "PeptideSequence", "ModifiedPeptideSequence", "PrecursorCharge",
        "ProductCharge", "FragmentType", "FragmentSeriesNumber",
        "transition_group_id", "transition_name", "Decoy", "PrecursorIonMobility"
  next
}

{
  mod = $(c["Modified.Sequence"])
  bare = mod
  gsub(/\([^)]*\)/, "", bare)          # (Carbamidomethyl), (UniMod:4)
  gsub(/\[[^]]*\]/, "", bare)          # [+57.0215]
  gsub(/[^A-Za-z]/, "", bare)          # any stray punctuation that is left
  print $(c["Precursor.Mz"]), $(c["Product.Mz"]), $(c["Relative.Intensity"]),
        $(c["RT"]), $(c["Protein.Group"]), bare, mod,
        $(c["Precursor.Charge"]), $(c["Fragment.Charge"]),
        $(c["Fragment.Type"]), $(c["Fragment.Series.Number"]),
        $(c["Precursor.Id"]), $(c["Precursor.Id"]) "_" NR, $(c["Decoy"]), k0()
}

# 1/K0 from the library's predicted CCS, for OpenSWATH's PASEF mode.
#
# Without this column OpenSwathWorkflow refuses a diaPASEF run outright:
#   "Transition 0 does not have a valid IM value, this must be set to use PASEF
#    mode (auto-detected or via -pasef flag)"
# -- and it auto-detects PASEF from the .d, so the column is not optional and
# omitting it is not a way to run without mobility.
#
# The library carries CCS in square angstroms, not the 1/K0 an instrument
# reports. ODIA refuses this conversion ON PURPOSE -- CCS and 1/K0 are related
# through the Mason-Schamp equation via the drift gas and the instrument
# calibration, so converting inside the library would write assumptions into a
# field consumers read as measured. Those assumptions are made HERE instead,
# in the benchmark harness, where they can be stated:
#
#   drift gas nitrogen (28.0134 Da), T = 305 K
#   1/K0 = CCS * sqrt(mu * T) / (18509 * z),  mu = M*28.0134/(M + 28.0134)
#
# Sanity check on a real row: CCS 527.49, z 2, m/z 464.568 gives 1/K0 = 1.30,
# which is in the normal band for a doubly-charged tryptic peptide.
#
# A wrong T or gas shifts every value by a common factor, which a per-run
# mobility calibration would absorb; a wrong FORM would not. Returns 0 when CCS
# is absent, which is what a non-PASEF run gets and does not use.
function k0(   ccs, z, m, mu) {
  ccs = ("CCS" in c) ? $(c["CCS"]) + 0 : 0
  z   = $(c["Precursor.Charge"]) + 0
  if (ccs <= 0 || z <= 0) { return 0 }
  m  = $(c["Precursor.Mz"]) * z - z * 1.00728      # neutral monoisotopic mass
  mu = (m * 28.0134) / (m + 28.0134)
  return ccs * sqrt(mu * 305.0) / (18509.0 * z)
}

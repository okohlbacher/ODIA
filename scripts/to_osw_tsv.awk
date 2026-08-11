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
        "transition_group_id", "transition_name", "Decoy"
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
        $(c["Precursor.Id"]), $(c["Precursor.Id"]) "_" NR, $(c["Decoy"])
}

#!/bin/bash

### Input start here ###
T=300         # Temperature in K
B=2000         # Number of bins
col=2         # Column in colvar_**.log containing ξ
input_file="../kappa-win.dat"  # Window config file
### Input finish here ###

shuffle_inplace() {
  local f="$1"
  if command -v shuf >/dev/null 2>&1; then
    shuf -o "$f" "$f"
  elif command -v gshuf >/dev/null 2>&1; then
    gshuf -o "$f" "$f"
  else
    perl -MList::Util=shuffle -e 'print shuffle(<STDIN>)' < "$f" > "${f}.tmp" && mv "${f}.tmp" "$f"
  fi
}


### Calculation starts here ###

# Process each line in kappa-win.dat, skipping header
mkdir -p tmp

tail -n +2 "$input_file" | while read index kappa center; do
    # Extract the ξ column from plumed output
    awk -v var="$col" 'NR > 1 {print $var}' ../plumed/colvar_${index}.log > colvar_${index}.log

    # Apply harmonic bias potential formatting
    awk -v var1="$kappa" -v var2="$center" '{ $2 = $2 var1; $3 = $3 var2; print }' colvar_${index}.log > window_${index}.log

    # Copy and shuffle data
    cat window_${index}.log > windowd_${index}.log
    shuffle_inplace "windowd_${index}.log"

done

# Clean up intermediate files
rm colvar_*
rm window_*

# Run umbrella integration
mv windowd_* tmp/
cd tmp
umbrella_integration.x -ui -d . -T $T -min -2.95 -max 0.10 -n $B -u 'kcal' -r -1 -seg 200 -v 1 > ../umbrella_integration.log
mv *xy ../
cd ..
rm -r tmp


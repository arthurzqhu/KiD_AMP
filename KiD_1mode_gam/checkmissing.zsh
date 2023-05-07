#!/bin/zsh

# Iterate through all directories matching the pattern fullmic_conftest_Na*_*w*
for top_dir in output/conftest_diffsp/fullmic_conftest_Na*w*; do
   [[ $top_dir =~ 'fullmic_conftest_Na([0-9]+)w([0-9]+)' ]]
   e="$match[1]"
   f="$match[2]"

   # Iterate through all subdirectories matching the pattern pmomxy*-*-*
   for sub_dir in $top_dir/BIN_TAU/pmomxy*-*; do
      [[ $sub_dir =~ 'pmomxy([0-9]+)-([0-9]+)' ]]
      a="$match[1]"
      b="$match[2]"

      # Iterate through all spcr directories
      for spcr_dir in $sub_dir/spcr*-*; do
         [[ $spcr_dir =~ 'spcr([0-9]{2})-([0-9]{2})' ]]
         c="$match[1]"
         d="$match[2]"

         # Check if the spcr directory is empty
         if [ -z "$(ls -A "$spcr_dir")" ]; then
            # Run the fullmic_conftest.zsh script with the extracted values
            ./fullmic_conftest.zsh $a $b $c $d $e $f
            # echo $a $b $c $d $e $f
         fi
      done
   done
done


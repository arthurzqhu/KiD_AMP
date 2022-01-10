#!/bin/zsh
  
case=$1

# make sure configured right in defines.inc
if grep -xq "^#CASE=$1" src/defines.inc; then
  sed -i "s/^CASE=/#CASE=/gm" src/defines.inc # replaces anything uncommented
  sed -i "s/^#CASE=$1/CASE=$1/gm" src/defines.inc
  make clean # probably need to clean in this case
fi

# make sure the right exe is compiled
if [ ! -f "./bin/KiD_$1.exe" ]; then
   make CASE=$1 all
fi

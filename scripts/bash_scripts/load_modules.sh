#!/bin/bash

##########################
# LOAD MODULES WITH LMOD #
##########################
containsElement () {
  local e
  for e in "${@:2}"; do [[ "$e" == "$1" ]] && return 0; done
  return 1
}
if [ $LMOD -eq 1 ]
then
    while read LINE
    do
        LINE_ARRAY=($LINE)
        PROGRAM_NAME=${LINE_ARRAY[0]}
        MODULE_NAME=${LINE_ARRAY[1]}
        if containsElement "$PROGRAM_NAME" "${REQUIRED_MODULES[@]}"
        then
            echo "Loading module $MODULE_NAME"
            eval "module load $MODULE_NAME"
        fi
    done < $resource_dir/OPTIONS_lmod 
fi

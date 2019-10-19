#!/bin/sh
################################################################################################
# Emerge code tools - consistent_trees_to_emerge.sh                                            #
# Call with <base_file_name> and <box_divisions>, e.g. for tree_0_0_0.dat to tree_1_1_1.dat:   #
# consistent_trees_to_emerge.sh tree 2                                                         #
################################################################################################

FILENAME="$1"
BOX_DIVISIONS=$2

A=0
for ((i=0; i<$BOX_DIVISIONS; i++)) do
  for ((j=0; j<$BOX_DIVISIONS; j++)) do
    for ((k=0; k<$BOX_DIVISIONS; k++)) do

      ./convert_CT_to_emerge tree_${i}_${j}_${k}.dat ${FILENAME}.${A}

      ((A++))
    done
  done
done

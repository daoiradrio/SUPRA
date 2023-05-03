#!/bin/bash

not_converged=0
total=0

for folder in optdir*
do
	if [ -f "$folder/not.uffconverged" ]
	then
		not_converged=$((not_converged+1))
	fi
	total=$((total+1))
done

echo $not_converged
echo $total

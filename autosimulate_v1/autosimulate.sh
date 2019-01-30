#!/bin/bash
python reset.py
echo "Press enter to continue"
read
var=1
time=0
while [[ 1 -eq 1 ]]; do
	# ((time++))
	# echo $time
	python set_input.py
	bin/ib
	# python collect_output.py
	var=$(python collect_output.py 2>&1)
	# echo "var=$var"
	if [ "$var" -gt "1" ]; then
		echo "var=$var"
		echo break!
		break;
		exit 0;
	fi
done



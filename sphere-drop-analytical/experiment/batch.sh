#! /bin/bash
for x in 4 3 2 1 0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.1 
do
for y in {1..5}
do
    ./sphere-drop-analytical 10 10 10 10 $x $x &
done
wait
done

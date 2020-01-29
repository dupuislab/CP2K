#!/bin/sh

cd gopt_state0
var1=$(awk '/Total FORCE_EVAL / {print $9}' output.out |tail -1)
echo "0.0", $var1
cd ..

for i in {1..5}
do
cd point_0.${i}
var1=$(awk '/Total FORCE_EVAL / {print $9}' output.out |tail -1)
echo 0.$i, $var1
cd ..
done

cd point_0.5_from1
var1=$(awk '/Total FORCE_EVAL / {print $9}' output.out |tail -1)
echo "0.5_from1", $var1
cd ..

for i in {6..9}
do
cd point_0.${i}
var1=$(awk '/Total FORCE_EVAL / {print $9}' output.out |tail -1)
echo 0.$i, $var1
cd ..
done

cd gopt_state1
var1=$(awk '/Total FORCE_EVAL / {print $9}' output.out |tail -1)
echo "1.0", $var1
cd ..


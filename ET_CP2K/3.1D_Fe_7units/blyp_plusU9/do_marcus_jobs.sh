#!/bin/sh

cd gopt_state0
echo "gopt_state0"
awk '/GEOMETRY OPTIMIZATION COMPLETED/' output.out
tail -79 p_fe_-pos-1.xyz |tee opt0.xyz
cd ../

cd gopt_state1
echo "gopt_state1"
awk '/GEOMETRY OPTIMIZATION COMPLETED/' output.out
tail -79 p_fe_-pos-1.xyz |tee opt1.xyz
cd ../

cp gopt_state0/opt0.xyz .
cp gopt_state1/opt1.xyz .

for i in {1..9}
do
mkdir point_0.$i
cp input.inp ./point_0.$i
cp slurmCP2K ./point_0.$i

done

mkdir point_0.5_from1
cp ./point_0.5/* ./point_0.5_from1/

for i in {1..9}
do
paste opt0.xyz opt1.xyz | awk -v point=0.$i '{if (FNR == 1) print $1;else if (FNR == 2) print "intermediate point",point;else printf "%s %.10f %.10f %.10f\n", $1, ($2*(1-point)+$6*point), ($3*(1-point)+$7*point), ($4*(1-point)+$8*point);}'|tee intermediate_point.xyz
cp intermediate_point.xyz ./point_0.$i/
done

cp ./point_0.5/intermediate_point.xyz ./point_0.5_from1/


for i in {1..5}
do
cp ./gopt_state0/p_fe_-RESTART.wfn ./point_0.$i
done

for i in {6..9}
do
cp ./gopt_state1/p_fe_-RESTART.wfn ./point_0.$i
done

cp ./gopt_state1/p_fe_-RESTART.wfn ./point_0.5_from1/

for i in {1..9}
do
cd ./point_0.$i
sed -i -e "s/job-name=/job-name=pont0.$i/g" slurmCP2K
scabatch slurmCP2K
cd ../
done

cd ./point_0.5_from1
sed -i -e "s/job-name=/job-name=pont0.5_from1/g" slurmCP2K
scabatch slurmCP2K
cd ../

cp point_0.1/* ./state1_at_point0_geo 
cp slurmCP2K ./state1_at_point0_geo
cp gopt_state0/opt0.xyz ./state1_at_point0_geo/ 
cp gopt_state1/p_fe_-RESTART.wfn ./state1_at_point0_geo/ 
cd state1_at_point0_geo/ 
sed -e -i "s/intermediate_point.xyz/opt0.xyz/g" input.inp
rm output.out 
sed -i -e "s/job-name=/job-name=state1_at_geo0/g" slurmCP2K
scabatch slurmCP2K
cd ../

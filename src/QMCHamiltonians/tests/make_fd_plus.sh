#!/usr/bin/bash

build_dir=/home/rclay/Software/qmcpack_worktree/zvzb_pseudo_stress/build_complex/src/QMCHamiltonians/tests

DELT="1e-7"
for i in 0 1 2; do 
  for j in `seq 0 1 $i`; do
    echo $i $j
    sed -e "s/int si=1,sj=2/int si=$j,sj=$i/g" -e "s/sdelta=1e-7/sdelta=$DELT/g" test_ion_derivs_REF.cpp > test_ion_derivs.cpp
    cd $build_dir
    make
    ./test_hamiltonian_force --turn-on-printout > e_${j}_${i}.p.d$DELT.txt
    cd -
  done
done
for i in 0 1 2; do 
  for j in `seq 0 1 $i`; do
    echo $i $j
    sed -e "s/int si=1,sj=2/int si=$j,sj=$i/g" -e "s/sdelta=1e-7/sdelta=-$DELT/g" test_ion_derivs_REF.cpp > test_ion_derivs.cpp
    cd $build_dir
    make
    ./test_hamiltonian_force --turn-on-printout > e_${j}_${i}.m.d$DELT.txt
    cd -
  done
done
#sed "s/int si=1,sj=2/int si="3",sj=0/g" test_ion_derivs.cpp
#int si=1,sj=2

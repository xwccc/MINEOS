#!/bin/bash
l=001
n=2
model='premANIC.card'
dir='tmp'
typ='sph'
./read_mineos << EOF
./${dir}/${model}_${typ}_0000${l}.fun
mode_fun.asc
$n $l
EOF

#!/bin/sh

truchas=/home/nnc/Telluride/truchas/remove-old-ht/bin/t-linux.x86_64.nag.serial.opt-2.7.dev
parser=/home/nnc/Telluride/truchas/remove-old-ht/tools/scripts/TBrookParse.py

$truchas runA.inp
$truchas runB.inp
python $parser < restart.mac
$truchas -r:restart.bin runAB.inp

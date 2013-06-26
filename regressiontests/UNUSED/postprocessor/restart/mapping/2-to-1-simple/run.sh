#!/bin/sh

truchas=/home/nnc/Telluride/truchas/remove-old-ht/bin/t-linux.x86_64.nag.serial.opt-2.7.dev
parser=/home/nnc/Telluride/truchas/remove-old-ht/tools/scripts/TBrookParse.py

$truchas boxA.inp
$truchas boxB.inp
python $parser < restart.mac
$truchas -r:restart.bin boxAB.inp

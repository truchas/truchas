#!/bin/sh

truchas -f test-init.inp
./add-solid-void.py test-init_output/test-init.h5 restart.bin
truchas -f -r:restart.bin test.inp

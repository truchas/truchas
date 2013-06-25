#!/bin/sh

### TEST-2A ####################################################################

command="addgaps mesh2.exo gap2a.exo 10 20"
echo $command
$command << EOF > /dev/null
1
2
EOF
cat addgaps.log

### TEST-2B ####################################################################

command="addgaps mesh2.exo gap2b.exo 10 22 40"
echo $command
$command << EOF > /dev/null
1
2
4
EOF
cat addgaps.log

### TEST-2C ####################################################################

command="addgaps mesh2.exo gap2c.exo 10 22"
echo $command
$command << EOF > /dev/null
1
2
EOF
cat addgaps.log

### TEST-2D ####################################################################

command="addgaps mesh2.exo gap2d.exo 10 20 30"
echo $command
$command << EOF > /dev/null
1
2
3
EOF
cat addgaps.log

### TEST-2E ####################################################################

command="addgaps mesh2.exo gap2e.exo 10 20 30 40"
echo $command
$command << EOF > /dev/null
1
2
3
4
EOF
cat addgaps.log

### TEST-2F ####################################################################

command="addgaps mesh2.exo gap2f.exo 10"
echo $command
$command << EOF > /dev/null
1
EOF
cat addgaps.log

### TEST-2G ####################################################################

command="addgaps -s mesh2.exo gap2g.exo 20"
echo $command
$command << EOF > /dev/null
2
EOF
cat addgaps.log

### TEST-2H ####################################################################

command="addgaps mesh2.exo junk.exo 20 22"
echo $command
$command

command="addgaps mesh2.exo junk.exo 21 22"
echo $command
$command


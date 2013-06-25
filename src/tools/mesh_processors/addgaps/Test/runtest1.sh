#!/bin/sh

### TEST-1A ####################################################################

addgaps mesh1.exo gap1a.exo 10 20 << EOF
1
2
EOF

cubit << EOF
import mesh geometry "gap1a.exo"
from 12.846238 -13.405936  12.452249
at    3.031168   2.849848   2.733799
up   -0.256370   0.377243   0.889923
draw block 1 color green
draw block 2 color red add
hardcopy "gap1a.tiff" tiff
quit
EOF
convert gap1a.tiff gap1a.png
rm gap1a.tiff

### TEST-1B ####################################################################

addgaps mesh1.exo gap1b.exo 10 22 40 << EOF
1
2
4
EOF

cubit << EOF
import mesh geometry "gap1b.exo"
from 12.846238 -13.405936  12.452249
at    3.031168   2.849848   2.733799
up   -0.256370   0.377243   0.889923
draw block 1 color green
draw block 2 color red add
draw block 4 color orange add
hardcopy "gap1b.tiff" tiff
quit
EOF
convert gap1b.tiff gap1b.png
rm gap1b.tiff

### TEST-1C ####################################################################

addgaps mesh1.exo gap1c.exo 10 22 << EOF
1
2
EOF

cubit << EOF
import mesh geometry "gap1c.exo"
from 12.846238 -13.405936  12.452249
at    3.031168   2.849848   2.733799
up   -0.256370   0.377243   0.889923
draw block 1 color green
draw block 2 color red add
hardcopy "gap1c.tiff" tiff
quit
EOF
convert gap1c.tiff gap1c.png
rm gap1c.tiff

### TEST-1D ####################################################################

addgaps mesh1.exo gap1d.exo 10 20 30 << EOF
1
2
3
EOF

cubit << EOF
import mesh geometry "gap1d.exo"
from 12.846238 -13.405936  12.452249
at    3.031168   2.849848   2.733799
up   -0.256370   0.377243   0.889923
draw block 1 color green
draw block 2 color red add
draw block 3 color blue add
hardcopy "gap1d.tiff" tiff
quit
EOF
convert gap1d.tiff gap1d.png
rm gap1d.tiff

### TEST-1E ####################################################################

addgaps mesh1.exo gap1e.exo 10 20 30 40 << EOF
1
2
3
4
EOF

cubit << EOF
import mesh geometry "gap1e.exo"
from 12.846238 -13.405936  12.452249
at    3.031168   2.849848   2.733799
up   -0.256370   0.377243   0.889923
draw block 1 color green
draw block 2 color red add
draw block 3 color blue add
draw block 4 color orange add
hardcopy "gap1e.tiff" tiff
quit
EOF
convert gap1e.tiff gap1e.png
rm gap1e.tiff

### TEST-1F ####################################################################

addgaps mesh1.exo gap1f.exo 10 << EOF
1
EOF

cubit << EOF
import mesh geometry "gap1f.exo"
from 12.846238 -13.405936  12.452249
at    3.031168   2.849848   2.733799
up   -0.256370   0.377243   0.889923
draw block 1 color green
hardcopy "gap1f.tiff" tiff
quit
EOF
convert gap1f.tiff gap1f.png
rm gap1f.tiff

### TEST-1G ####################################################################

addgaps -s mesh1.exo gap1g.exo 20 << EOF
2
EOF

cubit << EOF
import mesh geometry "gap1g.exo"
from 12.846238 -13.405936  12.452249
at    3.031168   2.849848   2.733799
up   -0.256370   0.377243   0.889923
draw block 2 color red
hardcopy "gap1g.tiff" tiff
quit
EOF
convert gap1g.tiff gap1g.png
rm gap1g.tiff


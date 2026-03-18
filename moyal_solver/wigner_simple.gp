#!/usr/bin/gnuplot

# Simple Wigner function GIF generator
# Modify the file list below to match your actual files

set terminal gif animate delay 10 size 800,600
set output 'wigner_evolution.gif'

# Plot settings
set pm3d map
set view map
set palette rgbformulae 33,13,10

# Axis limits as requested
set xrange [-5:5]
set yrange [-20:20]

set xlabel 'x'
set ylabel 'p' 
set colorbox

# File pattern - modify this list based on your actual files wigner_00000000000
do for [step=0:3000:10] {
    filename = sprintf("build/output/wigner_%011d.dat", step)

#    filename = "build/output/".word(files,i)
    set title sprintf('Wigner Function - Frame %d', i)
    splot filename using 1:2:3 with image
}

unset output
print "Done! Check wigner_evolution.gif"

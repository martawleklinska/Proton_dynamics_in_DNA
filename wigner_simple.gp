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
set yrange [-15:15]

set xlabel 'x'
set ylabel 'p' 
set colorbox

# Automatic loop through all files from wigner_00000000000.dat to wigner_00000001500.dat
# Files increment by 5: 0, 5, 10, 15, ..., 1500

# Animation loop using automatic file generation
do for [step=0:1500:5] {
    filename = sprintf("build/output/wigner_%011d.dat", step)
    set title sprintf('Wigner Function - Step %d (t = %.4f)', step, step*0.0001)
    
    # Check if we can plot (gnuplot will skip if file doesn't exist)
    splot filename using 1:2:3 with pm3d notitle
}

unset output
print "Done! Check wigner_evolution.gif"

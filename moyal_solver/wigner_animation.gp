#!/usr/bin/gnuplot

# Gnuplot script to create animated GIF of Wigner function evolution
# Data format: x p wigner_value irrelevant_column
# Files: wigner_00000000000.dat, wigner_00000000005.dat, etc.

reset

# Set terminal for GIF animation
set terminal gif animate delay 10 size 800,600
set output 'wigner_evolution.gif'

# Set plot parameters
set xlabel 'x' font ",14"
set ylabel 'p' font ",14" 
set title 'Wigner Function Evolution' font ",16"

# Set axis limits
set xrange [-5:5]
set yrange [-15:15]

# Set up color palette and contour settings
set palette rgbformulae 33,13,10  # Blue-white-red colormap
set pm3d map                      # Use surface mapping
set view map                      # Top-down view

# Optional: Add colorbox/colorbar
set colorbox
set cbrange [*:*]  # Auto-scale colorbar

# Optional: Set contour lines
# set contour base
# set cntrparam levels auto 10

# Function to extract frame number from filename for title
frame_number(filename) = int(substr(filename, 8, 11))

# Define data directory
datadir = "build/output/"

# Get list of files (you may need to adjust this based on your system)
# For this script, we'll manually define the file pattern

# Animation loop - adjust the range based on your actual files
do for [i=0:100:5] {  # Assuming files increment by 5, adjust as needed
    # Format filename with zero padding
    filename = sprintf("%swigner_%011d.dat", datadir, i)
    
    # Check if file exists (optional, may cause issues on some systems)
    # You can comment this out if it causes problems
    
    # Set title with frame information
    set title sprintf('Wigner Function Evolution - Step %d', i) font ",16"
    
    # Plot the data
    # Using columns 1 (x), 2 (p), and 3 (wigner value)
    splot filename using 1:2:3 with pm3d notitle
}

# Close output
unset output

print "Animation complete! Output saved as wigner_evolution.gif"

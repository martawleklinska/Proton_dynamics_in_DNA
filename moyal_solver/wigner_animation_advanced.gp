#!/usr/bin/gnuplot

# Advanced Wigner function animation script
# Handles variable file patterns and includes error checking

reset

# Configuration parameters
INPUT_DIR = "build/output/"
OUTPUT_FILE = "wigner_evolution.gif"
DELAY = 8          # Animation delay (1/100 seconds)
SIZE_X = 1000
SIZE_Y = 700

# Set terminal for GIF animation
set terminal gif animate delay DELAY size SIZE_X,SIZE_Y enhanced
set output OUTPUT_FILE

# Plot aesthetics
set xlabel 'Position x (a.u.)' font "Arial,14"
set ylabel 'Momentum p (a.u.)' font "Arial,14"
set cblabel 'Wigner Function W(x,p)' font "Arial,12"

# Set axis limits as requested
set xrange [-5:5]
set yrange [-15:15]

# Improved color palette (blue-white-red)
set palette defined (0 "blue", 0.5 "white", 1 "red")
set pm3d map
set view map

# Colorbar settings
set colorbox
set colorbox vertical user origin 0.85,0.15 size 0.05,0.7
set format cb "%.2e"

# Grid and tics
set grid
set mxtics 5
set mytics 5

# Function to check if file exists (system dependent)
file_exists(file) = system("test -f ".file) == 0 ? 1 : 0

# Main animation loop
# First, let's try common patterns
do for [step=0:200] {
    # Try different file naming patterns
    filename1 = sprintf("%swigner_%011d.dat", INPUT_DIR, step*5)
    filename2 = sprintf("%swigner_%08d.dat", INPUT_DIR, step*5)
    filename3 = sprintf("%swigner_%d.dat", INPUT_DIR, step*5)
    
    # Use the first pattern by default (adjust if needed)
    filename = filename1
    
    # Update title
    time_value = step * 5 * 0.0001  # Assuming dt=0.0001, adjust as needed
    set title sprintf('Wigner Function Evolution - t = %.4f a.u.', time_value) font "Arial,16"
    
    # Plot if file exists
    if (file_exists(filename)) {
        splot filename using 1:2:3 with pm3d notitle
    } else {
        # If file doesn't exist, break the loop
        break
    }
}

unset output
set terminal x11  # Reset terminal

print sprintf("Animation saved as %s", OUTPUT_FILE)
print "To view: open wigner_evolution.gif with your preferred image viewer"

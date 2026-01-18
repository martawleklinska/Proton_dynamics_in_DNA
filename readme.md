# Double harmonic oscilator dynamics
Harmonic model of the single proton transfer. Key aspects in the dynamics were the constant continuity of both forms with the barrier.
The project consists of two parts: polynomial derivation and calculation of analytical Wigner function.
Calculation of WDF via Moyal equation (directly) - set an initial condition for the system to be a gaussian made out of canonical states (for for example a superposition  of ground state and 1sth excited); solve the moyal equation; calculate the energy phase integral. 

Initialization via CMake and shell
```bash
cd moyal_solver/
./build.sh
cd build/
./main
'''

## TODO
Check what kind of a initial condition to use and how to caption the tautomerisation.
set terminal pngcairo size 1300,850 enhanced font "Arial,12"
set datafile commentschars "#"
set grid
set key outside

set output "bead_open_chain_reservoir_temperatures.png"
set xlabel "time (ps)"
set ylabel "temperature (K)"
set title "Coarse-grained open chain: Ar reservoir temperatures"
plot "bead_open_chain_observables.dat" using 2:3 with lines lw 2 title "Ar 1", \
     "bead_open_chain_observables.dat" using 2:4 with lines lw 2 title "Ar 2", \
     "bead_open_chain_observables.dat" using 2:5 with lines lw 2 title "Ar 3"

set output "bead_open_chain_bead_temperatures.png"
set xlabel "time (ps)"
set ylabel "temperature (K)"
set title "Coarse-grained open chain: bead kinetic temperatures"
plot "bead_open_chain_observables.dat" using 2:6 with lines lw 2 title "bead 1", \
     "bead_open_chain_observables.dat" using 2:7 with lines lw 2 title "bead 2", \
     "bead_open_chain_observables.dat" using 2:8 with lines lw 2 title "bead 3"

set output "bead_open_chain_x_displacements.png"
set xlabel "time (ps)"
set ylabel "x COM displacement (A)"
set title "Coarse-grained open chain: bead x displacements"
plot "bead_open_chain_observables.dat" using 2:12 with lines lw 2 title "u1", \
     "bead_open_chain_observables.dat" using 2:13 with lines lw 2 title "u2", \
     "bead_open_chain_observables.dat" using 2:14 with lines lw 2 title "u3"

set output "bead_open_chain_feedback_forces.png"
set xlabel "time (ps)"
set ylabel "feedback force (eV/A)"
set title "Coarse-grained open chain: controller forces"
plot "bead_open_chain_observables.dat" using 2:15 with lines lw 2 title "F1", \
     "bead_open_chain_observables.dat" using 2:16 with lines lw 2 title "F2", \
     "bead_open_chain_observables.dat" using 2:17 with lines lw 2 title "F3"

set output "bead_open_chain_total_energy.png"
set xlabel "time (ps)"
set ylabel "total energy (eV)"
set title "Coarse-grained open chain: total energy diagnostic"
plot "bead_open_chain_observables.dat" using 2:18 with lines lw 2 title "Etot"

set output "bead_open_chain_transverse_displacements.png"
set xlabel "time (ps)"
set ylabel "transverse displacement (A)"
set title "Coarse-grained open chain: bead transverse confinement"
plot "bead_open_chain_observables.dat" using 2:19 with lines lw 2 title "uy1", \
     "bead_open_chain_observables.dat" using 2:20 with lines lw 2 title "uy2", \
     "bead_open_chain_observables.dat" using 2:21 with lines lw 2 title "uy3", \
     "bead_open_chain_observables.dat" using 2:22 with lines lw 2 title "uz1", \
     "bead_open_chain_observables.dat" using 2:23 with lines lw 2 title "uz2", \
     "bead_open_chain_observables.dat" using 2:24 with lines lw 2 title "uz3"

set output "bead_open_chain_feedback_power.png"
set xlabel "time (ps)"
set ylabel "feedback power (eV/ps)"
set title "Coarse-grained open chain: instantaneous controller power"
plot "bead_open_chain_observables.dat" using 2:25 with lines lw 2 title "P1", \
     "bead_open_chain_observables.dat" using 2:26 with lines lw 2 title "P2", \
     "bead_open_chain_observables.dat" using 2:27 with lines lw 2 title "P3", \
     "bead_open_chain_observables.dat" using 2:28 with lines lw 2 title "P total"

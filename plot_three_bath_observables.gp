set terminal pngcairo size 1300,850 enhanced font "Arial,12"
set datafile commentschars "#"
set grid
set key outside

set output "three_bath_temperatures.png"
set xlabel "time (ps)"
set ylabel "temperature (K)"
set title "Self-consistent Ar reservoir temperatures"
plot "three_bath_observables.dat" using 2:3 with lines lw 1.5 title "ambient leak T", \
     "three_bath_observables.dat" using 2:5 with lines lw 2 title "Ar 1", \
     "three_bath_observables.dat" using 2:6 with lines lw 2 title "Ar 2", \
     "three_bath_observables.dat" using 2:7 with lines lw 2 title "Ar 3"

set output "three_bath_nanoparticle_temperatures.png"
set xlabel "time (ps)"
set ylabel "temperature (K)"
set title "Cu nanoparticle kinetic temperatures"
plot "three_bath_observables.dat" using 2:8 with lines lw 2 title "Cu NP 1", \
     "three_bath_observables.dat" using 2:9 with lines lw 2 title "Cu NP 2", \
     "three_bath_observables.dat" using 2:10 with lines lw 2 title "Cu NP 3"

set output "three_bath_com_displacements.png"
set xlabel "time (ps)"
set ylabel "COM displacement x - x0 (A)"
set title "Nanoparticle COM displacements"
plot "three_bath_observables.dat" using 2:14 with lines lw 2 title "u1", \
     "three_bath_observables.dat" using 2:15 with lines lw 2 title "u2", \
     "three_bath_observables.dat" using 2:16 with lines lw 2 title "u3"

set output "three_bath_feedback_forces.png"
set xlabel "time (ps)"
set ylabel "COM feedback force (eV/A)"
set title "Nonreciprocal controller forces"
plot "three_bath_observables.dat" using 2:17 with lines lw 2 title "Fcom1", \
     "three_bath_observables.dat" using 2:18 with lines lw 2 title "Fcom2", \
     "three_bath_observables.dat" using 2:19 with lines lw 2 title "Fcom3"

set output "three_bath_leakage_energy.png"
set xlabel "time (ps)"
set ylabel "cumulative Langevin leakage energy tally (eV)"
set title "Energy exchanged with weak external leakage baths"
plot "three_bath_observables.dat" using 2:20 with lines lw 2 title "Qleak1", \
     "three_bath_observables.dat" using 2:21 with lines lw 2 title "Qleak2", \
     "three_bath_observables.dat" using 2:22 with lines lw 2 title "Qleak3"

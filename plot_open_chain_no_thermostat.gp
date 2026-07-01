set terminal pngcairo size 1300,850 enhanced font "Arial,12"
set datafile commentschars "#"
set grid
set key outside

set output "open_chain_reservoir_temperatures.png"
set xlabel "time (ps)"
set ylabel "temperature (K)"
set title "Open-chain control: Ar reservoir temperatures, no production thermostat"
plot "open_chain_no_thermostat_observables.dat" using 2:3 with lines lw 2 title "Ar 1", \
     "open_chain_no_thermostat_observables.dat" using 2:4 with lines lw 2 title "Ar 2", \
     "open_chain_no_thermostat_observables.dat" using 2:5 with lines lw 2 title "Ar 3"

set output "open_chain_nanoparticle_temperatures.png"
set xlabel "time (ps)"
set ylabel "temperature (K)"
set title "Open-chain control: Cu nanoparticle kinetic temperatures"
plot "open_chain_no_thermostat_observables.dat" using 2:6 with lines lw 2 title "Cu NP 1", \
     "open_chain_no_thermostat_observables.dat" using 2:7 with lines lw 2 title "Cu NP 2", \
     "open_chain_no_thermostat_observables.dat" using 2:8 with lines lw 2 title "Cu NP 3"

set output "open_chain_com_displacements.png"
set xlabel "time (ps)"
set ylabel "COM displacement x - x0 (A)"
set title "Open-chain control: nanoparticle COM displacements"
plot "open_chain_no_thermostat_observables.dat" using 2:12 with lines lw 2 title "u1", \
     "open_chain_no_thermostat_observables.dat" using 2:13 with lines lw 2 title "u2", \
     "open_chain_no_thermostat_observables.dat" using 2:14 with lines lw 2 title "u3"

set output "open_chain_feedback_forces.png"
set xlabel "time (ps)"
set ylabel "COM feedback force (eV/A)"
set title "Open-chain control: nonreciprocal controller forces"
plot "open_chain_no_thermostat_observables.dat" using 2:15 with lines lw 2 title "Fcom1", \
     "open_chain_no_thermostat_observables.dat" using 2:16 with lines lw 2 title "Fcom2", \
     "open_chain_no_thermostat_observables.dat" using 2:17 with lines lw 2 title "Fcom3"

set output "open_chain_total_energy.png"
set xlabel "time (ps)"
set ylabel "total energy (eV)"
set title "Open-chain control: total energy diagnostic"
plot "open_chain_no_thermostat_observables.dat" using 2:18 with lines lw 2 title "Etot"

set output "open_chain_transverse_displacements.png"
set xlabel "time (ps)"
set ylabel "transverse COM displacement (A)"
set title "Open-chain control: transverse nanoparticle COM confinement"
plot "open_chain_no_thermostat_observables.dat" using 2:19 with lines lw 2 title "uy1", \
     "open_chain_no_thermostat_observables.dat" using 2:20 with lines lw 2 title "uy2", \
     "open_chain_no_thermostat_observables.dat" using 2:21 with lines lw 2 title "uy3", \
     "open_chain_no_thermostat_observables.dat" using 2:22 with lines lw 2 dt 2 title "uz1", \
     "open_chain_no_thermostat_observables.dat" using 2:23 with lines lw 2 dt 2 title "uz2", \
     "open_chain_no_thermostat_observables.dat" using 2:24 with lines lw 2 dt 2 title "uz3"

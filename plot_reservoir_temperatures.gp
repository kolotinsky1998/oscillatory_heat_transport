set terminal pngcairo size 1200,800 enhanced font "Arial,12"
set datafile commentschars "#"
set key outside
set grid

set output "reservoir_temperatures.png"
set xlabel "time"
set ylabel "reservoir target temperature"
set title "Oscillatory heat-reservoir temperatures"
plot "reservoir_temperatures.dat" using 2:3 with lines lw 2 title "T1", \
     "reservoir_temperatures.dat" using 2:4 with lines lw 2 title "T2", \
     "reservoir_temperatures.dat" using 2:5 with lines lw 2 title "T3"

set output "oscillator_positions.png"
set xlabel "time"
set ylabel "x"
set title "Coarse-grained nanoparticle coordinates"
plot "reservoir_temperatures.dat" using 2:6 with lines lw 2 title "x1", \
     "reservoir_temperatures.dat" using 2:7 with lines lw 2 title "x2", \
     "reservoir_temperatures.dat" using 2:8 with lines lw 2 title "x3"

set output "feedback_forces.png"
set xlabel "time"
set ylabel "nonreciprocal feedback force"
set title "Applied nonreciprocal forces"
plot "reservoir_temperatures.dat" using 2:9 with lines lw 2 title "F1", \
     "reservoir_temperatures.dat" using 2:10 with lines lw 2 title "F2", \
     "reservoir_temperatures.dat" using 2:11 with lines lw 2 title "F3"


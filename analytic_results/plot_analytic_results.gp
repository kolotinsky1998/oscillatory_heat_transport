set terminal pngcairo size 1300,850 enhanced font "Arial,12"
set datafile separator comma
set key outside
set grid

set output "current_lammps_reservoirs_analytic.png"
set xlabel "time (ps)"
set ylabel "temperature (K)"
set title "Current LAMMPS mapping: analytic reservoir temperatures"
plot "current_lammps_analytic.csv" using 1:2 with lines lw 2 title "Tar1 analytic", \
     "current_lammps_analytic.csv" using 1:3 with lines lw 2 title "Tar2 analytic", \
     "current_lammps_analytic.csv" using 1:4 with lines lw 2 title "Tar3 analytic"

set output "current_lammps_nanoparticles_analytic.png"
set xlabel "time (ps)"
set ylabel "temperature (K)"
set title "Current LAMMPS mapping: analytic Cu kinetic temperatures"
plot "current_lammps_analytic.csv" using 1:5 with lines lw 2 title "Tnp1 analytic", \
     "current_lammps_analytic.csv" using 1:6 with lines lw 2 title "Tnp2 analytic", \
     "current_lammps_analytic.csv" using 1:7 with lines lw 2 title "Tnp3 analytic"

set output "centered_excitation_reservoirs_analytic.png"
set xlabel "time (ps)"
set ylabel "temperature (K)"
set title "Centered excitation: analytic reservoir temperatures"
plot "centered_excitation_analytic.csv" using 1:2 with lines lw 2 title "Tar1 analytic", \
     "centered_excitation_analytic.csv" using 1:3 with lines lw 2 title "Tar2 analytic", \
     "centered_excitation_analytic.csv" using 1:4 with lines lw 2 title "Tar3 analytic"

set output "centered_excitation_nanoparticles_analytic.png"
set xlabel "time (ps)"
set ylabel "temperature (K)"
set title "Centered excitation: analytic Cu kinetic temperatures"
plot "centered_excitation_analytic.csv" using 1:5 with lines lw 2 title "Tnp1 analytic", \
     "centered_excitation_analytic.csv" using 1:6 with lines lw 2 title "Tnp2 analytic", \
     "centered_excitation_analytic.csv" using 1:7 with lines lw 2 title "Tnp3 analytic"

set output "current_lammps_alpha_scan_mean_reservoir.png"
set xlabel "time (ps)"
set ylabel "mean reservoir temperature (K)"
set title "Current LAMMPS initial condition: alphaHeat scan"
plot "current_lammps_alpha_scan.csv" using 2:(($3+$4+$5)/3.0) every ::0::360 with lines lw 2 title "alpha=1", \
     "current_lammps_alpha_scan.csv" using 2:(($3+$4+$5)/3.0) every ::361::721 with lines lw 2 title "alpha=0.3", \
     "current_lammps_alpha_scan.csv" using 2:(($3+$4+$5)/3.0) every ::722::1082 with lines lw 2 title "alpha=0.1", \
     "current_lammps_alpha_scan.csv" using 2:(($3+$4+$5)/3.0) every ::1083::1443 with lines lw 2 title "alpha=0.03", \
     "current_lammps_alpha_scan.csv" using 2:(($3+$4+$5)/3.0) every ::1444::1804 with lines lw 2 title "alpha=0.01"


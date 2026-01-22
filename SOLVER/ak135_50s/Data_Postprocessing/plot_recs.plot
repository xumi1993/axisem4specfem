set term png linewidth 1  
set output "GRAPHICS/A43_XZ_disp_post_mij_conv0000_Z.png"
set title "colat,lon: 45.5882, 238.87, epidist: 51.454009122882248"
plot "SEISMOGRAMS/A43_XZ_disp_post_mij_conv0000_Z.dat" with lines
set xrange [ 0: 6.6013495E+02];set xlabel 'time [s]';set ylabel 'displacement [m]' 

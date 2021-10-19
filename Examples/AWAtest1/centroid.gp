set term qt font "sans,14"
set grid
set xlabel 'zcentroid(m)'
set ylabel 'xcentroid(m)'
plot 'amom.out' u 6:2 w l lw 2 lt 8

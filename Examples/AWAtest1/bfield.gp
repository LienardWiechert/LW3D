set term qt font "sans,14"
set grid
set xlabel 'z(m)'
set ylabel "By(Tesla)
plot 'extern.out' u 3:5 w l lw 2 lt 8

set term qt font "sans,14"
set grid
set xlabel 'z(micron)'
set ylabel "Ezrad(V/m)
plot 'field00100' using ($22*1.e6):9 w l lw 2 lt 8

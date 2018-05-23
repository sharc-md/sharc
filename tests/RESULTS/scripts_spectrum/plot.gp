set title "Absorption spectrum (Gaussian, FWHM=0.100000 eV)\n1 Initial conditions, diagonal representation"

set xrange [0.000000:10.000000]
set yrange [0.000000:1.000000]
set xlabel 'Energy (eV)'
set ylabel 'Absorption spectrum (normalized)'

set style fill transparent solid 0.25 border
set term pngcairo size 640,480
set out 'spectrum.out.png'

N=0.529855

plot "spectrum.out" u 1:(($1)/N) w filledcu tit "State 1" lw 0.5 lc rgbcolor "#FF0000", \
     ""             u 1:(($2)/N) w filledcu tit "State 2" lw 0.5 lc rgbcolor "#51FF00", \
     ""             u 1:(($3)/N) w filledcu tit "State 3" lw 0.5 lc rgbcolor "#00FFC1", \
     ""             u 1:(($4)/N) w filledcu tit "State 4" lw 0.5 lc rgbcolor "#0028FF", \
     ""             u 1:(($5)/N) w filledcu tit "State 5" lw 0.5 lc rgbcolor "#EA00FF", \
     ""             u 1:(($7)/N)   w l        tit "Sum"       lw 3.0 lc rgbcolor "black"


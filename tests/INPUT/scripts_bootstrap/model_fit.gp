#>>
# +------------------------------------------------------------+
# |                       Fitting script                       |
# +------------------------------------------------------------+
#
#
# *** Definition of the kinetic model: ***
#     |  S1   T1   T2 
#-----+---------------
#  S1 |   .    .   k1 
#  T1 |   .    .    . 
#  T2 |   .   k2    . 
#
#(Initial species: rows; Final species: columns)
#
#
# *** Species and initial value: ***
# S1        S1__0
# T1        0
# T2        0
#
#
# *** Reaction rates: ***
# k2
# k1
#
#
# *** Species population function definitions: ***
S1(t) = exp(-k1*t)*S1__0

T1(t) = k1*exp(-k2*t)*S1__0/(k2-k1)-k2*exp(-k1*t)*S1__0/(k2-k1)+S1__0

T2(t) = k1*exp(-k1*t)*S1__0/(k2-k1)-k1*exp(-k2*t)*S1__0/(k2-k1)



# ========================================================
# *** Reaction rates initial guesses: ***
# These are given in units of inverse fs (e.g., time constant of 100 fs is written as 1./100.).
# TODO: Please change to some suitable values!
k2 = 1./100.
k1 = 1./101.

# *** Species population initial guesses: ***
# TODO: Please change to some suitable values!
S1__0 = 1./1.

#<<
# *** Gnuplot general options: ***
set xlabel "Time (fs)"
set ylabel "Population"
set xrange [0:1000.00]
set yrange [0:1]
set key at 1000.00,1.00 top right

# *** Label with time constants: ***
set label 1 sprintf("t_k2 = %7.1f fs\nt_k1 = %7.1f fs\nS1__0 = %7.1f\n",1./k2,1./k1,S1__0) at 300.00,0.95

# *** Initial plot before fit: ***
set title "Initial plot before fit\nPress ENTER to setup global fit."
p \
  S1(x) t "S1" w l lw 2 lc rgbcolor "#FF0000", \
  T1(x) t "T1" w l lw 2 lc rgbcolor "#00FF66", \
  T2(x) t "T2" w l lw 2 lc rgbcolor "#3200FF", \
  "model_fit.dat" u 1:($3) t "$3" w l lw 0.5 lc rgbcolor "#FF0000", \
  "model_fit.dat" u 1:($4) t "$4" w l lw 0.5 lc rgbcolor "#00FF66", \
  "model_fit.dat" u 1:($5) t "$5" w l lw 0.5 lc rgbcolor "#3200FF"

pause -1

#>>
# *** Global fit setup: ***
set xrange [0:3000.00]
unset label
set yrange [0:1]
set key at 3000.00,1.00 top right
set fit logfile "model_fit.log"
F(x)= x<1000.00 ? S1(x-0.00) : ( x<2000.00 ? T1(x-1000.00) : T2(x-2000.00) ) 

#<<
# *** Global plot before fitting: ***
set title "Global plot before fitting\nPress ENTER to perform fit."
p "model_fit.dat" u 1:( $1<1000.00 ? $3 : ( $1<2000.00 ? $4 : $5 )  ) t "Data" w p pt 6 ps 0.5 lc rgbcolor "black", F(x) t "Fitting function" w l lw 2 lc rgbcolor "red"
pause -1

#>>
# *** Execute global fit: ***
# TODO: remove the "#" on the next line to also fit the initial populations.
fit F(x) "model_fit.dat" u 1:( $1<1000.00 ? $3 : ( $1<2000.00 ? $4 : $5 )  ) via k2,k1,S1__0

#<<
# *** Global plot after fitting: ***
set title "Global plot after fitting\nPress ENTER to display final results."
p "model_fit.dat" u 1:( $1<1000.00 ? $3 : ( $1<2000.00 ? $4 : $5 )  ) t "Data" w p pt 6 ps 0.5 lc rgbcolor "black", F(x) t "Fitting function" w l lw 2 lc rgbcolor "red"
pause -1

#>>
# *** Gnuplot general options: ***
set xlabel "Time (fs)"
set ylabel "Population"
set xrange [0:1000.00]
set yrange [0:1]
set key at 1000.00,1.00 top right

# *** Label with time constants: ***
set label 1 sprintf("t_k2 = %7.1f fs\nt_k1 = %7.1f fs\nS1__0 = %7.1f\n",1./k2,1./k1,S1__0) at 300.00,0.95

#<><> 1
# *** Final plot after fit: ***
set title "Final plot after fitting\nPress ENTER to save as PNG and TXT."
p \
  S1(x) t "S1" w l lw 2 lc rgbcolor "#FF0000", \
  T1(x) t "T1" w l lw 2 lc rgbcolor "#00FF66", \
  T2(x) t "T2" w l lw 2 lc rgbcolor "#3200FF", \
  "model_fit.dat" u 1:($3) t "$3" w l lw 0.5 lc rgbcolor "#FF0000", \
  "model_fit.dat" u 1:($4) t "$4" w l lw 0.5 lc rgbcolor "#00FF66", \
  "model_fit.dat" u 1:($5) t "$5" w l lw 0.5 lc rgbcolor "#3200FF"
#<<

pause -1

#>>
# *** File output: ***
#<<
set title "mai@iron.itc.univie.ac.at, 2017-11-29 12:58:14.846512\nFile:/user/mai/Documents/NewSHARC/SHARC_2.0/bootstrap/pop.out"
set term pngcairo size 800,480
set out "model_fit.png"
replot
#>>
set table "model_fit.txt"
replot

print "&&& k2 : ", 1./k2
print "&&& k1 : ", 1./k1
print "&&& S1__0 : ", S1__0
# *** Infos: ***
# mai@iron.itc.univie.ac.at
# Date: 2017-11-29 12:58:14.846551
# Current directory: /user/mai/Documents/NewSHARC/SHARC_2.0/bootstrap


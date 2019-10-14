unset key
unset colorbox
set title "Energies in diagonal basis"
set xlabel "Time t in fs"
set ylabel "Energy in eV"
set cbrange [0:17]
set palette defined (0.0 "gray90", 1e-5 "gray60", 1e-4 "gray30", 1e-3 "orange", 1e-2 "red", 1e-1 "magenta", 1e-0 "blue", 10 "blue", 11 "green", 12 "red", 13 "turquoise", 14 "orange", 15 "cyan", 16 "brown", 17 "skyblue")

plot "output_data/expec.out" u 1:($4)               title "Total Energy" lw   0.50 lc rgbcolor "#000000" w l, \
""               u 1:($5):(abs($23)+10) title "State 1"     lw   4.50 pal w l, \
""               u 1:($5):(abs($41))    title "State 1"     lw   3.50 pal w l, \
""               u 1:($6):(abs($24)+10) title "State 2"     lw   4.50 pal w l, \
""               u 1:($6):(abs($42))    title "State 2"     lw   3.50 pal w l, \
""               u 1:($7):(abs($25)+10) title "State 3"     lw   4.50 pal w l, \
""               u 1:($7):(abs($43))    title "State 3"     lw   3.50 pal w l, \
""               u 1:($8):(abs($26)+10) title "State 4"     lw   4.50 pal w l, \
""               u 1:($8):(abs($44))    title "State 4"     lw   3.50 pal w l, \
""               u 1:($9):(abs($27)+10) title "State 5"     lw   4.50 pal w l, \
""               u 1:($9):(abs($45))    title "State 5"     lw   3.50 pal w l, \
""               u 1:($10):(abs($28)+10) title "State 6"     lw   4.50 pal w l, \
""               u 1:($10):(abs($46))    title "State 6"     lw   3.50 pal w l, \
""               u 1:($11):(abs($29)+10) title "State 7"     lw   4.50 pal w l, \
""               u 1:($11):(abs($47))    title "State 7"     lw   3.50 pal w l, \
""               u 1:($12):(abs($30)+10) title "State 8"     lw   4.50 pal w l, \
""               u 1:($12):(abs($48))    title "State 8"     lw   3.50 pal w l, \
""               u 1:($13):(abs($31)+10) title "State 9"     lw   4.50 pal w l, \
""               u 1:($13):(abs($49))    title "State 9"     lw   3.50 pal w l, \
""               u 1:($14):(abs($32)+10) title "State 10"     lw   4.50 pal w l, \
""               u 1:($14):(abs($50))    title "State 10"     lw   3.50 pal w l, \
""               u 1:($15):(abs($33)+10) title "State 11"     lw   4.50 pal w l, \
""               u 1:($15):(abs($51))    title "State 11"     lw   3.50 pal w l, \
""               u 1:($16):(abs($34)+10) title "State 12"     lw   4.50 pal w l, \
""               u 1:($16):(abs($52))    title "State 12"     lw   3.50 pal w l, \
""               u 1:($17):(abs($35)+10) title "State 13"     lw   4.50 pal w l, \
""               u 1:($17):(abs($53))    title "State 13"     lw   3.50 pal w l, \
""               u 1:($18):(abs($36)+10) title "State 14"     lw   4.50 pal w l, \
""               u 1:($18):(abs($54))    title "State 14"     lw   3.50 pal w l, \
""               u 1:($19):(abs($37)+10) title "State 15"     lw   4.50 pal w l, \
""               u 1:($19):(abs($55))    title "State 15"     lw   3.50 pal w l, \
""               u 1:($20):(abs($38)+10) title "State 16"     lw   4.50 pal w l, \
""               u 1:($20):(abs($56))    title "State 16"     lw   3.50 pal w l, \
""               u 1:($21):(abs($39)+10) title "State 17"     lw   4.50 pal w l, \
""               u 1:($21):(abs($57))    title "State 17"     lw   3.50 pal w l, \
""               u 1:($22):(abs($40)+10) title "State 18"     lw   4.50 pal w l, \
""               u 1:($22):(abs($58))    title "State 18"     lw   3.50 pal w l, \
""               u 1:($3)               title "Trajectory"   lw   1.00 lc rgbcolor "#000000" pt 6 w p

pause -1

unset key
unset colorbox
set title "Energies in diagonal basis"
set xlabel "Time t in fs"
set ylabel "Energy in eV"
set cbrange [0:17]
set palette defined (0.0 "gray90", 1e-5 "gray60", 1e-4 "gray30", 1e-3 "orange", 1e-2 "red", 1e-1 "magenta", 1e-0 "blue", 10 "blue", 11 "green", 12 "red", 13 "turquoise", 14 "orange", 15 "cyan", 16 "brown", 17 "skyblue")

plot "output_data/expec.out" u 1:($4-$5)               title "Total Energy" lw   0.50 lc rgbcolor "#000000" w l, \
""               u 1:($5-$5):(abs($23)+10) title "State 1"     lw   4.50 pal w l, \
""               u 1:($5-$5):(abs($41))    title "State 1"     lw   3.50 pal w l, \
""               u 1:($6-$5):(abs($24)+10) title "State 2"     lw   4.50 pal w l, \
""               u 1:($6-$5):(abs($42))    title "State 2"     lw   3.50 pal w l, \
""               u 1:($7-$5):(abs($25)+10) title "State 3"     lw   4.50 pal w l, \
""               u 1:($7-$5):(abs($43))    title "State 3"     lw   3.50 pal w l, \
""               u 1:($8-$5):(abs($26)+10) title "State 4"     lw   4.50 pal w l, \
""               u 1:($8-$5):(abs($44))    title "State 4"     lw   3.50 pal w l, \
""               u 1:($9-$5):(abs($27)+10) title "State 5"     lw   4.50 pal w l, \
""               u 1:($9-$5):(abs($45))    title "State 5"     lw   3.50 pal w l, \
""               u 1:($10-$5):(abs($28)+10) title "State 6"     lw   4.50 pal w l, \
""               u 1:($10-$5):(abs($46))    title "State 6"     lw   3.50 pal w l, \
""               u 1:($11-$5):(abs($29)+10) title "State 7"     lw   4.50 pal w l, \
""               u 1:($11-$5):(abs($47))    title "State 7"     lw   3.50 pal w l, \
""               u 1:($12-$5):(abs($30)+10) title "State 8"     lw   4.50 pal w l, \
""               u 1:($12-$5):(abs($48))    title "State 8"     lw   3.50 pal w l, \
""               u 1:($13-$5):(abs($31)+10) title "State 9"     lw   4.50 pal w l, \
""               u 1:($13-$5):(abs($49))    title "State 9"     lw   3.50 pal w l, \
""               u 1:($14-$5):(abs($32)+10) title "State 10"     lw   4.50 pal w l, \
""               u 1:($14-$5):(abs($50))    title "State 10"     lw   3.50 pal w l, \
""               u 1:($15-$5):(abs($33)+10) title "State 11"     lw   4.50 pal w l, \
""               u 1:($15-$5):(abs($51))    title "State 11"     lw   3.50 pal w l, \
""               u 1:($16-$5):(abs($34)+10) title "State 12"     lw   4.50 pal w l, \
""               u 1:($16-$5):(abs($52))    title "State 12"     lw   3.50 pal w l, \
""               u 1:($17-$5):(abs($35)+10) title "State 13"     lw   4.50 pal w l, \
""               u 1:($17-$5):(abs($53))    title "State 13"     lw   3.50 pal w l, \
""               u 1:($18-$5):(abs($36)+10) title "State 14"     lw   4.50 pal w l, \
""               u 1:($18-$5):(abs($54))    title "State 14"     lw   3.50 pal w l, \
""               u 1:($19-$5):(abs($37)+10) title "State 15"     lw   4.50 pal w l, \
""               u 1:($19-$5):(abs($55))    title "State 15"     lw   3.50 pal w l, \
""               u 1:($20-$5):(abs($38)+10) title "State 16"     lw   4.50 pal w l, \
""               u 1:($20-$5):(abs($56))    title "State 16"     lw   3.50 pal w l, \
""               u 1:($21-$5):(abs($39)+10) title "State 17"     lw   4.50 pal w l, \
""               u 1:($21-$5):(abs($57))    title "State 17"     lw   3.50 pal w l, \
""               u 1:($22-$5):(abs($40)+10) title "State 18"     lw   4.50 pal w l, \
""               u 1:($22-$5):(abs($58))    title "State 18"     lw   3.50 pal w l, \
""               u 1:($3-$5)               title "Trajectory"   lw   1.00 lc rgbcolor "#000000" pt 6 w p

pause -1

set key
set yrange [0:1]
set ylabel "Wavefunction Amplitude"
set title "MCH Quantum Amplitudes"
plot "output_data/coeff_MCH.out" 	u 1:2  	title "Sum of Amplitudes" 	lw   3.00 	lc rgbcolor "#000000" 	w l, \
"" 			u 1:($3**2+$4**2)  	title "S 0, 0" 	lw   1.00 	lc rgbcolor "#FF0000" 	w l, \
"" 			u 1:($5**2+$6**2)  	title "S 1, 0" 	lw   1.00 	lc rgbcolor "#FF7200" 	w l, \
"" 			u 1:($7**2+$8**2)  	title "S 2, 0" 	lw   1.00 	lc rgbcolor "#7FFF00" 	w l, \
"" 			u 1:($9**2+$10**2)  	title "S 3, 0" 	lw   1.00 	lc rgbcolor "#0CFF00" 	w l, \
"" 			u 1:($11**2+$12**2)  	title "T 1, -1" 	lw   3.00 	lc rgbcolor "#00FF66" 	w l, \
"" 			u 1:($13**2+$14**2)  	title "T 2, -1" 	lw   3.00 	lc rgbcolor "#00FFFF" 	w l, \
"" 			u 1:($15**2+$16**2)  	title "T 3, -1" 	lw   3.00 	lc rgbcolor "#0066FF" 	w l, \
"" 			u 1:($17**2+$18**2)  	title "T 1, 0" 	lw   3.00 	lc rgbcolor "#00FF66" 	w l, \
"" 			u 1:($19**2+$20**2)  	title "T 2, 0" 	lw   3.00 	lc rgbcolor "#00FFFF" 	w l, \
"" 			u 1:($21**2+$22**2)  	title "T 3, 0" 	lw   3.00 	lc rgbcolor "#0066FF" 	w l, \
"" 			u 1:($23**2+$24**2)  	title "T 1, 1" 	lw   3.00 	lc rgbcolor "#00FF66" 	w l, \
"" 			u 1:($25**2+$26**2)  	title "T 2, 1" 	lw   3.00 	lc rgbcolor "#00FFFF" 	w l, \
"" 			u 1:($27**2+$28**2)  	title "T 3, 1" 	lw   3.00 	lc rgbcolor "#0066FF" 	w l, \
"" 			u 1:($29**2+$30**2)  	title "Q 1, -2" 	lw   5.00 	lc rgbcolor "#3200FF" 	w l, \
"" 			u 1:($31**2+$32**2)  	title "Q 1, -1" 	lw   5.00 	lc rgbcolor "#3200FF" 	w l, \
"" 			u 1:($33**2+$34**2)  	title "Q 1, 0" 	lw   5.00 	lc rgbcolor "#3200FF" 	w l, \
"" 			u 1:($35**2+$36**2)  	title "Q 1, 1" 	lw   5.00 	lc rgbcolor "#3200FF" 	w l, \
"" 			u 1:($37**2+$38**2)  	title "Q 1, 2" 	lw   5.00 	lc rgbcolor "#3200FF" 	w l

pause -1

set yrange [0:1]
set ylabel "Wavefunction Amplitude"
set title "Diag Quantum Amplitudes"
plot "output_data/coeff_diag.out" 	u 1:2  	title "Sum of Amplitudes" 	lw   3.00 	lc rgbcolor "#000000" 	w l, \
"" 			u 1:($3**2+$4**2)  	title "State 1" 	lw   2.00 	lc rgbcolor "#FF0000" 	w l, \
"" 			u 1:($5**2+$6**2)  	title "State 2" 	lw   2.00 	lc rgbcolor "#FF4C00" 	w l, \
"" 			u 1:($7**2+$8**2)  	title "State 3" 	lw   2.00 	lc rgbcolor "#FF9900" 	w l, \
"" 			u 1:($9**2+$10**2)  	title "State 4" 	lw   2.00 	lc rgbcolor "#7FFF00" 	w l, \
"" 			u 1:($11**2+$12**2)  	title "State 5" 	lw   2.00 	lc rgbcolor "#32FF00" 	w l, \
"" 			u 1:($13**2+$14**2)  	title "State 6" 	lw   2.00 	lc rgbcolor "#00FF19" 	w l, \
"" 			u 1:($15**2+$16**2)  	title "State 7" 	lw   2.00 	lc rgbcolor "#00FF66" 	w l, \
"" 			u 1:($17**2+$18**2)  	title "State 8" 	lw   2.00 	lc rgbcolor "#00FFB2" 	w l, \
"" 			u 1:($19**2+$20**2)  	title "State 9" 	lw   2.00 	lc rgbcolor "#00FFFF" 	w l, \
"" 			u 1:($21**2+$22**2)  	title "State 10" 	lw   2.00 	lc rgbcolor "#00B2FF" 	w l, \
"" 			u 1:($23**2+$24**2)  	title "State 11" 	lw   2.00 	lc rgbcolor "#0066FF" 	w l, \
"" 			u 1:($25**2+$26**2)  	title "State 12" 	lw   2.00 	lc rgbcolor "#0019FF" 	w l, \
"" 			u 1:($27**2+$28**2)  	title "State 13" 	lw   2.00 	lc rgbcolor "#3300FF" 	w l, \
"" 			u 1:($29**2+$30**2)  	title "State 14" 	lw   2.00 	lc rgbcolor "#7F00FF" 	w l, \
"" 			u 1:($31**2+$32**2)  	title "State 15" 	lw   2.00 	lc rgbcolor "#CC00FF" 	w l, \
"" 			u 1:($33**2+$34**2)  	title "State 16" 	lw   2.00 	lc rgbcolor "#FF00E5" 	w l, \
"" 			u 1:($35**2+$36**2)  	title "State 17" 	lw   2.00 	lc rgbcolor "#FF0098" 	w l, \
"" 			u 1:($37**2+$38**2)  	title "State 18" 	lw   2.00 	lc rgbcolor "#FF004C" 	w l

pause -1

set xrange[0:*]
set ylabel "Hopping Probability"
set title "Diag Hopping Probablities and Random Number"
set style fill solid 0.25 border
plot "output_data/prob.out"	u 1:20  	title "State 18" 	lw   1.00 	lc rgbcolor "#FF004C" 	w boxes, \
""		u 1:19  	title "State 17" 	lw   1.00 	lc rgbcolor "#FF0098" 	w boxes, \
""		u 1:18  	title "State 16" 	lw   1.00 	lc rgbcolor "#FF00E5" 	w boxes, \
""		u 1:17  	title "State 15" 	lw   1.00 	lc rgbcolor "#CC00FF" 	w boxes, \
""		u 1:16  	title "State 14" 	lw   1.00 	lc rgbcolor "#7F00FF" 	w boxes, \
""		u 1:15  	title "State 13" 	lw   1.00 	lc rgbcolor "#3300FF" 	w boxes, \
""		u 1:14  	title "State 12" 	lw   1.00 	lc rgbcolor "#0019FF" 	w boxes, \
""		u 1:13  	title "State 11" 	lw   1.00 	lc rgbcolor "#0066FF" 	w boxes, \
""		u 1:12  	title "State 10" 	lw   1.00 	lc rgbcolor "#00B2FF" 	w boxes, \
""		u 1:11  	title "State 9" 	lw   1.00 	lc rgbcolor "#00FFFF" 	w boxes, \
""		u 1:10  	title "State 8" 	lw   1.00 	lc rgbcolor "#00FFB2" 	w boxes, \
""		u 1:9  	title "State 7" 	lw   1.00 	lc rgbcolor "#00FF66" 	w boxes, \
""		u 1:8  	title "State 6" 	lw   1.00 	lc rgbcolor "#00FF19" 	w boxes, \
""		u 1:7  	title "State 5" 	lw   1.00 	lc rgbcolor "#32FF00" 	w boxes, \
""		u 1:6  	title "State 4" 	lw   1.00 	lc rgbcolor "#7FFF00" 	w boxes, \
""		u 1:5  	title "State 3" 	lw   1.00 	lc rgbcolor "#FF9900" 	w boxes, \
""		u 1:4  	title "State 2" 	lw   1.00 	lc rgbcolor "#FF4C00" 	w boxes, \
""		u 1:3  	title "State 1" 	lw   1.00 	lc rgbcolor "#FF0000" 	w boxes, \
""		u 1:($1>0 ? $2 : 1/0)	title "Random number" 			lc rgbcolor "black"	w lp

pause -1



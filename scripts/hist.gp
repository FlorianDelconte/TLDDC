reset
max=30 #max value
min=-30 #min value
set size ratio -1
#width=0.01 #interval width
#n=(max - min)/width
#function used to map a value to the intervals
#hist(x,width)=width*floor(x/width)+width/2.0
set encoding utf8
#set terminal postscript enhanced color font ",14"
set terminal postscript size 4.5,3.5 enhanced color fontfile "/usr/share/texmf-dist/fonts/type1/public/cm-super/sfss1200.pfb" "SFSS1200" 10
set output '| ps2pdf - "dataFile.pdf"'

set xrange [min:max]
set yrange [0:]
#to put an empty boundary around the
#data inside an autoscaled graph.
set offset graph 0.05,0.05,0.05,0.0
#set xtics min,(max-min)/5,max
#set boxwidth width*0.9

set style line 103 lc rgb '#ff8800' pt 7 ps 1
set style line 101 lc rgb '#4499cc' lt 1 lw 1.3
set style line 102 lc rgb '#000000' lt 1 lw 1
set style line 100 lc rgb '#3322ee' lt 1 lw 2

set style fill solid 0.5 #fillstyle
set tics out nomirror
#set xlabel "Difference between the distance to the center line and the reference distance to the center line"
set xlabel "Distance"
set ylabel "Frequency"


#count and plot
plot "dataFile" smooth freq w boxes lc rgb "#555555" notitle

set terminal png enhanced
set terminal postscript enhanced
set terminal postscript eps enhanced
set term postscript enhanced color
set output "colorindex.ps"
set style line 1 lt 1 lw 2
set style line 2 lt 2 lw 2
set style line 3 lt 3 lw 2
plot 	"../data/RDF.dat" using 1:2 title "" with line ls 1
set xlabel "r (Å)"
set xrange [0:5]
set title "TiO_2"
set ylabel "RDF G_{TiO_2}(r)" 
set terminal png font arial 23 size 1024,768
set output "../graph/Total RDF.png"
replot


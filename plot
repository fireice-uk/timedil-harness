#!/usr/bin/env bash

plot_script=" 
    #set terminal x11
    set multiplot layout 2,1 rowsfirst
    set lmargin at screen 0.15
    set rmargin at screen 0.85
    set grid x y y2
    set xlabel 'step'
    set ylabel 'normalized difficulty'
    set y2label 'normalized hash rate'
   # set autoscale y
    #set autoscale y2
    set y2tics
    set y2tics nomirror
    set tics out

"

plot_line="plot \"$1\" u 1:15 t 'hash rate' axes x1y2"
plot_line_time="
    unset y2tics
    unset y2label
    set ylabel 'delta blocktime'
    delta_t(x) = ( tD = x - old_t, old_t = x, tD)
    old_t = NaN
    plot \"$1\" u 1:(delta_t(\$5)) t 'delta block timestamp'
    unset multiplot
"

for i in $*
do
    plot_line+=", \"$i\" u 1:9 t \"$i\""
done

echo "$plot_script $plot_line $plot_line_time" | gnuplot -persist

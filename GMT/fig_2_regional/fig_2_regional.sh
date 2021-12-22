# Script to generate Figure 2.3
ps=fig_2_regional.ps
gmt set PS_MEDIA=Custom_16.0cx20.0c
gmt set MAP_LABEL_OFFSET=3p
gmt set MAP_ANNOT_OFFSET=3p
gmt set MAP_FRAME_WIDTH=0.06c
gmt set MAP_TICK_LENGTH_PRIMARY=0.06c
gmt set MAP_GRID_PEN_PRIMARY=default,lightgrey
gmt set MAP_TITLE_OFFSET=0.05c
gmt set MAP_FRAME_PEN=0.5p,40/40/40
gmt set FONT_ANNOT_PRIMARY             = 8p,Hind-Light,40/40/40
gmt set FONT_ANNOT_SECONDARY           = 8p,Hind-Light,40/40/40
gmt set FONT_LABEL                     = 8p,Hind-Light,40/40/40
gmt set FONT_LOGO                      = 8p,Hind-Light,40/40/40
gmt set FONT_TITLE                     = 8p,Hind-Light,40/40/40
gmt set PS_PAGE_ORIENTATION=portrait
gmt set PS_LINE_CAP=butt
gmt set PS_LINE_JOIN=round

R=1970/2050/-0.20/0.80
J=X7c/4.0c
Bx=x20g20
By=y0.2g0.2
color0=60/60/60
color1=#4e79a7
color2=#f28e2b
color3=#e15759
color4=#76b7b2
color5=#59a14f
color6=#edc948
color7=#b07aa1
color8=#ff9da7
color9=#9c755f
color10=#bab0ac

colorlow=#330a5f
colorintlow=#781c6d
colorint=#bb3755
colorinthigh=#ed6925
colorhigh=#fcb519

print_lines () {
gmt psxy  -R -J -L+bD -t50 -N -G$color5 -O -K ${region}_Trajectory.txt >> $ps
gmt psxy  -R -J -L+bD -t90 -N -G$colorlow -O -K ${region}_Low.txt   >> $ps
gmt psxy  -R -J -L+bD -t90 -N -G$colorintlow -O -K ${region}_IntLow.txt   >> $ps
gmt psxy  -R -J -L+bD -t90 -N -G$colorint -O -K ${region}_Int.txt   >> $ps
gmt psxy  -R -J -L+bD -t90 -N -G$colorinthigh -O -K ${region}_IntHigh.txt   >> $ps
gmt psxy  -R -J -L+bD -t90 -N -G$colorhigh -O -K ${region}_High.txt   >> $ps
gmt psxy  -R -J -N -W1.0p,$colorlow -O -K ${region}_Low.txt   >> $ps
gmt psxy  -R -J -N -W1.0p,$colorintlow -O -K ${region}_IntLow.txt   >> $ps
gmt psxy  -R -J -N -W1.0p,$colorint -O -K ${region}_Int.txt   >> $ps
gmt psxy  -R -J -N -W1.0p,$colorinthigh -O -K ${region}_IntHigh.txt   >> $ps
gmt psxy  -R -J -N -W1.0p,$colorhigh -O -K ${region}_High.txt   >> $ps
gmt psxy  -R -J -N -W1.5p,$color5 -O -K ${region}_Trajectory.txt   >> $ps
gmt psxy  -R -J -N -W1.5p,$color0 -O -K ${region}_20c.txt   >> $ps
gmt psxy  -R -J -N -W1.5p,$color1 -O -K ${region}_Altimetry.txt   >> $ps
}

region="EC"
gmt psbasemap -K -R$R -J$J -Y15.6c -X1.0c -BWesn+t'Northeast' -B$Bx -B$By+l'Sea level (m)'  > $ps
print_lines
echo "1970 0.8 a)" | gmt pstext -D0.125c/-0.125c -Gwhite -R -J -F+f9p,Hind-SemiBold+jLT -O -K >> $ps

region="SE"
gmt psbasemap -O -K -R$R -J$J -X7.8c -Bwesn+t'Southeast' -B$Bx -B$By  >> $ps
print_lines
echo "1970 0.8 b)" | gmt pstext -D0.125c/-0.125c -Gwhite -R -J -F+f9p,Hind-SemiBold+jLT -O -K >> $ps

region="GCE"
gmt psbasemap -O -K -R$R -J$J -X-7.8c -Y-4.5c -BWesn+t'Eastern Gulf' -B$Bx -B$By+l'Sea level (m)'  >> $ps
print_lines
echo "1970 0.8 c)" | gmt pstext -D0.125c/-0.125c -Gwhite -R -J -F+f9p,Hind-SemiBold+jLT -O -K >> $ps

region="GCW"
gmt psbasemap -O -K -R$R -J$J -X7.8c -Bwesn+t'Western Gulf' -B$Bx -B$By  >> $ps
print_lines
echo "1970 0.8 d)" | gmt pstext -D0.125c/-0.125c -Gwhite -R -J -F+f9p,Hind-SemiBold+jLT -O -K >> $ps

region="SWC"
gmt psbasemap -O -K -R$R -J$J -X-7.8c -Y-4.5c -BWesn+t'Southwest' -B$Bx -B$By+l'Sea level (m)'  >> $ps
print_lines
echo "1970 0.8 e)" | gmt pstext -D0.125c/-0.125c -Gwhite -R -J -F+f9p,Hind-SemiBold+jLT -O -K >> $ps

region="NWC"
gmt psbasemap -O -K -R$R -J$J -X7.8c -Bwesn+t'Northwest' -B$Bx -B$By  >> $ps
print_lines
echo "1970 0.8 f)" | gmt pstext -D0.125c/-0.125c -Gwhite -R -J -F+f9p,Hind-SemiBold+jLT -O -K >> $ps

region="PAC"
gmt psbasemap -O -K -R$R -J$J -X-7.8c -Y-4.5c -BWeSn+t'Hawaiian Islands' -B$Bx -B$By+l'Sea level (m)'  >> $ps
print_lines
echo "1970 0.8 g)" | gmt pstext -D0.125c/-0.125c -Gwhite -R -J -F+f9p,Hind-SemiBold+jLT -O -K >> $ps

region="CAR"
gmt psbasemap -O -K -R$R -J$J -X7.8c -BweSn+t'Caribbean Islands' -B$Bx -B$By  >> $ps
print_lines
echo "1970 0.8 h)" | gmt pstext -D0.125c/-0.125c -Gwhite -R -J -F+f9p,Hind-SemiBold+jLT -O -K >> $ps

gmt psbasemap -O -K -R0/0.85/0.5/5.5 -JX6.8c/1.5c -X-3.8c -Y-2.05c -T100/100/1  >> $ps
echo -e "0.0 5 \n 0.05 5" | gmt psxy -R -J -O -K -W1.5p,$color1 >> $ps
echo "0.055 5 Altimetry observations " | gmt pstext -R -J -F+f8+jLM -O -K -N >> $ps

echo -e "0.0 4 \n 0.05 4" | gmt psxy -R -J -O -K -W1.5p,$color0 >> $ps
echo "0.055 4 Tide-gauge observations" | gmt pstext -R -J -F+f8+jLM -O -K -N >> $ps

echo -e "0.0 3 \n 0.05 3" | gmt psxy -R -J -O -K -t60 -W6p,$color5 >> $ps
echo -e "0.0 3 \n 0.05 3" | gmt psxy -R -J -O -K -W1.5p,$color5 >> $ps
echo "0.055 3 Present trajectory" | gmt pstext -R -J -F+f8+jLM -O -K -N >> $ps

echo -e "0.5 5 \n 0.55 5" | gmt psxy -R -J -O -K -t90 -W6p,$color0 >> $ps
echo -e "0.5 5 \n 0.55 5" | gmt psxy -R -J -O -K -W1.5p,$color0 >> $ps
echo "0.555 5 Low" | gmt pstext -R -J -F+f8+jLM -O -K -N >> $ps

echo -e "0.5 4 \n 0.55 4" | gmt psxy -R -J -O -K -t90 -W6p,$colorintlow >> $ps
echo -e "0.5 4 \n 0.55 4" | gmt psxy -R -J -O -K -W1.5p,$colorintlow >> $ps
echo "0.555 4 Intermediate-low" | gmt pstext -R -J -F+f8+jLM -O -K -N >> $ps

echo -e "0.5 3 \n 0.55 3" | gmt psxy -R -J -O -K -t90 -W6p,$colorint >> $ps
echo -e "0.5 3 \n 0.55 3" | gmt psxy -R -J -O -K -W1.5p,$colorint >> $ps
echo "0.555 3 Intermediate" | gmt pstext -R -J -F+f8+jLM -O -K -N >> $ps

echo -e "0.5 2 \n 0.55 2" | gmt psxy -R -J -O -K -t90 -W6p,$colorinthigh >> $ps
echo -e "0.5 2 \n 0.55 2" | gmt psxy -R -J -O -K -W1.5p,$colorinthigh >> $ps
echo "0.555 2 Intermediate-high" | gmt pstext -R -J -F+f8+jLM -O -K -N >> $ps

echo -e "0.5 1 \n 0.55 1" | gmt psxy -R -J -O -K -t90 -W6p,$colorhigh >> $ps
echo -e "0.5 1 \n 0.55 1" | gmt psxy -R -J -O -K -W1.5p,$colorhigh >> $ps
echo "0.555 1 High" | gmt pstext -R -J -F+f8+jLM -O -N >> $ps

gmt psconvert -Tf -C-sFONTPATH=/home/thomas/.local/share/fonts/ $ps
gmt psconvert -Tg -C-sFONTPATH=/home/thomas/.local/share/fonts/ $ps



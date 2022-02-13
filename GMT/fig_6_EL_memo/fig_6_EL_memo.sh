print_lines () {
gmt psxy  -R -J -N -W1.0p,$colorlow -O -K ${region}_Low.txt   >> $ps
gmt psxy  -R -J -N -W1.0p,$colorintlow -O -K ${region}_IntLow.txt   >> $ps
gmt psxy  -R -J -N -W1.0p,$colorint -O -K ${region}_Int.txt   >> $ps
gmt psxy  -R -J -N -W1.0p,$colorinthigh -O -K ${region}_IntHigh.txt   >> $ps
gmt psxy  -R -J -N -W1.0p,$colorhigh -O -K ${region}_High.txt   >> $ps
gmt psxy  -R -J -N -W2p,$color5 -O -K ${region}_Trajectory.txt   >> $ps
gmt psxy  -R -J -N -W1.5p,$color0 -O -K ${region}_20c.txt   >> $ps
}

ps=fig_6_EL_memo.ps
gmt set PS_MEDIA=Custom_16.0cx13.5c
gmt set MAP_LABEL_OFFSET=3p
gmt set MAP_ANNOT_OFFSET=3p
gmt set MAP_FRAME_WIDTH=0.06c
gmt set MAP_TICK_LENGTH_PRIMARY=0.06c
gmt set MAP_GRID_PEN_PRIMARY=default,lightgrey
gmt set MAP_TITLE_OFFSET=0.1c
gmt set MAP_FRAME_PEN=0.5p,40/40/40
gmt set FONT_ANNOT_PRIMARY             = 8p,Hind-Light,40/40/40
gmt set FONT_ANNOT_SECONDARY           = 8p,Hind-Light,40/40/40
gmt set FONT_LABEL                     = 8p,Hind-Light,40/40/40
gmt set FONT_LOGO                      = 8p,Hind-Light,40/40/40
gmt set FONT_TITLE                     = 8p,Hind-Medium,40/40/40
gmt set PS_PAGE_ORIENTATION=portrait
gmt set PS_LINE_CAP=butt
gmt set PS_LINE_JOIN=round

R=1970/2050/-0.12/0.54
J=X7c/4.5c
Bx=x20g20
By=y0.1g0.1
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
region="GMSL"
gmt psbasemap -K -R$R -J$J -Y8.5c -X1.0c -BWeSn+t'Global-mean sea level' -B$Bx -B$By+l'Sea level (m)'  > $ps
print_lines

region="USA"
gmt psbasemap -O -K -R$R -J$J -X7.8c -BweSn+t'Contiguous United States' -B$Bx -B$By  >> $ps
print_lines

region="NE"
gmt psbasemap -O -K -R$R -J$J -Y-5.8c -X-7.8c -BWeSn+t'Northeast Coast' -B$Bx -B$By+l'Sea level (m)'  >> $ps
print_lines

region="SW"
gmt psbasemap -O -K -R$R -J$J -X7.8c -BweSn+t'Southwest Coast' -B$Bx -B$By  >> $ps
print_lines


gmt psbasemap -O -K -R0/0.85/0.5/6.5 -JX6.8c/1.8c -X-3.8c -Y-2.55c -T100/100/1  >> $ps
echo "0 6.2 Observations:" | gmt pstext -R -J -F+fHind-Medium,8p+jLM -O -K -N >> $ps

echo -e "0.0 5 \n 0.05 5" | gmt psxy -R -J -O -K -W1.5p,$color0 >> $ps
echo "0.055 5 Tide-gauge observations" | gmt pstext -R -J -F+f8+jLM -O -K -N >> $ps

echo -e "0.0 4 \n 0.05 4" | gmt psxy -R -J -O -K -W2p,$color5 >> $ps
echo "0.055 4 Present trajectory" | gmt pstext -R -J -F+f8+jLM -O -K -N >> $ps

echo -e "0.5 5 \n 0.55 5" | gmt psxy -R -J -O -K -W1p,$colorlow >> $ps
echo "0.555 5 Low" | gmt pstext -R -J -F+f8+jLM -O -K -N >> $ps

echo "0.5 6.2 Scenarios:" | gmt pstext -R -J -F+fHind-Medium,8p+jLM -O -K -N >> $ps
echo -e "0.5 4 \n 0.55 4" | gmt psxy -R -J -O -K -W1p,$colorintlow >> $ps
echo "0.555 4 Intermediate-low" | gmt pstext -R -J -F+f8+jLM -O -K -N >> $ps

echo -e "0.5 3 \n 0.55 3" | gmt psxy -R -J -O -K -W1p,$colorint >> $ps
echo "0.555 3 Intermediate" | gmt pstext -R -J -F+f8+jLM -O -K -N >> $ps

echo -e "0.5 2 \n 0.55 2" | gmt psxy -R -J -O -K -W1p,$colorinthigh >> $ps
echo "0.555 2 Intermediate-high" | gmt pstext -R -J -F+f8+jLM -O -K -N >> $ps

echo -e "0.5 1 \n 0.55 1" | gmt psxy -R -J -O -K -W1p,$colorhigh >> $ps
echo "0.555 1 High" | gmt pstext -R -J -F+f8+jLM -O -N >> $ps

gmt psconvert -Tf -C-sFONTPATH=/home/thomas/.local/share/fonts/ $ps
gmt psconvert -Tg -C-sFONTPATH=/home/thomas/.local/share/fonts/ $ps


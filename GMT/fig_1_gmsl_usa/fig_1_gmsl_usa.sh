ps=fig_1_gmsl_usa.ps
gmt set PS_MEDIA=Custom_16.0cx7c
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
gmt set FONT_TITLE                     = 8p,Hind-Light,40/40/40
gmt set PS_PAGE_ORIENTATION=portrait
gmt set PS_LINE_CAP=butt
gmt set PS_LINE_JOIN=round

R=1970/2050/-0.1/0.7
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

gmt psbasemap -K -R$R -J$J -Y2.1c -X1.0c -BWeSn+t'Global average' -B$Bx -B$By+l'Sea level (m)'  > $ps
gmt psxy  -R -J -L+bD -t50 -N -G$color0 -O -K GMSL_20c.txt >> $ps
gmt psxy  -R -J -L+bD -t50 -N -G$color5 -O -K GMSL_Trajectory.txt >> $ps
gmt psxy  -R -J -L+bD -t80 -N -G$colorlow -O -K GMSL_Low.txt   >> $ps
gmt psxy  -R -J -L+bD -t80 -N -G$colorintlow -O -K GMSL_IntLow.txt   >> $ps
gmt psxy  -R -J -L+bD -t80 -N -G$colorint -O -K GMSL_Int.txt   >> $ps
gmt psxy  -R -J -L+bD -t80 -N -G$colorinthigh -O -K GMSL_IntHigh.txt   >> $ps
gmt psxy  -R -J -L+bD -t80 -N -G$colorhigh -O -K GMSL_High.txt   >> $ps

gmt psxy  -R -J -N -W1.0p,$colorlow -O -K GMSL_Low.txt   >> $ps
gmt psxy  -R -J -N -W1.0p,$colorintlow -O -K GMSL_IntLow.txt   >> $ps
gmt psxy  -R -J -N -W1.0p,$colorint -O -K GMSL_Int.txt   >> $ps
gmt psxy  -R -J -N -W1.0p,$colorinthigh -O -K GMSL_IntHigh.txt   >> $ps
gmt psxy  -R -J -N -W1.0p,$colorhigh -O -K GMSL_High.txt   >> $ps
gmt psxy  -R -J -N -W1.5p,$color5 -O -K GMSL_Trajectory.txt   >> $ps
gmt psxy  -R -J -N -W1.5p,$color0 -O -K GMSL_20c.txt   >> $ps
gmt psxy  -R -J -N -W1.5p,$color1 -O -K GMSL_Altimetry.txt   >> $ps

gmt psbasemap -O -K -R$R -J$J -X7.8c -BweSn+t'Contiguous United States' -B$Bx -B$By  >> $ps
gmt psxy  -R -J -L+bD -t50 -N -G$color5 -O -K USA_Trajectory.txt >> $ps
gmt psxy  -R -J -L+bD -t80 -N -G$colorlow -O -K USA_Low.txt   >> $ps
gmt psxy  -R -J -L+bD -t80 -N -G$colorintlow -O -K USA_IntLow.txt   >> $ps
gmt psxy  -R -J -L+bD -t80 -N -G$colorint -O -K USA_Int.txt   >> $ps
gmt psxy  -R -J -L+bD -t80 -N -G$colorinthigh -O -K USA_IntHigh.txt   >> $ps
gmt psxy  -R -J -L+bD -t80 -N -G$colorhigh -O -K USA_High.txt   >> $ps

gmt psxy  -R -J -N -W1.0p,$colorlow -O -K USA_Low.txt   >> $ps
gmt psxy  -R -J -N -W1.0p,$colorintlow -O -K USA_IntLow.txt   >> $ps
gmt psxy  -R -J -N -W1.0p,$colorint -O -K USA_Int.txt   >> $ps
gmt psxy  -R -J -N -W1.0p,$colorinthigh -O -K USA_IntHigh.txt   >> $ps
gmt psxy  -R -J -N -W1.0p,$colorhigh -O -K USA_High.txt   >> $ps
gmt psxy  -R -J -N -W1.5p,$color5 -O -K USA_Trajectory.txt   >> $ps
gmt psxy  -R -J -N -W1.5p,$color0 -O -K USA_20c.txt   >> $ps
gmt psxy  -R -J -N -W1.5p,$color1 -O -K USA_Altimetry.txt   >> $ps

gmt psbasemap -O -K -R0/0.85/0.5/5.5 -JX6.8c/1.5c -X-3.8c -Y-2.05c -T100/100/1  >> $ps
echo -e "0.0 5 \n 0.05 5" | gmt psxy -R -J -O -K -W1.5p,$color1 >> $ps
echo "0.055 5 Altimetry observations " | gmt pstext -R -J -F+f8+jLM -O -K -N >> $ps

echo -e "0.0 4 \n 0.05 4" | gmt psxy -R -J -O -K -t60 -W6p,$color0 >> $ps
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



# This script creates Figure 2.2
ps=fig_7_US.ps
gmt set PS_MEDIA=Custom_8.1cx5.3c
gmt set MAP_LABEL_OFFSET=4p
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

R=2020/2150/-0.1/4.5
J=X7c/4.5c
Bx=x25g25
By=y0.5g0.5

colorlow=#330a5f
colorintlow=#781c6d
colorint=#bb3755
colorinthigh=#ed6925
colorhigh=#fcb519

plot_unc="false"

gmt psbasemap -K -R$R -J$J -Y0.4c -X1.0c -BWeSn+t'Contiguous United States' -B$Bx -B$By+l'Sea level (m)'  > $ps

if [ "$plot_unc" = true ] ; then 
gmt psxy  -R -J -L+bD -t80 -N -G$colorlow -O -K USA_Low.txt   >> $ps
gmt psxy  -R -J -L+bD -t80 -N -G$colorintlow -O -K USA_IntLow.txt   >> $ps
gmt psxy  -R -J -L+bD -t80 -N -G$colorint -O -K USA_Int.txt   >> $ps
gmt psxy  -R -J -L+bD -t80 -N -G$colorinthigh -O -K USA_IntHigh.txt   >> $ps
gmt psxy  -R -J -L+bD -t80 -N -G$colorhigh -O -K USA_High.txt   >> $ps
fi

gmt psxy  -R -J -N -W1.0p,$colorlow -O -K USA_Low.txt   >> $ps
gmt psxy  -R -J -N -W1.0p,$colorintlow -O -K USA_IntLow.txt   >> $ps
gmt psxy  -R -J -N -W1.0p,$colorint -O -K USA_Int.txt   >> $ps
gmt psxy  -R -J -N -W1.0p,$colorinthigh -O -K USA_IntHigh.txt   >> $ps
gmt psxy  -R -J -N -W1.0p,$colorhigh -O -K USA_High.txt   >> $ps

echo "2024 4.3 Scenario:" | gmt pstext -R -J -F+f8,Hind-SemiBold+jLM -O -K -N >> $ps

if [ "$plot_unc" = true ] ; then 
echo -e "2024 3.95 \n 2034 3.95" | gmt psxy -R -J -O -K -t90 -W6p,$colorlow >> $ps
echo -e "2024 3.63 \n 2034 3.63" | gmt psxy -R -J -O -K -t90 -W6p,$colorintlow >> $ps
echo -e "2024 3.31 \n 2034 3.31" | gmt psxy -R -J -O -K -t90 -W6p,$colorint >> $ps
echo -e "2024 2.99 \n 2034 2.99" | gmt psxy -R -J -O -K -t90 -W6p,$colorinthigh >> $ps
echo -e "2024 2.67 \n 2034 2.67" | gmt psxy -R -J -O -K -t90 -W6p,$colorhigh >> $ps
fi

echo -e "2024 3.95 \n 2034 3.95" | gmt psxy -R -J -O -K -W1.5p,$colorlow >> $ps
echo -e "2024 3.63 \n 2034 3.63" | gmt psxy -R -J -O -K -W1.5p,$colorintlow >> $ps
echo -e "2024 3.31 \n 2034 3.31" | gmt psxy -R -J -O -K -W1.5p,$colorint >> $ps
echo -e "2024 2.99 \n 2034 2.99" | gmt psxy -R -J -O -K -W1.5p,$colorinthigh >> $ps
echo -e "2024 2.67 \n 2034 2.67" | gmt psxy -R -J -O -K -W1.5p,$colorhigh >> $ps

echo "2035 3.95 Low" | gmt pstext -R -J -F+f8+jLM -O -K -N >> $ps
echo "2035 3.63 Intermediate-low" | gmt pstext -R -J -F+f8+jLM -O -K -N >> $ps
echo "2035 3.31 Intermediate" | gmt pstext -R -J -F+f8+jLM -O -K -N >> $ps
echo "2035 2.99 Intermediate-high" | gmt pstext -R -J -F+f8+jLM -O -K -N >> $ps
echo "2035 2.67 High" | gmt pstext -R -J -F+f8+jLM -O -N >> $ps

# gmt psbasemap -O -K -R0/0.45/0.5/5.5 -JX3c/1.5c -X2c -Y-2.05c -T100/100/1  >> $ps

gmt psconvert -Tf -C-sFONTPATH=/home/thomas/.local/share/fonts/ $ps
gmt psconvert -Tg -C-sFONTPATH=/home/thomas/.local/share/fonts/ $ps



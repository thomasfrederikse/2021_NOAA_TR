ps=fig_4_gmsl.ps
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

R=1920/2018/-0.22/0.12
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


gmt psbasemap -K -R$R -J$J -Y2.1c -X1.0c -BWeSn+t'Causes of global-mean sea-level rise' -B$Bx -B$By+l'Sea level (m)'  > $ps
gmt psxy  -R -J -L+bD -t80 -N -G$color3 -O -K glac.txt >> $ps
gmt psxy  -R -J -L+bD -t80 -N -G$color4 -O -K GrIS.txt >> $ps
gmt psxy  -R -J -L+bD -t80 -N -G$color5 -O -K AIS.txt   >> $ps
gmt psxy  -R -J -L+bD -t80 -N -G$color6 -O -K tws.txt   >> $ps
gmt psxy  -R -J -L+bD -t80 -N -G$color7 -O -K steric.txt   >> $ps
gmt psxy  -R -J -L+bD -t60 -N -G$color2 -O -K Budget.txt   >> $ps
gmt psxy  -R -J -L+bD -t60 -N -G$color1 -O -K GMSL.txt   >> $ps

gmt psxy  -R -J -N -W1.5p,$color3 -O -K glac.txt   >> $ps
gmt psxy  -R -J -N -W1.5p,$color4 -O -K GrIS.txt   >> $ps
gmt psxy  -R -J -N -W1.5p,$color5 -O -K AIS.txt   >> $ps
gmt psxy  -R -J -N -W1.5p,$color6 -O -K tws.txt   >> $ps
gmt psxy  -R -J -N -W1.5p,$color7 -O -K steric.txt   >> $ps
gmt psxy  -R -J -N -W2p,$color2 -O -K Budget.txt   >> $ps
gmt psxy  -R -J -N -W2p,$color1 -O -K GMSL.txt   >> $ps
echo "1920 0.12 a)" | gmt pstext -D0.125c/-0.125c -Gwhite -R -J -F+f9p,Hind-SemiBold+jLT -O -K >> $ps

gmt psbasemap -O -K -R$R -J$J -X7.8c -BweSn+t'Contiguous United States versus global-mean sea level' -B$Bx -B$By  >> $ps
gmt psxy  -R -J -N -W2p,$color1 -O -K GMSL.txt   >> $ps
gmt psxy  -R -J -N -W2p,$color0 -O -K USA_MSL.txt   >> $ps
echo "1920 0.12 b)" | gmt pstext -D0.125c/-0.125c -Gwhite -R -J -F+f9p,Hind-SemiBold+jLT -O -K >> $ps

gmt psbasemap -O -K -R0/0.85/0.5/5.5 -JX6.8c/1.5c -X-3.8c -Y-2.05c -T100/100/1  >> $ps
echo -e "0.0 5 \n 0.05 5" | gmt psxy -R -J -O -K -t60 -W6p,$color3 >> $ps
echo -e "0.0 5 \n 0.05 5" | gmt psxy -R -J -O -K -W1.5p,$color3 >> $ps
echo "0.055 5 Glaciers " | gmt pstext -R -J -F+f8+jLM -O -K -N >> $ps

echo -e "0.0 4 \n 0.05 4" | gmt psxy -R -J -O -K -t60 -W6p,$color4 >> $ps
echo -e "0.0 4 \n 0.05 4" | gmt psxy -R -J -O -K -W1.5p,$color4 >> $ps
echo "0.055 4 Greenland Ice Sheet" | gmt pstext -R -J -F+f8+jLM -O -K -N >> $ps

echo -e "0.0 3 \n 0.05 3" | gmt psxy -R -J -O -K -t60 -W6p,$color5 >> $ps
echo -e "0.0 3 \n 0.05 3" | gmt psxy -R -J -O -K -W1.5p,$color5 >> $ps
echo "0.055 3 Antarctic Ice Sheet" | gmt pstext -R -J -F+f8+jLM -O -K -N >> $ps

echo -e "0.0 2 \n 0.05 2" | gmt psxy -R -J -O -K -t60 -W6p,$color6 >> $ps
echo -e "0.0 2 \n 0.05 2" | gmt psxy -R -J -O -K -W1.5p,$color6 >> $ps
echo "0.055 2 Terrestrial Water Storage" | gmt pstext -R -J -F+f8+jLM -O -K -N >> $ps

echo -e "0.0 1 \n 0.05 1" | gmt psxy -R -J -O -K -t60 -W6p,$color7 >> $ps
echo -e "0.0 1 \n 0.05 1" | gmt psxy -R -J -O -K -W1.5p,$color7 >> $ps
echo "0.055 1 Thermosteric expansion" | gmt pstext -R -J -F+f8+jLM -O -K -N >> $ps

echo -e "0.5 5 \n 0.55 5" | gmt psxy -R -J -O -K -t60 -W6p,$color2 >> $ps
echo -e "0.5 5 \n 0.55 5" | gmt psxy -R -J -O -K -W1.5p,$color2 >> $ps
echo "0.555 5 Sum of contributing processes" | gmt pstext -R -J -F+f8+jLM -O -K -N >> $ps

echo -e "0.5 4 \n 0.55 4" | gmt psxy -R -J -O -K -t60 -W6p,$color1 >> $ps
echo -e "0.5 4 \n 0.55 4" | gmt psxy -R -J -O -K -W1.5p,$color1 >> $ps
echo "0.555 4 Observed global-mean sea level" | gmt pstext -R -J -F+f8+jLM -O -K -N >> $ps

echo -e "0.5 3 \n 0.55 3" | gmt psxy -R -J -O -K -W1.5p,$color0 >> $ps
echo "0.555 3 Sea level along contiguous US" | gmt pstext -R -J -F+f8+jLM -O -N >> $ps


gmt psconvert -Tf -C-sFONTPATH=/home/thomas/.local/share/fonts/ $ps
gmt psconvert -Tg -C-sFONTPATH=/home/thomas/.local/share/fonts/ $ps



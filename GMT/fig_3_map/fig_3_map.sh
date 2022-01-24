# Script to generate map with TG locations and regions (in appendix)
ps=fig_3_map.ps
gmt set PS_MEDIA=Custom_12.0cx10.3c
gmt set MAP_LABEL_OFFSET=3p
gmt set MAP_ANNOT_OFFSET=3p
gmt set MAP_FRAME_WIDTH=0.06c
gmt set MAP_TICK_LENGTH_PRIMARY=0.06c
gmt set MAP_FRAME_PEN=0.5p,40/40/40
gmt set FONT_ANNOT_PRIMARY             = 8p,Hind-Light,40/40/40
gmt set FONT_ANNOT_SECONDARY           = 8p,Hind-Light,40/40/40
gmt set FONT_LABEL                     = 8p,Hind-Light,40/40/40
gmt set FONT_LOGO                      = 8p,Hind-Light,40/40/40
gmt set FONT_TITLE                     = 8p,Hind-Light,40/40/40
gmt set PS_PAGE_ORIENTATION=portrait
gmt set MAP_FRAME_TYPE=plain
gmt set COLOR_NAN=200/200/200
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
color10=#439894
gmt makecpt -Cwhite,$color1,$color2,$color3,$color4,$color5,$color6,$color7,$color8,$color9,$color10 -W >regioncolor.cpt
gmt grdmath ../../Data/region_mask.nc?mask_num 0 NAN = region_mask.nc
gmt psbasemap -K -R170/310/5/75 -JJ180/11c -BWeSn -X0.7c -Y3.2c -By20f10 -Bx30f15 > $ps
gmt grdimage -R -J -O -K region_mask.nc?mask_num -Cregioncolor.cpt -nn -t40 >> $ps
gmt pscoast -R -J -Dl -A10/0/1 -G120/120/120 -O -K >> $ps
gmt psxy -R -J -O -K EC.txt -Sc0.15c -G$color1 >> $ps
gmt psxy -R -J -O -K SE.txt -Sc0.15c -G$color2 >> $ps
gmt psxy -R -J -O -K GCE.txt -Sc0.15c -G$color3 >> $ps
gmt psxy -R -J -O -K GCW.txt -Sc0.15c -G$color4 >> $ps
gmt psxy -R -J -O -K SWC.txt -Sc0.15c -G$color5 >> $ps
gmt psxy -R -J -O -K NWC.txt -Sc0.15c -G$color6 >> $ps
gmt psxy -R -J -O -K PAC.txt -Sc0.15c -G$color7 >> $ps
gmt psxy -R -J -O -K CAR.txt -Sc0.15c -G$color8 >> $ps
gmt psxy -R -J -O -K ALN.txt -Sc0.15c -G$color9 >> $ps
gmt psxy -R -J -O -K ALS.txt -Sc0.15c -G$color10 >> $ps

gmt psbasemap -R0/2/0.5/5.5 -JX7c/2.5c -O -K -Y-3.1c -X2c -T100/100/1 >> $ps
echo -e "0 5 \n 0.2 5" | gmt psxy -R -J -O -K -t40 -W7p,$color1 >> $ps
echo -e "0.1 5 1" | gmt psxy -R -J -O -K -Sc0.15c -G$color1 >> $ps
echo "0.21 5 Northeast " | gmt pstext -R -J -F+f8+jLM -O -K -N >> $ps

echo -e "0 4 \n 0.2 4" | gmt psxy -R -J -O -K -t40 -W7p,$color2 >> $ps
echo -e "0.1 4 1" | gmt psxy -R -J -O -K -Sc0.15c -G$color2 >> $ps
echo "0.21 4 Southeast " | gmt pstext -R -J -F+f8+jLM -O -K -N >> $ps

echo -e "0 3 \n 0.2 3" | gmt psxy -R -J -O -K -t40 -W7p,$color3 >> $ps
echo -e "0.1 3 1" | gmt psxy -R -J -O -K -Sc0.15c -G$color3 >> $ps
echo "0.21 3 Eastern Gulf Coast" | gmt pstext -R -J -F+f8+jLM -O -K -N >> $ps

echo -e "0 2 \n 0.2 2" | gmt psxy -R -J -O -K -t40 -W7p,$color4 >> $ps
echo -e "0.1 2 1" | gmt psxy -R -J -O -K -Sc0.15c -G$color4 >> $ps
echo "0.21 2 Western Gulf Coast" | gmt pstext -R -J -F+f8+jLM -O -K -N >> $ps

echo -e "0 1 \n 0.2 1" | gmt psxy -R -J -O -K -t40 -W7p,$color5 >> $ps
echo -e "0.1 1 1" | gmt psxy -R -J -O -K -Sc0.15c -G$color5 >> $ps
echo "0.21 1 Southwest Coast " | gmt pstext -R -J -F+f8+jLM -O -K -N >> $ps
# # Line 2
echo -e "1 5 \n 1.2 5" | gmt psxy -R -J -O -K -t40 -W7p,$color6 >> $ps
echo -e "1.1 5 1" | gmt psxy -R -J -O -K -Sc0.15c -G$color6 >> $ps
echo "1.21 5 Northwest Coast " | gmt pstext -R -J -F+f8+jLM -O -K -N >> $ps

echo -e "1 4 \n 1.2 4" | gmt psxy -R -J -O -K -t40 -W7p,$color7 >> $ps
echo -e "1.1 4 1" | gmt psxy -R -J -O -K -Sc0.15c -G$color7 >> $ps
echo "1.21 4 Hawaiian Islands " | gmt pstext -R -J -F+f8+jLM -O -K -N >> $ps

echo -e "1 3 \n 1.2 3" | gmt psxy -R -J -O -K -t40 -W7p,$color8 >> $ps
echo -e "1.1 3 1" | gmt psxy -R -J -O -K -Sc0.15c -G$color8 >> $ps
echo "1.21 3 Caribbean Islands" | gmt pstext -R -J -F+f8+jLM -O -K -N >> $ps

echo -e "1 2 \n 1.2 2" | gmt psxy -R -J -O -K -t40 -W7p,$color9 >> $ps
echo -e "1.1 2 1" | gmt psxy -R -J -O -K -Sc0.15c -G$color9 >> $ps
echo "1.21 2 Alaska North" | gmt pstext -R -J -F+f8+jLM -O -K -N >> $ps

echo -e "1 1 \n 1.2 1" | gmt psxy -R -J -O -K -t40 -W7p,$color10 >> $ps
echo -e "1.1 1 1" | gmt psxy -R -J -O -K -Sc0.15c -G$color10 >> $ps
echo "1.21 1 Alaska South " | gmt pstext -R -J -F+f8+jLM -O -N >> $ps


gmt psconvert -Tf -C-sFONTPATH=/Users/tfrederi/Library/Fonts/ $ps
gmt psconvert -Tg -C-sFONTPATH=/Users/tfrederi/Library/Fonts/ $ps

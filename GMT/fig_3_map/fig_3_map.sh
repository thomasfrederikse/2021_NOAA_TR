ps=fig_3_map.ps
gmt set PS_MEDIA=Custom_8.0cx6.5c
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

gmt psbasemap -K -R170/310/5/75 -JJ180/7c -BWeSn -X0.7c -Y1c -By20f10 -Bx30f15 > $ps
# gmt grdimage -R -J -O -K ../../Data/region_mask.nc?mask_num -Cregioncolor.cpt -nn -t60 >> $ps
gmt pscoast -R -J -Dl -A5000/0/1 -G120/120/120 -O -K >> $ps
gmt psxy -R -J -O -K EC.txt -Sc0.15c -G$color1 >> $ps
gmt psxy -R -J -O -K SE.txt -Sc0.15c -G$color2 >> $ps
gmt psxy -R -J -O -K GCE.txt -Sc0.15c -G$color3 >> $ps
gmt psxy -R -J -O -K GCW.txt -Sc0.15c -G$color4 >> $ps
gmt psxy -R -J -O -K SWC.txt -Sc0.15c -G$color5 >> $ps
gmt psxy -R -J -O -K NWC.txt -Sc0.15c -G$color6 >> $ps
gmt psxy -R -J -O -K PAC.txt -Sc0.15c -G$color7 >> $ps
gmt psxy -R -J -O -K CAR.txt -Sc0.15c -G$color8 >> $ps
gmt psxy -R -J -O -K ALN.txt -Sc0.15c -G$color9 >> $ps
gmt psxy -R -J -O  ALS.txt -Sc0.15c -G$color10 >> $ps

gmt psconvert -Tf -C-sFONTPATH=/Users/tfrederi/Library/Fonts/ $ps
gmt psconvert -Tg -C-sFONTPATH=/Users/tfrederi/Library/Fonts/ $ps

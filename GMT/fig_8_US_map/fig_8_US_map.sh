# Script to generate map with TG locations and regions (in appendix)
ps=fig_8_US_map.ps
gmt set PS_MEDIA=Custom_12.0cx7.1c
gmt set MAP_LABEL_OFFSET=3p
gmt set MAP_TITLE_OFFSET=-0.2c
gmt set MAP_ANNOT_OFFSET=3p
gmt set MAP_FRAME_WIDTH=0.06c
gmt set MAP_TICK_LENGTH_PRIMARY=0.06c
gmt set MAP_FRAME_PEN=0.5p,40/40/40
gmt set FONT_ANNOT_PRIMARY             = 8p,Hind-Light,40/40/40
gmt set FONT_ANNOT_SECONDARY           = 8p,Hind-Light,40/40/40
gmt set FONT_LABEL                     = 8p,Hind-Light,40/40/40
gmt set FONT_LOGO                      = 8p,Hind-Light,40/40/40
gmt set FONT_TITLE                     = 8p,Hind-SemiBold,40/40/40
gmt set PS_PAGE_ORIENTATION=portrait
gmt set MAP_FRAME_TYPE=plain
gmt set COLOR_NAN=200/200/200

gmt makecpt -Clajolla -T2/6/41+n -M -D > rsl.cpt

gmt psbasemap -K -R220/310/16/51 -JM11c -Bwesn+t"Sea level change in 2100 for the Intermediate scenario" -X0.7c -Y1.2c -By10f10 -Bx30f15 > $ps
gmt grdimage -R -J -O -K fig_8_US_map.nc?rsl -Crsl.cpt -E400 -nl >> $ps
gmt pscoast -R -J -Di -A10/0/1 -G120/120/120 -EUS+p0.5p,80/80/80 -O -K >> $ps
gmt psscale -Dx1c/-0.75c+w9c/0.25ch+e -R -J  -By+l'ft' -Crsl.cpt -B0.5 -O -S >> $ps

gmt psconvert -Tg -C-sFONTPATH=/Users/tfrederi/Library/Fonts/ $ps
gmt psconvert -Tf -C-sFONTPATH=/Users/tfrederi/Library/Fonts/ $ps


#!/bin/sh

# converts the Inkscape SVG files to TIFF files formatted in the way PLOS want
convert="convert -font Arial -density 600 -format tiff -colorspace rgb -compress lzw -depth 8 -alpha off"
$convert Fig1_recap/recap.svg Fig1.tiff
$convert Fig2_simdiffpatts/simdiffpatts.svg Fig2.tiff
$convert Fig3_pattern/pattern.svg Fig3.tiff
$convert Fig4_elazorsi/elazorsi.svg Fig4.tiff
$convert Fig5_naturalstim/naturalstim.svg Fig5.tiff
$convert Fig6_avkernels/avkernels.svg Fig6.tiff


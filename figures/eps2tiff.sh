#!/bin/sh

# converts the EPS files to TIFF files formatted in the way PLOS want
mogrify -font Arial -density 600 -format tiff -colorspace rgb -compress lzw -depth 8 -alpha off *.eps

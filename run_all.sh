#!/bin/bash
CURDIR=`pwd`

$CURDIR/deblur -m 60 -n 60 -p 60 -i denoised_array_3D.bin -o deblur.out -b 0
$CURDIR/deblur -m 60 -n 60 -p 60 -i denoised_array_3D.bin -o deblur.out -b 1
$CURDIR/deblur -m 60 -n 60 -p 60 -i denoised_array_3D.bin -o deblur.out -b 2
$CURDIR/deblur -m 60 -n 60 -p 60 -i denoised_array_3D.bin -o deblur.out -b 3
$CURDIR/deblur -m 60 -n 60 -p 60 -i denoised_array_3D.bin -o deblur.out -b 4
$CURDIR/deblur -m 60 -n 60 -p 60 -i denoised_array_3D.bin -o deblur.out -b 5
$CURDIR/deblur -m 60 -n 60 -p 60 -i denoised_array_3D.bin -o deblur.out -b 6
$CURDIR/deblur -m 60 -n 60 -p 60 -i denoised_array_3D.bin -o deblur.out -b 7
$CURDIR/deblur -m 60 -n 60 -p 60 -i denoised_array_3D.bin -o deblur.out -b 8
$CURDIR/deblur -m 60 -n 60 -p 60 -i denoised_array_3D.bin -o deblur.out -b 9

#! /bin/csh

set spec="1H1D"
set date="030415"
set sample="1mA16_dsA6RNA_WajS"
set expt="Ann955_pH6.8"
set temp="25C"

mv out.ft2 $sample"_"$spec"_"$expt"_"$temp"_"$date".ft2"

#!/bin/bash

message() {
    echo "Usage: $0 [algorithm]"
    echo "[algorithm] must be 1, 2, or 3:"
    echo "   1  to run Accelerated DF-SANE in R"
    echo "   2  to run Accelerated DF-SANE in Fortran"
    echo "   3  to run NITSOL"
    exit 1
}

if [ -z $1 ]; then
    message
fi

if [ "$1" != "1" -a "$1" != "2" -a  "$1" != "3" ]; then
    message
fi

problems=("BOOTH" "CLUSTER" "CUBENE" "DENSCHNCNE" "DENSCHNFNE"
	  "FREURONE" "GOTTFR" "HIMMELBA" "HIMMELBC" "HIMMELBD" "HS8"
	  "HYPCIR" "POWELLBS" "POWELLSQ" "PRICE3NE" "PRICE4NE"
	  "RSNBRNE" "SINVALNE" "WAYSEA1NE" "WAYSEA2NE" "DENSCHNDNE"
	  "DENSCHNENE" "HATFLDF" "HATFLDFLNE" "HELIXNE" "HIMMELBE"
	  "RECIPE" "ZANGWIL3" "POWELLSE" "POWERSUMNE" "HEART6"
	  "HEART8" "COOLHANS" "MOREBVNE" "OSCIPANE" "TRIGON1NE"
	  "INTEQNE" "HATFLDG" "HYDCAR6" "METHANB8" "METHANL8"
	  "HYDCAR20" "LUKSAN21" "MANCINONE" "QINGNE" "ARGTRIG"
	  "BROWNALE" "CHANDHEU" "10FOLDTR" "KSS" "MSQRTA" "MSQRTB"
	  "EIGENAU" "EIGENB" "EIGENC" "NONMSQRTNE" "BROYDN3D"
	  "BROYDNBD" "BRYBNDNE" "NONDIANE" "SBRYBNDNE" "SROSENBRNE"
	  "SSBRYBNDNE" "TQUARTICNE" "OSCIGRNE" "CYCLIC3" "YATP1CNE"
	  "YATP1NE" "YATP2CNE" "YATP2SQ")

for prob in ${problems[@]}; do
    if [ "$1" == "1" ]; then
	echo "$prob" | Rscript dfsaneacc-forcutest.R
	cat tabline.txt >> table-dfsaneacc-r-cutest.txt
	rm tabline.txt
    elif [ "$1" == "2" ]; then
	./run-dfsaneaccma-forcutest.sh $prob
    elif [ "$1" == "3" ]; then
	./run-nitsolma-forcutest.sh $prob
    fi
done

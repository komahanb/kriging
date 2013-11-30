	#!/bin/bash		       
	# This script export the tecplot layout into eps and then to PDF

	if [ -f list ];then
	    rm list
	fi
	counter=0
	for i in $( ls *.lay ); do

	    let counter=counter+1

	    tec360 -mesa -b $i export.mcr  #opens a layout and then exports

	    NAME=`echo "$i" | cut -d'.' -f1` #get a filename to the figure based on the layout file
	    mv temp.eps $NAME.eps
	    epstopdf $NAME.eps
#	    rm $NAME.eps
	    if [$counter -eq 1];then
	    echo "$NAME.pdf"  > list
	    else
		echo "$NAME.pdf"  >> list
	    fi
	done

	if [ -f batch.log ];then
	    rm batch.log
	fi
	cat list
	

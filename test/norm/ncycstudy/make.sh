	#!/bin/bash		       
	# This script export the all tecplot layouts (layout.lay) in the folder into eps and then to PDF

#INPUT:
# 1) layout file (*.lay)
# 2) tecplot export macro (export.mcr) in the same directory
	
#OUTPUT:
#Example: kriging.lay is exported in to kriging.eps and then to kriging.pdf
	
#tec360 -mesa -b layout.lay export.mcr 
#b stands for tecplot batch mode (command line export)
	

	for i in $( ls *.lay ); do

	    let counter=counter+1

	    tec360 -mesa -b $i export.mcr  #opens a layout and then exports as eps

	    NAME=`echo "$i" | cut -d'.' -f1` #get a filename to the figure based on the layout file

	    mv temp.eps $NAME.eps

	    epstopdf $NAME.eps

	    echo "====================================="
	    echo "$NAME.eps created succesfully"
	    echo "====================================="

	done
	

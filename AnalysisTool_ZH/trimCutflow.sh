#!/bin/bash
set -e
# this script needs to be run with "./" instead of "source" or else errors will close the ssh session

INFILE=$1
CUTFLOWFILE=$INFILE\_cutflow.txt
if [ -f "$CUTFLOWFILE" ]; then
    rm $CUTFLOWFILE
fi

NUMLINES=$(wc -l < $INFILE)
IFS=''
rswitch=false
k=1
echo "0% done"
while read line ; do
    if [[ "$line" =~ "===========" ]]; then rswitch=true
    fi
    if $rswitch ; then echo "$line" >> $CUTFLOWFILE
    fi
    if [[ "$k % 100000" -eq "0" ]]; then echo "$(($k*100/$NUMLINES))% done"
    fi
    ((++k))
done < $INFILE
echo "100% done"

#j=1
#while read line ; do
#    ((++j))
#    if [[ "$j" -le "$k" ]]; then continue
#    fi
#    echo "$line" >> $CUTFLOWFILE
#done < $INFILE
echo "$INFILE" >> $CUTFLOWFILE

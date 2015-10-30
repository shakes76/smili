#!/bin/sh
 appname=`basename $0 | sed s,\.sh$,,`

 dirname=`dirname $0`
 tmp="${dirname#?}"
 echo "Using SMILI from:"
 echo $dirname

 if [ "${dirname%$tmp}" != "/" ]; then
 dirname=$PWD/$dirname
 fi
 LD_LIBRARY_PATH="${dirname}/lib:${LD_LIBRARY_PATH}"
 export LD_LIBRARY_PATH
 $dirname/bin/$appname $*
 
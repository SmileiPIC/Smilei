#!/bin/sh

set -u
set -e

H=$(pwd)
SCRIPTPATH=$( cd "$(dirname "$0")" ; pwd -P )

function error_exit
{
    echo " ${1:-"Unknown Error"}" 1>&2
    cd $H
    exit 1
}

function title {
    echo "::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"
    echo $1
    echo "::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"
}

function runTest {
    cd $SCRIPTPATH/../..
    for F in "$@"
    do
      D=`basename $F .py`
      rm -fr tests/collisions/$D
      ./smilei.sh 4 tests/collisions/$F || error_exit "Failed in $F"
    done
    cd $H
}

title "Starting Maxwellinization"
runTest Maxwellianization1.py
title "Ending Maxwellinization"


title "Starting StoppingPower"
runTest `ls Stopping_power?.py`
title "Ending StoppingPower"


title "Starting beam_relaxation"
runTest `ls beam_relaxation?.py`
title "Ending beam_relaxation"


title "Starting temperature_isotropization"
runTest temperature_isotropization1.py
title "Ending temperature_isotropization"


title "Starting thermalisation_ei"
runTest `ls thermalisation_ei?.py`


title "Starting conductivity"
runTest `ls conductivity?.py`
title "Ending conductivity"


title "Starting ionization_equilibrium"
runTest `ls ionization_equilibrium[A-Z]*.py`
title "Ending ionization_equilibrium"


title "Starting ionization_multiple"
runTest `ls ionization_multiple[A-Z]*.py`
title "Ending ionization_multiple"


title "Starting ionization_rate"
runTest `ls ionization_rate?.py`
title "Ending ionization_rate"


title "Starting ionization_stopping_power"
runTest `ls ionization_stopping_power?.py`
title "Ending ionization_stopping_power"

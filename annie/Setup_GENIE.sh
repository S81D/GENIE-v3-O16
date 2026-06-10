#! /usr/bin/env bash

#Application Paths

Gapp=/Genie
BUILD_PATH=/exp/annie/app/users/doran/GENIE_production/GENIE-v3-O16
# ^ this can be replaced with the usual path: ${Gapp}/GENIE-v3-master

export LIBGL_ALWAYS_INDIRECT=1

export DISPLAY=:0

export ROOTSYS=${Gapp}/root-6.24.06/install/

export GENIE=${BUILD_PATH}

export LD_LIBRARY_PATH=/lib:.:${Gapp}/log4cpp/install/lib:${Gapp}/Pythia6Support/v6_424/lib:${ROOTSYS}/lib:${Gapp}/LHAPDF-6.3.0/install/lib:${GENIE}/install/lib:${LD_LIBRARY_PATH}

export PYTHIA6_DIR=${Gapp}/Pythia6Support/v6_424/
export PYTHIA6_INCLUDE_DIR=${Gapp}/Pythia6Support/v6_424/inc/
export PYTHIA6_LIBRARY=${Gapp}/Pythia6Support/v6_424/lib/

export LHAPATH=${Gapp}/LHAPDF-6.3.0/install/share/LHAPDF:${LHAPATH}

export PATH=.:${Gapp}/fsplit/:${ROOTSYS}/bin:${Gapp}/LHAPDF-6.3.0/install/bin:${GENIE}/install/bin:${PATH}

export MANPATH=${ROOTSYS}/bin:${MANPATH}

export PYTHONPATH=${ROOTSYS}/lib:${PYTHONPATH}

#cannot include file in GitHub, too large
#set the path via the wrapper script
#export GENIEXSECFILE=${GENIE}/annie/gxspl-FNALsmall.xml

#!/usr/bin/env bash

# if building from a modified GENIE repo

source /exp/annie/app/users/doran/GENIE_production/Setup_GENIE.sh
cd $GENIE

./configure \
  --prefix=${GENIE}/install \
  --enable-lhapdf6 \
  --disable-lhapdf5 \
  --enable-fnal \
  --enable-flux-drivers \
  --enable-geom-drivers \
  --enable-dylibversion \
  --with-optimiz-level=O2 \
  --with-compiler=gcc \
  --with-pythia6-lib=/Genie/Pythia6Support/v6_424/lib/ \
  --with-lhapdf6-lib=/Genie/LHAPDF-6.3.0/install/lib \
  --with-lhapdf6-inc=/Genie/LHAPDF-6.3.0/install/include \
  --with-libxml2-inc=/usr/include/libxml2 \
  --with-libxml2-lib=/usr/lib64 \
  --with-log4cpp-inc=/Genie/log4cpp/install/include/ \
  --with-log4cpp-lib=/Genie/log4cpp/install/lib/

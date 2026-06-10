#!/usr/bin/env bash
source /exp/annie/app/users/doran/GENIE_production/Setup_GENIE.sh

RUN=0
SEED=42

chmod +x run_annie_genie.sh
./run_annie_genie.sh \
  -o . \
  -r ${RUN} \
  -n 20000 \
  -f gsimple_flux_${RUN}.root \
  -x gxspl-FNALsmall.xml \
  --seed ${SEED} \
  --topvol WORLD_LV \
  --tune G18_10a_02_11a

echo "Your wish is granted"

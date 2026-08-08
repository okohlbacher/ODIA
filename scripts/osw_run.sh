#!/bin/bash
# $1=tag  $2=mzML basename  $3=osw-format library  $4=threads
source /ceph/ibmi/abi/oliver/AI/OpenDIAlyzer/ODIA/scripts/env.sh
export OPENMS_TMPDIR=/dev/shm/$USER-osw; mkdir -p $OPENMS_TMPDIR
R=/ceph/ibmi/abi/oliver/AI/OpenDIAlyzer
D=/scratch/kohlbach/odia/data; O=/scratch/kohlbach/odia/out
tag=$1; mz=$2; lib=$3; th=${4:-32}
[ -f $D/$mz ] || cp $R/data/$mz $D/
/usr/bin/time -v -o $O/$tag.time \
  $R/opt/openms-3.6.0/bin/OpenSwathWorkflow \
  -in $D/$mz -tr $R/shared/lib/$lib \
  -out_features $O/$tag.osw -threads $th \
  > $O/$tag.log 2>&1
echo "exit=$?" >> $O/$tag.log
echo "$tag FINISHED"

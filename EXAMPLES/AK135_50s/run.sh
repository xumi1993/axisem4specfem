#!/bin/bash

# Create Diags directory
mkdir -p Diags

datapath=`grep "^DATA_DIR" inparam_advanced  |awk '{print $2}'| sed 's/\"//g'`
infopath=`grep "^INFO_DIR" inparam_advanced  |awk '{print $2}'| sed 's/\"//g'`
meshdir=./MESHES/
rm Mesh
rm -rf $meshdir ${datapath} ${infopath}
rm timestamp*.txt
rm simulation.info
rm OUTPUT_*
rm external_model.bm


# Run xmesh
../../bin/xmesh > OUTPUT_MESH

# move mesh files to MESHES & prepare for axisem run
mkdir -p $meshdir
mv meshdb.dat* $meshdir
cp inparam_mesh ak135.smooth.bm mesh_params.h  $meshdir
cp ak135.smooth.bm ./external_model.bm
ln -s $meshdir/ Mesh
mkdir -p $datapath
mkdir -p $infopath


# run on local machine
mpirun -n 4 ../../bin/axisem > OUTPUT_AXISEM

# post processing
cp ../param_post_processing ./
outdir=`grep "DATA_DIR" ./param_post_processing |awk '{print $2}' |sed 's/"/ /g'`
mkdir -p $outdir/SEISMOGRAMS
../../bin/xpost_processing > OUTPUT_POSTPROC
#!/bin/bash
set -e
change_par(){
  local param=$1
  local value=$2
  local file=$3
  
  local oldstr=`grep "^$param " $file`
  local newstr="$param        $value"
  sed -i "s?$oldstr?$newstr?g" $file
}

# set your study region lonmin/max latmin/max
lon0=114
lon1=132
lat0=41
lat1=46
maxdep=450

### STOP HERE ##########

cat CMT_DIR/injection_time |while read line;
do
  name=`echo $line | awk  '{print $1}'`
  t0=`echo $line |awk '{print $2-50.}'` 
  tend=`echo "$line" |awk '{print $2 + 300.}'`
  echo $name $t0 $tend

  # substitute dump_t0 and SEISMO_LENGTH
  change_par SEISMOGRAM_LENGTH $tend inparam_basic
  change_par DUMP_T0 $t0 inparam_advanced

  # set dump region
  evla=`grep ^latitude CMT_DIR/CMTSOLUTION_${name} |cut -d':' -f2`
  evlo=`grep ^longitude CMT_DIR/CMTSOLUTION_${name} |cut -d':' -f2`
  result=$(python << EOF
import numpy as np

def loc2degs(slati, sloni, rlati, rloni):
    deg2rad = np.pi / 180.0
    rad2deg = 180.0 / np.pi

    # Convert to radians
    slat = np.radians(slati)
    slon = np.radians(sloni)
    rlat = np.radians(rlati)
    rlon = np.radians(rloni)

    # Compute cosine of arc
    cosarc = np.sin(rlat) * np.sin(slat) + np.cos(rlat) * np.cos(slat) * np.cos(rlon - slon)
    cosarc = np.clip(cosarc, -1.0, 1.0)  # avoid numerical errors

    gcarc_rad = np.arccos(cosarc)
    gcarc = np.degrees(gcarc_rad)         # in degrees

    return gcarc

evla=$evla
evlo=$evlo
lon = np.linspace($lon0,$lon1,100)
lat = np.linspace($lat0,$lat1,100)
pairs = np.array([[xi, yi] for xi in lon for yi in lat])
degs = loc2degs(evla,evlo,pairs[:,1],pairs[:,0])
degmin = np.min(degs)
degmax = np.max(degs)
print("%f %f" %(degmin-10,degmax+10))
EOF
)
  delmin=`echo $result |awk '{print $1}'`
  delmax=`echo $result |awk '{print $2}'`
  echo $delmin $delmax
  Rmin=`echo $maxdep | awk '{print 6371-$1}'` 

  change_par KERNEL_COLAT_MIN $delmin inparam_advanced
  change_par KERNEL_COLAT_MAX $delmax inparam_advanced
  change_par KERNEL_RMIN $Rmin inparam_advanced
  change_par KERNEL_RMAX 7000 inparam_advanced

  # copy stations/cmtsolution
  \cp CMT_DIR/CMTSOLUTION_${name} CMTSOLUTION
  \cp CMT_DIR/STATIONS_${name} STATIONS

  ./submit.csh ak135.${name} -q slurm
   #exit 1
done 

#!/bin/sh

# Read in the arguements
WD=`pwd`
S_OPTIONS="-o out -e err"

NPROC=1
NODES=0
PPN=0
SERIAL=0


while [ $# -gt 0 ]
do
    case "$1" in
    (-*np) NPROC="$2"; shift;;
    (-*nodes) NODES="$2"; shift;;
    (-*ppn) PPN="$2"; shift;;
    (-*walltime) TIME="$2"; shift;;
    (-*serial) SERIAL=1; NPROC=1; NODES=1; PPN=8;;
    (-*) echo "$0: error - unrecognized option $1" 1>&2; exit 1;;
    (*)  break;;
    esac
    shift
done

#check for (exclusive) serial jobs
if [ $SERIAL -eq 0 ]
then
  if [ $NPROC -gt 1 ]
  then
    #in absense of any other information, go for 8 proc/node
    NODES=$(($NPROC/16))
    if [ $(($NPROC % 16)) -ne 0 ]; then NODES=$(($NODES+1)); fi
    PPN=$((NPROC/$NODES))
  elif [ $NODES -gt 0 ]
  then
    #if nodes and/or ppn are specified, use that
    if [ $PPN -eq 0 ]
    then
      PPN=16
    fi
    NPROC=$(($NODES*$PPN))
  else
    #default to (non-exclusive) serial
    NODES=1
    PPN=1
  fi
fi
RESOURCES="-N $NODES -n $NPROC -t $TIME"

while [ "$1" ]
do

INPUTFILE=$1
NAME=${INPUTFILE}

# Create the input file
TARGET=run.$$.sh

cat >> $TARGET << EOF
#!/bin/bash
#SBATCH --export=ALL
#SBATCH --job-name=$NAME
#SBATCH -N $NODES -n $NPROC -t $TIME --mem 5000mb
#SBATCH --tasks-per-node=$PPN

cd  \$SLURM_SUBMIT_DIR
./optimize.py $NAME

sleep 3
 
EOF

# Submit the input file to PBS and clean up
sbatch $S_OPTIONS $TARGET
rm $TARGET
sleep 1

shift
done

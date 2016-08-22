WD=`pwd`
cd $WORKDIR

module purge
module load openmpi/intel-13.0/1.6.3/64
module load intel-mkl/11.0/1/64
module load fftw/intel-13.0/3.3.3

HOSTS=`echo $SLURM_JOB_NODELIST | uniq`
for HOST in $HOSTS
do
  ssh $HOST mkdir -p $SCRDIR
done

export VIADEV_USE_AFFINITY=0
srun -n $NPROC /scratch/gpfs/kuangy/compile/Nwchem_emb/Nwchem-6.5-nlemb/bin/LINUX64/nwchem ${JOBNAME}.inp >& ${JOBNAME}.log
sleep 3

for HOST in $HOSTS
do
  ssh $HOST rm -rf $SCRDIR
done

cd $WD

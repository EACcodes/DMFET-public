# Installation

### 1. Embedding NWchem Installation

* First create your installation path (/embnwchem/path/):

		mkdir /embnwchem/path/
		cd /embnwchem/path/

 * The current patch is only compatible with Nwchem-6.5, the source code of which can be found in tigress-hsm. Obtain the code and unpack it:

 		cp /tigress-hsm/kuangy/compile/Nwchem-6.5.revision26243-src.2014-09-10.tar.gz ./
		tar xfvz Nwchem-6.5.revision26243-src.2014-09-10.tar.gz

* Patch the original Nwchem source code (assume you already have the embedding package downloaded in /dmfet_codes/), and enter the src folder:

		cp -r /dmfet_codes/nwchem_src/* Nwchem-6.5/src/
		cd Nwchem-6.5/src

* Open compile.sh file and configure it using your own settings. Especially, update the installation path in compile.sh:

		export NWCHEM_TOP=/embnwchem/path/Nwchem-6.5
		
	The existing MPI/MKL/FFTW settings in compile.sh were configured in compatible with the following modules in Tiger, so make sure you set them in .bashrc before compiling:

		module load openmpi/intel-13.0/1.6.3/64
		module load intel-mkl/11.1/0/64
		module load fftw/intel-13.0/3.3.3

* Compile embedding Nwchem:

		./compile.sh

	If success, the executable (**/embnwchem/path/Nwchem-6.5/bin/LINUX64/nwchem**) should be generated.

* NOTE: in this package, we only dicuss the nonlocal DMFET scheme in GTO function space. However, this version of embedding Nwchem also supports local embedding scheme based on realspace uniform grid, which will be discussed later.

### 2. Configure Wrappers

* All python wrappers are located in /dmfet_codes/optimizer_wrappers. The submission code is **submit.sh**, which generates the run script that calls **optimize.py** when submitted via SLURM scheduler. **optimize.py** then interfaces with embedding Nwchem to perform OEP optimization. 

	**optimize.py** calls the l-bfgs function from scipy, so make sure both **python2.7** and **scipy** packages are available.

	Before use, please configure the following settings in the wrappers, according to your environment.

* In **nwchem_runfile_template.sh**, make sure the correct mpi/mkl/fftw modules are loaded: they have to match the embedding Nwchem  configurations. So update the following lines:

		module load openmpi/intel-13.0/1.6.3/64
		module load intel-mkl/11.0/1/64
		module load fftw/intel-13.0/3.3.3

	Update the path of embedding Nwchem:

		srun -n $NPROC /embnwchem/path/Nwchem-6.5/bin/LINUX64/nwchem ${JOBNAME}.inp >& ${JOBNAME}.log

* In **wrappers.py**, update the scratch directory in line 9:

		scratch_dir = '/your/scratch'

	( You can leave the **embutil_path** variable unchanged. This setting was used in the molpro implementation, which was based on uniform grids. This tag is obsolete in the Nwchem implementation)

<br><br>
# Example & Introduction to Inputs

### 1. Howto Run the Example

* Enter the example folder:

		cd /dmfet_codes/examples/dmfet

* Copy all the wrapper files:

		cp /dmfet_codes/wrappers/* ./

* Run example:

		./submit.sh -np 1 -walltime 1:00:00 ch3ch2br

	Here, you can use **-np** and **-walltime** options to set the number of processors and walltime. And **ch3ch2br** is the job name, which has to match the input file name **ch3ch2br.xyz** and **ch3ch2br.input**

* The reference output is in **/dmfet_codes/examples/dmfet/reference/**.

### 2. Introduction to Inputs

* Typically, seven separate input files are needed: 

	* **jobname(ch3ch2br).input**: specifies the xyz file name, and the nonlocal grid dimension, which is simply the number of AO basis functions. If it is a fresh start with a reference calculation (see optimize.py input section), ***the nonlocal grid dimension can be simply set to 0***, as it will be automatically updated and overwriten by **optimize.py** later. However, if you are performing a continued job, this setting has to be correct.

	* **jobname.xyz**: specifies the coordinates of the system. The format is fairly self-explanory.

	* **cluster/environ.input**: specifies the xyz file name for cluster/environment calculations. Other settings are the same to **jobname.input**

	* **cluster/environ.xyz**: xyz files for cluster/environment. NOTE: since typically we use full dimer basis sets, so these two xyz files are identical to **jobname.xyz**. Different atoms will be set to be the ghost atoms later in **optimize.py**.

	* **basis(lanl2dz).bas**: the basis set file in Nwchem format. NOTE: do not forget Nwchem uses prefix "bq" to indicate ghost atoms, whose basis sets have to be set too (e.g., bqC, bqH, bqBr...) !

* In **optimize.py**, line 10:

		i_restart = 0

	This specifies the restart option: 

	* 0 means a fresh start. In this case, the code will perform a reference calculation on the entire system first to find out the reference density matrix, as well as the number of basis functions. The optimization will start from zero embedding potentials.

	* 1 means a restart. In this case, the code assumes the reference already exists and the basis function numbers in **input** files have already been updated. The code will read the initial guess file, whose name can be specified in line 11:

			restart_embpot = 'embpot.dat.final'

* In **optimize.py**, line 14-59: 

	* Specify the basis file name:

			bas='lanl2dz.bas'

	* The QM method (in Nwchem input style, see Nwchem manual for more info):

			dft_b3lyp="""\                               
			 dft                                         
			 iterations 150                              
			 xc b3lyp                                    
			# convergence energy 1.0e-5                  
			# convergence density 1.0e-4                 
			 grid fine                                   
			 smear 0.05                                  
			 vectors input INIT.MOVECS output OUT.MOVECS\
			""" 
			...
			method = dft_b3lyp

		***NOTE***: if **"vectors input INIT.MOVECS output OUT.MOVECS"** is set, then the MO vectors generated in the current OEP iteration will be used as the initial guess in the next OEP iteration.

		***NOTE***: currently, the wrapper only supports DFT and HF.

		***TIPS***: for open-shell subsystems, smearing is typically necessary for SCF convergence. BTW, ***pay special attentions to the SCF convergence of the subsystems during the OEP optimization***!

	* Define the cluster and the environment. The atoms do not belong to the system will be set as ghosts.

			cluster_atoms = [1,2,3,4]
			environ_atoms = [5,6,7,8]

		***NOTE***: Different to the python convention, the atom index starts from 1 (instead of 0) here!

	* Set the ECPs in Nwchem format:

			ecp_setting="""
			...
			"""
 
* In **optimize.py**, line 71 & 91 & 94:

		calc_ref = nwchem_calculator( '%s'%jobname, '%s.xyz'%jobname, basis='read:%s'%bas\
		           , method=method, lembed=0, lnlembed=0 )

		calc_cluster = nwchem_calculator( 'cluster', 'cluster.xyz', basis='read:%s'%bas \
                                , method=method, lembed=0, lnlembed=1, charge= 0, spin=1 )

		calc_environ = nwchem_calculator( 'environ', 'environ.xyz', basis='read:%s'%bas \
                                , method=method, lembed=0, lnlembed=1, charge= 0, spin=1 )

	These define the charge & spin (i.e., the number of unpaired electrons) in the total/cluster/environment calculations.

### 3. Outputs

* **run.log**: the log file for the OEP optimization.

* **nlembpot.dat**: current nonlocal embedding potential, in matrix format and atomic unit.

* **embpot.dat.final**: converged nonlocal embedding potential, in vector format (reshaped from matrix format) and atomic unit.

<br>
<br>
# Guides for Embedding Nwchem

The embedding Nwchem can be used separately without the python wrappers, which were designed for the purpose of nonlocal DMFET. This section introduces the new features in the embedding Nwchem, compared to the original Nwchem code. 

### 1. Local Embedding Calculations

* Embedding Nwchem supports local embedding potentials represented on a uniform grid. A simple example can be found in:

		cd /dmfet_codes/examples/loc_emb/

	To run, first updat the path in **snwchem.emb**:

		NWCHEMHOME=/embnwchem/path/

	Then run:

		./snwchem.emb -np 8 -walltime 1:00:00 -mem 16000mb cluster.inp

	**snwchem.emb** allows you to specify number of processors, walltime and memory via command line.

* When running local embedding calculations, the code requires two extra files:

	* **emb.in**: this file specifies the grid dimension, the periodic boundary condition box size of the potential. Tag **shift** is never tested, keep it zeros. The format is the same as embedding molcas.

	* **embpot.dat**: defines the 3d local potential, in the same format with embedding molcas.

* In the Nwchem input file, add:

		set embed logical true

	which will tell the code to enable the local embedding calculation.

### 2. Nonlocal Embedding Calculations

* Embedding Nwchem supports the nonlocal embedding calculations, which enables us to perform DMFET introduced above. The nonlocal embedding potential has to be represented in GTO space, as an 1-e integral matrix. This function can also be used to do local embedding calculations, if the provided 1-e integral matrix was computed using Caroline's **EmbeddingIntegralGenerator** code.

* Usually, the python wrappers will automatically take care of the input setups, but the nonlocal embedding Nwchem can be run independently.

* An example can be found in:

		cd /dmfet_codes/examples/nonloc_emb/

	Again, first update the path in **snwchem.emb**

		NWCHEMHOME=/embnwchem/path/

	Then run:

		./snwchem.emb -np 1 -walltime 1:00:00 cluster.inp

	Again, **snwchem.emb** allows you to specify number of processors, walltime and memory via command line.

* Here, no **emb.in** is needed, but the nonlocal embedding matrix has to be provided, in file ** nlembpot.dat**.

* The key inputs in the Nwchem input file are the following: 

		set nlembed logical true 

	This tells the code to read the nonlocal embedding potential.

		property
		 densmat
		end

		task dft property

	This tells the code to print the converged 1-e density matrix.


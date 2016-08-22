#!/usr/bin/python2.7
import sys
import os
import commands
import re
from numpy import *

embutil_path='/scratch/gpfs/kuangy/working_fall2015/dmfet_develop/EmbeddingIntegralGenerator/'
scratch_dir ='/scratch/kuangy'
exe_int1e = embutil_path + 'emb_ints.exe'
exe_dmconvert = embutil_path + 'convert_dm.exe'

def ipt2r( ipt, n1, n2, n3, box ):
    ix = ipt%n1
    iy = (ipt/n1)%n2
    iz = ipt/(n1*n2)
    return box[0,:]/n1*ix + box[1,:]/n2*iy + box[2,:]/n3*iz

def print_embedding_pot( data, ntot, embpotfn='embpot.dat' ):
    # print embpot.dat
    ofile = file( embpotfn, 'w' )
    for i in range(ntot):
        print >>ofile, '%20.8e'%data[i]
    ofile.close()

def read_embedding_pot( ntot, embpotfn='embpot.dat' ):
    extpot = array([0.0 for i in range(ntot)])
    # read local embedding potential
    ifile = file( embpotfn, 'r' )
    for i in range(ntot):
        line = ifile.readline()
        extpot[i] = float(line)
    ifile.close()
    return extpot

def print_nlembedding_pot( data, nbas, nlembpotfn='nlembpot.dat' ):
    ofile = file( nlembpotfn, 'w' )
    nprint = 5
    for ibas in range(nbas):
        iprint = 0
        for jbas in range(nbas):
            ielem = ibas*nbas + jbas
            if iprint%nprint == 0:
                print >> ofile, '%18.8e'%data[ielem],
            else:
                print >> ofile, '%17.8e'%data[ielem],
            iprint += 1
            if iprint%nprint == 0:
                print >> ofile, ''
        print >> ofile, ''
    ofile.close()

def read_nlembedding_pot( nbas, nlembpotfn='nlembpot.dat' ):
    nlembpot = []
    ifile = file( nlembpotfn, 'r' )
    for line in ifile:
        words = line.split()
        for word in words:
            nlembpot.append(float(word))
    if len(nlembpot) != nbas**2:
        sys.exit('Wrong dimension in nlembpot file: %s %d'%(nlembpotfn, nbas))
    ifile.close()
    return array(nlembpot)

class emb_input_file_class:
    def __init__(self, ifn):
        self.ifn = ifn
        self.n1, self.n2, self.n3, self.ntot, self.n_nonloc, self.box, self.xyzfn = self.read_info()

    def read_info(self):
        ifile = file( self.ifn, 'r' )
        iread = 0
        flag = 1
        n1 = 0
        n2 = 0
        n3 = 0
        ntot = 0
        n_nonloc = 0
        box = array([[0.0 for i in range(3)] for j in range(3)])
        xyzfn = ''
        for line in ifile:
            words = line.split()
            if iread == 0 and 'XYZ FILE' in line:
                iread = 1
                continue
            if iread == 0 and 'GRID DIMENSIONS' in line:
                iread = 2
                continue
            if iread == 0 and 'EMBEDDING INTS REAL' in line:
                flag = 0
            if iread == 0 and 'NONLOCAL GRID DIMENSION' in line:
                iread = 3
                continue
            if iread == 0 and 'LATTICE VECTORS' in line:
                iread = 4
                continue
            if iread == 1:
                xyzfn = line.split()[0]
                iread = 0
                continue
            elif iread == 2:
                [ n1, n2, n3 ] = [ int(words[i]) for i in range(3) ]
                ntot = n1 * n2 * n3
                iread = 0
                continue
            elif iread == 3:
                n_nonloc = int(line)
                iread = 0
                continue
            elif iread == 4:
                box = [[ float(words[i]) for i in range(3) ]]
                iread += 1
                continue
            elif iread == 5:
                box.append( [ float(words[i]) for i in range(3) ] )
                iread += 1
                continue
            elif iread == 6:
                box.append( [ float(words[i]) for i in range(3) ] )
                box = array(box)
                iread = 0
                continue
        return n1, n2, n3, ntot, n_nonloc, box, xyzfn

    def update( self, key, value ):
        # update the file in disk
        ofile = file( self.ifn+'.tmp', 'w' )
        ifile = file( self.ifn, 'r' )
        iread = 0
        for line in ifile:
            if iread == 0 and key in line:
                print >> ofile, line,
                iread = 1
            elif iread == 1:
                print >> ofile, value
                iread = 0
            else:
                print >> ofile, line,
        ifile.close()
        ofile.close()
        os.system('mv %s %s'%(self.ifn+'.tmp', self.ifn))
        # update the info in data structure
        self.n1, self.n2, self.n3, self.ntot, self.n_nonloc, self.box, self.xyzfn = self.read_info()

class molpro_calculator():
    def __init__( self, jobname, xyzfn, basis, method, lembed, lnlembed, charge=0, spin=0, nproc=1, mem='500 mb' ):
        self.jobname = jobname
        self.xyzfn = xyzfn
        self.basis = basis
        self.method = method.upper()
        self.lembed = lembed
        self.lnlembed = lnlembed
        self.nproc = nproc
        self.mem = mem
        self.comfn = self.jobname + '.com'
        self.charge = charge
        self.spin = spin
        self.extra_settings = []
        if self.lembed or self.lnlembed: # embedding input file
            self.embinputfile = emb_input_file_class( self.jobname + '.input' )
            if self.embinputfile.xyzfn != self.xyzfn:
                sys.exit('xyz file name in input file does not conform.')

    def set_files(self, path):
        # build com file
        self.set_path( path )
        if os.path.isdir(path): # ensures a clean start
            os.system('rm -r %s'%path)
        os.system('mkdir %s'%path)
        comfn = path + '/' + self.jobname + '.com'
        self.runpath = path
        ofile = file( comfn, 'w' )
        print >> ofile, "***,", self.jobname
        # restart from the wavefunction from last step
#        print >> ofile, 'file,2,%s.wfu'%self.jobname
        print >> ofile, "\nsymmetry,nosym"
        print >> ofile, "geometry={"
        xyzfile = file( self.xyzfn, 'r' )
        for line in xyzfile:
            words = line.split()
            if len(words) > 0:
                print >> ofile, line,
        xyzfile.close()
        print >> ofile, "}\n"
        if self.basis.split(':')[0] != 'read':
            print >> ofile, self.basis
        else:
            basfile = file( self.basis.split(':')[1] )
            for line in basfile:
                print >> ofile, line,
            basfile.close()
        print >> ofile, 'set,charge=%d'%self.charge
        print >> ofile, 'set,spin=%d'%self.spin
        for line in self.extra_settings:
            print >>ofile, line
        print >> ofile, ''
        if self.lembed:
            print >> ofile, "{matrop\nload,one_int,H0\nread,emb_int,FILE=emb_ints.molpro"
            print >> ofile, "add,H0,one_int,emb_int\nsave,H0,1200.1,h0\n}"
        print >> ofile, self.method
#       a couple of more times
#        for i in range(3):
#            print >> ofile, self.method
        print >> ofile, "en_emb = energy"
        print >> ofile, "{matrop;\nload,charged,den,,type=charge;print,charged\n}\n"
        print >> ofile, "show en_emb"
        ofile.close()
        if self.lembed:
            self.set_embfiles(path)

    def set_embfiles(self, path):
        # set up the files related to embedding
        inputfn = '%s.input'%self.jobname
        os.system('cp %s %s'%(inputfn, path))
        os.system('cp embpot.dat %s'%path)
        os.system('cp %s.xyz %s'%(self.jobname, path))
       
    # Run a complete single point calculation: 1e integral, SCF, and compute density matrix 
    def calculate(self):
        if not hasattr( self, 'runpath' ):
            sys.exit('runpath not set in molpro object, need to set files before calculate')
        self.calc_int1e()
        self.run_molpro()
        #self.calc_densmat()

    def calc_int1e(self):
        if not hasattr( self, 'runpath' ):
            sys.exit('runpath not set in molpro object, need to set files before calculate')
        # reprint the 1e integrals in the molpro format
        ifile = file( '%s/embpot.dat'%self.runpath )
        ofile = file( '%s/emb_ints.molpro'%self.runpath, 'w' )
        nbas = self.embinputfile.n_nonloc
        print >> ofile, 'BEGIN_DATA,'
        print >> ofile, '# MATRIX EMB_INT            EMB    SYMMETRY=1'
        i_elem = 0
        i_print = 0
        for line in ifile:
            words = line.split()
            for word in words:
                integral = float(word)
                print >> ofile, '%20.13e,'%integral,
                i_elem += 1
                i_print += 1
                if i_print%5 == 0:
                    print >> ofile, ''
                if i_elem%nbas == 0:
                    if i_print%5 != 0:
                        print >> ofile, ''
                    i_print = 0
        print >> ofile, 'END_DATA,'
        ifile.close()
        ofile.close()

    def run_molpro(self):
        if not hasattr( self, 'runpath' ):
            sys.exit('runpath not set in molpro object, need to set files before calculate')
        runfile = file( 'run.sh', 'w' )
        print >> runfile, '#!/bin/bash\n'
        print >> runfile, 'WORKDIR=%s'%self.runpath
        print >> runfile, 'NPROC=%d'%self.nproc
        print >> runfile, 'MEM=%s'%self.mem
        print >> runfile, 'JOBNAME=%s'%self.jobname
        ifile = file( 'molpro_runfile_template.sh', 'r' )
        for line in ifile:
            print >> runfile, line,
        ifile.close()
        runfile.close()
        os.system('bash run.sh')
        os.system('rm run.sh')

    #def calc_densmat(self):
    #    if not hasattr( self, 'runpath' ):
    #        sys.exit('runpath not set in molpro object, need to set files before calculate')
    #    print >> runfile, '#!/bin/bash'
    #    print >> runfile, 'WD=`pwd`'
    #    print >> runfile, 'module load intel/14.0/64/14.0.0.080 intel-mkl/11.1/0/64'
    #    print >> runfile, 'cd %s'%self.runpath
    #    print >> runfile, '%s %s.input %s.out > densmat_convert.log'%(exe_dmconvert, self.jobname, self.jobname)
    #    print >> runfile, 'cd $WD'
    #    runfile.close()
    #    os.system('bash run.sh')
    #    os.system('rm run.sh')

    def read_energy(self):
        if not hasattr( self, 'runpath' ):
            sys.exit('runpath not set in molpro object, need to set files before calculate')
        tag = commands.getoutput( 'grep CONVERGENCE\ REACHED %s/%s.out | wc'%(self.runpath, self.jobname) )
        if int(tag.split()[0]) == 0:
            sys.exit('ERROR: SCF convergence not reached in %s'%self.runpath)
        line = commands.getoutput( 'grep EN_EMB %s/%s.out | tail -n1 '%(self.runpath, self.jobname) )
        words = line.split()
        return float(words[2])

    #def read_gradient(self):
    #    ifile = file( '%s/dm_grids.dat'%self.runpath, 'r' )
    #    grad = []
    #    for line in ifile:
    #        grad.append(float(line))
    #    if len(grad) != self.embinputfile.ntot + self.embinputfile.n_nonloc:
    #        sys.exit('The length of dm_grids.dat does not conform with input file.')
    #    return array(grad)
    def read_gradient(self):
        ifile = file( '%s/%s.out'%(self.runpath,self.jobname), 'r')
        grad = []
        iread = 0
        for line in ifile:
            words = line.split()
            if 'MATRIX CHARGED' in line:
                iread = 1
                continue
            if iread == 1 and 'SYMMETRY BLOCK' in line:
                iread = 2
                continue
            if iread == 2 and len(words) == 0:
                iread = 0
                continue
            if iread == 2:
                for word in words:
                    grad.append(float(word))
        if len(grad) != self.embinputfile.n_nonloc**2:
            sys.exit('The length of dm in out file does not conform with input file.')
        return array(grad)

    def clean_files(self, path=''):
        if path == '':
            path = self.runpath
        os.system('rm -r %s'%path)

    def set_path(self, path):
        self.runpath = path

############################################################################################
#------------------------------------ Wrappers for NWchem ---------------------------------#
############################################################################################

class nwchem_calculator(molpro_calculator):
    def set_ghosts( self, ghlist ): # ghost atom list, index start from 1 !
        self.ghlist = list(ghlist)
    def set_files(self, path):
        # build inp file
        self.set_path( path )
        if os.path.isdir(path): # ensures a clean start
            os.system('rm -r %s'%path)
        os.system('mkdir %s'%path)
        inpfn = path + '/' + self.jobname + '.inp'
        self.runpath = path
        ofile = file( inpfn, 'w' )
        jobid = commands.getoutput('echo $SLURM_JOB_ID')
        print >> ofile, 'scratch_dir %s/%s'%(scratch_dir,jobid)
        print >> ofile, 'memory total %s'%self.mem
        print >> ofile, 'echo'
        print >> ofile, 'title \"%s\"'%self.jobname
        print >> ofile, ''
        print >> ofile, 'charge %s'%self.charge
        if self.lnlembed:
            print >> ofile, 'set nlembed logical true'
        elif self.lembed:
            print >> ofile, 'set embed logical true'
        # print geometry
        print >> ofile, 'geometry units angstroms nocenter noautoz\nsymmetry c1'
        xyzfile = file( self.xyzfn, 'r' )
        i_atom = 1
        for line in xyzfile:
            words = line.split()
            if len(words) == 4: # atom information 
                atype = words[0]
                [x,y,z] = [float(words[i]) for i in range(1,4)]
                # the input xyz may have atom numbers, remove it for nwchem input
                atype = re.search("[A-Za-z]+",atype).group(0) 
                # ghost atoms?
                if hasattr(self, 'ghlist') and i_atom in self.ghlist:
                    atype = 'bq' + atype
                print >> ofile, '%5s%15.8f%15.8f%15.8f'%(atype,x,y,z)
                i_atom += 1
        print >> ofile, 'end\n' # end of geometry section
        # print basis funciton
        print >> ofile, 'basis cartesian segment'
        if self.basis.split(':')[0] != 'read':
            sys.exit('Only support user-defined basis in NWchem, sorry :-)')
        else:
            basfile = file( self.basis.split(':')[1] )
            for line in basfile:
                print >> ofile, line,
            basfile.close()
        print >> ofile, 'end\n'
        # print computational method
        self.flag = 0 # 0: hf ; 1: dft
        if re.search('dft(?i)',self.method):
            self.flag = 1
        if self.flag:
            print >> ofile, 'dft'
            print >> ofile, self.method
            print >> ofile, ' mult %d'%(self.spin+1)
            print >> ofile, 'end\n'
        else:
            print >> ofile, 'scf'
            if self.spin == 0:
                print >> ofile, 'singlet'
            elif self.spin == 1:
                print >> ofile, 'doublet'
            elif self.spin == 2:
                print >> ofile, 'triplet'
            elif self.spin == 3:
                print >> ofile, 'quartet'
            else:
                sys.exit('You must be insane to do embedding with such high spins...')
            print >> ofile, self.method
            print >> ofile, 'end\n'
        # print densmat output
        print >> ofile, 'property'
        print >> ofile, ' densmat'
        print >> ofile, 'end\n' 
        for line in self.extra_settings:
            print >> ofile, line
        print >> ofile, ''
        if self.flag:
            print >> ofile, 'task dft property'
        else:
            print >> ofile, 'task scf property'
        ofile.close()
        if self.lembed:
            self.set_embfiles(path)
        elif self.lnlembed:
            self.set_nlembfiles(path)
        # initial guesses
        if os.path.isfile('%s.movecs'%self.jobname):
            os.system('cp %s.movecs %s/INIT.MOVECS'%(self.jobname, path))
    def set_nlembfiles(self, path):
        # set up the files related to embedding
        os.system('cp nlembpot.dat %s'%path)
    def calculate(self):
        if not hasattr( self, 'runpath' ):
            sys.exit('runpath not set in molpro object, need to set files before calculate')
        self.run_nwchem()
    def run_nwchem(self):
        if not hasattr( self, 'runpath' ):
            sys.exit('runpath not set in nwchem object, need to set files before calculate')
        runfile = file( 'run.sh', 'w' )
        print >> runfile, '#!/bin/bash\n'
        print >> runfile, 'WORKDIR=%s'%self.runpath
        print >> runfile, 'NPROC=%d'%self.nproc
        print >> runfile, 'JOBNAME=%s'%self.jobname
        print >> runfile, 'SCRDIR=%s/$SLURM_JOB_ID'%scratch_dir
        ifile = file( 'nwchem_runfile_template.sh', 'r' )
        for line in ifile:
            print >> runfile, line,
        ifile.close()
        runfile.close()
        os.system('bash run.sh')
        os.system('rm run.sh')
    def read_energy(self):
        if not hasattr( self, 'runpath' ):
            sys.exit('runpath not set in molpro object, need to set files before calculate')
        if self.flag:
            line = commands.getoutput('grep Total\ DFT\ energy %s/%s.log '%(self.runpath, self.jobname)) 
        else:
            line = commands.getoutput('grep Total\ SCF\ energy %s/%s.log '%(self.runpath, self.jobname))
        words = line.split('=')
        return float(words[1])
    def read_gradient(self):
        ifn =  '%s/densmat.dat'%self.runpath
        nbas = self.embinputfile.n_nonloc
        grad = read_nlembedding_pot( nbas, ifn )
        return array(grad)
    def read_nbas(self):
        ifn = '%s/%s.log'%(self.runpath, self.jobname)
        if self.flag:
            line = commands.getoutput( 'grep \'AO basis - number of functions\' %s'%ifn )
            nbas = int( line.split(':')[1] )
        else:
            line = commands.getoutput( 'grep functions %s'%ifn )
            nbas = int( line.split('=')[1] )
        return nbas

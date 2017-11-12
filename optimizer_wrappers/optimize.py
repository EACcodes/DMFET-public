#!/usr/bin/python2.7
import sys
from numpy import *
from scipy import optimize
from wrappers import *
from functions import *

# set up nonlocal potential
jobname = sys.argv[1]
i_restart = 0
restart_embpot = 'embpot.dat.final'

# set up run parameters
bas='lanl2dz.bas'
dft_b3lyp="""\
 dft
 iterations 200
 xc b3lyp
 #tolerances tight
 #convergence energy 1.0e-8
 #convergence density 1.0e-7
 #convergence gradient 1.0e-6
 grid fine
 smear 0.05
 vectors input INIT.MOVECS output OUT.MOVECS\
"""
hf="""\
 hf
 maxiter 150\
"""
method = dft_b3lyp
cluster_atoms = [1,2,3,4]
environ_atoms = [5,6,7,8]

pc_setting_cluster="""\
bq
      0.11216457      1.08420644     -0.00393857 1
end
"""
pc_setting_environ="""\
bq
      0.00000000      0.00000000      0.00000000 -1
end
"""

ecp_setting="""\
ECP
Br nelec 28
Br ul
1    213.6143969            -28.0000000        
2     41.0585380           -134.9268852        
2      8.7086530            -41.9271913        
2      2.6074661             -5.9336420        
Br S
0     54.1980682              3.0000000        
1     32.9053558             27.3430642        
2     13.6744890            118.8028847        
2      3.0341152             43.4354876        
Br P
0     54.2563340              5.0000000        
1     26.0095593             25.0504252        
2     28.2012995             92.6157463        
2      9.4341061             95.8249016        
2      2.5321764             26.2684983        
Br D
0     87.6328721              3.0000000        
1     61.7373377             22.5533557        
2     32.4385104            178.1241988        
2      8.7537199             76.9924162        
2      1.6633189              9.4818270        
END\
"""

# read embedding input files
ifile_ref = emb_input_file_class( '%s.input'%jobname )
ifile_cluster = emb_input_file_class( 'cluster.input' )
ifile_environ = emb_input_file_class( 'environ.input' )
[ n1, n2, n3 ] = [ ifile_ref.n1, ifile_ref.n2, ifile_ref.n3 ]
nbas = ifile_ref.n_nonloc
ndim = nbas*nbas # fully nonlocal density matrix
box = ifile_ref.box
xyzfn = ifile_ref.xyzfn

# print embpot.dat and nlembpot.dat for reference calculation
calc_ref = nwchem_calculator( '%s'%jobname, '%s.xyz'%jobname, basis='read:%s'%bas\
           , method=method, lembed=0, lnlembed=0, charge=0, spin=0, nproc=16, mem="5000 mb" )
calc_ref.extra_settings.append(ecp_setting)
calc_ref.set_path('ref')
if i_restart == 0: # a fresh start, no reference exists
    # do reference calculation
    calc_ref.clean_files('ref') 
    calc_ref.set_files('ref')   
    calc_ref.calculate()        
    nbas = calc_ref.read_nbas()
    ndim = nbas**2
    ifile_ref.update('NONLOCAL GRID DIMENSION', nbas)
    ifile_cluster.update('NONLOCAL GRID DIMENSION', nbas)
    ifile_environ.update('NONLOCAL GRID DIMENSION', nbas)
    extpot0 = array([ 0.0 for i in range(ndim) ])
    print_nlembedding_pot( extpot0, nbas, 'nlembpot.dat' )
calc_ref.embinputfile = ifile_ref
rho_ref = calc_ref.read_gradient()

# setup cluster and environ calculator
calc_cluster = nwchem_calculator( 'cluster', 'cluster.xyz', basis='read:%s'%bas \
                                , method=method, lembed=0, lnlembed=1, charge=-1, spin=0, nproc=16, mem="5000 mb" )
calc_cluster.extra_settings.append(ecp_setting)
calc_cluster.extra_settings.append(pc_setting_cluster)
calc_environ = nwchem_calculator( 'environ', 'environ.xyz', basis='read:%s'%bas \
                                , method=method, lembed=0, lnlembed=1, charge= 1, spin=0, nproc=16, mem="5000 mb" )
calc_environ.extra_settings.append(ecp_setting)
calc_environ.extra_settings.append(pc_setting_environ)
# set dummy atoms
calc_cluster.set_ghosts( environ_atoms )
calc_environ.set_ghosts( cluster_atoms )

calc_cluster.clean_files('cluster')
calc_environ.clean_files('environ')
calc_cluster.set_files('cluster')
calc_environ.set_files('environ')

# Definition of Lagrangian
eval_counter = 0
logfile = file( 'run.log', 'w' )
def Lagrangian( params, *args ):
    global eval_counter
    # prepare embedding file:
    print_nlembedding_pot( params, nbas, 'nlembpot.dat' )
    # prepare for cluster & environ run
    calc_cluster.clean_files('cluster')
    calc_environ.clean_files('environ')
    calc_cluster.set_files('cluster')
    calc_environ.set_files('environ')
    # run cluster & environ
    calc_cluster.calculate()
    calc_environ.calculate()
    # read result
    W = -( calc_cluster.read_energy() + calc_environ.read_energy() - dot(params, rho_ref) )
    grad = -( calc_cluster.read_gradient() + calc_environ.read_gradient() - rho_ref )
    if eval_counter == 0 and i_restart == 0:
        os.system('cp -r cluster cluster.0')
        os.system('cp -r environ environ.0')
    # output
    eval_counter += 1
    print >>logfile, '--------------------'
    print >>logfile, ' Evaluation No: %d'%eval_counter
    print >>logfile, ' Maximum of gradient: %20.8e'%max(abs(grad))
    print >>logfile, ' Norm of gradient: %20.8e'%sqrt(dot(grad,grad))
    print >>logfile, ' Lagrangian value: %20.8e'%W
    # add penalty
    #pen, pen_grad = penalty( params[:ntot], n1, n2, n3 )
    #W -= pen
    #grad[:ntot] -= pen_grad[:]
    #print >>logfile, ' Maximum of gradient with penalty: %20.8e'%max(abs(grad))
    #print >>logfile, ' Norm of gradient with penalty: %20.8e'%sqrt(dot(grad,grad))
    #print >>logfile, ' Lagrangian value with penalty: %20.8e'%W
    logfile.flush()
    os.system('cp cluster/OUT.MOVECS cluster.movecs')
    os.system('cp environ/OUT.MOVECS environ.movecs')
    return W, grad

# read initial guesses
if i_restart:
    extpot0 = read_nlembedding_pot( nbas )

## test gradient
#delta = 1.0e-3
#dn = 2
#results = [ 0 for i in range(2*dn+1) ]
#W, grad = Lagrangian( extpot0, () )
#direct = grad
##results[dn] = W
##direct = array([0.0 for i in range(1000001)])
##direct[1000000] = 1.0
#for i in range(-dn, dn+1):
##    if i == 0:
##        continue
#    extpot = extpot0 + i*delta*direct
#    W, g = Lagrangian( extpot, () )
#    os.system('mv cluster cluster.%s'%eval_counter)
#    os.system('mv environ environ.%s'%eval_counter)
#    results[i+dn] = W
#print '# %20.10e'%dot(direct,grad)
#for i in range(-dn, dn+1):
#    print '%15.6e%20.8e'%(i*delta,results[i+dn])

# the real optimization
x,f,d = optimize.fmin_l_bfgs_b( Lagrangian, x0=extpot0, args=(), factr=1e4, pgtol=1e-05 )
print_embedding_pot( x, ndim, 'embpot.dat.final' )
print >>logfile, d['warnflag']
if d['warnflag'] == 2:
    print d['task']
W,grad = Lagrangian( x )
logfile.close()

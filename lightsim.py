import pynetlogo
import numpy as np
import networkx as nx
import os
from itertools import product
from scipy.special import kl_div
from scipy.stats import relfreq, norm

from tqdm import tqdm

from multiprocessing import Pool

## VARIABLES ##
# Resolution
res = 2
# Pruning parameters
beta = np.linspace(0.001, 10, res)
# Separation parameter
dist = np.linspace(0, 10, res)
# Set uncertainty on new observations
varc = 10
# Get dimensional d
d_norm = dist[:]/np.sqrt(varc)
# Set number of agents
N = 100
# Set number of iterations (NEW!!!)
its=100
# Set number of repeat runs
rep = 5
# Where to store values
sv_pth = './b[' + str(beta[0]) + ',' + str(beta[-1]) + ']d[' + str(dist[0]) + ',' + str(dist[-1]) + ']/'

################

# Netlogo installation and model
NL_PATH = '/opt/netlogo/' if os.environ.get('USER')=='daniele' else '/home/lunis/Programs/NetLogo-6.4.0-64'
NL_MODEL = "./EC3.0.nlogo"
NL_GUI=False

#Initial distribution parameters
MUTRUE = 0.
VARTRUE = 1.

# Multiprocessing parameters
workers = 4#os.cpu_count() - 6
csize = 6


# Create iterable single simulation function for multiprocessing
# takes netlogo as a parameter, so that it doesn't start again at each iteration
# in multiprocessing, each thread must have a different netlogo instance
def single_sim_nonl(netlogo, N : int, beta : float, dist : float, var_c : float, iters : int):
    '''Executes a single simulation.
    returns kullback-leibler of initial and final mu distros wrt initial distro'''

    # Set values for simulation
    global_vars = {
        'N': N,
        'beta': beta,
        'var-c': var_c,
        'mutrue': MUTRUE,
        'vartrue': VARTRUE,
        'update-type': 2,
        'network-type': "\"scale-free\"",
        'p': 0.01,
        'pref': 1,
        'initial-sampling' : "\"bivariate\"",
        'dist' : 0.
    }
    global_vars['N'] = N
    global_vars['beta'] = beta
    global_vars['var-c'] = var_c
    global_vars['dist'] = dist

    # Simple function for retrieving nodes variables
    def values(var: str):
        c = netlogo.report(f"map [s -> [{var}] of s] sort nodes")
        return c
    
    # Prepare NetLogo enviroment
    netlogo.command('clear-all')

    for name in global_vars:
        netlogo.command(f'set {name} {global_vars[name]}')

    netlogo.command('setup')
 
    # Initialize results arrays
    mus = np.empty((0,N))
    sigma2s = mus
    nets = []

    mus = np.concatenate((mus, values('mu0')[np.newaxis, :]), axis = 0)
    sigma2s = np.concatenate((sigma2s, values('var0')[np.newaxis, :]), axis = 0)

    nets.append(netlogo.report("[list ([label] of end1) ([label] of end2)] of edges").astype(int))

    # Run the simulation
    for n in range(1,iters+1):
        netlogo.command('go')
        ### save only initial and final mus,sigmas and networks.
        # mus = np.concatenate((mus, values('mu')[np.newaxis, :]), axis = 0)
        # sigma2s = np.concatenate((sigma2s, values('var')[np.newaxis, :]), axis = 0)
        # if (n%t_rep)==0:
        #     nets.append(netlogo.report("[list ([label] of end1) ([label] of end2)] of edges").astype(int))
    mus = np.concatenate((mus, values('mu')[np.newaxis, :]), axis = 0)
    sigma2s = np.concatenate((sigma2s, values('var')[np.newaxis, :]), axis = 0)
    nets.append(netlogo.report("[list ([label] of end1) ([label] of end2)] of edges").astype(int))

    #Create final network to calculate connected components
    G = nx.from_edgelist(nets[-1],create_using=nx.DiGraph)
    for i in sorted(G.nodes):
        (G.nodes)[i]['mu'] = mus[-1,i]
        (G.nodes)[i]['sigma2']= sigma2s[-1,i]
    #Remove edges to distant nodes
    toremove = []
    for edge in G.edges:
        if abs(G.nodes[edge[0]]['mu']-G.nodes[edge[1]]['mu'])>=beta*np.sqrt(G.nodes[edge[0]]['sigma2']):
            toremove.append(edge)
    for edge in toremove:
        G.remove_edge(*edge)
    #Calculate dimensions of strongly and weakly connected components
    scc = sorted([len(i) for i in nx.strongly_connected_components(G)],reverse=True)
    wcc = sorted([len(i) for i in nx.weakly_connected_components(G)],reverse=True)


    #Get KL divergence at start and end
    dv0 = gauss_DV(mus[0,:], int(N/10))
    dv = gauss_DV(mus[-1,:], int(N/10))

    return dv0, dv, scc, wcc

# Divergence function conveniently wrapped
def gauss_DV(x, nbins: int):

    blind_mean = MUTRUE #x.mean()
    blind_dev = VARTRUE # np.sqrt(x.var()) this method creates problems and changes too much per dataset, results are hard to make sense of

    # Get relative frequencies of the data

    x_freq, llim, wbin, __ = relfreq(x, numbins = nbins)

    # Get gaussian frequencies from the 

    y_freq = np.empty((0))

    for i in range(0, nbins):

        y_freq = np.concatenate((y_freq, norm.pdf( (llim + wbin/2 + i*wbin), loc = blind_mean, scale = blind_dev)[np.newaxis]))

    dv = kl_div(x_freq, y_freq).sum()

    return dv

# Create iterable function
def it_single_sim(bd : tuple):
    b, d = bd

    dv0 = np.empty((0))
    dv = np.empty((0))
    scc = []
    wcc= []

    for i in range(0,rep):#,desc='subiter:',leave=False):
        #'netlogo' gets initialized in init_worker 
        d0_, d_, scc_, wcc_  = single_sim_nonl(netlogo,N = N, beta = b, dist = d, var_c = varc, iters=its)

        dv0 = np.concatenate((dv0, d0_[np.newaxis]))
        dv = np.concatenate((dv, d_[np.newaxis]))
        scc.append(scc_)
        wcc.append(wcc_)
    
    dv0_var = dv0.var()
    dv0 = dv0.mean()

    dv_var = dv.var()
    dv = dv.mean()
    
    avg_scc_n = np.mean([len(scc_) for scc_ in scc])
    avg_wcc_n = np.mean([len(wcc_) for wcc_ in wcc])

    return dv0, dv0_var, dv, dv_var, avg_scc_n, avg_wcc_n

def init_worker():
    '''Initialize process by instantiating a netlogolink.
    This way, it gets opened only once per process, saving (some) time.'''
    global netlogo
    # global here is misleading: since multiprocess spawns completely separate python interpreters,
    # it means that the process will be able to access it until it's not closed,
    # not that all processes can see it. So each process gets its own netlogo.
    
    netlogo = pynetlogo.NetLogoLink(
        gui = NL_GUI,
        netlogo_home = NL_PATH,
    )
    netlogo.load_model(NL_MODEL)
    #print('One netlogo started')

    #stop properly netlogo when the process is terminated
    def kill_netlogo():
        if 'netlogo' in globals():
            netlogo.kill_workspace()
            #print('killed netlogo!')
        #else:
        #    print('no netlogo global found')
    import atexit
    atexit.register(kill_netlogo) 

def multiple_sim():
    bdpairs = list(product(beta,dist)) # list of (b,d) tuples
    with Pool(workers,initializer=init_worker) as p:
        ps = list(tqdm(p.imap(it_single_sim, bdpairs),# chunksize = csize),
                        total=len(bdpairs),
                        leave=True)
                    )
    
    psarr = np.array(ps)
    return bdpairs, psarr


if __name__ == '__main__':
    #### EVALUATE THINGS
    bdpairs,psarr = multiple_sim()
    # Save results to local
    os.makedirs(sv_pth, exist_ok = True)
    print('output shape:', psarr.shape)
    np.save(sv_pth + str(res) + 'ps.npy', psarr)
    np.save(sv_pth + str(res) + 'bdpairs.npy', bdpairs)

import numpy as np
import sys
import os
import argparse
from sklearn.linear_model import LinearRegression
from scipy.optimize import curve_fit

bme_dir='/lindorffgrp-isilon/thomasen/software/BME2_original'
sys.path.append(bme_dir)
import BME as BME

blocking_dir='/storage1/thomasen/software/BLOCKING/MonoCV'
sys.path.append(blocking_dir)
import block as block

# Parsing command line arguments

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--iterations", type=int, default=1000, help="Number of iterations")
parser.add_argument("-sr", "--iteration_start_reweighting", type=int, default=500, help="Number of iterations before reweighting begins")
parser.add_argument("-os", "--oligomersize", type=int, help="Largest oligomer size (nr subunits)")
parser.add_argument("-s", "--oligomersizestep", type=int, help="Nr of subunits added in each oligomer")
parser.add_argument("-cp", "--calcpath", type=str, help="Path to calc_data files for BME")
parser.add_argument("-fs", "--fitscaleglobaloffsetlocal", action="store_true", help="Set this flag to fit SAXS scale globally across protein concentrations and offset locally for each protein concentration")
parser.add_argument("-fsl", "--fitscaleoffsetlocal", action="store_true", help="Set this flag to fit scale and offset of calculated SAXS data to experimental locally for each protein concentration")
parser.add_argument("-fc", "--fitctot", action="store_true", help="Set this flag to fit ctot")
parser.add_argument("-fk", "--fitKa", action="store_true", help="Set this flag to fit Ka")
parser.add_argument("-c", "--ctots", nargs="+", type=float, help="Total protein concentrations for isodesmic model. Used as intial value if fitting total protein concentration.")
parser.add_argument("-k", "--ka", type=float, help="Association constant Ka for isodesmic model")
parser.add_argument("-t", "--theta", type=float, default=1000000, help="Initial theta value for BME reweighting")
parser.add_argument("-tf", "--thetafinal", type=float, default=0.1, help="Final theta value for BME reweighting")
parser.add_argument("-td", "--thetadecrease", type=float, default=0.025, help="How fast theta is decreased. Theta is decreased by currentvalue*thetadecrease each iteration.")
parser.add_argument("-p", "--philimit", type=float, default=0.5, help="When phi effective drops below this value, theta stops decreasing.")
parser.add_argument("-ss", "--mcstepsize", type=float, default=0.02, help="Maximum stepsize for simulated annealing of total protein concentration")
parser.add_argument("-tinit", "--temperatureinit", type=float, default=0.001, help="Initial temperature for simulated annealing of total protein concentration")
parser.add_argument("-tstop", "--temperaturestop", type=float, default=0.000001, help="Final temperature for simulated annealing of total protein concentration")
parser.add_argument("-tdecr", "--temperaturedecrease", type=float, default=0.0000001, help="Step-size of temperature decrease for simulated annealing of total protein concentration")
args = parser.parse_args()

#Settings
iterations=args.iterations

largest_oligomer = args.oligomersize
oligomersizestep = args.oligomersizestep
oligomersizes = range(oligomersizestep, largest_oligomer+1, oligomersizestep)

calcpath = args.calcpath

fit_scale_global_offset_local = args.fitscaleglobaloffsetlocal
fit_scale_offset_local = args.fitscaleoffsetlocal
if fit_scale_global_offset_local == True:
    assert fit_scale_offset_local != True, 'Cannot fit scale both locally and globally'

fit_ctot_factor = args.fitctot
fit_Ka_global = args.fitKa

Ka_init = args.ka
ctots = args.ctots

theta=args.theta
theta_min=args.thetafinal
theta_decrease=args.thetadecrease
phi_limit=args.philimit
iteration_start_reweighting=args.iteration_start_reweighting

#Simulated annealing settings
T_init = args.temperatureinit
T_stop = args.temperaturestop
T_decrease = args.temperaturedecrease
factor_MC_step = args.mcstepsize
Ka_MC_step = args.mcstepsize

# Functions for loading files


def load_weights(weightsfile):
    weights = np.genfromtxt(weightsfile, usecols=(1))
    return np.array(weights)

#Loads SAXS intensities calculated from simulations - in BME format
def load_sim_SAXS(calc_file):
    calc_data = np.array([[float(x) for x in line.split()[1:]] for line in open(calc_file) if("#" not in line)])
    return np.array(calc_data)

#Loads experimental SAXS q, I and err - in BME format
def load_exp_SAXS(exp_file):
    q, I, err = np.genfromtxt(exp_file, usecols=(0,1,2), skip_header=0, unpack=True)
    return np.array(q), np.array(I), np.array(err)

# Function to calculate concentrations with isodesmic model

#Calculate concentrations with isodesmic model
#Returns array of oligomer concentrations
def isodesmic(Ka, Ctot, oligomersizes):

    BTBdimer_conc = Ctot/2.0

    C = []

    #Calculate monomer concentration
    conc = (2*BTBdimer_conc*Ka+1-np.sqrt(4*BTBdimer_conc*Ka+1))/(2*BTBdimer_conc*(Ka**2))
    C.append(conc)

    #Calculate oligomer concentrations
    for i in range(1,len(oligomersizes)):
        conc = Ka*C[i-1]*C[0]
        C.append(conc)

    return np.array(C)

# Function for generating uniform weights

#Make initial uniform weights
def uniform_weights(calc_data):
    
    weights = []
    for i in range(len(oligomersizes)):
        weights.append([1/len(calc_data[i])]*len(calc_data[i]))
    
    return np.array(weights)

# Functions for SAXS and getting chi2

#Takes list of Iaverages from the different species and calculates volume fraction weighted average based on conc and oligomersizes
def volume_average_all_SAXS(calc_data, weights, concs, oligomersizes):
    Iaverages = np.sum(calc_data*weights[:,np.newaxis].reshape(calc_data.shape[0],calc_data.shape[1],1),axis=1)
    
    volume_weights = concs*oligomersizes
    volume_weights = volume_weights/np.sum(volume_weights)

    Iaverage_all = np.sum(Iaverages*volume_weights[:,np.newaxis],axis=0)
        
    return Iaverage_all

def get_chi2(Iaverage_sim, I_exp, err_exp):
    
    chi2 = np.average(np.square((Iaverage_sim-I_exp)/err_exp))
    
    return chi2

def block_error_SAXS(calc_data):

    err_sim = []

    for i,size in enumerate(oligomersizes):

        calc_data_Iq_vs_frames = calc_data[i].transpose()

        Iq_err = []
        
        for Iq_vs_frames in calc_data_Iq_vs_frames:
            dataset = block.check(Iq_vs_frames, multi=1)
            blocks_stat = block.blocking(dataset)
            corr_len, err = block.optimal_block(len(dataset),blocks_stat,method="hline")

            Iq_err.append(err)

        err_sim.append(Iq_err)

    err_sim = np.array(err_sim)

    return err_sim

def decompose_SAXS(oligomerindex, concs, calc_data, weights, oligomersizes, q, Iexp, err_exp, err_sim):    
    #Calculate normalized volume fraction for all species
    volume_weights = concs*oligomersizes
    volume_weights = volume_weights/np.sum(volume_weights)
    
    #Remove oligomer species i from all arrays
    oligomersizes_rest = np.delete(oligomersizes, oligomerindex, 0)
    concs_rest = np.delete(concs, oligomerindex, 0)
    calc_data_rest = np.delete(calc_data, oligomerindex, 0)
    weights_rest = np.delete(weights, oligomerindex, 0)
    volume_weights_rest = np.delete(volume_weights, oligomerindex, 0)
    err_sim_rest = np.delete(err_sim, oligomerindex, 0)
    
    #Calculate average SAXS of each oligomer species except oligomer species i
    Iaverages_rest = np.sum(calc_data_rest*weights_rest[:,np.newaxis].reshape(calc_data_rest.shape[0],calc_data_rest.shape[1],1),axis=1)
    
    #Calculate average SAXS profile of all species minus contribution from species i
    Iaverage_rest = np.sum(Iaverages_rest*volume_weights_rest[:,np.newaxis],axis=0)
    
    #Find the contribution from species i
    Iaverage_sim = Iexp - Iaverage_rest
    
    #Scale Iaverage back from volume fraction contribution
    Iaverage_sim = Iaverage_sim/volume_weights[oligomerindex]
    
    #Get error on Iaverage_sim from propagating simulation SAXS errors and experimental SAXS errors
    err_propagated = np.sqrt(np.square(err_exp) + np.sum(np.square(err_sim_rest*volume_weights_rest[:,np.newaxis])))
    dI = err_propagated/volume_weights[oligomerindex]
    
    #Write input file for BME reweighting
    with open("exp_file_sim.dat", 'w') as f:
        f.write("# DATA=SAXS PRIOR=GAUSS\n")
        for i in range(len(Iaverage_sim)):
            f.write("%e \t %e \t %e \n" % (q[i],Iaverage_sim[i],dI[i]))

# Functions for SAXS scale and offset

#Fit calc data profiles scale and offset with linear regression
def fit_global_scale_expdata(Iaverage_vs_ctots, Iexp_vs_ctots, err_exp_vs_ctots):
    #Get weight for each point based on exp error
    sample_weight=1.0/np.square(err_exp_vs_ctots.flatten())
    
    #Linear regression
    reg = LinearRegression(fit_intercept=False).fit(Iaverage_vs_ctots.reshape(-1,1),Iexp_vs_ctots.reshape(-1,1),sample_weight=sample_weight)
    r_value = reg.score(Iaverage_vs_ctots.reshape(-1,1),Iexp_vs_ctots.reshape(-1,1),sample_weight=sample_weight)
    slope = reg.coef_[0]
    
    slope = (1.0/slope)
    
    #Transform exp profiles with new scale and offset
    Iexp_vs_ctots *= slope
    err_exp_vs_ctots *= abs(slope)
    
    return Iexp_vs_ctots, err_exp_vs_ctots, slope, r_value

def fit_offset_expdata(Iaverage, Iexp, err_exp):
    
    #Offset function
    def offset_func(x, a):
        y=x+a
        return y
    
    #Fit offset
    offset = curve_fit(offset_func, Iaverage, Iexp, sigma=err_exp)
    
    #Modify exp data
    offset = -float(offset[0][0])
    Iexp += offset
    
    return Iexp, offset

def evaluate_scale_offset(Iexp, Iexp_start):
    
    #Linear regression
    reg = LinearRegression(fit_intercept=True).fit(Iexp.reshape(-1,1),Iexp_start.reshape(-1,1))
    r_value = reg.score(Iexp.reshape(-1,1),Iexp_start.reshape(-1,1))
    scale, offset = reg.coef_[0],reg.intercept_
    
    return scale, offset

#Fit calc data profiles scale and offset with linear regression
def fit_scale_offset_expdata(Iaverage_sim, I_exp, err_exp):
    #Get weight for each point based on exp error
    sample_weight=1.0/(err_exp**2)
    
    #Linear regression
    reg = LinearRegression(fit_intercept=True).fit(Iaverage_sim.reshape(-1,1),I_exp.reshape(-1,1),sample_weight=sample_weight)
    r_value = reg.score(Iaverage_sim.reshape(-1,1),I_exp.reshape(-1,1),sample_weight=sample_weight)
    slope,intercept = reg.coef_[0],reg.intercept_
    
    slope = (1.0/slope)
    intercept = -intercept
    #Transform exp profiles with new scale and offset
    I_exp_fitted = I_exp+intercept
    I_exp_fitted *= slope
    err_exp_fitted = abs(slope)*err_exp
    
    return I_exp_fitted, err_exp_fitted, slope, intercept, r_value

# Function for BME reweighting

#Runs BME and returns weights
def run_BME(exp_file, calc_file, theta, initial_weights='uniform'):
    
    if initial_weights=='uniform':
        rew = BME.Reweight("BMEreweighting")
    else:
        rew = BME.Reweight("BMEreweighting", w0=initial_weights)
    
    rew.load(exp_file,calc_file)

    chi2_before, chi2_after, phi = rew.fit(theta=theta)
    new_weights = rew.get_weights()
    
    return new_weights, chi2_before, chi2_after, phi


# Functions for optimization of isodesmic parameters

def nudge(param, step):
    
    #Get new param (draw from gaussian centered at old param)
    new = np.random.normal(loc=param, scale=step)
    
    #If new param is <= 0, try again
    while new <= 0:
        new = np.random.normal(loc=param, scale=step)
    
    return new

#Optimizes global Ctot-factor using monte-carlo 
def montecarlo_opt_ctotfactor(T_init, T_stop, T_decrease, ctots, ctot_factor_init, factor_MC_step, Ka, weights_vs_ctots, oligomersizes, Iexp_vs_ctots, exp_err_vs_ctots, calc_data):
    
    #Get initial average chi2
    chi2_old = 0
    for i,ctot in enumerate(ctots):
        #Get initial chi2 
        concs = isodesmic(Ka, ctot*ctot_factor_init, oligomersizes)
        Iaverage = volume_average_all_SAXS(calc_data, weights_vs_ctots[i], concs, oligomersizes)
        chi2_old += get_chi2(Iaverage, Iexp_vs_ctots[i], exp_err_vs_ctots[i])
    chi2_old /= len(ctots)
    
    #Get initial ctot_factor
    ctot_factor_old = ctot_factor_init
    
    T = T_init

    accept = 0
    reject = 0
    j = 0
    
    chi2_lst = []
    ctot_factor_lst = []
    
    chi2_lst.append(chi2_old)
    ctot_factor_lst.append(ctot_factor_old)
    
    #Monte carlo loop
    while T > T_stop:
        
        #Nudge parameters
        ctot_factor_new = nudge(ctot_factor_old, factor_MC_step)

        #Get SAXS and average chi2
        chi2_new = 0
        for i,ctot in enumerate(ctots):
            concs = isodesmic(Ka, ctot*ctot_factor_new, oligomersizes)
            Iaverage = volume_average_all_SAXS(calc_data, weights_vs_ctots[i], concs, oligomersizes)
            chi2_new += get_chi2(Iaverage, Iexp_vs_ctots[i], exp_err_vs_ctots[i])
        chi2_new /= len(ctots)
    
        #Calculate boltzmann factor for metropolis criterion
        criterion = np.exp(-(chi2_new-chi2_old)/T)
    
        #Accept
        if criterion > 1:
            chi2_old = chi2_new
            ctot_factor_old = ctot_factor_new
            accept += 1
    
        #Accept with rng probability or reject
        if criterion <= 1:
            rng = np.random.random()
            if criterion > rng:
                chi2_old = chi2_new
                ctot_factor_old = ctot_factor_new
                accept += 1
            else:
                reject += 1
    
        chi2_lst.append(chi2_old)
        ctot_factor_lst.append(ctot_factor_old)
    
        #Decrease temperature
        T -= T*T_decrease
        j += 1
        
    acceptance = accept/j
    rejection = reject/j

    #Find minimum chi2
    chi2_min = np.amin(chi2_lst)

    #Find index of minimum chi2 in chi2_lst
    min_index = chi2_lst.index(chi2_min)

    #Find ctot and Ka giving rise to minimum chi2 using index from above
    assert len(ctot_factor_lst) == len(chi2_lst)

    ctot_factor_opt = ctot_factor_lst[min_index]
    
    print("Optimal ctot_factor: %f " % ctot_factor_opt)

    return ctot_factor_opt, acceptance, rejection

#Optimizes global Ka using monte-carlo 
def montecarlo_opt_globalKa(T_init, T_stop, T_decrease, ctots, ctot_factor, Ka_MC_step, Ka_init, weights_vs_ctots, oligomersizes, Iexp_vs_ctots, exp_err_vs_ctots, calc_data):
    
    #Get initial average chi2
    chi2_old = 0
    for i,ctot in enumerate(ctots):
        #Get initial chi2 
        concs = isodesmic(Ka_init, ctot*ctot_factor, oligomersizes)
        Iaverage = volume_average_all_SAXS(calc_data, weights_vs_ctots[i], concs, oligomersizes)
        chi2_old += get_chi2(Iaverage, Iexp_vs_ctots[i], exp_err_vs_ctots[i])
    chi2_old /= len(ctots)
    
    #Get initial ctot_factor
    Ka_old = Ka_init
    
    T = T_init

    accept = 0
    reject = 0
    j = 0
    
    chi2_lst = []
    Ka_lst = []
    
    chi2_lst.append(chi2_old)
    Ka_lst.append(Ka_old)
    
    #Monte carlo loop
    while T > T_stop:
        
        #Nudge parameters
        Ka_new = nudge(Ka_old, Ka_MC_step)

        #Get SAXS and average chi2
        chi2_new = 0
        for i,ctot in enumerate(ctots):
            concs = isodesmic(Ka_new, ctot*ctot_factor, oligomersizes)
            Iaverage = volume_average_all_SAXS(calc_data, weights_vs_ctots[i], concs, oligomersizes)
            chi2_new += get_chi2(Iaverage, Iexp_vs_ctots[i], exp_err_vs_ctots[i])
        chi2_new /= len(ctots)
    
        #Calculate boltzmann factor for metropolis criterion
        criterion = np.exp(-(chi2_new-chi2_old)/T)
    
        #Accept
        if criterion > 1:
            chi2_old = chi2_new
            Ka_old = Ka_new
            accept += 1
    
        #Accept with rng probability or reject
        if criterion <= 1:
            rng = np.random.random()
            if criterion > rng:
                chi2_old = chi2_new
                Ka_old = Ka_new
                accept += 1
            else:
                reject += 1
    
        chi2_lst.append(chi2_old)
        Ka_lst.append(Ka_old)
    
        #Decrease temperature
        T -= T*T_decrease
        j += 1

    acceptance = accept/j
    rejection = reject/j

    #Find minimum chi2
    chi2_min = np.amin(chi2_lst)

    #Find index of minimum chi2 in chi2_lst
    min_index = chi2_lst.index(chi2_min)

    #Find ctot and Ka giving rise to minimum chi2 using index from above
    assert len(Ka_lst) == len(chi2_lst)

    Ka_opt = Ka_lst[min_index]
    
    print("Optimal Ka: %f " % Ka_opt)

    return Ka_opt, acceptance, rejection

# Fitting

Iexp_start_vs_ctots = []
err_exp_start_vs_ctots = []
Iexp_vs_ctots = []
err_exp_vs_ctots = []
weights_vs_ctots = []
phi_limit_reached_vs_ctots = []
chi2_start_vs_ctots = []
theta_vs_ctots = []
Iaverage_vs_ctots = []
offset_vs_ctots = []
scale_vs_ctots = []

#Make array of simulation SAXS data
calc_data = []
for oligomersize in oligomersizes:
    calc_data.append(load_sim_SAXS(f'{calcpath}/calc_data_{oligomersize}mer.dat'))
calc_data = np.array(calc_data)

#Get errors on simulation SAXS data from blocking
err_sim = block_error_SAXS(calc_data)

#Set initial isodesmic and SAXS params
ctot_factor = 1.0
Ka = Ka_init
scale = 1.0
Ka_acceptance = 0
ctot_factor_acceptance = 0

#For getting average chi2 over all experiments
chi2_global=0

#Loop over experiments at different protein concentrations
for ctot in ctots:
    
    #Load experimental data
    exp_file = f'exp_data_{ctot}uM.dat'
    q, Iexp_start, err_exp_start = load_exp_SAXS(exp_file)
    #Normalize for total protein concentration
    Iexp=Iexp_start/ctot
    err_exp=err_exp_start/ctot

    #Generate uniform weights
    weights = uniform_weights(calc_data)

    #Get initial chi2
    concs = isodesmic(Ka, ctot*ctot_factor, oligomersizes)
    Iaverage = volume_average_all_SAXS(calc_data, weights, concs, oligomersizes)
    chi2 = get_chi2(Iaverage, Iexp, err_exp)
    chi2_global += chi2
    print(f"Initial chi2 for {ctot}uM: {chi2}")

    #Start concentration-specific log file and write initial params
    os.system(f'mkdir {ctot}_uM_outfiles')
    with open(f'{ctot}_uM_outfiles/optimization.log', 'w') as f:
        f.write('# Iteration \tChi2 \tKa \tctot \tTheta_min \tScale \tOffset \tPhi_min \n')
        f.write('0 \t%f \t%f \t%f \t%f \t%f \t%f \t%f \n' % (chi2, Ka, ctot*ctot_factor, theta, 1.0, 0.0, 1.0))
    
    #Start warning file for failed reweightings
    with open(f'{ctot}_uM_outfiles/warnings.log', 'w') as f:
        f.write('# List of failed BME reweightings \n')
    
    #Append parameters to lists
    Iexp_start_vs_ctots.append(Iexp_start)
    err_exp_start_vs_ctots.append(err_exp_start)
    Iexp_vs_ctots.append(Iexp)
    err_exp_vs_ctots.append(err_exp)
    weights_vs_ctots.append(weights)
    chi2_start_vs_ctots.append(chi2)
    theta_vs_ctots.append([theta]*len(oligomersizes))
    Iaverage_vs_ctots.append(Iaverage)
    offset_vs_ctots.append(0.0)
    scale_vs_ctots.append(1.0)

#Get average chi2 over all ctot experiments
chi2_global /= len(ctots)

#Make lists numpy arrays
Iexp_start_vs_ctots = np.array(Iexp_start_vs_ctots)
err_exp_start_vs_ctots = np.array(err_exp_start_vs_ctots)
Iexp_vs_ctots = np.array(Iexp_vs_ctots)
err_exp_vs_ctots = np.array(err_exp_vs_ctots)
weights_vs_ctots = np.array(weights_vs_ctots)
chi2_start_vs_ctots = np.array(chi2_start_vs_ctots)
theta_vs_ctots = np.array(theta_vs_ctots)
Iaverage_vs_ctots = np.array(Iaverage_vs_ctots)
offset_vs_Ctots = np.array(offset_vs_ctots)
scale_vs_ctots = np.array(scale_vs_ctots)

#Start global log file and write initial params
with open('optimization_global.log', 'w') as f:
    f.write('#Iteration \tGlobal_chi2 \tKa \tctot_factor \tScale_avg \tKa_accpt \tctot_accpt \n')
    f.write('0 \t%f \t%f \t%f \t%f \t%f \t%f \n' % (chi2_global, Ka, ctot_factor, scale, Ka_acceptance, ctot_factor_acceptance))

#Loop over iterations
for iteration in range(iterations):
    
    print(f'Starting iteration {iteration+1}')  
    print('Fitting isodesmic model and SAXS parameters')
    
    #Fit SAXS scale and offset
    if fit_scale_global_offset_local == True:
        #Fit offset for each ctot
        for i,ctot in enumerate(ctots):
            Iexp_vs_ctots[i], offset_new  = fit_offset_expdata(Iaverage_vs_ctots[i], Iexp_vs_ctots[i], err_exp_vs_ctots[i])
        
        #Fit global scale
        Iexp_vs_ctots, err_exp_vs_ctots, scale_new, r_value = fit_global_scale_expdata(Iaverage_vs_ctots, Iexp_vs_ctots, err_exp_vs_ctots)
        
    if fit_scale_offset_local == True:
        
        for i,ctot in enumerate(ctots):
            #Fit scale and offset, get new exp intensities and errors
            Iexp_vs_ctots[i], err_exp_vs_ctots[i], scale_new, offset_new, r_value = fit_scale_offset_expdata(Iaverage_vs_ctots[i], Iexp_start_vs_ctots[i], err_exp_start_vs_ctots[i])
    
    #Evaluate scale and offset for each ctot
    for i,ctot in enumerate(ctots):
        scale_vs_ctots[i], offset_vs_ctots[i] = evaluate_scale_offset(Iexp_vs_ctots[i], Iexp_start_vs_ctots[i])
    
    #Fit global ctot factor
    if fit_ctot_factor == True:
        ctot_factor, ctot_factor_acceptance, ctot_factor_rejection = montecarlo_opt_ctotfactor(T_init, T_stop, T_decrease, ctots, ctot_factor, factor_MC_step, Ka, weights_vs_ctots, oligomersizes, Iexp_vs_ctots, err_exp_vs_ctots, calc_data)
    
    #Fit global Ka
    if fit_Ka_global == True:
        Ka, Ka_acceptance, Ka_rejection = montecarlo_opt_globalKa(T_init, T_stop, T_decrease, ctots, ctot_factor, Ka_MC_step, Ka, weights_vs_ctots, oligomersizes, Iexp_vs_ctots, err_exp_vs_ctots, calc_data)
    
    #Loop over total protein concentrations (different experiments)
    chi2_global = 0
    for i,ctot in enumerate(ctots):
        
        #Get concentrations from isodesmic model
        concs = isodesmic(Ka, ctot*ctot_factor, oligomersizes)
        
        #Make lists to save phi_eff and chi2 before and after reweighting if last iteration
        if iteration == iterations-1:
            phieff_lst = []
            chi2_before_lst = []
            chi2_after_lst = []

        #Reset to check whether any oligomers go below phi_eff limit
        phi_vs_oligomersizes = [1.0]*len(oligomersizes)
        
        #Check if it's time for reweighting
        if iteration >= iteration_start_reweighting:
            
            print(f'Starting reweighting for {ctot} uM')
            
            #Loop over oligomer species
            for j,oligomersize in enumerate(oligomersizes):

                #Calculate SAXS profile of species
                decompose_SAXS(j, concs, calc_data, weights_vs_ctots[i], oligomersizes, q, Iexp_vs_ctots[i], err_exp_vs_ctots[i], err_sim)

                #Optimize weights with BME
                exp_file = 'exp_file_sim.dat'
                calc_file = f'{calcpath}/calc_data_{oligomersize}mer.dat'
                
                reweighting_check = False
                while reweighting_check == False: 
                    new_weights, chi2_before, chi2_after, phi = run_BME(exp_file, calc_file, theta=theta_vs_ctots[i][j])
                    phi_vs_oligomersizes[j] = phi
                    
                    #Check whether reweighting worked (new_weights not empty)
                    if len(new_weights) == len(calc_data[j]):
                        reweighting_check = True
                    else:
                        #Increase theta slightly
                        theta_vs_ctots[i][j] += theta_vs_ctots[i][j]*theta_decrease*0.1

                        #Write to warning file
                        with open(f'{ctot}_uM_outfiles/warnings.log', 'a') as f:
                            f.write(f'{oligomersize}mer, iteration {iteration}, theta {theta_vs_ctots[i][j]}, phi {phi} \n')

                #Replace old weights with new
                weights_vs_ctots[i][j] = new_weights

                #Save phi_eff and chi2 before and after if last iteration
                if iteration == iterations-1:
                    phieff_lst.append(phi)
                    chi2_before_lst.append(chi2_before)
                    chi2_after_lst.append(chi2_after)

        #Get average intensities with new weights
        Iaverage = volume_average_all_SAXS(calc_data, weights_vs_ctots[i], concs, oligomersizes)
        Iaverage_vs_ctots[i] = Iaverage

        #Get chi2
        chi2 = get_chi2(Iaverage, Iexp_vs_ctots[i], err_exp_vs_ctots[i])
        chi2_global += chi2
        print(f'Chi2 for iteration {iteration+1} at {ctot}uM = {chi2}')
        
        #Write to log file
        with open(f'{ctot}_uM_outfiles/optimization.log', 'a') as f:
            f.write('%i \t%f \t%f \t%f \t%f \t%f \t%f \t%f \n' % (iteration+1, chi2, Ka, ctot*ctot_factor, np.amin(theta_vs_ctots[i]), scale_vs_ctots[i], offset_vs_ctots[i], np.amin(phi_vs_oligomersizes)))
        
        #Check if it's time for reweighting
        if iteration >= iteration_start_reweighting:
            #Update theta if it's higher than the set minimum 
            for j,oligomersize in enumerate(oligomersizes):
                #Decrease theta if minimum theta has not been reached and if phi is not below the limit
                if theta_vs_ctots[i][j] > theta_min and phi_vs_oligomersizes[j] > phi_limit:
                    theta_vs_ctots[i][j] -= theta_vs_ctots[i][j]*theta_decrease
                #Increase theta if phi is below the limit
                elif phi_vs_oligomersizes[j] < phi_limit:
                    theta_vs_ctots[i][j] += theta_vs_ctots[i][j]*theta_decrease
                
                #Check whether theta is below the minimum theta, and if so set it to minimum theta
                if theta_vs_ctots[i][j] < theta_min:
                    theta_vs_ctots[i][j] = theta_min

        #Write output files if it's the last iteration
        if iteration == iterations-1:

            #Write weights files
            for j,oligomersize in enumerate(oligomersizes):
                with open(f'{ctot}_uM_outfiles/weights_%imer.dat' % oligomersize, 'w') as f:
                    for k in range(len(weights_vs_ctots[i][j])):
                        f.write('frame_%i \t%f \n' % (k+1, weights_vs_ctots[i][j][k]))

            #Write summary file
            with open(f'{ctot}_uM_outfiles/summary.log', 'w') as f:
                #Overall stats
                f.write('#Chi2_before \tChi2_after \tKa_before \tKa_after \tCtot_before \tCtot_after \tScale \tOffset \n')
                f.write('%f \t%f \t%f \t%f \t%f \t%f \t%f \t%f \n \n' % (chi2_start_vs_ctots[i], chi2, Ka_init, Ka, ctot, ctot*ctot_factor, scale_vs_ctots[i], offset_vs_ctots[i]))
            
                #Check whether reweighting was performed in final round
                if len(chi2_before_lst) > 0:
                    #BME stats for each oligomer
                    f.write('# Oligomer_size \tChi2_before \tChi2_after \tPhi_eff \tTheta \n')
                    for j,oligomersize in enumerate(oligomersizes):
                        f.write('%i \t%f \t%f \t%f \t%f \n' % (oligomersize, chi2_before_lst[j], chi2_after_lst[j], phieff_lst[j], theta_vs_ctots[i][j]))

    #Get average chi2 over all ctot experiments
    chi2_global /= len(ctots)
    
    print(f'Global chi2: {chi2_global}')
    
    #Write global params to log file
    with open('optimization_global.log', 'a') as f:
        f.write('%i \t%f \t%f \t%f \t%f \t%f \t%f \n' % (iteration+1, chi2_global, Ka, ctot_factor, np.average(scale_vs_ctots), Ka_acceptance, ctot_factor_acceptance))

#rm temperorary exp_file
os.system('rm exp_file_sim.dat')

#Write final average SAXS file
uniform_weights = uniform_weights(calc_data)

for i,ctot in enumerate(ctots):
    
    #Modify calc_data by fitted scale and offset
    scale = scale_vs_ctots[i]
    offset = offset_vs_ctots[i]
    calc_data_scale_offset = calc_data*scale+offset
    
    #Calculate concentrations
    concs = isodesmic(Ka, ctot*ctot_factor, oligomersizes)

    #Calculate Iaverage with and without uniform weights - note that new scale, offset, ctot, Ka fit is used for both    
    Iaverage_noweights = volume_average_all_SAXS(calc_data_scale_offset, uniform_weights, concs, oligomersizes)
    Iaverage = volume_average_all_SAXS(calc_data_scale_offset, weights_vs_ctots[i], concs, oligomersizes)

    #Calculate residuals before and after reweighting
    residuals_noweights = (Iexp_start_vs_ctots[i]-Iaverage_noweights)/err_exp_start_vs_ctots[i]
    residuals_weights = (Iexp_start_vs_ctots[i]-Iaverage)/err_exp_start_vs_ctots[i]

    #Write SAXS file
    with open(f'{ctot}_uM_outfiles/SAXS_Iaverage.dat', 'w') as f:
        f.write('#q \t I_exp \t I_exp_err \t I_avg_noweights \t Residuals_noweights \t I_avg_rew \t Residuals_rew \n')
        for j,q_value in enumerate(q):
            f.write('%e \t%e \t%e \t%e \t%e \t%e \t%e \n' % (q_value, Iexp_start_vs_ctots[i][j], err_exp_start_vs_ctots[i][j], Iaverage_noweights[j], residuals_noweights[j], Iaverage[j], residuals_weights[j]))
          

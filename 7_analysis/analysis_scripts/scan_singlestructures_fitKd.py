import numpy as np
from sklearn.linear_model import LinearRegression
import pickle as pkl
from scipy import stats

iterations = 10000

#Kd values to scan
Kd_scan = np.logspace(-2, 2, num=10000)

MATH_type='MATHfree'
oligomersizes = np.arange(2, 61, 2)
calcpath = f'calc_data_{MATH_type}'
concs = [5.0, 10.0, 20.0, 30.0, 40.0]

def save_pickle(filename, pickle_obj):
    with open(filename, 'wb') as f:
        pkl.dump(pickle_obj, f)

def load_pickle(filename):
    with open(filename, 'rb') as f:
        loaded_obj = pkl.load(f)
        
    return loaded_obj

#Loads experimental SAXS q, I and err - in BME format
def load_exp_SAXS(exp_file):
    q, I, err = np.genfromtxt(exp_file, usecols=(0,1,2), skip_header=0, unpack=True)
    return np.array(q), np.array(I), np.array(err)

#Loads SAXS intensities calculated from simulations - in BME format
def load_sim_SAXS(calc_file):
    calc_data = np.array([[float(x) for x in line.split()[1:]] for line in open(calc_file) if("#" not in line)])
    return np.array(calc_data)

#Fit calc data profiles scale and offset with linear regression
def fit_scale_offset(Iaverage_sim, I_exp, err_exp):
    #Get weight for each point based on exp error
    sample_weight=1.0/(err_exp**2)
    
    #Linear regression
    reg = LinearRegression(fit_intercept=True).fit(Iaverage_sim.reshape(-1,1),I_exp.reshape(-1,1),sample_weight=sample_weight)
    r_value = reg.score(Iaverage_sim.reshape(-1,1),I_exp.reshape(-1,1),sample_weight=sample_weight)
    slope,intercept = reg.coef_[0],reg.intercept_
    
    Iaverage_sim_fit = Iaverage_sim*slope+intercept
    
    return Iaverage_sim_fit, slope, intercept, r_value

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

def get_chi2(Iaverage_sim, I_exp, err_exp):
    
    chi2 = np.average(np.square((Iaverage_sim-I_exp)/err_exp))
    
    return chi2

#Make array of simulation SAXS data
calc_data = []
for oligomersize in oligomersizes:
    calc_data.append(load_sim_SAXS(f'{calcpath}/calc_data_{oligomersize}mer.dat'))
calc_data = np.array(calc_data)
nr_frames = len(calc_data[0])

#Get experimental SAXS
Iexp_vs_concs = []
err_exp_vs_concs = []
for conc in concs:
    exp_file = f'WT_{MATH_type}_scaleoffsetonly/exp_data_{conc}uM.dat'
    q, Iexp, err_exp = load_exp_SAXS(exp_file)
    Iexp_vs_concs.append(Iexp)
    err_exp_vs_concs.append(err_exp)

chi2_global_vs_iterations = []
Ka_fitted_vs_iterations = []

#Iterate over single structure selection
for iteration in range(iterations):

    print(f'Iteration {iteration}')

    #Select a single random frame for each oligomer and get SAXS profile
    sel_frames = np.random.randint(0, high=nr_frames, size=len(oligomersizes))
    sel_SAXS = calc_data[range(len(oligomersizes)), sel_frames]

    #Scan Kd values for given set of single structures
    chi2_global_vs_Kd_scan = []
    for Kd in Kd_scan:
        Ka = 1/Kd

        chi2_global = 0
        for i,conc in enumerate(concs):

            #Get concentrations from isodesmic model
            isodesmic_concs = isodesmic(Ka, conc, oligomersizes)
            volume_weights = isodesmic_concs*oligomersizes
            volume_weights = volume_weights/np.sum(volume_weights)

            #Calculate average SAXS
            Iaverage = np.sum(sel_SAXS*volume_weights[:,np.newaxis],axis=0)

            #Fit SAXS scale and cst
            Iaverage_fit, slope, intercept, r_value = fit_scale_offset(Iaverage, Iexp_vs_concs[i], err_exp_vs_concs[i])

            #Calculate chi2 and add to running sum
            chi2 = get_chi2(Iaverage_fit, Iexp_vs_concs[i], err_exp_vs_concs[i])
            chi2_global += chi2

        #Average global chi2 and append
        chi2_global /= len(concs)
        chi2_global_vs_Kd_scan.append(chi2_global)

    chi2_global_vs_Kd_scan = np.array(chi2_global_vs_Kd_scan)

    #Select best fit Kd and corresponding chi2
    min_index = np.argmin(chi2_global_vs_Kd_scan)
    Kd_fitted = Kd_scan[min_index]
    Ka_fitted = 1/Kd_fitted
    chi2_min = chi2_global_vs_Kd_scan[min_index]

    chi2_global_vs_iterations.append(chi2_min)
    Ka_fitted_vs_iterations.append(Ka_fitted)

    if iteration == 0 or chi2_global < chi2_global_min:
        chi2_global_min = chi2_global
        frames_min = sel_frames

    print(chi2_min)

chi2_global_vs_iterations = np.array(chi2_global_vs_iterations)
Ka_fitted_vs_iterations = np.array(Ka_fitted_vs_iterations)

scan_singlestructure_results = {'chi2':chi2_global_vs_iterations, 'frames_min':frames_min, 'Ka_fit':Ka_fitted_vs_iterations}
save_pickle(f'pickles/scan_singlestructure_results_{MATH_type}_fitKa.pkl', scan_singlestructure_results)


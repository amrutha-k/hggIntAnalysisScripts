##### import ROOT
import os
import numpy as np
import math
import uproot as up
import pandas as pd
import matplotlib.pyplot as plt
import mplhep as mpl
from scipy.optimize import minimize
from iminuit import Minuit
import numdifftools as nd
import argparse as ap
from collections import OrderedDict as od

from python.model import dscb, splusb
import python.plotting as plottool

mpl.style.use("ROOT")
#mpl.style.use("CMS")

parser = ap.ArgumentParser()

parser.add_argument("-fsig", "--inputSigFile", nargs='+', required=True,
                    help='input signal file ')
parser.add_argument("-fintgg", "--inputIntggFile", nargs='+', required=True,
                    help='input gg interference file ')
parser.add_argument("-fintqg", "--inputIntqgFile", nargs='+', required=True,
                    help='input qg interference file ')
parser.add_argument("-o", "--outFile", required=True,
                    help='output txt file name')
parser.add_argument("-outdir", "--outDir", required=True,
                    help='output directory path')
parser.add_argument("-y", "--year", required=True,
                    help='year')
#parser.add_argument("-m", "--mass", required=True,
#                    help='mass')
parser.add_argument("-c", "--categ", nargs='+', required=True,
                    help='category')
args = parser.parse_args()

## GLOBAL VARIABLES       
GammaRatio = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
#GammaRatio = [0.5, 1.0, 2.0, 5.0, 10.0]
#GammaRatio = [0.1, 0.3, 0.5, 0.7, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 3.0, 4.0, 5.0, 7.0, 9.0]
n = len(GammaRatio)
ncats = 11

def get_fixed_param(categ_):
    '''
    GammaRatio: aL, aR, nL, nR, norm, mean, sigma
    '''
    with open(f'./fixed_params_txts/{args.year}/category{categ_}.txt', 'r') as file:
        fp_dict = {}
        for line in file:
            if line.strip():
                k, values = line.strip().split(':')
                k = k.strip()
                values_list = [float(v.strip()) for v in values.split(',')]
                if k !='noint' and k !='alpha':
                    k = float(k)
                    if k not in fp_dict:
                        fp_dict[k] = values_list
                elif(k == 'noint'):
                    fp_dict['noint'] = values_list
                else:
                    fp_dict['alpha'] = values_list
    return fp_dict

mh = 125.0
initial_params = []

# parameters look-up-table         
pLUT = od()
pLUT['noint'] = od()
pLUT['int'] = od()
pLUT['noint']['alpha_p1'] = [0.0,-0.08,0.08]
#pLUT['noint']['alpha_p1_m125'] = [0.0,-0.05,0.05]
#pLUT['noint']['alpha_p1_m130'] = [0.0,-0.05,0.05]
pLUT['noint']['m_p1'] = [0.0,-0.1,0.1]
#pLUT['noint']['m_p1_m125'] = [0.0,-0.08,0.08]
#pLUT['noint']['m_p1_m130'] = [0.0,-0.08,0.08]                                                                                                                                               
pLUT['noint']['sigma_p1'] = [0.0,-0.1,0.1]
#pLUT['noint']['sigma_p1_m125'] = [0.0,-0.1,0.1]
#pLUT['noint']['sigma_p1_m130'] = [0.0,-0.1,0.1]
pLUT['noint']['n1_p1'] = [0.0,-0.1,0.1]
#pLUT['noint']['n1_p1_m125'] = [0.0,-1.5,1.5]
#pLUT['noint']['n1_p1_m130'] = [0.0,-1.5,1.5]
pLUT['noint']['n2_p1'] = [0.0,-0.1,0.1]
#pLUT['noint']['n2_p1_m125'] = [0.0,-1.5,1.5]
#pLUT['noint']['n2_p1_m130'] = [0.0,-1.5,1.5]
pLUT['noint']['a1_p1'] = [0.0,-0.1,0.1]
#pLUT['noint']['a1_p1_m125'] = [0.0,-1.5,1.5]
#pLUT['noint']['a1_p1_m130'] = [0.0,-1.5,1.5]
pLUT['noint']['a2_p1'] = [0.0,-0.1,0.1]
#pLUT['noint']['a2_p1_m125'] = [0.0,-1.5,1.5]
#pLUT['noint']['a2_p1_m130'] = [0.0,-1.5,1.5]
pLUT['noint']['N_p1'] = [0.0,-0.1,0.1]
#pLUT['noint']['N_p1_m125'] = [0.0,-0.1,0.1]
#pLUT['noint']['N_p1_m130'] = [0.0,-0.1,0.1]
pLUT['int']['sigma_int_p1'] = [0.0,-0.1,0.1]
#pLUT['int']['sigma_int_p1_m125'] = [0.0,-0.1,0.1]
#pLUT['int']['sigma_int_p1_m130'] = [0.0,-0.1,0.1]
pLUT['int']['n1_int_p1'] = [0.0,-1.0,1.0]
#pLUT['int']['n1_int_p1_m125'] = [0.0,-2.5,2.5]
#pLUT['int']['n1_int_p1_m130'] = [0.0,-2.5,2.5]
pLUT['int']['n2_int_p1'] = [0.0,-0.01,0.01]
#pLUT['int']['n2_int_p1_m125'] = [0.0,-2.5,2.5]
#pLUT['int']['n2_int_p1_m130'] = [0.0,-2.5,2.5]
pLUT['int']['a1_int_p1'] = [0.0,-0.5,0.5]
#pLUT['int']['a1_int_p1_m125'] = [0.0,-2.5,2.5]
#pLUT['int']['a1_int_p1_m130'] = [0.0,-2.5,2.5]
pLUT['int']['a2_int_p1'] = [0.0,-0.01,0.01]
#pLUT['int']['a2_int_p1_m125'] = [0.0,-2.5,2.5]
#pLUT['int']['a2_int_p1_m130'] = [0.0,-2.5,2.5]
pLUT['int']['N_int_p1'] = [0.0,-0.1,0.1]
#pLUT['int']['N_int_p1_m125'] = [0.0,-0.5,0.5]
#pLUT['int']['N_int_p1_m130'] = [0.0,-0.5,0.5]


# Luminosity map in fb^-1                                                                                                                                                                                  
lumiMap = {'2016':36.33, '2016preVFP': 19.52, '2016postVFP': 16.81, '2017':41.48, '2018':63.67, 'combined':138}
lumi = lumiMap[args.year]
# 3 mass points
hint_array = np.empty((3, n, ncats), dtype=object)
hint_err_array = np.empty((3, n, ncats), dtype=object)
filtered_hint_array = np.empty((3, n, ncats), dtype=object)
filtered_hint_err_array = np.empty((3, n, ncats), dtype=object)
filtered_bin_centers_array = np.empty((3, n, ncats), dtype=object)
hnoint_array = np.empty((3, ncats), dtype=object)
hnoint_err_array = np.empty((3, ncats), dtype=object)
filtered_hnoint_array = np.empty((3, ncats), dtype=object)
filtered_hnoint_err_array = np.empty((3, ncats), dtype=object)
bin_centers = []

# to initialize the normalization parameter
int_norm = np.empty((3, ncats))
noint_norm = np.empty((3, ncats))
#int_norm = np.zeros(11)
#noint_norm = np.zeros(11)
#s_int_array = [0.95, 1.1, 0.7, 0.7, 0.9, 1.0, 1.0, 1.2, 0.9, 1.2, 1.22]

ar_hint = None
ar_hint_err = None
ar_hnoint = None
ar_hnoint_err = None
x_int = None
x_noint = None
bkg_func = None

alphaVars = []
hlow = 118.0
hhigh = 130.0
nbins = 120
catlabel = [f'UntaggedTag_{i}_{args.year}' for i in range(ncats - 1)] + [f'VBFTag_0_{args.year}']
m_dict = {
    120.0:(52.22, 0.002218, 0.000517, 1.35, -2.168e-04),
    125.0:(48.58, 0.00227, 0.000517, 1.26, -2.051e-04),
    130.0:(45.31, 0.002238, 0.000517, 1.14, -2.067e-04)
    }

fixed_params_dict = {}


def get_alpha_err(best_fit_params):
    chi2alphap = []
    chi2alpham = []
    variations = [0.0, 0.04, 0.08, 0.12, 0.16, 0.2, 0.25, 0.3, .35, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7]
    global alphaVars
    #alpha = fixed_params_dict['alpha'][0] + mref * best_fit_params[im]
    best_fit_list = list(best_fit_params)
    best_fit_list_copy = best_fit_list.copy()
    alphaVars = [(best_fit_params[0]*v) for v in variations]
    for i in alphaVars:
        best_fit_list_copy[0] = best_fit_list[0] + i
        chi2alphap.append(global_chi2(best_fit_list_copy)[0])
        best_fit_list_copy[0] = best_fit_list[0] - i
        chi2alpham.append(global_chi2(best_fit_list_copy)[0])
    return chi2alphap, chi2alpham

def get_chi2_mass(h, herr, x_, a1_, a2_, n1_, n2_, N_, M, sigma_):
    #y_fit = splusb(x_, bfunc, bscale, a1_, a2_, n1_, n2_, N_, M, sigma_)
    y_fit = dscb(x_, a1_, a2_, n1_, n2_, N_, M, sigma_)
    if np.isnan(y_fit).any() or np.isinf(y_fit).any():
        indices = np.where(np.isnan(y_fit) | np.isinf(y_fit))
        #print(indices)
        for i in indices[0]:
            global initial_params
            y_fit[i] = initial_params[i]
    #print(herr)
    res = (h - np.array(y_fit)) / (herr**0.5)
    #if np.isnan(res).any() or np.isinf(res).any():
    #    res = 99999.99
    chi2 = np.sum(res**2)
    return chi2

def get_params(im, mp, params):
    mref = mp - mh
    #print(mref)
    global fixed_params_dict

    '''
    alpha = fixed_params_dict['alpha'][0] + mref * params[im]
    Mnoint_fit = fixed_params_dict['noint'][5] + mref * params[im + 3]
    sigma_noint = fixed_params_dict['noint'][6] + mref * params[im + 6]
    n1_noint = fixed_params_dict['noint'][2] + mref * params[im + 9]
    n2_noint = fixed_params_dict['noint'][3] + mref * params[im + 12]
    a1_noint = fixed_params_dict['noint'][0] + mref * params[im + 15]
    a2_noint = fixed_params_dict['noint'][1] + mref * params[im + 18]
    N_noint = fixed_params_dict['noint'][4] + mref * params[im + 21]

    sigma_int = np.array([fixed_params_dict[gam][6] for gam in GammaRatio]) + mref * np.array(params[n*im+24: n + n*im+24])
    n1_int = np.array([fixed_params_dict[gam][2] for gam in GammaRatio]) + mref * np.array(params[3*n + n*im+24: 4*n + n*im+24])
    n2_int = np.array([fixed_params_dict[gam][3] for gam in GammaRatio]) + mref * np.array(params[6*n + n*im+24: 7*n + n*im+24])
    a1_int = np.array([fixed_params_dict[gam][0] for gam in GammaRatio]) + mref * np.array(params[9*n + n*im+24: 10*n + n*im+24])
    a2_int = np.array([fixed_params_dict[gam][1] for gam in GammaRatio]) + mref * np.array(params[12*n + n*im+24: 13*n + n*im+24])
    N_int = np.array([fixed_params_dict[gam][4] for gam in GammaRatio]) + mref * np.array(params[15*n + n*im+24: 16*n + n*im+24])
    '''
    alpha = fixed_params_dict['alpha'][0] + mref * params[0]
    Mnoint_fit = fixed_params_dict['noint'][5] + mref * params[1]
    sigma_noint = fixed_params_dict['noint'][6] + mref * params[2]
    n1_noint = fixed_params_dict['noint'][2] + mref * params[3]
    n2_noint = fixed_params_dict['noint'][3] + mref * params[4]
    a1_noint = fixed_params_dict['noint'][0] + mref * params[5]
    a2_noint = fixed_params_dict['noint'][1] + mref * params[6]
    N_noint = fixed_params_dict['noint'][4] + mref * params[7]

    sigma_int = np.array([fixed_params_dict[gam][6] for gam in GammaRatio]) + mref * np.array(params[8:n+8])
    n1_int = np.array([fixed_params_dict[gam][2] for gam in GammaRatio]) + mref * np.array(params[n+8: 2*n+8])
    n2_int = np.array([fixed_params_dict[gam][3] for gam in GammaRatio]) + mref * np.array(params[2*n+8: 3*n+8])
    a1_int = np.array([fixed_params_dict[gam][0] for gam in GammaRatio]) + mref * np.array(params[3*n+8: 4*n+8])
    #print('a1_int', fixed_params_dict[0.1][0])
    a2_int = np.array([fixed_params_dict[gam][1] for gam in GammaRatio]) + mref * np.array(params[4*n+8: 5*n+8])
    N_int = np.array([fixed_params_dict[gam][4] for gam in GammaRatio]) + mref * np.array(params[5*n+8: 6*n+8])
    return a1_noint, a2_noint, n1_noint, n2_noint, N_noint, Mnoint_fit, sigma_noint, a1_int, a2_int, n1_int, n2_int, N_int, sigma_int, alpha
    #return alpha, Mnoint_fit, sigma_noint, n1_noint, n2_noint, a1_noint, a2_noint, N_noint, sigma_int, n1_int, n2_int, a1_int, a2_int, N_int 


def global_chi2(params):
    chi2_int_total = 0
    chi2_noint_total = 0
    chi2_total = 0
    dof_int = 0
    for im, mp in enumerate(m_dict.keys()):
        global ar_hint, ar_hint_err, ar_hnoint, ar_hnoint_err, x_int, x_noint
        a1_noint, a2_noint, n1_noint, n2_noint, N_noint, Mnoint_fit, sigma_noint, a1_int, a2_int, n1_int, n2_int, N_int, sigma_int, alpha = get_params(im,mp,params)
        for g, gamma in enumerate(GammaRatio):
            chi2_int = get_chi2_mass(ar_hint[im][g], ar_hint_err[im][g], x_int[im][g], a1_int[g], a2_int[g], n1_int[g], n2_int[g], N_int[g], Mnoint_fit+alpha*np.sqrt(gamma), sigma_int[g])
            chi2_int_total += chi2_int
            dof_int += len(x_int[im][g])
        
        chi2_noint = get_chi2_mass(ar_hnoint[im], ar_hnoint_err[im], x_noint, a1_noint, a2_noint, n1_noint, n2_noint, N_noint, Mnoint_fit, sigma_noint)
        chi2_noint_total += chi2_noint
    chi2_total = chi2_int_total + chi2_noint_total
    dof_noint = 3*len(x_noint) - 7
    dof_total = (dof_int - n*6) + dof_noint - 1
    return chi2_total, dof_total

## MAIN
def main():
    files_intgg = args.inputIntggFile
    files_intqg = args.inputIntqgFile
    files_sig = args.inputSigFile

    gg_int_events = []
    qg_int_events = []
    no_int_events = []

    for i in range(ncats): 
        gg_int_dfs = []
        qg_int_dfs = []
        no_int_dfs = []
        for f in range(len(files_intgg)):
            if i == (ncats - 1):
                df_vbf = up.open(files_intgg[f])['tagsDumper/trees/ggh_125_13TeV_VBFTag_0'].arrays(['CMS_hgg_mass','weight'], library='pd')
                df_vbf['weight'] = df_vbf['weight']*(lumi)*(0.3215174129)
                df_noInt_vbf = up.open(files_sig[f])['tagsDumper/trees/ggh_125_13TeV_VBFTag_0'].arrays(['CMS_hgg_mass','weight'], library='pd')
                df_noInt_vbf['weight'] = df_noInt_vbf['weight']*(lumi)
                gg_int_dfs.append(df_vbf)
                no_int_dfs.append(df_noInt_vbf)
                if f == 0:
                    df_qg_vbf = up.open(files_intqg[f])['tagsDumper/trees/ggh_125_13TeV_VBFTag_0'].arrays(['CMS_hgg_mass','weight'], library='pd')
                    df_qg_vbf['weight'] = df_qg_vbf['weight']*(lumi)
                    qg_int_dfs.append(df_qg_vbf)
            else:
                df = up.open(files_intgg[f])[f'tagsDumper/trees/ggh_125_13TeV_UntaggedTag_%d'%(i)].arrays(['CMS_hgg_mass','weight'], library='pd')
                df['weight'] = df['weight']*(lumi)*(0.3215174129)
                gg_int_dfs.append(df)
                df_noInt = up.open(files_sig[f])[f'tagsDumper/trees/ggh_125_13TeV_UntaggedTag_%d'%(i)].arrays(['CMS_hgg_mass','weight'], library='pd')
                df_noInt['weight'] = df_noInt['weight']*(lumi)
                no_int_dfs.append(df_noInt)
                if f == 0:
                    df_qg = up.open(files_intqg[f])[f'tagsDumper/trees/ggh_125_13TeV_UntaggedTag_%d'%(i)].arrays(['CMS_hgg_mass','weight'], library='pd')
                    df_qg['weight'] = df_qg['weight']*(lumi)
                    qg_int_dfs.append(df_qg)

        
        #merge categories if there are multiple files    
        gg_int_events.append(pd.concat(gg_int_dfs, ignore_index=True))
        no_int_events.append(pd.concat(no_int_dfs, ignore_index=True))
        qg_int_events.append(pd.concat(qg_int_dfs, ignore_index=True))
        #remove 0 weights (present in samples Ruben produced)
        no_int_events[i] = no_int_events[i].loc[no_int_events[i]['weight'] != 0]
        
    for m,mp in enumerate(m_dict.keys()):
        m_ref = m_dict[125.0]
        gg_int_dfs_mp = [df.copy() for df in gg_int_events]
        qg_int_dfs_mp = [df.copy() for df in qg_int_events]
        no_int_dfs_mp = [df.copy() for df in no_int_events]
        for d in range(len(gg_int_dfs_mp)):
            gg_int_dfs_mp[d]['weight'] = gg_int_dfs_mp[d]['weight']*(m_dict[mp][2]/m_ref[2])*(m_dict[mp][3]/m_ref[3])
            qg_int_dfs_mp[d]['weight'] = qg_int_dfs_mp[d]['weight']*(m_dict[mp][4]/m_ref[4])
            no_int_dfs_mp[d]['weight'] = no_int_dfs_mp[d]['weight']*(m_dict[mp][0]/m_ref[0])*(m_dict[mp][1]/m_ref[1])
            
        for j in range(len(GammaRatio)):
            gg_int_dfs_scaled = [df.copy() for df in gg_int_dfs_mp]
            qg_int_dfs_scaled = [df.copy() for df in qg_int_dfs_mp]
            for k in range(len(gg_int_dfs_scaled)): 
                gg_int_dfs_scaled[k]['weight'] = gg_int_dfs_scaled[k]['weight']*np.sqrt(GammaRatio[j])
                qg_int_dfs_scaled[k]['weight'] = qg_int_dfs_scaled[k]['weight']*np.sqrt(GammaRatio[j])
            df_cats = [pd.concat([df_int, df_intqg, df_noint], ignore_index=True) for df_int, df_intqg, df_noint in zip(gg_int_dfs_scaled, qg_int_dfs_scaled, no_int_dfs_mp)]
            #df_cats = [pd.concat([df_int, df_noint], ignore_index=True) for df_int, df_noint in zip(gg_int_dfs_scaled, no_int_dfs)]
            for c, df_cat in enumerate(df_cats):
                if j==0:
                    h_, bin_edges = np.histogram(no_int_events[c]['CMS_hgg_mass'], bins=nbins, range=(hlow,hhigh), weights=no_int_events[c]['weight'])
                    herr_, bin_edges = np.histogram(no_int_events[c]['CMS_hgg_mass'], bins=nbins, range=(hlow,hhigh), weights=no_int_events[c]['weight']**2)
                    global bin_centers
                    bin_centers = 0.5*(bin_edges[1:] + bin_edges[:-1])
                    hnoint_array[m,c] = h_ 
                    hnoint_err_array[m,c] = np.sqrt(herr_) 
                    noint_norm[m,c] = float(h_.max())     


                h_, bin_edges = np.histogram(df_cat['CMS_hgg_mass'], bins=nbins, range=(hlow, hhigh), weights=df_cat['weight'])
                herr_, bin_edges = np.histogram(df_cat['CMS_hgg_mass'], bins=nbins, range=(hlow, hhigh), weights=df_cat['weight']**2)

                hint_array[m,j,c] = h_
                hint_err_array[m,j,c] = np.sqrt(herr_) 
                last_positive_bin = nbins - 1
                first_positive_bin = 0
                ##fit only the positive part of the spectrum
                midpoint_h_ = int(0.5*len(h_))
                for ibin in range(midpoint_h_):
                    #if (h_[ibin]<0 and h_[ibin+1]<0):
                    #    first_positive_bin = ibin + 1
                    if(ibin + midpoint_h_ + 1 == len(h_)):
                        continue
                    else:
                        if (h_[ibin + midpoint_h_]<0 and h_[ibin + midpoint_h_ + 1]<0):
                            last_positive_bin = (ibin + midpoint_h_) - 1
                            break

                
                mask = (bin_centers >= bin_centers[first_positive_bin]) & (bin_centers <= bin_centers[last_positive_bin])
                filtered_bin_centers = bin_centers[mask]
                filtered_h_ = (h_)[mask]
                #filtered_h_ = (h_ + bkg_h[c])[mask]
                #filtered_herr_ = (np.sqrt(herr_ + bkg_h[c]))[mask]   
                filtered_herr_ = (np.sqrt(herr_))[mask] 
                int_norm[m,c] = float(h_.max())

                filtered_hint_array[m,j,c] = filtered_h_
                filtered_hint_err_array[m,j,c] = filtered_herr_
                filtered_bin_centers_array[m,j,c] = filtered_bin_centers

    ## MINIMIZER
    if (args.categ[0] == 'all'):
        categories = [icat for icat in range(ncats)]
    else:
        categories = [int(icat) for icat in args.categ]

    for ct in categories:
        print(f'CATEGORY {ct}')
        outfile = open(f'{args.outFile}_cat{ct}.txt', 'w')

        global n, fixed_params_dict
        fixed_params_dict = get_fixed_param(ct)

        # Set up parameter names
        param_names = list(pLUT['noint'].keys()) + [f'{name}[{i}]' for name in pLUT['int'].keys() for i in range(n)]
        #print(param_names)
        # Initial guesses for parameters
        global initial_params
        initial_params = [vnoint[0] for vnoint in pLUT['noint'].values()] + [vint[0] for vint in pLUT['int'].values() for _ in range(n)]
        #print(initial_params)
        #parameter bounds
        lim = [(vnoint[1],vnoint[2]) for vnoint in pLUT['noint'].values()]  + [(vint[1],vint[2]) for vint in pLUT['int'].values() for _ in range(n)] 
                

        global ar_hint, ar_hint_err, ar_hnoint, ar_hnoint_err, x_int, x_noint
        ar_hint = filtered_hint_array[:,:,ct]
        ar_hint_err = filtered_hint_err_array[:,:,ct]
        ar_hnoint = hnoint_array[:,ct]
        ar_hnoint_err = hnoint_err_array[:,ct]
        x_int=filtered_bin_centers_array[:,:,ct]
        x_noint = bin_centers

        InitialFit = minimize(lambda parameters: global_chi2(parameters)[0], initial_params, bounds=lim, method='Powell')
        #for im, mp in enumerate(m_dict.keys()):
        #    initial_fit_results = get_params(im, mp, InitialFit.x)
        #    for x in initial_fit_results:
        #        if x.any()<0.0:
        #           print(x)

        # Initialize Minuit
        print("initializing Minuit...")
        if np.isnan(InitialFit.x).any() or np.isinf(InitialFit.x).any():
            m = Minuit(lambda parameters: global_chi2(parameters)[0], initial_params, name=param_names)
        else:
            m = Minuit(lambda parameters: global_chi2(parameters)[0], InitialFit.x, name=param_names)
        #m1 = Minuit(global_chi2, initial_params, name=param_names)
        m.errordef = Minuit.LEAST_SQUARES
        m.print_level = 1
        #m.tol = 0.1
        for ip,p in enumerate(param_names):
            m.limits[p] = lim[ip]

        # Perform minimization
        print("performing minimization...")
        m.simplex().migrad(ncall=100000)
        print("Fitted parameters:")
        for name, value, error in zip(m.parameters, m.values, m.errors):
            print(f"{name} = {value:.4f} ± {error:.4f}")   

        ## plot
        
        dof = global_chi2(m.values)[1]
        fit_params_err = [0.00 for _ in range(7)]
        for im, mp in enumerate(m_dict.keys()):
            fit_params = get_params(im, mp, m.values)
            fit_params_int = np.empty((6,n))
            for g,gam in enumerate(GammaRatio):
                #fit_params_err = [0.00 for _ in range(7)]
                Mint = fit_params[5]+fit_params[-1]*np.sqrt(gam)
                fit_params_int = np.vstack(fit_params[7:-1])
                fit_params_int = np.insert(fit_params_int,-1, Mint, axis=0)
                label = f'mgg_int_gamma{gam}_cat{ct}'
                

                #plt_hfit = plottool.plot_data(hint_array[:,ct][g], hint_err_array[:,ct][g], bin_centers, x_int, fit_params, fit_params_err, bkg_func, bscale, m.fval/dof, args.outDir, label)
                plt_hfit = plottool.plot_data(hint_array[im,:,ct][g], hint_err_array[im,:,ct][g], x_noint, x_int[im,g], fit_params_int[:,g], fit_params_err, m.fval/dof, args.outDir, mp, label)
            label = f'mgg_noint_cat{ct}'
            plt_hfit_noint = plottool.plot_data(ar_hnoint[im], ar_hnoint_err[im], x_noint, x_noint, fit_params[0:7], fit_params_err, m.fval/dof, args.outDir, mp, label)
            #plt_dm_vs_gr = plottool.plot_dm(fit_params[-1], GammaRatio, args.outDir, f'dm_vs_gr_cat{ct}')
            #plt_dm_vs_gr = plottool.plot_dm(best_fit_params, param_errors, GammaRatio, args.outDir, f'dm_vs_gr_cat{ct}')

        # get alpha error
        chi2_ap, chi2_am = get_alpha_err(m.values)
        a1err = plottool.plot_a_err(m.values[0], alphaVars, chi2_ap, chi2_am, args.outDir, f'alphaErrCat{ct}')
        #best_alpha = fixed_params_dict['alpha'][0] + (mp-mh)*m.values[im]
        #if mp-mh == 0:
        #    alphaErr = fixed_params_dict['alpha'][1]
        #    abErr = 0.0
        #else:
        #   abErr = plottool.plot_a_err(m.values[im], alphaVars, chi2_ap, chi2_am, args.outDir, mp, f'alphaErrCat{ct}')
        #    alphaErr = np.sqrt(fixed_params_dict['alpha'][1]**2 + (mp-mh)**2 * abErr**2)
        #outfile.write(f'M{int(mp)}: {best_alpha:.5f}\n')
        outfile.write(f"{fixed_params_dict['alpha'][0]:.5f}, {m.values[0]:.7f}, {fixed_params_dict['alpha'][1]:.7f}, {a1err:.7f}\n")
        outfile.close()

if __name__ == "__main__":
    main()     

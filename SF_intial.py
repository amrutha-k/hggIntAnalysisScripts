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

from python.model import dscb, splusb
import python.plotting_ini as plottool

mpl.style.use("ROOT")


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
parser.add_argument("-m", "--mass", required=True,
                    help='mass')
parser.add_argument("-c", "--categ", nargs='+', required=True,
                    help='category')
args = parser.parse_args()

## GLOBAL VARIABLES       
#GammaRatio = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 15.0, 20.0, 30.0, 40.0]
GammaRatio = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.1, 1.2, 1.4, 1.6, 1.8, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
#GammaRatio = [0.1, 0.3, 0.5, 0.7, 0.9, 1.0, 1.1, 1.3, 1.5, 1.7, 1.9, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]
n = len(GammaRatio)
hlow = 118.0
hhigh = 130.0
nbins = 120
#bscale = (hhigh - hlow) / nbins

ncats = 11
# Luminosity map in fb^-1                                                                                                                                                                                  
lumiMap = {'2016':36.33, '2016preVFP': 19.52, '2016postVFP': 16.81, '2017':41.48, '2018':63.67, 'combined':138}
lumi = lumiMap[args.year]

hint_array = np.empty((n, ncats), dtype=object)
hint_err_array = np.zeros((n, ncats), dtype=object)
filtered_hint_array = np.zeros((n, ncats), dtype=object)
filtered_hint_err_array = np.zeros((n, ncats), dtype=object)
filtered_bin_centers_array = np.empty((n, ncats), dtype=object)
bkg_curve_array = np.zeros(ncats, dtype=object)
bkg_h = np.zeros(ncats, dtype=object)
hnoint_array = np.zeros(ncats, dtype=object)
hnoint_err_array = np.zeros(ncats, dtype=object)
bin_centers = []

# to initialize the normalization parameter
int_norm = np.zeros(ncats)
noint_norm = np.zeros(ncats)
#int_norm = np.zeros(11)
#noint_norm = np.zeros(11)
#s_int_array = [0.94, 1.1, 1.0, 1.2, 1.0, 1.24, 1.1, 1.2, 1.2, 1.35, 1.25]
s_int_array = [0.95, 1.1, 0.7, 0.7, 0.9, 1.0, 1.0, 1.2, 0.9, 1.2, 1.22]

ar_hint = None
ar_hint_err = None
ar_hnoint = None
ar_hnoint_err = None
x_int = None
x_noint = None
bkg_func = None

alphaVars = []


def get_alpha_err(best_fit_params):
    chi2alphap = []
    chi2alpham = []
    variations = [0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.12]
    global alphaVars
    #if abs(best_fit_params[0]) < 10e-4:
    #    best_fit_params[0] = -0.010
        
    alphaVars = [(best_fit_params[0]*v) for v in variations]
    #print(alphaVars)
    best_mnoint = best_fit_params[6 * n + 1]
    best_mint = best_mnoint + best_fit_params[0]*np.sqrt(GammaRatio)
    best_fit_list = list(best_fit_params)
    best_fit_list_copy = best_fit_list.copy()
    for i in alphaVars:
        best_fit_list_copy[0] = best_fit_list[0] + i
        chi2alphap.append(global_chi2(best_fit_list_copy)[0])
        best_fit_list_copy[0] = best_fit_list[0] - i
        chi2alpham.append(global_chi2(best_fit_list_copy)[0])
    return chi2alphap, chi2alpham

def get_chi2_mass(h, herr, x_, a1_, a2_, n1_, n2_, N_, M, sigma_):
    y_fit = dscb(x_, a1_, a2_, n1_, n2_, N_, M, sigma_)
    if np.any(np.isinf(y_fit)):
        return np.inf
    res = (h - np.array(y_fit)) / (herr**0.5)
    chi2 = np.sum(res**2)
    return chi2

def global_chi2(params):
    alpha = params[0]
    # Extract DSCB parameters for int and noint histograms
    global n
    #Mint_fit = np.array(params[1:n + 1])  # Mint for each gamma
    sigma_int = np.array(params[1:n + 1])
    a1_int = np.array(params[n + 1:2 * n + 1])
    n1_int = np.array(params[2 * n + 1:3 * n + 1])
    a2_int = np.array(params[3 * n + 1:4 * n + 1])
    n2_int = np.array(params[4 * n + 1:5 * n + 1])
    N_int = np.array(params[5 * n + 1:6 * n + 1])

    Mnoint_fit = params[6 * n + 1]
    sigma_noint = params[6 * n + 2]
    a1_noint = params[6 * n + 3]
    n1_noint = params[6 * n + 4]
    a2_noint = params[6 * n + 5]
    n2_noint = params[6 * n + 6]
    N_noint = params[6 * n + 7]

    global ar_hint, ar_hint_err, ar_hnoint, ar_hnoint_err, x_int, x_noint, bkg_func

    chi2_int_total = 0
    dof_int = 0
    for g, gamma in enumerate(GammaRatio):
        #chi2_int = get_chi2_mass(ar_hint[g], ar_hint_err[g], x_int, bkg_func, a1_int[g], a2_int[g], n1_int[g], n2_int[g], N_int[g], Mint_fit[g], sigma_int[g])
        chi2_int = get_chi2_mass(ar_hint[g], ar_hint_err[g], x_int[g], a1_int[g], a2_int[g], n1_int[g], n2_int[g], N_int[g], Mnoint_fit + alpha*np.sqrt(gamma), sigma_int[g])
        chi2_int_total += chi2_int
        dof_int += (len(x_int[g]) - 6)

    chi2_noint = get_chi2_mass(ar_hnoint, ar_hnoint_err, x_noint, a1_noint, a2_noint, n1_noint, n2_noint, N_noint, Mnoint_fit, sigma_noint)
    chi2_total = chi2_int_total + chi2_noint
    dof_noint = len(x_noint) - 7
    dof_total = dof_int + dof_noint - 1
    #chi2_ndf = chi2_total / dof_total
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
                if f==0:
                    df_qg_vbf = up.open(files_intqg[f])['tagsDumper/trees/ggh_125_13TeV_VBFTag_0'].arrays(['CMS_hgg_mass','weight'], library='pd')
                    df_qg_vbf['weight'] = df_qg_vbf['weight']*(lumi)
                    qg_int_dfs.append(df_qg_vbf)

                df_vbf = up.open(files_intgg[f])['tagsDumper/trees/ggh_125_13TeV_VBFTag_0'].arrays(['CMS_hgg_mass','weight'], library='pd')
                df_vbf['weight'] = df_vbf['weight']*(lumi)
                df_noInt_vbf = up.open(files_sig[f])['tagsDumper/trees/ggh_125_13TeV_VBFTag_0'].arrays(['CMS_hgg_mass','weight'], library='pd')
                df_noInt_vbf['weight'] = df_noInt_vbf['weight']*(lumi)
                gg_int_dfs.append(df_vbf)
                no_int_dfs.append(df_noInt_vbf)
                

            else:
                if f==0:
                    df_qg = up.open(files_intqg[f])['tagsDumper/trees/ggh_125_13TeV_UntaggedTag_%d'%(i)].arrays(['CMS_hgg_mass','weight'], library='pd')
                    df_qg['weight'] = df_qg['weight']*(lumi)
                    qg_int_dfs.append(df_qg)
                    
                df = up.open(files_intgg[f])['tagsDumper/trees/ggh_125_13TeV_UntaggedTag_%d'%(i)].arrays(['CMS_hgg_mass','weight'], library='pd')
                df['weight'] = df['weight']*(lumi)
                gg_int_dfs.append(df)
                df_noInt = up.open(files_sig[f])['tagsDumper/trees/ggh_125_13TeV_UntaggedTag_%d'%(i)].arrays(['CMS_hgg_mass','weight'], library='pd')
                df_noInt['weight'] = df_noInt['weight']*(lumi)
                no_int_dfs.append(df_noInt)
                    
        
        #merge categories if there are multiple files    
        gg_int_events.append(pd.concat(gg_int_dfs, ignore_index=True))
        no_int_events.append(pd.concat(no_int_dfs, ignore_index=True))
        qg_int_events.append(pd.concat(qg_int_dfs, ignore_index=True))
        #remove 0 weights (present in samples Ruben produced)
        no_int_events[i] = no_int_events[i].loc[no_int_events[i]['weight'] != 0]

    #catlabel = [f'UntaggedTag_{i}_{args.year}' for i in range(ncats - 1)] + [f'VBFTag_0_{args.year}']

    for j in range(len(GammaRatio)):
        gg_int_dfs_scaled = [df.copy() for df in gg_int_events]
        qg_int_dfs_scaled = [df.copy() for df in qg_int_events]
        for k in range(len(gg_int_dfs_scaled)): 
            gg_int_dfs_scaled[k]['weight'] = gg_int_dfs_scaled[k]['weight']*np.sqrt(GammaRatio[j])
            qg_int_dfs_scaled[k]['weight'] = qg_int_dfs_scaled[k]['weight']*np.sqrt(GammaRatio[j])
        df_cats = [pd.concat([df_int, df_intqg, df_noint], ignore_index=True) for df_int, df_intqg, df_noint in zip(gg_int_dfs_scaled, qg_int_dfs_scaled, no_int_events)]
        #df_cats = [pd.concat([df_int, df_noint], ignore_index=True) for df_int, df_noint in zip(gg_int_dfs_scaled, no_int_dfs)]
        for c, df_cat in enumerate(df_cats):
            if j==0:
                h_, bin_edges = np.histogram(no_int_events[c]['CMS_hgg_mass'], bins=nbins, range=(hlow,hhigh), weights=no_int_events[c]['weight'])
                herr_, bin_edges = np.histogram(no_int_events[c]['CMS_hgg_mass'], bins=nbins, range=(hlow,hhigh), weights=no_int_events[c]['weight']**2)
                global bin_centers
                bin_centers = 0.5*(bin_edges[1:] + bin_edges[:-1])
                hnoint_array[c] = h_
                hnoint_err_array[c] = np.sqrt(herr_) 
                noint_norm[c] = float(h_.max())

            h_, bin_edges = np.histogram(df_cat['CMS_hgg_mass'], bins=nbins, range=(hlow, hhigh), weights=df_cat['weight'])
            herr_, bin_edges = np.histogram(df_cat['CMS_hgg_mass'], bins=nbins, range=(hlow, hhigh), weights=df_cat['weight']**2)
            hint_array[j,c] = h_
            hint_err_array[j,c] = np.sqrt(herr_) 
            ##fit only the positive part of the spectrum
            last_positive_bin = nbins - 1
            first_positive_bin = 0
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
            filtered_herr_ = (np.sqrt(herr_))[mask] 
            int_norm[c] = float(h_.max())

            filtered_bin_centers_array[j,c] = filtered_bin_centers
            filtered_hint_array[j,c] = filtered_h_
            filtered_hint_err_array[j,c] = filtered_herr_

    ## MINIMIZER
    if (args.categ[0] == 'all'):
        categories = [icat for icat in range(ncats)]
    else:
        categories = [int(icat) for icat in args.categ]

    for ct in categories:
        print(f'CATEGORY {ct}')
        outfile = open(f'category{ct}.txt', 'w')

        # Initial guesses for parameters
        global n
        initial_alpha = -0.02
        #initial_Mint = [124.6] * n
        initial_sigma_int = [s_int_array[ct]] * n
        #initial_sigma_int = [1.2] * n
        initial_a1_int = [1.0] * n
        initial_n1_int = [5.0] * n
        initial_a2_int = [2.0] * n
        initial_n2_int = [10.0] * n
        initial_N_int = [float(int_norm[ct])] * n
        initial_Mnoint = [124.8]
        initial_sigma_noint = [1.2]
        initial_a1_noint = [1.0]
        initial_n1_noint = [5.0]
        initial_a2_noint = [2.0]
        initial_n2_noint = [10.0]
        initial_N_noint = [noint_norm[ct]]

        # Combine all parameters
        initial_params = (
            [initial_alpha]
            #+ initial_Mint
            + initial_sigma_int
            + initial_a1_int
            + initial_n1_int
            + initial_a2_int
            + initial_n2_int
            + initial_N_int
            + initial_Mnoint
            + initial_sigma_noint
            + initial_a1_noint
            + initial_n1_noint
            + initial_a2_noint
            + initial_n2_noint
            + initial_N_noint
        )

        initial_params = np.array(initial_params)

        # Set up parameter names
        param_names = (
            ['alpha']
            #+ [f'Mint[{i}]' for i in range(n)]
            + [f'sigma_int[{i}]' for i in range(n)]
            + [f'a1_int[{i}]' for i in range(n)]
            + [f'n1_int[{i}]' for i in range(n)]
            + [f'a2_int[{i}]' for i in range(n)]
            + [f'n2_int[{i}]' for i in range(n)]
            + [f'N_int[{i}]' for i in range(n)]
            + ['Mnoint']
            + ['sigma_noint']
            + ['a1_noint']
            + ['n1_noint']
            + ['a2_noint']
            + ['n2_noint']
            + ['N_noint']
        )

        lim = []
        ## scipy minimize
        lim.extend([(-0.4,0.04)])
        #lim.extend([(123.8,125.5)] * n)
        lim.extend([(0.5,5)] * n)
        lim.extend([(0.001,5)] * n)
        lim.extend([(0.001,20)] * n)
        lim.extend([(0.001,7)] * n)
        lim.extend([(0.001,20)] * n)
        lim.extend([(int_norm[ct]-2.0, int_norm[ct]+2.0)] * n)
        lim.extend([(123.8,125.5)])
        lim.extend([(0.5,5)])
        lim.extend([(0.001,5)])
        lim.extend([(0.001,20)])
        lim.extend([(0.05,7)])
        lim.extend([(0.001,20)])
        lim.extend([(noint_norm[ct]-2.0, noint_norm[ct]+2.0)])


        global ar_hint, ar_hint_err, ar_hnoint, ar_hnoint_err, x_int, x_noint, bkg_func
        ar_hint = filtered_hint_array[:,ct]
        ar_hint_err = filtered_hint_err_array[:,ct]
        ar_hnoint = hnoint_array[ct]
        ar_hnoint_err = hnoint_err_array[ct]
        x_int=filtered_bin_centers_array[:,ct]
        x_noint = bin_centers
        #bkg_func = bkg_curve_array[ct]

        InitialFit = minimize(lambda parameters: global_chi2(parameters)[0], initial_params, bounds=lim, method='L-BFGS-B')

        # Initialize Minuit
        print("initializing Minuit...")
        m = Minuit(lambda parameters: global_chi2(parameters)[0], InitialFit.x, name=param_names)
        m.errordef = Minuit.LEAST_SQUARES
        m.print_level = 1
        #m.tol = 0.1
        m.limits['alpha'] = (-0.4,0.04)
        m.limits['sigma_noint'] = (0.5,2)
        m.limits['a1_noint'] = (0.001,50)
        m.limits['n1_noint'] = (0.001,300)
        m.limits['a2_noint'] = (0.01,70)
        m.limits['n2_noint'] = (0.001,300)
        m.limits['N_noint'] = (noint_norm[ct]-2.0, noint_norm[ct]+2.0)
        m.limits['Mnoint'] = (123.8,125.5)
        for i in range(n):
            #m.limits[f'Mint[{i}]'] = (123.8,125.5)
            m.limits[f'sigma_int[{i}]'] = (0.5,2)
            m.limits[f'a1_int[{i}]'] = (0.001,50)
            m.limits[f'n1_int[{i}]'] = (0.001,300)
            m.limits[f'a2_int[{i}]'] = (0.01,70)
            m.limits[f'n2_int[{i}]'] = (0.001,300)
            m.limits[f'N_int[{i}]'] = ((int_norm[ct]-2.0, int_norm[ct]+2.0))
            

        # Perform minimization
        print("performing minimization...")
        m.simplex().migrad(ncall=100000)
        #m.migrad(ncall=100000)
        m.hesse()
        print("Done!")

        # Print results
        print("Fitted parameters:")
        for name, value, error in zip(m.parameters, m.values, m.errors):
            print(f"{name} = {value:.4f} ± {error:.4f}")   

        chi2_ap, chi2_am = get_alpha_err(m.values)
        #print(chi2_ap, chi2_am)
        alphaErr = plottool.plot_a_err(m.values[0], alphaVars, chi2_ap, chi2_am, args.outDir, f'alphaErrCat{ct}')

        ## plot
        Mint = []
        Mint_err = []
        dof = global_chi2(m.values)[1]
        for g,gam in enumerate(GammaRatio):
            Mint.append(m.values[f'Mnoint'] + m.values[f'alpha']*np.sqrt(gam))
            Mint_err.append(np.sqrt(m.errors[f'Mnoint']**2 + gam*m.values[f'alpha']**2))
            fit_params = (m.values[f'a1_int[{g}]'], m.values[f'a2_int[{g}]'], m.values[f'n1_int[{g}]'], 
                            m.values[f'n2_int[{g}]'], m.values[f'N_int[{g}]'], Mint[g], m.values[f'sigma_int[{g}]'])
            fit_params_err = (m.errors[f'a1_int[{g}]'], m.errors[f'a2_int[{g}]'], m.errors[f'n1_int[{g}]'], 
                            m.errors[f'n2_int[{g}]'], m.errors[f'N_int[{g}]'], Mint_err[g], m.errors[f'sigma_int[{g}]'])
            label = f'mgg_int_gamma{gam}_cat{ct}'
            plt_hfit = plottool.plot_data(hint_array[:,ct][g], hint_err_array[:,ct][g], bin_centers, x_int[g], fit_params, fit_params_err, m.fval/dof, args.outDir, label)
            #plt_hfit = plottool.plot_data(hint_array[:,ct][g], hint_err_array[:,ct][g], x_noint, x_int[g], fit_params, fit_params_err, FitResult.fun, args.outDir, label)
            outfile.write(f"{gam:.1f} : {m.values[f'a1_int[{g}]']:.3f}, {m.values[f'a2_int[{g}]']:.3f}, {m.values[f'n1_int[{g}]']:.3f}, {m.values[f'n2_int[{g}]']:.3f}, {m.values[f'N_int[{g}]']:.3f}, {Mint[g]:.3f}, {m.values[f'sigma_int[{g}]']:.3f}\n")

        ## plot no int
        fit_params = (m.values[f'a1_noint'], m.values[f'a2_noint'], m.values[f'n1_noint'], 
                            m.values[f'n2_noint'], m.values[f'N_noint'], m.values[f'Mnoint'], m.values[f'sigma_noint'])
        fit_params_err = (m.errors[f'a1_noint'], m.errors[f'a2_noint'], m.errors[f'n1_noint'], 
                            m.errors[f'n2_noint'], m.errors[f'N_noint'], m.errors[f'Mnoint'], m.errors[f'sigma_noint'])
        label = f'mgg_noint_cat{ct}'
        plt_hfit_noint = plottool.plot_data(hnoint_array[ct], hnoint_err_array[ct], bin_centers, x_noint, fit_params, fit_params_err, m.fval/dof, args.outDir, label)
        #plt_hfit_noint = plottool.plot_data(ar_hnoint, ar_hnoint_err, x_noint, x_noint, fit_params, fit_params_err, FitResult.fun, args.outDir, label)
        plt_dm_vs_gr = plottool.plot_dm(m.values, m.errors, GammaRatio, args.outDir, f'dm_vs_gr_cat{ct}')
        #plt_dm_vs_gr = plottool.plot_dm(best_fit_params, param_errors, GammaRatio, args.outDir, f'dm_vs_gr_cat{ct}')
        outfile.write(f"noint : {m.values[f'a1_noint']:.3f}, {m.values[f'a2_noint']:.3f}, {m.values[f'n1_noint']:.3f}, {m.values[f'n2_noint']:.3f}, {m.values[f'N_noint']:.3f}, {m.values[f'Mnoint']:.3f}, {m.values[f'sigma_noint']:.3f}\n")
        outfile.write(f"alpha : {m.values[f'alpha']:.5f}, {alphaErr:.7f}\n")
        outfile.close()        

if __name__ == "__main__":
    main()

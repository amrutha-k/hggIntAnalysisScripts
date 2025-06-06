##### import ROOT
import os
import numpy as np
import math
import uproot as up
import pandas as pd
import matplotlib.pyplot as plt
import mplhep as mpl
from scipy.optimize import curve_fit

from python.model import dscb
import python.plotting_ini as plottool

plt.style.use(mpl.style.ROOT)
plt.style.use(mpl.style.CMS)
plt.rcParams['text.usetex'] = True
plt.rcParams.update({'lines.markeredgewidth': 1})
plt.rcParams['axes.labelweight'] = 'bold'


### open all files as pandas dataframe 

#GammaRatio = [0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0]
GammaRatio = [0.1, 0.3, 0.5, 0.7, 0.9, 1.0, 1.1, 1.3, 1.5, 1.7, 1.9, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]

#infile_gg ='/eos/user/a/amkrishn/hggWidth/mcNtuples/condor_output/2018/intggSummer20UL18oldCats/output_GluGluHToGG_int_M125_13TeV-sherpa.root'
#infile_nonInt = '/eos/user/a/amkrishn/hggWidth/mcNtuples/condor_output/2018/sigSummer20UL18_newFNUF/hadded_trees/output_GluGluHToGG_M125_TuneCP5_13TeV-amcatnloFXFX-pythia8.root'
#infile_gg ='/eos/user/a/amkrishn/hggWidth/mcNtuples/fromRuben/output_GluGluHToGG_int_M125_13TeV-sherpa.root'
#infile_nonInt = '/eos/user/r/rgargiul/dataHggWidth/trees/trees_postVBFcat_sig/output_GluGluHToGG_M125_TuneCP5_13TeV-amcatnloFXFX-pythia8.root'


#UL18_sig = '/eos/user/a/amkrishn/hggWidth/mcNtuples/condor_output/2018/sigNewBoundaries/hadded_tree/output_GluGluHToGG_M125_TuneCP5_13TeV-amcatnloFXFX-pythia8.root'
#UL18_intgg = '/eos/user/a/amkrishn/hggWidth/mcNtuples/condor_output/2018/intNewBoundaries/output_GluGluHToGG_int_M125_13TeV-sherpa.root'

UL17_sig = '/eos/cms/store/group/phys_higgs/cmshgg/rgargiul/trees/trees_sig_UL17/hadded/output_GluGluHToGG_M125_TuneCP5_13TeV-amcatnloFXFX-pythia8_corr.root'
UL17_intgg = '/eos/cms/store/group/phys_higgs/cmshgg/rgargiul/trees/trees_int_UL17/hadded/output_GluGluHToGG_int_M125_13TeV-sherpa_corr.root'
UL17_intqg = '/eos/cms/store/group/phys_higgs/cmshgg/rgargiul/trees/trees_intqg_UL17/hadded/output_GluGluHToGG_intqg_M125_13TeV-sherpa_corr.root'
#files_intgg = [UL16pre_intgg, UL16post_intgg, UL17_intgg, UL18_intgg]
#files_sig = [UL16pre_sig, UL16post_sig, UL17_sig, UL18_sig]

files_intgg = [UL17_intgg]
files_intqg = [UL17_intqg]
files_sig = [UL17_sig]

'''
for i in range(10): 
    gg_int_dfs = []
    no_int_dfs = []

    for f in range(len(files_intgg)):
        df = up.open(files_intgg[f])['tagsDumper/trees/ggh_125_13TeV_UntaggedTag_%d'%(i)].arrays(['CMS_hgg_mass','weight'], library='pd')
        #df['weight'] = df['weight']*(138.0)
        gg_int_dfs.append(df)
        df_noInt = up.open(files_sig[f])['tagsDumper/trees/ggh_125_13TeV_UntaggedTag_%d'%(i)].arrays(['CMS_hgg_mass','weight'], library='pd')
        #df_noInt['weight'] = df_noInt['weight']*(138.0)
        no_int_dfs.append(df_noInt)
        
    gg_int_events.append(pd.concat(gg_int_dfs))
    no_int_events.append(pd.concat(no_int_dfs))

'''
gg_int_events = []
qg_int_events = []
no_int_events = []

for i in range(11):
    qg_int_dfs = []
    gg_int_dfs = []
    no_int_dfs = []

    for f in range(len(files_intgg)):
        if i == 10:
            df_vbf = up.open(files_intgg[f])['tagsDumper/trees/ggh_125_13TeV_VBFTag_0'].arrays(['CMS_hgg_mass','weight'], library='pd')
            df_vbf['weight'] = df_vbf['weight']*(41.48)
            df_intqg_vbf = up.open(files_intqg[f])['tagsDumper/trees/ggh_125_13TeV_VBFTag_0'].arrays(['CMS_hgg_mass','weight'], library='pd')
            df_intqg_vbf['weight'] = df_intqg_vbf['weight']*(41.48)
            df_noInt_vbf = up.open(files_sig[f])['tagsDumper/trees/ggh_125_13TeV_VBFTag_0'].arrays(['CMS_hgg_mass','weight'], library='pd')
            df_noInt_vbf['weight'] = df_noInt_vbf['weight']*(41.48)
            gg_int_dfs.append(df_vbf)
            no_int_dfs.append(df_noInt_vbf)
            qg_int_dfs.append(df_intqg_vbf)
        else:
            df = up.open(files_intgg[f])['tagsDumper/trees/ggh_125_13TeV_UntaggedTag_%d'%(i)].arrays(['CMS_hgg_mass','weight'], library='pd')
            df['weight'] = df['weight']*(41.48)
            gg_int_dfs.append(df)
            df_intqg = up.open(files_intqg[f])['tagsDumper/trees/ggh_125_13TeV_UntaggedTag_%d'%(i)].arrays(['CMS_hgg_mass','weight'], library='pd')
            df_intqg['weight'] = df_intqg['weight']*(41.48)
            qg_int_dfs.append(df_intqg)
            df_noInt = up.open(files_sig[f])['tagsDumper/trees/ggh_125_13TeV_UntaggedTag_%d'%(i)].arrays(['CMS_hgg_mass','weight'], library='pd')
            df_noInt['weight'] = df_noInt['weight']*(41.48)
            no_int_dfs.append(df_noInt)
        
            
    gg_int_events.append(pd.concat(gg_int_dfs))
    no_int_events.append(pd.concat(no_int_dfs))
    qg_int_events.append(pd.concat(qg_int_dfs))
    
gg_int_merged = pd.concat(gg_int_events)
qg_int_merged = pd.concat(qg_int_events)
no_int_merged = pd.concat(no_int_events)


plotDir = '/eos/user/a/amkrishn/www/hggWidth/AdditionalPlotsReview/'
#no_int_events['weight'] = no_int_events['weight']*(138.0)
#gg_int_events['weight'] = gg_int_events['weight']*(138.0) # modify gg int weights

#qg_int_events = up.open('/eos/user/a/amkrishn/hggWidth/mcNtuples/condor_output/UntaggedTag_dump_UL17_sigMC_intqg_correctNormalisation_v2/output_GluGluHToGG_intqg_M125_13TeV-sherpa_full.root')\
#['tagsDumper/trees/GluGluHToGG_intqg_M125_13TeV_sherpa_13TeV_UntaggedTag'].arrays(['CMS_hgg_mass','weight','diphoton_pt','diphoton_mva'], library='pd')

    
def make_histo_fit(df_, hname):
    
    nbins = 20
    
    h_, bin_edges = np.histogram(df_['CMS_hgg_mass'], bins=nbins, range=(118,130), weights=df_['weight'])
    h_err, bin_edges = np.histogram(df_['CMS_hgg_mass'], bins=nbins, range=(118,130), weights=df_['weight']**2)

    bin_centers = 0.5*(bin_edges[1:] + bin_edges[:-1])
    last_positive_bin = nbins-1
    first_positive_bin = 0
    #fit only the positive part of the spectrum
    midpoint_h_ = int(0.5*len(h_))
    for ibin in range(midpoint_h_):
        if h_[ibin]<0:
            first_positive_bin = ibin + 1
        if h_[ibin + midpoint_h_]<0:
            #last_positive_bin = (ibin + midpoint_h_) - 1
            last_positive_bin = ibin + midpoint_h_ 
            break
        
    chi2_ndf = 10.0
    while (chi2_ndf > 6.0 or np.isnan(chi2_ndf)) and last_positive_bin > first_positive_bin:
        last_positive_bin -= 1
        mask = (bin_centers >= bin_centers[first_positive_bin]) & (bin_centers <= bin_centers[last_positive_bin])
        filtered_bin_centers = bin_centers[mask]
        filtered_h_ = h_[mask]
        filtered_h_err = h_err[mask]   
        
        initial_mean = np.average(filtered_bin_centers, weights=filtered_h_)
        initial_var = np.average((filtered_bin_centers - initial_mean)**2, weights=filtered_h_)
        

        try:
            popt, pcov = curve_fit(
                dscb, 
                filtered_bin_centers, 
                filtered_h_, 
                p0=([1, 1, 1, 1, filtered_h_.max(), initial_mean, np.sqrt(initial_var)]), 
                sigma=filtered_h_err**0.5, 
                method='trf', 
                maxfev=30000
            )
        except RuntimeError as e:
            print(f"Fit failed: {e}")
            popt, pcov = None, None  # Assign default values or handle as needed

        # Continue with your script
        if popt is not None:
            ##chi2 per dof
            y_fit_skim = dscb(filtered_bin_centers,*popt)
            res = (filtered_h_ - y_fit_skim) / (filtered_h_err**0.5)
            chi2 = np.sum(res**2)
            chi2_ndf = chi2 / (len(filtered_h_) - len(popt))

        #popt, pcov = curve_fit(dscb, filtered_bin_centers, filtered_h_, p0=([1,1,1,1,filtered_h_.max(),
        #                    initial_mean,np.sqrt(initial_var)]),sigma=filtered_h_err**0.5, method='trf', maxfev=40000)
        
        
    


    #plot
    fig, ax = plt.subplots(figsize=(10,8))
    fig.set_dpi(75)
    plt.errorbar(bin_centers, h_, xerr=0.5*(bin_centers[1]-bin_centers[0]), yerr=h_err**0.5, fmt='o', markersize=6, c='black', elinewidth=0.7, capsize=1)
    plottool.common_plot_options_mass()
    x_fit = np.linspace(bin_centers[first_positive_bin],bin_centers[last_positive_bin],100)
    y_fit = dscb(x_fit,*popt)
    plt.plot(x_fit, y_fit, color='red')
    perr = [np.sqrt(pcov[j,j]) for j, _ in enumerate(popt)]
    plottool.add_text_box(popt, perr)
    textbox = r'$\chi^2/\mathrm{ndf}=%.2f$' % chi2_ndf

    plt.text(0.75, 0.95, textbox, transform=plt.gca().transAxes, fontsize=19, fontweight='bold', verticalalignment='top', bbox=dict(facecolor='white', edgecolor='none'))
    
    plt.savefig('%s/%s.png'%(plotDir,hname))
    plt.savefig('%s/%s.pdf'%(plotDir,hname))
    plt.clf()
    
    #save data
    np.savez_compressed('%s/%s_data.npz'%(plotDir, hname),
                    fit_params=popt,
                    fit_param_err=perr
                    )
        
    


# gg interference
h_int, bin_edges = np.histogram(gg_int_merged['CMS_hgg_mass'], bins=60, range=(118,130), weights=gg_int_merged['weight'])
weights_sq = (gg_int_merged['weight'])**2
h_int_err, bin_edges = np.histogram(gg_int_merged['CMS_hgg_mass'], bins=60, range=(118,130), weights=weights_sq)
bin_centers = 0.5*(bin_edges[1:] + bin_edges[:-1])
fig, ax = plt.subplots(figsize=(10,8))
fig.set_dpi(75)
plt.errorbar(bin_centers, h_int, yerr=h_int_err**0.5, fmt='o', markersize=6, c='black')
plottool.common_plot_options_int()
plt.rcParams['text.usetex'] = False
plt.savefig('%s/gg_interference.png'%plotDir)
plt.savefig('%s/gg_interference.pdf'%plotDir)
plt.clf()


#qg interference
h_int_qg, bin_edges = np.histogram(qg_int_merged['CMS_hgg_mass'], bins=60, range=(118,130), weights=qg_int_merged['weight'])
h_int_qg_err, bin_edges = np.histogram(qg_int_merged['CMS_hgg_mass'], bins=60, range=(118,130), weights=qg_int_merged['weight']**2)
bin_centers = 0.5*(bin_edges[1:] + bin_edges[:-1])
plt.errorbar(bin_centers, h_int_qg, yerr=h_int_qg_err**0.5, fmt='o', markersize=6, c='black')
plottool.common_plot_options_int()
plt.savefig('%s/qg_interference.png'%plotDir)
plt.savefig('%s/qg_interference.pdf'%plotDir)
plt.clf()


# Mgg with interference (inclusive pT)                                                                                                                                                                    
#df_with_int = pd.concat([nonInt_events, gg_int_events, qg_int_events], ignore_index=True)  ##df with interference added                                                                                  
#df_with_int = pd.concat([no_int_merged, gg_int_merged, qg_int_merged], ignore_index=True)  ##df with interference added                                                                                  

# no interference
h_noint, bin_edges = np.histogram(no_int_merged['CMS_hgg_mass'], bins=60, range=(118,130), weights=no_int_merged['weight'])
h_noint_err, bin_edges = np.histogram(no_int_merged['CMS_hgg_mass'], bins=60, range=(118,130), weights=no_int_merged['weight']**2)
bin_centers = 0.5*(bin_edges[1:] + bin_edges[:-1])
plt.errorbar(bin_centers, h_noint, yerr=h_noint_err**0.5, fmt='o', markersize=6, c='black')
plottool.common_plot_options_int()
plt.savefig('%s/no_interference.png'%plotDir)
plt.savefig('%s/no_interference.pdf'%plotDir)
plt.clf()

h_ratio_gg=h_int/h_noint
h_ratio_gg_err=np.abs(h_ratio_gg)*np.sqrt((h_int_err**0.5/h_int)**2 + (h_noint_err**0.5/h_noint)**2)

h_ratio_qg=h_int_qg/h_noint
h_ratio_qg_err=np.abs(h_ratio_qg)*np.sqrt((h_int_qg_err**0.5/h_int_qg)**2 + (h_noint_err**0.5/h_noint)**2)
plt.errorbar(bin_centers, h_ratio_gg, yerr=h_ratio_gg_err, fmt='o', markersize=6, c='red', label='gg_int')
plt.errorbar(bin_centers, h_ratio_qg, yerr=h_ratio_qg_err, fmt='o', markersize=6, c='blue', label='qg_int')
plt.xlabel(r'$M_{\gamma\gamma}\ (GeV)$')
plt.ylabel('ratio')
plt.xlim(118,130)
plt.xticks(np.arange(118,130,2))
plt.legend()
plt.savefig('%s/ratio_gg_qg.png'%plotDir)
plt.savefig('%s/ratio_gg_qg.pdf'%plotDir)
plt.clf()

plt.errorbar(bin_centers, h_ratio_qg, yerr=h_ratio_qg_err, fmt='o', markersize=6, c='blue')
#plt.xlabel(r'$M_{\gamma\gamma}\ (GeV)$')
#plt.ylabel('ratio')
#plt.xlim(118,130)
#plt.xticks(np.arange(118,130,2))                                                                                                           
#plottool.common_plot_options_int()
#plt.savefig('%s/ratio_qg.png'%plotDir)
#plt.savefig('%s/ratio_qg.pdf'%plotDir)
#plt.clf()


print("***********DONE**************")

#make_histo_fit(no_int_merged, 'mgg_noint')

# Mgg with interference (inclusive pT)
#df_with_int = pd.concat([nonInt_events, gg_int_events, qg_int_events], ignore_index=True)  ##df with interference added
#df_with_int = pd.concat([no_int_merged, gg_int_merged], ignore_index=True)  ##df with interference added
#make_histo_fit(df_with_int, 'mgg_int')  #duplicate in the loop below


# Mgg inclusive
#for gR in range(len(GammaRatio)):
#    gg_int_merged_copy = gg_int_merged.copy()
#    gg_int_merged_copy['weight'] = gg_int_merged_copy['weight']*np.sqrt(GammaRatio[gR])
#    df_int_g = pd.concat([no_int_merged, gg_int_merged_copy], ignore_index=True)  ##df with interference added
#    make_histo_fit(df_int_g, 'mgg_int_gamma%0.1f'%(GammaRatio[gR]))
    
"""    
# Mgg in pT bins
#remove the low MVA cats from ptbin1 and 2
#gg_int_dfs_trimmed = gg_int_dfs.copy()
#del gg_int_dfs_trimmed[5]
#del gg_int_dfs_trimmed[3]
#paired_int_dfs = zip(gg_int_dfs_trimmed[0::2],gg_int_dfs_trimmed[1::2])
paired_int_dfs = zip(gg_int_dfs[0::2],gg_int_dfs[1::2])
pt_bins_int = [pd.concat([df1,df2], ignore_index=True) for df1, df2 in paired_int_dfs]

#remove the low MVA cats from ptbin1 and 2
#no_int_dfs_trimmed = no_int_dfs.copy()
#del no_int_dfs_trimmed[5]
#del no_int_dfs_trimmed[3]
#paired_noint_dfs = zip(no_int_dfs_trimmed[0::2],no_int_dfs_trimmed[1::2])
paired_noint_dfs = zip(no_int_dfs[0::2],no_int_dfs[1::2])
pt_bins_noint = [pd.concat([df1,df2], ignore_index=True) for df1, df2 in paired_noint_dfs]

for j in range(len(GammaRatio)):
    pt_bins_int_copy = [df.copy() for df in pt_bins_int]
    for k in range(len(pt_bins_int_copy)):
        pt_bins_int_copy[k]['weight'] = pt_bins_int_copy[k]['weight']*np.sqrt(GammaRatio[j]) 
    df_pts = [pd.concat([df_int, df_noint], ignore_index=True) for df_int, df_noint in zip(pt_bins_int_copy,pt_bins_noint)]
    #df_pt12_merged = [df_pts[0], pd.concat([df_pts[1], df_pts[2]]), df_pts[3], df_pts[4]]  #merge pT bin 1 and 2
    for nbin, df_pt in enumerate(df_pts):
        make_histo_fit(df_pt, 'mgg_int_gamma%0.1f_ptbin%d'%(GammaRatio[j],nbin))
        if j==0:
            make_histo_fit(pt_bins_noint[nbin], 'mgg_noint_ptbin%d'%nbin)
"""

for j in range(len(GammaRatio)):
    gg_int_dfs_scaled = [df.copy() for df in gg_int_events]
    for k in range(len(gg_int_dfs_scaled)):
        gg_int_dfs_scaled[k]['weight'] = gg_int_dfs_scaled[k]['weight']*np.sqrt(GammaRatio[j])
        #if j==0:
            #no_int_events[k]['weight'] = no_int_events[k]['weight']*138
    df_cats = [pd.concat([df_int, df_noint], ignore_index=True) for df_int, df_noint in zip(gg_int_dfs_scaled,no_int_events)]
    #df_pt12_merged = [df_pts[0], pd.concat([df_pts[1], df_pts[2]]), df_pts[3], df_pts[4]]  #merge pT bin 1 and 2
    for nbin, df_cat in enumerate(df_cats):
        make_histo_fit(df_cat, 'mgg_int_gamma%0.1f_cat%d'%(GammaRatio[j],nbin))
        if j==0:
            #print(no_int_dfs[nbin])
            make_histo_fit(no_int_events[nbin], 'mgg_noint_cat%d'%nbin)


# dM vs Gamma ratio per pT bin
#plotDir = '/eos/user/a/amkrishn/www/hggWidth/dMvsWidthXsecScaling/newFNUF_VBFTag0_UL17'

##model definition
def model_func(x,a):
    return a*np.sqrt(x)

## plot settings
fig, ax = plt.subplots(figsize=(9,6))
fig.set_dpi(100)
#cats = [1]
alpha_list = []
alpha_err_list =[]

for cat in range(6):
    M_int = [np.load(f'{plotDir}/mgg_int_gamma{GammaRatio[gr]:0.1f}_cat{cat}_data.npz')['fit_params'][5] for gr in range(len(GammaRatio))]
    M_noint = np.load(f'{plotDir}/mgg_noint_cat{cat}_data.npz')['fit_params'][5]
    #dM = [((np.load(f'{plotDir}/mgg_int_gamma{GammaRatio[gr]:0.1f}_cat{cat}_data.npz')['fit_params'][5] - np.load(f'{plotDir}/mgg_noint_cat{cat}_data.npz')['fit_params'][5])) for gr in range(len(GammaRatio))]
    dM = [M_int[i] - M_noint for i in range(len(M_int))]
    #print(dM)
    noint_err_sq = pow(np.load(f'{plotDir}/mgg_noint_cat{cat}_data.npz')['fit_param_err'][5], 2)
    int_err_sq = [pow(np.load(f'{plotDir}/mgg_int_gamma{GammaRatio[gr]:0.1f}_cat{cat}_data.npz')['fit_param_err'][5], 2) for gr in range(len(GammaRatio))]
    #print(noint_err_sq)
    #print(int_err_sq)
    dMerr = [np.sqrt(int_err_sq[gr] + noint_err_sq) for gr in range(len(GammaRatio))]
    #print(dMerr)

    x = []
    y = []
    y_err = []
    count = 0
    
    #if (cat != 10 and cat != 0):
    for ie, ele in enumerate(M_int):
        if (ie == len(M_int)-1):
            break
        if (M_int[ie+1] - M_int[ie]) >= 0.00:
            count += 1
        if (count >= 2): break
        if (dMerr[ie] > 0.2): break
            
        if M_int[ie] == None:
            print(M_int[ie])
            continue
                
        x.append(GammaRatio[ie])
        y.append(dM[ie]*1000)
        y_err.append(dMerr[ie]*1000)
    #else:
    #    for ie, ele in enumerate(M_int):
    #        x.append(GammaRatio[ie])
    #        y.append(dM[ie]*1000)
    #        y_err.append(dMerr[ie]*1000)
        
        
    #print("x: ", x)
    #print("y: ", y)
    #print("y_err: ", y_err)
    
            
    plt.errorbar(x, y, yerr=y_err, fmt='o', markersize=4, c='black', elinewidth=0.7, capsize=1)
    #plt.xscale('log')
    plt.xlabel(r'${\Gamma_H}/{\Gamma_H^{SM}}$', fontsize=22)
    plt.ylabel(r'$\Delta M$ (MeV)', fontsize=22)
    plt.xticks(fontsize=16, fontweight='bold')
    plt.yticks(fontsize=16, fontweight='bold')

    popt, pcov = curve_fit(model_func, x, y, sigma=y_err, p0=[-0.02], maxfev = 5000, absolute_sigma=True)
    #print(popt)
    #print(pcov)
    perr = np.sqrt(np.diag(pcov)) #perr
    #print(perr)
    print(f'cat{cat}: alpha{popt[0]}, error{perr}')
    alpha_list.append(popt[0])
    alpha_err_list.append(perr)
    
    x_fit = np.linspace(0,max(GammaRatio)+2.0,1000)
    y_fit = model_func(x_fit, *popt)
    y_fit_u = model_func(x_fit, *(popt+perr))
    y_fit_l = model_func(x_fit, *(popt-perr))

    ##plot fit function with error envelopes
    plt.plot(x_fit,y_fit,'r')
    plt.fill_between(x_fit, y_fit_u, y_fit_l, color='orange', alpha=0.5)
    plt.tick_params(axis='both',which='minor',labelsize=9)

    ##chi2 per dof
    res = (y - model_func(x, *popt)) / y_err
    #print(res)
    chi2 = np.sum(res**2)
    #print(chi2)
    chi2_ndf = chi2 / (len(dM) - len(popt))
    #print(chi2_ndf)

    textbox = r'$\chi^2/\mathrm{ndf}=%.2f$' % chi2_ndf
    textbox += '\n'
    textbox += r'$\alpha=%.2f \pm %.2f$' % (popt[0], perr[0])

    plt.text(0.75, 0.9, textbox, transform=plt.gca().transAxes, fontsize=14, verticalalignment='top', bbox=dict(facecolor='white', alpha=0.2))
    plt.savefig('%s/dM_vs_gR_cat%d.png'%(plotDir,cat))
    plt.savefig('%s/dM_vs_gR_cat%d.pdf'%(plotDir,cat))
    plt.clf()
    
    ### Second Plot (Logarithmic x-axis)
    fig, ax = plt.subplots(figsize=(10, 6))
    fig.set_dpi(100)
    
    ax.set_xscale('log')  # Set x-axis to log scale
    plt.errorbar(x, y, yerr=y_err, fmt='o', markersize=4, c='black', elinewidth=0.7, capsize=1)
    plt.plot(x_fit, y_fit, 'r')
    plt.fill_between(x_fit, y_fit_u, y_fit_l, color='orange', alpha=0.5)
    plt.xlabel(r'$\Gamma_H / \Gamma_H^{SM}$', fontsize=22)
    plt.ylabel(r'$\Delta M$ (MeV)', fontsize=22)
    plt.xticks(fontsize=16, fontweight='bold')
    plt.yticks(fontsize=16, fontweight='bold')
    
    # Add chi-squared text box
    plt.text(0.75, 0.9, textbox, transform=plt.gca().transAxes, fontsize=14, verticalalignment='top', bbox=dict(facecolor='white', alpha=0.4))

    # Save log-scale plot
    plt.savefig(f'{plotDir}/dM_vs_gR_cat{cat}_logx.png')
    plt.savefig(f'{plotDir}/dM_vs_gR_cat{cat}_logx.pdf')
    plt.clf()

print(alpha_list)
alpha_err_list_flat = np.concatenate(alpha_err_list).ravel().tolist()
print(alpha_err_list_flat)

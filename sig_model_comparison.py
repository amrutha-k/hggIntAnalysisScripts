#!/usr/bin/env bash
#!/usr/bin/env python3

import ROOT
from ROOT import TMath
from ROOT import TLatex
import numpy as np
import math
import uproot as up
import pandas as pd
import python.models as model

### open all files as pandas dataframe (all files are scaled to 137/fb)

infile_gg ='/eos/user/a/amkrishn/hggWidth/mcNtuples/condor_output/UntaggedTag_dump_hggInt_UL17_sigMC_genweight0/output_GluGluHToGG_int_M125_13TeV-sherpa_run2UL17_0.root'
gg_int_files = []
for i in range(4):
    df = up.open(infile_gg)['tagsDumper/trees/GluGluHToGG_int_M125_13TeV_sherpa_13TeV_UntaggedTag_%d'%(i)].arrays(['CMS_hgg_mass','weight','diphoton_pt','diphoton_mva'], library='pd')
    gg_int_files.append(df)
    #print(df)
gg_int_events = pd.concat(gg_int_files)
gg_int_events['weight'] = gg_int_events['weight']*(1.583e-03/(43.92*0.00227*1.06)) # modify gg int weights

qg_int_events = up.open('/eos/user/a/amkrishn/hggWidth/mcNtuples/condor_output/UntaggedTag_dump_UL17_sigMC_intqg_correctNormalisation_v2/output_GluGluHToGG_intqg_M125_13TeV-sherpa_full.root')\
['tagsDumper/trees/GluGluHToGG_intqg_M125_13TeV_sherpa_13TeV_UntaggedTag'].arrays(['CMS_hgg_mass','weight','diphoton_pt','diphoton_mva'], library='pd')

nonInt_events = up.open('/eos/user/a/amkrishn/hggWidth/mcNtuples/condor_output/UntaggedTag_dump_UL17_sigMC_nonInt_genweight0/output_GluGluHToGG_M125_TuneCP5_13TeV-amcatnloFXFX-pythia8_storeWeights.root')['tagsDumper/trees/GluGluHToGG_M125_TuneCP5_13TeV_amcatnloFXFX_pythia8_13TeV_UntaggedTag'].arrays(['CMS_hgg_mass','weight','diphoton_pt','diphoton_mva'], library='pd')

def make_histo(df_, hname):
    hist = ROOT.TH1F(hname,";M^{#gamma#gamma} (GeV); entries",28,116,130)
    # pt and mva bins
    #pt_bins = [15, 25, 40, 64, 114, 3000]
    #mva_bins = [0.72, 0.86, 0.74, 0.86, 0.76, 0.86, 0.77, 0.88, 0.87, 0.95]

    # mappings between categories and their corresponding pt and mva masks
    category = {
    'cat0': {'pt': (114, 3000), 'mva': (0.95, 1)},
    'cat1': {'pt': (114, 3000), 'mva': (0.87, 0.95)},
    'cat2': {'pt': (64, 114), 'mva': (0.88, 1)},
    'cat3': {'pt': (64, 114), 'mva': (0.77, 0.88)},
    'cat4': {'pt': (40, 64), 'mva': (0.86, 1)},
    'cat5': {'pt': (40, 64), 'mva': (0.76, 0.86)},
    'cat6': {'pt': (25, 40), 'mva': (0.86, 1)},
    'cat7': {'pt': (25, 40), 'mva': (0.74, 0.86)},
    'cat8': {'pt': (15, 25), 'mva': (0.86, 1)},
    'cat9': {'pt': (15, 25), 'mva': (0.72, 0.86)},
    'inclusive': {'pt': (15, 3000), 'mva': (0.72, 1)},    
    }

    # Iterate over the category_masks dictionary and apply conditions for each category
    for cat, ranges in category.items():
        if cat in hname:
            pt_range = ranges['pt']
            mva_range = ranges['mva']
            pt_mask = np.logical_and(df_['diphoton_pt'] >= pt_range[0], df_['diphoton_pt'] < pt_range[1])
            mva_mask = np.logical_and(df_['diphoton_mva'] > mva_range[0], df_['diphoton_mva'] <= mva_range[1])
            skimmed_df = df_[pt_mask & mva_mask]
            for index, row in skimmed_df.iterrows():
                ev = row['CMS_hgg_mass']
                wei = row['weight']
                hist.Fill(ev,wei)
            break
    ## initial fit to estimate parameter values
    #func = ROOT.TF1("func",dscb,118,130,7)
    #func.SetParameters(1,1,1,1,hist.GetMaximum(), hist.GetMean(), hist.GetRMS())
    #hist.Fit(func,"QR")
    
    #c = ROOT.TCanvas()
    #hist.Draw()
    #c.SaveAs('%s/%s.png'%(plotDir,hname))
    #c.SaveAs('%s/%s.pdf'%(plotDir,hname))
    return hist

def fitnGaus(h, hname, n):
    ### fit ###
    # convert TH1F to RooDataHist
    print("test1")
    m = ROOT.RooRealVar("m", "M_{#gamma#gamma} (GeV)", 116.0, 130.0)
    data = ROOT.RooDataHist("h","h",ROOT.RooArgSet(m), h)
    print("test2")
    # Define the model parameters
    # Get initial values from the curve_fit results
    #mean_initial = np.load(f'{hname}_data.npz')['fit_params'][5]
    #sigma_initial = np.load(f'{hname}_data.npz')['fit_params'][6]
    mean_initial = 124.8
    sigma_initial = 1.19
    model = []
    frac = []
    print("test3")
    for i in range(n):
        mean = ROOT.RooRealVar("mean", "mean", mean_initial, mean_initial - 5.0, mean_initial + 5.0)
        sigma = ROOT.RooRealVar("sigma", "sigma", sigma_initial, 0.01, sigma_initial + 2)
        model.append(ROOT.RooGaussian("model_%d"%i, "Gaussian", m, mean, sigma))
        
    
    # fit
    _pdfs, _coeffs = ROOT.RooArgList(), ROOT.RooArgList()
    for pdf in model: _pdfs.add(pdf)
    for g in range(n-1):
        frac.append(ROOT.RooRealVar("frac","frac",0.5-0.05*g,0.01,0.99))
    #frac_constrained = ROOT.RooFormulaVar("frac_constrained","frac_constrained","(@0>0)*(@0<1)*@0+ (@0>1.0)*0.9999",ROOT.RooArgList(frac))
    for coeff in frac: _coeffs.add(coeff)
    #print(_coeffs)
    fmodel = ROOT.RooAddPdf("fmodel", "nGaussians", _pdfs, _coeffs, recursiveFraction=True)
    #print(going to fit...)
    fit = fmodel.fitTo(data, ROOT.RooFit.Save(), ROOT.RooFit.SumW2Error(False))
    #print(fitting...)
    #fit.Print("V")
    # plot
    #c.Clear()
    c = ROOT.TCanvas()
    frame = m.frame()
    data.plotOn(frame)
    fmodel.plotOn(frame)
    fmodel.paramOn(frame)
    frame.Draw()
    c.SaveAs('%s/%s_ngausfit.png'%(plotDir,hname))
    c.SaveAs('%s/%s_ngausfit.pdf'%(plotDir,hname))

    ## main
plotDir = '/eos/user/a/amkrishn/hggWidth/mgg_sig_model_comparison/nGaus'
# no interference
#hnoint = make_histo(nonInt_events, 'mgg_noint_inclusive')
#fitCB(hnoint, 'mgg_noint_inclusive')

# Mgg with interference (inclusive pT)
df_with_int = pd.concat([nonInt_events, gg_int_events, qg_int_events], ignore_index=True)  ##df with interference added
#hint = make_histo(df_with_int, 'mgg_int_inclusive')
#fitCB(hint, 'mgg_int_inclusive')

# Make histograms per category and fit
hcat0 = make_histo(df_with_int, 'mgg_int_cat0')
fitnGaus(hcat0, 'mgg_int_cat0', 2)
    


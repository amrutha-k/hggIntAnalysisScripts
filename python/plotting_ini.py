import matplotlib.pyplot as plt
import numpy as np
from python.model import dscb,splusb,parabola
import mplhep as mpl
mpl.style.use("ROOT")
#mpl.style.use("CMS")
from scipy.optimize import curve_fit

#plot setup
plt.rcParams.update(
    {
        "text.usetex": True,
        "text.latex.preamble": r"\usepackage{bm}",

        # Enforce default LaTeX font.
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica"],
    }
)

def common_plot_options_int():
    plt.xlabel(r'$M_{\gamma\gamma}\ (GeV)$')
    plt.ylabel('entries')
    #mpl.cms.text('Internal')
    plt.xlim(115,135)
    plt.xticks(np.arange(115,135,2))
    
def common_plot_options_mass(hl, hh):
    plt.xlabel(r'$\bm{M_{\gamma\gamma}}$ (GeV)', fontsize=22, weight='bold')
    plt.ylabel('entries', fontsize=22, weight='bold')
    #mpl.cms.text('Internal')
    plt.xlim(hl,hh)
    plt.xticks(np.arange(hl,hh,2), fontsize=18, fontweight='bold')
    plt.yticks(fontsize=18, fontweight='bold')
    
    
def add_text_box(params, param_err):
    param_names = ['aL', 'aR', 'nL', 'nR', 'norm', 'mean', 'sigma']
    param_dict = {name: (value, error) for name, value, error in zip(param_names, params, param_err)}
    params_str = [f'{name} = {val:.2f} ± {err:.2f}' for name, (val, err) in param_dict.items()]
    textstr = '\n'.join(params_str)
    props = dict(boxstyle='square', facecolor='white', alpha=0.8)
    plt.text(0.05, 0.95, textstr, transform=plt.gca().transAxes, fontsize=14, verticalalignment='top', bbox=props)

#def plot_data(hy, hyerr, hx, fhx, fpars, fparerrs, bfunc_, bscale_, chi2ndf, plotdir, hname):
def plot_data(hy, hyerr, hx, fhx, fpars, fparerrs, chi2ndf, plotdir, hname):
    #fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10,8), gridspec_kw={'height_ratios': [3, 1]}, sharex=True)
    fig, ax = plt.subplots(figsize=(10,8))
    fig.set_dpi(75)
    plt.errorbar(hx, hy, xerr=0.5*(hx[1]-hx[0]), yerr=hyerr**0.5, fmt='o', markersize=6, c='black', elinewidth=0.7, capsize=1)
    common_plot_options_mass(hx[0]-0.05,hx[-1]-0.05)
    xfit = np.linspace(fhx[0],fhx[-1],100)
    #yfit = splusb(xfit, bfunc_, bscale_, *fpars)
    yfit = dscb(xfit, *fpars)
    plt.plot(xfit, yfit, color='blue')
    add_text_box(fpars, fparerrs)
    textbox = r'$\bm{\chi^2/\mathrm{ndf}=%.2f}$' % chi2ndf
    plt.text(0.74, 0.95, textbox, transform=plt.gca().transAxes, fontsize=19, fontweight='bold', verticalalignment='top', bbox=dict(facecolor='white', edgecolor='none'))

    # --- Residuals Plot ---
    #residuals = (hy - dscb(hx, *fpars)) / np.sqrt(hyerr)  # Residuals = (data - fit) / error
    #ax2.errorbar(hx, residuals, xerr=0.5*(hx[1]-hx[0]), yerr=1, fmt='o', markersize=6, c='black', elinewidth=0.7, capsize=1)
    #ax2.axhline(0, color='red', linestyle='--', linewidth=1)  # Zero residual line

    # Formatting residuals plot
    #ax2.set_ylabel('Residuals', fontsize=22)
    #ax2.set_xlabel('Mass')
    #ax2.set_ylim(-5, 5)  # Adjust depending on your residuals' range
    #ax2.grid(True, linestyle='dotted', alpha=0.7)
 
    plt.savefig('%s/%s.png'%(plotdir,hname))
    plt.savefig('%s/%s.pdf'%(plotdir,hname))
    plt.close()

def plot_dm(pvalues, pverrs, gr, plotdir, pname):
    fig, ax = plt.subplots(figsize=(10,8))
    fig.set_dpi(75)
    
    a_ = pvalues[0]
    Mnoint_ = pvalues[6 * len(gr) + 1]
    #Mnointerr_ = pverrs[6 * len(gr) + 1]
    Mint_ = Mnoint_ + a_*np.sqrt(gr)
    #Minterr_ = pverrs[1:len(gr) + 1]

    y_val = [(Mint_[j] - Mnoint_) for j in range(len(gr))]
    #print(y_val)
    #y_val_err = [np.sqrt(Minterr_[k]**2 + Mnointerr_**2) for k in range(len(gr))]

    plt.errorbar(gr, y_val, fmt='o', markersize=4, c='black', elinewidth=0.7, capsize=1)
    x_fit = np.linspace(0,max(gr)+1.0,100)
    y_fit = pvalues[0]*np.sqrt(x_fit)
    plt.plot(x_fit,y_fit,'r')
    plt.tick_params(axis='both',which='minor',labelsize=9)
    plt.xlabel(r'$\bm{\Gamma_H / \Gamma^{H}_{SM}}$')
    plt.ylabel(r'$\bm{\Delta\ M\ (GeV)}$')

    plt.savefig('%s/%s.png'%(plotdir,pname))
    plt.savefig('%s/%s.pdf'%(plotdir,pname))
    plt.close()


def plot_a_err(b, aVars, chi2aplus, chi2aminus, plotdir, plabel):
    x1 = [b+i for i in aVars] 
    plt.plot(x1, chi2aplus, '*', color='black')
    x2 = [b-i for i in aVars]
    plt.plot(x2, chi2aminus, '*', color='black')

    if b < 0:
        x1.reverse()
        chi2aplus.reverse()
        x_ = x1 + x2
        y_ = chi2aplus+chi2aminus
    else:
        x2.reverse()
        chi2aminus.reverse()
        x_ = x2 + x1
        y_ = chi2aminus+chi2aplus
    
    #popt, pcov = curve_fit(parabola, x_, y_, p0 = (1000, a+mref*b, min(y_)) )
    #print((10000, b, min(y_)))
    dx = x_[1] - x_[0]
    dy = y_[1] - y_[0]
    approx_a = dy / (dx**2) if dx != 0 else 1e-3
    popt, pcov = curve_fit(parabola, x_, y_, p0 = (approx_a, b, min(y_)) )
    fit_ = parabola(x_, *popt)
    plt.plot(x_, fit_, color='blue')
    
    plt.xlabel(r'$\bm{\alpha_0}$ (GeV)')
    plt.ylabel(r'$\bm{\chi^2}$')

    plt.savefig(f"{plotdir}/{plabel}.png")
    plt.savefig(f"{plotdir}/{plabel}.pdf")
    plt.close()

    return 1/np.sqrt(popt[0])


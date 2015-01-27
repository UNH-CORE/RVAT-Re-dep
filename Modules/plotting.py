# -*- coding: utf-8 -*-
"""
Created on Wed Dec  3 09:15:20 2014

@author: Pete
"""
from Modules.processing import *
import matplotlib
import os

style_file = "Config/ggplot-mod.mplstyle"

try:
    matplotlib.style.use(style_file)
except ValueError:
    matplotlib.style.use(os.path.join("../", style_file))

ylabels = {"mean_u" : r"$U/U_\infty$",
           "std_u" : r"$\sigma_u/U_\infty$",
           "mean_v" : r"$V/U_\infty$",
           "mean_w" : r"$W/U_\infty$",
           "mean_upvp" : r"$\overline{u'v'}/U_\infty^2$",
           "mean_u_diff" : r"$\Delta U$ (\%)",
           "mean_v_diff" : r"$\Delta V$ (\%)",
           "mean_w_diff" : r"$\Delta W$ (\%)"}
           
def plot_trans_wake_profile(quantity, U_infty=0.4, z_H=0.0, save=False, savedir="Figures", 
                            savetype=".pdf", newfig=True, marker="-ok",
                            fill="none", oldwake=False, figsize=(10, 5)):
    """Plots the transverse wake profile of some quantity. These can be
      * mean_u
      * mean_v
      * mean_w
      * std_u
    """
    Re_D = U_infty*D/nu
    label = "{:.1f}e6".format(Re_D/1e6)
    section = "Wake-" + str(U_infty)
    df = pd.read_csv(os.path.join("Data", "Processed", section+".csv"))
    df = df[df.z_H==z_H]
    q = df[quantity]
    y_R = df.y_R
    if newfig:
        plt.figure(figsize=figsize)
    if oldwake:
        plot_old_wake(quantity, y_R)
    if quantity in ["mean_upvp"]:
        unorm = U_infty**2
    else:
        unorm = U_infty
    plt.plot(y_R, q/unorm, marker, markerfacecolor=fill, label=label)
    plt.xlabel(r"$y/R$")
    plt.ylabel(ylabels[quantity])
    plt.tight_layout()
    
def plot_perf_re_dep(save=False, savedir="Figures", savetype=".pdf", errorbars=False,
                     cfd=False, normalize_by="default", dual_xaxes=False, show=True):
    speeds = np.array([0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3])
    cp = np.zeros(len(speeds))
    std_cp = np.zeros(len(speeds))
    exp_unc_cp = np.zeros(len(speeds))
    cd = np.zeros(len(speeds))
    std_cd = np.zeros(len(speeds))
    exp_unc_cd = np.zeros(len(speeds))
    Re_D = speeds*D/1e-6
    for n in range(len(speeds)):
        if speeds[n] in [0.3, 0.5, 0.7, 0.9, 1.1, 1.3]:
            section = "Perf-"+str(speeds[n])
            df = pd.read_csv(os.path.join("Data", "Processed", section+".csv"))
            cp_s = df.mean_cp
            dcp_s = df.exp_unc_cp
            cd_s = df.mean_cd
            dcd_s = df.exp_unc_cd
            cp[n] = np.mean(cp_s)
            exp_unc_cp[n] = np.mean(dcp_s)
            cd[n] = np.mean(cd_s)
            exp_unc_cd[n] = np.mean(dcd_s)
        else:
            section = "Wake-"+str(speeds[n])
            df = pd.read_csv(os.path.join("Data", "Processed", section+".csv"))
            cp_s = df.mean_cp
            dcp_s = df.exp_unc_cp
            cd_s = df.mean_cd
            dcd_s = df.exp_unc_cd
            cp[n], std_cp[n] = np.mean(cp_s), np.std(cp_s)
            cd[n], std_cd[n] = np.mean(cd_s), np.std(cd_s)
            exp_unc_cp[n] = np.mean(dcp_s)
            exp_unc_cd[n] = np.mean(dcd_s)
    plt.figure()
    if normalize_by == "default":
        norm_cp = cp[-4]
        norm_cd = cd[-4]
        norm_cfd = "CFD"
    else:
        norm_cp = normalize_by
        norm_cd = normalize_by
        norm_cfd = normalize_by
    if errorbars:    
        plt.errorbar(Re_D, cp/norm_cp, yerr=exp_unc_cp/2/cp[-4], fmt="-ok",
                     markerfacecolor="none", label="Experiment")
    else:
        plt.plot(Re_D, cp/norm_cp, '-ok', markerfacecolor="none", label="Experiment")
    if cfd:
        plot_cfd_perf("cp", normalize_by=norm_cfd)
    plt.xlabel(r"$Re_D$")
    if normalize_by == "default":
        plt.ylabel(r"$C_P/C_{P0}$")
    else:
        plt.ylabel(r"$C_P$")
#    plt.ylim((0.4, 1.2))
    ax = plt.gca()
    plt.grid(True)
    if dual_xaxes:
        plt.text(1.27e6, 1.11, r"$\times 10^5$", color=r"#555555")
        ax2 = ax.twiny()
        ax.xaxis.get_majorticklocs()
        ticklabs = np.arange(0.2e6, 1.6e6, 0.2e6)
        ticklabs = ticklabs/D*1.9*0.14/1e5
        ticklabs = [str(np.round(ticklab, decimals=1)) for ticklab in ticklabs]
        ax2.set_xticks(ax.xaxis.get_ticklocs())
        ax2.set_xlim((0.2e6, 1.4e6))
        ax2.set_xticklabels(ticklabs)
#        ax2.xaxis.major.formatter.set_powerlimits((0,0)) 
        ax2.set_xlabel(r"$Re_{c, \mathrm{ave}}$")
    if cfd:
        plt.legend(loc=4)
    ax.xaxis.major.formatter.set_powerlimits((0,0)) 
    plt.tight_layout()
    if save:
        plt.savefig(savedir + "/re_dep_cp" + savetype)
    plt.figure()
    if errorbars:
        plt.errorbar(Re_D, cd/norm_cd, yerr=exp_unc_cd/cd[-4]/2, fmt="-ok",
                     markerfacecolor="none", label="Experiment")
    else:
        plt.plot(Re_D, cd/cd[-4], '-ok', markerfacecolor="none", label="Experiment")
    plt.xlabel(r"$Re_D$")
    if normalize_by == "default":
        plt.ylabel(r"$C_D/C_{D0}$")
    else:
        plt.ylabel(r"$C_D$")
    plt.hold(True)
    if cfd:
        plot_cfd_perf("cd", normalize_by=norm_cfd)
#    plt.ylim((0.5,1.1))
    plt.grid(True)
    if cfd:
        plt.legend(loc=4)
    ax = plt.gca()
    ax.xaxis.major.formatter.set_powerlimits((0,0)) 
    plt.tight_layout()
    if save:
        plt.savefig(savedir + "/re_dep_cd" + savetype)
    if show:
        plt.show()
    
def plot_old_wake(quantity, y_R):
    plt.hold(True)
    runs = range(32, 77)
    ind = [run-1 for run in runs]
    f = "../2013.03 VAT/Processed/"+quantity+".npy"
    q = np.load(f)[ind]
    plt.plot(y_R, q, 'xr', label=r"$Re_D=1.0 \times 10^6$", 
             markerfacecolor="none")
             
def plot_cfd_perf(quantity="cp", normalize_by="CFD"):
    Re_D = np.load(cfd_path + "/processed/Re_D.npy")
    q = np.load(cfd_path + "/processed/" + quantity + ".npy")
    if normalize_by=="CFD":
        normval = q[-3]
    else:
        normval = normalize_by
    plt.plot(Re_D, q/normval, "--^k", label="Simulation")
    
def plot_settling(tow_speed):
    """Plot data from the settling experiments."""
    testplan = pd.read_csv("Config/Test plan/Settling.csv")
    nrun = testplan["Run"][testplan["U"] == tow_speed].iloc[0]
    fpath = "Data/Raw/Settling/{}/vecdata.dat".format(nrun)
    data = np.loadtxt(fpath, unpack=True)
    u = data[2] # 2 for x velocity
    t = data[0]*0.005
    uf = u.copy()
    uf[t>80] = ts.sigmafilter(uf[t>80], 4, 1)
    t_std, u_std = ts.runningstd(t, uf, 1000)
    u = ts.smooth(u, 200)
    plt.figure()
    plt.plot(t, u, "k")
    plt.xlabel("t (s)")
    plt.ylabel("$u$ (m/s)")
    plt.tight_layout()
    plt.figure()
    plt.plot(t_std, u_std)
    plt.xlabel("t (s)")
    plt.ylabel(r"$\sigma_u$")
    plt.tight_layout()
    
def plot_perf_curves(subplots=True, save=False, savedir="Figures", savetype=".pdf"):
    """Plots all performance curves."""
    if subplots:
        plt.figure(figsize=(12,5))
        plt.subplot(121)
    PerfCurve(0.4).plotcp(newfig=not subplots, show=False, marker=">")
    PerfCurve(0.6).plotcp(newfig=False, show=False, marker="s")
    PerfCurve(0.8).plotcp(newfig=False, show=False, marker="<")
    PerfCurve(1.0).plotcp(newfig=False, show=False, marker="o")
    PerfCurve(1.2).plotcp(newfig=False, show=False, marker="^")
    if subplots:
        plt.subplot(122)
    PerfCurve(0.4).plotcd(newfig=not subplots, show=False, marker=">")
    PerfCurve(0.6).plotcd(newfig=False, show=False, marker="s")
    PerfCurve(0.8).plotcd(newfig=False, show=False, marker="<")
    PerfCurve(1.0).plotcd(newfig=False, show=False, marker="o")
    PerfCurve(1.2).plotcd(newfig=False, show=False, marker="^")
    plt.legend(("0.4e6", "0.6e6", "0.8e6", "1.0e6", "1.2e6"), 
               loc="lower right", ncol=2)
    plt.show()
    if save:
        if savedir != "":
            savedir += "/"
        plt.savefig(savedir + "perf_curves" + savetype)
    
def plot_wake_profiles(z_H=0.25, save=False, savedir="Figures", savetype=".pdf"):
    """Plots all wake profiles of interest."""
    legendlocs = {"mean_u" : 4,
                  "std_u" : 1,
                  "mean_upvp" : 1}
    for q in ["mean_u", "std_u", "mean_upvp"]:
        plot_trans_wake_profile(q, U_infty=0.4, z_H=z_H, newfig=True, marker="--vb",
                                fill="blue")
        plot_trans_wake_profile(q, U_infty=0.6, z_H=z_H, newfig=False, marker="sk",
                                fill="lightblue")
        plot_trans_wake_profile(q, U_infty=0.8, z_H=z_H, newfig=False, marker="<k",
                                fill="gray")
        plot_trans_wake_profile(q, U_infty=1.0, z_H=z_H, newfig=False, marker="-ok",
                                fill="orange")
        plot_trans_wake_profile(q, U_infty=1.2, z_H=z_H, newfig=False, marker="^k",
                                fill="red")
        plt.legend(loc=legendlocs[q])
        if q == "mean_upvp":
            plt.ylim((-0.015, 0.025))
        if save:
            plt.savefig(os.path.join(savedir, q+savetype))
    plt.show()
    
def plot_meancomboquiv(U_infty=1.0):
    wm = WakeMap(U_infty)
    wm.plot_meancomboquiv()

if __name__ == "__main__":
    pass
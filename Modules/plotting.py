# -*- coding: utf-8 -*-
"""
This module contains classes and functions for plotting data.

"""
from Modules.processing import *
import os

plt.style.use("Config/plotstyle.mplstyle")

ylabels = {"mean_u" : r"$U/U_\infty$",
           "std_u" : r"$\sigma_u/U_\infty$",
           "mean_v" : r"$V/U_\infty$",
           "mean_w" : r"$W/U_\infty$",
           "mean_upvp" : r"$\overline{u'v'}/U_\infty^2$",
           "mean_u_diff" : r"$\Delta U$ (\%)",
           "mean_v_diff" : r"$\Delta V$ (\%)",
           "mean_w_diff" : r"$\Delta W$ (\%)"}
           
           
class PerfCurve(object):
    """Object that represents a performance curve."""
    def __init__(self, tow_speed):
        self.tow_speed = tow_speed
        self.Re_D = tow_speed*D/nu
        self.section = "Perf-{}".format(tow_speed)
        self.raw_data_dir = os.path.join("Data", "Raw", self.section)
        self.df = pd.read_csv(os.path.join("Data", "Processed", self.section+".csv"))
        self.testplan = pd.read_csv(os.path.join("Config", "Test plan", self.section+".csv")) 
        
    def plotcp(self, newfig=True, show=True, save=False, savedir="Figures",
               savetype=".pdf", splinefit=False, marker="o"):
        """Generates power coefficient curve plot."""
        # Check to see if processed data exists and if not, process it
        label = "$Re_D = {:0.1e}$".format(self.Re_D)
        self.tsr = self.df.mean_tsr
        self.cp = self.df.mean_cp
        if newfig:
            plt.figure()
        if splinefit and not True in np.isnan(self.tsr):
            plt.plot(self.tsr, self.cp, marker+"k", markerfacecolor="None", 
                     label=label)
            plt.hold(True)
            tsr_fit = np.linspace(np.min(self.tsr), np.max(self.tsr), 200)
            tck = interpolate.splrep(self.tsr[::-1], self.cp[::-1], s=1e-3)
            cp_fit = interpolate.splev(tsr_fit, tck)
            plt.plot(tsr_fit, cp_fit, "k")
        else:
            if splinefit:
                print("Cannot fit spline. NaN present in array.")
            plt.plot(self.tsr, self.cp, "-"+marker+"k", markerfacecolor="None",
                     label=label)
        plt.xlabel(r"$\lambda$")
        plt.ylabel(r"$C_P$")
        plt.grid(True)
        plt.tight_layout()
        if save:
            plt.savefig(os.path.join(savedir, "cp_vs_tsr" + savetype))
        if show:
            plt.show()
            
    def plotcd(self, newfig=True, show=True, save=False, savedir="Figures",
               savetype=".pdf", splinefit=False, marker="o"):
        """Generates power coefficient curve plot."""
        # Check to see if processed data exists and if not, process it
        label = "$Re_D = {:0.1e}$".format(self.Re_D)
        self.tsr = self.df.mean_tsr
        self.cd = self.df.mean_cd
        if newfig:
            plt.figure()
        if splinefit and not True in np.isnan(self.tsr):
            plt.plot(self.tsr, self.cd, marker+"k", markerfacecolor="None", 
                     label=label)
            plt.hold(True)
            tsr_fit = np.linspace(np.min(self.tsr), np.max(self.tsr), 200)
            tck = interpolate.splrep(self.tsr[::-1], self.cd[::-1], s=1e-3)
            cd_fit = interpolate.splev(tsr_fit, tck)
            plt.plot(tsr_fit, cd_fit, "k")
        else:
            if splinefit:
                print("Cannot fit spline. NaN present in array.")
            plt.plot(self.tsr, self.cd, "-"+marker+"k", markerfacecolor="None",
                     label=label)
        plt.xlabel(r"$\lambda$")
        plt.ylabel(r"$C_D$")
        plt.ylim((0, 1.2))
        plt.grid(True)
        plt.tight_layout()
        if save:
            plt.savefig(os.path.join(savedir, "cd_vs_tsr" + savetype))
        if show:
            plt.show()
        
class WakeProfile(object):
    def __init__(self, tow_speed, z_H, quantity, orientation="horizontal"):
        self.tow_speed = tow_speed
        self.z_H = z_H
        self.section = "Wake-" + str(tow_speed)
        self.testplan = pd.read_csv(os.path.join("Config", "Test plan", 
                                                 self.section+".csv"))
        self.runs = self.testplan.Run[self.testplan["z/H"]==z_H]
        self.quantity = quantity
        self.load()
        
    def load(self):
        """Loads the processed data"""
        self.df = pd.read_csv(os.path.join(processed_data_dir, 
                                           self.section+".csv"))
        self.df = self.df[self.df.z_H==self.z_H]
        self.y_R = self.df["y_R"]
        
    def plot(self, quantity, newfig=True, show=True, save=False, 
             savedir="Figures", savetype=".pdf", linetype='--ok'):
        """Plots some quantity"""
        y_R = self.df["y_R"]
        q = self.df[quantity]
        loc = 1
        if quantity == "mean_u":
            q = q/self.tow_speed
            ylab = r"$U/U_\infty$"
            loc = 3
        if quantity == "mean_w":
            q = q/self.tow_speed
            ylab = r"$U/U_\infty$"
            loc = 4
        if quantity == "mean_v":
            q = q/self.tow_speed
            ylab = r"$V/U_\infty$"
            loc=4
        if quantity == "std_u":
            q = q/self.tow_speed
            ylab = r"$\sigma_u/U_\infty$"
        if quantity is "mean_upvp":
            q = q/(self.tow_speed**2)
            ylab = r"$\overline{u'v'}/U_\infty^2$" 
        if newfig:
            if quantity == "mean_u":
                plt.figure(figsize=(10,5))
            else: plt.figure()
            plt.ylabel(ylab)
            plt.xlabel(r"$y/R$")
        plt.plot(y_R, q, "-.^k", label=r"$Re_D=0.4 \times 10^6$")
        plt.legend(loc=loc)
        plt.tight_layout()
        if show:
            plt.show()
        if save:
            plt.savefig(savedir+quantity+"_Re_dep_exp"+savetype)
    
class WakeMap(object):
    def __init__(self, U_infty):
        self.U_infty = U_infty
        self.z_H = np.array([0.0, 0.125, 0.25, 0.375, 0.5, 0.625])
        self.load()
        self.calc_transport()
        
    def load(self):
        self.df = pd.DataFrame() 
        self.y_R = WakeProfile(self.U_infty, 0, "mean_u").y_R.values
        for z_H in self.z_H:
            wp = WakeProfile(self.U_infty, z_H, "mean_u")
            self.df = self.df.append(wp.df, ignore_index=True)
        self.mean_u = self.df.mean_u
        self.mean_v = self.df.mean_v
        self.mean_w = self.df.mean_w
        self.df["mean_k"] = \
                0.5*(self.df.mean_u**2 + self.df.mean_v**2 + self.df.mean_w**2)
        self.mean_k = self.df.mean_k
        self.grdims = (len(self.z_H), len(self.y_R))
        self.df = self.df.pivot(index="z_H", columns="y_R")
        
    def calc_transport(self):
        """
        Calculates wake tranport terms similar to Bachant and Wosnik (2015)
        "Characterising the near wake of a cross-flow turbine."
        """
        self.calc_mom_transport()
        self.calc_mean_k_grad()
        self.calc_k_prod_mean_diss()
        self.calc_mean_k_turb_trans()
    
    def calc_mean_k_turb_trans(self):
        """Calculates the transport of $K$ by turbulent fluctuations."""
        y, z  = self.y_R*R, self.z_H*H
        self.ddy_uvU = np.zeros(self.grdims)
        self.ddz_uwU = np.zeros(self.grdims)
        self.ddy_vvV = np.zeros(self.grdims)
        self.ddz_vwV = np.zeros(self.grdims)
        self.ddy_vwW = np.zeros(self.grdims)
        self.ddz_wwW = np.zeros(self.grdims)
        for n in range(len(z)):
            self.ddy_uvU[n,:] = \
                fdiff.second_order_diff((self.df.mean_upvp*self.df.mean_u)\
                .iloc[n,:], y)
            self.ddy_vvV[n,:] = \
                fdiff.second_order_diff((self.df.mean_vpvp*self.df.mean_v)\
                .iloc[n,:], y)
            self.ddy_vwW[n,:] = \
                fdiff.second_order_diff((self.df.mean_vpwp*self.df.mean_w)\
                .iloc[n,:], y)
        for n in range(len(y)):
            self.ddz_uwU[:,n] = \
                fdiff.second_order_diff((self.df.mean_upwp*self.df.mean_u)\
                .iloc[:,n], z)
            self.ddz_vwV[:,n] = \
                fdiff.second_order_diff((self.df.mean_vpwp*self.df.mean_v)\
                .iloc[:,n], z)
            self.ddz_wwW[:,n] = \
                fdiff.second_order_diff((self.df.mean_wpwp*self.df.mean_w)\
                .iloc[:,n], z)
        self.mean_k_turb_trans = -0.5*(self.ddy_uvU + \
                                       self.ddz_uwU + \
                                       self.ddy_vvV + \
                                       self.ddz_vwV + \
                                       self.ddy_vwW + \
                                       self.ddz_wwW)
        self.mean_k_turb_trans_y = -0.5*(self.ddy_uvU + \
                                         self.ddy_vvV + \
                                         self.ddy_vwW) # Only ddy terms
        self.mean_k_turb_trans_z = -0.5*(self.ddz_uwU + \
                                         self.ddz_vwV + \
                                         self.ddz_wwW) # Only ddz terms
        
    def calc_k_prod_mean_diss(self):
        """
        Calculates the production of turbulent kinetic energy and dissipation
        from mean shear. Note that the mean streamwise velocity derivatives
        have already been calculated by this point.
        """
        y, z = self.y_R*R, self.z_H*H
        self.dVdy = np.zeros(self.grdims)
        self.dVdz = np.zeros(self.grdims)
        self.dWdy = np.zeros(self.grdims)
        self.dWdz = np.zeros(self.grdims)
        for n in range(len(z)):
            self.dVdy[n,:] = \
                fdiff.second_order_diff(self.df.mean_v.iloc[n,:], y)
            self.dWdy[n,:] = \
                fdiff.second_order_diff(self.df.mean_w.iloc[n,:], y)
        for n in range(len(y)):
            self.dVdz[:,n] = \
                fdiff.second_order_diff(self.df.mean_v.iloc[:,n], z)
            self.dWdz[:,n] = \
                fdiff.second_order_diff(self.df.mean_w.iloc[:,n], z)
        self.k_prod = self.df.mean_upvp*self.dUdy + \
                      self.df.mean_upwp*self.dUdz + \
                      self.df.mean_vpwp*self.dVdz + \
                      self.df.mean_vpwp*self.dWdy + \
                      self.df.mean_vpvp*self.dVdy + \
                      self.df.mean_wpwp*self.dWdz
        self.mean_diss = -2.0*nu*(self.dUdy**2 + self.dUdz**2 + self.dVdy**2 +\
                                  self.dVdz**2 + self.dWdy**2 + self.dWdz**2)
        
    def calc_mean_k_grad(self):
        """Calulates $y$- and $z$-derivatives of $K$."""
        z = self.z_H*H
        y = self.y_R*R
        self.dKdy = np.zeros(self.grdims)
        self.dKdz = np.zeros(self.grdims)
        for n in range(len(z)):
            self.dKdy[n,:] = \
                fdiff.second_order_diff(self.df.mean_k.iloc[n,:], y)
        for n in range(len(y)):
            self.dKdz[:,n] = \
                fdiff.second_order_diff(self.df.mean_k.iloc[:,n], z)

    def calc_mom_transport(self):
        """
        Calculates relevant (and available) momentum transport terms in the 
        RANS equations.
        """
        y = self.y_R*R
        z = self.z_H*H
        self.ddy_upvp = np.zeros(self.grdims)
        self.ddz_upwp = np.zeros(self.grdims)
        self.d2Udy2 = np.zeros(self.grdims)
        self.d2Udz2 = np.zeros(self.grdims)
        self.dUdy = np.zeros(self.grdims)
        self.dUdz = np.zeros(self.grdims)
        for n in range(len(z)):
            self.ddy_upvp[n, :] = \
                fdiff.second_order_diff(self.df.mean_upvp.iloc[n, :], y)
            self.dUdy[n, :] = \
                fdiff.second_order_diff(self.df.mean_u.iloc[n, :], y)
            self.d2Udy2[n, :] = fdiff.second_order_diff(self.dUdy[n, :], y)
        for n in range(len(y)):
            self.ddz_upwp[:, n] = \
                fdiff.second_order_diff(self.df.mean_upwp.iloc[:, n], z)
            self.dUdz[:, n] = \
                fdiff.second_order_diff(self.df.mean_u.iloc[:, n], z)
            self.d2Udz2[:, n] = fdiff.second_order_diff(self.dUdz[:, n], z)
        
    def turb_lines(self, linestyles="solid", linewidth=3, colors="gray"):
        plt.hlines(0.5, -1, 1, linestyles=linestyles, colors="gray",
                   linewidth=linewidth)
        plt.vlines(-1, -0.2, 0.5, linestyles=linestyles, colors="gray",
                   linewidth=linewidth)
        plt.vlines(1, -0.2, 0.5, linestyles=linestyles, colors="gray",
                   linewidth=linewidth)
                   
    def plot_contours(self, quantity, label="", cb_orientation="vertical"):
        """Plots contours of given quantity."""
        plt.figure(figsize=(10,5))
        cs = plt.contourf(self.y_R, self.z_H, quantity, 20,
                          cmap=plt.cm.coolwarm)
        plt.xlabel(r"$y/R$")
        plt.ylabel(r"$z/H$")
        if cb_orientation == "horizontal":
            cb = plt.colorbar(cs, shrink=1, extend="both", 
                              orientation="horizontal", pad=0.3)
        elif cb_orientation == "vertical":
            cb = plt.colorbar(cs, shrink=0.401, extend="both", 
                              orientation="vertical", pad=0.02)
        cb.set_label(label)
        self.turb_lines()
        plt.ylim((0, 0.63))
        ax = plt.axes()
        ax.set_aspect(2)
        plt.yticks([0,0.13,0.25,0.38,0.5,0.63])
        plt.tight_layout()
        
    def plot_mean_u(self, save=False, show=False, savedir="Figures", 
                    savetype=".pdf"):
        """Plot contours of mean streamwise velocity."""
        self.plot_contours(self.df.mean_u/self.U_infty, 
                           label=r"$U/U_\infty$")
        if save:
            plt.savefig(savedir+"/mean_u_cont"+savetype)
        if show:
            self.show()
            
    def plot_k(self, save=False, savetype=".pdf", show=False):
        """Plots contours of turbulence kinetic energy."""
        self.plot_contours(self.df.k/(0.5*self.U_infty**2), 
                           label=r"$k/\frac{1}{2}U_\infty^2$")
        if save:
            plt.savefig("Figures/k_contours_{}{}".format(self.U_infty, savetype))
        if show:
            plt.show()
    
    def plot_meancontquiv(self, save=False, show=False, savedir="Figures",
                          savetype=".pdf", cb_orientation="vertical"):
        """
        Plot contours of mean velocity and vector arrows showing mean
        cross-stream and vertical velocity.
        """
        print("Plotting mean velocity at U_infty = {} m/s".format(self.U_infty))
        plt.figure(figsize=(10,6))
        # Add contours of mean velocity
        cs = plt.contourf(self.y_R, self.z_H, self.df.mean_u/self.U_infty,
                          np.arange(0.15, 1.25, 0.05), cmap=plt.cm.coolwarm)
        if cb_orientation == "horizontal":
            cb = plt.colorbar(cs, shrink=1, extend="both",
                              orientation="horizontal", pad=0.14)
        elif cb_orientation == "vertical":
            cb = plt.colorbar(cs, shrink=0.401, extend="both", 
                              orientation="vertical", pad=0.02)
        cb.set_label(r"$U/U_{\infty}$")
        plt.hold(True)
        # Make quiver plot of v and w velocities
        Q = plt.quiver(self.y_R, self.z_H, self.df.mean_v/self.U_infty, 
                       self.df.mean_w/self.U_infty, width=0.0022,
                       edgecolor="none")
        plt.xlabel(r"$y/R$")
        plt.ylabel(r"$z/H$")
        plt.ylim(-0.2, 0.78)
        plt.xlim(-3.2, 3.2)
        if cb_orientation == "horizontal":
            plt.quiverkey(Q, 0.65, 0.26, 0.1, r"$0.1 U_\infty$",
                          labelpos="E",
                          coordinates="figure",
                          fontproperties={"size": "small"})
        elif cb_orientation == "vertical":
            plt.quiverkey(Q, 0.65, 0.23, 0.1, r"$0.1 U_\infty$",
                          labelpos="E",
                          coordinates="figure",
                          fontproperties={"size": "small"})
        self.turb_lines()
        ax = plt.axes()
        ax.set_aspect(2)
        plt.yticks([0,0.13,0.25,0.38,0.5,0.63])
        plt.grid(True)
        plt.tight_layout()
        if show:
            self.show()
        if save:
            plt.savefig(savedir+"/meancontquiv"+savetype)
    
    def plot_xvorticity(self):
        pass
    
    def plot_diff(self, quantity="mean_u", U_infty_diff=1.0, save=False, 
                  show=False, savedir="Figures", savetype=""):
        wm_diff = WakeMap(U_infty_diff)
        q_ref, q_diff = None, None
        if quantity in ["mean_u", "mean_v", "mean_w"]:
            exec("q_ref = self." + quantity)
            exec("q_diff = wm_diff." + quantity)
            print(q_ref)
        else:
            print("Not a valid quantity")
            return None
        a_diff = (q_ref/self.U_infty - \
                  q_diff/wm_diff.U_infty)#/q_ref/self.U_infty*100
        plt.figure(figsize=(12,3.75))
        cs = plt.contourf(self.y_R, self.z_H, a_diff, 20,
                          cmap=plt.cm.coolwarm)
        cb = plt.colorbar(cs, shrink=1, fraction=0.15,
                          orientation="vertical", pad=0.05)
        cb.set_label(ylabels[quantity+"_diff"])
        plt.xlabel(r"$y/R$")
        plt.ylabel(r"$z/H$")
        plt.axes().set_aspect(2)
        plt.yticks([0,0.13,0.25,0.38,0.5,0.63])
        plt.tight_layout()
        if show:
            self.show()
        if save:
            if savedir: savedir += "/"
            plt.savefig(savedir+"/"+quantity+"_diff"+savetype)
    
    def plot_meancontquiv_diff(self, U_infty_diff, save=False, show=False,
                               savedir="Figures", savetype="", percent=True):
        wm_diff = WakeMap(U_infty_diff)
        mean_u_diff = (self.mean_u/self.U_infty - \
                wm_diff.mean_u/wm_diff.U_infty)
        mean_v_diff = (self.mean_v/self.U_infty - \
                wm_diff.mean_v/wm_diff.U_infty)
        mean_w_diff = (self.mean_w/self.U_infty - \
                wm_diff.mean_w/wm_diff.U_infty)
        if percent:
            mean_u_diff = mean_u_diff/self.mean_u/self.U_infty*100
            mean_v_diff = mean_v_diff/self.mean_v/self.U_infty*100
            mean_w_diff = mean_w_diff/self.mean_w/self.U_infty*100
        plt.figure(figsize=(12,4))
        cs = plt.contourf(self.y_R, self.z_H, mean_u_diff, 20,
                          cmap=plt.cm.coolwarm)
        cb = plt.colorbar(cs, shrink=1, fraction=0.15,
                          orientation="vertical", pad=0.05)
        cb.set_label(r"$\Delta U$ (\%)")
        plt.hold(True)
        # Make quiver plot of v and w velocities
        Q = plt.quiver(self.y_R, self.z_H, mean_v_diff, 
                       mean_w_diff, width=0.0022)
        plt.xlabel(r"$y/R$")
        plt.ylabel(r"$z/H$")
        plt.ylim(-0.2, 0.78)
        plt.xlim(-3.2, 3.2)
        if percent:
            keylen = 100
        else:
            keylen = 0.05
        plt.quiverkey(Q, 0.75, 0.05, keylen, str(keylen),
                      labelpos="E",
                      coordinates="figure",
                      fontproperties={"size": "small"})
        plt.axes().set_aspect(2)
        plt.yticks([0,0.13,0.25,0.38,0.5,0.63])
        plt.tight_layout()
        if show:
            self.show()
        if save:
            if savedir: savedir += "/"
            plt.savefig(savedir+"/meancontquiv_diff"+savetype)
            
    def plot_mean_u_diff_std(self):
        u_ref = 1.0
        mean_u_ref = WakeMap(u_ref).mean_u/u_ref
        std = []
        u_array = np.arange(0.4, 1.4, 0.2)
        for u in u_array:
            wm = WakeMap(u)
            mean_u = wm.mean_u/wm.U_infty
            std.append(np.std((mean_u - mean_u_ref)/mean_u_ref))
        std = np.asarray(std)
        plt.figure()
        plt.plot(u_array, std)
        plt.show()
        
    def show(self):
        plt.show()

           
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
    
def plot_perf_re_dep(save=False, savedir="Figures", savetype=".pdf", 
                     errorbars=False, cfd=False, normalize_by="default", 
                     dual_xaxes=False, show=False):
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
        plt.text(1.31e6, 1.11, "1e5", color=r"#555555")
        ax2 = ax.twiny()
        ax.xaxis.get_majorticklocs()
        ticklabs = np.arange(0.2e6, 1.6e6, 0.2e6)
        ticklabs = ticklabs/D*1.9*0.14/1e5
        ticklabs = [str(np.round(ticklab, decimals=1)) for ticklab in ticklabs]
        ax2.set_xticks(ax.xaxis.get_ticklocs())
        ax2.set_xlim((0.2e6, 1.4e6))
        ax2.set_xticklabels(ticklabs)
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
    ax = plt.gca()
    if dual_xaxes:
        plt.text(1.31e6, 1.035, "1e5", color=r"#555555")
        ax2 = ax.twiny()
        ax.xaxis.get_majorticklocs()
        ticklabs = np.arange(0.2e6, 1.6e6, 0.2e6)
        ticklabs = ticklabs/D*1.9*0.14/1e5
        ticklabs = [str(np.round(ticklab, decimals=1)) for ticklab in ticklabs]
        ax2.set_xticks(ax.xaxis.get_ticklocs())
        ax2.set_xlim((0.2e6, 1.4e6))
        ax2.set_xticklabels(ticklabs)
        ax2.set_xlabel(r"$Re_{c, \mathrm{ave}}$")
    plt.hold(True)
    if cfd:
        plot_cfd_perf("cd", normalize_by=norm_cfd)
#    plt.ylim((0.5,1.1))
    plt.grid(True)
    if cfd:
        plt.legend(loc=4)
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
    
def plot_tare_drag():
    df = pd.read_csv("Data/Processed/Tare drag.csv")
    plt.figure()
    plt.plot(df.tow_speed, df.tare_drag, "-ok")
    plt.xlabel("Tow speed (m/s)")
    plt.ylabel("Tare drag (N)")
    plt.show()
    
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
    
def plot_cp_curve(u_infty, save=False, show=False, savedir="Figures",
                  savetype=".pdf"):
    pc = PerfCurve(u_infty)
    pc.plotcp(save=False, show=False)
    if save:
        savepath = os.path.join(savedir, "cp_vs_tsr_{}".format(u_infty) + savetype)
        plt.savefig(savepath)
    if show:
        plt.show()
    
def plot_perf_curves(subplots=True, save=False, savedir="Figures", 
                     show=False, savetype=".pdf"):
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
    if save:
        if savedir != "":
            savedir += "/"
        plt.savefig(savedir + "perf_curves" + savetype)
    if show:
        plt.show()
    
def plot_wake_profiles(z_H=0.25, save=False, show=False, savedir="Figures", 
                       savetype=".pdf"):
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
    if show:
        plt.show()
    
def plot_meancontquiv(U_infty=1.0, save=False, savetype=".pdf", show=False, 
                      cb_orientation="vertical"):
    wm = WakeMap(U_infty)
    wm.plot_meancontquiv(show=False, cb_orientation=cb_orientation)
    if save:
        plt.savefig("Figures/meancontquiv_{}{}".format(U_infty, savetype))
    if show:
        plt.show()

if __name__ == "__main__":
    pass
# -*- coding: utf-8 -*-
"""This module contains classes and functions for plotting data."""

from .processing import *
import os
from scipy.optimize import curve_fit
import string


ylabels = {"mean_u": r"$U/U_\infty$",
           "std_u": r"$\sigma_u/U_\infty$",
           "mean_v": r"$V/U_\infty$",
           "mean_w": r"$W/U_\infty$",
           "mean_upvp": r"$\overline{u^\prime v^\prime}/U_\infty^2$",
           "mean_u_diff": r"$\Delta U$ (\%)",
           "mean_v_diff": r"$\Delta V$ (\%)",
           "mean_w_diff": r"$\Delta W$ (\%)",
           "k": r"$k/U_\infty^2$"}


class PerfCurve(object):
    """Object that represents a performance curve."""
    def __init__(self, tow_speed):
        self.tow_speed = tow_speed
        self.Re_D = tow_speed*D/nu
        self.section = "Perf-{}".format(tow_speed)
        self.raw_data_dir = os.path.join("Data", "Raw", self.section)
        self.df = pd.read_csv(os.path.join("Data", "Processed", self.section+".csv"))
        self.testplan = pd.read_csv(os.path.join("Config", "Test plan", self.section+".csv"))
        self.label = r"$Re_D = {:.1f} \times 10^6$".format(self.Re_D/1e6)

    def plotcp(self, ax=None, fig=None, save=False, savedir="Figures",
               savetype=".pdf", splinefit=False, **kwargs):
        """Plot mean power coefficient versus tip speed ratio."""
        label = self.label
        self.tsr = self.df.mean_tsr
        self.cp = self.df.mean_cp
        if ax is None:
            fig, ax = plt.subplots()
        if not "marker" in kwargs.keys():
            kwargs["marker"] = "o"
        if splinefit and not True in np.isnan(self.tsr):
            ax.plot(self.tsr, self.cp, label=label, **kwargs)
            tsr_fit = np.linspace(np.min(self.tsr), np.max(self.tsr), 200)
            tck = interpolate.splrep(self.tsr[::-1], self.cp[::-1], s=1e-3)
            cp_fit = interpolate.splev(tsr_fit, tck)
            ax.plot(tsr_fit, cp_fit, color=kwargs["color"])
        else:
            if splinefit:
                print("Cannot fit spline. NaN present in array.")
            ax.plot(self.tsr, self.cp, label=label, **kwargs)
        ax.set_xlabel(r"$\lambda$")
        ax.set_ylabel(r"$C_P$")
        if fig is not None:
            fig.tight_layout()
        if save:
            plt.savefig(os.path.join(savedir, "cp_vs_tsr" + savetype))

    def plotcd(self, ax=None, fig=None, save=False, savedir="Figures",
               savetype=".pdf", splinefit=False, **kwargs):
        """Plot mean power coefficient versus tip speed ratio."""
        label = self.label
        self.tsr = self.df.mean_tsr
        self.cd = self.df.mean_cd
        if ax is None:
            fig, ax = plt.subplots()
        if not "marker" in kwargs.keys():
            kwargs["marker"] = "o"
        if splinefit and not True in np.isnan(self.tsr):
            ax.plot(self.tsr, self.cd, label=label, **kwargs)
            tsr_fit = np.linspace(np.min(self.tsr), np.max(self.tsr), 200)
            tck = interpolate.splrep(self.tsr[::-1], self.cd[::-1], s=1e-3)
            cd_fit = interpolate.splev(tsr_fit, tck)
            ax.plot(tsr_fit, cd_fit, color=kwargs["color"])
        else:
            if splinefit:
                print("Cannot fit spline. NaN present in array.")
            ax.plot(self.tsr, self.cd, label=label, **kwargs)
        ax.set_xlabel(r"$\lambda$")
        ax.set_ylabel(r"$C_D$")
        if fig is not None:
            fig.tight_layout()
        if save:
            plt.savefig(os.path.join(savedir, "cd_vs_tsr" + savetype))


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
                plt.figure(figsize=(7.5, 3.75))
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
        """Calculate the production of turbulent kinetic energy and dissipation
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
        self.dUdx = -self.dVdy - self.dWdz
        self.k_prod = self.df.mean_upvp*self.dUdy + \
                      self.df.mean_upwp*self.dUdz + \
                      self.df.mean_vpwp*self.dVdz + \
                      self.df.mean_vpwp*self.dWdy + \
                      self.df.mean_vpvp*self.dVdy + \
                      self.df.mean_wpwp*self.dWdz
        self.mean_diss = -2.0*nu*(self.dUdy**2 + self.dUdz**2 + self.dVdy**2 +\
                                  self.dVdz**2 + self.dWdy**2 + self.dWdz**2)

    def calc_mean_k_grad(self):
        """Calulate $y$- and $z$-derivatives of $K$."""
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
        """Calculate relevant (and available) momentum transport terms in the
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

    def turb_lines(self, linestyles="solid", linewidth=2, color="gray"):
        plt.hlines(0.5, -1, 1, linestyles=linestyles, colors=color,
                   linewidth=linewidth)
        plt.vlines(-1, -0.2, 0.5, linestyles=linestyles, colors=color,
                   linewidth=linewidth)
        plt.vlines(1, -0.2, 0.5, linestyles=linestyles, colors=color,
                   linewidth=linewidth)

    def plot_contours(self, quantity, label="", cb_orientation="vertical",
                      newfig=True, levels=None):
        """Plots contours of given quantity."""
        if newfig:
            plt.figure(figsize=(7.5, 2.0))
        cs = plt.contourf(self.y_R, self.z_H, quantity, 20,
                          cmap=plt.cm.coolwarm, levels=levels)
        plt.xlabel(r"$y/R$")
        plt.ylabel(r"$z/H$")
        if cb_orientation == "horizontal":
            cb = plt.colorbar(cs, shrink=1, extend="both",
                              orientation="horizontal", pad=0.3)
        elif cb_orientation == "vertical":
            cb = plt.colorbar(cs, shrink=1, extend="both",
                              orientation="vertical", pad=0.02)
        cb.set_label(label)
        self.turb_lines(color="black")
        plt.ylim((0, 0.63))
        ax = plt.axes()
        ax.set_aspect(2)
        plt.yticks([0,0.13,0.25,0.38,0.5,0.63])
        plt.tight_layout()

    def plot_mean_u(self, save=False, show=False, savedir="Figures",
                    savetype=".pdf", newfig=True):
        """Plot contours of mean streamwise velocity."""
        self.plot_contours(self.df.mean_u/self.U_infty,
                           label=r"$U/U_\infty$", newfig=newfig)
        if save:
            plt.savefig(savedir+"/mean_u_cont"+savetype)
        if show:
            self.show()

    def plot_k(self, newfig=True, save=False, savetype=".pdf", show=False):
        """Plots contours of turbulence kinetic energy."""
        self.plot_contours(self.df.k/(self.U_infty**2),
                           newfig=newfig,
                           label=r"$k/U_\infty^2$",
                           levels=np.linspace(0, 0.09, num=19))
        if save:
            label = str(self.U_infty).replace(".", "")
            plt.savefig("Figures/k_contours_{}{}".format(label, savetype))
        if show:
            plt.show()

    def plot_meancontquiv(self, save=False, show=False, savedir="Figures",
                          savetype=".pdf", cb_orientation="vertical",
                          newfig=True):
        """Plot contours of mean velocity and vector arrows showing mean
        cross-stream and vertical velocity.
        """
        if newfig:
            plt.figure(figsize=(7.5, 2.625))
        # Add contours of mean velocity
        cs = plt.contourf(self.y_R, self.z_H, self.df.mean_u/self.U_infty,
                          np.arange(0.15, 1.25, 0.05), cmap=plt.cm.coolwarm)
        if cb_orientation == "horizontal":
            cb = plt.colorbar(cs, shrink=1, extend="both",
                              orientation="horizontal", pad=0.14)
        elif cb_orientation == "vertical":
            cb = plt.colorbar(cs, shrink=0.83, extend="both",
                              orientation="vertical", pad=0.02)
        cb.set_label(r"$U/U_{\infty}$")
        # Make quiver plot of v and w velocities
        Q = plt.quiver(self.y_R, self.z_H, self.df.mean_v/self.U_infty,
                       self.df.mean_w/self.U_infty, width=0.0022, scale=3,
                       edgecolor="none")
        plt.xlabel(r"$y/R$")
        plt.ylabel(r"$z/H$")
        plt.ylim(-0.2, 0.78)
        plt.xlim(-3.2, 3.2)
        if cb_orientation == "horizontal":
            plt.quiverkey(Q, 0.65, 0.26, 0.1, r"$0.1 U_\infty$",
                          labelpos="E", coordinates="figure")
        elif cb_orientation == "vertical":
            plt.quiverkey(Q, 0.65, 0.088, 0.1, r"$0.1 U_\infty$",
                          labelpos="E", coordinates="figure",
                          fontproperties={"size": "small"})
        self.turb_lines()
        ax = plt.axes()
        ax.set_aspect(2)
        plt.yticks([0, 0.13, 0.25, 0.38, 0.5, 0.63])
        plt.tight_layout()
        if save:
            label = str(self.U_infty).replace(".", "")
            plt.savefig(savedir+"/meancontquiv_{}{}".format(label, savetype))
        if show:
            self.show()

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
        plt.figure(figsize=(7.5, 2.34))
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
                               savedir="Figures", savetype="", percent=True,
                               cb_orientation="vertical"):
        wm_diff = WakeMap(U_infty_diff)
        mean_u_diff = (self.df.mean_u/self.U_infty - \
                wm_diff.df.mean_u/wm_diff.U_infty)
        mean_v_diff = (self.df.mean_v/self.U_infty - \
                wm_diff.df.mean_v/wm_diff.U_infty)
        mean_w_diff = (self.df.mean_w/self.U_infty - \
                wm_diff.df.mean_w/wm_diff.U_infty)
        if percent:
            mean_u_diff = mean_u_diff/self.df.mean_u/self.U_infty*100
            mean_v_diff = mean_v_diff/self.df.mean_v/self.U_infty*100
            mean_w_diff = mean_w_diff/self.df.mean_w/self.U_infty*100
        plt.figure(figsize=(7.5, 2.625))
        cs = plt.contourf(self.y_R, self.z_H, mean_u_diff, 20,
                          cmap=plt.cm.coolwarm)
        if cb_orientation == "horizontal":
            cb = plt.colorbar(cs, shrink=1, extend="both",
                              orientation="horizontal", pad=0.14)
        elif cb_orientation == "vertical":
            cb = plt.colorbar(cs, shrink=0.785, extend="both",
                              orientation="vertical", pad=0.02)
        cb.set_label(r"$\Delta U$ (\%)")
        # Make quiver plot of v and w velocities
        Q = plt.quiver(self.y_R, self.z_H, mean_v_diff,
                       mean_w_diff, width=0.0022, edgecolor="none", scale=3)
        if cb_orientation == "horizontal":
            plt.quiverkey(Q, 0.65, 0.26, 0.1, r"$0.1 U_\infty$",
                          labelpos="E", coordinates="figure")
        elif cb_orientation == "vertical":
            plt.quiverkey(Q, 0.65, 0.08, 0.1, r"$0.1 U_\infty$",
                          labelpos="E", coordinates="figure")
        plt.xlabel(r"$y/R$")
        plt.ylabel(r"$z/H$")
        plt.ylim(-0.2, 0.78)
        plt.xlim(-3.2, 3.2)
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


class WakeMapDiff(WakeMap):
    """Object representing the difference between two wake maps. Quantities are
    calculated as `wm1 - wm2`.
    """
    def __init__(self, U1, U2):
        WakeMap.__init__(self, 1.0)
        self.U1 = U1
        self.U2 = U2
        self.wm1 = WakeMap(U1)
        self.wm2 = WakeMap(U2)
        self.df = self.wm1.df - self.wm2.df
        self.mean_u = self.wm1.df.mean_u/self.wm1.U_infty \
                    - self.wm2.df.mean_u/self.wm2.U_infty

    def plot_mean_u(self):
        self.plot_contours(self.mean_u, label="$U_{\mathrm{diff}}$")


def plot_trans_wake_profile(ax=None, quantity="mean_u", U_infty=0.4, z_H=0.0,
                            save=False, marker="-ok", color="black",
                            oldwake=False, figsize=(7.5, 3.75)):
    """Plot the transverse wake profile of some quantity.

    These can be
      * mean_u
      * mean_v
      * mean_w
      * std_u
    """
    Re_D = U_infty*D/nu
    label = "{:.1f}e6".format(Re_D/1e6)
    section = f"Wake-{U_infty:.1f}"
    df = pd.read_csv(os.path.join("Data", "Processed", section+".csv"))
    df = df[df.z_H==z_H]
    q = df[quantity]
    y_R = df.y_R
    if quantity in ["mean_upvp", "k"]:
        unorm = U_infty**2
    else:
        unorm = U_infty
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    ax.plot(y_R, q/unorm, marker, color=color, label=label)
    ax.set_xlabel(r"$y/R$")
    ax.set_ylabel(ylabels[quantity])
    ax.grid(True)
    try:
        fig.tight_layout()
    except UnboundLocalError:
        pass


def plot_perf_re_dep(ax1=None, ax2=None, save=False, savedir="Figures",
                     savetype=".pdf", errorbars=False, subplots=True,
                     normalize_by=1.0, dual_xaxes=True, power_law=False,
                     label_subplots=True, **kwargs):
    """Plot the Reynolds number dependence of power and drag coefficients."""
    if not subplots:
        label_subplots = False
    if not "marker" in kwargs.keys():
        kwargs["marker"] = "o"
    speeds = np.array([0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3])
    cp = np.zeros(len(speeds))
    std_cp = np.zeros(len(speeds))
    exp_unc_cp = np.zeros(len(speeds))
    cd = np.zeros(len(speeds))
    std_cd = np.zeros(len(speeds))
    exp_unc_cd = np.zeros(len(speeds))
    Re_D = speeds*D/1e-6
    Re_c = Re_D/D*chord*1.9
    for n in range(len(speeds)):
        if speeds[n] in [0.3, 0.5, 0.7, 0.9, 1.1, 1.3]:
            section = "Perf-"+str(speeds[n])
        else:
            section = "Wake-"+str(speeds[n])
        df = pd.read_csv(os.path.join("Data", "Processed", section+".csv"))
        cp[n] = df.mean_cp.mean()
        cd[n] = df.mean_cd.mean()
        if errorbars:
            exp_unc_cp[n] = ts.calc_multi_exp_unc(df.sys_unc_cp, df.n_revs,
                                                  df.mean_cp,
                                                  df.std_cp_per_rev, df.dof_cp,
                                                  confidence=0.95)
            exp_unc_cd[n] = ts.calc_multi_exp_unc(df.sys_unc_cd, df.n_revs,
                                                  df.mean_cd,
                                                  df.std_cd_per_rev, df.dof_cd,
                                                  confidence=0.95)
    df = pd.DataFrame()
    df["Re_D"] = Re_D
    df["Re_c_ave"] = Re_c
    df["mean_cp"] = cp
    df["mean_cd"] = cd
    df.to_csv("Data/Processed/Perf-tsr_0.csv", index=False)
    if ax1 is None and ax2 is None and subplots:
        fig1, (ax1, ax2) = plt.subplots(figsize=(7.5, 3.5), nrows=1, ncols=2)
    elif ax1 is None and ax2 is None and not subplots:
        fig1, ax1 = plt.subplots()
        fig2, ax2 = plt.subplots()
    if normalize_by == "C_P_0":
        norm_cp = cp[-4]
        norm_cd = cd[-4]
    else:
        norm_cp = normalize_by
        norm_cd = normalize_by
    if errorbars:
        ax1.errorbar(Re_D, cp/norm_cp, yerr=exp_unc_cp/norm_cp,
                     label="Experiment", **kwargs)
    else:
        ax1.plot(Re_D, cp/norm_cp, label="Experiment", **kwargs)
    ax1.set_xlabel(r"$Re_D$")
    if normalize_by == "default":
        ax1.set_ylabel(r"$C_P/C_{P_0}$")
    else:
        ax1.set_ylabel(r"$C_P$")
    if dual_xaxes:
        twiny_sci_label(ax=ax1, power=5, subplots=subplots)
        ax12 = ax1.twiny()
        ticklabs = ax1.xaxis.get_majorticklocs()
        ticklabs = ticklabs/D*1.9*0.14/1e5
        ticklabs = [str(np.round(ticklab, decimals=1)) for ticklab in ticklabs]
        ax12.set_xticks(ax1.xaxis.get_ticklocs())
        ax12.set_xlim(ax1.get_xlim())
        ax12.set_xticklabels(ticklabs)
        ax12.set_xlabel(r"$Re_{c, \mathrm{ave}}$")
        ax12.grid(False)
    if power_law:
        # Calculate power law fits for quantities
        def func(x, a, b):
            return a*x**b
        coeffs_cd, covar_cd = curve_fit(func, Re_c, cd)
        coeffs_cp, covar_cp = curve_fit(func, Re_c, cp)
        print("Power law fits:")
        print("C_P = {:.3f}*Re_c**{:.3f}".format(coeffs_cp[0], coeffs_cp[1]))
        print("C_D = {:.3f}*Re_c**{:.3f}".format(coeffs_cd[0], coeffs_cd[1]))
        Re_D_curve = np.linspace(0.3e6, 1.3e6)
        cp_power_law = coeffs_cp[0]*(Re_D_curve/D*chord*1.9)**coeffs_cp[1]
        ax1.plot(Re_D_curve, cp_power_law, "--k",
                 label=r"${:.3f}Re_c^{{ {:.3f} }}$".format(coeffs_cp[0],
                 coeffs_cp[1]))
        ax1.legend(loc="lower right")
    ax1.set_ylim((0.14/normalize_by, 0.28/normalize_by))
    ax1.xaxis.major.formatter.set_powerlimits((0, 0))
    ax1.grid(True)
    if label_subplots:
        label_subplot(ax1, text="(a)")
    try:
        fig1.tight_layout()
    except UnboundLocalError:
        pass
    if save and not subplots:
        fig1.savefig(savedir + "/re_dep_cp" + savetype)
    if errorbars:
        ax2.errorbar(Re_D, cd/norm_cd, yerr=exp_unc_cd/norm_cd,
                     label="Experiment", **kwargs)
    else:
        ax2.plot(Re_D, cd/norm_cd, label="Experiment", **kwargs)
    ax2.set_xlabel(r"$Re_D$")
    if normalize_by == "default":
        ax2.set_ylabel(r"$C_D/C_{D_0}$")
    else:
        ax2.set_ylabel(r"$C_D$")
    if dual_xaxes:
        twiny_sci_label(ax=ax2, power=5, subplots=subplots)
        ax22 = ax2.twiny()
        ticklabs = ax2.xaxis.get_majorticklocs()
        ticklabs = ticklabs/D*1.9*0.14/1e5
        ticklabs = [str(np.round(ticklab, decimals=1)) for ticklab in ticklabs]
        ax22.set_xticks(ax2.xaxis.get_ticklocs())
        ax22.set_xlim(ax2.get_xlim())
        ax22.set_xticklabels(ticklabs)
        ax22.set_xlabel(r"$Re_{c, \mathrm{ave}}$")
        ax22.grid(False)
    ax2.set_ylim((0.82/norm_cd, 0.96/norm_cd))
    ax2.xaxis.major.formatter.set_powerlimits((0,0))
    if label_subplots:
        label_subplot(text="(b)")
    ax2.grid(True)
    if power_law:
        cd_power_law = coeffs_cd[0]*(Re_D_curve/D*chord*1.9)**coeffs_cd[1]
        ax2.plot(Re_D_curve, cd_power_law, "--k",
                 label=r"${:.3f}Re_c^{{ {:.3f} }}$".format(coeffs_cd[0],
                 coeffs_cd[1]))
        ax2.legend(loc="lower right")
    try:
        fig1.tight_layout()
        fig2.tight_layout()
    except UnboundLocalError:
        pass
    if save:
        if subplots:
            fig1.savefig(savedir + "/perf_re_dep" + savetype)
        else:
            fig2.savefig(savedir + "/re_dep_cd" + savetype)


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


def plot_perf_curves(ax1=None, ax2=None, subplots=True, save=False,
                     savedir="Figures", savetype=".pdf", **kwargs):
    """Plot all performance curves."""
    if ax1 is None and ax2 is None and subplots:
        fig1, (ax1, ax2) = plt.subplots(figsize=(7.5, 3.07), nrows=1, ncols=2)
    elif ax1 is None and ax2 is None and not subplots:
        fig1, ax1 = plt.subplots()
        fig2, ax2 = plt.subplots()
    speeds = np.round(np.arange(0.4, 1.3, 0.2), decimals=1)
    cm = plt.cm.coolwarm
    colors = [cm(int(n/4*256)) for n in range(len(speeds))]
    markers = [">", "s", "<", "o", "^"]
    for speed, color, marker, in zip(speeds, colors, markers):
        PerfCurve(speed).plotcp(ax=ax1, marker=marker, color=color, **kwargs)
    ax1.grid(True)
    if not subplots:
        ax1.legend(loc="best", ncol=2)
    try:
        fig1.tight_layout()
    except UnboundLocalError:
        pass
    if save and not subplots:
        fig1.savefig(os.path.join(savedir, "cp_curves" + savetype))
    for speed, color, marker, in zip(speeds, colors, markers):
        PerfCurve(speed).plotcd(ax=ax2, marker=marker, color=color, **kwargs)
    ax2.legend(loc="lower right", ncol=2)
    ax2.set_ylim((0, 1.2))
    ax2.grid(True)
    try:
        fig1.tight_layout()
        fig2.tight_layout()
    except UnboundLocalError:
        pass
    if save:
        if subplots:
            fig1.savefig(os.path.join(savedir, "perf_curves" + savetype))
        else:
            fig2.savefig(os.path.join(savedir, "cd_curves" + savetype))


def plot_wake_profiles(z_H=0.0, save=False, show=False, savedir="Figures",
                       quantities=["mean_u", "k"], figsize=(7.5, 3.25),
                       savetype=".pdf", subplots=True, label_subplots=True):
    """Plot wake profiles for all Re."""
    tow_speeds = np.arange(0.4, 1.3, 0.2)
    cm = plt.cm.coolwarm
    colors = [cm(int(n/4*256)) for n in range(len(tow_speeds))]
    markers = ["--v", "s", "<", "-o", "^"]
    letters = list(string.ascii_lowercase)[:len(quantities)]
    if subplots:
        fig, ax = plt.subplots(figsize=figsize, nrows=1, ncols=len(quantities))
    else:
        ax = [None]*len(quantities)
        label_subplots = False
    for a, q, letter in zip(ax, quantities, letters):
        if not subplots:
            fig, a = plt.subplots(figsize=figsize)
        for U, marker, color in zip(tow_speeds, markers, colors):
            plot_trans_wake_profile(ax=a, quantity=q, U_infty=U, z_H=z_H,
                                    marker=marker, color=color)
        if q == quantities[0] or not subplots:
            a.legend(loc="lower left")
        if q == "mean_upvp":
            a.set_ylim((-0.015, 0.025))
        fig.tight_layout()
        if label_subplots:
            label_subplot(ax=a, text="({})".format(letter))
        if save and not subplots:
            fig.savefig(os.path.join(savedir, q + "_profiles" + savetype))
    if save and subplots:
        fig.savefig(os.path.join(savedir,
                                 "_".join(quantities) + "_profiles" + savetype))


def plot_meancontquiv(U_infty=1.0, save=False, savetype=".pdf", show=False,
                      cb_orientation="vertical"):
    wm = WakeMap(U_infty)
    wm.plot_meancontquiv(save=save, show=show, savetype=savetype,
                         cb_orientation=cb_orientation)


def plot_all_meancontquiv(save=False, savetype=".pdf", show=False):
    """Plot all mean velocity contours/quivers."""
    for n, U in enumerate([0.4, 0.6, 0.8, 1.0, 1.2]):
        WakeMap(U).plot_meancontquiv(newfig=True, save=save, savetype=savetype)
    if show:
        plt.show()


def plot_all_kcont(save=False, savetype=".pdf"):
    """Plot contours of turbulence kinetic energy for all Reynolds numbers
    tested.
    """
    for n, U in enumerate([0.4, 0.6, 0.8, 1.0, 1.2]):
        WakeMap(U).plot_k(save=save, savetype=savetype)


def make_k_bar_graph(save=False, savetype=".pdf", show=False,
                     print_analysis=True):
    """Makes a bar graph from the mean kinetic energy transport terms for four
    Reynolds numbers.
    """
    names = [r"$y$-adv.", r"$z$-adv.", r"$y$-turb.", r"$z$-turb.",
             r"$k$-prod.", r"Mean diss. $(\times 10^3)$"]
    plt.figure(figsize=(7.5, 3.2))
    cm = plt.cm.coolwarm
    for n, U in enumerate([0.4, 0.6, 0.8, 1.0, 1.2]):
        Re_D = U*D/nu
        wm = WakeMap(U)
        tty, ttz = wm.mean_k_turb_trans_y, wm.mean_k_turb_trans_z
        kprod, meandiss = wm.k_prod, wm.mean_diss
        dKdy, dKdz = wm.dKdy, wm.dKdz
        y_R, z_H = wm.y_R, wm.z_H
        meanu, meanv, meanw = wm.df.mean_u, wm.df.mean_v, wm.df.mean_w
        quantities = [ts.average_over_area(-2*meanv/meanu*dKdy/(0.5*U**2)*D, y_R, z_H),
                      ts.average_over_area(-2*meanw/meanu*dKdz/(0.5*U**2)*D, y_R, z_H),
                      ts.average_over_area(2*tty/meanu/(0.5*U**2)*D, y_R, z_H),
                      ts.average_over_area(2*ttz/meanu/(0.5*U**2)*D, y_R, z_H),
                      ts.average_over_area(2*kprod/meanu/(0.5*U**2)*D, y_R, z_H),
                      ts.average_over_area(2*meandiss/meanu/(0.5*U**2)*D*1e3, y_R, z_H)]
        ax = plt.gca()
        color = cm(int(n/4*256))
        ax.bar(np.arange(len(names))+n*0.15, quantities, color=color, edgecolor="black",
               hatch=None, width=0.15,
               label=r"$Re_D={:.1f}\times 10^6$".format(Re_D/1e6))
        if print_analysis:
            quantities[-1] /= 1e3
            print("K recovery rate at {:.1f} m/s: {:.2f} (%/D)".format(U,
                  np.sum(quantities)*100))
    ax.set_xticks(np.arange(len(names)) + 5*.15/2)
    ax.set_xticklabels(names)
    plt.hlines(0, 0, len(names), color="black")
    plt.ylabel(r"$\frac{K \, \mathrm{ transport}}{UK_\infty D^{-1}}$")
    plt.legend(loc="upper right", ncol=2)
    plt.tight_layout()
    if save:
        plt.savefig("Figures/K_trans_bar_graph" + savetype)
    if show:
        plt.show()


def make_mom_bar_graph(ax=None, save=False, savetype=".pdf",
                       print_analysis=True):
    """Create a bar graph of terms contributing to dU/dx:
      * Cross-stream advection
      * Vertical advection
      * Cross-stream Re stress gradient
      * Vertical Re stress gradient
      * Cross-steam diffusion
      * Vertical diffusion
    """
    names = [r"$-V \frac{\partial U}{\partial y}$",
             r"$-W \frac{\partial U}{\partial z}$",
             r"$-\frac{\partial}{\partial y} \overline{u^\prime v^\prime}$",
             r"$-\frac{\partial}{\partial z} \overline{u^\prime w^\prime}$",
             r"$\nu \frac{\partial^2 U}{\partial y^2} (\times 10^3)$",
             r"$\nu \frac{\partial^2 U}{\partial z^2} (\times 10^3)$"]
    if ax is None:
        fig, ax = plt.subplots(figsize=(7.5, 3.2))
    cm = plt.cm.coolwarm
    for n, U in enumerate([0.4, 0.6, 0.8, 1.0, 1.2]):
        wm = WakeMap(U)
        dUdy = wm.dUdy
        dUdz = wm.dUdz
        tty = wm.ddy_upvp
        ttz = wm.ddz_upwp
        d2Udy2 = wm.d2Udy2
        d2Udz2 = wm.d2Udz2
        meanu, meanv, meanw = wm.df.mean_u, wm.df.mean_v, wm.df.mean_w
        y_R, z_H = wm.y_R, wm.z_H
        quantities = [ts.average_over_area(-2*meanv*dUdy/meanu/U*D, y_R, z_H),
                      ts.average_over_area(-2*meanw*dUdz/meanu/U*D, y_R, z_H),
                      ts.average_over_area(-2*tty/meanu/U*D, y_R, z_H),
                      ts.average_over_area(-2*ttz/meanu/U*D, y_R, z_H),
                      ts.average_over_area(2*nu*d2Udy2/meanu/U*D*1e3, y_R, z_H),
                      ts.average_over_area(2*nu*d2Udz2/meanu/U*D*1e3, y_R, z_H)]
        dUdx = ts.average_over_area(2*wm.dUdx/U*D, y_R, z_H)
        color = cm(int(n/4*256))
        ax.bar(np.arange(len(names)) + n*.15, quantities, color=color,
               width=0.15, edgecolor="black",
               label=r"$Re_D={:.1f}\times 10^6$".format(U*D/nu/1e6))
        if print_analysis:
            quantities[4] /= 1e3
            quantities[5] /= 1e3
            print("U recovery rate at {:.1f} m/s: {:.2f} (%/D)".format(U,
                  np.sum(quantities)*100))
    ax.set_xticks(np.arange(len(names)) + 5*.15/2)
    ax.set_xticklabels(names)
    ax.hlines(0, 0, len(names), color="black")
    ax.set_ylabel(r"$\frac{U \, \mathrm{ transport}}{UU_\infty D^{-1}}$")
    ax.legend(loc="upper right", ncol=2)
    try:
        fig.tight_layout()
    except UnboundLocalError:
        pass
    if save:
        plt.savefig("Figures/mom_bar_graph"+savetype)


def plot_wake_trans_totals(ax=None, save=False, savetype=".pdf", **kwargs):
    """Plot totals for wake transport quantities for all Reynolds numbers
    tested, both for the momentum and kinetic energy.
    """
    momentum_totals = []
    energy_totals = []
    speeds = [0.4, 0.6, 0.8, 1.0, 1.2]
    Re_D = np.array(speeds)*D/nu
    for U in speeds:
        wm = WakeMap(U)
        dUdy = wm.dUdy
        dUdz = wm.dUdz
        tty = wm.ddy_upvp
        ttz = wm.ddz_upwp
        d2Udy2 = wm.d2Udy2
        d2Udz2 = wm.d2Udz2
        meanu, meanv, meanw = wm.df.mean_u, wm.df.mean_v, wm.df.mean_w
        y_R, z_H = wm.y_R, wm.z_H
        quantities = [ts.average_over_area(-2*meanv*dUdy/meanu/U*D, y_R, z_H),
                      ts.average_over_area(-2*meanw*dUdz/meanu/U*D, y_R, z_H),
                      ts.average_over_area(-2*tty/meanu/U*D, y_R, z_H),
                      ts.average_over_area(-2*ttz/meanu/U*D, y_R, z_H),
                      ts.average_over_area(2*nu*d2Udy2/meanu/U*D, y_R, z_H),
                      ts.average_over_area(2*nu*d2Udz2/meanu/U*D, y_R, z_H)]
        momentum_totals.append(np.sum(quantities))
        tty, ttz = wm.mean_k_turb_trans_y, wm.mean_k_turb_trans_z
        kprod, meandiss = wm.k_prod, wm.mean_diss
        dKdy, dKdz = wm.dKdy, wm.dKdz
        y_R, z_H = wm.y_R, wm.z_H
        meanu, meanv, meanw = wm.df.mean_u, wm.df.mean_v, wm.df.mean_w
        quantities = [ts.average_over_area(-2*meanv/meanu*dKdy/(0.5*U**2)*D, y_R, z_H),
                      ts.average_over_area(-2*meanw/meanu*dKdz/(0.5*U**2)*D, y_R, z_H),
                      ts.average_over_area(2*tty/meanu/(0.5*U**2)*D, y_R, z_H),
                      ts.average_over_area(2*ttz/meanu/(0.5*U**2)*D, y_R, z_H),
                      ts.average_over_area(2*kprod/meanu/(0.5*U**2)*D, y_R, z_H),
                      ts.average_over_area(2*meandiss/meanu/(0.5*U**2)*D, y_R, z_H)]
        energy_totals.append(np.sum(quantities))
    if ax is None:
        fig, ax = plt.subplots()
    ax.plot(Re_D, momentum_totals, marker="o", label="$U$", **kwargs)
    ax.plot(Re_D, energy_totals, marker="s", label="$K$", **kwargs)
    ax.set_xlabel("$Re_D$")
    ax.set_ylabel("Normalized total transport")
    ax.legend(loc="best")
    ax.grid(True)
    try:
        fig.tight_layout()
    except UnboundLocalError:
        pass
    if save:
        plt.savefig("Figures/wake_trans_totals" + savetype)


def plot_vel_spec(U_infty, y_R, z_H, n_band_ave=4, plot_conf_int=False,
                  show=False, newfig=True, plot_lines=True, color="black"):
    """Plot the cross-stream velocity spectrum (normalized by the free stream
    velocity) for a single run.

    Any NaNs in the velocity data are replaced with the mean.

    Bachant and Wosnik (2015, JoT) used locations y/R = (-1, 1.5) and
    z/H = 0.25 to compare spectra with high and low levels of turbulence,
    respectively.
    """
    print("Plotting cross-stream velocity spectrum at ({}, {})".format(y_R,
                                                                       z_H))
    # Find index for the desired parameters
    s_name = "Wake-{:.1f}".format(U_infty)
    tp = Section(s_name).test_plan
    tp = tp[tp["y/R"] == y_R]
    n = int(tp[tp["z/H"] == z_H].iloc[0]["Run"])
    r = Run(s_name, n)
    v = r.v
    print("Replacing {} datapoints with mean".format(r.nbadv))
    v[np.isnan(v)] = r.mean_v
    f, spec = ts.psd(r.time_vec, v/U_infty, window="Hanning",
                     n_band_average=n_band_ave)
    f_turbine = r.mean_tsr*U_infty/R/(2*np.pi)
    # Find maximum frequency and its relative strength
    f_max = f[spec==spec.max()][0]
    strength = np.max(spec)/r.std_v**2*(f[1] - f[0])
    print("Strongest frequency f/f_turbine: {:.3f}".format(f_max/f_turbine))
    print("Spectral concentration: {:.3f}".format(strength))
    if newfig:
        plt.figure()
    plt.loglog(f/f_turbine, spec, color=color,
               label=r"$Re_D = {:.1f} \times 10^6$".format(U_infty))
    plt.xlim((0, 50))
    plt.xlabel(r"$f/f_{\mathrm{turbine}}$")
    plt.ylabel(r"Spectral density")
    # Should the spectrum be normalized?
    if plot_lines:
        f_line = np.linspace(10,40)
        spec_line = f_line**(-5./3)*0.5*1e-2
        plt.loglog(f_line, spec_line, "black")
        plt.ylim((1e-8, 1e-1))
        plot_vertical_lines([1, 3, 6, 9], color="lightgray")
    if plot_conf_int:
        dof = n_band_ave*2
        chi2 = scipy.stats.chi2.interval(alpha=0.95, df=dof)
        y1 = dof*spec/chi2[1]
        y2 = dof*spec/chi2[0]
        plt.fill_between(f/f_turbine, y1, y2, facecolor=color, alpha=0.3)
    plt.grid(True)
    plt.tight_layout()
    if show:
        plt.show()


def plot_multi_spec(n_band_ave=4, plot_conf_int=False, save=False, show=False,
                    savetype=".pdf"):
    """Plot the cross-stream velocity spectra for two cross-stream locations at
    all Reynolds numbers.
    """
    u_list = [0.4, 0.6, 0.8, 1.0, 1.2]
    cm = plt.cm.coolwarm
    y_R_a = -1.0
    y_R_b = 1.5
    z_H = 0.25
    plt.figure(figsize=(7.5, 3.75))
    plt.subplot(1, 2, 1)
    label_subplot(text="(a)")
    for n, u in enumerate(u_list):
        plot_vel_spec(u, y_R_a, z_H, n_band_ave=n_band_ave, newfig=False,
                      plot_conf_int=plot_conf_int, plot_lines=(u==1.2),
                      color=cm(int(n/4*256)))
    plt.legend(loc="best")
    plt.subplot(1, 2, 2)
    label_subplot(text="(b)")
    for n, u in enumerate(u_list):
        plot_vel_spec(u, y_R_b, z_H, n_band_ave=n_band_ave, newfig=False,
                      plot_conf_int=plot_conf_int, plot_lines=(u==1.2),
                      color=cm(int(n/4*256)))
    if save:
        plt.savefig("Figures/wake_spectra" + savetype)
    if show:
        plt.show()


def plot_vertical_lines(xlist, ymaxscale=1, color="gray"):
    if not isinstance(xlist, list):
        x = [x]
    ymin = plt.axis()[2]
    ymax = plt.axis()[3]*ymaxscale
    for x in xlist:
        plt.vlines(x, ymin, ymax,
                   color=color, linestyles="dashed")
    plt.ylim((ymin, ymax))


def plot_wake_re_dep(y_R=0.0, z_H=0.25, save=False):
    """Plot the Reynolds number dependence of the streamwise mean velocity and
    turbulence intensity at y/R = 0, z/H = 0.25. Averages are taken from two
    runs for each speed.
    """
    mean_u = []
    std_u = []
    wake_speeds = np.arange(0.4, 1.4, 0.2)
    if y_R == 0.0 and z_H == 0.25:
        speeds = np.arange(0.3, 1.4, 0.1)
    else:
        speeds = wake_speeds
    for speed in speeds:
        if y_R == 0.0 and z_H == 0.25:
            s = Section("Perf-{:.1f}".format(speed))
            df = s.data
            df = df[np.round(df.mean_tsr, decimals=1)==1.9]
            df = df[df.y_R == y_R]
            df = df[df.z_H == z_H]
            mean_u.append(df.mean_u.mean())
            std_u.append(df.std_u.mean())
            if speed in wake_speeds:
                s = Section("Wake-{:.1f}".format(speed))
                df = s.data
                df = df[np.round(df.mean_tsr, decimals=1)==1.9]
                df = df[df.y_R == y_R]
                df = df[df.z_H == z_H]
                mean_u[-1] = (mean_u[-1] + df.mean_u.mean())/2
                std_u[-1] = (std_u[-1] + df.std_u.mean())/2
        else:
            s = Section("Wake-{:.1f}".format(speed))
            df = s.data
            df = df[np.round(df.mean_tsr, decimals=1)==1.9]
            df = df[df.y_R == y_R]
            df = df[df.z_H == z_H]
            mean_u.append(df.mean_u.mean())
            std_u.append(df.std_u.mean())
    mean_u = np.asarray(mean_u)
    std_u = np.asarray(std_u)
    plt.plot(speeds, mean_u/speeds, "-o")
#    plt.ylim((0.15, 0.35))
    plt.figure()
    plt.plot(speeds, std_u/speeds, "-o")
#    plt.ylim((0.05, 0.2))


def label_subplot(ax=None, x=0.5, y=-0.25, text="(a)", **kwargs):
    """Create a subplot label."""
    if ax is None:
        ax = plt.gca()
    ax.text(x=x, y=y, s=text, transform=ax.transAxes,
            horizontalalignment="center", verticalalignment="top", **kwargs)


def twiny_sci_label(ax=None, power=5, subplots=True):
    """Put scientific notation label on secondary x-axis."""
    if ax is None:
        ax = plt.gca()
    use_mathtext = plt.rcParams["axes.formatter.use_mathtext"]
    if use_mathtext:
        x, y = 0.90, 1.1
        if subplots:
            x, y = x*0.955, y*1.03
        text = r"$\times\mathregular{{10^{}}}$".format(power)
    else:
        x, y = 0.95, 1.08
        if subplots:
            x, y = x*0.955, y*1.03
        text = "1e{}".format(power)
    ax.text(x=x, y=y, s=text, transform=ax.transAxes)


def make_velocity_unc_table(save=False):
    """Create table with average velocity uncertainties."""
    tow_speeds = np.arange(0.4, 1.3, 0.2)
    mean_unc_u = []
    mean_unc_v = []
    mean_unc_w = []
    for u in tow_speeds:
        s = Section("Wake-{:.1f}".format(u))
        mean_unc_u.append(s.data.exp_unc_mean_u.mean())
        mean_unc_v.append(s.data.exp_unc_mean_v.mean())
        mean_unc_w.append(s.data.exp_unc_mean_w.mean())
    df = pd.DataFrame()
    df[r"$U_\infty$ (m/s)"] = tow_speeds
    df["$U$"] = mean_unc_u
    df["$V$"] = mean_unc_v
    df["$W$"] = mean_unc_w
    def speed_format(speed):
        return "${:.1f}$".format(speed)
    def unc_format(unc):
        un = "{:.0e}".format(unc).split("e")
        return r"${} \times 10^{{{}}}$".format(un[0], int(un[1]))
    fmt = [speed_format, unc_format, unc_format, unc_format]
    if save:
        if not os.path.isdir("Tables"):
            os.mkdir("Tables")
        df.to_latex(buf="Tables/mean_vel_unc.tex", index=False,
                    column_format="cccc", escape=False, formatters=fmt)
        df.to_csv("Tables/mean_vel_unc.csv", index=False)
    print("\nAverage wake velocity uncertainties (LaTeX formatted):\n")
    print(df.to_latex(index=False, column_format="cccc", escape=False,
                      formatters=fmt))

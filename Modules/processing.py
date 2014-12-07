# -*- coding: utf-8 -*-
"""
This code processes the data from the 2014 Spring RVAT Re dependence 
experiment.

@author: Pete Bachant

"""
from __future__ import division, print_function
import numpy as np
from pxl import timeseries as ts
import matplotlib.pyplot as plt
from scipy.io import loadmat
from scipy import interpolate
import multiprocessing as mp
import scipy.stats
from scipy.stats import nanmean, nanstd
from pxl import fdiff
import json
import os
import sys
import pandas as pd
try:
    import Modules.plotting
    from Modules.plotting import ylabels
except ImportError:
    import plotting
    from plotting import ylabels

# Dict for runs corresponding to each height
wakeruns = {0.0 : np.arange(0, 45),
            0.125 : np.arange(45, 90),
            0.25 : np.arange(90, 135),
            0.375 : np.arange(135, 180),
            0.5 : np.arange(180, 225),
            0.625 : np.arange(225, 270)}
           
# Constants
H = 1.0
D = 1.0
A = D*H
R = D/2
rho = 1000.0
nu = 1e-6

# Directory constants
raw_data_dir = os.path.join("Data", "Raw")
processed_data_dir = os.path.join("Data", "Processed")

def calc_b_vec(vel):
    """Calculates the systematic error of a Vectrino measurement (in m/s)
    from their published specs. Returns half the +/- value as b."""
    return 0.5*(0.005*np.abs(vel) + 0.001)
    
def calc_uncertainty(quantity, b):
    return np.sqrt(nanstd(quantity)**2 + b**2)
             
def calc_tare_torque(rpm):
    """Returns tare torque array given RPM array."""
    return 0.000474675989476*rpm + 0.876750155952
             
times = {0.3 : (20.0, 80.0),
         0.4 : (20.0, 60.0),
         0.5 : (20.0, 50.0),
         0.6 : (20.0, 45.0),
         0.7 : (20.0, 38.0),
         0.8 : (18.0, 34.0),
         0.9 : (16.0, 32.0),
         1.0 : (15.0, 30.0),
         1.1 : (15.0, 28.0),
         1.2 : (14.0, 25.0),
         1.3 : (13.0, 23.0),
         1.4 : (12.0, 20.0)}

if "linux" in sys.platform:
    cfd_path = "/media/pete/BigPocket/OpenFOAM/pete-2.3.0/run/unh-rvat-2d_re-dep_2"
elif "win" in sys.platform:
    cfd_path = "G:/OpenFOAM/pete-2.3.0/run/unh-rvat-2d_re-dep_2"

class Run(object):
    """Object that represents a single turbine tow"""
    def __init__(self, section, nrun):
        self.section = section
        section_raw_dir = os.path.join("Data", "Raw", section)
        if nrun < 0:
            runs = []
            for f in os.listdir(section_raw_dir):
                try: 
                    runs.append(int(f))
                except ValueError:
                    pass
            self.nrun = sorted(runs)[nrun]
        else:
            self.nrun = nrun
        self.raw_dir = os.path.join(section_raw_dir, str(self.nrun))
        self.loaded = False
        self.t2found = False
        self.not_loadable = False
        self.wake_calculated = False
        # Do all processing
        self.load()
        self.subtract_tare_drag()
        self.add_tare_torque()
        self.calc_perf_instantaneous()
        self.make_trimmed()
        self.filter_wake()
        self.calc_wake_instantaneous()
        self.calc_perf_per_rev()
        self.calc_perf_stats()
        self.calc_wake_stats()
        self.calc_perf_uncertainty()
        self.calc_perf_exp_uncertainty()
        
    def load(self):
        """Loads the data from the run into memory"""
        self.load_metadata()
        self.load_nidata()
        self.load_acsdata()
        self.load_vecdata()
        self.loaded = True
        
    def load_metadata(self):
        """Loads run metadata."""
        try: 
            with open(os.path.join(self.raw_dir, "metadata.json")) as f:
                self.metadata = json.load(f)
        except IOError:
            if self.nrun < 0:
                self.nrun -= 1
                self.raw_dir = os.path.join("Data", "Raw", section, str(self.nrun))
            else:
                pass
        try:
            with open(os.path.join(self.raw_dir, "metadata.json")) as f:
                self.metadata = json.load(f)
        except IOError:
            self.not_loadable = True
            return None
        self.tow_speed_nom = np.round(self.metadata["Tow speed (m/s)"], decimals=1)
        self.y_R = self.metadata["Vectrino y/R"]
        self.z_H = self.metadata["Vectrino z/H"]
        
    def load_nidata(self):
        nidata = loadmat(os.path.join(self.raw_dir, "nidata.mat"), squeeze_me=True)
        self.time_ni = nidata["t"]
        self.sr_ni = (1.0/(self.time_ni[1] - self.time_ni[0]))
        if "carriage_pos" in nidata:
            self.lin_enc = True
            self.carriage_pos = nidata["carriage_pos"]
            self.tow_speed_ni = fdiff.second_order_diff(self.carriage_pos, self.time_ni)
            self.tow_speed_ni = ts.smooth(self.tow_speed_ni, 8)
            self.tow_speed_ref = self.tow_speed_ni
        else:
            self.lin_enc = False
            self.tow_speed_ref = self.tow_speed_nom
        self.torque = nidata["torque_trans"]
        self.torque_arm = nidata["torque_arm"]
        self.drag = nidata["drag_left"] + nidata["drag_right"]
        # Remove offsets from drag, not torque
        t0 = 2
        self.drag = self.drag - np.mean(self.drag[0:self.sr_ni*t0])
        # Compute RPM and omega
        self.angle = nidata["turbine_angle"]
        self.rpm_ni = fdiff.second_order_diff(self.angle, self.time_ni)/6.0
        self.rpm_ni = ts.smooth(self.rpm_ni, 8)
        self.omega_ni = self.rpm_ni*2*np.pi/60.0
        self.omega = self.omega_ni
        self.tow_speed = self.tow_speed_ref
        
    def load_acsdata(self):
        fpath = os.path.join(self.raw_dir, "acsdata.mat")
        acsdata = loadmat(fpath, squeeze_me=True)
        self.tow_speed_acs = acsdata["carriage_vel"]
        self.rpm_acs = acsdata["turbine_rpm"]
        self.rpm_acs = ts.sigmafilter(self.rpm_acs, 3, 3)
        self.omega_acs = self.rpm_acs*2*np.pi/60.0
        self.time_acs = acsdata["t"]
        if len(self.time_acs) != len(self.omega_acs):
            newlen = np.min((len(self.time_acs), len(self.omega_acs)))
            self.time_acs = self.time_acs[:newlen]
            self.omega_acs = self.omega_acs[:newlen]
        self.omega_acs_interp = np.interp(self.time_ni, self.time_acs, self.omega_acs)
        self.rpm_acs_interp = self.omega_acs_interp*60.0/(2*np.pi)
        
    def load_vecdata(self):
        try:
            vecdata = loadmat(self.raw_dir + "/" + "vecdata.mat", 
                              squeeze_me=True)
            self.sr_vec = 200.0
            self.time_vec = vecdata["t"]
            self.u = vecdata["u"]
            self.v = vecdata["v"]
            self.w = vecdata["w"]
        except IOError:
            self.vecdata = None
        
    def subtract_tare_drag(self):
        df = pd.read_csv(os.path.join("Data", "Processed", "Tare drag.csv"))
        self.tare_drag = df.tare_drag[df.tow_speed==self.tow_speed_nom].values[0]
        self.drag = self.drag - self.tare_drag
        
    def add_tare_torque(self):
        # Choose reference RPM, using NI for all except Perf-0.4
        if self.section == "Perf-0.4":
            rpm_ref = self.rpm_acs_interp
        else:
            rpm_ref = self.rpm_ni
        # Add tare torque
        self.tare_torque = calc_tare_torque(rpm_ref)
        self.torque += self.tare_torque
        
    def calc_perf_instantaneous(self):
        if self.section == "Perf-0.4":
            omega_ref = self.omega_acs_interp
        else:
            omega_ref = self.omega_ni
        # Compute power
        self.power = self.torque*omega_ref
        self.tsr = omega_ref*R/self.tow_speed_ref
        # Compute power, drag, and torque coefficients
        self.cp = self.power/(0.5*rho*A*self.tow_speed_ref**3)
        self.cd = self.drag/(0.5*rho*A*self.tow_speed_ref**2)
        self.ct = self.torque/(0.5*rho*A*R*self.tow_speed_ref**2)
        # Remove datapoints for coefficients where tow speed is small
        self.cp[np.abs(self.tow_speed_ref < 0.01)] = np.nan
        self.cd[np.abs(self.tow_speed_ref < 0.01)] = np.nan
        
    def load_vectxt(self):
        """Loads Vectrino data from text (*.dat) file."""
        data = np.loadtxt(self.raw_dir + "/vecdata.dat", unpack=True)
        self.time_vec_txt = data[0]
        self.u_txt = data[3]
        
    def make_trimmed(self):
        """Trim all time series and replace the full run names with names with
        the '_all' suffix."""
        # Put in some guesses for t1 and t2
        self.t1, self.t2 = times[self.tow_speed_nom]
        self.find_t2()
        # Trim performance quantities
        self.time_ni_all = self.time_ni
        self.time_perf_all = self.time_ni
        self.time_ni = self.time_ni_all[self.t1*self.sr_ni:self.t2*self.sr_ni]
        self.time_perf = self.time_ni
        self.angle_all = self.angle
        self.angle = self.angle_all[self.t1*self.sr_ni:self.t2*self.sr_ni]
        self.torque_all = self.torque
        self.torque = self.torque_all[self.t1*self.sr_ni:self.t2*self.sr_ni]
        self.torque_arm_all = self.torque_arm
        self.torque_arm = self.torque_arm_all[self.t1*self.sr_ni:self.t2*self.sr_ni]
        self.omega_all = self.omega
        self.omega = self.omega_all[self.t1*self.sr_ni:self.t2*self.sr_ni]
        self.tow_speed_all = self.tow_speed
        if self.lin_enc:
            self.tow_speed = self.tow_speed_all[self.t1*self.sr_ni:self.t2*self.sr_ni]
        self.tsr_all = self.tsr
        self.tsr = self.tsr_all[self.t1*self.sr_ni:self.t2*self.sr_ni]
        self.cp_all = self.cp
        self.cp = self.cp_all[self.t1*self.sr_ni:self.t2*self.sr_ni]
        self.ct_all = self.ct
        self.ct = self.ct_all[self.t1*self.sr_ni:self.t2*self.sr_ni]
        self.cd_all = self.cd
        self.cd = self.cd_all[self.t1*self.sr_ni:self.t2*self.sr_ni]
        self.rpm_ni_all = self.rpm_ni
        self.rpm_ni = self.rpm_ni_all[self.t1*self.sr_ni:self.t2*self.sr_ni]
        self.rpm = self.rpm_ni
        self.rpm_all = self.rpm_ni_all
        self.drag_all = self.drag
        self.drag = self.drag_all[self.t1*self.sr_ni:self.t2*self.sr_ni]
        # Trim wake quantities
        self.u_all = self.u
        self.u = self.u_all[self.t1*self.sr_vec:self.t2*self.sr_vec]
        self.v_all = self.v
        self.v = self.v_all[self.t1*self.sr_vec:self.t2*self.sr_vec]
        self.w_all = self.w
        self.w = self.w_all[self.t1*self.sr_vec:self.t2*self.sr_vec]
        
    def find_t2(self):
        sr = self.sr_ni
        angle1 = self.angle[sr*self.t1]
        angle2 = self.angle[sr*self.t2]
        n3rdrevs = np.floor((angle2-angle1)/120.0)
        self.n_revs = int(np.floor((angle2-angle1)/360.0))
        self.n_blade_pass = int(n3rdrevs)
        angle2 = angle1 + n3rdrevs*120
        t2i = np.where(np.round(self.angle)==np.round(angle2))[0][0]
        t2 = self.time_ni[t2i]
        self.t2 = np.round(t2, decimals=2)
        self.t2found = True
        self.t1_wake = self.t1
        self.t2_wake = self.t2
        
    def calc_perf_stats(self):
        """Calculates mean performance based on trimmed time series."""
        self.mean_tsr, self.std_tsr = nanmean(self.tsr), nanstd(self.tsr)
        self.mean_cp, self.std_cp = nanmean(self.cp), nanstd(self.cp)
        self.mean_cd, self.std_cd = nanmean(self.cd), nanstd(self.cd)
        self.mean_ct, self.std_ct = nanmean(self.ct), nanstd(self.ct)
        if self.lin_enc:
            self.mean_u_enc = nanmean(self.tow_speed)
            self.std_u_enc = nanstd(self.tow_speed)
        else:
            self.mean_u_enc = np.nan
            self.std_u_enc = np.nan

    def print_perf_stats(self):
        print("tow_speed_nom =", self.tow_speed_nom)
        if self.lin_enc:
            print("tow_speed_enc =", self.mean_u_enc, "std =", self.std_u_enc)
        print("TSR = {:.2f} +/- {:.2f}".format(self.mean_tsr, self.exp_unc_tsr))
        print("C_P = {:.2f} +/- {:.2f}".format(self.mean_cp, self.exp_unc_cp))
        print("C_D = {:.2f} +/- {:.2f}".format(self.mean_cd, self.exp_unc_cd))
        
    def calc_perf_uncertainty(self):
        """See uncertainty IPython notebook for equations."""
        # Systematic uncertainty estimates
        b_torque = 0.5/2
        b_angle = 3.14e-5/2
        b_car_pos = 0.5e-5/2
        b_force = 0.28/2
        # Uncertainty of C_P
        omega = self.omega.mean()
        torque = self.torque.mean()
        u_infty = np.mean(self.tow_speed)
        const = 0.5*rho*A
        b_cp = np.sqrt((omega/(const*u_infty**3))**2*b_torque**2 + \
                       (torque/(const*u_infty**3))**2*b_angle**2 + \
                       (-3*torque*omega/(const*u_infty**4))**2*b_car_pos**2)
        self.b_cp = b_cp
        self.unc_cp = calc_uncertainty(self.cp_per_rev, b_cp)
        # Drag coefficient
        drag = self.drag.mean()
        b_cd = np.sqrt((1/(const*u_infty**2))**2*b_force**2 + \
                       (1/(const*u_infty**2))**2*b_force**2 +
                       (-2*drag/(const*u_infty**3))**2*b_car_pos**2)
        self.unc_cd = calc_uncertainty(self.cd_per_rev, b_cd)
        self.b_cd = b_cd
        # Tip speed ratio
        b_tsr = np.sqrt((R/(u_infty))**2*b_angle**2 + \
                        (-omega*R/(u_infty**2))**2*b_car_pos**2)
        self.unc_tsr = calc_uncertainty(self.tsr_per_rev, b_tsr)
        self.b_tsr = b_tsr
        
    def calc_perf_exp_uncertainty(self):
        """See uncertainty IPython notebook for equations."""
        # Power coefficient
        s_cp = self.std_cp_per_rev
        nu_s_cp = len(self.cp_per_rev) - 1
        b_cp = self.b_cp
        b_cp_rel_unc = 0.25 # A guess
        nu_b_cp = 0.5*b_cp_rel_unc**(-2)
        nu_cp = ((s_cp**2 + b_cp**2)**2)/(s_cp**4/nu_s_cp + b_cp**4/nu_b_cp)
        t = scipy.stats.t.interval(alpha=0.95, df=nu_cp)[-1]
        self.exp_unc_cp = t*self.unc_cp
        self.dof_cp = nu_cp
        # Drag coefficient
        s_cd = self.std_cd_per_rev
        nu_s_cd = len(self.cd_per_rev) - 1
        b_cd = self.b_cd
        b_cd_rel_unc = 0.25 # A guess
        nu_b_cd = 0.5*b_cd_rel_unc**(-2)
        nu_cd = ((s_cd**2 + b_cd**2)**2)/(s_cd**4/nu_s_cd + b_cd**4/nu_b_cd)
        t = scipy.stats.t.interval(alpha=0.95, df=nu_cd)[-1]
        self.exp_unc_cd = t*self.unc_cd
        self.dof_cd = nu_cd
        # Tip speed ratio
        s_tsr = self.std_tsr_per_rev
        nu_s_tsr = len(self.tsr_per_rev) - 1
        b_tsr = self.b_tsr
        b_tsr_rel_unc = 0.25 # A guess
        nu_b_tsr = 0.5*b_tsr_rel_unc**(-2)
        nu_tsr = ((s_tsr**2 + b_tsr**2)**2)/(s_tsr**4/nu_s_tsr + b_tsr**4/nu_b_tsr)
        t = scipy.stats.t.interval(alpha=0.95, df=nu_tsr)[-1]
        self.exp_unc_tsr = t*self.unc_tsr
        self.dof_tsr = nu_tsr
        
    def calc_wake_instantaneous(self):
        """Creates fluctuating and Reynolds stress time series. Note that
        time series must be trimmed first, or else subtracting the mean makes
        no sense. Prime variables are denoted by a `p` e.g., $u'$ is `up`."""
        self.up = self.u - nanmean(self.u)
        self.vp = self.v - nanmean(self.v)
        self.wp = self.w - nanmean(self.w)
        self.upup = self.up**2
        self.upvp = self.up*self.vp
        self.upwp = self.up*self.wp
        self.vpvp = self.vp**2
        self.vpwp = self.vp*self.wp
        self.wpwp = self.wp**2
        
    def filter_wake(self, stdfilt=True, threshfilt=True):
        """Applies filtering to wake velocity data with a standard deviation
        filter, threshold filter, or both. Renames unfiltered time series with
        the '_unf' suffix. Time series are already trimmed before they reach 
        this point, so no slicing is necessary"""
        std = 8
        passes = 1
        fthresh = 0.9
        # Calculate means
        mean_u = self.u.mean()
        mean_v = self.v.mean()
        mean_w = self.w.mean()
        # Create new unfiltered arrays
        self.u_unf = self.u*1
        self.v_unf = self.v*1
        self.w_unf = self.w*1
        if stdfilt:
        # Do standard deviation filters
            self.u = ts.sigmafilter(self.u, std, passes)
            self.v = ts.sigmafilter(self.v, std, passes)
            self.w = ts.sigmafilter(self.w, std, passes)
        if threshfilt:
            # Do threshold filter on u
            ibad = np.where(self.u > mean_u + fthresh)[0]
            ibad = np.append(ibad, np.where(self.u < mean_u - fthresh)[0])
            i = np.where(np.logical_and(ibad > self.t1*200, 
                                        ibad < self.t2*200))[0]
            self.u[ibad[i]] = np.nan
            # Do threshold filter on v
            ibad = np.where(self.v > mean_v + fthresh)[0]
            ibad = np.append(ibad, np.where(self.v < mean_v - fthresh)[0])
            i = np.where(np.logical_and(ibad > self.t1*200, 
                                        ibad < self.t2*200))[0]
            self.v[ibad[i]] = np.nan
            # Do threshold filter on w
            ibad = np.where(self.w > mean_w + fthresh)[0]
            ibad = np.append(ibad, np.where(self.w < mean_w - fthresh)[0])
            i = np.where(np.logical_and(ibad > self.t1*200, 
                                        ibad < self.t2*200))[0]
            self.w[ibad[i]] = np.nan
        # Count up bad datapoints
        self.nbadu = len(np.where(np.isnan(self.u)==True)[0])
        self.nbadv = len(np.where(np.isnan(self.v)==True)[0])
        self.nbadw = len(np.where(np.isnan(self.w)==True)[0])
        self.nbad = self.nbadu + self.nbadv + self.nbadw
        
    def calc_wake_stats(self):
        if self.not_loadable:
            self.mean_u = np.nan
            self.mean_v = np.nan
            self.mean_w = np.nan
            return None
        if not self.t2found:
            self.find_t2()
        self.filter_wake()
        self.mean_u, self.std_u = nanmean(self.u), nanstd(self.u)
        self.mean_v, self.std_v = nanmean(self.v), nanstd(self.v)
        self.mean_w, self.std_w = nanmean(self.w), nanstd(self.w)
        self.mean_upup, self.std_upup = nanmean(self.upup), nanstd(self.upup)
        self.mean_upvp, self.std_upvp = nanmean(self.upvp), nanstd(self.upvp)
        self.mean_upwp, self.std_upwp = nanmean(self.upwp), nanstd(self.upwp)
        self.mean_vpvp, self.std_vpvp = nanmean(self.vpvp), nanstd(self.vpvp)
        self.mean_vpwp, self.std_vpwp = nanmean(self.vpwp), nanstd(self.vpwp)
        self.mean_wpwp, self.std_wpwp = nanmean(self.wpwp), nanstd(self.wpwp)
        self.k = 0.5*(self.mean_upup + self.mean_vpvp + self.mean_wpwp)
            
    def print_wake_stats(self):
        ntotal = int((self.t2 - self.t1)*self.sr_vec*3)
        print("y/R =", self.y_R)
        print("z/H =", self.z_H)
        print("tow_speed_vec/tow_speed_nom =", self.mean_u/self.tow_speed_nom, "+/-", 
              self.delta_mean_u/2/self.tow_speed_nom)
        print("std_u/tow_speed_nom =", self.std_u/self.tow_speed_nom, "+/-",
              self.delta_std_u/2/self.tow_speed_nom)
        print(str(self.nbad)+"/"+str(ntotal), "data points omitted")
        
    def calc_wake_uncertainty(self):
        """Computes delta values for wake measurements from Vectrino accuracy
        specs, not statistical uncertainties."""
        self.unc_mean_u = np.nan
        self.unc_std_u = np.nan
        
    def calc_perf_per_rev(self):
        """Computes mean power coefficient over each revolution."""
        angle = self.angle*1
        angle -= angle[0]
        cp = np.zeros(self.n_revs)
        cd = np.zeros(self.n_revs)
        tsr = np.zeros(self.n_revs)
        torque = np.zeros(self.n_revs)
        omega = np.zeros(self.n_revs)
        u_infty3 = np.zeros(self.n_revs)
        start_angle = 0.0
        for n in range(self.n_revs):
            end_angle = start_angle + 360
            ind = np.logical_and(angle >= start_angle, end_angle > angle)
            cp[n] = self.cp[ind].mean()
            cd[n] = self.cd[ind].mean()
            tsr[n] = self.tsr[ind].mean()
            torque[n] = self.torque[ind].mean()
            omega[n] = self.omega[ind].mean()
            try:
                u_infty3[n] = (self.tow_speed_ni**3)[ind].mean()
            except AttributeError:
                u_infty3[n] = self.tow_speed_nom**3
            start_angle += 360
        self.cp_per_rev = cp
        self.std_cp_per_rev = cp.std()
        self.cd_per_rev = cd
        self.std_cd_per_rev = cd.std()
        self.tsr_per_rev = tsr
        self.std_tsr_per_rev = tsr.std()
        self.torque_per_rev = torque
        self.std_torque_per_rev = torque.std()
        self.u_infty3_per_rev = u_infty3
        
    @property
    def cp_conf_interval(self, alpha=0.95):
        self.calc_perf_per_rev()
        t_val = scipy.stats.t.interval(alpha=alpha, df=self.n_revs-1)[1]
        std = self.std_cp_per_rev
        return t_val*std/np.sqrt(self.n_revs)
        
    def detect_badvec(self):
        """Detects if Vectrino data is bad by looking at first 2 seconds of
        data, and checking if there are many datapoints."""
        nbad = len(np.where(np.abs(self.u[:400]) > 0.5)[0])
        print(nbad, "bad Vectrino datapoints in first 2 seconds")
        if nbad > 50:
            self.badvec = True
            print("Vectrino data bad")
        else:
            self.badvec = False
            print("Vectrino data okay")
        
    def plot_perf(self, quantity="power coefficient"):
        """Plots the run's data"""
        if not self.loaded:
            self.load()
        if quantity == "drag":
            quantity = self.drag
            ylabel = "Drag (N)"
            ylim = None
        elif quantity == "torque":
            quantity = self.torque
            ylabel = "Torque (Nm)"
            ylim = None
        elif quantity.lower == "power coefficient" or "cp" or "c_p":
            quantity = self.cp
            ylabel = "$C_P$"
            ylim = (-1, 1)
        plt.figure()
        plt.plot(self.time_ni, quantity, 'k')
        plt.xlabel("Time (s)")
        plt.ylabel(ylabel)
        plt.ylim(ylim)
        plt.tight_layout()
        
    def plot_wake(self):
        """Plot streamwise velocity over experiment."""
        if not self.loaded:
            self.load()
        plt.figure()
        self.filter_wake()
        plt.plot(self.time_vec, self.u, 'k')
        plt.xlabel("Time (s)")
        plt.ylabel("$u$ (m/s)")
        
    def plot_acs(self):
        if not self.loaded:
            self.load()
        plt.figure()
        plt.plot(self.time_acs, self.rpm_acs)
        plt.hold(True)
        plt.plot(self.time_ni, self.rpm_ni)
        plt.figure()
        plt.plot(self.time_ni, self.tow_speed_ni)
        plt.hold(True)
        plt.plot(self.time_acs, self.tow_speed_acs)
        plt.show()
        
    def plot_carriage_vel(self):
        if not self.loaded:
            self.load()
        plt.figure()
        plt.plot(self.time_ni, self.tow_speed_ni)
        plt.tight_layout()
        plt.show()

        
class Section(object):
    def __init__(self, name):
        self.name = name
        self.processed_path = os.path.join(processed_data_dir, name+".csv")
        self.test_plan_path = os.path.join("Config", "Test plan", name+".csv")
        self.load()    
    def load(self):
        self.data = pd.read_csv(self.processed_path)
        self.test_plan = pd.read_csv(self.test_plan_path)
    @property
    def cp(self):
        return self.data.cp
    def process(self):
        """To-do: Process an entire section of data."""
        pass
    def process_parallel(self, nproc=8, nruns=24):
        """To-do: Process a section in parallel."""
        output = mp.Queue()
        processes = []
        runs_per_proc = int(np.ceil(nruns/nproc))
        for n in range(nproc):
            n1 = n*runs_per_proc
            n2 = (n+1)*runs_per_proc
            if n2 > nruns:
                n2 = nruns + 1
            runs = list(range(n1, n2))
            processes.append(mp.Process(target=process_runs_parallel, 
                                        args=(self.name, runs, output)))
        for p in processes:
            p.start()
        for p in processes:
            p.join()
        results = [output.get() for p in processes]
        results.sort()
        all_results = []
        for r in results:
            all_results += r[-1]
        print(all_results)
            

def process_runs_parallel(section, nrunlist, output):
    mean_cp = []
    for nrun in nrunlist:
        run = Run(section, nrun)
        mean_cp.append(run.mean_cp)
    output.put((nrunlist[0], mean_cp))


class PerfCurve(object):
    """Object that represents a performance curve."""
    def __init__(self, tow_speed):
        self.tow_speed = tow_speed
        self.Re_D = tow_speed*D/nu
        self.section = "Perf-{}".format(tow_speed)
        self.raw_data_dir = os.path.join("Data", "Raw", self.section)
        self.df = pd.read_csv(os.path.join("Data", "Processed", self.section+".csv"))
        self.testplan = pd.read_csv(os.path.join("Config", "Test plan", self.section+".csv")) 
        
    def process(self, reprocess=True):
        """Calculates power and drag coefficients for each run"""
        print("Processing", self.section)
        if not reprocess:
            print("Leaving processed runs as they are")
            try:
                runsdone = np.load(self.folder+"/Processed/runs.npy")
                tsr_old = np.load(self.folder+"/Processed/tsr.npy")
                cp_old = np.load(self.folder+"/Processed/cp.npy")
                cd_old = np.load(self.folder+"/Processed/cp.npy")
            except IOError:
                reprocess = False
                runsdone = []
        tsr = np.zeros(len(self.runs))
        cp = np.zeros(len(self.runs))
        cd = np.zeros(len(self.runs))
        for n in range(len(self.runs)):
            nrun = self.runs[n]
            if reprocess or nrun not in runsdone or np.isnan(tsr_old[n]):
                print("Processing run " + str(self.runs[nrun]) + "...")
                run = Run("Perf-{:0.1f}".format(self.tow_speed), nrun)
                run.calc_perf()
                tsr[n] = run.meantsr
                cp[n] = run.meancp
                cd[n] = run.meancd
            else:
                tsr[n] = tsr_old[np.where(runsdone==nrun)[0]]
                cp[n] = cp_old[np.where(runsdone==nrun)[0]]
                cd[n] = cd_old[np.where(runsdone==nrun)[0]]
        # Save updated DataFrame
        
    def plotcp(self, newfig=True, show=True, save=False, savepath="",
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
            plt.savefig(os.path.join(savepath, "cp_vs_tsr" + savetype))
        if show:
            plt.show()
            
    def plotcd(self, newfig=True, show=True, save=False, savepath="",
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
            plt.savefig(os.path.join(savepath, "cd_vs_tsr" + savetype))
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
             savepath="", savetype=".pdf", linetype='--ok'):
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
            plt.grid()
        plt.plot(y_R, q, "-.^k", label=r"$Re_D=0.4 \times 10^6$")
#        plot_old_wake(quantity, y_R)
        plt.legend(loc=loc)
        plt.tight_layout()
        if show:
            plt.show()
        if save:
            plt.savefig(savepath+quantity+"_Re_dep_exp"+savetype)
    
class WakeMap(object):
    def __init__(self, U_infty):
        self.U_infty = U_infty
        self.z_H = [0.0, 0.125, 0.25, 0.375, 0.5, 0.625]
        self.loaded = False
        self.load()
        
    def load(self):
        self.y_R = WakeProfile(self.U_infty, 0, "mean_u").y_R
        self.mean_u = np.zeros((len(self.z_H), len(self.y_R)))
        self.mean_v = self.mean_u*1
        self.mean_w = self.mean_u*1
        for z_H in self.z_H:
            wp = WakeProfile(self.U_infty, z_H, "mean_u")
            self.mean_u[self.z_H.index(z_H)] = wp.df.mean_u
            self.mean_v[self.z_H.index(z_H)] = wp.df.mean_v
            self.mean_w[self.z_H.index(z_H)] = wp.df.mean_w
        self.loaded = True
        
    def turb_lines(self, linestyles="solid", linewidth=3, colors="gray"):
        plt.hlines(0.5, -1, 1, linestyles=linestyles, colors="gray",
                   linewidth=linewidth)
        plt.vlines(-1, -0.2, 0.5, linestyles=linestyles, colors="gray",
                   linewidth=linewidth)
        plt.vlines(1, -0.2, 0.5, linestyles=linestyles, colors="gray",
                   linewidth=linewidth)
        
    def plot_mean_u(self, save=False, show=False, savepath="", savetype=".pdf"):
        # Plot contours of mean streamwise velocity
        plt.figure(figsize=(10,5))
        cs = plt.contourf(self.y_R, self.z_H, self.mean_u, 20,
                          cmap=plt.cm.coolwarm)
        plt.xlabel(r"$y/R$")
        plt.ylabel(r"$z/H$")
        cb = plt.colorbar(cs, shrink=1, extend="both", 
                          orientation="horizontal", pad=0.3)
        cb.set_label(r"$U/U_{\infty}$")
        self.turb_lines()
        ax = plt.axes()
        ax.set_aspect(2)
        plt.yticks([0,0.13,0.25,0.38,0.5,0.63])
        plt.tight_layout()
        if save:
            plt.savefig(savepath+"/mean_u_cont"+savetype)
        if show:
            self.show()
    
    def plot_meancomboquiv(self, save=False, show=False, savepath="",
                           savetype=".pdf"):
        plt.figure(figsize=(10,6))
        # Add contours of mean velocity
        cs = plt.contourf(self.y_R, self.z_H, self.mean_u/self.U_infty, 20, 
                          cmap=plt.cm.coolwarm)
        cb = plt.colorbar(cs, shrink=1, extend="both", 
                          orientation="horizontal", pad=0.2)
        cb.set_label(r"$U/U_{\infty}$")
        plt.hold(True)
        # Make quiver plot of v and w velocities
        Q = plt.quiver(self.y_R, self.z_H, self.mean_v/self.U_infty, 
                       self.mean_w/self.U_infty, width=0.0022)
        plt.xlabel(r"$y/R$")
        plt.ylabel(r"$z/H$")
        plt.ylim(-0.2, 0.78)
        plt.xlim(-3.2, 3.2)
        plt.quiverkey(Q, 0.75, 0.3, 0.1, r"$0.1 U_\infty$",
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
            plt.savefig(savepath+"/meancomboquiv"+savetype)
    
    def plot_xvorticity(self):
        pass
    
    def plot_diff(self, quantity="mean_u", U_infty_diff=1.0, save=False, 
                  show=False, savepath="", savetype=""):
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
            if savepath: savepath += "/"
            plt.savefig(savepath+"/"+quantity+"_diff"+savetype)
    
    def plot_meancomboquiv_diff(self, U_infty_diff, save=False, show=False,
                                savepath="", savetype="", percent=True):
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
            if savepath: savepath += "/"
            plt.savefig(savepath+"/meancomboquiv_diff"+savetype)
            
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
        
def load_test_plan_section(section):
    df = pd.read_csv(os.path.join("Config", "Test plan", section+".csv"))
    df = df.dropna(how="all", axis=1).dropna(how="all", axis=0)
    if "Run" in df:
        df["Run"] = df["Run"].astype(int)
    return df
            
def batch_process_section(section, reprocess=True):
    """Processes all data in a section. Will skip already processed
    runs if `reprocess = False`. Something is up with this algorithm, as it
    sometimes will save the wrong y_R value."""
    test_plan = load_test_plan_section(section)
    if not reprocess:
        try:
            data = pd.read_csv(os.path.join("Data", "Processed", section+".csv"))
            runs = data.run
        except:
            reprocess = True
            data = pd.DataFrame()
            data["run"] = runs = test_plan["Run"]
    else:
        print("Reprocessing", section)
        data = pd.DataFrame()
        data["run"] = runs = test_plan["Run"]
    # Create a empty arrays for all quantities
    if reprocess:
        for q in ["mean_tow_speed", "std_tow_speed", "t1", "t2", "n_blade_pass", 
                  "n_revs", "mean_tsr", "mean_cp", "mean_cd", "std_tsr", 
                  "std_cp", "std_cd", "std_tsr_per_rev", "std_cp_per_rev",
                  "std_cd_per_rev", "exp_unc_tsr", "exp_unc_cp", "exp_unc_cd", 
                  "dof_tsr", "dof_cp", "dof_cd", "t1_wake", "t2_wake", "y_R", 
                  "z_H", "mean_u", "mean_v", "mean_w", "std_u", "std_v", 
                  "std_w", "mean_upvp", "mean_upwp", "mean_vpwp", "k"]:
            data[q] = np.zeros(len(data.run))*np.nan
    for n in runs:
        if reprocess or True in np.isnan(data.iloc[n,:]).values:
            r = Run(section, n)
            if r.not_loadable:
                pass
            else:
                print("Processing data from {} run {}".format(section, n))
                data.mean_tow_speed.iloc[n] = r.mean_u_enc
                data.std_tow_speed.iloc[n] = r.std_u_enc
                data.t1.iloc[n] = r.t1
                data.t1_wake.iloc[n] = r.t1_wake
                data.t2.iloc[n] = r.t2
                data.t2_wake.iloc[n] = r.t2_wake
                data.n_revs.iloc[n] = r.n_revs
                data.n_blade_pass.iloc[n] = r.n_blade_pass
                data.y_R.iloc[n] = r.y_R
                data.z_H.iloc[n] = r.z_H
                data.mean_tsr.iloc[n] = r.mean_tsr
                data.mean_cp.iloc[n] = r.mean_cp
                data.mean_cd.iloc[n] = r.mean_cd
                data.std_tsr.iloc[n] = r.std_tsr
                data.std_cp.iloc[n] = r.std_cp
                data.std_cd.iloc[n] = r.std_cd
                data.std_tsr_per_rev.iloc[n] = r.std_tsr_per_rev
                data.std_cp_per_rev.iloc[n] = r.std_cp_per_rev
                data.std_cd_per_rev.iloc[n] = r.std_cd_per_rev
                data.exp_unc_tsr.iloc[n] = r.exp_unc_tsr
                data.exp_unc_cp.iloc[n] = r.exp_unc_cp
                data.exp_unc_cd.iloc[n] = r.exp_unc_cd
                data.dof_tsr.iloc[n] = r.dof_tsr
                data.dof_cp.iloc[n] = r.dof_cp
                data.dof_cd.iloc[n] = r.dof_cd
                data.mean_u.iloc[n] = r.mean_u
                data.mean_v.iloc[n] = r.mean_v
                data.mean_w.iloc[n] = r.mean_w
                data.mean_upvp.iloc[n] = r.mean_upvp
                data.mean_upwp.iloc[n] = r.mean_upwp
                data.mean_vpwp.iloc[n] = r.mean_vpwp
                data.std_u.iloc[n] = r.std_u
                data.std_v.iloc[n] = r.std_v
                data.std_w.iloc[n] = r.std_w
                data.k.iloc[n] = r.k
    data.to_csv(os.path.join(processed_data_dir, section+".csv"), index=False,
                na_rep="NaN")
    
def batch_process_all():
    """Batch processes all sections."""
    sections = ["Perf-0.3", "Perf-0.4", "Perf-0.5",
                "Perf-0.6", "Perf-0.7", "Perf-0.8",
                "Perf-0.9", "Perf-1.0", "Perf-1.1",
                "Perf-1.2", "Perf-1.3", "Wake-0.4",
                "Wake-0.6", "Wake-0.8", "Wake-1.0",
                "Wake-1.2"]
    for section in sections:
        batch_process_section(section)
    
def plot_trans_wake_profile(quantity, U_infty=0.4, z_H=0.0, save=False, savepath="", 
                            savetype=".pdf", newfig=True, marker="-ok",
                            fill="none", oldwake=False, figsize=(10, 5)):
    """Plots the transverse wake profile of some quantity. These can be
      * mean_u
      * mean_v
      * mean_w
      * std_u
    """
    Re_D = U_infty*D/nu
    label = str((Re_D/1e6)) + "e6"
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
    plt.grid(True)
    plt.tight_layout()
    
def plot_perf_re_dep(save=False, savepath="", savetype=".pdf", errorbars=False,
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
        plt.savefig(savepath + "/re_dep_cp" + savetype)
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
        plt.savefig(savepath + "/re_dep_cd" + savetype)
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
    data = np.loadtxt("Settling/tow_speed_" + str(tow_speed) + "/vecdata.dat", unpack=True)
    u = data[2] # 2 for x velocity
    t = data[0]*0.005
#    i = np.where(np.round(t)==88)[0][0]
    i = 0
    u = u[i:]
    t = t[i:]
#    u_f = ts.sigmafilter(u, 3, 1)
    t_std, u_std = ts.runningstd(t, u, 400)
#    u_std = ts.smooth(u_std, 500)
    u = ts.smooth(u, 400)
    plt.figure()
    plt.plot(t-t[0], u, "k")
    plt.xlabel("t (s)")
    plt.ylabel("$u$ (m/s)")
    plt.tight_layout()
    plt.figure()
    plt.plot(t_std, u_std)
    plt.xlabel("t (s)")
    plt.ylabel(r"$\sigma_u$")
    plt.tight_layout()
    
def process_tare_drag(nrun, plot=False):
    """Processes a single tare drag run."""
    print("Processing tare drag run", str(nrun)+"...")
    times = {0.3 : (10, 77),
             0.4 : (8, 60),
             0.5 : (8, 47),
             0.6 : (10, 38),
             0.7 : (8, 33),
             0.8 : (7, 30),
             0.9 : (8, 27),
             1.0 : (6, 24),
             1.1 : (6, 22),
             1.2 : (7, 21),
             1.3 : (7, 19),
             1.4 : (6, 18)}
    rdpath = os.path.join(raw_data_dir, "Tare drag", str(nrun))
    with open(os.path.join(rdpath, "metadata.json")) as f:
        metadata = json.load(f)
    speed = float(metadata["Tow speed (m/s)"])
    nidata = loadmat(os.path.join(rdpath, "nidata.mat"), squeeze_me=True)
    time_ni  = nidata["t"]
    drag = nidata["drag_left"] + nidata["drag_right"]
    drag = drag - np.mean(drag[:2000])
    t1, t2 = times[speed]
    meandrag, x = ts.calcstats(drag, t1, t2, 2000) 
    print("Tare drag =", meandrag, "N at", speed, "m/s")
    if plot:
        plt.figure()
        plt.plot(time_ni, drag, 'k')
        plt.show()
    return speed, meandrag
        
def batch_process_tare_drag(plot=False):
    """Processes all tare drag data."""
    runs = os.listdir("Raw/Tare drag")
    runs = sorted([int(run) for run in runs])
    speed = np.zeros(len(runs))
    taredrag = np.zeros(len(runs))
    for n in range(len(runs)):
        speed[n], taredrag[n] = process_tare_drag(runs[n])
    data = pd.DataFrame()
    data["run"] = runs
    data["tow_speed"] = speed
    data["tare_drag"] = taredrag
    data.to_csv("Data/Processed/Tare drag.csv", index=False)
    if plot:
        plt.figure()
        plt.plot(speed, taredrag, "-ok", markerfacecolor="None")
        plt.xlabel("Tow speed (m/s)")
        plt.ylabel("Tare drag (N)")
        plt.tight_layout()
        plt.show()
        
def plot_tare_drag():
    df = pd.read_csv("Data/Processed/Tare drag.csv")
    plt.figure()
    plt.plot(df.tow_speed, df.tare_drag, "-ok")
    plt.xlabel("Tow speed (m/s)")
    plt.ylabel("Tare drag (N)")
    plt.show()
    
def process_tare_torque(nrun, plot=False):
    """Processes a single tare torque run."""
    print("Processing tare torque run", str(nrun)+"...")
    times = {0 : (35, 86),
             1 : (12, 52),
             2 : (11, 32),
             3 : (7, 30)}
    nidata = loadmat("Data/Raw/Tare torque/" + str(nrun) + "/nidata.mat", 
                     squeeze_me=True)
    # Compute RPM
    time_ni  = nidata["t"]
    angle = nidata["turbine_angle"]
    rpm_ni = fdiff.second_order_diff(angle, time_ni)/6.0
    rpm_ni = ts.smooth(rpm_ni, 8)
    try:
        t1, t2 = times[nrun]
    except KeyError:
        t1, t2 = times[3]
    meanrpm, x = ts.calcstats(rpm_ni, t1, t2, 2000)
    torque = nidata["torque_trans"]
#    torque = torque - np.mean(torque[:2000]) # 2000 samples of zero torque
    meantorque, x = ts.calcstats(torque, t1, t2, 2000)
    print("Tare torque =", meantorque, "Nm at", meanrpm, "RPM")
    if plot:
        plt.figure()
        plt.plot(time_ni, torque)
        plt.xlabel("Time (s)")
        plt.ylabel("Torque (Nm)")
        plt.tight_layout()
        plt.show()
    return meanrpm, -meantorque
    
def batch_process_tare_torque(plot=False):
    """Processes all tare torque data."""
    runs = os.listdir("Data/Raw/Tare torque")
    runs = sorted([int(run) for run in runs])
    rpm = np.zeros(len(runs))
    taretorque = np.zeros(len(runs))
    for n in range(len(runs)):
        rpm[n], taretorque[n] = process_tare_torque(runs[n])
    df = pd.DataFrame()
    df["run"] = runs
    df["rpm"] = rpm
    df["tare_torque"] = taretorque
    df.to_csv("Data/Processed/Tare torque.csv", index=False)
    m, b = np.polyfit(rpm, taretorque, 1)
    print("tare_torque = "+str(m)+"*rpm +", b)
    if plot:
        plt.figure()
        plt.plot(rpm, taretorque, "-ok", markerfacecolor="None")
        plt.plot(rpm, m*rpm + b)
        plt.xlabel("RPM")
        plt.ylabel("Tare torque (Nm)")
        plt.ylim((0, 1))
        plt.tight_layout()
        plt.show()
        
def plot_perf_curves(subplots=True, save=False, savepath="", savetype=".pdf"):
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
        if savepath != "":
            savepath += "/"
        plt.savefig(savepath + "perf_curves" + savetype)
    
def plot_wake_profiles(z_H=0.25, save=False, savepath="", savetype=".pdf"):
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
            plt.savefig(os.path.join(savepath, q+savetype))
    plt.show()
    
def main():
    plot_perf_re_dep(dual_xaxes=True)
        
if __name__ == "__main__":
    if os.getcwd()[-7:] == "Modules":
        print("Changing working directory to experiment root directory")
        os.chdir("../")
    if len(sys.argv) == 3:
        section = sys.argv[1]
        nrun = int(sys.argv[2])
        run = Run(section, nrun)
        run.calc_perf()
        run.calc_wake()
    else:
        main()

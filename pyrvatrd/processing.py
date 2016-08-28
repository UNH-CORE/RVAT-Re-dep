# -*- coding: utf-8 -*-
"""This module contains classes and functions for processing data."""

from __future__ import division, print_function
import numpy as np
from pxl import timeseries as ts
from pxl.timeseries import calc_uncertainty, calc_exp_uncertainty
import matplotlib.pyplot as plt
from scipy.io import loadmat
import multiprocessing as mp
import scipy.stats
from numpy import nanmean, nanstd
from scipy.signal import decimate
from pxl import fdiff
import progressbar
import json
import os
import sys
import pandas as pd

if sys.version_info[0] == 3:
    from urllib.request import urlretrieve
else:
    from urllib import urlretrieve


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
chord = 0.14

# Directory constants
raw_data_dir = os.path.join("Data", "Raw")
processed_data_dir = os.path.join("Data", "Processed")


def calc_b_vec(vel):
    """Calculates the systematic error of a Vectrino measurement (in m/s)
    from their published specs. Returns half the +/- value as b."""
    return 0.5*(0.005*np.abs(vel) + 0.001)


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
         1.2 : (14.0, 27.0),
         1.3 : (13.0, 23.0),
         1.4 : (12.0, 20.0)}


class Run(object):
    """Object that represents a single turbine tow"""
    def __init__(self, section, nrun):
        self.section = section
        nrun = int(nrun)
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
        if self.loaded:
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
            self.calc_wake_per_rev()
            self.calc_wake_uncertainty()
            self.calc_wake_exp_uncertainty()

    def load(self):
        """Loads the data from the run into memory."""
        self.loaded = True
        try:
            with open("Config/raw_data_urls.json") as f:
                raw_data_urls = json.load(f)
        except IOError:
            raw_data_urls = {}
        # Load metadata if it exists
        fpath_metadata = os.path.join(self.raw_dir, "metadata.json")
        if os.path.isfile(fpath_metadata):
            self.load_metadata()
        elif make_remote_name(fpath_metadata) in raw_data_urls:
            self.download_raw("metadata.json")
            self.load_metadata()
        else:
            self.loaded = False
        # Load NI data if it exists
        fpath_nidata = os.path.join(self.raw_dir, "nidata.mat")
        if os.path.isfile(fpath_nidata):
            self.load_nidata()
        elif make_remote_name(fpath_nidata) in raw_data_urls:
            self.download_raw("nidata.mat")
            self.load_nidata()
        else:
            self.loaded = False
        # Load ACS data if it exists
        fpath_acsdata = os.path.join(self.raw_dir, "acsdata.mat")
        if os.path.isfile(fpath_acsdata):
            self.load_acsdata()
        elif make_remote_name(fpath_acsdata) in raw_data_urls:
            self.download_raw("acsdata.mat")
            self.load_acsdata()
        else:
            self.loaded = False
        # Load Vectrino data if it exists
        fpath_vecdata = os.path.join(self.raw_dir, "vecdata.mat")
        if os.path.isfile(fpath_vecdata):
            self.load_vecdata()
        elif make_remote_name(fpath_vecdata) in raw_data_urls:
            self.download_raw("vecdata.mat")
            self.load_vecdata()
        else:
            self.loaded = False

    def load_metadata(self):
        """Loads run metadata."""
        with open(os.path.join(self.raw_dir, "metadata.json")) as f:
            self.metadata = json.load(f)
        self.tow_speed_nom = np.round(self.metadata["Tow speed (m/s)"], decimals=1)
        self.tsr_nom = self.metadata["Tip speed ratio"]
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

    def download_raw(self, name):
        download_raw(self.section, self.nrun, name)

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
        if self.tow_speed_nom == 1.2:
            if self.tsr_nom == 2.9:
                self.t1 = 19
            elif self.tsr_nom > 2.9:
                self.t1 = 23
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
        self.time_vec_all = self.time_vec
        self.time_vec = self.time_vec_all[self.t1*self.sr_vec:self.t2*self.sr_vec]
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
        # Calculate cosine fits for performance
        self.amp_tsr, self.phase_tsr = np.nan, np.nan
        self.amp_cp, self.phase_cp = np.nan, np.nan
        if self.n_revs >= 1:
            angle_seg = np.deg2rad(self.angle)
            self.amp_tsr, self.phase_tsr = ts.find_amp_phase(angle_seg,
                    self.tsr, min_phase=np.deg2rad(10))
            self.amp_cp, self.phase_cp = ts.find_amp_phase(angle_seg, self.cp,
                    min_phase=np.deg2rad(10))

    def print_perf_stats(self):
        print("tow_speed_nom =", self.tow_speed_nom)
        if self.lin_enc:
            print("mean_tow_speed_enc =", self.mean_u_enc)
            print("std_tow_speed_enc =", self.std_u_enc)
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
        self.exp_unc_cp, self.dof_cp = calc_exp_uncertainty(self.n_revs,
                self.std_cp_per_rev, self.unc_cp, self.b_cp)
        # Drag coefficient
        self.exp_unc_cd, self.dof_cd = calc_exp_uncertainty(self.n_revs,
                self.std_cd_per_rev, self.unc_cd, self.b_cd)
        # Tip speed ratio
        self.exp_unc_tsr, self.dof_tsr = calc_exp_uncertainty(self.n_revs,
                self.std_tsr_per_rev, self.unc_tsr, self.b_tsr)

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
        self.u_unf = self.u.copy()
        self.v_unf = self.v.copy()
        self.w_unf = self.w.copy()
        if stdfilt:
        # Do standard deviation filters
            self.u = ts.sigmafilter(self.u, std, passes)
            self.v = ts.sigmafilter(self.v, std, passes)
            self.w = ts.sigmafilter(self.w, std, passes)
        if threshfilt:
            # Do threshold filter on u
            ibad = np.where(self.u > mean_u + fthresh)[0]
            ibad = np.append(ibad, np.where(self.u < mean_u - fthresh)[0])
            self.u[ibad] = np.nan
            # Do threshold filter on v
            ibad = np.where(self.v > mean_v + fthresh)[0]
            ibad = np.append(ibad, np.where(self.v < mean_v - fthresh)[0])
            self.v[ibad] = np.nan
            # Do threshold filter on w
            ibad = np.where(self.w > mean_w + fthresh)[0]
            ibad = np.append(ibad, np.where(self.w < mean_w - fthresh)[0])
            self.w[ibad] = np.nan
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
        print("mean_u/tow_speed_nom =", self.mean_u/self.tow_speed_nom)
        print("std_u/tow_speed_nom =", self.std_u/self.tow_speed_nom)
        print(str(self.nbad)+"/"+str(ntotal), "data points omitted")

    def calc_wake_uncertainty(self):
        """Compute uncertainty for wake statistics."""
        # Mean u
        self.unc_mean_u = np.sqrt(np.nanmean(calc_b_vec(self.u))**2 \
                          + (self.std_u_per_rev/np.sqrt(self.n_revs))**2)
        # Mean v
        self.unc_mean_v = np.sqrt(np.nanmean(calc_b_vec(self.v))**2 \
                          + (self.std_v_per_rev/np.sqrt(self.n_revs))**2)
        # Mean w
        self.unc_mean_w = np.sqrt(np.nanmean(calc_b_vec(self.w))**2 \
                          + (self.std_w_per_rev/np.sqrt(self.n_revs))**2)
        # TODO: Standard deviations
        self.unc_std_u = np.nan

    def calc_wake_exp_uncertainty(self):
        """Calculate expanded uncertainty of wake statistics."""
        # Mean u
        self.exp_unc_mean_u, self.dof_mean_u = calc_exp_uncertainty(
                self.n_revs, self.std_u_per_rev, self.unc_mean_u,
                np.nanmean(calc_b_vec(self.u)))
        # Mean v
        self.exp_unc_mean_v, self.dof_mean_v = calc_exp_uncertainty(
                self.n_revs, self.std_v_per_rev, self.unc_mean_v,
                np.nanmean(calc_b_vec(self.v)))
        # Mean w
        self.exp_unc_mean_w, self.dof_mean_w = calc_exp_uncertainty(
                self.n_revs, self.std_w_per_rev, self.unc_mean_w,
                np.nanmean(calc_b_vec(self.w)))

    def calc_perf_per_rev(self):
        """Computes mean power coefficient over each revolution."""
        angle = self.angle*1
        angle -= angle[0]
        cp = np.zeros(self.n_revs)
        cd = np.zeros(self.n_revs)
        tsr = np.zeros(self.n_revs)
        torque = np.zeros(self.n_revs)
        omega = np.zeros(self.n_revs)
        start_angle = 0.0
        for n in range(self.n_revs):
            end_angle = start_angle + 360
            ind = np.logical_and(angle >= start_angle, end_angle > angle)
            cp[n] = self.cp[ind].mean()
            cd[n] = self.cd[ind].mean()
            tsr[n] = self.tsr[ind].mean()
            torque[n] = self.torque[ind].mean()
            omega[n] = self.omega[ind].mean()
            start_angle += 360
        self.cp_per_rev = cp
        self.std_cp_per_rev = cp.std()
        self.cd_per_rev = cd
        self.std_cd_per_rev = cd.std()
        self.tsr_per_rev = tsr
        self.std_tsr_per_rev = tsr.std()
        self.torque_per_rev = torque
        self.std_torque_per_rev = torque.std()

    def calc_wake_per_rev(self):
        """Computes wake stats per revolution."""
        # Downsample angle measurements to match Vectrino sample rate
        angle = self.angle.copy()
        angle -= angle[0]
        angle = decimate(angle, 10)
        angle = angle[:len(self.u)]
        mean_u = np.zeros(self.n_revs)*np.nan
        mean_v = np.zeros(self.n_revs)*np.nan
        mean_w = np.zeros(self.n_revs)*np.nan
        start_angle = 0.0
        for n in range(self.n_revs):
            end_angle = start_angle + 360
            ind = np.logical_and(angle >= start_angle, angle < end_angle)
            if np.any(ind):
                mean_u[n] = np.nanmean(self.u[ind])
                mean_v[n] = np.nanmean(self.v[ind])
                mean_w[n] = np.nanmean(self.w[ind])
            start_angle += 360
        self.std_u_per_rev = mean_u.std()
        self.std_v_per_rev = mean_v.std()
        self.std_w_per_rev = mean_w.std()

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

    @property
    def summary(self):
        s = pd.Series()
        s["run"] = self.nrun
        s["mean_tow_speed"] = self.mean_u_enc
        s["std_tow_speed"] = self.std_u_enc
        s["t1"] = self.t1
        s["t2"] = self.t2
        s["n_blade_pass"] = self.n_blade_pass
        s["n_revs"] = self.n_revs
        s["mean_tsr"] = self.mean_tsr
        s["mean_cp"] = self.mean_cp
        s["mean_cd"] = self.mean_cd
        s["amp_tsr"] = self.amp_tsr
        s["phase_tsr"] = self.phase_tsr
        s["amp_cp"] = self.amp_cp
        s["phase_cp"] = self.phase_cp
        s["std_tsr"] = self.std_tsr
        s["std_cp"] = self.std_cp
        s["std_cd"] = self.std_cd
        s["std_tsr_per_rev"] = self.std_tsr_per_rev
        s["std_cp_per_rev"] = self.std_cp_per_rev
        s["std_cd_per_rev"] = self.std_cd_per_rev
        s["sys_unc_tsr"] = self.b_tsr
        s["sys_unc_cp"]  = self.b_cp
        s["sys_unc_cd"] = self.b_cd
        s["exp_unc_tsr"] = self.exp_unc_tsr
        s["exp_unc_cp"] = self.exp_unc_cp
        s["exp_unc_cd"] = self.exp_unc_cd
        s["dof_tsr"] = self.dof_tsr
        s["dof_cp"] = self.dof_cp
        s["dof_cd"] = self.dof_cd
        s["t1_wake"] = self.t1_wake
        s["t2_wake"] = self.t2_wake
        s["y_R"] = self.y_R
        s["z_H"] = self.z_H
        s["mean_u"] = self.mean_u
        s["mean_v"] = self.mean_v
        s["mean_w"] = self.mean_w
        s["std_u"] = self.std_u
        s["std_v"] = self.std_v
        s["std_w"] = self.std_w
        s["mean_upup"] = self.mean_upup
        s["mean_upvp"] = self.mean_upvp
        s["mean_upwp"] = self.mean_upwp
        s["mean_vpvp"] = self.mean_vpvp
        s["mean_vpwp"] = self.mean_vpwp
        s["mean_wpwp"] = self.mean_wpwp
        s["k"] = self.k
        s["exp_unc_mean_u"] = self.exp_unc_mean_u
        s["exp_unc_mean_v"] = self.exp_unc_mean_v
        s["exp_unc_mean_w"] = self.exp_unc_mean_w
        return s

    def plot_perf(self, quantity="power coefficient", verbose=True):
        """Plot the run's performance data."""
        qname = quantity
        if verbose:
            print("Plotting {} from {} run {}".format(quantity, self.section,
                  self.nrun))
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
            ylim = None
        if verbose:
            print("Mean TSR: {:.3f}".format(self.mean_tsr))
            print(qname.capitalize(), "statistics:")
            print("Min: {:.3f}, Max: {:.3f}, Mean: {:.3f}".format(np.min(quantity),
                  np.max(quantity), nanmean(quantity)))
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
    def mean_cp(self):
        return self.data.mean_cp
    def process(self, nproc=4, save=True):
        """Process an entire section of data."""
        if nproc > 1:
            self.process_parallel(nproc=nproc)
        else:
            self.process_serial()
        self.data.run = [int(run) for run in self.data.run]
        if save:
            self.data.to_csv(self.processed_path, na_rep="NaN", index=False)
    def process_parallel(self, nproc=4, nruns="all"):
        s = self.name
        runs = self.test_plan["Run"]
        if nruns != "all":
            runs = runs[:nruns]
        pool = mp.Pool(processes=nproc)
        results = [pool.apply_async(process_run, args=(s,n)) for n in runs]
        output = [p.get() for p in results]
        self.data = pd.DataFrame(output)
        pool.close()
    def process_serial(self):
        s = self.name
        runs = self.test_plan["Run"]
        summaries = []
        for nrun in runs:
            r = Run(s, int(nrun))
            summaries.append(r.summary)
        self.data = pd.DataFrame(summaries)


def process_run(section, nrun):
    run = Run(section, nrun)
    return run.summary


def process_latest_run(section):
    """Automatically detects the most recently acquired run and processes it,
    printing a summary to the shell.
    """
    print("Processing latest run in", section)
    raw_dir = os.path.join("Data", "Raw", section)
    dirlist = [os.path.join(raw_dir, d) for d in os.listdir(raw_dir) \
               if os.path.isdir(os.path.join(raw_dir, d))]
    dirlist = sorted(dirlist, key=os.path.getmtime, reverse=True)
    for d in dirlist:
        try:
            nrun = int(os.path.split(d)[-1])
            break
        except ValueError:
            print(d, "is not a properly formatted directory")
    print("\nSummary for {} run {}:".format(section, nrun))
    print(Run(section, nrun).summary)


def load_test_plan_section(section):
    df = pd.read_csv(os.path.join("Config", "Test plan", section+".csv"))
    df = df.dropna(how="all", axis=1).dropna(how="all", axis=0)
    if "Run" in df:
        df["Run"] = df["Run"].astype(int)
    return df


def process_section(name):
    s = Section(name)
    s.process()


def batch_process_all():
    """Batch processes all sections."""
    sections = ["Perf-0.3", "Perf-0.4", "Perf-0.5",
                "Perf-0.6", "Perf-0.7", "Perf-0.8",
                "Perf-0.9", "Perf-1.0", "Perf-1.1",
                "Perf-1.2", "Perf-1.3", "Wake-0.4",
                "Wake-0.6", "Wake-0.8", "Wake-1.0",
                "Wake-1.2"]
    for section in sections:
        print("Processing {}".format(section))
        process_section(section)


def process_tare_drag(nrun, plot=False):
    """Process a single tare drag run."""
    print("Processing tare drag run", nrun)
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


def process_tare_torque(nrun, plot=False):
    """Processes a single tare torque run."""
    print("Processing tare torque run", nrun)
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


def make_remote_name(local_path):
    return "_".join(local_path.split("\\")[-3:])


def download_raw(section, nrun, name):
    """Download a run's raw data.

    `name` can be either the file name with extension, or
      * `"metadata"` -- Metadata in JSON format
      * `"nidata"` -- Data from the NI DAQ system
      * `"acsdata"` -- Data from the tow tank's motion controller
      * `"vecdata"` -- Data from the Nortek Vectrino
    """
    if name == "metadata":
        filename = "metadata.json"
    elif name in ["vecdata", "nidata", "acsdata"]:
        filename = name + ".mat"
    else:
        filename = name
    print("Downloading", filename, "from", section, "run", nrun)
    local_dir = os.path.join("Data", "Raw", section, str(nrun))
    if not os.path.isdir(local_dir):
        os.makedirs(local_dir)
    local_path = os.path.join(local_dir, filename)
    remote_name = make_remote_name(local_path)
    with open("Config/raw_data_urls.json") as f:
        urls = json.load(f)
    url = urls[remote_name]
    pbar = progressbar.ProgressBar()
    def download_progress(blocks_transferred, block_size, total_size):
        percent = int(blocks_transferred*block_size*100/total_size)
        try:
            pbar.update(percent)
        except ValueError:
            pass
        except AssertionError:
            pass
    pbar.start()
    urlretrieve(url, local_path, reporthook=download_progress)
    pbar.finish()

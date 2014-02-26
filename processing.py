# -*- coding: utf-8 -*-
"""
Created on Sun Sep 29 12:51:12 2013

This code processes the data from the 2013 Fall Re dependence experiment

Tare drag and torque will have to be done after APS...

@author: Pete
"""
from __future__ import division, print_function
import numpy as np
import timeseries as ts
import matplotlib.pyplot as plt
from scipy.io import loadmat
from styleplot import styleplot
import json
import os
import fdiff

folders = {"Perf-0.4" : "Performance/U_0.4",
           "Perf-0.6" : "Performance/U_0.6",
           "Perf-0.8" : "Performance/U_0.8",
           "Perf-1.0" : "Performance/U_1.0",
           "Perf-1.2" : "Performance/U_1.2",
           "Perf-1.4" : "Performance/U_1.4",
           "Wake-0.5" : "Wake/U_0.5",
           "Shakedown" : "Shakedown"}
           
# Constants
H = 1.0
D = 1.0
A = D*H
R = D/2
rho = 1000.0
tare_torque = 1.0 # Approximate

# Tare drag in Newtons for each speed
tare_drag = {0.5 : 8.0110528524,
             1.0 : 45.462475976,
             1.5 : 107.981,
             2.0 : 193.399,
             0.3 : 4.0,
             0.4 : 5.0,
             0.7 : 30.0,
             0.8 : 35.0,
             0.9 : 40.0,
             1.1 : 55,
             0.6 : 20.0}
             
times = {0.4 : (20.0, 60.0),
         0.6 : (18.0, 42.0),
         0.8 : (18.0, 34.0),
         1.0 : (15.0, 29.0),
         1.2 : (14.0, 24.0),
         1.4 : (12.0, 20.0)}

cfd_path = "C:/Users/Pete/Google Drive/OpenFOAM/pete-2.2.2/run/unh-rvat-2d_Re-dep/processed/"

class Run(object):
    """Object that represents a single turbine tow"""
    def __init__(self, section, nrun):
        self.section = section
        self.nrun = nrun
        self.folder = folders[section] + "/" + str(nrun)
        self.loaded = False
        self.t2found = False
        
    def load(self):
        """Loads the data from the run into memory"""
        # Load metadata
        with open(self.folder + "/" + "metadata.json") as f:
            self.metadata = json.load(f)
        self.U_nom = np.round(self.metadata["Tow speed (m/s)"], decimals=1)
        self.y_R = self.metadata["Vectrino y/R"]
        # Load NI data
        nidata = loadmat(self.folder + "/" + "nidata.mat", squeeze_me=True)
        self.t_ni = nidata["t"]
        if "carriage_pos" in nidata:
            self.lin_enc = True
            self.carriage_pos = nidata["carriage_pos"]
            self.U_ni = fdiff.second_order_diff(self.carriage_pos, self.t_ni)
            self.U_ref = self.U_ni
        else:
            self.lin_enc = False
            self.U_ref = self.U_nom
        self.sr_ni = (1.0/(self.t_ni[1] - self.t_ni[0]))
        self.torque = nidata["torque_trans"]
        self.drag = nidata["drag_left"] + nidata["drag_right"]
        # Load ACS data
        acsdata = loadmat(self.folder + "/" + "acsdata.mat", squeeze_me=True)
        self.U_acs = acsdata["carriage_vel"]
        self.rpm_acs = acsdata["turbine_rpm"]
        self.omega_acs = self.rpm_acs*2*np.pi/60.0
        self.t_acs = acsdata["t"]
        self.omega_acs_interp = np.interp(self.t_ni, self.t_acs, self.omega_acs)
        self.rpm_acs_interp = self.omega_acs_interp*60.0/(2*np.pi)
        # Remove offsets from torque and drag
        t0 = 2
        self.torque = self.torque - np.mean(self.torque[:self.sr_ni*t0])
        self.drag = self.drag - np.mean(self.drag[0:self.sr_ni*t0])
        # Subtract tare drag
        self.drag = self.drag - tare_drag[self.U_nom]
        # Compute RPM and omega
        self.angle = nidata["turbine_angle"]
        self.rpm_ni = fdiff.second_order_diff(self.angle, self.t_ni)/6.0
        self.rpm_ni = ts.smooth(self.rpm_ni, 50)
        self.omega_ni = self.rpm_ni*2*np.pi/60.0
        # Add tare torque
        tare_torque = 0.00174094659759*self.rpm_acs_interp + 0.465846267394
        self.torque = self.torque + tare_torque
        # Compute power
        self.power = self.torque*self.omega_acs_interp
        self.tsr = self.omega_acs_interp*R/self.U_ref
        # Compute power and drag coefficients
        self.cp = self.power/(0.5*rho*A*self.U_ref**3)
        self.cd = self.drag/(0.5*rho*A*self.U_ref**2)
        # Load Vectrino data
        try:
            vecdata = loadmat(self.folder + "/" + "vecdata.mat", 
                              squeeze_me=True)
            self.sr_vec = 200.0
            self.t_vec = vecdata["t"]
            self.u = vecdata["u"]
            self.v = vecdata["v"]
            self.w = vecdata["w"]
        except IOError:
            self.vecdata = None
        # Put in some guesses for t1 and t2
        self.t1, self.t2 = times[self.U_nom]
        self.loaded = True
        
    def find_t2(self):
        sr = self.sr_ni
        angle1 = self.angle[sr*self.t1]
        angle2 = self.angle[sr*self.t2]
        n3rdrevs = np.floor((angle2-angle1)/120.0)
        self.nrevs = np.floor((angle2-angle1)/360.0)
        angle2 = angle1 + n3rdrevs*120
        t2i = np.where(np.round(self.angle)==np.round(angle2))[0][0]
        t2 = self.t_ni[t2i]
        self.t2 = np.round(t2, decimals=2)
        self.t2found = True
        
    def calcperf(self):
        """Calculates mean performance based on data between t0 and t1"""
        if not self.loaded:
            self.load()
        self.find_t2()
        self.meantsr, x = ts.calcstats(self.tsr, self.t1, self.t2, self.sr_ni)
        self.meancd, x = ts.calcstats(self.cd, self.t1, self.t2, self.sr_ni)
        self.meancp, x = ts.calcstats(self.cp, self.t1, self.t2, self.sr_ni)
        print("U_nom =", self.U_nom)
        if self.lin_enc:
            self.meanu_enc, x = ts.calcstats(self.U_ni, self.t1, self.t2, self.sr_ni)
            print("U_enc =", self.meanu_enc)
        print("tsr =", self.meantsr)
        print("C_P =", self.meancp)
        print("C_D =", self.meancd)
        
    def filter_wake(self):
        std = 3
        passes = 2
        self.u_f = self.u*1
        self.u_f[200*self.t1:200*self.t2] = \
                ts.sigmafilter(self.u[200*self.t1:200*self.t2], std, passes)
        meanu, x = ts.calcstats(self.u, self.t1, self.t2, self.sr_vec)
        ibad = np.where(self.u > 1.5*meanu)[0]
        ibad = np.append(ibad, np.where(self.u < 0.5*meanu)[0])
        i = np.where(np.logical_and(ibad > self.t1*200, ibad < self.t2*200))[0]
#        self.u_f[ibad[i]] = np.nan
        self.nbad = len(i)
        self.v_f = self.v*1
        self.v_f[200*self.t1:200*self.t2] = \
                ts.sigmafilter(self.v[200*self.t1:200*self.t2], std, passes)
        self.w_f = self.w*1
        self.w_f[200*self.t1:200*self.t2] = \
                ts.sigmafilter(self.w[200*self.t1:200*self.t2], std, passes)
        
    def calcwake(self):
        if not self.loaded:
            self.load()
        if not self.t2found:
            self.find_t2()
        self.filter_wake()
        self.meanu, self.stdu = ts.calcstats(self.u_f, self.t1, self.t2, self.sr_vec)
        self.meanv, self.stdv = ts.calcstats(self.v_f, self.t1, self.t2, self.sr_vec)
        self.meanw, self.stdw = ts.calcstats(self.w_f, self.t1, self.t2, self.sr_vec)
        uv = (self.u_f - self.meanu)*(self.v_f - self.meanv)
        self.meanuv, self.stduv = ts.calcstats(uv, self.t1, self.t2, self.sr_vec)
        print("y/R =", self.y_R)
        print("U_vec/U =", self.meanu/self.U)
        print("std_u/U =", self.stdu/self.U)
        print(self.nbad, "data points omitted")
        
    def plotperf(self, quantity="torque"):
        """Plots the run's data"""
        if not self.loaded:
            self.load()
        if quantity == "drag":
            quantity = self.drag
        elif quantity == "torque":
            quantity = self.torque
        plt.figure()
        plt.plot(self.t_ni, quantity, 'k')
        styleplot()
        
    def plotwake(self):
        """Plot mean and standard deviation for each velocity component"""
        if not self.loaded:
            self.load()
        plt.figure()
        self.filter_wake()
        plt.plot(self.t_vec, self.u_f, 'k')
        
    def plotacs(self):
        if not self.loaded:
            self.load()
        plt.figure()
        plt.plot(self.t_acs, self.rpm_acs)
        plt.hold(True)
        plt.plot(self.t_ni, self.rpm_ni)
        plt.show()
        
    def plotvel(self):
        if not self.loaded:
            self.load()
        plt.figure()
        plt.plot(self.t_ni, self.U_ni)
        styleplot()
        plt.show()

class PerfCurve(object):
    """Object that represents a performance curve."""
    def __init__(self, U):
        self.U = U
        self.folder = folders["Perf-" + str(U)]
        dirconts = os.listdir(self.folder)
        self.runs = []
        for item in dirconts:
            try:
                self.runs.append(int(item))
            except ValueError:
                pass
        self.runs.sort()
        
    def process(self, reprocess=True):
        """Calculates power and drag coefficients for each run"""
        if not reprocess:
            runsdone = np.load(self.folder+"/Processed/runs.npy")
            tsr_old = np.load(self.folder+"/Processed/tsr.npy")
            cp_old = np.load(self.folder+"/Processed/cp.npy")
            cd_old = np.load(self.folder+"/Processed/cp.npy")
        tsr = np.zeros(len(self.runs))
        cp = np.zeros(len(self.runs))
        cd = np.zeros(len(self.runs))
        for n in range(len(self.runs)):
            nrun = self.runs[n]
            if reprocess or nrun not in runsdone:
                print("Processing run " + str(self.runs[nrun]) + "...")
                run = Run("Perf-{:0.1f}".format(self.U), nrun)
                run.calcperf()
                tsr[n] = run.meantsr
                cp[n] = run.meancp
                cd[n] = run.meancd
            else:
                tsr[n] = tsr_old[np.where(runsdone==nrun)[0]]
                cp[n] = cp_old[np.where(runsdone==nrun)[0]]
                cd[n] = cd_old[np.where(runsdone==nrun)[0]]
        # Save processed data in the folder
        if not os.path.isdir(self.folder+"/Processed"):
            os.mkdir(self.folder+"/Processed")
        np.save(self.folder+"/Processed/runs.npy", self.runs)
        np.save(self.folder+"/Processed/tsr.npy", tsr)
        np.save(self.folder+"/Processed/cp.npy", cp)
        np.save(self.folder+"/Processed/cd.npy", cd)
        
    def plotcp(self, newfig=True, show=True, save=False, figname="test.pdf"):
        """Generates power coefficient curve plot."""
        # Check to see if processed data exists and if not, process it
        try:
            self.tsr = np.load(self.folder+"/Processed/tsr.npy")
            self.cp = np.load(self.folder+"/Processed/cp.npy")
        except IOError:
            self.process()
        if newfig:
            plt.figure()
        plt.plot(self.tsr, self.cp, "-ok", markerfacecolor="None")
        plt.xlabel(r"$\lambda$")
        plt.ylabel(r"$C_P$")
        plt.grid(True)
        styleplot()
        if show:
            plt.show()
        if save:
            plt.savefig(figname)
        
class WakeProfile(object):
    def __init__(self, U):
        self.U = U
        self.section = "Wake-U_" + str(U)
        self.folder = folders[self.section]
        self.runs = range(45)
        
    def process(self):
        """Runs through data to calculate statistics"""
        meanu = np.zeros(len(self.runs))
        meanv = np.zeros(len(self.runs))
        meanw = np.zeros(len(self.runs))
        stdu = np.zeros(len(self.runs))
        meanuv = np.zeros(len(self.runs))
        y_R = np.zeros(len(self.runs))
        for n in self.runs:
            run = Run(self.section, n)
            run.calcwake()
            meanu[n] = run.meanu
            meanv[n] = run.meanv
            meanw[n] = run.meanw
            stdu[n] = run.stdu
            meanuv[n] = run.meanuv
            y_R[n] = run.y_R
        # Save processed data
        if not os.path.isdir(self.folder+"/Processed"):
            os.mkdir(self.folder+"/Processed")
        np.save(self.folder+"/Processed/meanu.npy", meanu)
        np.save(self.folder+"/Processed/meanv.npy", meanv)
        np.save(self.folder+"/Processed/meanw.npy", meanw)
        np.save(self.folder+"/Processed/stdu.npy", stdu)
        np.save(self.folder+"/Processed/meanuv.npy", meanuv)
        np.save(self.folder+"/Processed/y_R.npy", y_R)
        
    def load(self):
        """Loads the processed data"""
        self.y_R = np.load(self.folder+"/Processed/y_R.npy")
        self.meanu = np.load(self.folder+"/Processed/meanu.npy")
        self.meanuv = np.load(self.folder+"/Processed/meanuv.npy")
        self.stdu = np.load(self.folder+"/Processed/stdu.npy")
        
    def plot(self, quantity, newfig=True, show=True, save=False, 
             savepath="", savetype=".pdf", linetype='--ok'):
        """Plots some quantity"""
        y_R = np.load(self.folder+"/Processed/y_R.npy")
        q = np.load(self.folder+"/Processed/"+quantity+".npy")
        loc = 1
        if quantity == "meanu":
            q = q/self.U
            ylab = r"$U/U_\infty$"
            loc = 3
        if quantity == "meanw":
            q = q/self.U
            ylab = r"$U/U_\infty$"
            loc = 4
        if quantity == "meanv":
            q = q/self.U
            ylab = r"$V/U_\infty$"
            loc=4
        if quantity == "stdu":
            q = q/self.U
            ylab = r"$\sigma_u/U_\infty$"
        if quantity is "meanuv":
            q = q/(self.U**2)
            ylab = r"$\overline{u'v'}/U_\infty^2$" 
        if newfig:
            if quantity == "meanu":
                plt.figure(figsize=(11,6))
            else: plt.figure()
            plt.ylabel(ylab)
            plt.xlabel(r"$y/R$")
            plt.grid()
        plt.plot(y_R, q, "-.^k", label=r"$Re_D=0.5 \times 10^6$")
        plot_old_wake(quantity, y_R)
        plt.legend(loc=loc)
        styleplot()
        if show:
            plt.show()
        if save:
            plt.savefig(savepath+quantity+"_Re_dep_exp"+savetype)

def calc_re_dep():
    runs = [3, 4, 5, 8] # this will need to be fixed later...
    cp = np.zeros(len(runs))
    cd = np.zeros(len(runs))
    U = np.zeros(len(runs))
    for n in range(len(runs)):
        run = Run("Shakedown", runs[n])
        run.calcperf()
        cp[n] = run.meancp
        cd[n] = run.meancd
        U[n] = run.U
    if not os.path.isdir("Shakedown/Processed"):
        os.mkdir("Shakedown/Processed")
    np.save("Shakedown/Processed/cp.npy", cp)
    np.save("Shakedown/Processed/cd.npy", cd)
    np.save("Shakedown/Processed/U.npy", U)
    
def plot_re_dep(save=False, savepath=""):
    cp = np.load("Shakedown/Processed/cp.npy")
    cd = np.load("Shakedown/Processed/cd.npy")
    U = np.load("Shakedown/Processed/U.npy")
    Re_D = U*D/1e-6
    plt.figure()
    plt.plot(Re_D, cp/cp[1], '-ok', markerfacecolor="none", label="Experiment")
    plt.hold(True)
    plot_cfd_perf("cp")
    plt.xlabel(r"$Re_D$")
    plt.ylabel(r"$C_P/C_{P0}$")
    plt.ylim((0.6,1.6))
    ax = plt.gca()
    ax.xaxis.major.formatter.set_powerlimits((0,0)) 
    plt.grid()
    plt.legend(loc=2)
    styleplot()
    if save:
        plt.savefig(savepath+"re_dep_cp.pdf")
    plt.figure()
    plt.plot(Re_D, cd/cd[1], '-ok', markerfacecolor="none", label="Experiment")
    plt.xlabel(r"$Re_D$")
    plt.ylabel(r"$C_D/C_{D0}$")
    plt.hold(True)
    plot_cfd_perf("cd")
    plt.ylim((0.8,1.1))
    plt.grid()
    plt.legend(loc=3)
    ax = plt.gca()
    ax.xaxis.major.formatter.set_powerlimits((0,0)) 
    styleplot()
    plt.show()
    if save:
        plt.savefig(savepath+"re_dep_cd.pdf")
    
def plot_old_wake(quantity, y_R):
    plt.hold(True)
    runs = range(32, 77)
    ind = [run-1 for run in runs]
    f = "../2013.03 VAT/Processed/"+quantity+".npy"
    q = np.load(f)[ind]
    plt.plot(y_R, q, '--ok', label=r"$Re_D=1.0 \times 10^6$", 
             markerfacecolor="none")
             
def plot_cfd_perf(quantity="cp"):
    Re_D = np.load(cfd_path+"Re_D.npy")
    q = np.load(cfd_path+quantity+".npy")
    plt.plot(Re_D, q/q[1], "--^k", label="Simulation")
    
        
if __name__ == "__main__":
    import sys
    if len(sys.argv) == 3:
        section = sys.argv[1]
        nrun = int(sys.argv[2])
        run = Run(section, nrun)
    else:
        run = Run("Perf-0.4", 7)
    plt.close("all")
    p = "C:/Users/Pete/Google Drive/Research/Presentations/2013.11.24 APS-DFD/Figures/"
#    run = Run("Perf-0.4", 23)
#    run.plotperf("torque")
#    run.calcperf()
#    run.plotacs()
    pc = PerfCurve(0.4)
    pc.process(reprocess=False)
    pc.plotcp()
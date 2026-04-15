import json
import matplotlib.pyplot as plt
import scipy.io
import pandas as pd
from scipy.optimize import curve_fit
import numpy as np

# Enable true LaTeX rendering
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Computer Modern Roman"],
})

class TGAnalysis:
    def __init__(self, json_path):
        self.json_path = json_path
        self.load_data()

    def filter_time(self, x, y):
        # 1. Find where the difference is negative
        decreasing_mask = np.diff(x) < 0
    
        # 2. Find the index of the first True
        # np.argmax returns the index of the first maximum value (True)
        first_decrease_idx = np.argmax(decreasing_mask)
    
        first_part_time = x[:first_decrease_idx + 1]
        first_part_signal = y[:first_decrease_idx + 1]
    
        return first_part_time, first_part_signal

    def load_data(self):
        # open and load JSON file
        with open(self.json_path, "r") as f:
            self.json_data = json.load(f)

        self.data_path = self.json_data["main_path"]
        scans_only = [value for key, value in self.json_data.items() if key.startswith("Scan")]

        self.df_scans = pd.DataFrame(scans_only)

    def get_data_E_const(self, energy_eV):
        df_const_E = self.df_scans[(self.df_scans["flag_secure"] == True) & (self.df_scans["Energy_eV"] == energy_eV) & (self.df_scans["repeated"] == False)]
        df_const_E = df_const_E.sort_values(by="XUV_intensity_uJ")
        scans_const_E = self.df_const_E["Scan"].to_numpy()
        I_const_E = self.df_const_E["XUV_intensity_uJ"].to_numpy()
        E_const_E = self.df_const_E["Energy_eV"].to_numpy()
        filter = self.df_const_E["filter"].to_numpy()

        return scans_const_E, E_const_E, I_const_E, filter

    def get_data_I_const(self, intensity_uJ):
        df_const_I = self.df_scans[(self.df_scans["flag_secure"] == True) & (self.df_scans["XUV_intensity_uJ"] == intensity_uJ) & (self.df_scans["repeated"] == False)]
        df_const_I = df_const_I.sort_values(by="Energy_eV")
        scans_const_I = df_const_I["Scan"].to_numpy()
        E_const_I = df_const_I["Energy_eV"].to_numpy()
        I_const_I = df_const_I["XUV_intensity_uJ"].to_numpy()
        filter = df_const_I["filter"].to_numpy()

        return scans_const_I, E_const_I, I_const_I, filter

    def get_data_scan(self, params_scan):
        if(params_scan["E"] == "all" and params_scan["I"] != "all"):
            self.scans, self.energies, self.intensities, filter = self.get_data_I_const(params_scan["I"])
        elif(params_scan["E"] != "all" and params_scan["I"] == "all"):
            self.scans, self.energies, self.intensities, filter = self.get_data_E_const(params_scan["E"])
        else:
            raise ValueError("Invalid parameters for scan")

        self.time_scans = []
        self.tgsignal_scans = []
        self.iomsh_ave_scans = []
        self.iomsh_ave2_scans = []

        for i, scan in enumerate(self.scans):
           data_scan = scipy.io.loadmat(self.data_path + f"{scan}.mat")
       
           self.time = data_scan["R"][0][0][0][:,0]
           #tgsum = data_scan["R"][0][0][1][:,0]
           #tgsumnorm = data_scan["R"][0][0][2][:,0]
           #fluxfacsum = data_scan["R"][0][0][3][:,0]
           #cbar = data_scan["R"][0][0][4][:,0]
           #rbar = data_scan["R"][0][0][5][:,0]
           #cg = data_scan["R"][0][0][6][:,0]
           #rg = data_scan["R"][0][0][7][:,0]
           self.iomsh_ave = data_scan["R"][0][0][8][:,0]
           self.iomsh_ave2 = data_scan["R"][0][0][9][:,0]
           self.tgsignal = data_scan["R"][0][0][10][:,0]
       
           #mask = time < 10
           #time = time[mask]
           #tgsignal = tgsignal[mask]
           
           if filter[i]:
               self.time, self.tgsignal = self.filter_time(self.time, self.tgsignal)

           self.tgsignal = self.tgsignal - np.min(self.tgsignal)
           self.tgsignal = (self.tgsignal - np.min(self.tgsignal)) / (np.max(self.tgsignal) - np.min(self.tgsignal))
           self.time = self.time - self.time[np.argmax(self.tgsignal)]

           self.time_scans.append(self.time)
           self.tgsignal_scans.append(self.tgsignal)
           self.iomsh_ave_scans.append(self.iomsh_ave)
           self.iomsh_ave2_scans.append(self.iomsh_ave2)

    def plot_phase_space(self):
        df_phase_space = self.df_scans[(self.df_scans["flag_secure"] == True) & (self.df_scans["repeated"] == False)]
        scans_phase_space = df_phase_space["Scan"].to_numpy()
        E_phase_space = df_phase_space["Energy_eV"].to_numpy()
        I_phase_space = df_phase_space["XUV_intensity_uJ"].to_numpy()

        plt.figure(figsize=(6, 5))
        plt.scatter(E_phase_space, I_phase_space, marker='o', color="olive")
        for i, (E, I) in enumerate(zip(E_phase_space, I_phase_space)):
            plt.annotate("S" + scans_phase_space[i][4:], (E, I), xytext=(5, 5), textcoords='offset points', fontsize=8)
        plt.xlabel("Energy (eV)")
        plt.ylabel("Intensity (uJ)")
        plt.title("Phase Space")
        plt.grid(linestyle='--', alpha=0.5)
        plt.xlim(E_phase_space.min()-2, E_phase_space.max()+2)
        plt.show()

    def plot_stacked_signals(self, limits_time=None, limits_signal=None, ylog_scale = False, plot_ind=None):
        plt.figure(figsize=(6, 5))
        for i, (time, tgsignal) in enumerate(zip(self.time_scans, self.tgsignal_scans)):
            if plot_ind is None:
                plt.plot(time, tgsignal, label=f"S{self.scans[i][4:]}, E={self.energies[i]:.1f} eV, I={self.intensities[i]:.1f} $\\mu$J")
            else:
                if plot_ind == "all":
                    plt.plot(time, tgsignal, label=f"S{self.scans[i][4:]}, E={self.energies[i]:.1f} eV, I={self.intensities[i]:.1f} $\\mu$J")
                else:
                    if i in plot_ind:
                        plt.plot(time, tgsignal, label=f"S{self.scans[i][4:]}, E={self.energies[i]:.1f} eV, I={self.intensities[i]:.1f} $\\mu$J")
        plt.xlabel("Time (ps)")
        plt.ylabel("TG Signal")
        plt.title("Stacked Signals")
        if limits_time is not None:
            plt.xlim(limits_time[0], limits_time[1])
        if limits_signal is not None:
            plt.ylim(limits_signal[0], limits_signal[1])
        plt.grid(linestyle='--', alpha=0.5)
        plt.legend(loc="upper center", bbox_to_anchor=(0.5, -0.15), ncol=3, fontsize=7)
        plt.tight_layout(rect=[0, 0.08, 1, 1])
        if ylog_scale:
            plt.yscale("log")
        plt.show()
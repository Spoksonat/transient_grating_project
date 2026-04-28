"""Core tools for transient grating (TG) analysis.

This module provides the :class:`TGAnalysis` class to load scan metadata from
JSON files, extract TG traces from MATLAB ``.mat`` files, fit different model
functions, and generate diagnostic plots.
"""

import json
import matplotlib.pyplot as plt
import scipy.io
import pandas as pd
from scipy.optimize import curve_fit
import numpy as np
import scipy.special as special

# Enable true LaTeX rendering
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Computer Modern Roman"],
})

class TGAnalysis:
    """Transient grating analysis pipeline.

    Parameters
    ----------
    json_path : str
        Path to the metadata JSON file. The file must contain a ``main_path``
        entry with the folder of ``.mat`` scans and multiple ``ScanXXX`` items.
    """

    def __init__(self, json_path):
        """Initialize the analysis object and load metadata."""
        self.json_path = json_path
        self.load_data()

    def filter_time(self, x, y):
        """Trim arrays before the first detected time-axis reset.

        Parameters
        ----------
        x : numpy.ndarray
            Time axis.
        y : numpy.ndarray
            TG signal values associated with ``x``.

        Returns
        -------
        tuple[numpy.ndarray, numpy.ndarray]
            Filtered time and signal arrays up to the first negative time step.
        """
        # 1. Find where the difference is negative
        decreasing_mask = np.diff(x) < 0
    
        # 2. Find the index of the first True
        # np.argmax returns the index of the first maximum value (True)
        first_decrease_idx = np.argmax(decreasing_mask)
    
        first_part_time = x[:first_decrease_idx + 1]
        first_part_signal = y[:first_decrease_idx + 1]
    
        return first_part_time, first_part_signal

    def load_data(self):
        """Read metadata JSON and build the scans dataframe."""
        # open and load JSON file
        with open(self.json_path, "r") as f:
            self.json_data = json.load(f)

        self.data_path = self.json_data["main_path"]
        scans_only = [value for key, value in self.json_data.items() if key.startswith("Scan")]

        self.df_scans = pd.DataFrame(scans_only)

    def get_data_E_const(self, energy_eV):
        """Select scans at fixed energy and variable intensity.

        Parameters
        ----------
        energy_eV : float
            Target energy in eV.

        Returns
        -------
        tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray, numpy.ndarray]
            Scan labels, energies, intensities, and per-scan filter flags.
        """
        df_const_E = self.df_scans[(self.df_scans["flag_secure"] == True) & (self.df_scans["Energy_eV"] == energy_eV) & (self.df_scans["repeated"] == False)]
        df_const_E = df_const_E.sort_values(by="XUV_intensity_uJ")
        scans_const_E = self.df_const_E["Scan"].to_numpy()
        I_const_E = self.df_const_E["XUV_intensity_uJ"].to_numpy()
        E_const_E = self.df_const_E["Energy_eV"].to_numpy()
        filter = self.df_const_E["filter"].to_numpy()

        return scans_const_E, E_const_E, I_const_E, filter

    def get_data_I_const(self, intensity_uJ):
        """Select scans at fixed intensity and variable energy.

        Parameters
        ----------
        intensity_uJ : float
            Target XUV intensity in microjoules.

        Returns
        -------
        tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray, numpy.ndarray]
            Scan labels, energies, intensities, and per-scan filter flags.
        """
        df_const_I = self.df_scans[(self.df_scans["flag_secure"] == True) & (self.df_scans["XUV_intensity_uJ"] == intensity_uJ) & (self.df_scans["repeated"] == False)]
        df_const_I = df_const_I.sort_values(by="Energy_eV")
        scans_const_I = df_const_I["Scan"].to_numpy()
        E_const_I = df_const_I["Energy_eV"].to_numpy()
        I_const_I = df_const_I["XUV_intensity_uJ"].to_numpy()
        filter = df_const_I["filter"].to_numpy()

        return scans_const_I, E_const_I, I_const_I, filter

    def get_data_scan(self, params_scan):
        """Load and preprocess TG traces for a scan selection.

        Parameters
        ----------
        params_scan : dict
            Scan-selection configuration. Supported patterns are:
            ``{"E": "all", "I": value}`` or ``{"E": value, "I": "all"}``.

        Raises
        ------
        ValueError
            If both ``E`` and ``I`` are fixed (or both are ``"all"``).
        """
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
           self.iomsh_ave = np.mean(data_scan["R"][0][0][8][:,0])
           self.iomsh_ave2 = np.mean(data_scan["R"][0][0][9][:,0])
           self.tgsignal = data_scan["R"][0][0][10][:,0]
           
           if filter[i]:
               self.time, self.tgsignal = self.filter_time(self.time, self.tgsignal)

           self.tgsignal = self.tgsignal - np.min(self.tgsignal)
           self.tgsignal = self.tgsignal / np.sum(self.tgsignal)
           #self.time = self.time - self.time[np.argmax(self.tgsignal)]

           mask = (self.time < 2.5) & (self.time > -1)
           self.time = self.time[mask]
           self.tgsignal = self.tgsignal[mask]

           self.time_scans.append(self.time)
           self.tgsignal_scans.append(self.tgsignal)
           self.iomsh_ave_scans.append(self.iomsh_ave)
           self.iomsh_ave2_scans.append(self.iomsh_ave2)

    def model1(self, t, amp1, t0, k_tau, sigma, amp2):
        """Model 1: mono-exponential decay convolved with Gaussian response."""
        exp1 = np.exp(-(t - t0)*k_tau)
        exp2 = np.exp(0.5*(k_tau*sigma)**2)
        erf1 = special.erf((t - t0 - k_tau*sigma**2)/(np.sqrt(2)*sigma))
        erf2 = special.erf((t-t0)/(np.sqrt(2)*sigma))
        function = amp1*exp1*exp2*(1 + erf1) + amp2*(1 + erf2)
        return function

    def model2(self, t, amp1, t0, k1, sigma, amp2, omega, phi, amp3, k2):
        """Model 2: bi-exponential decay: one decay plus one damped oscillatory contribution."""
        exp1 = np.exp(-(t - t0)*k1)
        exp2 = np.exp(0.5*(k1*sigma)**2)
        exp3 = np.exp(-k2*(t - t0))
        exp4 = np.exp(-0.5*(k2*sigma)**2)
        erf1 = special.erf((t - t0 - k1*sigma**2)/(np.sqrt(2)*sigma))
        erf2 = special.erf((t-t0)/(np.sqrt(2)*sigma))
        erf3 = special.erf((t-t0-k2*sigma**2)/(np.sqrt(2)*sigma))
        osc_factor = np.sin(omega*t - phi) - np.cos(omega*t - phi)
        function = amp1*exp1*exp2*(1 + erf1) + amp2*(1 + erf2) + amp3*exp3*exp4*(1 + erf3)*osc_factor 
        return function

    def model3(self, t, amp1, amp2, amp3, amp4, t0, k1, k2, sigma, omega, phi):
        """Model 3: multi-exponential decay: two decay plus one damped oscillatory contribution."""
        k10 = k1 + 1/t0
        k20 = k2 + 1/t0
        exp1 = np.exp(-(t - t0)*k1)
        exp2 = np.exp(-(t - t0)*k2)
        exp3 = np.exp(-(t - t0)*k10)
        exp4 = np.exp(-(t - t0)*k20)
        exp5 = np.exp(-0.5*(k1*sigma)**2)
        exp6 = np.exp(-0.5*(k2*sigma)**2)
        exp7 = np.exp(-0.5*(k10*sigma)**2)
        exp8 = np.exp(-0.5*(k20*sigma)**2)
        erf1 = special.erf((t - t0 - k1*sigma**2)/(np.sqrt(2)*sigma))
        erf2 = special.erf((t - t0 - k2*sigma**2)/(np.sqrt(2)*sigma))
        erf3 = special.erf((t-t0-k10*sigma**2)/(np.sqrt(2)*sigma))
        erf4 = special.erf((t-t0-k20*sigma**2)/(np.sqrt(2)*sigma))
        erf5 = special.erf((t-t0)/(np.sqrt(2)*sigma))
        osc_factor = np.sin(omega*t - phi) - np.cos(omega*t - phi)
        function = amp1*exp1*exp5*(1 + erf1) + amp2*exp2*exp6*(1 + erf2) - amp3*exp3*exp7*osc_factor*(1 + erf3) + amp3*exp4*exp8*(1 + erf4)*osc_factor + amp4*(1 + erf5)
        return function

    def get_fit_parameters(self, model_idxs, initial_guess_bool=False, bounds=False):
        """Fit selected model(s) to every loaded scan.

        Parameters
        ----------
        model_idxs : int or list[int]
            Model index (1, 2, or 3) for all scans, or one index per scan.
        initial_guess_bool : bool, optional
            If ``True``, pass predefined initial guesses to ``curve_fit``.
        bounds : bool, optional
            If ``True``, use predefined lower/upper bounds in the optimization.
        """

        self.params_fit = []
        self.times_fit = []
        self.tgsignals_fit = []
        self.tgsignals_sampled_fit = []
        self.taus_fit = np.array([])
        self.errors_tau_fit = np.array([])
        self.chi2_fit = np.array([])
        self.r2_fit = np.array([])


        for i in range(len(self.time_scans)):
            time = self.time_scans[i]
            tgsignal = self.tgsignal_scans[i]

            if isinstance(model_idxs, int):
                model_idx = model_idxs
            else:
                if len(model_idxs) != len(self.time_scans):
                    raise ValueError("Length of model_idxs must be equal to the number of time scans")
                else:
                    model_idx = model_idxs[i]
            
            if model_idx == 1:
                # Params: amp1, t0, tau, sigma, amp2
                model = self.model1
                initial_guess = [0.05, 0, 0.1, 0.01, 0.05] 

                lower_bounds = [0, -np.inf, 0.1, 0, 0]
                upper_bounds = [1.0, np.inf, np.inf, np.inf, np.inf]

            elif model_idx == 2:
                # Params: amp1, t0, k_tau, sigma, amp2, omega, phi, amp3, k
                model = self.model2

                omega_dom = 0.1
                phi_dom = 1.1
                initial_guess = [0.05, 0, 2, 0.01, 0.05, omega_dom, phi_dom, 0.05, 0.001] 

                lower_bounds = [0, -np.inf, 0, 0, 0, 0, 0, 0, 0]
                upper_bounds = [1.0, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf]

            elif model_idx == 3:
                # Params: amp1, amp2, amp3, amp4, t0, k1, k2, sigma, omega, phi

                omega_dom = 1.1
                phi_dom = 1.0  
                model = self.model3
                initial_guess = [0.05, 0.05, 0.05, 0.05, 0, 2, 2, 0.01, omega_dom, phi_dom] 

                lower_bounds = [0, 0, 0, 0, -np.inf, 0, 0, 0, 0, 0]
                upper_bounds = [1.0, 1.0, 1.0, 1.0, np.inf, 10, 10, np.inf, np.inf, 2*np.pi]

            else:
                raise ValueError("Invalid model index")

            # Perform the fit
            if initial_guess_bool:
                if bounds:
                    popt, pcov = curve_fit(model, time, tgsignal, p0=initial_guess, bounds=(lower_bounds, upper_bounds), maxfev=1000000)
                else:   
                    popt, pcov = curve_fit(model, time, tgsignal, p0=initial_guess, maxfev=1000000)
            else:
                if bounds:
                    popt, pcov = curve_fit(model, time, tgsignal, bounds=(lower_bounds, upper_bounds), maxfev=1000000)
                else:
                    popt, pcov = curve_fit(model, time, tgsignal, maxfev=1000000)
            perr = np.sqrt(np.diag(pcov))
            new_time = np.linspace(time.min(), time.max(), 1000)
            new_signal = model(new_time, *popt)
            new_signal_sampled = model(time, *popt)

            self.times_fit.append(new_time)
            self.tgsignals_fit.append(new_signal)
            self.tgsignals_sampled_fit.append(new_signal_sampled)

            residuals = tgsignal - new_signal_sampled
            n_params = len(popt)
            dof = len(tgsignal) - n_params        # degrees of freedom
            errors = np.full_like(tgsignal, np.std(residuals))
            chi2 = np.sum((residuals / errors)**2)
            chi2_reduced = chi2 / dof
            self.chi2_fit = np.append(self.chi2_fit, chi2_reduced)

            SS_residual = np.sum((tgsignal - new_signal_sampled)**2)
            SS_total = np.sum((tgsignal - np.mean(tgsignal))**2)
            r2 = 1 - SS_residual / SS_total
            self.r2_fit = np.append(self.r2_fit, r2)

            if model_idx == 1:
                params = {
                    "amp1": [popt[0], perr[0], "Amplitude 1 (a.u.)", "Amplitude 1"],
                    "t0": [popt[1], perr[1], r"$t_0$(ps)", "Time 0"],
                    "tau": [(1/popt[2])*1000, (perr[2]/popt[2]**2)*1000, "Decay time (fs)", "Decay time"],
                    "sigma": [popt[3], perr[3], r"$\sigma$(ps)", "Sigma"],
                    "amp2": [popt[4], perr[4], "Amplitude 2 (a.u.)", "Amplitude 2"]
                }

            elif model_idx == 2:
                params = {
                    "amp1": [popt[0], perr[0], "Amplitude 1 (a.u.)", "Amplitude 1"],
                    "t0": [popt[1], perr[1], r"$t_0$(ps)", "Time 0"],
                    "tau": [(1/popt[2])*1000, (perr[2]/popt[2]**2)*1000, "Decay time (fs)", "Decay time"],
                    "sigma": [popt[3], perr[3], r"$\sigma$(ps)", "Sigma"],
                    "amp2": [popt[4], perr[4], "Amplitude 2 (a.u.)", "Amplitude 2"],
                    "omega": [popt[5], perr[5], r"$\omega$(rad/ps)", "Omega"],
                    "phi": [popt[6], perr[6], r"$\phi$(rad)", "Phi"],
                    "amp3": [popt[7], perr[7], "Amplitude 3 (a.u.)", "Amplitude 3"],
                    "tau2": [(1/popt[8])*1000, (perr[8]/popt[8]**2)*1000, "Decay time 2 (fs)", "Decay time 2"]
                }

            elif model_idx == 3:
                params = {
                    "amp1": [popt[0], perr[0], "Amplitude 1 (a.u.)", "Amplitude 1"],
                    "amp2": [popt[1], perr[1], "Amplitude 2 (a.u.)", "Amplitude 2"],
                    "amp3": [popt[2], perr[2], "Amplitude 3 (a.u.)", "Amplitude 3"],
                    "amp4": [popt[3], perr[3], "Amplitude 4 (a.u.)", "Amplitude 4"],
                    "t0": [popt[4], perr[4], r"$t_0$(ps)", "Time 0"],
                    "tau": [(1/popt[5])*1000, (perr[5]/popt[5]**2)*1000, "Decay time (fs)", "Decay time"],
                    "tau2": [(1/popt[6])*1000, (perr[6]/popt[6]**2)*1000, "Decay time 2 (fs)", "Decay time 2"],
                    "sigma": [popt[7], perr[7], r"$\sigma$(ps)", "Sigma"],
                    "omega": [popt[8], perr[8], r"$\omega$(rad/ps)", "Omega"],
                    "phi": [popt[9], perr[9], r"$\phi$(rad)", "Phi"]
                }

            self.taus_fit = np.append(self.taus_fit, params["tau"][0])
            self.errors_tau_fit = np.append(self.errors_tau_fit, params["tau"][1])
            self.params_fit.append(params)

    def plot_phase_space(self):
        """Plot the phase-space coverage of secure, non-repeated scans."""
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

    def plot_fits(self):
        """Plot data vs fitted curves and relative error per scan."""
        if not hasattr(self, "times_fit") or not hasattr(self, "tgsignals_fit"):
            raise ValueError("Fit data not found. Run get_fit_parameters() first.")
        if not hasattr(self, "tgsignals_sampled_fit"):
            raise ValueError("Sampled fit not found. Run get_fit_parameters() first.")

        n_scans = len(self.scans)
        n_cols = int(np.ceil(np.sqrt(n_scans)))
        n_rows = int(np.ceil(n_scans / n_cols))

        fig = plt.figure(figsize=(4 * n_cols, 4 * n_rows))
        outer_gs = fig.add_gridspec(n_rows, n_cols, hspace=0.45, wspace=0.35)

        eps = 1e-12
        for i in range(n_scans):
            row, col = divmod(i, n_cols)
            inner = outer_gs[row, col].subgridspec(
                2, 1, height_ratios=[3, 1], hspace=0.08
            )
            ax_main = fig.add_subplot(inner[0, 0])
            ax_res = fig.add_subplot(inner[1, 0], sharex=ax_main)

            time = self.time_scans[i]
            y_data = self.tgsignal_scans[i]
            y_fit_s = self.tgsignals_sampled_fit[i]

            ax_main.scatter(time, y_data, color="blue", s=12)
            ax_main.plot(
                self.times_fit[i],
                self.tgsignals_fit[i],
                color="red",
                linestyle="-",
                linewidth=1.5,
            )
            ax_main.set_title(self.scans[i] + " - E: " + f"{self.energies[i]:.1f} eV, I: {self.intensities[i]:.1f} $\\mu$J" + "\n" + rf"$\chi^2$/dof = {self.chi2_fit[i]:.4f}, $R^2$ = {self.r2_fit[i]:.4f}")
            ax_main.set_ylabel("TG Signal")
            ax_main.grid(linestyle="--", alpha=0.5)
            ax_main.tick_params(labelbottom=False)

            rel_pct = y_data - y_fit_s
        
            ax_res.plot(time, rel_pct, color="green", linewidth=1.0)
            ax_res.set_xlabel("Time (ps)")
            ax_res.set_ylabel(r"Rel.\ Diff.\ (a.u.)", fontsize=8)
            ax_res.tick_params(axis="both", labelsize=8)
            ax_res.grid(linestyle="--", alpha=0.5)

        fig.tight_layout()
        plt.show()

    def plot_params_vs_energy(self, param_name, errors_bool=False):
        """Plot fitted decay times vs energy with absorbance references."""
        
        data_Co2plus = np.genfromtxt("/Users/manuelfernandosanchezalarcon/Desktop/Trieste_Project/Transient_Grating/transient_grating_project/external_files/Co2plus_absorbance.txt")
        data_Co3plus = np.genfromtxt("/Users/manuelfernandosanchezalarcon/Desktop/Trieste_Project/Transient_Grating/transient_grating_project/external_files/Co3plus_absorbance.txt")
        data_Co3O4 = np.genfromtxt("/Users/manuelfernandosanchezalarcon/Desktop/Trieste_Project/Transient_Grating/transient_grating_project/external_files/Co3O4_absorbance.txt")
        
        energy_Co2, Co2_absorbance = data_Co2plus[:,0]-0.4, data_Co2plus[:,1]
        energy_Co3, Co3_absorbance = data_Co3plus[:,0]-0.4, data_Co3plus[:,1]
        energy_Co3O4, Co3O4_absorbance = data_Co3O4[:,0]-0.4, data_Co3O4[:,1]
        
        fig, ax1 = plt.subplots(figsize=(8, 5))
        
        # --- Left Y axis: Absorbance ---
        ax1.fill_between(energy_Co2, Co2_absorbance, alpha=0.2, color='gray', label=r'Co$^{2+}$ absorbance')
        ax1.fill_between(energy_Co3, Co3_absorbance, alpha=0.2, color='green', label=r'Co$^{3+}$ absorbance')
        ax1.plot(energy_Co3O4, Co3O4_absorbance, color='black', label=r'Co$_3$O$_4$ absorbance')
        
        ax1.set_xlabel("Energy (eV)")
        ax1.set_ylabel("Absorbance (a.u.)")
        ax1.set_ylim(0, 0.35)
        ax1.set_xlim(energy_Co2.min(), energy_Co2.max())
        
        # --- Right Y axis: Decay Time ---
        ax2 = ax1.twinx()

        param = np.array([self.params_fit[i][param_name][0] for i in range(len(self.params_fit))])
        label_param = self.params_fit[0][param_name][2]
        label_title = self.params_fit[0][param_name][3]
        if errors_bool:
            errors = np.array([self.params_fit[i][param_name][1] for i in range(len(self.params_fit))])
            ax2.errorbar(self.energies, param, yerr=errors, fmt='o-', color='blue', ecolor='blue', capsize=5)
        else:
            errors = 0
            ax2.plot(self.energies, param, 'o-', color='blue')
            
            
        ax2.set_ylabel(label_param, color='blue')
        ax2.tick_params(axis='y', colors='blue')
        ax2.spines['right'].set_color('blue')  # Right spine also in blue

        # --- Legend (combine both axes) ---
        lines1, labels1 = ax1.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax1.legend(lines1 + lines2, labels1 + labels2, loc='best')
        
        plt.title(f"{label_title} vs Energy")
        plt.tight_layout()
        plt.show()

    def plot_stacked_signals(self, limits_time=None, limits_signal=None, ylog_scale = False, plot_ind=None, data_over_fit=False):
        """Plot stacked TG signals and optionally stacked fits.

        Parameters
        ----------
        limits_time : tuple[float, float] or None, optional
            Time-axis limits in ps.
        limits_signal : tuple[float, float] or None, optional
            Signal-axis limits.
        ylog_scale : bool, optional
            If ``True``, display y-axis in logarithmic scale.
        plot_ind : list[int] or str or None, optional
            Subset of scan indices to plot, or ``"all"`` for every scan.
        data_over_fit : bool, optional
            If ``True``, overlay sampled data in the fit panel.
        """
        fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
        ax_data, ax_fit = axes
        for i, (time, tgsignal) in enumerate(zip(self.time_scans, self.tgsignal_scans)):
            if plot_ind is None:
                ax_data.plot(time, tgsignal, label=f"S{self.scans[i][4:]}, E={self.energies[i]:.1f} eV, I={self.intensities[i]:.1f} $\\mu$J")
                if hasattr(self, "times_fit") and hasattr(self, "tgsignals_fit"):
                    ax_fit.plot(
                        self.times_fit[i],
                        self.tgsignals_fit[i],
                        label=f"S{self.scans[i][4:]}, E={self.energies[i]:.1f} eV, I={self.intensities[i]:.1f} $\\mu$J",
                    )
            else:
                if plot_ind == "all":
                    ax_data.plot(time, tgsignal, label=f"S{self.scans[i][4:]}, E={self.energies[i]:.1f} eV, I={self.intensities[i]:.1f} $\\mu$J")
                    if hasattr(self, "times_fit") and hasattr(self, "tgsignals_fit"):
                        ax_fit.plot(
                            self.times_fit[i],
                            self.tgsignals_fit[i],
                            label=f"S{self.scans[i][4:]}, E={self.energies[i]:.1f} eV, I={self.intensities[i]:.1f} $\\mu$J",
                        )
                else:
                    if i in plot_ind:
                        ax_data.plot(time, tgsignal, label=f"S{self.scans[i][4:]}, E={self.energies[i]:.1f} eV, I={self.intensities[i]:.1f} $\\mu$J")
                        if hasattr(self, "times_fit") and hasattr(self, "tgsignals_fit"):
                            ax_fit.plot(
                                self.times_fit[i],
                                self.tgsignals_fit[i],
                                label=f"S{self.scans[i][4:]}, E={self.energies[i]:.1f} eV, I={self.intensities[i]:.1f} $\\mu$J",
                            )
                            if data_over_fit:
                                ax_fit.scatter(time, tgsignal, label=f"S{self.scans[i][4:]}, Experimental data")
                                ax_fit.plot(self.times_fit[i], np.exp(-(self.times_fit[i] - self.params_fit[i]["t0"])/self.params_fit[i]["tau"]) + np.mean(tgsignal[time > 2.5]), color="red", linestyle="--")

        ax_data.set_xlabel("Time (ps)")
        ax_data.set_ylabel("TG Signal")
        ax_data.set_title("Stacked Signals (Data)")
        ax_fit.set_xlabel("Time (ps)")
        ax_fit.set_ylabel("TG Signal")
        ax_fit.set_title("Stacked Signals (Fit)")

        if limits_time is not None:
            ax_data.set_xlim(limits_time[0], limits_time[1])
            ax_fit.set_xlim(limits_time[0], limits_time[1])
        if limits_signal is not None:
            ax_data.set_ylim(limits_signal[0], limits_signal[1])
            ax_fit.set_ylim(limits_signal[0], limits_signal[1])
        ax_data.grid(linestyle='--', alpha=0.5)
        ax_fit.grid(linestyle='--', alpha=0.5)
        handles, labels = ax_data.get_legend_handles_labels()
        fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, -0.02), ncol=3, fontsize=7)

        if ylog_scale:
            ax_data.set_yscale("log")
            ax_fit.set_yscale("log")
        fig.tight_layout(rect=[0, 0.08, 1, 1])
        plt.show()

    def plot_params_all_models(self, models_config, param_name, errors_bool=False, y_limits=None):
        """Compare fitted decay times for multiple model configurations.

        Parameters
        ----------
        models_config : list[dict]
            List of fitting configurations passed to ``get_fit_parameters``.
        errors_bool : bool, optional
            If ``True``, include error bars for fitted ``tau`` values.
        """

        plt.figure(figsize=(8, 5))
        for model_config in models_config:
            self.get_fit_parameters(model_idxs=model_config["model_idxs"], initial_guess_bool=model_config["initial_guess_bool"], bounds=model_config["bounds"])
            param = np.array([self.params_fit[i][param_name][0] for i in range(len(self.params_fit))])
            if errors_bool:
                errors = np.array([self.params_fit[i][param_name][1] for i in range(len(self.params_fit))])
                plt.errorbar(self.energies, param, yerr=errors, fmt='o-', capsize=5, label=model_config["label_model"])
            else:
                errors = 0
                plt.plot(self.energies, param, 'o-', label=model_config["label_model"])

        plt.grid(linestyle='--', alpha=0.5)
        if y_limits is not None:
            plt.ylim(y_limits[0], y_limits[1])
        plt.legend()
        plt.show()
            



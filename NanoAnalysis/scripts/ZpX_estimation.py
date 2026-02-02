import numpy as np

import matplotlib.pyplot as plt
import mplhep as hep

import ROOT

class ZpX:
    def __init__(self, zpx_procs = ["DY","TT","WZ"], max_sip = 4):
        self.zpx_procs = zpx_procs

        self.max_sip = max_sip
        self.max_idx = lambda edges: int(np.argwhere(edges==self.max_sip)) + 1

        self.fstates = ["fs_4l", "fs_2x2e", "fs_2x2mu"]

        self.fstate_map = dict(
            fs_4l    = r"$4l$",
            fs_2x2e  = r"$2X2e$",
            fs_2x2mu = r"$2X2mu$"
        )

        self.plot_info = dict(
            N_ZpX_MidMass = dict(
                ylabel = r"$N_{Z+X}^{SR}$",
            ),
            N_ZPP_SS = dict(
                ylabel = r"$N_{Z+X}^{CR}$",
            ),
            r_OS_SS_MidMass = dict(
                ylabel = r"$r_{OS/SS}$",
            ),
            r_OS_SS_LowMass = dict(
                ylabel = r"$r_{OS/SS}$",
            ),
            N_ZpX_LowMass = dict(
                ylabel = r"$N_{Z+X}^{SR}$",
            ),
        )

        self.lumi_2022 = 34.6532e3
        self.lumi_2023 = 27.245e3

    def get_count(self, hist, errs, var):
        """Return histogram count and error for variable var.
        If checking sip, count+err for sip<=4."""
        counts, edges = hist

        # Array index corresponding to SIP=4
        max_idx = len(counts) if var != "Lepton_sip3d_Z2" else self.max_idx(edges)
        count = counts[:max_idx].sum()
        err   = np.sqrt(np.square(errs[:max_idx]).sum())

        return count, err          

    def nZPP(self, reg, fs, all_hists, all_errors, var = "mass", max_edges = 4):
        """Get count+err of Data - MC_{!ZpX}."""
        data_hist = all_hists[var][reg][fs]["Data"]
        data_err  = all_errors[var][reg][fs]["Data"]
        
        count, n_err = self.get_count(data_hist, data_err, var)
        
        tot_errs = [n_err]
        for proc in all_hists[var][reg][fs]["MC"].keys():

            if proc in self.zpx_procs:
                continue
            
            mc_hist = all_hists[var][reg][fs]["MC"][proc]
            mc_err  = all_errors[var][reg][fs]["MC"][proc]

            mc_count, mc_err = self.get_count(mc_hist, mc_err, var)

            count -= mc_count
            tot_errs.append(mc_err)

        tot_err = np.sqrt(np.square(tot_errs).sum())
        
        return count, tot_err

    def get_nZPPSS(self, fs, all_hists, all_errors):
        """Get _event_ count for N_{Z+X}^{CR}."""
        lep_count, lep_err = self.nZPP("HighMassSSSIP", fs, all_hists, all_errors, var = "Lepton_sip3d_Z2")
        return lep_count/2, lep_err/2

    def get_r(self, fs, mass_reg, all_hists, all_errors):
        """Calculate os_ss transfer function and associated error
        in a given mass region."""

        os_reg = mass_reg + "OSSIP"
        ss_reg = mass_reg + "SSSIP"

        os_count, os_err = self.nZPP(os_reg, fs, all_hists, all_errors)
        ss_count, ss_err = self.nZPP(ss_reg, fs, all_hists, all_errors)

        ratio = os_count/ss_count
        
        ratio_err = ratio*np.sqrt((os_err/os_count)**2 + (ss_err/ss_count)**2)
        if ratio < 0:
            ratio = 0
            ratio_err *= -1 # If ratio < 0, ratio_err would be negative without this

        return ratio, ratio_err

    def get_zpx(self, all_hists, all_errors):
        """Calcualte and save ZpX^SR estimate using the
        Mid and Low mass transfer functions."""

        zpx_info = dict(
            N_ZPP_SS        = {},
            r_OS_SS_MidMass = {},
            r_OS_SS_LowMass = {},
            N_ZpX_MidMass   = {},
            N_ZpX_LowMass   = {}
        )

        for fs in self.fstates:
            # N_ZpX^CR
            n, n_err = self.get_nZPPSS(fs, all_hists, all_errors)

            zpx_info["N_ZPP_SS"][fs] = (n, n_err)

            # Transfer function r_os_ss for Mid (Low) Mass
            r_mm, r_mm_err = self.get_r(fs, "MidMass", all_hists, all_errors)
            r_lm, r_lm_err = self.get_r(fs, "LowMass", all_hists, all_errors)

            zpx_info["r_OS_SS_MidMass"][fs] = (r_mm, r_mm_err)
            zpx_info["r_OS_SS_LowMass"][fs] = (r_lm, r_lm_err)

            # N_ZpX^SR = N_ZpX^CR x r_os_ss
            zpx_mm = n*r_mm
            zpx_mm_err = np.sqrt(n_err**2 + r_mm_err**2)

            zpx_info["N_ZpX_MidMass"][fs] = (zpx_mm, zpx_mm_err)

            zpx_lm = n*r_lm
            zpx_lm_err = np.sqrt(n_err**2 + r_lm_err**2)

            zpx_info["N_ZpX_LowMass"][fs] = (zpx_lm, zpx_lm_err)
        
        return zpx_info

    def plot_zpx(self, zpx_info, year, step="N_ZpX_MidMass", tag=""):
        """Quick plotter for ZpX estimates for a given year. By default
        plots 4l state on the left, and 2x2e and 2x2mu final states on the right,
        seperated by a dotted black line. Steps correspond to self.plot_info.keys()."""
        
        hep.style.use("CMS")
        
        zpx_dict = zpx_info[step]
        
        vals = dict(
            Old = [zpx_dict[fs][0] for fs in zpx_dict.keys()]
        )
        errs = dict(
            Old = [zpx_dict[fs][1] for fs in zpx_dict.keys()]
        )
        
        final_states = [self.fstate_map[key] for key in zpx_dict.keys()]

        x = np.arange(len(final_states))
        group_width = 0.5
        offset_step = group_width/len(vals)

        fig, ax = plt.subplots(layout='constrained')

        for i, (key, val) in enumerate(vals.items()):
            offset = x - (group_width - offset_step)/2 + i*offset_step
            ax.errorbar(
                offset,
                vals[key],
                yerr=errs[key],
                fmt="o",
            )

        plt.axvline(x=0.5, linestyle="--", color="black")

        ax.set_ylabel(self.plot_info[step]["ylabel"], fontsize=40)

        ax.set_xticks(x)
        ax.set_xticklabels(final_states)
        ax.tick_params(axis='both', labelsize=40)

        filename = f"{step}_{year}{tag}.pdf"
        fig.savefig(filename, dpi=600, format="pdf")

    def plot_zpx_years(self, zpx_22, zpx_23, step="N_ZpX_MidMass", tag=""):
        """Quick plotter for ZpX estimates for 2022 and 2023, shown /fb^-1.
        Otherwise identical to self.plot_zpx()."""

        hep.style.use("CMS")
        
        dict_22 = zpx_22[step]
        dict_23 = zpx_23[step]

        vals = {
            "2022": [dict_22[fs][0] for fs in dict_22.keys()],
            "2023": [dict_23[fs][0] for fs in dict_22.keys()],
        }

        errs = {
            "2022": [dict_22[fs][1] for fs in dict_22.keys()],
            "2023": [dict_23[fs][1] for fs in dict_22.keys()],
        }

        norm_counts = {
            "2022": [cnt/(self.lumi_2022*1e-3) for cnt in vals["2022"]],
            "2023": [cnt/(self.lumi_2023*1e-3) for cnt in vals["2023"]]
        }

        # Err. in Lumi for 2022 (2023) is 1.4% (1.3%)
        norm_errs = {
            "2022": [abs(nm_cnt)*np.sqrt((err/cnt)**2 + (0.014)**2) for nm_cnt, err, cnt in zip(norm_counts["2022"], errs["2022"], vals["2022"])],
            "2023": [abs(nm_cnt)*np.sqrt((err/cnt)**2 + (0.013)**2) for nm_cnt, err, cnt in zip(norm_counts["2023"], errs["2023"], vals["2023"])]
        }
        
        final_states = [self.fstate_map[key] for key in dict_22.keys()]

        x = np.arange(len(final_states))
        group_width = 0.5
        offset_step = group_width/len(vals)

        fig, ax = plt.subplots(layout='constrained')
        for i, (key, val) in enumerate(vals.items()):
            offset = x - (group_width - offset_step)/2 + i*offset_step
            ax.errorbar(
                offset,
                norm_counts[key],
                yerr=norm_errs[key],
                label=key,
                fmt="o",
            )

        ax.legend()

        plt.axvline(x=0.5, linestyle="--", color="black")

        ylabel = r"$N_{Z+X}^{SR} / fb^{-1}$"

        ax.set_ylabel(ylabel, fontsize=40)

        ax.set_xticks(x)
        ax.set_xticklabels(final_states)
        ax.tick_params(axis='both', labelsize=40)

        filename = f"{step}_2022_2023{tag}.pdf"
        fig.savefig(filename, dpi=600, format="pdf")      
import numpy as np
import matplotlib.pyplot as plt

class ZpX:
    def __init__(self, zpx_procs = ["DY","TT"], max_sip = 4):
        self.zpx_procs = zpx_procs

        self.max_sip = max_sip
        self.max_idx = lambda edges: int(np.argwhere(edges==self.max_sip)) + 1

        self.plot_info = dict(
            N_ZPP_SS        = dict(
                y_label = r"$N_{ZPP_{SS}}$",
                title   = "Passing Same Sign High Mass",
            ),
            r_OS_SS_MidMass = dict(
                y_label = r"$r_{OS/SS}$",
                title   = "Mid Mass Transfer Function"
            ),
            r_OS_SS_LowMass = dict(
                y_label = r"$r_{OS/SS}$",
                title   = "Low Mass Transfer Function"
            ),
            N_ZpX_MidMass   = dict(
                y_label = r"$N_{ZpX}$",
                title   = "ZpX Estimate (Mid Mass Transfer Function)"
            ),
            N_ZpX_LowMass   = dict(
                y_label = r"$N_{ZpX}$",
                title   = "ZpX Estimate (Low Mass Transfer Function)"
            )
        )

    def get_count(self, hist, err, var):
        counts, edges = hist
        max_idx = len(counts) if var != "Lepton_sip3d_Z2" else self.max_idx(edges)
        count = counts[:max_idx].sum()
        err   = np.sqrt(np.square(err[:max_idx]).sum())

        return count, err          

    def nZPP(self, reg, fs, all_hists, all_errors, var = "mass", max_edges = 4):
        data_hist = all_hists[reg][var][fs]["Data"]
        data_err  = all_errors[reg][var][fs]["Data"]
        
        n, n_err = self.get_count(data_hist, data_err, var)
        
        tot_errs = [n_err]
        for proc in all_hists[reg][var][fs]["MC"].keys():
            
            if proc in self.zpx_procs:
                continue
            
            mc_hist = all_hists[reg][var][fs]["MC"][proc]
            mc_err  = all_errors[reg][var][fs]["MC"][proc]

            mc_count, mc_err = self.get_count(mc_hist, mc_err, var)

            n -= mc_count
            tot_errs.append(mc_err)

        tot_err = np.sqrt(np.square(tot_errs).sum())
        
        return n, tot_err

    def get_nZPPSS(self, fs, all_hists, all_errors):
        lep_count, lep_err = self.nZPP("SS_NoSIP_HighMass", fs, all_hists, all_errors, var = "Lepton_sip3d_Z2")
        return lep_count/2, lep_err/2

    def get_r(self, fs, mass_reg, all_hists, all_errors):
        os_reg = "OS_NoSIP_" + mass_reg
        ss_reg = "SS_NoSIP_" + mass_reg

        os_count, os_err = self.nZPP(os_reg, fs, all_hists, all_errors)
        ss_count, ss_err = self.nZPP(ss_reg, fs, all_hists, all_errors)

        ratio = os_count/ss_count
        ratio_err = ratio*np.sqrt((os_err/os_count)**2 + (ss_err/ss_count)**2)

        return ratio, ratio_err

    def get_zpx(self, all_hists, all_errors, fstates):

        zpx_info = dict(
            N_ZPP_SS        = {},
            r_OS_SS_MidMass = {},
            r_OS_SS_LowMass = {},
            N_ZpX_MidMass   = {},
            N_ZpX_LowMass   = {}
        )

        for fs in fstates:
            n, n_err = self.get_nZPPSS(fs, all_hists, all_errors)

            zpx_info["N_ZPP_SS"][fs] = (n, n_err)

            r_mm, r_mm_err = self.get_r(fs, "MidMass", all_hists, all_errors)
            r_lm, r_lm_err = self.get_r(fs, "LowMass", all_hists, all_errors)

            zpx_info["r_OS_SS_MidMass"][fs] = (r_mm, r_mm_err)
            zpx_info["r_OS_SS_LowMass"][fs] = (r_lm, r_lm_err)
            
            zpx_mm = n*r_mm
            zpx_mm_err = np.sqrt(n_err**2 + r_mm_err**2)

            zpx_info["N_ZpX_MidMass"][fs] = (zpx_mm, zpx_mm_err)

            zpx_lm = n*r_lm
            zpx_lm_err = np.sqrt(n_err**2 + r_lm_err**2)

            zpx_info["N_ZpX_LowMass"][fs] = (zpx_lm, zpx_lm_err)
        
        return zpx_info

    def plot_zpx(self, zpx_info, step):
        fstates = zpx_info[step].keys()

        counts = [zpx_info[step][fs][0] for fs in fstates]
        errs   = [zpx_info[step][fs][1] for fs in fstates]

        y_label = self.plot_info[step]["y_label"]
        title   = self.plot_info[step]["title"]

        fig, ax = plt.subplots()
        ax.errorbar(fstates, counts, yerr=errs, linestyle="None", marker = "o", color="black")
        ax.set_ylabel(y_label, rotation="horizontal")

        ax.set_title(title)

        fig.savefig(step+".png")
        



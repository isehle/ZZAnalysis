import numpy as np
import matplotlib.pyplot as plt

import ROOT

# Note: originally done considering ZpX = DY+TT,
# should really include also WZto3LNu
class ZpX:
    def __init__(self, zpx_procs = ["DY","TT","WZto3LNu"], max_sip = 4):
        self.zpx_procs = zpx_procs

        self.max_sip = max_sip
        self.max_idx = lambda edges: int(np.argwhere(edges==self.max_sip)) + 1

        self.plot_info = dict(
            N_ZPP_SS        = dict(
                # y_label = r"$N_{ZPP_{SS}}/{fb^{-1}}$",
                # title   = r"$N_{ZPP_{SS}}/{fb^{-1}}$"
                y_label = r"$N_{ZPP_{SS}}$",
                title   = r"$N_{ZPP_{SS}}$"
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
        data_hist = all_hists[var][reg][fs]["Data"]
        data_err  = all_errors[var][reg][fs]["Data"]
        
        n, n_err = self.get_count(data_hist, data_err, var)
        
        tot_errs = [n_err]
        for proc in all_hists[var][reg][fs]["MC"].keys():

            if proc in self.zpx_procs:
                continue
            
            mc_hist = all_hists[var][reg][fs]["MC"][proc]
            mc_err  = all_errors[var][reg][fs]["MC"][proc]

            mc_count, mc_err = self.get_count(mc_hist, mc_err, var)

            n -= mc_count
            tot_errs.append(mc_err)

        tot_err = np.sqrt(np.square(tot_errs).sum())
        
        return n, tot_err

    def get_nZPPSS(self, fs, all_hists, all_errors):
        #lep_count, lep_err = self.nZPP("SS_NoSIP_HighMass", fs, all_hists, all_errors, var = "Lepton_sip3d_Z2")
        lep_count, lep_err = self.nZPP("HighMassSSSIP", fs, all_hists, all_errors, var = "Lepton_sip3d_Z2")
        return lep_count/2, lep_err/2

    def get_r(self, fs, mass_reg, all_hists, all_errors):
        # os_reg = "OS_NoSIP_" + mass_reg
        # ss_reg = "SS_NoSIP_" + mass_reg

        os_reg = mass_reg + "OSSIP"
        ss_reg = mass_reg + "SSSIP"

        os_count, os_err = self.nZPP(os_reg, fs, all_hists, all_errors)
        ss_count, ss_err = self.nZPP(ss_reg, fs, all_hists, all_errors)

        if os_count < 0: os_count = 1e-10
        if ss_count < 0: ss_count = 1e-10

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

    def get_yields(self, all_hists, all_errors, fstates):
        sip_less4_count = lambda sipHist: sipHist[0][1:5].sum() #1st bin is underflow
        get_err = lambda errArr: np.sqrt(np.sum(np.square(errArr[1:5])))

        # sip3d_z2_highMass_SS = all_hists["SS_NoSIP_HighMass"]["Lepton_sip3d_Z2"]
        # err_arrs = all_errors["SS_NoSIP_HighMass"]["Lepton_sip3d_Z2"]

        sip3d_z2_highMass_SS = all_hists["Lepton_sip3d_Z2"]["HighMassSSSIP"]
        err_arrs = all_errors["Lepton_sip3d_Z2"]["HighMassSSSIP"]

        counts = {}
        for fs in ["fs_4e", "fs_4mu", "fs_2e2mu", "fs_2mu2e"]:
            data_z2_leps = sip3d_z2_highMass_SS[fs]["Data"]["Data"]
            data_err_arr = err_arrs[fs]["Data"]["Data"]

            # Arbitrary normalization for Z+X which will be fit by combine
            # Need to divide by 2 to get event counts since these leptons from Z2-->ll, and Data must be integer (CHECK THAT THIS IS OKAY)
            counts[fs] = {"data_obs": (round(sip_less4_count(data_z2_leps)/2), get_err(data_err_arr)/2),
                          "ZpX":  (1., 0.)}

            mc_hists = sip3d_z2_highMass_SS[fs]["MC"]
            mc_errs  = err_arrs[fs]["MC"]

            for proc in ["ZZ_NLO", "ggZZ", "H", "VVV"]:
                count, err = sip_less4_count(mc_hists[proc])/2, get_err(mc_errs[proc])

                counts[fs][proc] = (count, err)

        return counts

    def plot_zpx(self, zpx_info, step, *args):
        r_os_ss_y_lim  = (0, 7)
        n_zpp_ss_y_lim = (-2, 4)
        fstates = zpx_info[step].keys()

        if "N_ZpX" in step:
            categories = fstates
            counts = dict(
                MidMass = [zpx_info["N_ZpX_MidMass"][fs][0] for fs in fstates],
                LowMass = [zpx_info["N_ZpX_LowMass"][fs][0] for fs in fstates]
            )
            errs = dict(
                MidMass = [zpx_info["N_ZpX_MidMass"][fs][1] for fs in fstates],
                LowMass = [zpx_info["N_ZpX_LowMass"][fs][1] for fs in fstates]
            )

            x = np.arange(len(fstates))
            group_width = 0.5
            offset_step = group_width/len(counts)

            fig, ax = plt.subplots(layout='constrained')
            for i, (key, val) in enumerate(counts.items()):
                offset = x - (group_width - offset_step)/2 + i*offset_step
                ax.errorbar(
                    offset,
                    counts[key],
                    yerr=errs[key],
                    fmt="o",
                    label=key
                )
            
            ax.legend()

            ax.set_xticks(x)
            ax.set_xticklabels(fstates)

            # ax.set_ylim()

            ax.set_ylabel(r"$N_{Z+X}$", rotation="horizontal")

            #title = "N_ZpX", "N_ZpX"
            title = "N_ZpX"
            for arg in args:
                title += " {}".format(arg)
            
            ax.set_title(title)
        
        else:
            counts = [zpx_info[step][fs][0] for fs in fstates]
            errs   = [zpx_info[step][fs][1] for fs in fstates]

            # if step == "N_ZPP_SS":
            #     lumi = 34.6521 if int(args[0]) == 2022 else 27.245

            #     norm_counts = [cnt/lumi for cnt in counts]
            #     norm_errs   = [abs(nm_cnt)*np.sqrt((err/cnt)**2 + (0.015)**2) for nm_cnt, err, cnt in zip(norm_counts, errs, counts)]

            #     counts, errs = norm_counts, norm_errs

            y_label = self.plot_info[step]["y_label"]
            title   = self.plot_info[step]["title"]

            fig, ax = plt.subplots()
            ax.errorbar(fstates, counts, yerr=errs, linestyle="None", marker = "o", color="black")
            ax.set_ylabel(y_label)

            # if "r_OS" in step:
            #     ax.set_ylim(*r_os_ss_y_lim)
            # else:
            #     ax.set_ylim(*n_zpp_ss_y_lim)
            #     #ax.set_ylim(-0.1, 0.1)

            ax.set_title(title)

        return fig
            
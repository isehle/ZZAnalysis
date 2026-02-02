## Histogramming and Plotting Scripts

### General — histManager

These are a collection of Python scripts and YAML configuration files to process the ntuples output from CJLST. The histogramming is run in `histWriter.py` and the plotting in `histPlotter.py`, with file management handled by `fileHandler.py`. Systematics (up/down variations) are handled through `histSystematics.py`. Typically any call would be handled through `histManager.py`. Currently, the script accepts `--year` and `--era` arguments for the data taking eras in 2022 and 2023, though it should be simple to extend to further years. The script can be run in four different "modes" (controlled by the `--mode` command line argument):

- `plot_hists` : Will first look for histograms stored in the base directory in `files_cfg.yaml` under `output/histograms/{year}/{era}` under the name `hists{tag}.root` where `--tag` is an additional argument which can be passed to `histManager.py` to add this name changing tag. If the file is not found, it will first write the hists, store them in this directory, and then make the plots and store them in the directory under `output/plots` in `files_cfg.yaml`

- `plot_zpx` : Will create plots for the estimation of the Z+X background, using the SIP method. By default it uses data stored in the `.json` files whose path is given in `gen_cfg.yaml` under `zpx/year_{year}`. If `--year -1` is given, it will plot the 2022 and 2023 years side by side on the same graph.

- `combine_eras` : Adds the histograms in the files in  `gen_cfg.yaml` under `combine_eras/{year_tag}`.

### General Configuration
- The writing configuration can be used to set which region and final states are run on, and for technical reasons also links to which histogram configuration file is used.
- `zpx` points to the `json` files which store the data for the Z+X estimation.
- `combine_eras` lists which two input files (of histograms) should be combined.

### File Configuration
- The paths to the ntuples from CJLST should be listed in `year_{year}/{Era}`, where `MC` has subdirectories for each process. The names of these subdirectories correspond to the values under `MC_Procs/{proc_key}` where `proc_key` is the name you wish the histogram title to be. If `proc_key` corresponds to several different processes, they will be combined into one histogram. The corresponding `lumi` should also be stored here. By default the output of CJLST names each file `ZZ4lAnalysis.root`, which is saved in the config file as `mc_file_name`. All paths are assumed to be relative to the `eos_base`.
- Output histogram and output file directories listed in `output`.

### Histogram Configuration
- Stored in `hist_cfg.yaml`
- First key should always be the basename of the variable (`ZZCand_<basename>` for CJLST).
    - Exception for lepton sampels which begin with `Lepton` in the CJLST output. This should be kept, and the following additional options can be appended: `[_Z1, _Z2, _Z1l1, _Z1l2, _Z2l1, _Z2l2]`. By default we store `_Z2` since this is needed for the Z+X estimation.
- Optional: can set different binnings under `SR`, `MidMass` and `LowMass` regions. The high mass same sign region of the Z+X SIP method for polarized ZZ uses the same binning as under `SR`.
- Optional: Can include `systematics` for specific variables, with the variations calculated for specific procsses governed by `vars` and `procs`, respectively.

### Systematics Configuration
- Flag to run or not
- Which regions to calculate systematics
- Which final states to run on (must be subset or set of `fstates` in `gen_cfg.yaml`)
- Which variables to calculate variations for
- Which variations to calculate for which processes

### Plotting Configuration
- First several sub-dictionaries control plot aesthetics. The `hatch_style` is used for MC statistical error.
- `extra` allows for a few more options:
    - `path` : if given, will plot histograms in this file
    - `rebin` : Boolean, if True will rebin the variables that are given with the new `bins` that are given.
    - `group` : Will combine specific histgorams into one for the plot. The title will be given by the key, with the histograms to be added given by those in the corresponding list.

### Combine Hists
The one standalone script (not managed by `histManager.py`) is `write_combine_hists.py`. This was simply a lack of time and can and should be implemented into the more global framework. The file reformats an input histogram file and outputs a new file which stores histograms designed to be used by Combine. This requires A `.json` file which stores the expected yields of `ZZ NLO` (unpolarized) nominally, and the expected ratios for each systematic variation. 

    Example:

        "2022": {
        "fs_2x2e": {
            "nom": 304.98520074802946,
            "LHEScaleWeightUp": 1.0179943641421512,
            "LHEScaleWeightDown": 0.9836792958822939,
            "LHEPdfWeightUp": 1.017567420756459,
            "LHEPdfWeightDown": 0.9824325792435409,
            "puWeightUp": 0.9864548913254596,
            "puWeightDown": 1.0132319813679047,
            "lepIDRecUp": 1.0451611,
            "lepIDRecDown": 0.95483891
        },
        }
It additionally requires `json` files for the expected Z+X estimation.
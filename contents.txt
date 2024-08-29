**MATLAB codes for regressing model for oxalic acid yields**
Kinetic model for numerical simulations of oxalic acid yields during benzyl alcohol oxidation, and determining benzyl alcohol -> hydroxyoxopentenal rates by regression. Includes equations for acid-base equilibrium.

path:
"/matlab_code_regression_analysis/regress_HOP_r_Ox"

contents:
acid_base.m
model_HOP.m
run_reg.m
rxn_network.m

---------------
**generate thermochemistry**
Codes used to generate thermochemistry from DFT output implemented in python with example. Contained in the "python_script_thermochemical_calculations" folder.

jupyter notebook to read DFT outputs and calculate thermochemistry
	thermochemistry_analysis.ipynb

SMILES string reprentation of reaction set for sample benzene fragmentaiton network
	sample_calculation/benzene_frag_rxns.csv

python scripts used to analyze energy and thermochemistry from DFT outputs called on by "thermochemistry_analysis.ipynb"
	custom_functions/thermo_parameters_py

---------------
**reaction sequences in SMILES representations**
plots and tables for reaction thermochemistry were generated from DFT outputs using SMILES representations of reaction steps. These SMILES representations are contained in .csv files in the "/reaction_networks_SMILES" folder

benzyl alcohol ortho fragmenation (Fig. 3b and Table S9):
	BA_ortho_fragmentation.csv

5hydroxy4oxo2pentenal fragmentation initiated by superoxide (Scheme 4, Scheme S3, Table S10, Eq. S25-26)
	5hydroxy4oxo2pentenal_frag.csv

benzene fragmentation (Fig. S1 and Table S8)
	benzene_frag.csv

ion hydration (table S15)
	ion_association.csv

fragmented products activation (Table S4)
	frag_activation.csv

benzyl alcohol initial reactions (Table S6-S7, Schemes 1-2)
	BA_initial_rxns.csv

----------------
**DFT output files**
All DFT output files and post-processing data needed to calculate kinetic and thermodynamic parameters are contained in "/DFT_datasets".

the species and TS referred to by "BA_ortho_fragmentation.csv" are listed in:
	species_benzyl_alcohol.csv
	TS_benzyl_alcohol.csv

by "BA_initial_rxns.csv" in:
	species_benzyl_alcohol.csv
	TS_benzyl_alcohol.csv

by "benzene_frag.csv" in:
	species_benzene.csv
	TS_benzene.csv

by "frag_activate.csv" in:
	species_fragments_activation.csv
	TS_fragments_activation.csv

by "5hydroxy4oxo2pentenal_frag.csv" in:
	species_pentenal_frag.csv
	TS_pentenal_frag.csv
	
by "ion_association.csv" in:
	species_ion_association.csv
	TS_ion_association.csv

additional DFT calculations to support results included in "misc_DFT_outputs" directory

----------------
**molecular volumes**
volumes calculated from DFT structures reported in "volume_outputs" directory and used for thermochemical and kinetics calculations







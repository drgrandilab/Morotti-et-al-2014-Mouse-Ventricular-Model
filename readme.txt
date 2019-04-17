Matlab code for Morotti et al. mouse ventricular model.

______________________________________________________________________________________
Contents:

readme.txt			this file

morotti_et_al_masterCompute.m	loads initial conditions and runs the simulation
morotti_et_al_masterODEfile.m	integrates the following model components
morotti_et_al_barODEfile.m	beta-adrenergic (PKA) phosphorylation module
morotti_et_al_camkiiODEfile.m	CaMKII phosphorylation module
morotti_et_al_camODEfile.m	CaM module
morotti_et_al_eccODEfile.m	excitation-contraction coupling module

.mat files			initial conditions (obtained at 1 Hz pacing)
- yfin_WT_1Hz			WT model
- yfin_WT_1Hz_120sISO		WT model + 120 s ISO administration (0.1 uM)
- yfin_WT_NaGain_1Hz		WT model, with Na loading parameters ON
- yfin_OE_1Hz			CaMKII-OE model
- yfin_OE_loop_1Hz		CaMKII-OE model, with close CaMKII-Na-Ca_CaMKII loop
- yfin_OE_NoNaGain_1Hz		CaMKII-OE model, with Na loading parameters OFF
______________________________________________________________________________________


Reference:

S. Morotti, A.G. Edwards, A.D. McCulloch, D.M. Bers, E. Grandi (2014).
A novel computational model of mouse myocyte electrophysiology to assess the synergy
between Na+ loading and CaMKII. J Physiol, doi: 10.1113/jphysiol.2013.266676.

Please cite the above paper when using this model.
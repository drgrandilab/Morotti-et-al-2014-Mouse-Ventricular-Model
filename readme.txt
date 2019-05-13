Matlab code for Morotti et al. mouse ventricular model.

This model describes excitation-contraction coupling in the mouse ventricular myocyte
with integrated descriptions of Ca, Ca/calmodulin-dependent protein kinase II (CaMKII)
and protein kinase A signaling pathways.
This model was built upon the Soltis and Saucerman rabbit ventricular model (Biophys J.
2010 Oct 6;99(7):2038-47) and incorporates the experimentally known differences in 
mouse vs. rabbit electrophysiology, Ca handling, and kinase signaling and target 
phosphorylation.
Model parameters were also tuned to recapitulate the electrophysiologic changes during
chronic (transgenic) CaMKII overexpression.

______________________________________________________________________________________
Contents:

readme.txt			        this file

morotti_et_al_masterCompute.m	  loads initial conditions and runs the simulation
morotti_et_al_masterODEfile.m	  integrates the following model components
morotti_et_al_barODEfile.m	    beta-adrenergic (PKA) phosphorylation module
morotti_et_al_camkiiODEfile.m	  CaMKII phosphorylation module
morotti_et_al_camODEfile.m	    CaM module
morotti_et_al_eccODEfile.m	    excitation-contraction coupling module

.mat files			        initial conditions (obtained at 1 Hz pacing)
- yfin_WT_1Hz			      WT model
- yfin_WT_1Hz_120sISO		WT model + 120 s ISO administration (0.1 uM)
- yfin_WT_NaGain_1Hz		WT model, with Na loading parameters ON
- yfin_OE_1Hz			      CaMKII-OE model
- yfin_OE_loop_1Hz		  CaMKII-OE model, with close CaMKII-Na-Ca_CaMKII loop
- yfin_OE_NoNaGain_1Hz	CaMKII-OE model, with Na loading parameters OFF
______________________________________________________________________________________


Reference:

S. Morotti, A.G. Edwards, A.D. McCulloch, D.M. Bers, E. Grandi.
A novel computational model of mouse myocyte electrophysiology to assess the synergy
between Na+ loading and CaMKII.
J Physiol. 2014 Mar 15;592(6):1181-97.
doi: https://doi.org/10.1113/jphysiol.2013.266676

Please cite the above paper when using this model.

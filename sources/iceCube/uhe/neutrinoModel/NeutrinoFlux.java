package iceCube.uhe.neutrinoModel;

import java.io.*;
import java.util.*;

import numRecipes.*;

/**
<pre>
UHE Neutrino fluxes based on the following moels
are calculated from the table.

</pre>
<UL>
  <DT> <cite>Extremely High Energy Neutrinos and their Detection</cite>
       <TABLE CELLPADDING=2>
       <TR>
         <TD ALIGN="LEFT"><it>S.Yoshida et al.,</it></TD>
         <TD ALIGN="LEFT">
           <a href="http://www.journals.uchicago.edu/ApJ/journal/issues/ApJ/v479n2/35099/sc0.html">
            The Astrophysical Journal <b>479</b> 547-559 (1997)</a></TD>
       </TR>
       </TABLE><BR>
  <DT> <cite>Extremely High Energy Neutrinos, Neutrino Hot Dark Matter, and 
        the Highest Energy Cosmic Rays</cite> 
       <TABLE CELLPADDING=2>
       <TR>
         <TD ALIGN="LEFT"><it> S.Yoshida, G. Sigl and S. Lee,</it></TD>
         <TD ALIGN="LEFT">
            <a href="http://publish.aps.org/abstract/PRL/v81/p5505">
            Phys.Rev.Lett. <b>81</b> 5505 (1998)</a></TD>
       </TR>
       </TABLE><BR>
  <DT> <cite>Probing grand unified theories with cosmic-ray, gamma-ray,
  and neutrino astrophysics </cite> 
       <TABLE CELLPADDING=2>
       <TR>
         <TD ALIGN="LEFT"><it> G.Sigl et al ,</it></TD>
         <TD ALIGN="LEFT">
            <a href="http://publish.aps.org/abstract/PRD/v59/p043504">
            Phys.Rev.Rev. <b>D 59</b> 043504 (1998)</a></TD>
       </TR>
       </TABLE><BR>
  <DT> <cite>Ultrahigh-energy Neutrino Fluxes and Their Constraints</cite>
       <TABLE CELLPADDING=2>
       <TR>
         <TD ALIGN="LEFT"><it> O.E..Kalahsev, G. Sigl et al ,</it></TD>
         <TD ALIGN="LEFT">
            <a href="http://publish.aps.org/abstract/PRD/v66/p063004">
            Phys.Rev.Rev. <b>D 66</b> 063004 (2002)</a></TD>
       </TR>
       </TABLE><BR>
  <DT> <cite>GZK neutrinos after the Fermi-LAT diffuse photon flux measurment</cite>
       <TABLE CELLPADDING=2>
       <TR>
         <TD ALIGN="LEFT"><it> M.Ahlers et al ,</it></TD>
         <TD ALIGN="LEFT">
            <a href="http://dx.doi.org/10.1016/j.astropartphys.2010.06.003">
            Astropart. Phys. <b>34</b> 106-115 (2010)</a></TD>
       </TR>
       </TABLE><BR>
  <DT> <cite>Cosmogenic neutrinos: parameter space and detectability from PeV to ZeV</cite>
       <TABLE CELLPADDING=2>
       <TR>
         <TD ALIGN="LEFT"><it> K. Kotera, D.Allard, A.V.Olinto</it></TD>
         <TD ALIGN="LEFT">
            <a href="http://dx.doi.org/10.1088/1475-7516/2010/10/013">
            JCAP <b>10</b> 013 (2010)</a></TD>
       </TR>
       </TABLE><BR>
  <DT> <cite>Diffuse neutrino intensity from the inner jets of active galactic nuclei: 
  Impacts of external photon fields and the blazar</cite>
       <TABLE CELLPADDING=2>
       <TR>
         <TD ALIGN="LEFT"><it>Kohta Murase, Yoshiyuki Inoue, and Charles D. Dermer </it></TD>
         <TD ALIGN="LEFT">
            <a href="http://dx.doi.org/10.1103/PhysRevD.90.023007">
            Phys Rev D <b>90</b> 023007 (2014)</a></TD>
       </TR>
       </TABLE><BR>
  <DT> <cite>Testing the newborn pulsar origin of ultrahigh energy cosmic rays with EeV neutrinos</cite>
       <TABLE CELLPADDING=2>
       <TR>
         <TD ALIGN="LEFT"><it>Ke Fang, Kumiko Kotera, Kohta Murase, and Angela V. Olinto</it></TD>
         <TD ALIGN="LEFT">
            <a href="http://dx.doi.org/10.1103/PhysRevD.90.103005">
            Phys Rev D <b>90</b> 103005 (2014)</a></TD>
       </TR>
       </TABLE><BR>
</UL>

<pre>

Note that NO NEUTRINO OSCILATION was considered in these models.
One can access the original flux with getEFlux and/or getDFDLogE method.

The flux at the earth surface after propagation in space with neutrino
oscillation effect are also calculated based on the recent result of
the oscillation parameters.
One can access the flux after having the oscillation effect with getEFluxwzOsci
and/or getDFDLogEwzOsci.

The tables to contain the flux data is stored in $JAVADIR/classes/icecube/uhe/neutrinoModel
which are read out by this class.

Usage:
NeutrinoFlux  model-parameter

model-parameter 
                1  the GZK neutrinos m=0 Zmax = 2
                2  the GZK neutrinos m=2 Zmax = 2
                3  the GZK neutrinos m=2 Zmax = 4
                4  the GZK neutrinos m=4 Zmax = 4
		5  the GZK neutrinos m=4 Zmax = 5 gamma=1.5
		6  the GZK neutrinos m=7 Zmax = 5 gamma=1.5
                7  the Z-burst
                8  the Top Down (SUSY)
                9  the Top Down (QCD)
		10 the GZK neutrinos (by Sigl) m=3 Zmax=2 gamma=1 Emax=100ZeV
		11 the GZK neutrinos (by Sigl) m=5 Zmax=3 gamma=2 Emax=10ZeV
		12 the GZK neutrinos (by Ahles, Francis et al) Emin 10EeV best fit with the Fermi constraint
                   m=4.6 Zmax=2 gamma=2.5 Emax = 1ZeV
		13 the GZK neutrinos (by Ahles, Francis et al) Emax 10EeV max fit with the Fermi constraint
                   m=4.4 Zmax=2 gamma=2.1 Emax = 1ZeV
		14 the GZK neutrinos (by Ahles, Francis et al) Emax 10EeV min fit with the Fermi constraint
                   m=2.0 Zmax=2 gamma=2.88 Emax = 1ZeV
                15 the GZK neutrinos (by Kotera) Uniform, Emax 316 EeV
                16 the GZK neutrinos (by Kotera) SFR1, Emax 316 EeV
                17 the GZK neutrinos (by Kotera) GRB1, Emax 316 EeV
                18 the GZK neutrinos (by Kotera) FR2, Emax 316 EeV
		19 the GZK neutrinos (by Ahles, Francis et al) Emin 3EeV max fit with the Fermi constraint
                   m=5.35 Zmax=2 gamma=2.28 Emax = 1ZeV
		20 the GZK neutrinos (by Ahles, Francis et al) Emin 3EeV min fit with the Fermi constraint
                   m=2.0 Zmax=2 gamma=2.63 Emax = 1ZeV
		21 the GZK neutrinos (by Ahles, Francis et al) Emin 3EeV best fit with the Fermi constraint
                   m=4.05 Zmax=2 gamma=2.47 Emax = 1ZeV

                22 the GZK neutrinos (by Kotera) no evolution, CMB only, Emax 316 EeV
		23 the GZK neutrinos (by Ahles, Francis et al) maximal with the Fermi constraint
                                  		m=7 Zmax=4 gamma=2.0 Emax=3ZeV
                24 The baseline E<sup>-2</sup> flux - J(E)E^2 = 10^-8 [GeV/cm^2 sec sr] per neutrino flavor.
		25 GZK by simProp (Alosio et al) SFR, KneiskeEBL
		26 GZK by simProp (Alosio et al) SFR, SteckerEBL
		26 GZK by simProp (Alosio et al) SFR, CMB only
		28 GZK by simProp (Alosio et al) FR2, KneiskeEBL
		29 GZK by simProp (Alosio et al) FR2, SteckerEBL
		30 GZK by simProp (Alosio et al) FR2, CMB only
		31 AGN neutrinos (Murase et al) gamma = 2.0  CR loading factor 10
		32 AGN neutrinos (Murase et al) gamma = 2.3  CR loading factor 10
		33 Newborn pulsar neutrinos (Ke, Kumiko, Angela) SFR
		34 Newborn pulsar neutrinos (Ke, Kumiko, Angela) uniform
		35 AGN Blazars  (Padovani et al) HBL only Y_nu_gamma=0.8
		36 AGN Blazars  (Padovani et al) all Y_nu_gamma=0.8
		37 the GZK neutrinos (by Ahles, Francis et al) Emin 1EeV best fit with the Fermi constraint
                   m=3.2 Zmax=2 gamma=2.52 Emax = 1ZeV

particleID
                1  nu-e
                2  nu-mu
                3  nu-tau (No nu-tau in the GZK model) 
</pre>

*/


public class NeutrinoFlux implements Function{

    private String dataFile = null;
    private int numberOfFlavor = 2;
    private double[ ][ ] logEArray = new double[3][170];
    private double[ ][ ] EFluxArray = new double[3][170];
    private static final double ln10 = Math.log(10.0);
    private int dataNumber = 1; // Initial value. Should be changed later.
    private boolean isPowerLaw = false;
    /** all these parameters are concerned with the power law model (model ID 24) */
    private double  efluxPowerLaw = 1.0e-8;  // [GeV/cm2 sec sr] for E^-2 
    private double  energyBase = 1.0e5;  // [GeV] the energy at which efluxPowerLaw.
    private double  powerLawIndex = 2.0;
    private boolean alreadyOscillated = false;
    private double logEnergyMax = 12.0; // Maximum log(Energy/GeV) for the power-law model

    private static final double sinSq_2theta12 = 0.8700;
    // The value of sin^2(2*theta12) taken from KamLAND result 
    // (S. Abe et al., PRL, 100, 221803 (2008))


    /** Constructor: Reads out from the table stored in iceCube/uhe/neutrinoModel
    logE [eV]  dF/dE E^2 [eV/cm^2 sec sr]. Energy unit is transfered
    to GeV. */

    public NeutrinoFlux(int model) throws IOException {

        if(model==1){ // GZK Neutrino m=0 Zmax=2
            dataFile = "iceCube/uhe/neutrinoModel/gzk0_2.data";
	    dataNumber = 68;
        }else if(model==2){ // GZK Neutrino m=2 Zmax=2
            dataFile = "iceCube/uhe/neutrinoModel/gzk2_2.data";
	    dataNumber = 68;
        }else if(model==3){ // GZK Neutrino m=2 Zmax=4
            dataFile = "iceCube/uhe/neutrinoModel/gzk2_4.data";
	    dataNumber = 68;
        }else if(model==4){ // GZK Neutrino m=4 Zmax=4
            dataFile = "iceCube/uhe/neutrinoModel/gzk4_4.data";
	    dataNumber = 68;
        }else if(model==5){ // GZK Neutrino m=4 Zmax=5 gamma=1.5
            dataFile = "iceCube/uhe/neutrinoModel/gzk4_5.data";
	    dataNumber = 68;
        }else if(model==6){ // GZK Neutrino m=7 Zmax=5 gamma=1.5
            dataFile = "iceCube/uhe/neutrinoModel/gzk7_5.data";
	    dataNumber = 68;
        }else if(model==7){ // Z-burst model
            dataFile = "iceCube/uhe/neutrinoModel/z-burst.data";
	    dataNumber = 170;
	    numberOfFlavor = 3;
        }else if(model==8){ // Top Down SUZY model
            dataFile = "iceCube/uhe/neutrinoModel/td_susy.data";
	    dataNumber = 170;
	    numberOfFlavor = 3;
        }else if(model==9){ // Top Down QCD model
            dataFile = "iceCube/uhe/neutrinoModel/td_qcd.data";
	    dataNumber = 170;
	    numberOfFlavor = 3;
        }else if(model==10){ // GZK Neutrino (by Sigl) m=3 Zmax=2 gamma=1 Emax=100ZeV
            dataFile = "iceCube/uhe/neutrinoModel/gzk_siglMax.data";
	    dataNumber = 75;
	    numberOfFlavor = 1;
        }else if(model==11){ // GZK Neutrino (by Sigl) m=5 Zmax=3 gamma=2 Emax=10ZeV
            dataFile = "iceCube/uhe/neutrinoModel/gzk_siglLowE.data";
	    dataNumber = 80;
	    numberOfFlavor = 1;
        }else if(model==12){ // GZK Neutrino (by Francis) m=4.6 Zmax=2 gamma=2.5 Emax=1ZeV
            dataFile = "iceCube/uhe/neutrinoModel/gzk_fermi19_best.data";
	    dataNumber = 73;
	    numberOfFlavor = 1;
        }else if(model==13){ // GZK Neutrino (by Francis) m=4.4 Zmax=2 gamma=2.1 Emax=1ZeV
            dataFile = "iceCube/uhe/neutrinoModel/gzk_fermi19_max.data";
	    dataNumber = 73;
	    numberOfFlavor = 1;
        }else if(model==14){ // GZK Neutrino (by Francis) m=2.0 Zmax=2 gamma=2.8 Emax=1ZeV
            dataFile = "iceCube/uhe/neutrinoModel/gzk_fermi19_min.data";
	    dataNumber = 73;
	    numberOfFlavor = 1;
        }else if(model==15){ // GZK Neutrino (by Kotera) Uniform Emax=316EeV
            dataFile = "iceCube/uhe/neutrinoModel/gzk_kotera_uniform.data";
	    dataNumber = 100;
        }else if(model==16){ // GZK Neutrino (by Kotera) SFR1 Emax=316EeV
            dataFile = "iceCube/uhe/neutrinoModel/gzk_kotera_SFR1.data";
	    dataNumber = 100;
        }else if(model==17){ // GZK Neutrino (by Kotera) GRB1 Emax=316EeV
            dataFile = "iceCube/uhe/neutrinoModel/gzk_kotera_GRB1.data";
	    dataNumber = 100;
        }else if(model==18){ // GZK Neutrino (by Kotera) FR II Emax=316EeV
            dataFile = "iceCube/uhe/neutrinoModel/gzk_kotera_FR2.data";
	    dataNumber = 100;
        }else if(model==19){ // GZK Neutrino (by Francis) m=5.35 Zmax=2 gamma=2.28 Emax=1ZeV
            dataFile = "iceCube/uhe/neutrinoModel/gzk_fermi18_5_max.data";
	    dataNumber = 73;
	    numberOfFlavor = 1;
        }else if(model==19){ // GZK Neutrino (by Francis) m=5.35 Zmax=2 gamma=2.28 Emax=1ZeV
            dataFile = "iceCube/uhe/neutrinoModel/gzk_fermi18_5_max.data";
	    dataNumber = 73;
	    numberOfFlavor = 1;
        }else if(model==20){ // GZK Neutrino (by Francis) m=2.0 Zmax=2 gamma=2.63 Emax=1ZeV
            dataFile = "iceCube/uhe/neutrinoModel/gzk_fermi18_5_min.data";
	    dataNumber = 73;
	    numberOfFlavor = 1;
        }else if(model==21){ // GZK Neutrino (by Francis) m=4.05 Zmax=2 gamma=2.52 Emax=1ZeV
            dataFile = "iceCube/uhe/neutrinoModel/gzk_fermi18_5_best.data";
	    dataNumber = 73;
	    numberOfFlavor = 1;
        }else if(model==22){ // GZK Neutrino (by Kotera) uniform Emax=316EeV CMB only
            dataFile = "iceCube/uhe/neutrinoModel/gzk_kotera_uniform.data";
	    dataNumber = 100;
	    numberOfFlavor = 1;
        }else if(model==23){ // GZK Neutrino (by Ahlers) m=7 Zmax=4 gamma=2.0 Emax=3ZeV
            dataFile = "iceCube/uhe/neutrinoModel/gzk_fermi20_superMax.data";
	    dataNumber = 73;
	    numberOfFlavor = 1;
        }else if(model==24){ // E^-2 10^-8 [GeV/cm^2 sec sr]
	    numberOfFlavor = 1;
	    isPowerLaw = true;
        }else if(model==25){ // GZK Neutrino (simProp) SFR Kneiske EBL
            dataFile = "iceCube/uhe/neutrinoModel/gzk_simProp_SFR_KneiskeEBL.data";
	    dataNumber = 45;
	    numberOfFlavor = 1;
        }else if(model==26){ // GZK Neutrino (simProp) SFR Stecker EBL
            dataFile = "iceCube/uhe/neutrinoModel/gzk_simProp_SFR_SteckerEBL.data";
	    dataNumber = 45;
	    numberOfFlavor = 1;
        }else if(model==27){ // GZK Neutrino (simProp) CMB only
            dataFile = "iceCube/uhe/neutrinoModel/gzk_simProp_SFR_CMBonly.data";
	    dataNumber = 45;
	    numberOfFlavor = 1;
        }else if(model==28){ // GZK Neutrino (simProp) FR2 Kneiske EBL
            dataFile = "iceCube/uhe/neutrinoModel/gzk_simProp_AGN_KneiskeEBL.data";
	    dataNumber = 45;
	    numberOfFlavor = 1;
        }else if(model==29){ // GZK Neutrino (simProp) FR2 Stecker EBL
            dataFile = "iceCube/uhe/neutrinoModel/gzk_simProp_AGN_SteckerEBL.data";
	    dataNumber = 45;
	    numberOfFlavor = 1;
        }else if(model==30){ // GZK Neutrino (simProp) CMB only
            dataFile = "iceCube/uhe/neutrinoModel/gzk_simProp_AGN_CMBonly.data";
	    dataNumber = 45;
	    numberOfFlavor = 1;
        }else if(model==31){ // AGN neutrinos (Murase) gamma=2.0 CR loading factor 10
            dataFile = "iceCube/uhe/neutrinoModel/agn_murase_2_0.data";
	    dataNumber = 120;
	    numberOfFlavor = 3;
	    alreadyOscillated = true;
        }else if(model==32){ // AGN neutrinos (Murase) gamma=2.3 CR loading factor 10
            dataFile = "iceCube/uhe/neutrinoModel/agn_murase_2_3.data";
	    dataNumber = 120;
	    numberOfFlavor = 3;
	    alreadyOscillated = true;
        }else if(model==33){ // Newborn pulsars (by Ke et al ) SFR
            dataFile = "iceCube/uhe/neutrinoModel/pulsar_SFR.data";
	    dataNumber = 79;
	    numberOfFlavor = 1;
        }else if(model==34){ // Newborn pulsars (by Ke et al ) UNIFORM
            dataFile = "iceCube/uhe/neutrinoModel/pulsar_uniform.data";
	    dataNumber = 79;
	    numberOfFlavor = 1;
        }else if(model==35){ // Blazars (by Padovani et al ) HBL only
            dataFile = "iceCube/uhe/neutrinoModel/agn_padovani_hbl.data";
	    dataNumber = 10;
	    numberOfFlavor = 1;
        }else if(model==36){ // Blazars (by Padovani et al ) HBL only
            dataFile = "iceCube/uhe/neutrinoModel/agn_padovani_all.data";
	    dataNumber = 10;
	    numberOfFlavor = 1;
        }else if(model==37){ // GZK Neutrino (by Francis) m=3.2 Zmax=2 gamma=2.52 Emax=1ZeV
            dataFile = "iceCube/uhe/neutrinoModel/gzk_fermi18_0_best.data";
	    dataNumber = 73;
	    numberOfFlavor = 1;
        }else{
	    System.err.println("Illiegal parameters. mode(" + model + ")");
	    System.exit(0);
	}


	if(isPowerLaw) return;

	BufferedReader in =  
	    new BufferedReader(new InputStreamReader(ClassLoader.getSystemResourceAsStream(dataFile)));
        char separator = ' ';
	int n = 0; int k = 0;
	int sepstart = 0; int sep = 0;
        String buffer;
   
	for(n=0;n<numberOfFlavor;n++){

	    for(k=0;k<dataNumber;k++){
		buffer=in.readLine( );

		sepstart = 0;
		sep = buffer.indexOf(separator,sepstart+1);
		logEArray[n][k] =
		    Double.valueOf(buffer.substring(sepstart,sep)).doubleValue( )-9.0;
		// [GeV]

		sepstart = sep;
		sep = buffer.indexOf(separator,sepstart+1);
		EFluxArray[n][k] =
		    Double.valueOf(buffer.substring(sepstart,sep)).doubleValue( );

	    }
	}

	in.close( );

    }

    /** 
	Set the logEmax[GeV] for the power law model (model ID 24)
     */
    public void setLogEmax(double logEnergy){
	if(isPowerLaw) logEnergyMax = logEnergy;
    }
    /** 
	Return the logEmax[GeV] setted for the power law model (model ID 24)
     */
    public double getLogEmax(){
	if(isPowerLaw) return(logEnergyMax);
	else return(-1.0);
    }

    /** 
	Set the normalization E^2dF/dE [GeV/cm2 sec sr] for the power law model (model ID 24)
     */
    public void setEflux(double eFlux){
	if(isPowerLaw) efluxPowerLaw = eFlux;
    }

    /** 
	Get the normalization E^2dF/dE [GeV/cm2 sec sr] for the power law model (model ID 24)
     */
    public double getEflux(){
	if(isPowerLaw) return efluxPowerLaw;
	else return(-1.0);
    }

    /** 
	Set the powerlaw index for the power law model (model ID 24)
     */
    public void setPowerLawIndex(double gamma){
	if(isPowerLaw) powerLawIndex = gamma;
    }

    /** 
	Get the powerlaw index for the power law model (model ID 24)
     */
    public double getPowerLawIndex(){
	if(isPowerLaw) return powerLawIndex;
	else return(Double.POSITIVE_INFINITY);
    }

    /**
       tell if the setted model is a powerlaw model or not
     */
    public boolean isPowerLaw(){
	return isPowerLaw;
    }


    /** 
	Set the base energy [GeV] to define normalization E^2dF/dE [GeV/cm2 sec sr] 
	for the power law model (model ID 24)
     */
    public void setEnergyBase(double eBase){
	if(isPowerLaw) energyBase = eBase;
    }

    /** 
	Set the base energy [GeV] to define normalization E^2dF/dE [GeV/cm2 sec sr] 
	for the power law model (model ID 24)
     */
    public double getEnergyBase(){
	if(isPowerLaw) return energyBase;
	else return(Double.POSITIVE_INFINITY);
    }


    /** 
	<pre>
	calculate the differential Energy Flux [GeV /cm^2 sec sr] 

	logEnergy [GeV]

	particleID
                1  nu-e
                2  nu-mu
                3  nu-tau (No nu-tau in the GZK model) 

    */
    public double getEFlux(double logEnergy, int particleID){

	if(isPowerLaw){
	    double neutrinoEnergy = Math.pow(10.0,logEnergy);
	    double cutOffEnergy = Math.pow(10.0,logEnergyMax);
	    double xE = neutrinoEnergy/energyBase;
	    double powerLawTerm = Math.pow(xE,-powerLawIndex);
	    double xEcut = cutOffEnergy/neutrinoEnergy;
	    double bzCutOffTerm = Math.exp(-xEcut)*(xEcut-1.0)+1.0;
	    return (efluxPowerLaw*powerLawTerm*bzCutOffTerm*xE*xE);
	}

	if(numberOfFlavor!=1 && numberOfFlavor<3 && particleID == 3) return 0.0;

	double EFlux;
	if(numberOfFlavor==1){ // Already assuming oscillation
	    EFlux = 
	    Interpolation.mThPolynominalInterpolate(logEArray[0],
	    EFluxArray[0],dataNumber,logEnergy,4);
	}else{
	    EFlux = 
	    Interpolation.mThPolynominalInterpolate(logEArray[particleID-1],
	    EFluxArray[particleID-1],dataNumber,logEnergy,4);
	}
	
	double flux = EFlux*1.0e-9;

	return(flux);
    }


    /** 
	<pre>
	calculate the log differential Flux dF/dLogE [/cm^2 sec sr] 

	logEnergy [GeV]

	particleID
                1  nu-e
                2  nu-mu
                3  nu-tau (No nu-tau in the GZK model) 

    */
    public double getDFDLogE(double logEnergy, int particleID){

	if(numberOfFlavor!=1 && numberOfFlavor<3 && particleID == 3) return 0.0;

	if(isPowerLaw){
	    double neutrinoEnergy = Math.pow(10.0,logEnergy); // [GeV]
	    double cutOffEnergy = Math.pow(10.0,logEnergyMax);
	    double xE = neutrinoEnergy/energyBase;
	    double powerLawTerm = Math.pow(xE,-(powerLawIndex-1.0));
	    double xEcut = cutOffEnergy/neutrinoEnergy;
	    double bzCutOffTerm = Math.exp(-xEcut)*(xEcut-1.0)+1.0;
	    return (efluxPowerLaw*powerLawTerm*bzCutOffTerm/energyBase*ln10);
	}

	double EFlux;
	if(numberOfFlavor==1){ // Already assuming oscillation
	    EFlux = 
	    Interpolation.mThPolynominalInterpolate(logEArray[0],
	    EFluxArray[0],dataNumber,logEnergy,4);
	}else{
	    EFlux = 
	    Interpolation.mThPolynominalInterpolate(logEArray[particleID-1],
	    EFluxArray[particleID-1],dataNumber,logEnergy,4);
	}


	double energy = Math.pow(10.0,logEnergy+9.0); // [eV]
	double flux = EFlux/energy*ln10;

	return(flux);
    }


    /** 
	<pre>
	Calculate the neutrino flux at the earth surface after propgation in space
	with neutrino oscillation.

	Calculation is based on a paper of J. Jones et al., PRD, 69, 033004 (2004)
	Note that the mixing between nu_e and nu_tau is known to be very small, so
	the mixing is ignored.

	input: logEnergy [GeV],
	       nu_e flux before neutrino oscillation,
	       nu_mu flux before neutrino oscillation,
	       nu_tau flux before neutrino oscillation

	output: neutrino fluxes after neutrino oscillation 
	        (0: nu_e, 1:nu_mu, 2:nu_tau)
    */
    protected double[] getFluxwzOsci(double logEnergy, double nueflux,
				     double numuflux, double nutauflux){

	double[] flux = {nueflux, numuflux, nutauflux};
	if(!alreadyOscillated){
	    double nueflux_osci = nueflux - sinSq_2theta12*(2.*nueflux - numuflux - nutauflux)/4.;
	    double numuflux_osci = (numuflux + nutauflux)/2. 
		+ sinSq_2theta12*(2.*nueflux - numuflux - nutauflux)/8.;
	    double nutauflux_osci = numuflux_osci;

	    flux[0] = nueflux_osci;
	    flux[1] = numuflux_osci;
	    flux[2] = nutauflux_osci;
	}

	
	return flux;
    }


    /** 
	<pre>
	Calculate the differential Energy Flux [GeV /cm^2 sec sr] 
	after the propagation in the universe (with taking into account 
	the neutrino oscillation)

	input: logEnergy [GeV]
	
	output: neutrino fluxes after neutrino oscillation 
	        (0: nu_e, 1:nu_mu, 2:nu_tau)
    */
    public double[] getEFluxwzOsci(double logEnergy){

	double nueEflux = getEFlux(logEnergy, 1);
	double numuEflux = getEFlux(logEnergy, 2);
	double nutauEflux = getEFlux(logEnergy, 3);

	return getFluxwzOsci(logEnergy, nueEflux, numuEflux, nutauEflux);
    }


    /** 
	<pre>
	Calculate the differential Energy Flux [GeV /cm^2 sec sr] 
	after the propagation in the universe (with taking into account 
	the neutrino oscillation)

	input: logEnergy [GeV]
	
	particleID
                1  nu-e
                2  nu-mu
                3  nu-tau (No nu-tau in the GZK model) 
    */
    public double getEFluxwzOsci(double logEnergy, int particleID){

	double[] flux = getEFluxwzOsci(logEnergy);
	
	return flux[particleID-1];
    }


    /** 
	<pre>
	Calculate the log differential Flux dF/dLogE [/cm^2 sec sr] 
	after the propagation in the universe (with taking into account 
	the neutrino oscillation)

	input: logEnergy [GeV]
	
	output: neutrino fluxes after neutrino oscillation 
	        (0: nu_e, 1:nu_mu, 2:nu_tau)
    */
    public double[] getDFDLogEwzOsci(double logEnergy){

	double nueflux = getDFDLogE(logEnergy, 1);
	double numuflux = getDFDLogE(logEnergy, 2);
	double nutauflux = getDFDLogE(logEnergy, 3);
	
	return getFluxwzOsci(logEnergy, nueflux, numuflux, nutauflux);
    }


    /** 
	<pre>
	Calculate the log differential Flux dF/dLogE [/cm^2 sec sr] 
	after the propagation in the universe (with taking into account 
	the neutrino oscillation)

	Calculation is based on a paper of J. Jones et al., PRD, 69, 033004 (2004)
	Note that the mixing between nu_e and nu_tau is known to be very small, so
	the mixing is ignored.

	input: logEnergy [GeV]
	
	particleID
                1  nu-e
                2  nu-mu
                3  nu-tau (No nu-tau in the GZK model) 
    */
    public double getDFDLogEwzOsci(double logEnergy, int particleID){

	double[] flux = getDFDLogEwzOsci(logEnergy);
	
	return flux[particleID-1];
    }

    /**
       Interface to Integration class. This is for calculating
       the integfral flux.

       function index 0 :  dF/dLogE before the neutrino oscillation
       function index 1 :  dF/dLogE taking into account the neutrino oscillation
     */
    public double getFunction(int functionIndex, double[] parameters, 
			      double x){

	int particleID = (int )parameters[0];
	if(functionIndex == 0) return getDFDLogE(x,particleID);
	else return  getDFDLogEwzOsci(x, particleID);
    }



}


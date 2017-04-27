package  iceCube.uhe.muonModel;

import iceCube.uhe.particles.*;
import iceCube.uhe.interactions.*;
import iceCube.uhe.points.*;
import iceCube.uhe.propagation.*;
import numRecipes.*;

import java.io.*;

/**
   This class calculates differential flux dF/dLogE [/cm^2 sec sr]
   of neutrinos and charged leptons after propagation in the earth
   for a given model of primary cosmic neutrino production in the Universe.
   The primary flux of UHE cosmic neutrinos is given by the AtmMuonFlux
   class. The transfer(propagation) matrix of particles during the propagation
   in the earth, for example, FnuMuToTau = dN(nuMu -> tau)/dLogE_tau (E_nuMu, E_tau)
   is read out from the file generated a priori by PropagationMatrix.java 
   in the propagation package.

   The argument for the constructor, "model" is for the AtmMuonFlux class.
   Consult the details to the API document of the AtmMuonFlux.java
   in this package.
*/

public class PropagatingAtmMuonFlux implements Function{

    static int dimension = Particle.getDimensionOfLogEnergyMatrix();
    int inputParticle, outputParticle;
    int inLogE = -1;

    /** Flag to choose matrix with or without Glashow Resonance. */
    boolean includeGlashowResonance = true;

    /** Glashow Resonance **/
    double logYmin,logYmax;
    double[][] FnuEToNuE,FnuMuToNuE,FnuTauToNuE,FmuToNuE,FtauToNuE;
    double[][] FnuEToNuMu,FnuMuToNuMu,FnuTauToNuMu,FmuToNuMu,FtauToNuMu;
    double[][] FnuEToNuTau, FnuMuToNuTau,FnuTauToNuTau,FmuToNuTau,FtauToNuTau;
    double[][] FnuEToE,FnuMuToE,FnuTauToE,FmuToE,FtauToE;
    double[][] FnuEToMu,FnuMuToMu,FnuTauToMu,FmuToMu,FtauToMu;
    double[][] FnuEToTau,FnuMuToTau,FnuTauToTau,FmuToTau,FtauToTau;
    double[][] FnuEToHadron,FnuMuToHadron,FnuTauToHadron,FmuToHadron,FtauToHadron;
    AtmMuonBundleFlux atmMuFlux = null;
    double muonThresholdEnergy;
    //AtmMuonFluxCorsika atmMuFlux;
    //AtmMuonFlux atmMuFlux;
    private static final double ln10 = Math.log(10.0);

    /** Constructor */
    public PropagatingAtmMuonFlux( ) {
	generatePropagationMatrixArray();

	// Generate NeurinoFlux class
        atmMuFlux = new AtmMuonBundleFlux( );
	muonThresholdEnergy = atmMuFlux.getMuonThresholdEnergy();
        //atmMuFlux = new AtmMuonFluxCorsika( );
        //atmMuFlux = new AtmMuonFlux( );
    }

    /** Constructor */
    public PropagatingAtmMuonFlux(double muonThresholdEnergy, boolean cutoffExists) {
	generatePropagationMatrixArray();

	// Generate NeurinoFlux class
        atmMuFlux = new AtmMuonBundleFlux( );
	atmMuFlux.setMuonThresholdEnergy(muonThresholdEnergy);
	atmMuFlux.setCutOffFeature(cutoffExists);
	this.muonThresholdEnergy = muonThresholdEnergy;
    }


    /** Constructor */
    public PropagatingAtmMuonFlux(double alpha, double muonThresholdEnergy, 
				  boolean cutoffExists) {
	generatePropagationMatrixArray();

	// Generate NeurinoFlux class
        atmMuFlux = new AtmMuonBundleFlux( );
	atmMuFlux.setAlpha(alpha);
	atmMuFlux.setMuonThresholdEnergy(muonThresholdEnergy);
	atmMuFlux.setCutOffFeature(cutoffExists);
	if(cutoffExists){
	    System.err.println("alpha = " + alpha + " Eth = " + muonThresholdEnergy +
			       " with GZK cutoff");
	}else{
	    System.err.println("alpha = " + alpha + " Eth = " + muonThresholdEnergy +
			       " with NO GZK cutoff");
	}
	this.muonThresholdEnergy = muonThresholdEnergy;
    }

    /** Allocate memory for the propagation matrix array */
    private void generatePropagationMatrixArray(){
	/** For Glashow Resonance -begin **/
        FnuEToNuE= new double[dimension][dimension];
        FnuEToNuMu= new double[dimension][dimension];
        FnuEToNuTau= new double[dimension][dimension];
        FnuEToE= new double[dimension][dimension];
        FnuEToMu= new double[dimension][dimension];
        FnuEToTau= new double[dimension][dimension];
        FnuEToHadron= new double[dimension][dimension];
	/** For Glashow Resonance -end **/
        FnuMuToNuE= new double[dimension][dimension];
        FnuMuToNuMu= new double[dimension][dimension];
        FnuMuToNuTau= new double[dimension][dimension];
        FnuMuToE= new double[dimension][dimension];
        FnuMuToMu= new double[dimension][dimension];
        FnuMuToTau= new double[dimension][dimension];
        FnuMuToHadron= new double[dimension][dimension];
        FnuTauToNuE= new double[dimension][dimension];
        FnuTauToNuMu= new double[dimension][dimension];
        FnuTauToNuTau= new double[dimension][dimension];
        FnuTauToE= new double[dimension][dimension];
        FnuTauToMu= new double[dimension][dimension];
        FnuTauToTau= new double[dimension][dimension];
        FnuTauToHadron= new double[dimension][dimension];
        FmuToNuE= new double[dimension][dimension];
        FmuToNuMu= new double[dimension][dimension];
        FmuToNuTau= new double[dimension][dimension];
        FmuToE= new double[dimension][dimension];
        FmuToMu= new double[dimension][dimension];
        FmuToTau= new double[dimension][dimension];
        FmuToHadron= new double[dimension][dimension];
        FtauToNuE= new double[dimension][dimension];
        FtauToNuMu= new double[dimension][dimension];
        FtauToNuTau= new double[dimension][dimension];
        FtauToE= new double[dimension][dimension];
        FtauToMu= new double[dimension][dimension];
        FtauToTau= new double[dimension][dimension];
        FtauToHadron= new double[dimension][dimension];
    }


    /** Read the calculated propagation matrix */
    public void readMatrix(DataInputStream in) throws IOException {
        int iLogE;
        for(iLogE=0;iLogE<dimension;iLogE++){
            int jLogE;
            for(jLogE=0;jLogE<=iLogE;jLogE++){

                FnuEToNuE[iLogE][jLogE]= in.readDouble( );
		if(includeGlashowResonance){
		    FnuEToNuMu[iLogE][jLogE]= in.readDouble( );
		    FnuEToNuTau[iLogE][jLogE]= in.readDouble( );
		}else{
		    FnuEToNuMu[iLogE][jLogE] = 0.0;
		    FnuEToNuTau[iLogE][jLogE]= 0.0;
		}
                FnuEToE[iLogE][jLogE]= in.readDouble( );
                if(includeGlashowResonance){
		    FnuEToMu[iLogE][jLogE]= in.readDouble( );
		    FnuEToTau[iLogE][jLogE]= in.readDouble( );
		}else{
		    FnuEToMu[iLogE][jLogE]= 0.0;
		    FnuEToTau[iLogE][jLogE]= 0.0;
		}
                FnuEToHadron[iLogE][jLogE]= in.readDouble( );
                FnuMuToNuE[iLogE][jLogE]= in.readDouble( );
                FnuMuToNuMu[iLogE][jLogE]= in.readDouble( );
                FnuMuToNuTau[iLogE][jLogE]= in.readDouble( );
                FnuMuToE[iLogE][jLogE]= in.readDouble( );
                FnuMuToMu[iLogE][jLogE]= in.readDouble( );
                FnuMuToTau[iLogE][jLogE]= in.readDouble( );
                FnuMuToHadron[iLogE][jLogE]= in.readDouble( );
                FnuTauToNuE[iLogE][jLogE]= in.readDouble( );
                FnuTauToNuMu[iLogE][jLogE]= in.readDouble( );
                FnuTauToNuTau[iLogE][jLogE]= in.readDouble( );
                FnuTauToE[iLogE][jLogE]= in.readDouble( );
                FnuTauToMu[iLogE][jLogE]= in.readDouble( );
                FnuTauToTau[iLogE][jLogE]= in.readDouble( );
                FnuTauToHadron[iLogE][jLogE]= in.readDouble( );
                FmuToNuE[iLogE][jLogE]= in.readDouble( );
                FmuToNuMu[iLogE][jLogE]= in.readDouble( );
                FmuToNuTau[iLogE][jLogE]= in.readDouble( );
                FmuToE[iLogE][jLogE]= in.readDouble( );
                FmuToMu[iLogE][jLogE]= in.readDouble( );
                FmuToTau[iLogE][jLogE]= in.readDouble( );
                FmuToHadron[iLogE][jLogE]= in.readDouble( );
                FtauToNuE[iLogE][jLogE]= in.readDouble( );
                FtauToNuMu[iLogE][jLogE]= in.readDouble( );
                FtauToNuTau[iLogE][jLogE]= in.readDouble( );
                FtauToE[iLogE][jLogE]= in.readDouble( );
                FtauToMu[iLogE][jLogE]= in.readDouble( );
                FtauToTau[iLogE][jLogE]= in.readDouble( );
                FtauToHadron[iLogE][jLogE]= in.readDouble( );

            }
        }

        in.close( );

	// Set the threshold energy of muons in a bundle 
	// at the Earth surface
	if(atmMuFlux!=null){ // Atmospheric muon object exists
	    double thresholdE = (muonThresholdEnergy+atmMuFlux.criticalEnergy)/
		getAverageMuonEnergyLossAfterPropagation()-atmMuFlux.criticalEnergy;
	    atmMuFlux.setMuonThresholdEnergy(thresholdE);
	    System.err.println("JULIET: Muon E threshold = " + 
	    	       atmMuFlux.getMuonThresholdEnergy());
	}else{
	    System.err.println("No object of the atm muon flux can be found");
	    System.exit(0);
	}
	    

    }

    /** set the mode on  the flux caulation 
	<pre>
	boolean includeFluctuationEffects : 
                    true  default. calculate the flux taking into account flucuation of 
                          muon energies due to the EAS cascading. 
                          CascadeFluctuationFactory class does this part of calculation.
                    false calculate the bare flux given by the Elbert formula.
                          an energy of muon bundles is associated with the energy of primary
                          cosmic ray by one-on-one relation.

	boolean fluctuateEventByEvent : Valid when includeFluctuationEffects = true

	   true :  Flux is calculated such that an event by event flucuation is included.
                   Useful to evaluate flucuation of number of background events
           false: default. the flucuation effects is included in average way. i.e.,
                  The systematic factor of the flux shift due the cascade flucuation
                  is calculated and multiplied in the flux calculation.
                  Useful to the MC-data fitting, calculation of average number of events
	</pre>
     **/
    public void setFluxCalculationMode(boolean includeFluctuationEffects,
				       boolean fluctuateEventByEvent){
	atmMuFlux.setFluxCalculationMode(includeFluctuationEffects,fluctuateEventByEvent);
    }


    /** Set Confidence level of the Energy rato parameter 
	R = (E_muon/E0)/ Bar(E_muon/E0). If this value is non zero (zero in its default)
	and (includeFluctuationEffect fluctuateEventByEvent) = (true, true)
	then all the relevant fluxes are calculated based on a fixed value of R 
        given by this confidence level. If zero (default), the fluxes are given
	by the interal over R which correponds to taking the central value of R-distribution.
    */
    public void setConfidenceLevelOfFluctuation(double prob){
	atmMuFlux.setConfidenceLevelOfFluctuation(prob);
    }

    /** set the GZK cutoff feature or not */
    public void setCutOffFeature(boolean cutoffExists){
	atmMuFlux.setCutOffFeature(cutoffExists);
    }


    /** Calculate dF/dLogE [/cm^2 sec sr] for nu-mu */
    public double getDFMuDLogE(double logE, double cos_th){

	int jLogE = (int)((logE - Particle.getLogEnergyMinimum())
			  /Particle.getDeltaLogEnergy());
	double count = 0.0;
	int iLogE;

	for(iLogE=jLogE;iLogE<dimension;iLogE++){

	    double logEprimary = Particle.getLogEnergyMinimum( ) + 
		Particle.getDeltaLogEnergy( )*(double )iLogE;

	    if(FmuToMu[iLogE][jLogE]>0.0){
		count += atmMuFlux.getDFDLogE(logEprimary,cos_th)*FmuToMu[iLogE][jLogE]; // Mu to Mu
	    }
	}
	return(count);
    }


    /** Calculate dF^2/dLogEcrDLogE [/cm^2 sec sr] for nu-mu 
	Because the cosmic ray follows the rapidly falling spectrum, the flux is 
	averaged over a given log(cosmic ray energy [GeV]) in a half width given
	by halfWidthOfLogE. Because the bin width of the propagation matrix is
	Particle.getDeltaLogEnergy() (presently 0.01), the width should be set
	as equal to or larger than this value.
	<pre>
	double logCosmicRayEnergy     : log10(Cosmic Ray Energy [GeV]) producing muon bundles with logE
	double halfWidthOfLogE        : a half width of integration for averaging the cosmic ray flux
	double logE                   : log10(Muon Bundle Energy [GeV])
	double cos_th                 : cos(zenith angle)
	</pre>
     */
    public double getDFMuDLogCREDLogE(double logCosmicRayEnergy, double halfWidthOfLogE,
				      double logE, double cos_th){

	// impossible - CR energy should be higher than mu energies.
	if(logCosmicRayEnergy+halfWidthOfLogE < logE ) return 0.0;

	int jLogE = (int)((logE - Particle.getLogEnergyMinimum())
			  /Particle.getDeltaLogEnergy());
	double count = 0.0;
	int iLogE;

	for(iLogE=jLogE;iLogE<dimension;iLogE++){

	    double logEprimary = Particle.getLogEnergyMinimum( ) + 
		Particle.getDeltaLogEnergy( )*(double )iLogE;

	    if(FmuToMu[iLogE][jLogE]>0.0){
		count += 
		    atmMuFlux.integralDFDLogCREDLogEOverCREnergy(logCosmicRayEnergy-halfWidthOfLogE,
								 logCosmicRayEnergy+halfWidthOfLogE,
								 logEprimary,cos_th)/
		    (2.0*halfWidthOfLogE)*FmuToMu[iLogE][jLogE]; // Mu to Mu
	    }
	}
	return(count);
    }

    /** Calculate and return the averaged relative energy loss of a muon 
      after its propagation in the Earth. This value is used as 
      the threshold energy of muons in a bundle
      in the calulcation of the parent Cosmic Ray energy with the
      <pre>
      AtmMuonBundle.getEffectiveEnergyOfCRs(double muonBundleEnergy, 
                                           double cosTheta_ice3)
      </pre>
    */
    
    protected double getAverageMuonEnergyLossAfterPropagation(){
	double logReferenceEnergy = 7.0; // 10^7 GeV

	int jLogE = (int)((logReferenceEnergy - Particle.getLogEnergyMinimum())
			  /Particle.getDeltaLogEnergy());
	double averagedLogEprimary = 0.0; double count = 0;
	int iLogE;
	for(iLogE=jLogE;iLogE<dimension;iLogE++){
	    double logEprimary = Particle.getLogEnergyMinimum( ) + 
		Particle.getDeltaLogEnergy( )*(double )iLogE;
	    if(FmuToMu[iLogE][jLogE]>0.0){
		averagedLogEprimary += logEprimary*FmuToMu[iLogE][jLogE];
		count += FmuToMu[iLogE][jLogE];
	    }
	}
	averagedLogEprimary /= count;

	double deltaLogE = averagedLogEprimary - logReferenceEnergy;
	double energyRatio = Math.pow(10.0,-deltaLogE);

	return energyRatio;
    }

    /** or Glashow Resonance */
    // Switch includeGlashowResonance flag.
    public void whetherPropagationMatrixWithGlashowResonance(boolean flag){
	includeGlashowResonance = flag;
    }

    /** integrate getDFDLogCREDLogE(logEcr,logEmu,cosTheta) over logEcr
	and compare with getDFDLogE(logEmu,cosTheta) for the consistency check
    */
    protected void listFluxes(double logEnergy, double cosTheta){

	double[] param = new double[2];
	param[0] = logEnergy;
	param[1]= cosTheta;
	// integral from 100TeV to 1000 EeV
	double dFdLogEbyIntegral = 
	    Integration.RombergIntegral(this,1,param, 4.0,12.0); 
	System.err.println(" muon bundle flux by integral(logE=" + logEnergy + ", cosZen=" + cosTheta +
			   ") = " + dFdLogEbyIntegral);

	double dFDLogE = getDFMuDLogE(logEnergy, cosTheta);
	System.err.println(" muon bundle flux(logE=" + logEnergy + ", cosZen=" + cosTheta +
			   ") = " + dFDLogE);
    }


    /** Implementation of interface numRecipes.Function
	used for numerical integration of the d^F/dLogEcrDlogEmu over
	logEcr, the cosmic ray energy. This should be equal to
	dF/dLogEmu, the return value of the method getDFMuDLogE(logE, cosZenith).
	This implementation is only for debugging purpose.

	<pre>
	double x = logCosmicRayEnergy

	parameter[0]   : log10(Muon Energy [GeV])
	parameter[1]   : cosine(Zenith angle)

	</pre>

    */
    public double getFunction(int functionIndex, double[] parameters, 
	double x){

	double logEnergy = parameters[0];
	double cosZenith = parameters[1];
	double halfWidthOfLogE = Particle.getDeltaLogEnergy()*0.5;
	return getDFMuDLogCREDLogE(x,halfWidthOfLogE,logEnergy,cosZenith);
    }
}

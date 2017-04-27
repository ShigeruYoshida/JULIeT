package  iceCube.uhe.neutrinoModel;

import iceCube.uhe.particles.*;
import iceCube.uhe.interactions.*;
import iceCube.uhe.points.*;
import iceCube.uhe.propagation.*;
import java.io.*;
import java.util.*;

/**
   This class calculates differential flux dF/dLogE [/cm^2 sec sr]
   of neutrinos and charged leptons after propagation in the earth
   for a given model of primary cosmic neutrino production in the Universe.
   The primary flux of UHE cosmic neutrinos is given by the NeutrinoFlux
   class with taking into account the neutrino oscillation effect. 
   The transfer(propagation) matrix of particles during the propagation
   in the earth, for example, FnuMuToTau = dN(nuMu -> tau)/dLogE_tau (E_nuMu, E_tau)
   is read out from the file generated a priori by PropagationMatrix.java 
   in the propagation package.

   The argument for the constructor, "model" is for the NeutrinoFlux class.
   Consult the detais to the API document of the NeutrinoFlux.java
   in this package.
*/

public class PropagatingNeutrinoFlux {

    static int dimension = Particle.getDimensionOfLogEnergyMatrix();
    int inputParticle, outputParticle;
    int inLogE = -1;
    double logYmin,logYmax;

    /** Flag to choose matrix with or without Glashow Resonance. */
    boolean includeGlashowResonance = true;

    /** For Glashow Resonance FnuETo(Mu/Tau flavor) are added. **/
    double[][] FnuEToNuE,FnuMuToNuE,FnuTauToNuE,FmuToNuE,FtauToNuE;
    double[][] FnuEToNuMu,FnuMuToNuMu,FnuTauToNuMu,FmuToNuMu,FtauToNuMu;
    double[][] FnuEToNuTau,FnuMuToNuTau,FnuTauToNuTau,FmuToNuTau,FtauToNuTau;
    double[][] FnuEToE,FnuMuToE,FnuTauToE,FmuToE,FtauToE;
    double[][] FnuEToMu,FnuMuToMu,FnuTauToMu,FmuToMu,FtauToMu;
    double[][] FnuEToTau,FnuMuToTau,FnuTauToTau,FmuToTau,FtauToTau;
    double[][] FnuEToHadron,FnuMuToHadron,FnuTauToHadron,FmuToHadron,FtauToHadron;
    private List neutFluxList = null;
    private List modelNumberList = null;

    private static final double ln10 = Math.log(10.0);

    /** Constructor */
    public PropagatingNeutrinoFlux(int model) throws IOException {

	generateMatrix();
	// Generate the list to accomodate neutrinoFlux object(s)
	neutFluxList = new LinkedList();
	modelNumberList = new LinkedList();

	// Generate the neutrinoFlux object with the given model number
	generateNeutrinoFluxObject(model);
    }

    /** This constroctor is for the subclass QuickPropagatingNeutrinoFlux */
    protected PropagatingNeutrinoFlux(int model, boolean nomatrix) throws IOException {

	// Generate the list to accomodate neutrinoFlux object(s)
	neutFluxList = new LinkedList();
	modelNumberList = new LinkedList();

	// Generate the neutrinoFlux object with the given model number
	generateNeutrinoFluxObject(model);
    }

    private void generateMatrix() {
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

    /** Generate Neutrino flux object and add it to the list */
    protected void generateNeutrinoFluxObject(int model) throws IOException {
        neutFluxList.add(new NeutrinoFlux(model));
	modelNumberList.add(new Integer(model));
    }


    /** Read the calculated propagation matrix */
    public void readMatrix(DataInputStream in) throws IOException {
        int iLogE;
        for(iLogE=0;iLogE<dimension;iLogE++){
            int jLogE;
            for(jLogE=0;jLogE<=iLogE;jLogE++){

		/** For Glashow Resonance -begin **/
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
		    FnuEToMu[iLogE][jLogE] = 0.0;
		    FnuEToTau[iLogE][jLogE]= 0.0;
		}
                FnuEToHadron[iLogE][jLogE]= in.readDouble( );
		/** For Glashow Resonance -end **/
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

    }


    /**
       Search for the neutrinoFlux objects stored in the list
       and return it if exists.
    */
    protected NeutrinoFlux getNeutrinoFluxFromList(int model){

	ListIterator neutFluxIterator = neutFluxList.listIterator();
	ListIterator modelNumberIterator = modelNumberList.listIterator();
	NeutrinoFlux neutFlux = null;
	boolean found = false;
	while(modelNumberIterator.hasNext()){
	    int neutModel = ((Integer )modelNumberIterator.next()).intValue();
	    //System.err.println(" given model=" + model + " listed model = " + neutModel);
	    neutFlux = (NeutrinoFlux)neutFluxIterator.next();
	    if(neutModel == model){ // The object has been already generated
		found = true;
		break;
	    }
	}
	if(found) return neutFlux;
	else return null;
    }



    /** Calculate dF/dLogE [/cm^2 sec sr] for nu-e 
	The neutrino flux is calculated for the model
        you specified in calling the constructor.
     */
    public double getDFNuEDLogE(double logE) throws IOException {

	int jLogE = (int)((logE - Particle.getLogEnergyMinimum())
			  /Particle.getDeltaLogEnergy());
	return getDFNuEDLogE(jLogE);
    }

    /** Calculate dF/dLogE [/cm^2 sec sr] for nu-e 
	The neutrino flux is calculated for the model
        you specified in calling the constructor.
     */
    public double getDFNuEDLogE(int jLogE) throws IOException {
	
	ListIterator modelNumberIterator = modelNumberList.listIterator();
	int model = ((Integer )modelNumberIterator.next()).intValue();
	return getDFNuEDLogE(model,jLogE);

    }

    /** Calculate dF/dLogE [/cm^2 sec sr] for nu-e 
	The neutrino flux is calculated for the model
        you specified in the argument.
     */
    public double getDFNuEDLogE(int model, double logE) throws IOException {

	int jLogE = (int)((logE - Particle.getLogEnergyMinimum())
			  /Particle.getDeltaLogEnergy());
	return getDFNuEDLogE(model,jLogE);
    }

    /** Calculate dF/dLogE [/cm^2 sec sr] for nu-e 
	The neutrino flux is calculated for the model
        you specified in the argument.
     */
    public double getDFNuEDLogE(int model, int jLogE) throws IOException {

	NeutrinoFlux neutFlux = getNeutrinoFluxFromList(model);
	if(neutFlux == null){ // This is a new model. Needs to generate the object
	    generateNeutrinoFluxObject(model);
	    neutFlux = getNeutrinoFluxFromList(model);
	    System.err.println("JULIeT: new NeutrinoFlux object generated with model "
			       + model);
	}


	double logE;
	double count = 0.0;
	int iLogE;

	for(iLogE=jLogE;iLogE<dimension;iLogE++){

	    logE = Particle.getLogEnergyMinimum( ) + 
		Particle.getDeltaLogEnergy( )*(double )iLogE;

	    // get nuetrino flux at the earth surface after propagation in space
	    // with neutrino oscillation
	    double[] nuflux_osci = neutFlux.getDFDLogEwzOsci(logE);
	    
	    if(FnuEToNuE[iLogE][jLogE]>0.0){
		count += nuflux_osci[0]*FnuEToNuE[iLogE][jLogE]; // nuE to nuE
	    }
	    if(FnuMuToNuE[iLogE][jLogE]>0.0){
		count += nuflux_osci[1]*FnuMuToNuE[iLogE][jLogE]; // nuMu to nuE
	    }
	    if(FnuTauToNuE[iLogE][jLogE]>0.0){
		count += nuflux_osci[2]*FnuTauToNuE[iLogE][jLogE]; // nuTau to nuE
	    }

	}
	return(count);
    }



    /** Calculate dF/dLogE [/cm^2 sec sr] for nu-mu
	The neutrino flux is calculated for the model
        you specified in calling the constructor.
    */
    public double getDFNuMuDLogE(double logE) throws IOException {

	int jLogE = (int)((logE - Particle.getLogEnergyMinimum())
			  /Particle.getDeltaLogEnergy());
	return getDFNuMuDLogE(jLogE);
    }

    /** Calculate dF/dLogE [/cm^2 sec sr] for nu-mu 
	The neutrino flux is calculated for the model
        you specified in calling the constructor.
     */
    public double getDFNuMuDLogE(int jLogE) throws IOException {
	
	ListIterator modelNumberIterator = modelNumberList.listIterator();
	int model = ((Integer )modelNumberIterator.next()).intValue();
	return getDFNuMuDLogE(model,jLogE);

    }

    /** Calculate dF/dLogE [/cm^2 sec sr] for nu-mu 
	The neutrino flux is calculated for the model
        you specified in the argument.
     */
    public double getDFNuMuDLogE(int model, double logE) throws IOException {

	int jLogE = (int)((logE - Particle.getLogEnergyMinimum())
			  /Particle.getDeltaLogEnergy());
	return getDFNuMuDLogE(model,jLogE);
    }

    /** Calculate dF/dLogE [/cm^2 sec sr] for nu-mu 
	The neutrino flux is calculated for the model
        you specified in the argument.
     */
    public double getDFNuMuDLogE(int model, int jLogE) throws IOException {

	NeutrinoFlux neutFlux = getNeutrinoFluxFromList(model);
	if(neutFlux == null){ // This is a new model. Needs to generate the object
	    generateNeutrinoFluxObject(model);
	    neutFlux = getNeutrinoFluxFromList(model);
	    System.err.println("JULIeT: new NeutrinoFlux object generated with model "
			       + model);
	}

	double logE;
	double count = 0.0;
	int iLogE;

	for(iLogE=jLogE;iLogE<dimension;iLogE++){

	    logE = Particle.getLogEnergyMinimum( ) + 
		Particle.getDeltaLogEnergy( )*(double )iLogE;

	    // get nuetrino flux at surface after propagation in space with
	    // neutrino oscillation
	    double[] nuflux_osci = neutFlux.getDFDLogEwzOsci(logE);

	    /** For Glashow Resonance -begin **/
	    if(FnuEToNuMu[iLogE][jLogE]>0.0){
		count += nuflux_osci[0]*FnuEToNuMu[iLogE][jLogE]; // nuE to nuMu
	    }
	    /** For Glashow Resonance -end **/
	    if(FnuMuToNuMu[iLogE][jLogE]>0.0){
		count += nuflux_osci[1]*FnuMuToNuMu[iLogE][jLogE]; // nuMu to nuMu
	    }
	    if(FnuTauToNuMu[iLogE][jLogE]>0.0){
		count += nuflux_osci[2]*FnuTauToNuMu[iLogE][jLogE]; // nuTau to nuMu
	    }
	}
	return(count);
    }
 


    /** Calculate dF/dLogE [/cm^2 sec sr] for nu-tau 
	The neutrino flux is calculated for the model
        you specified in calling the constructor.
     */
    public double getDFNuTauDLogE(double logE) throws IOException {

	int jLogE = (int)((logE - Particle.getLogEnergyMinimum())
			  /Particle.getDeltaLogEnergy());
	return getDFNuTauDLogE(jLogE);
    }

    /** Calculate dF/dLogE [/cm^2 sec sr] for nu-tau 
	The neutrino flux is calculated for the model
        you specified in calling the constructor.
     */
    public double getDFNuTauDLogE(int jLogE) throws IOException {
	
	ListIterator modelNumberIterator = modelNumberList.listIterator();
	int model = ((Integer )modelNumberIterator.next()).intValue();
	return getDFNuTauDLogE(model,jLogE);

    }

    /** Calculate dF/dLogE [/cm^2 sec sr] for nu-tau 
	The neutrino flux is calculated for the model
        you specified in the argument.
     */
    public double getDFNuTauDLogE(int model, double logE) throws IOException {

	int jLogE = (int)((logE - Particle.getLogEnergyMinimum())
			  /Particle.getDeltaLogEnergy());
	return getDFNuTauDLogE(model,jLogE);
    }


    /** Calculate dF/dLogE [/cm^2 sec sr] for nu-tau 
	The neutrino flux is calculated for the model
        you specified in the argument.
     */
    public double getDFNuTauDLogE(int model, int jLogE) throws IOException {

	NeutrinoFlux neutFlux = getNeutrinoFluxFromList(model);
	if(neutFlux == null){ // This is a new model. Needs to generate the object
	    generateNeutrinoFluxObject(model);
	    neutFlux = getNeutrinoFluxFromList(model);
	    System.err.println("JULIeT: new NeutrinoFlux object generated with model "
			       + model);
	}

	double logE;
	double count = 0.0;
	int iLogE;

	for(iLogE=jLogE;iLogE<dimension;iLogE++){

	    logE = Particle.getLogEnergyMinimum( ) + 
		Particle.getDeltaLogEnergy( )*(double )iLogE;

	    // get nuetrino flux at surface after propagation in space with
	    // neutrino oscillation
	    double[] nuflux_osci = neutFlux.getDFDLogEwzOsci(logE);

	    /** For Glashow Resonance -begin **/
	    if(FnuEToNuTau[iLogE][jLogE]>0.0){
		count += nuflux_osci[0]*FnuEToNuTau[iLogE][jLogE]; // nuE to nuTau
	    }
	    /** For Glashow Resonance -end **/
	    if(FnuMuToNuTau[iLogE][jLogE]>0.0){
		count += nuflux_osci[1]*FnuMuToNuTau[iLogE][jLogE]; // nuMu to nuTau
	    }
	    if(FnuTauToNuTau[iLogE][jLogE]>0.0){
		count += nuflux_osci[2]*FnuTauToNuTau[iLogE][jLogE]; // nuTau to nuTau
	    }
	}
	return(count);
    }
 


    /** Calculate dF/dLogE [/cm^2 sec sr] for mu 
	The neutrino flux is calculated for the model
        you specified in calling the constructor.
     */
    public double getDFMuDLogE(double logE) throws IOException {

	int jLogE = (int)((logE - Particle.getLogEnergyMinimum())
			  /Particle.getDeltaLogEnergy());
	return getDFMuDLogE(jLogE);
    }

    /** Calculate dF/dLogE [/cm^2 sec sr] for mu
	The neutrino flux is calculated for the model
        you specified in calling the constructor.
     */
    public double getDFMuDLogE(int jLogE) throws IOException {
	
	ListIterator modelNumberIterator = modelNumberList.listIterator();
	int model = ((Integer )modelNumberIterator.next()).intValue();
	return getDFMuDLogE(model,jLogE);

    }

    /** Calculate dF/dLogE [/cm^2 sec sr] for mu 
	The neutrino flux is calculated for the model
        you specified in the argument.
     */
    public double getDFMuDLogE(int model, double logE) throws IOException {

	int jLogE = (int)((logE - Particle.getLogEnergyMinimum())
			  /Particle.getDeltaLogEnergy());
	return getDFMuDLogE(model,jLogE);
    }

    /** Calculate dF/dLogE [/cm^2 sec sr] for mu 
	The neutrino flux is calculated for the model
        you specified in the argument.
     */
    public double getDFMuDLogE(int model, int jLogE) throws IOException {

	NeutrinoFlux neutFlux = getNeutrinoFluxFromList(model);
	if(neutFlux == null){ // This is a new model. Needs to generate the object
	    generateNeutrinoFluxObject(model);
	    neutFlux = getNeutrinoFluxFromList(model);
	    System.err.println("JULIeT: new NeutrinoFlux object generated with model "
			       + model);
	}

	double logE;
	double count = 0.0;
	int iLogE;

	for(iLogE=jLogE;iLogE<dimension;iLogE++){

	    logE = Particle.getLogEnergyMinimum( ) + 
		Particle.getDeltaLogEnergy( )*(double )iLogE;

	    // get nuetrino flux at surface after propagation in space with
	    // neutrino oscillation
	    double[] nuflux_osci = neutFlux.getDFDLogEwzOsci(logE);

	    /** For Glashow Resonance -begin **/
	    if(FnuEToMu[iLogE][jLogE]>0.0){
		count += nuflux_osci[0]*FnuEToMu[iLogE][jLogE]; // nuE to Mu
	    }
	    /** For Glashow Resonance -end **/
	    if(FnuMuToNuE[iLogE][jLogE]>0.0){
		count += nuflux_osci[1]*FnuMuToMu[iLogE][jLogE]; // nuMu to Mu
	    }
	    if(FnuTauToNuE[iLogE][jLogE]>0.0){
		count += nuflux_osci[2]*FnuTauToMu[iLogE][jLogE]; // nuTau to Mu
	    }
	}
	return(count);
    }
 


    /** Calculate dF/dLogE [/cm^2 sec sr] for tau 
	The neutrino flux is calculated for the model
        you specified in calling the constructor.
     */
    public double getDFTauDLogE(double logE) throws IOException {

	int jLogE = (int)((logE - Particle.getLogEnergyMinimum())
			  /Particle.getDeltaLogEnergy());
	return getDFTauDLogE(jLogE);
    }

    /** Calculate dF/dLogE [/cm^2 sec sr] for tau 
	The neutrino flux is calculated for the model
        you specified in calling the constructor.
     */
    public double getDFTauDLogE(int jLogE) throws IOException {
	
	ListIterator modelNumberIterator = modelNumberList.listIterator();
	int model = ((Integer )modelNumberIterator.next()).intValue();
	return getDFTauDLogE(model,jLogE);

    }

    /** Calculate dF/dLogE [/cm^2 sec sr] for tau 
	The neutrino flux is calculated for the model
        you specified in the argument.
     */
    public double getDFTauDLogE(int model, double logE) throws IOException {

	int jLogE = (int)((logE - Particle.getLogEnergyMinimum())
			  /Particle.getDeltaLogEnergy());
	return getDFTauDLogE(model,jLogE);
    }


    /** Calculate dF/dLogE [/cm^2 sec sr] for tau 
	The neutrino flux is calculated for the model
        you specified in the argument.
     */
    public double getDFTauDLogE(int model, int jLogE) throws IOException {

	NeutrinoFlux neutFlux = getNeutrinoFluxFromList(model);
	if(neutFlux == null){ // This is a new model. Needs to generate the object
	    generateNeutrinoFluxObject(model);
	    neutFlux = getNeutrinoFluxFromList(model);
	    System.err.println("JULIeT: new NeutrinoFlux object generated with model "
			       + model);
	}

	double logE;
	double count = 0.0;
	int iLogE;

	for(iLogE=jLogE;iLogE<dimension;iLogE++){

	    logE = Particle.getLogEnergyMinimum( ) + 
		Particle.getDeltaLogEnergy( )*(double )iLogE;

	    // get nuetrino flux at surface after propagation in space with
	    // neutrino oscillation
	    double[] nuflux_osci = neutFlux.getDFDLogEwzOsci(logE);

	    /** For Glashow Resonance -begin **/
	    if(FnuEToTau[iLogE][jLogE]>0.0){
		count += nuflux_osci[0]*FnuEToTau[iLogE][jLogE]; // nuE to Mu
	    }
	    /** For Glashow Resonance -end **/
	    if(FnuMuToNuE[iLogE][jLogE]>0.0){
		count += nuflux_osci[1]*FnuMuToTau[iLogE][jLogE]; // nuMu to Tau
	    }
	    if(FnuTauToNuE[iLogE][jLogE]>0.0){
		count += nuflux_osci[2]*FnuTauToTau[iLogE][jLogE]; // nuTau to Tau
	    }
	}
	return(count);
    }
 
    /** or Glashow Resonance */
    // Switch includeGlashowResonance flag.
    public void whetherPropagationMatrixWithGlashowResonance(boolean flag){
	includeGlashowResonance = flag;
    }

}

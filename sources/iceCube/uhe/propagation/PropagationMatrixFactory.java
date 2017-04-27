package  iceCube.uhe.propagation;

import iceCube.uhe.particles.*;
import iceCube.uhe.propagation.*;
import java.io.*;

/**

   This class handles the generated propagation matrix data.
   It provides methods to return the particle flux differential
   dF(input particle, output particle), which is an element
   of the PropagationMatrix itself,
   for various combinations of input/output particles.

   This class has been primarily written for calculating
   the sensitivity/upper bound of EHE neutrino search
   with the neutrio flux model independent way.

   Written by S. Yoshida 2007/3/18
*/

public class PropagationMatrixFactory {

    static int dimension = Particle.getDimensionOfLogEnergyMatrix();
    /** For Glashow resonance NuETo(Mu & Tau flavor)'s added. **/
    double[][] FnuEToNuE,FnuMuToNuE,FnuTauToNuE,FmuToNuE,FtauToNuE;
    double[][] FnuEToNuMu,FnuMuToNuMu,FnuTauToNuMu,FmuToNuMu,FtauToNuMu;
    double[][] FnuEToNuTau,FnuMuToNuTau,FnuTauToNuTau,FmuToNuTau,FtauToNuTau;
    double[][] FnuEToE,FnuMuToE,FnuTauToE,FmuToE,FtauToE;
    double[][] FnuEToMu,FnuMuToMu,FnuTauToMu,FmuToMu,FtauToMu;
    double[][] FnuEToTau,FnuMuToTau,FnuTauToTau,FmuToTau,FtauToTau;
    double[][] FnuEToHadron,FnuMuToHadron,FnuTauToHadron,FmuToHadron,FtauToHadron;

    /** Flag to choose matrix with or without Glashow Resonance. */
    boolean includeGlashowResonance = true;

    /** Default constructor.
	Allocating the array memory to store the propagation matricis
     */
    public PropagationMatrixFactory(){
	generatePropagationMatrixArray();
    }

    /** This constructor is for subclass such as ExtendedMuonPropMatrixFactory */
    public PropagationMatrixFactory(boolean generateMatrix){
	if(generateMatrix) generatePropagationMatrixArray();
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
	/** For Glashow Resonance -begin **/
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
	Returns dF/dLogE*deltaLogE (inputParticle ---> outputParticle).

	The flavor/doublet/energy should be defined in each of the Particle object
        in the argument.

	delta LogE is the bin size of the propagation matrix. It is defined by
	Particle.getDimensionOfLogEnergyMatrix()
    */
    public double getDF(Particle inputParticle, Particle outputParticle){

	int inputFlavor = inputParticle.getFlavor();
	int inputDoublet = inputParticle.getDoublet();
	double logEinput = inputParticle.getLogEnergy();

	int outputFlavor = outputParticle.getFlavor();
	int outputDoublet = outputParticle.getDoublet();
	double logEoutput = outputParticle.getLogEnergy();


	return(getDF(inputFlavor,inputDoublet,logEinput,
		     outputFlavor,outputDoublet,logEoutput));
    }


    /** 
	Returns dF/dLogE * deltaLogE (inputParticle ---> outputParticle).

	<pre>
	int   inputFlavor   : flavor of the Particle object entering into the earth.
	int   inputDoublet  : doublet of the Particle object entering into the earth.
	double logEinput    : logE [GeV] of the particle entering into the earth.
	int   outputFlavor   : flavor of the Particle object after the propagation.
	int   outputDoublet  : doublet of the Particle object after the propagation.
	double logEoutput    : logE [GeV] of the particle after the propagation.
	</pre>

	delta LogE is the bin size of the propagation matrix. It is defined by
	Particle.getDimensionOfLogEnergyMatrix()
    */
    public double getDF(int inputFlavor, int inputDoublet, double logEinput,
			     int outputFlavor, int outputDoublet, double logEoutput){

	int iLogE = (int)((logEinput - Particle.getLogEnergyMinimum())
			  /Particle.getDeltaLogEnergy());
	int jLogE = (int)((logEoutput - Particle.getLogEnergyMinimum())
			  /Particle.getDeltaLogEnergy());

	if(((0 <= iLogE) && (iLogE<dimension)) && 
	   ((jLogE<= iLogE) && (0 <= jLogE))){ // In the valid energy range

	    double count = 0.0;

            if(inputFlavor==0 && inputDoublet == 0){ // nu_e
		/** For Glashow Resonance -begin **/
		if(outputFlavor==0 && outputDoublet ==0) count = FnuEToNuE[iLogE][jLogE];
                else if(outputFlavor==1 && outputDoublet ==0) count = FnuEToNuMu[iLogE][jLogE];
                else if(outputFlavor==2 && outputDoublet ==0) count = FnuEToNuTau[iLogE][jLogE];
                else if(outputFlavor==0 && outputDoublet ==1) count = FnuEToE[iLogE][jLogE];
                else if(outputFlavor==1 && outputDoublet ==1) count = FnuEToMu[iLogE][jLogE];
                else if(outputFlavor==2 && outputDoublet ==1) count = FnuEToTau[iLogE][jLogE];
                else if(outputFlavor==3) count = FnuEToHadron[iLogE][jLogE];
		/** For Glashow Resonance -end **/

	    }else if(inputFlavor==1 && inputDoublet ==0){ // nu_mu
                if(outputFlavor==0 && outputDoublet ==0) count = FnuMuToNuE[iLogE][jLogE];
                else if(outputFlavor==1 && outputDoublet ==0) count = FnuMuToNuMu[iLogE][jLogE];
                else if(outputFlavor==2 && outputDoublet ==0) count = FnuMuToNuTau[iLogE][jLogE];
                else if(outputFlavor==0 && outputDoublet ==1) count = FnuMuToE[iLogE][jLogE];
                else if(outputFlavor==1 && outputDoublet ==1) count = FnuMuToMu[iLogE][jLogE];
                else if(outputFlavor==2 && outputDoublet ==1) count = FnuMuToTau[iLogE][jLogE];
                else if(outputFlavor==3) count = FnuMuToHadron[iLogE][jLogE];

	    }else if(inputFlavor==2 && inputDoublet ==0){ // nu_tau
                if(outputFlavor==0 && outputDoublet ==0) count = FnuTauToNuE[iLogE][jLogE];
                else if(outputFlavor==1 && outputDoublet ==0) count = FnuTauToNuMu[iLogE][jLogE];
                else if(outputFlavor==2 && outputDoublet ==0) count = FnuTauToNuTau[iLogE][jLogE];
                else if(outputFlavor==0 && outputDoublet ==1) count = FnuTauToE[iLogE][jLogE];
                else if(outputFlavor==1 && outputDoublet ==1) count = FnuTauToMu[iLogE][jLogE];
                else if(outputFlavor==2 && outputDoublet ==1) count = FnuTauToTau[iLogE][jLogE];
                else if(outputFlavor==3) count = FnuTauToHadron[iLogE][jLogE];

	    }else if(inputFlavor==1 && inputDoublet ==1){ // mu
                if(outputFlavor==0 && outputDoublet ==0) count = FmuToNuE[iLogE][jLogE];
                else if(outputFlavor==1 && outputDoublet ==0) count = FmuToNuMu[iLogE][jLogE];
                else if(outputFlavor==2 && outputDoublet ==0) count = FmuToNuTau[iLogE][jLogE];
                else if(outputFlavor==0 && outputDoublet ==1) count = FmuToE[iLogE][jLogE];
                else if(outputFlavor==1 && outputDoublet ==1) count = FmuToMu[iLogE][jLogE];
                else if(outputFlavor==2 && outputDoublet ==1) count = FmuToTau[iLogE][jLogE];
                else if(outputFlavor==3) count = FmuToHadron[iLogE][jLogE];

	    }else if(inputFlavor==2 && inputDoublet ==1){ // tau
                if(outputFlavor==0 && outputDoublet ==0) count = FtauToNuE[iLogE][jLogE];
                else if(outputFlavor==1 && outputDoublet ==0) count = FtauToNuMu[iLogE][jLogE];
                else if(outputFlavor==2 && outputDoublet ==0) count = FtauToNuTau[iLogE][jLogE];
                else if(outputFlavor==0 && outputDoublet ==1) count = FtauToE[iLogE][jLogE];
                else if(outputFlavor==1 && outputDoublet ==1) count = FtauToMu[iLogE][jLogE];
                else if(outputFlavor==2 && outputDoublet ==1) count = FtauToTau[iLogE][jLogE];
                else if(outputFlavor==3) count = FtauToHadron[iLogE][jLogE];
	    }

	    return(count);

	}else{ // out of the vaid energy range
	    return(0.0);
	}

    }

    /** or Glashow Resonance */
    // Switch includeGlashowResonance flag.
    public void whetherPropagationMatrixWithGlashowResonance(boolean flag){
	includeGlashowResonance = flag;
    }


    /** Return average Epropagated/Eincoming  of propagating muon with Energy logEnergy [GeV] */
    public double getAverageMuonEnergyLossAfterPropagation(double logEnergy){

	int jLogE = (int)((logEnergy - Particle.getLogEnergyMinimum())
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

	double deltaLogE = averagedLogEprimary - logEnergy;
	double energyRatio = Math.pow(10.0,-deltaLogE);

	return energyRatio;
    }

}

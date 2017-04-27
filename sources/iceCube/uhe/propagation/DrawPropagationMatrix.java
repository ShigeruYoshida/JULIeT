package iceCube.uhe.propagation;

import iceCube.uhe.particles.*;
import iceCube.uhe.interactions.*;
import iceCube.uhe.points.*;
import iceCube.uhe.propagation.*;
import java.io.*;

/******
<pre>
       Draw the Prpgatation matrix to calculate the energy distribution
	and the flux of particles after propagation in the Earth. 

	The Incoming Nadir angle, the switches to control
	the interaction and decay channels involved, and the file name
	to record the calculated matrix are read thorugh the arguments.

inputParticle
bit7     bit6      bit5     bit4     bit3     bit2     bit1     bit0
SUM      Reserv.   tau      mu       Reserv.  nuTau    nuMu     nuE
outputParticle
bit7     bit6      bit5     bit4     bit3     bit2     bit1     bit0
Reserv.  hadron    tau      mu       e/gamma  nuTau    nuMu     nuE

</pre>

*/

public class DrawPropagationMatrix {

    int dimension = Particle.getDimensionOfLogEnergyMatrix();
    int inputParticle, outputParticle;
    int inLogE = -1;
    double logYmin,logYmax;

    /** Flag to choose matrix with or without Glashow Resonance. */
    boolean includeGlashowResonance = true;

    /** For Glashow Resonance, FnuETo(Mu/Tau flavor) are added. **/
    double[][] FnuEToNuE,FnuMuToNuE,FnuTauToNuE,FmuToNuE,FtauToNuE;
    double[][] FnuEToNuMu,FnuMuToNuMu,FnuTauToNuMu,FmuToNuMu,FtauToNuMu;
    double[][] FnuEToNuTau,FnuMuToNuTau,FnuTauToNuTau,FmuToNuTau,FtauToNuTau;
    double[][] FnuEToE,FnuMuToE,FnuTauToE,FmuToE,FtauToE;
    double[][] FnuEToMu,FnuMuToMu,FnuTauToMu,FmuToMu,FtauToMu;
    double[][] FnuEToTau,FnuMuToTau,FnuTauToTau,FmuToTau,FtauToTau;
    double[][] FnuEToHadron,FnuMuToHadron,FnuTauToHadron,FmuToHadron,FtauToHadron;
    private final static int NU_E = 1;
    private final static int NU_MU = 2;
    private final static int NU_TAU = 4;
    private final static int EL = 8;
    private final static int MU = 16;
    private final static int TAU = 32;
    private final static int HADRON = 64;
    private final static int SUM = 128;

    /** Constructor */
    public DrawPropagationMatrix(int inputParticle, int outputParticle,
				 double logYmin, double logYmax) throws IOException {

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

	this.inputParticle = inputParticle;
	this.outputParticle = outputParticle;
	this.logYmin = logYmin;
	this.logYmax = logYmax;

    }
    public DrawPropagationMatrix(int inputParticle, int outputParticle,
				 double logYmin, double logYmax, int inLogE) 
	throws IOException {

	this(inputParticle, outputParticle, logYmin, logYmax);

	if(0<=inLogE && inLogE<dimension){
	    this.inLogE = inLogE;
	}else{
	    System.err.println("inputLogE " + inLogE + " out of the matrix range");
	    System.exit(0);
	}

    }

    /** Read the calculated propagatin matrix */
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

    /** Draw the matrix data by gragig-Xfig */
    public void drawOnXfig( ){

	System.out.println("titx Log E[GeV]");
	System.out.println("tity Log Count");
	System.out.println("scal 6.0 12.0 " + logYmin + " " + logYmax);

	int iLogE;
	int jLogE;
	double count,logCount,logE;
	double[] sumF = new double[dimension];

	for(jLogE=0;jLogE<dimension;jLogE++) sumF[jLogE] = 0.0;

	if(((inputParticle & NU_E) == NU_E ) && ((outputParticle & NU_E) == NU_E) ){
            System.out.println("lnpt 0");
            System.out.println("lnth 1");
            System.out.println("lncl 0");
	    for(jLogE=0;jLogE<dimension;jLogE++){
		count = 0.0;
		if(inLogE==-1){
		    for(iLogE=jLogE;iLogE<dimension;iLogE++){
			if(FnuEToNuE[iLogE][jLogE]>0.0){
			    count += FnuEToNuE[iLogE][jLogE];
			    sumF[jLogE] += FnuEToNuE[iLogE][jLogE];
			}
		    }
		}else{
		    count += FnuEToNuE[inLogE][jLogE];
		    sumF[jLogE] += FnuEToNuE[inLogE][jLogE];
		}
		if(count<=0.0) logCount = logYmin;
		else logCount = Math.log(count)/Math.log(10.0);

		logE = Particle.getLogEnergyMinimum( ) + 
		    Particle.getDeltaLogEnergy( )*(double )jLogE;

		System.out.println("data " + logE + " 0.0 " + logCount + "0.0");

	    }
	    System.out.println("join");
	    System.out.println("disp");
	    System.out.println("cont");
	}

	/** For Glashow Resonance -begin **/
	if(((inputParticle & NU_E) == NU_E ) && ((outputParticle & NU_MU) == NU_MU) ){
            System.out.println("lnpt 0");
            System.out.println("lnth 1");
            System.out.println("lncl 0");
	    for(jLogE=0;jLogE<dimension;jLogE++){
		count = 0.0;
		if(inLogE==-1){
		    for(iLogE=jLogE;iLogE<dimension;iLogE++){
			if(FnuEToNuMu[iLogE][jLogE]>0.0){
			    count += FnuEToNuMu[iLogE][jLogE];
			    sumF[jLogE] += FnuEToNuMu[iLogE][jLogE];
			}
		    }
		}else{
		    count += FnuEToNuMu[inLogE][jLogE];
		    sumF[jLogE] += FnuEToNuMu[inLogE][jLogE];
		}
		if(count<=0.0) logCount = logYmin;
		else logCount = Math.log(count)/Math.log(10.0);

		logE = Particle.getLogEnergyMinimum( ) + 
		    Particle.getDeltaLogEnergy( )*(double )jLogE;

		System.out.println("data " + logE + " 0.0 " + logCount + "0.0");

	    }
	    System.out.println("join");
	    System.out.println("disp");
	    System.out.println("cont");
	}

	if(((inputParticle & NU_E) == NU_E ) && ((outputParticle & NU_TAU) == NU_TAU) ){
            System.out.println("lnpt 0");
            System.out.println("lnth 1");
            System.out.println("lncl 0");
	    for(jLogE=0;jLogE<dimension;jLogE++){
		count = 0.0;
		if(inLogE==-1){
		    for(iLogE=jLogE;iLogE<dimension;iLogE++){
			if(FnuEToNuTau[iLogE][jLogE]>0.0){
			    count += FnuEToNuTau[iLogE][jLogE];
			    sumF[jLogE] += FnuEToNuTau[iLogE][jLogE];
			}
		    }
		}else{
		    count += FnuEToNuTau[inLogE][jLogE];
		    sumF[jLogE] += FnuEToNuTau[inLogE][jLogE];
		}
		if(count<=0.0) logCount = logYmin;
		else logCount = Math.log(count)/Math.log(10.0);

		logE = Particle.getLogEnergyMinimum( ) + 
		    Particle.getDeltaLogEnergy( )*(double )jLogE;

		System.out.println("data " + logE + " 0.0 " + logCount + "0.0");

	    }
	    System.out.println("join");
	    System.out.println("disp");
	    System.out.println("cont");
	}
	/** For Glashow Resonance -end **/

	if(((inputParticle & NU_E) == NU_E ) && ((outputParticle & EL) == EL) ){
            System.out.println("lnpt 0");
            System.out.println("lnth 1");
            System.out.println("lncl 3");
	    for(jLogE=0;jLogE<dimension;jLogE++){
		count = 0.0;
		if(inLogE==-1){
		    for(iLogE=jLogE;iLogE<dimension;iLogE++){
			if(FnuEToE[iLogE][jLogE]>0.0){
			    count += FnuEToE[iLogE][jLogE];
			    sumF[jLogE] += FnuEToE[iLogE][jLogE];
			}
		    }
		}else{
		    count += FnuEToE[inLogE][jLogE];
		    sumF[jLogE] += FnuEToE[inLogE][jLogE];
		}
		if(count<=0.0) logCount = logYmin;
		else logCount = Math.log(count)/Math.log(10.0);

		logE = Particle.getLogEnergyMinimum( ) + 
		    Particle.getDeltaLogEnergy( )*(double )jLogE;

		System.out.println("data " + logE + " 0.0 " + logCount + "0.0");

	    }
	    System.out.println("join");
	    System.out.println("disp");
	    System.out.println("cont");
	}

	/** For Glashow Resonance -begin **/
	if(((inputParticle & NU_E) == NU_E ) && ((outputParticle & MU) == MU) ){
            System.out.println("lnpt 0");
            System.out.println("lnth 1");
            System.out.println("lncl 3");
	    for(jLogE=0;jLogE<dimension;jLogE++){
		count = 0.0;
		if(inLogE==-1){
		    for(iLogE=jLogE;iLogE<dimension;iLogE++){
			if(FnuEToMu[iLogE][jLogE]>0.0){
			    count += FnuEToMu[iLogE][jLogE];
			    sumF[jLogE] += FnuEToMu[iLogE][jLogE];
			}
		    }
		}else{
		    count += FnuEToMu[inLogE][jLogE];
		    sumF[jLogE] += FnuEToMu[inLogE][jLogE];
		}
		if(count<=0.0) logCount = logYmin;
		else logCount = Math.log(count)/Math.log(10.0);

		logE = Particle.getLogEnergyMinimum( ) + 
		    Particle.getDeltaLogEnergy( )*(double )jLogE;

		System.out.println("data " + logE + " 0.0 " + logCount + "0.0");

	    }
	    System.out.println("join");
	    System.out.println("disp");
	    System.out.println("cont");
	}

	if(((inputParticle & NU_E) == NU_E ) && ((outputParticle & TAU) == TAU) ){
            System.out.println("lnpt 0");
            System.out.println("lnth 1");
            System.out.println("lncl 3");
	    for(jLogE=0;jLogE<dimension;jLogE++){
		count = 0.0;
		if(inLogE==-1){
		    for(iLogE=jLogE;iLogE<dimension;iLogE++){
			if(FnuEToTau[iLogE][jLogE]>0.0){
			    count += FnuEToTau[iLogE][jLogE];
			    sumF[jLogE] += FnuEToTau[iLogE][jLogE];
			}
		    }
		}else{
		    count += FnuEToTau[inLogE][jLogE];
		    sumF[jLogE] += FnuEToTau[inLogE][jLogE];
		}
		if(count<=0.0) logCount = logYmin;
		else logCount = Math.log(count)/Math.log(10.0);

		logE = Particle.getLogEnergyMinimum( ) + 
		    Particle.getDeltaLogEnergy( )*(double )jLogE;

		System.out.println("data " + logE + " 0.0 " + logCount + "0.0");

	    }
	    System.out.println("join");
	    System.out.println("disp");
	    System.out.println("cont");
	}
	/** For Glashow Resonance -end **/

	if(((inputParticle & NU_E) == NU_E ) && ((outputParticle & HADRON) == HADRON) ){
            System.out.println("lnpt 0");
            System.out.println("lnth 1");
            System.out.println("lncl 6");
	    for(jLogE=0;jLogE<dimension;jLogE++){
		count = 0.0;
		if(inLogE==-1){
		    for(iLogE=jLogE;iLogE<dimension;iLogE++){
			if(FnuEToHadron[iLogE][jLogE]>0.0){
			    count += FnuEToHadron[iLogE][jLogE];
			    sumF[jLogE] += FnuEToHadron[iLogE][jLogE];
			}
		    }
		}else{
		    count += FnuEToHadron[inLogE][jLogE];
		    sumF[jLogE] += FnuEToHadron[inLogE][jLogE];
		}
		if(count<=0.0) logCount = logYmin;
		else logCount = Math.log(count)/Math.log(10.0);

		logE = Particle.getLogEnergyMinimum( ) + 
		    Particle.getDeltaLogEnergy( )*(double )jLogE;

		System.out.println("data " + logE + " 0.0 " + logCount + "0.0");

	    }
	    System.out.println("join");
	    System.out.println("disp");
	    System.out.println("cont");
	}



	if(((inputParticle & NU_MU) == NU_MU ) && ((outputParticle & NU_E) == NU_E) ){
            System.out.println("lnpt 1");
            System.out.println("lnth 1");
            System.out.println("lncl 0");
	    for(jLogE=0;jLogE<dimension;jLogE++){
		count = 0.0;
		if(inLogE==-1){
		    for(iLogE=jLogE;iLogE<dimension;iLogE++){
			if(FnuMuToNuE[iLogE][jLogE]>0.0){
			    count += FnuMuToNuE[iLogE][jLogE];
			    sumF[jLogE] += FnuMuToNuE[iLogE][jLogE];
			}
		    }
		}else{
		    count += FnuMuToNuE[inLogE][jLogE];
		    sumF[jLogE] += FnuMuToNuE[inLogE][jLogE];
		}
		if(count<=0.0) logCount = logYmin;
		else logCount = Math.log(count)/Math.log(10.0);

		logE = Particle.getLogEnergyMinimum( ) + 
		    Particle.getDeltaLogEnergy( )*(double )jLogE;

		System.out.println("data " + logE + " 0.0 " + logCount + "0.0");

	    }
	    System.out.println("join");
	    System.out.println("disp");
	    System.out.println("cont");
	}

	if(((inputParticle & NU_MU) == NU_MU ) && ((outputParticle & NU_MU) == NU_MU) ){
            System.out.println("lnpt 1");
            System.out.println("lnth 1");
            System.out.println("lncl 1");  
	    for(jLogE=0;jLogE<dimension;jLogE++){
		count = 0.0;
		if(inLogE==-1){
		    for(iLogE=jLogE;iLogE<dimension;iLogE++){
			if(FnuMuToNuMu[iLogE][jLogE]>0.0){
			    count += FnuMuToNuMu[iLogE][jLogE];
			    sumF[jLogE] += FnuMuToNuMu[iLogE][jLogE];
			}
		    }
		}else{
		    count += FnuMuToNuMu[inLogE][jLogE];
		    sumF[jLogE] += FnuMuToNuMu[inLogE][jLogE];
		}
		if(count<=0.0) logCount = logYmin;
		else logCount = Math.log(count)/Math.log(10.0);

		logE = Particle.getLogEnergyMinimum( ) + 
		    Particle.getDeltaLogEnergy( )*(double )jLogE;

		System.out.println("data " + logE + " 0.0 " + logCount + "0.0");

	    }
	    System.out.println("join");
	    System.out.println("disp");
	    System.out.println("cont");
	}

	if(((inputParticle & NU_MU) == NU_MU ) && ((outputParticle & NU_TAU) == NU_TAU) ){
            System.out.println("lnpt 1");
            System.out.println("lnth 1");
            System.out.println("lncl 2");
	    for(jLogE=0;jLogE<dimension;jLogE++){
		count = 0.0;
		if(inLogE==-1){
		    for(iLogE=jLogE;iLogE<dimension;iLogE++){
			if(FnuMuToNuTau[iLogE][jLogE]>0.0){
			    count += FnuMuToNuTau[iLogE][jLogE];
			    sumF[jLogE] += FnuMuToNuTau[iLogE][jLogE];
			}
		    }
		}else{
		    count += FnuMuToNuTau[inLogE][jLogE];
		    sumF[jLogE] += FnuMuToNuTau[inLogE][jLogE];
		}
		if(count<=0.0) logCount = logYmin;
		else logCount = Math.log(count)/Math.log(10.0);

		logE = Particle.getLogEnergyMinimum( ) + 
		    Particle.getDeltaLogEnergy( )*(double )jLogE;

		System.out.println("data " + logE + " 0.0 " + logCount + "0.0");

	    }
	    System.out.println("join");
	    System.out.println("disp");
	    System.out.println("cont");
	}

	if(((inputParticle & NU_MU) == NU_MU ) && ((outputParticle & EL) == EL) ){
            System.out.println("lnpt 1");
            System.out.println("lnth 1");
            System.out.println("lncl 3");
	    for(jLogE=0;jLogE<dimension;jLogE++){
		count = 0.0;
		if(inLogE==-1){
		    for(iLogE=jLogE;iLogE<dimension;iLogE++){
			if(FnuMuToE[iLogE][jLogE]>0.0){
			    count += FnuMuToE[iLogE][jLogE];
			    sumF[jLogE] += FnuMuToE[iLogE][jLogE];
			}
		    }
		}else{
		    count += FnuMuToE[inLogE][jLogE];
		    sumF[jLogE] += FnuMuToE[inLogE][jLogE];
		}
		if(count<=0.0) logCount = logYmin;
		else logCount = Math.log(count)/Math.log(10.0);

		logE = Particle.getLogEnergyMinimum( ) + 
		    Particle.getDeltaLogEnergy( )*(double )jLogE;

		System.out.println("data " + logE + " 0.0 " + logCount + "0.0");

	    }
	    System.out.println("join");
	    System.out.println("disp");
	    System.out.println("cont");
	}

	if(((inputParticle & NU_MU) == NU_MU ) && ((outputParticle & MU) == MU) ){ 
            System.out.println("lnpt 1");
            System.out.println("lnth 1");
            System.out.println("lncl 4"); 
	    for(jLogE=0;jLogE<dimension;jLogE++){
		count = 0.0;
		if(inLogE==-1){
		    for(iLogE=jLogE;iLogE<dimension;iLogE++){
			if(FnuMuToMu[iLogE][jLogE]>0.0){
			    count += FnuMuToMu[iLogE][jLogE];
			    sumF[jLogE] += FnuMuToMu[iLogE][jLogE];
			}
		    }
		}else{
		    count += FnuMuToMu[inLogE][jLogE];
		    sumF[jLogE] += FnuMuToMu[inLogE][jLogE];
		}
		if(count<=0.0) logCount = logYmin;
		else logCount = Math.log(count)/Math.log(10.0);

		logE = Particle.getLogEnergyMinimum( ) + 
		    Particle.getDeltaLogEnergy( )*(double )jLogE;

        	System.out.println("data " + logE + " 0.0 " + logCount + "0.0");

	    }
	    System.out.println("join");
	    System.out.println("disp");
	    System.out.println("cont");
	}

	if(((inputParticle & NU_MU) == NU_MU ) && ((outputParticle & TAU) == TAU) ){
            System.out.println("lnpt 1");
            System.out.println("lnth 1");
            System.out.println("lncl 5");
	    for(jLogE=0;jLogE<dimension;jLogE++){
		count = 0.0;
		if(inLogE==-1){
		    for(iLogE=jLogE;iLogE<dimension;iLogE++){
			if(FnuMuToTau[iLogE][jLogE]>0.0){
			    count += FnuMuToTau[iLogE][jLogE];
			    sumF[jLogE] += FnuMuToTau[iLogE][jLogE];
			}
		    }
		}else{
		    count += FnuMuToTau[inLogE][jLogE];
		    sumF[jLogE] += FnuMuToTau[inLogE][jLogE];
		}
		if(count<=0.0) logCount = logYmin;
		else logCount = Math.log(count)/Math.log(10.0);

		logE = Particle.getLogEnergyMinimum( ) + 
		    Particle.getDeltaLogEnergy( )*(double )jLogE;

	      	System.out.println("data " + logE + " 0.0 " + logCount + "0.0");

	    }
	    System.out.println("join");
	    System.out.println("disp");
	    System.out.println("cont");
	}

	if(((inputParticle & NU_MU) == NU_MU ) && ((outputParticle & HADRON) == HADRON) ){
            System.out.println("lnpt 1");
            System.out.println("lnth 1");
            System.out.println("lncl 6");
	    for(jLogE=0;jLogE<dimension;jLogE++){
		count = 0.0;
		if(inLogE==-1){
		    for(iLogE=jLogE;iLogE<dimension;iLogE++){
			if(FnuMuToHadron[iLogE][jLogE]>0.0){
			    count += FnuMuToHadron[iLogE][jLogE];
			    sumF[jLogE] += FnuMuToHadron[iLogE][jLogE];
			}
		    }
		}else{
		    count += FnuMuToHadron[inLogE][jLogE];
		    sumF[jLogE] += FnuMuToHadron[inLogE][jLogE];
		}
		if(count<=0.0) logCount = logYmin;
		else logCount = Math.log(count)/Math.log(10.0);

		logE = Particle.getLogEnergyMinimum( ) + 
		    Particle.getDeltaLogEnergy( )*(double )jLogE;

		System.out.println("data " + logE + " 0.0 " + logCount + "0.0");

	    }
	    System.out.println("join");
	    System.out.println("disp");
	    System.out.println("cont");
	}




	if(((inputParticle & NU_TAU) == NU_TAU ) && ((outputParticle & NU_E) == NU_E) ){
            System.out.println("lnpt 2");
            System.out.println("lnth 1");
            System.out.println("lncl 0"); 
	    for(jLogE=0;jLogE<dimension;jLogE++){
		count = 0.0;
		if(inLogE==-1){
		    for(iLogE=jLogE;iLogE<dimension;iLogE++){
			if(FnuTauToNuE[iLogE][jLogE]>0.0){
			    count += FnuTauToNuE[iLogE][jLogE];
			    sumF[jLogE] += FnuTauToNuE[iLogE][jLogE];
			}
		    }
		}else{
		    count += FnuTauToNuE[inLogE][jLogE];
		    sumF[jLogE] += FnuTauToNuE[inLogE][jLogE];
		}
		if(count<=0.0) logCount = logYmin;
		else logCount = Math.log(count)/Math.log(10.0);

		logE = Particle.getLogEnergyMinimum( ) + 
		    Particle.getDeltaLogEnergy( )*(double )jLogE;

		System.out.println("data " + logE + " 0.0 " + logCount + "0.0");

	    }
	    System.out.println("join");
	    System.out.println("disp");
	    System.out.println("cont");
	}

	if(((inputParticle & NU_TAU) == NU_TAU ) && ((outputParticle & NU_MU) == NU_MU) ){
            System.out.println("lnpt 2");
            System.out.println("lnth 1");
            System.out.println("lncl 1");
	    for(jLogE=0;jLogE<dimension;jLogE++){
		count = 0.0;
		if(inLogE==-1){
		    for(iLogE=jLogE;iLogE<dimension;iLogE++){
			if(FnuTauToNuMu[iLogE][jLogE]>0.0){
			    count += FnuTauToNuMu[iLogE][jLogE];
			    sumF[jLogE] += FnuTauToNuMu[iLogE][jLogE];
			}
		    }
		}else{
		    count += FnuTauToNuMu[inLogE][jLogE];
		    sumF[jLogE] += FnuTauToNuMu[inLogE][jLogE];
		}
		if(count<=0.0) logCount = logYmin;
		else logCount = Math.log(count)/Math.log(10.0);

		logE = Particle.getLogEnergyMinimum( ) + 
		    Particle.getDeltaLogEnergy( )*(double )jLogE;

		System.out.println("data " + logE + " 0.0 " + logCount + "0.0");

	    }
	    System.out.println("join");
	    System.out.println("disp");
	    System.out.println("cont");
	}

	if(((inputParticle & NU_TAU) == NU_TAU ) && ((outputParticle & NU_TAU) == NU_TAU) ){
            System.out.println("lnpt 2");
            System.out.println("lnth 1");
            System.out.println("lncl 2");
	    for(jLogE=0;jLogE<dimension;jLogE++){
		count = 0.0;
		if(inLogE==-1){
		    for(iLogE=jLogE;iLogE<dimension;iLogE++){
			if(FnuTauToNuTau[iLogE][jLogE]>0.0){
			    count += FnuTauToNuTau[iLogE][jLogE];
			    sumF[jLogE] += FnuTauToNuTau[iLogE][jLogE];
			}
		    }
		}else{
		    count += FnuTauToNuTau[inLogE][jLogE];
		    sumF[jLogE] += FnuTauToNuTau[inLogE][jLogE];
		}
		if(count<=0.0) logCount = logYmin;
		else logCount = Math.log(count)/Math.log(10.0);

		logE = Particle.getLogEnergyMinimum( ) + 
		    Particle.getDeltaLogEnergy( )*(double )jLogE;

		System.out.println("data " + logE + " 0.0 " + logCount + "0.0");

	    }
	    System.out.println("join");
	    System.out.println("disp");
	    System.out.println("cont");
	}

	if(((inputParticle & NU_TAU) == NU_TAU ) && ((outputParticle & EL) == EL) ){
	    System.out.println("lnpt 2");
            System.out.println("lnth 1");
            System.out.println("lncl 3");
	    for(jLogE=0;jLogE<dimension;jLogE++){
		count = 0.0;
		if(inLogE==-1){
		    for(iLogE=jLogE;iLogE<dimension;iLogE++){
			if(FnuTauToE[iLogE][jLogE]>0.0){
			    count += FnuTauToE[iLogE][jLogE];
			    sumF[jLogE] += FnuTauToE[iLogE][jLogE];
			}
		    }
		}else{
		    count += FnuTauToE[inLogE][jLogE];
		    sumF[jLogE] += FnuTauToE[inLogE][jLogE];
		}
		if(count<=0.0) logCount = logYmin;
		else logCount = Math.log(count)/Math.log(10.0);

		logE = Particle.getLogEnergyMinimum( ) + 
		    Particle.getDeltaLogEnergy( )*(double )jLogE;

		System.out.println("data " + logE + " 0.0 " + logCount + "0.0");

	    }
	    System.out.println("join");
	    System.out.println("disp");
	    System.out.println("cont");
	}

	if(((inputParticle & NU_TAU) == NU_TAU ) && ((outputParticle & MU) == MU) ){
            System.out.println("lnpt 2");
            System.out.println("lnth 1");
            System.out.println("lncl 4");
	    for(jLogE=0;jLogE<dimension;jLogE++){
		count = 0.0;
		if(inLogE==-1){
		    for(iLogE=jLogE;iLogE<dimension;iLogE++){
			if(FnuTauToMu[iLogE][jLogE]>0.0){
			    count += FnuTauToMu[iLogE][jLogE];
			    sumF[jLogE] += FnuTauToMu[iLogE][jLogE];
			}
		    }
		}else{
		    count += FnuTauToMu[inLogE][jLogE];
		    sumF[jLogE] += FnuTauToMu[inLogE][jLogE];
		}
		if(count<=0.0) logCount = logYmin;
		else logCount = Math.log(count)/Math.log(10.0);

		logE = Particle.getLogEnergyMinimum( ) + 
		    Particle.getDeltaLogEnergy( )*(double )jLogE;

	      	System.out.println("data " + logE + " 0.0 " + logCount + "0.0");

	    }
	    System.out.println("join");
	    System.out.println("disp");
	    System.out.println("cont");
	}

	if(((inputParticle & NU_TAU) == NU_TAU ) && ((outputParticle & TAU) == TAU) ){
            System.out.println("lnpt 2");
            System.out.println("lnth 1");
            System.out.println("lncl 5");
	    for(jLogE=0;jLogE<dimension;jLogE++){
		count = 0.0;
		if(inLogE==-1){
		    for(iLogE=jLogE;iLogE<dimension;iLogE++){
			if(FnuTauToTau[iLogE][jLogE]>0.0){
			    count += FnuTauToTau[iLogE][jLogE];
			    sumF[jLogE] += FnuTauToTau[iLogE][jLogE];
			}
		    }
		}else{
		    count += FnuTauToTau[inLogE][jLogE];
		    sumF[jLogE] += FnuTauToTau[inLogE][jLogE];
		}
		if(count<=0.0) logCount = logYmin;
		else logCount = Math.log(count)/Math.log(10.0);

		logE = Particle.getLogEnergyMinimum( ) + 
		    Particle.getDeltaLogEnergy( )*(double )jLogE;

	       	System.out.println("data " + logE + " 0.0 " + logCount + "0.0");

	    }
	    System.out.println("join");
	    System.out.println("disp");
	    System.out.println("cont");
	}

	if(((inputParticle & NU_TAU) == NU_TAU ) && ((outputParticle & HADRON) == HADRON) ){
            System.out.println("lnpt 2");
            System.out.println("lnth 1");
            System.out.println("lncl 6");
	    for(jLogE=0;jLogE<dimension;jLogE++){
		count = 0.0;
		if(inLogE==-1){
		    for(iLogE=jLogE;iLogE<dimension;iLogE++){
			if(FnuTauToHadron[iLogE][jLogE]>0.0){
			    count += FnuTauToHadron[iLogE][jLogE];
			    sumF[jLogE] += FnuTauToHadron[iLogE][jLogE];
			}
		    }
		}else{
		    count += FnuTauToHadron[inLogE][jLogE];
		    sumF[jLogE] += FnuTauToHadron[inLogE][jLogE];
		}
		if(count<=0.0) logCount = logYmin;
		else logCount = Math.log(count)/Math.log(10.0);

		logE = Particle.getLogEnergyMinimum( ) + 
		    Particle.getDeltaLogEnergy( )*(double )jLogE;

		System.out.println("data " + logE + " 0.0 " + logCount + "0.0");

	    }
	    System.out.println("join");
	    System.out.println("disp");
	    System.out.println("cont");
	}




	if(((inputParticle & MU) == MU ) && ((outputParticle & NU_E) == NU_E) ){
            System.out.println("lnpt 1");
            System.out.println("lnth 3");
            System.out.println("lncl 0");
	    for(jLogE=0;jLogE<dimension;jLogE++){
		count = 0.0;
		if(inLogE==-1){
		    for(iLogE=jLogE;iLogE<dimension;iLogE++){
			if(FmuToNuE[iLogE][jLogE]>0.0){
			    count += FmuToNuE[iLogE][jLogE];
			    sumF[jLogE] += FmuToNuE[iLogE][jLogE];
			}
		    }
		}else{
		    count += FmuToNuE[inLogE][jLogE];
		    sumF[jLogE] += FmuToNuE[inLogE][jLogE];
		}
		if(count<=0.0) logCount = logYmin;
		else logCount = Math.log(count)/Math.log(10.0);

		logE = Particle.getLogEnergyMinimum( ) + 
		    Particle.getDeltaLogEnergy( )*(double )jLogE;

		System.out.println("data " + logE + " 0.0 " + logCount + "0.0");

	    }
	    System.out.println("join");
	    System.out.println("disp");
	    System.out.println("cont");
	}

	if(((inputParticle & MU) == MU ) && ((outputParticle & NU_MU) == NU_MU) ){
            System.out.println("lnpt 1");
            System.out.println("lnth 3");
            System.out.println("lncl 1");
	    for(jLogE=0;jLogE<dimension;jLogE++){
		count = 0.0;
		if(inLogE==-1){
		    for(iLogE=jLogE;iLogE<dimension;iLogE++){
			if(FmuToNuMu[iLogE][jLogE]>0.0){
			    count += FmuToNuMu[iLogE][jLogE];
			    sumF[jLogE] += FmuToNuMu[iLogE][jLogE];
			}
		    }
		}else{
		    count += FmuToNuMu[inLogE][jLogE];
		    sumF[jLogE] += FmuToNuMu[inLogE][jLogE];
		}
		if(count<=0.0) logCount = logYmin;
		else logCount = Math.log(count)/Math.log(10.0);

		logE = Particle.getLogEnergyMinimum( ) + 
		    Particle.getDeltaLogEnergy( )*(double )jLogE;

		System.out.println("data " + logE + " 0.0 " + logCount + "0.0");

	    }
	    System.out.println("join");
	    System.out.println("disp");
	    System.out.println("cont");
	}

	if(((inputParticle & MU) == MU ) && ((outputParticle & NU_TAU) == NU_TAU) ){
            System.out.println("lnpt 1");
            System.out.println("lnth 3");
            System.out.println("lncl 2");
	    for(jLogE=0;jLogE<dimension;jLogE++){
		count = 0.0;
		if(inLogE==-1){
		    for(iLogE=jLogE;iLogE<dimension;iLogE++){
			if(FmuToNuTau[iLogE][jLogE]>0.0){
			    count += FmuToNuTau[iLogE][jLogE];
			    sumF[jLogE] += FmuToNuTau[iLogE][jLogE];
			}
		    }
		}else{
		    count += FmuToNuTau[inLogE][jLogE];
		    sumF[jLogE] += FmuToNuTau[inLogE][jLogE];
		}
		if(count<=0.0) logCount = logYmin;
		else logCount = Math.log(count)/Math.log(10.0);

		logE = Particle.getLogEnergyMinimum( ) + 
		    Particle.getDeltaLogEnergy( )*(double )jLogE;

		System.out.println("data " + logE + " 0.0 " + logCount + "0.0");

	    }
	    System.out.println("join");
	    System.out.println("disp");
	    System.out.println("cont");
	}

	if(((inputParticle & MU) == MU ) && ((outputParticle & EL) == EL) ){
            System.out.println("lnpt 1");
            System.out.println("lnth 3");
            System.out.println("lncl 3");
	    for(jLogE=0;jLogE<dimension;jLogE++){
		count = 0.0;
		if(inLogE==-1){
		    for(iLogE=jLogE;iLogE<dimension;iLogE++){
			if(FmuToE[iLogE][jLogE]>0.0){
			    count += FmuToE[iLogE][jLogE];
			    sumF[jLogE] += FmuToE[iLogE][jLogE];
			}
		    }
		}else{
		    count += FmuToE[inLogE][jLogE];
		    sumF[jLogE] += FmuToE[inLogE][jLogE];
		}
		if(count<=0.0) logCount = logYmin;
		else logCount = Math.log(count)/Math.log(10.0);

		logE = Particle.getLogEnergyMinimum( ) + 
		    Particle.getDeltaLogEnergy( )*(double )jLogE;

		System.out.println("data " + logE + " 0.0 " + logCount + "0.0");

	    }
	    System.out.println("join");
	    System.out.println("disp");
	    System.out.println("cont");
	}

	if(((inputParticle & MU) == MU ) && ((outputParticle & MU) == MU) ){               System.out.println("lnpt 1");
            System.out.println("lnth 3");
            System.out.println("lncl 4");
	    for(jLogE=0;jLogE<dimension;jLogE++){
		count = 0.0;
		if(inLogE==-1){
		    for(iLogE=jLogE;iLogE<dimension;iLogE++){
			if(FmuToMu[iLogE][jLogE]>0.0){
			    count += FmuToMu[iLogE][jLogE];
			    sumF[jLogE] += FmuToMu[iLogE][jLogE];
			}
		    }
		}else{
		    count += FmuToMu[inLogE][jLogE];
		    sumF[jLogE] += FmuToMu[inLogE][jLogE];
		}
		if(count<=0.0) logCount = logYmin;
		else logCount = Math.log(count)/Math.log(10.0);

		logE = Particle.getLogEnergyMinimum( ) + 
		    Particle.getDeltaLogEnergy( )*(double )jLogE;

		System.out.println("data " + logE + " 0.0 " + logCount + "0.0");

	    }
	    System.out.println("join");
	    System.out.println("disp");
	    System.out.println("cont");
	}

	if(((inputParticle & MU) == MU ) && ((outputParticle & TAU) == TAU) ){
            System.out.println("lnpt 1");
            System.out.println("lnth 3");
            System.out.println("lncl 5");
	    for(jLogE=0;jLogE<dimension;jLogE++){
		count = 0.0;
		if(inLogE==-1){
		    for(iLogE=jLogE;iLogE<dimension;iLogE++){
			if(FmuToTau[iLogE][jLogE]>0.0){
			    count += FmuToTau[iLogE][jLogE];
			    sumF[jLogE] += FmuToTau[iLogE][jLogE];
			}
		    }
		}else{
		    count += FmuToTau[inLogE][jLogE];
		    sumF[jLogE] += FmuToTau[inLogE][jLogE];
		}
		if(count<=0.0) logCount = logYmin;
		else logCount = Math.log(count)/Math.log(10.0);

		logE = Particle.getLogEnergyMinimum( ) + 
		    Particle.getDeltaLogEnergy( )*(double )jLogE;

		System.out.println("data " + logE + " 0.0 " + logCount + "0.0");

	    }
	    System.out.println("join");
	    System.out.println("disp");
	    System.out.println("cont");
	}

	if(((inputParticle & MU) == MU ) && ((outputParticle & HADRON) == HADRON) ){
            System.out.println("lnpt 1");
            System.out.println("lnth 3");
            System.out.println("lncl 6");
	    for(jLogE=0;jLogE<dimension;jLogE++){
		count = 0.0;
		if(inLogE==-1){
		    for(iLogE=jLogE;iLogE<dimension;iLogE++){
			if(FmuToHadron[iLogE][jLogE]>0.0){
			    count += FmuToHadron[iLogE][jLogE];
			    sumF[jLogE] += FmuToHadron[iLogE][jLogE];
			}
		    }
		}else{
		    count += FmuToHadron[inLogE][jLogE];
		    sumF[jLogE] += FmuToHadron[inLogE][jLogE];
		}
		if(count<=0.0) logCount = logYmin;
		else logCount = Math.log(count)/Math.log(10.0);

		logE = Particle.getLogEnergyMinimum( ) + 
		    Particle.getDeltaLogEnergy( )*(double )jLogE;

		System.out.println("data " + logE + " 0.0 " + logCount + "0.0");

	    }
	    System.out.println("join");
	    System.out.println("disp");
	    System.out.println("cont");
	}




	if(((inputParticle & TAU) == TAU ) && ((outputParticle & NU_E) == NU_E) ){
            System.out.println("lnpt 2");
            System.out.println("lnth 3");
            System.out.println("lncl 0");
	    for(jLogE=0;jLogE<dimension;jLogE++){
		count = 0.0;
		if(inLogE==-1){
		    for(iLogE=jLogE;iLogE<dimension;iLogE++){
			if(FtauToNuE[iLogE][jLogE]>0.0){
			    count += FtauToNuE[iLogE][jLogE];
			    sumF[jLogE] += FtauToNuE[iLogE][jLogE];
			}
		    }
		}else{
		    count += FtauToNuE[inLogE][jLogE];
		    sumF[jLogE] += FtauToNuE[inLogE][jLogE];
		}
		if(count<=0.0) logCount = logYmin;
		else logCount = Math.log(count)/Math.log(10.0);

		logE = Particle.getLogEnergyMinimum( ) + 
		    Particle.getDeltaLogEnergy( )*(double )jLogE;

		System.out.println("data " + logE + " 0.0 " + logCount + "0.0");

	    }
	    System.out.println("join");
	    System.out.println("disp");
	    System.out.println("cont");
	}

	if(((inputParticle & TAU) == TAU ) && ((outputParticle & NU_MU) == NU_MU) ){
            System.out.println("lnpt 2");
            System.out.println("lnth 3");
            System.out.println("lncl 1");
	    for(jLogE=0;jLogE<dimension;jLogE++){
		count = 0.0;
		if(inLogE==-1){
		    for(iLogE=jLogE;iLogE<dimension;iLogE++){
			if(FtauToNuMu[iLogE][jLogE]>0.0){
			    count += FtauToNuMu[iLogE][jLogE];
			    sumF[jLogE] += FtauToNuMu[iLogE][jLogE];
			}
		    }
		}else{
		    count += FtauToNuMu[inLogE][jLogE];
		    sumF[jLogE] += FtauToNuMu[inLogE][jLogE];
		}
		if(count<=0.0) logCount = logYmin;
		else logCount = Math.log(count)/Math.log(10.0);

		logE = Particle.getLogEnergyMinimum( ) + 
		    Particle.getDeltaLogEnergy( )*(double )jLogE;

		System.out.println("data " + logE + " 0.0 " + logCount + "0.0");

	    }
	    System.out.println("join");
	    System.out.println("disp");
	    System.out.println("cont");
	}

	if(((inputParticle & TAU) == TAU ) && ((outputParticle & NU_TAU) == NU_TAU) ){
            System.out.println("lnpt 2");
            System.out.println("lnth 3");
            System.out.println("lncl 2");
	    for(jLogE=0;jLogE<dimension;jLogE++){
		count = 0.0;
		if(inLogE==-1){
		    for(iLogE=jLogE;iLogE<dimension;iLogE++){
			if(FtauToNuTau[iLogE][jLogE]>0.0){
			    count += FtauToNuTau[iLogE][jLogE];
			    sumF[jLogE] += FtauToNuTau[iLogE][jLogE];
			}
		    }
		}else{
		    count += FtauToNuTau[inLogE][jLogE];
		    sumF[jLogE] += FtauToNuTau[inLogE][jLogE];
		}
		if(count<=0.0) logCount = logYmin;
		else logCount = Math.log(count)/Math.log(10.0);

		logE = Particle.getLogEnergyMinimum( ) + 
		    Particle.getDeltaLogEnergy( )*(double )jLogE;

		System.out.println("data " + logE + " 0.0 " + logCount + "0.0");

	    }
	    System.out.println("join");
	    System.out.println("disp");
	    System.out.println("cont");
	}

	if(((inputParticle & TAU) == TAU ) && ((outputParticle & EL) == EL) ){
            System.out.println("lnpt 2");
            System.out.println("lnth 3");
            System.out.println("lncl 3");
	    for(jLogE=0;jLogE<dimension;jLogE++){
		count = 0.0;
		if(inLogE==-1){
		    for(iLogE=jLogE;iLogE<dimension;iLogE++){
			if(FtauToE[iLogE][jLogE]>0.0){
			    count += FtauToE[iLogE][jLogE];
			    sumF[jLogE] += FtauToE[iLogE][jLogE];
			}
		    }
		}else{
		    count += FtauToE[inLogE][jLogE];
		    sumF[jLogE] += FtauToE[inLogE][jLogE];
		}
		if(count<=0.0) logCount = logYmin;
		else logCount = Math.log(count)/Math.log(10.0);

		logE = Particle.getLogEnergyMinimum( ) + 
		    Particle.getDeltaLogEnergy( )*(double )jLogE;

		System.out.println("data " + logE + " 0.0 " + logCount + "0.0");

	    }
	    System.out.println("join");
	    System.out.println("disp");
	    System.out.println("cont");
	}

	if(((inputParticle & TAU) == TAU ) && ((outputParticle & MU) == MU) ){
            System.out.println("lnpt 2");
            System.out.println("lnth 3");
            System.out.println("lncl 4");
	    for(jLogE=0;jLogE<dimension;jLogE++){
		count = 0.0;
		if(inLogE==-1){
		    for(iLogE=jLogE;iLogE<dimension;iLogE++){
			if(FtauToMu[iLogE][jLogE]>0.0){
			    count += FtauToMu[iLogE][jLogE];
			    sumF[jLogE] += FtauToMu[iLogE][jLogE];
			}
		    }
		}else{
		    count += FtauToMu[inLogE][jLogE];
		    sumF[jLogE] += FtauToMu[inLogE][jLogE];
		}
		if(count<=0.0) logCount = logYmin;
		else logCount = Math.log(count)/Math.log(10.0);

		logE = Particle.getLogEnergyMinimum( ) + 
		    Particle.getDeltaLogEnergy( )*(double )jLogE;

		System.out.println("data " + logE + " 0.0 " + logCount + "0.0");

	    }
	    System.out.println("join");
	    System.out.println("disp");
	    System.out.println("cont");
	}

	if(((inputParticle & TAU) == TAU ) && ((outputParticle & TAU) == TAU) ){
            System.out.println("lnpt 2");
            System.out.println("lnth 3");
            System.out.println("lncl 5");
	    for(jLogE=0;jLogE<dimension;jLogE++){
		count = 0.0;
		if(inLogE==-1){
		    for(iLogE=jLogE;iLogE<dimension;iLogE++){
			if(FtauToTau[iLogE][jLogE]>0.0){
			    count += FtauToTau[iLogE][jLogE];
			    sumF[jLogE] += FtauToTau[iLogE][jLogE];
			}
		    }
		}else{
		    count += FtauToTau[inLogE][jLogE];
		    sumF[jLogE] += FtauToTau[inLogE][jLogE];
		}
		if(count<=0.0) logCount = logYmin;
		else logCount = Math.log(count)/Math.log(10.0);

		logE = Particle.getLogEnergyMinimum( ) + 
		    Particle.getDeltaLogEnergy( )*(double )jLogE;

                System.out.println("data " + logE + " 0.0 " + logCount + "0.0");

	    }
	    System.out.println("join");
	    System.out.println("disp");
	    System.out.println("cont");
	}

	if(((inputParticle & TAU) == TAU ) && ((outputParticle & HADRON) == HADRON) ){
            System.out.println("lnpt 2");
            System.out.println("lnth 3");
            System.out.println("lncl 6");
	    for(jLogE=0;jLogE<dimension;jLogE++){
		count = 0.0;
		if(inLogE==-1){
		    for(iLogE=jLogE;iLogE<dimension;iLogE++){
			if(FtauToHadron[iLogE][jLogE]>0.0){
			    count += FtauToHadron[iLogE][jLogE];
			    sumF[jLogE] += FtauToHadron[iLogE][jLogE];
			}
		    }
		}else{
		    count += FtauToHadron[inLogE][jLogE];
		    sumF[jLogE] += FtauToHadron[inLogE][jLogE];
		}
		if(count<=0.0) logCount = logYmin;
		else logCount = Math.log(count)/Math.log(10.0);

		logE = Particle.getLogEnergyMinimum( ) + 
		    Particle.getDeltaLogEnergy( )*(double )jLogE;

		System.out.println("data " + logE + " 0.0 " + logCount + "0.0");

	    }
	    System.out.println("join");
	    System.out.println("disp");
	    System.out.println("cont");
	}


	// Summation of the all fluxes you selected.
        System.out.println("lnpt 0");
        System.out.println("lnth 3");
        System.out.println("lncl 0");
	for(jLogE=0;jLogE<dimension;jLogE++){
	    count = sumF[jLogE];
	    if(count<=0.0) logCount = logYmin;
	    else logCount = Math.log(count)/Math.log(10.0);

	    logE = Particle.getLogEnergyMinimum( ) + 
		Particle.getDeltaLogEnergy( )*(double )jLogE;

	    System.out.println("data " + logE + " 0.0 " + logCount + "0.0");
	}

	System.out.println("join");
	System.out.println("disp");
	System.out.println("endg");
    }


    public static void main(String[] args) throws IOException {

        String fileName = null;
	String matrixName = null;
	int inP = 255;
	int outP = 255;
	int inLogE = 100;
	double logYmin = -12.0;
	double logYmax = 1.0;
	DrawPropagationMatrix draw;


        if(args.length<5){
            System.out.println(
   "Usage: DrawPropagationMatrix file-name inputP outputP logYmin logYmax (inputenergy)");
	    System.exit(0);
        }else{
            fileName = args[0];
	    inP = Integer.valueOf(args[1]).intValue();
	    outP = Integer.valueOf(args[2]).intValue();
	    logYmin = Double.valueOf(args[3]).doubleValue();
	    logYmax = Double.valueOf(args[4]).doubleValue();
        }

	if(args.length == 6){
	    inLogE = Integer.valueOf(args[5]).intValue();
	    draw =
		new DrawPropagationMatrix(inP,outP,logYmin,logYmax,inLogE);
	}else{
	    draw = new DrawPropagationMatrix(inP,outP,logYmin,logYmax);
	}

	// Read the serialized object of the Neutrino Charged Interaction Matrix
	DataInputStream in = new DataInputStream(new FileInputStream(fileName));
	draw.readMatrix(in);
	in.close( );

	draw.drawOnXfig( );

    }

    /** or Glashow Resonance */
    // Switch includeGlashowResonance flag.
    public void whetherPropagationMatrixWithGlashowResonance(boolean flag){
	includeGlashowResonance = flag;
    }

}

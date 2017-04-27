package iceCube.uhe.analysis;

import iceCube.uhe.particles.*;
import iceCube.uhe.propagation.*;
import iceCube.uhe.analysis.*;
import geometry.*;

import java.io.*;
import java.util.*;

/** 
    This class provides the methods to fill I3Particles with 
    the pregenerated propagation matrix. PropagationMatrixFactory
    object in the iceCube.uhe.propagation package is a wrapper
    to read out the matrix.
    Every methods are static so that you can call them 
    directly from your object.

    Matrix(nuE -> I3Particle) + Matrix(nuMu -> I3Particle)
    + Matrix(nuTau -> I3Particle) is filled based upon
    the assumption that the primary cosmic neutrino flux
    at the earth surface is equal between the three flavors.

    Written by S. Yoshida 2007 March 21st
*/

public class I3ParticlePropMatrixFiller {

    Particle nuE;  // electron neutrino as a primary particle
    Particle nuMu; // muon neutrino as a primary particle
    Particle nuTau;// tau neutrino as a primary particle

    /** Constructor. Generate nu-e, nu-mu, and nu-tau objects
	as primary neutrinos at the earth surface.
    */
    public I3ParticlePropMatrixFiller(){
	nuE = new Particle(0,0);// electron neutrino as a primary particle
	nuMu = new Particle(1,0);// muon neutrino as a primary particle
	nuTau =new Particle(2,0);// tau neutrino as a primary particle
    }

    /** 
	Fill the propagation matrix with PropagatingMatrixFactory. 

	I3Particle iceParticle    : inIceParticle (IceCube MC event)
	PropagationMatrixFactory matrix : matrix from the surface to InIce
    */
    public boolean fillPropagationMatrixWeight(I3Particle iceParticle, 
					    PropagationMatrixFactory matrix)
	throws IOException{

	    // All the parameters should be MC truth 
	    iceParticle.switchToMCTruth();

	    J3UnitVector n_ice = iceParticle.getDirectionInIceCubeCoordinate();
	    double distance = 
		iceParticle.getDistanceFromEarthSurfaceToIceCube();
	    boolean isDownWard = true;
	    double cosZenith = -n_ice.getZ(); // Reversed vector
	    if(cosZenith<0.0) isDownWard = false; // this is an upward event

	    // This event is so small. Skip the matrix filling
	    if(iceParticle.getIceCubeData().getNDOMsLaunch()<
	       I3ParticleWeightFiller.minNDOMsToFill){
		return false;
	    }

	    readPropagationMatrix(n_ice,distance,isDownWard,matrix);
	    iceParticle.generateLogEnergyMatrix();

	    double logInIceEnergy = iceParticle.getLogEnergy();

	    double epsilon = 1.0e-4; // round-off margin for binning
	    double logSurfaceEnergy = Particle.getLogEnergyMinimum() + epsilon;
	    double logSurfaceMaxEnergy = Particle.getLogEnergyMinimum()
		+ Particle.getDeltaLogEnergy()*(double )(Particle.getDimensionOfLogEnergyMatrix());

	    int iLogE = 0;
	    while(logSurfaceEnergy < logSurfaceMaxEnergy){ // E < 10^12 GeV

		nuE.putLogEnergy(logSurfaceEnergy);
		double flux = matrix.getDF(nuE,iceParticle);

		nuMu.putLogEnergy(logSurfaceEnergy);
		flux += matrix.getDF(nuMu,iceParticle);

		nuTau.putLogEnergy(logSurfaceEnergy);
		flux += matrix.getDF(nuTau,iceParticle);

		iceParticle.putLogEnergyMatrix(iLogE,flux);

		iLogE++;
		logSurfaceEnergy += Particle.getDeltaLogEnergy();

		//System.err.println("logSurfaceEnergy= " + logSurfaceEnergy + " iLogE=" + iLogE);
	    }

	    return true;
    }

    /** 
	Reading the propagation matrix stored in the file corresponsing
	to the track geometry. 
	<pre>
	J3UnitVector n_ice    : Direction of I3Particle trajectory in the IceCube Coordinate
	double distance    : Propagation distance of I3Particle from the earth surface [cm]
	boolean isDownWard :   whether this event is downward-going or not
	</pre>
    */
    private static void readPropagationMatrix(J3UnitVector n_ice, double distance,
					      boolean isDownWard,
					      PropagationMatrixFactory matrix) 
	throws IOException {

	String matrixFileName = 
	    I3ParticleWeightFiller.getMatrixFileName(n_ice,distance,isDownWard);
	System.err.println(" --matrix filename " + matrixFileName);
	DataInputStream in = 
	    new DataInputStream(new FileInputStream(matrixFileName));
	matrix.readMatrix(in);
	in.close( );

    }

    /** Main method -- Reading out the stored I3Particles and fills 
	the propagation matrix */
    public static void main(String[] args) throws IOException{

	String inputFileName = null;
	String outputFileName = null;
	if(args.length!=2){
	    System.out.println("Usage: I3ParticlePropMatrixFiller input-file-name output-file-name");
	    System.exit(0);
        }else{
            inputFileName = args[0];
            outputFileName = args[1];
        }

	I3ParticlePropMatrixFiller filler = new I3ParticlePropMatrixFiller();

	// Open Input data stream
	InputStream in = ClassLoader.getSystemResourceAsStream(inputFileName);
	// Open Output data stream 
	FileOutputStream out = new FileOutputStream(outputFileName);

	// PropagatingMatrixFactory object
	PropagationMatrixFactory matrix = new PropagationMatrixFactory();

	I3Particle iceParticle = null; int numberOfData = 0; int numberOfFilledData = 0;
	while((iceParticle = I3ParticleInputStream.inputI3Particle(in))
	      !=null){

	    numberOfData++;

	    // All the parameters should be MC truth 
	    iceParticle.switchToMCTruth();

	    // Fills the propagation matrix weights
	    if(filler.fillPropagationMatrixWeight(iceParticle,
						  matrix)){
		numberOfFilledData++;
		System.err.println(numberOfData + " filled(" + numberOfFilledData + ") " +
				   iceParticle.particleName(iceParticle.getFlavor(),
							iceParticle.getDoublet()) +
			       " propMtx(" + iceParticle.getLogEnergy() + 
			       " from max)=" + iceParticle.getLogEnergyMatrix(
			     iceParticle.getDimensionOfLogEnergyMatrix()-1));
	    }
							
	    // Write out I3Particle object
	    I3ParticleOutputStream.outputI3Particle(iceParticle, out);
	}
	in.close();
	out.close();
    }

}

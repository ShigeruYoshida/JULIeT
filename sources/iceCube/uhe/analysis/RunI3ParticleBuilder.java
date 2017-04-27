package iceCube.uhe.analysis;

import iceCube.uhe.particles.*;
import iceCube.uhe.interactions.*;
import iceCube.uhe.analysis.*;
import geometry.*;

import java.io.*;
import java.util.*;

/** Run I3ParticleBulder to  generate a series of I3Particles.
    Generated particles are written out to file with name you
    specified.

    If you set isMCTruth = true in the main method,
    then I3ParticleBuilder fills the primary spectrum weight. 
*/
public class RunI3ParticleBuilder {

    public static void main(String[] args) throws IOException{

	boolean isMCTruth = true;
	boolean isFullData = true;
	double powerLaw = 1.0; // power law index of primary spectrum
	boolean isNeutrino = true;
	InteractionsMatrix nuCCmtx = null;   // CC cross section for weighting
	InteractionsMatrix nuNCmtx = null;  // NC cross section for weighting
	String ccMtxFileName = "iceCube/uhe/interactions/ice/MuNeutrinoChargeMtx";
	String ncMtxFileName = "iceCube/uhe/interactions/ice/MuNeutrinoNeutralMtx";

	String fileName = null;
	I3Particle iceTrack;

        if(args.length!=1){
            System.out.println("Usage: I3ParticleStreamDemo file-name");
	    System.exit(0);
        }else{
             fileName = args[0];
        }

	/** I3ParticleBuilder */
	I3ParticleBuilder builder =  new I3ParticleBuilder(isMCTruth);
	if(isNeutrino){
	    InputStream in = 
		ClassLoader.getSystemResourceAsStream(ccMtxFileName);
	    nuCCmtx = InteractionsMatrixInput.inputInteractionsMatrix(in);
	    in.close( );
	    in = ClassLoader.getSystemResourceAsStream(ncMtxFileName);
	    nuNCmtx = InteractionsMatrixInput.inputInteractionsMatrix(in);
	    in.close( );
	}

	// weight
	if(isMCTruth) builder.fillMCPrimarySpectrumWeight(powerLaw);
	if(isNeutrino) builder.fillNeutrinoWeight(nuCCmtx,nuNCmtx);

	/** Open Output data stream */
	FileOutputStream out = new FileOutputStream(fileName);
	DataInputStream input = new DataInputStream(System.in);

	/** Builds the I3Particle from the data via input.*/
	builder.process(input, out, isFullData);
	out.close( );
	System.err.println("Building The I3Particle  done.");

	/** Check the serialized I3Particle objects stored in disk */

	FileInputStream in = new FileInputStream(fileName);
	I3Particle iceParticle = null;
	while((iceParticle = I3ParticleInputStream.inputI3Particle(in))
	      !=null){

	    double npeATWD = iceParticle.getIceCubeData().getNpeATWD();
	    double npeFADC = iceParticle.getIceCubeData().getNpeFADC();
	    double npeBest = iceParticle.getIceCubeData().getBestNpe();
	    if(npeBest>1.0e4){
		System.out.print("- EventNumber is " + 
				 iceParticle.getIceCubeData().getEventNumber());
		System.out.println(" The I3Particle Name is " + 
				   iceParticle.particleName(iceParticle.getFlavor(), 
							    iceParticle.getDoublet()));

		System.out.print("- The I3Particle Energy (" + 
				 iceParticle.getEnergy( ) + ") [GeV] ");
		System.out.println("- The I3Particle Reco Energy (" + 
				   iceParticle.getRecoEnergy( ) + ") [GeV] ");
		J3UnitVector n = iceParticle.getDirectionInIceCubeCoordinate();
		System.out.println("- The I3Particle direction in ice3(" + 
				   n.getX() + ", " + n.getY() + ", " + n.getZ() +
				   ")");
		System.out.print("- Npe by ATWD = " + 
				   iceParticle.getIceCubeData().getNpeATWD());
		System.out.print(" Npe by FADC = " + 
				   iceParticle.getIceCubeData().getNpeFADC());
		System.out.print(" Best Npe = " + 
				   iceParticle.getIceCubeData().getBestNpe());
		double fitQuality = iceParticle.getFirstGuessQuality();
		System.out.println("  FG velocity(" + fitQuality + ")");
	    }
	}

	in.close();

	System.err.println(" everything done --");
    }
}

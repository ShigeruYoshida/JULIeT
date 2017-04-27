package iceCube.uhe.analysis;

import iceCube.uhe.particles.*;
import iceCube.uhe.points.*;
import iceCube.uhe.muonModel.*;
import iceCube.uhe.neutrinoModel.*;
import iceCube.uhe.geometry.*;
import iceCube.uhe.analysis.*;
import geometry.*;

import java.io.*;
import java.util.*;

/** 
    This class provides the methods to fill I3Particles with MC primary spectrum
    weights.

    Written by S. Yoshida 2007 February 11th
*/

public class I3ParticleMCPrimaryWeightFiller {


    /** Main method -- Reading out the stored I3Particles and fills the weight */
    public static void main(String[] args) throws IOException{

	String inputFileName = null;
	String outputFileName = null;
	double powerLaw = 1.0;
	double log_EnergyMinimum = 5.0;   // log[GeV]
	double log_EnergyMaximum = 11.0;   // log[GeV]
	if(args.length!=3){
	    System.out.println("Usage: I3ParticleMCPrimaryWeightFiller input-file-name output-file-name powerLaw");
	    System.exit(0);
        }else{
            inputFileName = args[0];
            outputFileName = args[1];
	    powerLaw = Double.valueOf(args[2]).doubleValue();
        }

	// Open Input data stream
	InputStream in = ClassLoader.getSystemResourceAsStream(inputFileName);
	// Open Output data stream 
	FileOutputStream out = new FileOutputStream(outputFileName);


	I3Particle iceParticle = null; 
	while((iceParticle = I3ParticleInputStream.inputI3Particle(in))
	      !=null){

	    // All the parameters should be MC truth 
	    iceParticle.switchToMCTruth();
	    double energy = iceParticle.getEnergy();
	    double flux = I3ParticleAnalysisFactory.getDNDLogE(powerLaw,energy,
							       log_EnergyMinimum,
							       log_EnergyMaximum);

	    iceParticle.setMCPrimarySpectrumWeight(flux);

	    // Write out I3Particle object
	    I3ParticleOutputStream.outputI3Particle(iceParticle, out);
	}
	in.close();
	out.close();
    }
}

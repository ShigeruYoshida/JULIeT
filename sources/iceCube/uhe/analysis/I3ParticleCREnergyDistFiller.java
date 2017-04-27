
package iceCube.uhe.analysis;

import iceCube.uhe.particles.*;
import iceCube.uhe.muonModel.*;
import iceCube.uhe.analysis.*;

import java.io.*;
import java.util.*;


public class I3ParticleCREnergyDistFiller {


    /** Main method -- Reading out the stored I3Particles and fills 
	the CR Energy spectrum */
    public static void main(String[] args) throws IOException{

	String inputFileName = null;
	String outputFileName = null;
	double muEthreshold = 1505.0; // IC22 data default
	double alpha = 1.97; // IC22 data default
	if(args.length!=2){
	    System.out.println("Usage: I3ParticleCREnergyDistFiller input-file-name output-file-name");
	    System.exit(0);
        }else{
            inputFileName = args[0];
            outputFileName = args[1];
        }

	// Open Input data stream
	InputStream in = ClassLoader.getSystemResourceAsStream(inputFileName);
	// Open Output data stream 
	FileOutputStream out = new FileOutputStream(outputFileName);

	PropagatingAtmMuonFlux muonFlux = new PropagatingAtmMuonFlux(alpha,muEthreshold, true);//GZK cutoff

	I3Particle iceParticle = null;
	while((iceParticle = I3ParticleInputStream.inputI3Particle(in))
	      !=null){

	    // All the parameters should be MC truth 
	    iceParticle.switchToMCTruth();

	    // Fills the propagation matrix weights
	    I3ParticleWeightFiller.fillCRFluxWeight(iceParticle,muonFlux);
							
	    // Write out I3Particle object
	    I3ParticleOutputStream.outputI3Particle(iceParticle, out);
	}
	in.close();
	out.close();
    }

}

package iceCube.uhe.analysis;

import iceCube.uhe.particles.*;
import iceCube.uhe.points.*;
import iceCube.uhe.muonModel.*;
import iceCube.uhe.neutrinoModel.*;
import iceCube.uhe.analysis.*;
import geometry.*;

import java.io.*;
import java.util.*;


public class I3ParticleCorsikaWeightFiller {

    /** Main method -- Reading out the stored I3Particles and fills the weight */
    public static void main(String[] args) throws IOException{

	String inputFileName = null;
	String outputFileName = null;
	String atmMuon_name = "corsika"; 


	if(args.length<2){
	    System.out.println("Usage: I3ParticleCorsikaWeightFiller input-file-name output-file-name");
	    System.exit(0);
        }else { 
            inputFileName = args[0];
            outputFileName = args[1];
        }


	System.out.println("Atm muon model Name " + atmMuon_name);

	// Open Input data stream
	InputStream in = ClassLoader.getSystemResourceAsStream(inputFileName);
	// Open Output data stream 
	FileOutputStream out = new FileOutputStream(outputFileName);


	// PropagatingNeutrinoFlux object
	CosmicRayFlux crFlux = new CosmicRayFlux();
	crFlux.setCutOffFeature(true); // with the GZK cutoff

	I3Particle iceParticle = null; 
	while((iceParticle = I3ParticleInputStream.inputI3Particle(in))
	      !=null){

	    // All the parameters should be MC truth 
	    iceParticle.switchToMCTruth();

	    // Remove the previousely stored weight by the same model name
	    // if exists
	    iceParticle.removeAtmosphericMuonFlux(atmMuon_name);

	    // Fills the flux weights

	    I3ParticleWeightFiller.fillPropagatingAtmMuonFluxWeight(iceParticle,
								     crFlux,
								     atmMuon_name);
	    // Write out I3Particle object
	    I3ParticleOutputStream.outputI3Particle(iceParticle, out);
	}
	in.close();
	out.close();
    }
}

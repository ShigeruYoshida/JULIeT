package iceCube.uhe.analysis;

import iceCube.uhe.particles.*;
import iceCube.uhe.points.*;
import iceCube.uhe.muonModel.*;
import iceCube.uhe.neutrinoModel.*;
import iceCube.uhe.analysis.*;
import geometry.*;

import java.io.*;
import java.util.*;


public class I3ParticleGZKWeightFiller {

    /** Main method -- Reading out the stored I3Particles and fills the weight */
    public static void main(String[] args) throws IOException{

	String inputFileName = null;
	String outputFileName = null;
	String gzk_YT_name = "GZK_YT_4_4"; int gzk_YT_number = 4;
	String gzk_sigl_name = "GZK_sigl"; int gzk_sigl_number = 11;
	String td_name = "TD_SUSY"; int td_number = 8;


	PropagatingNeutrinoFlux propNeutFlux = null;

	if(args.length<2){
	    System.out.println("Usage: I3ParticleGZKWeightFiller input-file-name output-file-name");
	    System.exit(0);
        }else { 
            inputFileName = args[0];
            outputFileName = args[1];
        }


	System.out.println("Neutrino model Name " + gzk_YT_name + 
			   " model number(" + gzk_YT_number + ")"); 
	System.out.println("Neutrino model Name " + gzk_sigl_name + 
			   " model number(" + gzk_sigl_number + ")"); 
	System.out.println("Neutrino model Name " + td_name + 
			   " model number(" + td_number + ")"); 

	// Open Input data stream
	InputStream in = ClassLoader.getSystemResourceAsStream(inputFileName);
	// Open Output data stream 
	FileOutputStream out = new FileOutputStream(outputFileName);


	// PropagatingNeutrinoFlux object
	propNeutFlux = new PropagatingNeutrinoFlux(gzk_YT_number);

	I3Particle iceParticle = null; 
	while((iceParticle = I3ParticleInputStream.inputI3Particle(in))
	      !=null){

	    // All the parameters should be MC truth 
	    iceParticle.switchToMCTruth();

	    // Remove the previousely stored weight by the same model name
	    // if exists
	    iceParticle.removeGZKNeutrinoFlux(gzk_YT_name);
	    iceParticle.removeGZKNeutrinoFlux(gzk_sigl_name);
	    iceParticle.removeGZKNeutrinoFlux(td_name);

	    // Fills the flux weights

	                   // This method also reads out the prop matix
	    I3ParticleWeightFiller.fillPropagatingNeutrinoFluxWeight(iceParticle,
								     propNeutFlux,
								     gzk_YT_name);

	    I3ParticleWeightFiller.fillPropagatingNeutrinoFluxWeight(iceParticle,
								     propNeutFlux,
								     gzk_sigl_number,
								     gzk_sigl_name);

	    I3ParticleWeightFiller.fillPropagatingNeutrinoFluxWeight(iceParticle,
								     propNeutFlux,
								     td_number,
								     td_name);
	    // Write out I3Particle object
	    I3ParticleOutputStream.outputI3Particle(iceParticle, out);
	}
	in.close();
	out.close();
    }
}

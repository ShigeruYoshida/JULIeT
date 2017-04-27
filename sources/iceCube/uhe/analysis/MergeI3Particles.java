package iceCube.uhe.analysis;

import iceCube.uhe.analysis.*;
import iceCube.uhe.particles.*;

import java.io.*;
import java.util.*;

public class MergeI3Particles {

    public static void main(String[] args) throws IOException{

	String mcDataFileName = null;
	String mcDataFileName2 = null;
	String mcDataOutputFileName = null;

        if(args.length!=3){
            System.out.println("Usage: MergeI3Particles inputFileName inputFileName2 outputFileName");
	    System.exit(0);
        }else{
	    mcDataFileName = args[0];
	    mcDataFileName2 = args[1];
	    mcDataOutputFileName = args[2];
        }

	/** Open Output data stream */
	FileOutputStream out = new FileOutputStream(mcDataOutputFileName);
	FileInputStream in = new FileInputStream(mcDataFileName); 

	I3Particle iceParticle = null; 
	while((iceParticle = I3ParticleInputStream.inputI3Particle(in))
	      !=null){
	    // Write out I3Particle object
	    I3ParticleOutputStream.outputI3Particle(iceParticle, out);

	}
	in.close();


	in = new FileInputStream(mcDataFileName2);
	iceParticle = null; 
	while((iceParticle = I3ParticleInputStream.inputI3Particle(in))
	      !=null){
	    // Write out I3Particle object
	    I3ParticleOutputStream.outputI3Particle(iceParticle, out);

	}
	in.close();
	out.close();
    }
}

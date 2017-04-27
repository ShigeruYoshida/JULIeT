package iceCube.uhe.analysis;

import iceCube.uhe.particles.*;
import iceCube.uhe.analysis.*;
import geometry.*;

import java.io.*;
import java.util.*;

/** I3Particle.class Demo program. */
public class I3ParticleBuilderDemo {

    public static void main(String[] args) throws IOException{

	String fileName = null;
	I3Particle iceTrack;

        if(args.length!=1){
            System.out.println("Usage: I3ParticleStreamDemo file-name");
	    System.exit(0);
        }else{
             fileName = args[0];
        }

	/** I3ParticleBuilder */
	boolean isMCTruth = true;
	I3ParticleBuilder builder =  new I3ParticleBuilder(isMCTruth);

	/** Open Output data stream */
	FileOutputStream out = new FileOutputStream(fileName);
	boolean isFullData = false;
	DataInputStream input = new DataInputStream(System.in);
	builder.process(input, out, isFullData);
	out.close( );
	System.out.println("Building The I3Particle  done.");

	/** Input the serialized I3Particle class */

	FileInputStream in = new FileInputStream(fileName);
	I3Particle iceParticle = null;
	while((iceParticle = I3ParticleInputStream.inputI3Particle(in))
	      !=null){


	    System.out.println("- The I3Particle Name is " + 
			       iceParticle.particleName(iceParticle.getFlavor(), 
							iceParticle.getDoublet()));

	    System.out.println("- The I3Particle Mass (" + 
			       iceParticle.getMass( )+ ") [GeV]");
	    System.out.println("- The I3Particle Energy (" + 
			       iceParticle.getEnergy( ) + ") [GeV] " + 
			       iceParticle.getLogEnergy( ));
	    System.out.println("- The I3Particle Reco Energy (" + 
			       iceParticle.getRecoEnergy( ) + ") [GeV] " + 
			       iceParticle.getLogRecoEnergy( ));
	    System.out.println("- The I3Particle lifetime (" + 
			       iceParticle.getLifeTime( ) + ") [sec]");

	    J3UnitVector n = iceParticle.getDirectionInIceCubeCoordinate();
	    System.out.println("- The I3Particle direction in ice3(" + 
			       n.getX() + ", " + n.getY() + ", " + n.getZ());

	    n = iceParticle.getDirectionInEarthCenterCoordinate();
	    System.out.println("- The I3Particle direction in EarthCenter(" + 
			       n.getX() + ", " + n.getY() + ", " + n.getZ());

	    J3Vector r = iceParticle.getR0InIceCubeCoordinate();
	    System.out.println("- The I3Particle location in ice3(" + 
			       r.getX() + ", " + r.getY() + ", " + r.getZ());
	    r = iceParticle.getR0InEarthCenterCoordinate();
	    System.out.println("- The I3Particle location in EarthCenter(" + 
			       r.getX() + ", " + r.getY() + ", " + r.getZ());

	    double fitQuality = iceParticle.getFirstGuessQuality();
	    System.out.println("- The I3Particle FG quality(" + 
			       fitQuality + ")");


	    System.out.println("- Npe by ATWD = " + 
			       iceParticle.getIceCubeData().getNpeATWD());
	    System.out.println("- Log Npe by ATWD = " + 
			       iceParticle.getIceCubeData().getLogNpeATWD());
	    System.out.println("- NDOMs with  ATWD = " + 
			       iceParticle.getIceCubeData().getNDOMsATWD());
	    System.out.println("- NDOMs with FADC = " + 
			       iceParticle.getIceCubeData().getNDOMsFADC());
	}

	in.close();

	System.out.println("Demo program ended --");
    }
}

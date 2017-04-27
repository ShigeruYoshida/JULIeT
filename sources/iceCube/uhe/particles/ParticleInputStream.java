package iceCube.uhe.particles;

import java.io.*;

/** The Object Particle is serialized and sent out to
    Outputstream */

public class ParticleInputStream {

    public static Particle inputParticle(InputStream in) 
    throws IOException {

        Particle particle = null;

	try{

	    ObjectInputStream objectIn = new ObjectInputStream(in);
	    particle = (Particle )objectIn.readObject();

	}catch(ClassNotFoundException e){
	    System.err.println("Caught ClassNotFoundException: " + 
			       e.getMessage( ));
	    System.exit(0);
	}catch (EOFException e){
	    System.err.println("Caught EOFException: " + e.getMessage());
	    particle = null;
	    return particle;
	}

	return particle;

    }

}


package iceCube.uhe.particles;

import java.io.EOFException;
import java.io.IOException;
import java.io.InputStream;
import java.io.ObjectInputStream;

/** The Object Particle is serialized and sent out to
    Outputstream */

public class I3ParticleInputStream {

    public static I3Particle inputI3Particle(InputStream in)
    throws IOException {

        I3Particle particle = null;
	try{

	    ObjectInputStream objectIn = new ObjectInputStream(in);
	    particle = (I3Particle )objectIn.readObject();

	}catch(ClassNotFoundException e){
	    System.err.println("Caught ClassNotFoundException: " + 
			       e.getMessage( ));
	    System.exit(0);
	}catch (EOFException e){
	    //System.err.println("Caught EOFException: " + e.getMessage());
	    particle = null;
	    return particle;
	}

	return particle;

    }

}


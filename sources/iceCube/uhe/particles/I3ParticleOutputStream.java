package iceCube.uhe.particles;

import java.io.*;

/** The Object I3Particle is serialized and sent out to
    Outputstream */

public class I3ParticleOutputStream {

    public static void outputI3Particle(I3Particle particles, OutputStream out) 
    throws IOException {

	ObjectOutputStream objectOut = new ObjectOutputStream(out);

	objectOut.writeObject(particles);
	objectOut.flush();
    }

}


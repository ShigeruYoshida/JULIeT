package iceCube.uhe.particles;

import java.io.*;

/** The Object Particle is serialized and sent out to
    Outputstream */

public class ParticleOutputStream {

    public static void outputParticle(Particle particles, OutputStream out) 
    throws IOException {

	ObjectOutputStream objectOut = new ObjectOutputStream(out);

	objectOut.writeObject(particles);
	objectOut.flush();
    }

}


package iceCube.uhe.interactions;

import java.io.*;

/** The Object InteractionsMatrix is serialized and sent out to
    Outputstream */

public class InteractionsMatrixOutput {

    public static void outputInteractionsMatrix(InteractionsMatrix intMtx, 
						OutputStream out) 
    throws IOException {

	ObjectOutputStream objectOut = new ObjectOutputStream(out);

	objectOut.writeObject(intMtx);
	objectOut.flush();
    }

}


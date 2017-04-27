package iceCube.uhe.interactions;

import java.io.*;

/** The Object InteractionsMatrix is read out from
    Inputstream */

public class InteractionsMatrixInput {

    public static InteractionsMatrix inputInteractionsMatrix(InputStream in) 
    throws IOException {

        InteractionsMatrix intMtx = null;
	ObjectInputStream objectIn = new ObjectInputStream(in);

	try{

	    intMtx = (InteractionsMatrix )objectIn.readObject();

	}catch(ClassNotFoundException e){
	    System.err.println("Caught ClassNotFoundException: " + 
			       e.getMessage( ));
	    System.exit(0);
	}

	return intMtx;

    }

}


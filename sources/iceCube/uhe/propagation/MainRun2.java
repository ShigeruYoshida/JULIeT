package iceCube.uhe.propagation;

import iceCube.uhe.particles.*;
import iceCube.uhe.interactions.*;
import iceCube.uhe.points.*;
import iceCube.uhe.propagation.*;
import java.io.*;


public class MainRun2 {
    public static void main(String[] args) throws IOException {

        String fileName = null;
	int intSwitch = 511;
	int decaySwitch = 255;
	int upDownFlag = 0;  // downgoing
	int mediumNumber = 0; //ice
	double nadirAngle = 0.0;
	double trajectoryLength = 1.0;
	double detectorDepth = 1.4e5;  // Detector Depth.. 1400m = 1.4e5 cm below sea level
	RunPropagationMatrix run = null;

        if(args.length!=2 ){
            System.out.println("Usage: MainRun2 log10(trajectoryLength[cm]) file-name");
	    System.exit(0);
        }else if(args.length==2){
	    double logLength = Double.valueOf(args[0]).doubleValue();
	    if(logLength>0.0) trajectoryLength = Math.pow(10.0,logLength);
            fileName = args[1];
	    System.err.println("Log(Distance)=" + logLength + " Distance=" + trajectoryLength + 
			       "[cm] filename " + fileName);
			       
	}

	System.err.println("Down-going propagation");
	mediumNumber = 0; //ice

	run = new RunPropagationMatrix(nadirAngle,intSwitch,decaySwitch,mediumNumber);
	    
	run.traceParticles(trajectoryLength);

	// Write out produced matrix
	DataOutputStream out = new DataOutputStream(new FileOutputStream(fileName));
	run.saveMatrix(out);

    }
}

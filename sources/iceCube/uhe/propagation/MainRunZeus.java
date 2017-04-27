package iceCube.uhe.propagation;

import iceCube.uhe.particles.*;
import iceCube.uhe.interactions.*;
import iceCube.uhe.points.*;
import iceCube.uhe.propagation.*;
import java.io.*;


public class MainRunZeus {
    public static void main(String[] args) throws IOException {

        String fileName = null;
	int intSwitch = 255;
	int decaySwitch = 255;
	int upDownFlag = 1;
	int mediumNumber = 1;
	double nadirAngle = 0.0;
	double trajectoryLength = 1.0;
	double detectorDepth = 1.4e5;  // Detector Depth.. 1400m = 1.4e5 cm below sea level
	RunPropagationMatrixZeus run = null;

        if(args.length!=5 && args.length!=6){
            System.out.println(
"Usage: MainRunPropagation NadirAngle[deg] log10(trajectoryLength[cm]) intSwitch decaySwitch upDownFlag file-name");
	    System.exit(0);
        }else if(args.length==5){
	    nadirAngle = Double.valueOf(args[0]).doubleValue();
            intSwitch = Integer.valueOf(args[1]).intValue();
            decaySwitch = Integer.valueOf(args[2]).intValue();
            upDownFlag = Integer.valueOf(args[3]).intValue();
            fileName = args[4];
        }else if(args.length==6){
	    nadirAngle = Double.valueOf(args[0]).doubleValue();
	    double logLength = Double.valueOf(args[1]).doubleValue();
            if(logLength>0.0) trajectoryLength = Math.pow(10.0,logLength);
             intSwitch = Integer.valueOf(args[2]).intValue();
            decaySwitch = Integer.valueOf(args[3]).intValue();
            upDownFlag = Integer.valueOf(args[4]).intValue();
            fileName = args[5];
            System.err.println("Log(Distance)=" + logLength + " Distance=" + trajectoryLength +
                               "[cm] filename " + fileName);
	}

	if(upDownFlag ==1){
	    System.err.println("Up-going propagation");
	    mediumNumber = 1; //Rock
	    run =
		new RunPropagationMatrixZeus(nadirAngle,intSwitch,decaySwitch,mediumNumber);

	    if(args.length==5){
		run.traceParticles( );
	    } else{
		run.traceParticles(trajectoryLength);
	    }

	}else if(upDownFlag == 0){
	    System.err.println("Down-going propagation");
	    mediumNumber = 0; //ice
	    double zenithAngle =  nadirAngle;  // Zenith Angle  0deg for vertical
	    double cos_zenith = Math.cos(zenithAngle*Math.PI/180.0);
	    double sq_term = Math.sqrt((ParticlePoint.REarth-detectorDepth)*(ParticlePoint.REarth-detectorDepth)
				       *cos_zenith*cos_zenith + 
				       2.0*ParticlePoint.REarth*detectorDepth-detectorDepth*detectorDepth);
	    double cos_nadir = sq_term/ParticlePoint.REarth;
	    nadirAngle = Math.acos(cos_nadir)*180.0/Math.PI;

	    run =
		new RunPropagationMatrixZeus(nadirAngle,intSwitch,decaySwitch,mediumNumber);
	    
	    if(args.length==5){
		trajectoryLength = sq_term - (ParticlePoint.REarth-detectorDepth)*cos_zenith;
		System.err.println("Zenith " + zenithAngle + " Nadir " + nadirAngle + 
				   " Trajectory " + trajectoryLength/100.0 + " [m]");
		run.traceParticles(trajectoryLength);
	    } else{
		run.traceParticles(trajectoryLength);
	    }

	}
	// Write out produced matrix
	DataOutputStream out = new DataOutputStream(new FileOutputStream(fileName));
	run.saveMatrix(out);

    }
}

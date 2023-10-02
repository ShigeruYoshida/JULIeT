package iceCube.uhe.propagation;

import iceCube.uhe.particles.*;
import iceCube.uhe.interactions.*;
import iceCube.uhe.points.*;
import iceCube.uhe.propagation.*;
import java.io.*;


public class MainRunGlashowNLO {
    public static void main(String[] args) throws IOException {

        String fileName = null;
        int intSwitch = 255;
        int decaySwitch = 255;
        int upDownFlag = 1;
        int mediumNumber = 1;
        double nadirAngle = 0.0;
        double trajectoryLength = 1.0;
        double logLength = 0.0;
        double detectorDepth = 1.4e5;  // Detector Depth.. 1400m = 1.4e5 cm below sea level
        double neutrinoFactor = 1.0;
        RunPropagationMatrixGlashowNLO run = null;
        boolean receivedTrajLength = true;
        int idxOffset = 1;
        boolean receivedCustomNuFiles = false;
        String customNuCCMtxFile = null;
        String customNuNCMtxFile = null;

        String cliUsage = (
            "Usage: MainRunPropagation NadirAngle[deg] log10(trajectoryLength[cm])[optional] " +
            "intSwitch decaySwitch upDownFlag file-name neutrinoFactor " +
            "customNuCCMtxFile[optional] customNuNCMtxFile[optional]");

        if(args.length < 6 || args.length > 9) {
            System.out.println(cliUsage);
            System.exit(0);
        } else {
            if (args.length % 2 == 0) {
                receivedTrajLength = false;
                idxOffset = 0;
            }
            else {
                logLength = Double.parseDouble(args[1]);
            }
            System.out.println("logLength: " + logLength);

            nadirAngle = Double.valueOf(args[0]).doubleValue();

            intSwitch = Integer.valueOf(args[1+idxOffset]).intValue();
            decaySwitch = Integer.valueOf(args[2+idxOffset]).intValue();
            upDownFlag = Integer.valueOf(args[3+idxOffset]).intValue();
            fileName = args[4+idxOffset];
            neutrinoFactor = Double.valueOf(args[5+idxOffset]).doubleValue();

            if (receivedTrajLength){
                if(logLength > 0.0) trajectoryLength = Math.pow(10., logLength);
                System.err.println(
                    "Log(Distance)= " + logLength +
                    " Distance= " + trajectoryLength +
                    "[cm] filename " + fileName);
            }

            if (args.length > 6 + idxOffset){
                if (args.length == 6 + idxOffset + 2){
                    customNuCCMtxFile = args[6+idxOffset];
                    customNuNCMtxFile = args[6+idxOffset+1];
                    receivedCustomNuFiles = true;
                }
                else {
                    System.out.println(cliUsage);
                    System.out.println(
                        "When providing custom neutrino cross section matrices, " +
                        "they have to be provided for both CC and NC!");
                    System.exit(0);
                }
            }
        }

    if (upDownFlag == 1) {
        System.err.println("Up-going propagation");
        mediumNumber = 1; //Rock

        if (receivedCustomNuFiles){
            run = new RunPropagationMatrixGlashowNLO(
                nadirAngle,
                intSwitch,
                decaySwitch,
                mediumNumber,
                neutrinoFactor,
                customNuCCMtxFile,
                customNuNCMtxFile);
        } else {
            run = new RunPropagationMatrixGlashowNLO(
                nadirAngle,
                intSwitch,
                decaySwitch,
                mediumNumber,
                neutrinoFactor);
        }

        if (receivedTrajLength) {
            run.traceParticles(trajectoryLength);
        } else {
            run.traceParticles();
        }

    } else if(upDownFlag == 0) {
        System.err.println("Down-going propagation");
        mediumNumber = 0; //ice
        double zenithAngle =  nadirAngle;  // Zenith Angle  0deg for vertical
        double cos_zenith = Math.cos(zenithAngle*Math.PI/180.0);
        double sq_term = Math.sqrt((ParticlePoint.REarth-detectorDepth)*(ParticlePoint.REarth-detectorDepth)
                       *cos_zenith*cos_zenith + 
                       2.0*ParticlePoint.REarth*detectorDepth-detectorDepth*detectorDepth);
        double cos_nadir = sq_term/ParticlePoint.REarth;
        nadirAngle = Math.acos(cos_nadir)*180.0/Math.PI;

        if (receivedCustomNuFiles){
            run = new RunPropagationMatrixGlashowNLO(
                nadirAngle,
                intSwitch,
                decaySwitch,
                mediumNumber,
                neutrinoFactor,
                customNuCCMtxFile,
                customNuNCMtxFile);
        } else {
            run = new RunPropagationMatrixGlashowNLO(
                nadirAngle,
                intSwitch,
                decaySwitch,
                mediumNumber,
                neutrinoFactor);
        }
        
        if (receivedTrajLength) {
            run.traceParticles(trajectoryLength);
        } else {
            trajectoryLength = sq_term - (ParticlePoint.REarth-detectorDepth)*cos_zenith;
            System.err.println("Zenith " + zenithAngle + " Nadir " + nadirAngle + 
                       " Trajectory " + trajectoryLength/100.0 + " [m]");
            run.traceParticles(trajectoryLength);
        }

    }
    // Write out produced matrix
    DataOutputStream out = new DataOutputStream(new FileOutputStream(fileName));
    run.saveMatrix(out);
    }
}

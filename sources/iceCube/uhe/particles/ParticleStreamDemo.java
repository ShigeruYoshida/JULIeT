package iceCube.uhe.particles;

import iceCube.uhe.particles.*;
import java.io.*;

/** Particle.class Demo program. */
public class ParticleStreamDemo {

    public static void main(String[] args) throws IOException{

	int doublet =1;
	String fileName = null;
	Particle[] uheParticles;

        if(args.length!=1){
            System.out.println("Usage: ParticleStreamDemo file-name");
	    System.exit(0);
        }else{
             fileName = args[0];
        }

	// Generate the particle class.
	uheParticles = new Particle[Particle.NumberOfFlavor];

	for(int flavor=0;flavor<Particle.NumberOfFlavor;flavor++){
	    uheParticles[flavor] = 
		new Particle(flavor, doublet, Math.pow(10.0,8.0+(double )flavor));

	    // Generate the logE matrix.
	    uheParticles[flavor].generateLogEnergyMatrix( );

	    // Put the values into the matrix.
	    for(int i=0;i<Particle.getDimensionOfLogEnergyMatrix();i++){
		double logEnergy = Particle.getLogEnergyMinimum()+
		    Particle.getDeltaLogEnergy()*(double )i;
		double energy = Math.pow(10.0,logEnergy);
		double sigma=6.9542e-34*Math.pow(energy/1.0e6,0.402);
		uheParticles[flavor].putLogEnergyMatrix(i,sigma);
	    }
	}

	for(int flavor=0;flavor<Particle.NumberOfFlavor;flavor++){
	    System.out.println("The Particle Name is " + 
            uheParticles[flavor].particleName(uheParticles[flavor].getFlavor(), 
			uheParticles[flavor].getDoublet()));

	    System.out.println("The Particle Mass (" + 
			       uheParticles[flavor].getMass( )+ ") [GeV]");
	    System.out.println("The Particle Energy (" + 
			       uheParticles[flavor].getEnergy( ) + ") [GeV] " + 
			       uheParticles[flavor].getLogEnergy( ));
	    System.out.println("The Particle lifetime (" + 
			       uheParticles[flavor].getLifeTime( ) + ") [sec]");

	    System.out.println("Matrix(" + 100 +")=" +
			       uheParticles[flavor].getLogEnergyMatrix(100));
	}

	/** Output the serialized Particle class */

	FileOutputStream out = new FileOutputStream(fileName);
	for(int flavor=0;flavor<Particle.NumberOfFlavor;flavor++){
	    ParticleOutputStream.outputParticle(uheParticles[flavor], out);
	}
	out.close( );

	System.out.println("The Particle output done.");




	/** Input the serialized Particle class */

	FileInputStream in = new FileInputStream(fileName);


	Particle[] vheParticles ;
	vheParticles = new Particle[Particle.NumberOfFlavor];
	for(int flavor=0;flavor<Particle.NumberOfFlavor;flavor++){
	    vheParticles[flavor] = ParticleInputStream.inputParticle(in);
	}
	in.close( );


	for(int flavor=0;flavor<Particle.NumberOfFlavor;flavor++){
	    System.out.println("The Particle Name is " + 
            vheParticles[flavor].particleName(vheParticles[flavor].getFlavor(), 
			vheParticles[flavor].getDoublet()));

	    System.out.println("The Particle Mass (" + 
			       vheParticles[flavor].getMass( )+ ") [GeV]");
	    System.out.println("The Particle Energy (" + 
			       vheParticles[flavor].getEnergy( ) + ") [GeV] " + 
			       vheParticles[flavor].getLogEnergy( ));
	    System.out.println("The Particle lifetime (" + 
			       vheParticles[flavor].getLifeTime( ) + ") [sec]");

	    System.out.println("Matrix(" + 100 +")=" +
			       vheParticles[flavor].getLogEnergyMatrix(100));
	}


    }
}

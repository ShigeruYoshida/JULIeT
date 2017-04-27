package iceCube.uhe.particles;

import iceCube.uhe.particles.*;

/** Particle.class Demo program. */

public class ParticleDemo {

    public static void main(String[] args){

	int flavor =0;
	int doublet =0;

        if(args.length!=2){
            System.out.println("Usage: ParticleDemo flavor doublet");
        }else{
            flavor = Integer.valueOf(args[0]).intValue();
            doublet = Integer.valueOf(args[1]).intValue();
        }

	// Generate the particle class.
	Particle uheParticle = new Particle(flavor, doublet, Math.pow(10.0,9.02));


	System.out.println("The Particle flavor (" + uheParticle.getFlavor( )
			   + ")");
	System.out.println("The Particle doublet (" + uheParticle.getDoublet( )
			   + ")");
	System.out.println("The Particle Mass (" + uheParticle.getMass( )
			   + ") [GeV]");
	System.out.println("The Particle Energy (" + uheParticle.getEnergy( )
			   + ") [GeV] " + uheParticle.getLogEnergy( ));
	System.out.println("The Particle lifetime (" + uheParticle.getLifeTime( )
			   + ") [sec]");

	uheParticle.putEnergy(1.0e6);
	uheParticle.putLogEnergy(6.0);

	System.out.println("The Particle Energy (" + uheParticle.getEnergy( )
			   + ") [GeV] " + uheParticle.getLogEnergy( ));

	System.out.println("The Particle Name is " + 
	uheParticle.particleName(uheParticle.getFlavor(), 
				 uheParticle.getDoublet()));

	// Generate the logE matrix.
	uheParticle.generateLogEnergyMatrix( );

	// Put the values into the matrix.
	for(int i=0;i<Particle.getDimensionOfLogEnergyMatrix();i++){
	    double logEnergy = Particle.getLogEnergyMinimum()+
		Particle.getDeltaLogEnergy()*(double )i;
	    double energy = Math.pow(10.0,logEnergy);
	    double sigma=6.9542e-34*Math.pow(energy/1.0e6,0.402);
	    //uheParticle.putLogEnergyMatrix(logEnergy,sigma);
	    uheParticle.putLogEnergyMatrix(i,sigma);
	}
	for(int i=0;i<Particle.getDimensionOfLogEnergyMatrix();i++){
	    double logEnergy = Particle.getLogEnergyMinimum()+
		Particle.getDeltaLogEnergy()*(double )i;
	    double energy = Math.pow(10.0,logEnergy);
	    System.out.println("logE(" + logEnergy + ") Matrix(" + i +")=" +
			       uheParticle.getLogEnergyMatrix(i));
	}

    }
}

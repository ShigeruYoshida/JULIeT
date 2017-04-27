package iceCube.uhe.event;

import iceCube.uhe.particles.*;
import iceCube.uhe.interactions.*;
import iceCube.uhe.event.*;
import iceCube.uhe.decay.*;
import iceCube.uhe.points.*;
import numRecipes.*;

import java.io.*;
import java.util.*;

/** Uses a SwingWorker to perform time-consuming tasks in JulietEventGenerator
such as generating the Interaction Matrix and running 
the propagating particle. These tasks are performed in a thred
for the SWING-based GUI application. It also generates the message
to indicate the progress of the invoked task for the ProgressBar
in the SWING. This is a sort of "interface" class between SWING and 
JULIeT. */

public class JulietTask {
    private int lengthOfTask;
    private int current = 0;
    private boolean done = false;
    private boolean propagationDone = false;
    private boolean canceled = false;
    private String statMessage;
    private SwingRun swingRun;

    /** JULIeTGenerator. Generate particles and let them run*/
    public JulietEventGenerator generator= null;


    public JulietTask(SwingRun swingRun) {
        //Compute length of task...
        //In a real program, this would figure out
        //the number of bytes to read or whatever.
        lengthOfTask = 1000;
	this.swingRun =  swingRun;
    }

    /**
     * start the task to generate a particle.
     */
    public void generate() {
        final SwingWorker worker = new SwingWorker() {
            public Object construct() {
                current = 0;
                done = false;
                propagationDone = false;
                canceled = false;
                statMessage = null;
                return new generateParticleAndMonteCarloBase();
            }
        };
        worker.start();
    }

    /**
     * start the task to run a particle.
     */
    public void propagate() {
        final SwingWorker workerRun = new SwingWorker() {
            public Object construct() {
                current = 0;
                done = false;
                propagationDone = false;
                canceled = false;
                statMessage = null;
                return new runParticle();
            }
        };
        workerRun.start();
    }

    /**
     * find out how much work needs
     * to be done.
     */
    public int getLengthOfTask() {
        return lengthOfTask;
    }

    /**
     * Called from ProgressBarDemo to find out how much has been done.
     */
    public int getCurrent() {
        return current;
    }

    public void stop() {
        canceled = true;
        statMessage = null;
    }

    /**
     * Called from ProgressBarDemo to find out if the task has completed.
     */
    public boolean isDone() {
        return done;
    }

    /** Tell outside objects whether runEvent is done or not.*/
    public boolean isPropagationDone() {
	return propagationDone;
    }

    public void setPropagationFlagToFalse() {
	propagationDone = false;
    }

    /**
     * Returns the most recent status message, or null
     * if there is no current status message.
     */
    public String getMessage() {
        return statMessage;
    }

    /**
     * The actual long running task to generate Particle.  
     * This runs in a SwingWorker thread.
     */
    class generateParticleAndMonteCarloBase {
        generateParticleAndMonteCarloBase() {
            while (!canceled && !done) {

	//
	// Select interactions and decay
	//

	// CC interaction
	int doCC = 0;
	if(swingRun.nuCcInt.getSelectedObjects()!=null) {
	    swingRun.messageField.setText("Select nu CC object");
	    doCC = 1;
	}
	current += 50;

	// NC interaction
	int doNC = 0;
	if(swingRun.nuNcInt.getSelectedObjects()!=null){
	    swingRun.messageField.setText("Select nu NC object");
	    doNC = 1;
	}
	current += 50;

	// muon Bremsstrahlung
	int doMuBrems = 0;
	if(swingRun.muBrems.getSelectedObjects()!=null){
	    swingRun.messageField.setText("Select mu Bremsstralung object");
	    doMuBrems = 1;
	}
	current += 50;

	// tau Bremsstrahlung
	int doTauBrems = 0;
	if(swingRun.tauBrems.getSelectedObjects()!=null){
	    swingRun.messageField.setText("Select tau Bremsstralung object");
	    doTauBrems = 1;
	}
	current += 50;

	// muon KnockOn
	int doMuKnock = 0;
	if(swingRun.muKnock.getSelectedObjects()!=null){
	    swingRun.messageField.setText("Select mu Knock-on electron object");
	    doMuKnock = 1;
	}
	current += 50;

	// tau KnockOn
	int doTauKnock = 0;
	if(swingRun.tauKnock.getSelectedObjects()!=null){
	    swingRun.messageField.setText("Select tau Knock-on electron object");
	    doTauKnock = 1;
	}
	current += 50;

	// mu e+e- pair creation
	int doMu2e = 0;
	if(swingRun.muToePair.getSelectedObjects()!=null){
	    swingRun.messageField.setText("Select mu e+e- object");
	    doMu2e = 1;
	}
	current += 50;

	//tau e+e- pair creation
	int doTau2e = 0;
	if(swingRun.tauToePair.getSelectedObjects()!=null){
	    swingRun.messageField.setText("Select tau e+e- object");
	    doTau2e = 1;
	}
	current += 50;

	// heavy lepton pair creations - Not selected in this class.
	int doMu2mu = 0;
	int doTau2mu = 0;
	int doMu2tau = 0;
	int doTau2tau = 0;

	//muon photonuclear reaction
	int doMuPN = 0;
	if(swingRun.muPhotoNucl.getSelectedObjects()!=null){
	    swingRun.messageField.setText("Select mu Photo Nuclear object");
	    doMuPN = 1;
	}
	current += 50;

	//tau photonuclear reaction
	int doTauPN = 0;
	if(swingRun.tauPhotoNucl.getSelectedObjects()!=null){
	    swingRun.messageField.setText("Select tau Photo Nuclear object");
	    doTauPN = 1;
	}

	// decay
	int doMuDecay = 0;
	int doTauDecay = 0;
	if(swingRun.decay.getSelectedObjects()!=null){
	    doMuDecay = 1;
	    doTauDecay = 1;
	}

	swingRun.messageField.setText("Now invoking the JulietEventGenerator...");
	try{
	generator = new JulietEventGenerator(swingRun.flavor, swingRun.doublet, swingRun.energy, 0,
                                doCC, doNC, doMuBrems, doTauBrems, 
                                doMuKnock, doTauKnock, doMu2e, doTau2e,
                                doMu2mu, doTau2mu, doMu2tau, doTau2tau,
                                doMuPN, doTauPN, doMuDecay, doTauDecay, 0);
	swingRun.messageField.setText("done!");
	}catch(IOException e) {
            swingRun.messageField.setText("IOException: " + e.getMessage());
	}


	current += 50;
        swingRun.particleMass.setText("Mass "+ 
				      Particle.particleMasses[swingRun.flavor][swingRun.doublet] 
				      + "[GeV]");
        swingRun.particleEnergy.setText("Energy "+ swingRun.energy + "[GeV]");
        swingRun.particleName.setText(Particle.particleName(swingRun.flavor,swingRun.doublet));

	current = lengthOfTask;
	if (current >= lengthOfTask) {
	    done = true;
	    current = lengthOfTask;
	}
	statMessage = "Completed " + current +
                                  " out of " + lengthOfTask + ".";
            }
        }
    }



    /**
     * The actual long running task to generate Particle.  
     * This runs in a SwingWorker thread.
     */
    class runParticle {
        runParticle() {
            while (!canceled && !done) {

		swingRun.julietIconLabel.setIcon(swingRun.julietRunIcon);

		swingRun.messageField.setText("Configure the geometry...");
		generator.definePropagationGeometry(swingRun.x_ice3,swingRun.y_ice3,
						    swingRun.z_ice3,swingRun.nadirAngle_ice3,
						    swingRun.azimuthAngle_ice3);
		generator.configurePropagationGeometry();
		current += 10;
		try{
		    Thread.sleep(2000); //sleep for a two second
		}catch(InterruptedException e){
		    swingRun.messageField.setText("RUN clashed");
		}

		// Run a single event
		swingRun.julietIconLabel.setIcon(swingRun.julietRunIcon2);
		swingRun.messageField.setText("A JULIeT RUN starts.");
		generator.runSingleEvent();
		swingRun.messageField.setText("The Journey ends as " +
 	        generator.event.propParticle.particleName(
			  generator.event.propParticle.getFlavor(),
			  generator.event.propParticle.getDoublet()));
		try{
		    Thread.sleep(2000); //sleep for a two second
		}catch(InterruptedException e){
		    swingRun.messageField.setText("RUN clashed");
		}

		swingRun.messageField.setText("Energy Deposit " + 
					     generator.event.getCascadeTotalEnergy() + "[GeV]");

		current = lengthOfTask;

		if (current >= lengthOfTask) {
		    done = true;
		    propagationDone = true;
		    current = lengthOfTask;
		}
		statMessage = "Completed " + current +
		    " out of " + lengthOfTask + ".";
            }
        }
    }

}


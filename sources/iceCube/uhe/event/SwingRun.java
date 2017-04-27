package iceCube.uhe.event;

import iceCube.uhe.particles.*;

import java.io.*;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import java.net.URL;


/** Run the JulietEventGenerator  with the Swing User Interface */
public class SwingRun implements ActionListener {

    final static int START_INDEX = 0;
    // Parameters related Propagation Particle and its track geometry
    int flavor = 0;
    int doublet = 0;
    double energy = 0.0;
    double logE = 7.0;
    double x_ice3 = 0.0;
    double y_ice3 = 0.0;
    double z_ice3 = 0.0;
    double nadirAngle_ice3 = 0.0;
    double azimuthAngle_ice3 = 0.0;


    // The Generator Panel
    JPanel particleGenerator;
    JLabel iceCubeIconLabel;
    TrackGeometryGraphics trackGraph;
    JLabel energyLabel;
    JLabel nadirLabel;
    JLabel azimuthLabel;
    JLabel xLabel;
    JLabel yLabel;
    JLabel zLabel;
    JTextField inputEnergy;
    JTextField inputNadirAngle;
    JTextField inputAzimuthAngle;
    JTextField inputX;
    JTextField inputY;
    JTextField inputZ;
    JComboBox particleFlavor = null;
    JComboBox particleDoublet = null;

    // The Particle Display Panel
    JPanel particleDisplay;
    JLabel particleMass, particleName, particleEnergy;

    // The Intetraction Panel
    JPanel interactionGenerator;
    JCheckBox nuCcInt; // Charged Current
    JCheckBox nuNcInt; // Neutral Current
    JCheckBox muBrems, tauBrems;        // Bremsstrahlug
    JCheckBox muKnock, tauKnock;        // Bremsstrahlug
    JCheckBox muToePair, tauToePair;    // Pair Creation
    JCheckBox muPhotoNucl, tauPhotoNucl;    // Photo Nuclear Interaction
    JCheckBox decay; // Decay
    JButton generateParticle;
    JButton runParticle;

    //The Message Panel
    JPanel messageDisplay;
    JTextField messageField;
    private JProgressBar progressBar;
    private Timer timer;
    public final static int HALF_SECOND = 500;

    // JulietTask: Interface to JulietEventGenerator
    public JulietTask  julietTask;


    // The Particle RUN Display
    JPanel particleRunDisplay;
    JLabel julietIconLabel;
    ImageIcon julietIcon;
    ImageIcon julietRunIcon;
    ImageIcon julietRunIcon2;
    JLabel messageLabel;
    ParticleTravelGraphics particleTravel;
    InteractionGraphics interactionFreq;

    // Main Panel
    JPanel mainPanel;

    //constructor
    public SwingRun( ) {

	// Create the particle generator and display panels.
	particleGenerator = new JPanel();
	particleGenerator.setBackground(new Color(0xECFFD9));
	particleDisplay = new JPanel();
	particleDisplay.setBackground(new Color(0xECFFD9));
	interactionGenerator = new JPanel();
	interactionGenerator.setBackground(new Color(0xECFFD9));
	messageDisplay = new JPanel();
	messageDisplay.setBackground(new Color(0xffffff));
	particleRunDisplay = new JPanel();
	particleRunDisplay.setBackground(new Color(0xffffff));

        // Add various widgets to the sub panels.
        addWidgets();
        // Create the main panel to contain the two sub panels.
        mainPanel = new JPanel();
        mainPanel.setBorder(BorderFactory.createEmptyBorder(5,5,5,5));
        GridBagLayout gridbagMain = new GridBagLayout();
        GridBagConstraints cMain = new GridBagConstraints();
	mainPanel.setLayout(gridbagMain);
	cMain.fill = GridBagConstraints.HORIZONTAL;

        // Add the select and display panels to the main panel.
        cMain.gridx = 0;
        cMain.gridy = 0;
        cMain.gridheight = 2;
        gridbagMain.setConstraints(particleGenerator, cMain);
        mainPanel.add(particleGenerator);

        cMain.gridx = 0;
        cMain.gridy = 2;
        cMain.gridheight = 1;
        gridbagMain.setConstraints(particleDisplay, cMain);
        mainPanel.add(particleDisplay);

        cMain.gridx = 0;
        cMain.gridy = 3;
        cMain.gridheight = 2;
        gridbagMain.setConstraints(interactionGenerator, cMain);
        mainPanel.add(interactionGenerator);

        cMain.gridx = 0;
        cMain.gridy = 5;
        cMain.gridheight = 1;
        gridbagMain.setConstraints(messageDisplay, cMain);
        mainPanel.add(messageDisplay);

        cMain.gridx = 0;
        cMain.gridy = 6;
        cMain.gridheight = 3;
        gridbagMain.setConstraints(particleRunDisplay, cMain);
        mainPanel.add(particleRunDisplay);
    }


    // Create and add the widgets to select particle spieces and
    // its propaty.
    private void addWidgets(){

	// Generate JulietTask
        julietTask = new JulietTask(this);

	// Get the images and put them into ImageIcon.
	String imageIceCube = "images/icecube-logo.jpg";
	URL iceCubeIconURL =ClassLoader.getSystemResource(imageIceCube);
	ImageIcon iceCubeIcon = new ImageIcon(iceCubeIconURL);

	// The widgets to select primary particle.
	iceCubeIconLabel = new JLabel();
	iceCubeIconLabel.setBackground(new Color(0x666699));
        iceCubeIconLabel.setHorizontalAlignment(JLabel.CENTER);
        iceCubeIconLabel.setVerticalAlignment(JLabel.CENTER);
        iceCubeIconLabel.setVerticalTextPosition(JLabel.CENTER);
        iceCubeIconLabel.setHorizontalTextPosition(JLabel.CENTER);
        iceCubeIconLabel.setBorder(BorderFactory.createCompoundBorder(
                            BorderFactory.createLoweredBevelBorder(),
                            BorderFactory.createEmptyBorder(5,5,5,5)));

        iceCubeIconLabel.setBorder(BorderFactory.createCompoundBorder(
                            BorderFactory.createEmptyBorder(0,0,10,0),
                            iceCubeIconLabel.getBorder()));
	iceCubeIconLabel.setIcon(iceCubeIcon);
	iceCubeIconLabel.setText("");

	String[] flavorString = 
	{"e   Flavor","mu  Flavor","tau Flavor","Hadron"};
	particleFlavor = new JComboBox(flavorString);
	particleFlavor.setBackground(new Color(0x666699));
	particleFlavor.setForeground(new Color(0x65d0f1));
	particleFlavor.setSelectedIndex(START_INDEX);

	String[] doubletString = {"Doblet 1","Doublet 2"};
	particleDoublet = new JComboBox(doubletString);
	particleDoublet.setBackground(new Color(0x666699));
	particleDoublet.setForeground(new Color(0x65d0f1));
	particleDoublet.setSelectedIndex(START_INDEX);

	energyLabel = new JLabel("log (Energy [GeV])");
	energyLabel.setForeground(new Color(0x000070));
	inputEnergy = new JTextField(6);

	nadirLabel = new JLabel("Nadir [Deg]");
	nadirLabel.setForeground(new Color(0x000070));
	inputNadirAngle = new JTextField(6);

	azimuthLabel = new JLabel("Azimuth [Deg]");
	azimuthLabel.setForeground(new Color(0x000070));
	inputAzimuthAngle = new JTextField(6);

	xLabel = new JLabel("X in ice3 [cm]");
	xLabel.setForeground(new Color(0x000070));
	inputX = new JTextField(10);

	yLabel = new JLabel("Y in ice3 [cm]");
	yLabel.setForeground(new Color(0x000070));
	inputY = new JTextField(10);

	zLabel = new JLabel("Z in ice3 [cm]");
	zLabel.setForeground(new Color(0x000070));
	inputZ = new JTextField(10);

	//Add the labels to the panal.
        particleGenerator.setBorder(BorderFactory.createCompoundBorder(
	      BorderFactory.createTitledBorder(
	      "Select Particle"),
                BorderFactory.createEmptyBorder(5,5,5,5)));
        GridBagLayout gridbag = new GridBagLayout();
        GridBagConstraints c = new GridBagConstraints();
	particleGenerator.setLayout(gridbag);
	c.fill = GridBagConstraints.HORIZONTAL; 
	c.weightx = 0.5;

        c.gridx = 0;
        c.gridy = 0;
        c.gridheight = 8;
        gridbag.setConstraints(iceCubeIconLabel, c);
	particleGenerator.add(iceCubeIconLabel);

	trackGraph = new TrackGeometryGraphics(julietTask);
        trackGraph.setPreferredSize(new Dimension(240,120));
        c.gridx = 1;
        c.gridy = 0;
        c.gridheight = 8;
        gridbag.setConstraints(trackGraph, c);
	particleGenerator.add(trackGraph);

        c.gridx = 3;
        c.gridy = 0;
        c.gridheight = 1;
        gridbag.setConstraints(particleFlavor, c);
	particleGenerator.add(particleFlavor);

        c.gridx = 4;
        c.gridy = 0;
        gridbag.setConstraints(particleDoublet, c);
	particleGenerator.add(particleDoublet);

        c.gridx = 3;
        c.gridy = 1;
        gridbag.setConstraints(energyLabel, c);
	particleGenerator.add(energyLabel);

        c.gridx = 4;
        c.gridy = 1;
        gridbag.setConstraints(inputEnergy, c);
	particleGenerator.add(inputEnergy);

        c.gridx = 3;
        c.gridy = 2;
        gridbag.setConstraints(nadirLabel, c);
	particleGenerator.add(nadirLabel);

        c.gridx = 4;
        c.gridy = 2;
        gridbag.setConstraints(inputNadirAngle, c);
	particleGenerator.add(inputNadirAngle);

        c.gridx = 3;
        c.gridy = 3;
        gridbag.setConstraints(azimuthLabel, c);
	particleGenerator.add(azimuthLabel);

        c.gridx = 4;
        c.gridy = 3;
        gridbag.setConstraints(inputAzimuthAngle, c);
	particleGenerator.add(inputAzimuthAngle);

        c.gridx = 3;
        c.gridy = 4;
        gridbag.setConstraints(xLabel, c);
	particleGenerator.add(xLabel);

        c.gridx = 4;
        c.gridy = 4;
        gridbag.setConstraints(inputX, c);
	particleGenerator.add(inputX);

        c.gridx = 3;
        c.gridy = 5;
        gridbag.setConstraints(yLabel, c);
	particleGenerator.add(yLabel);

        c.gridx = 4;
        c.gridy = 5;
        gridbag.setConstraints(inputY, c);
	particleGenerator.add(inputY);

        c.gridx = 3;
        c.gridy = 6;
        gridbag.setConstraints(zLabel, c);
	particleGenerator.add(zLabel);

        c.gridx = 4;
        c.gridy = 6;
        gridbag.setConstraints(inputZ, c);
	particleGenerator.add(inputZ);

	// The widgets to display particle propaties.
	particleName = new JLabel("Name",SwingConstants.LEFT);
	particleMass = new JLabel("Mass",SwingConstants.LEFT);
	particleEnergy = new JLabel("Energy",SwingConstants.LEFT);

	//Add the labels to the panal.
        particleDisplay.setBorder(BorderFactory.createCompoundBorder(
                BorderFactory.createTitledBorder("Particle Propaty"),
                BorderFactory.createEmptyBorder(5,5,5,5)));
	particleDisplay.setLayout(new GridLayout(1,0));
	particleDisplay.add(particleName);
	particleDisplay.add(particleMass);
	particleDisplay.add(particleEnergy);

	// The widgets to select Interactions
	nuCcInt = new JCheckBox("nu CC");
	nuCcInt.setBackground(new Color(0xECFFD9));
	nuNcInt = new JCheckBox("nu NC");
	nuNcInt.setBackground(new Color(0xECFFD9));
	muBrems = new JCheckBox("mu Bremss");
	muBrems.setBackground(new Color(0xECFFD9));
	tauBrems = new JCheckBox("tau Bremss");
	tauBrems.setBackground(new Color(0xECFFD9));
	muKnock = new JCheckBox("mu Knock-on");
	muKnock.setBackground(new Color(0xECFFD9));
	tauKnock = new JCheckBox("tau Knock-on");
	tauKnock.setBackground(new Color(0xECFFD9));
	muToePair = new JCheckBox("mu e+e-");
	muToePair.setBackground(new Color(0xECFFD9));
	tauToePair = new JCheckBox("tau e+e-");
	tauToePair.setBackground(new Color(0xECFFD9));
	muPhotoNucl = new JCheckBox("mu Photo-nucl");
	muPhotoNucl.setBackground(new Color(0xECFFD9));
	tauPhotoNucl = new JCheckBox("tau Photo-nucl");
	tauPhotoNucl.setBackground(new Color(0xECFFD9));
	decay = new JCheckBox("mu/tau decay");
	decay.setBackground(new Color(0xECFFD9));

	//Add the labels to the panal.
        interactionGenerator.setBorder(BorderFactory.createCompoundBorder(
                BorderFactory.createTitledBorder("Select Interaction"),
                BorderFactory.createEmptyBorder(5,5,5,5)));
        GridBagLayout gridbagInt = new GridBagLayout();
        GridBagConstraints cInt = new GridBagConstraints();
	interactionGenerator.setLayout(gridbagInt);
	cInt.fill = GridBagConstraints.HORIZONTAL; 
	cInt.weightx = 0.5;

        cInt.gridx = 0;
        cInt.gridy = 0;
        cInt.gridheight = 1;
        gridbagInt.setConstraints(nuCcInt, cInt);
	interactionGenerator.add(nuCcInt);
	nuCcInt.setSelected(false);

        cInt.gridx = 1;
        cInt.gridy = 0;
        gridbagInt.setConstraints(nuNcInt, cInt);
	interactionGenerator.add(nuNcInt);
	nuNcInt.setSelected(false);

        cInt.gridx = 0;
        cInt.gridy = 1;
        gridbagInt.setConstraints(muBrems, cInt);
	interactionGenerator.add(muBrems);
	muBrems.setSelected(true);

        cInt.gridx = 1;
        cInt.gridy = 1;
        gridbagInt.setConstraints(tauBrems, cInt);
	interactionGenerator.add(tauBrems);
	tauBrems.setSelected(true);

        cInt.gridx = 2;
        cInt.gridy = 1;
        gridbagInt.setConstraints(muKnock, cInt);
	interactionGenerator.add(muKnock);
	muKnock.setSelected(true);

        cInt.gridx = 3;
        cInt.gridy = 1;
        gridbagInt.setConstraints(tauKnock, cInt);
	interactionGenerator.add(tauKnock);
	tauKnock.setSelected(true);

        cInt.gridx = 0;
        cInt.gridy = 2;
        gridbagInt.setConstraints(muToePair, cInt);
	interactionGenerator.add(muToePair);
	muToePair.setSelected(true);

        cInt.gridx = 1;
        cInt.gridy = 2;
        gridbagInt.setConstraints(tauToePair, cInt);
	interactionGenerator.add(tauToePair);
	tauToePair.setSelected(true);

        cInt.gridx = 2;
        cInt.gridy = 2;
        gridbagInt.setConstraints(muPhotoNucl, cInt);
	interactionGenerator.add(muPhotoNucl);
	muPhotoNucl.setSelected(true);

        cInt.gridx = 3;
        cInt.gridy = 2;
        gridbagInt.setConstraints(tauPhotoNucl, cInt);
	interactionGenerator.add(tauPhotoNucl);
	tauPhotoNucl.setSelected(true);

        cInt.gridx = 3;
        cInt.gridy = 0;
        gridbagInt.setConstraints(decay, cInt);
	interactionGenerator.add(decay);
	decay.setSelected(true);

	generateParticle = new JButton("Step1: Generate Interactions");
	generateParticle.setBackground(new Color(0x72CE68));
	generateParticle.setActionCommand("GENERATE");
        cInt.gridx = 4;
        cInt.gridy = 0;
        cInt.gridheight = 1;
        gridbagInt.setConstraints(generateParticle, cInt);
	interactionGenerator.add(generateParticle);

	runParticle =  new JButton("Step2: Run the Particle ...");
	runParticle.setBackground(new Color(0x72CE68));
	runParticle.setActionCommand("RUN");
        cInt.gridx = 4;
        cInt.gridy = 2;
        cInt.gridheight = 1;
        gridbagInt.setConstraints(runParticle, cInt);
	interactionGenerator.add(runParticle);

	// The widgets to display messages.
        messageDisplay.setBorder(BorderFactory.createCompoundBorder(
                BorderFactory.createTitledBorder("Message Display"),
                BorderFactory.createEmptyBorder(5,5,5,5)));
	messageDisplay.setLayout(new GridLayout(0,1));
	messageField = new JTextField(64);

        progressBar = new JProgressBar(0, julietTask.getLengthOfTask());
	progressBar.setBackground(new Color(0xffff00));
        progressBar.setValue(0);
        progressBar.setStringPainted(true);


	//Add the labels to the panal.
	messageDisplay.add(messageField);
	messageDisplay.add(progressBar);

	// Get the images and put them into ImageIcon.
	String imageJULIeT = "images/romeoAndJulietSmall.jpg";
	URL julietIconURL =ClassLoader.getSystemResource(imageJULIeT);
	julietIcon = new ImageIcon(julietIconURL);
	String imageJULIeTRun = "images/runJuliet.png";
	URL julietRunIconURL =ClassLoader.getSystemResource(imageJULIeTRun);
	julietRunIcon = new ImageIcon(julietRunIconURL);
	String imageJULIeTRun2 = "images/runJuliet2.png";
	URL julietRunIconURL2 =ClassLoader.getSystemResource(imageJULIeTRun2);
	julietRunIcon2 = new ImageIcon(julietRunIconURL2);

	// The widgets to display an simulating event
	julietIconLabel = new JLabel();
	julietIconLabel.setIcon(julietIcon);
	julietIconLabel.setText("");
	particleTravel = new ParticleTravelGraphics(julietTask);
	interactionFreq = new InteractionGraphics(julietTask);

	//Add the labels to the panal.
        particleRunDisplay.setBorder(BorderFactory.createCompoundBorder(
                BorderFactory.createTitledBorder("Event Display"),
                BorderFactory.createEmptyBorder(5,5,5,5)));
	particleRunDisplay.setLayout(new GridLayout(1,0));
	particleRunDisplay.add(julietIconLabel);
	particleRunDisplay.add(particleTravel);
	particleRunDisplay.add(interactionFreq);

        //Create a timer.
        timer = new Timer(HALF_SECOND, new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
                progressBar.setValue(julietTask.getCurrent());
		trackGraph.repaint();
		particleTravel.repaint();
		interactionFreq.repaint();
                String s = julietTask.getMessage();
                if (s != null) {
		    if (progressBar.isIndeterminate()) {
			progressBar.setIndeterminate(false);
			//progressBar.setString(null); //display % string
		    }
		}
		if (julietTask.isDone()) {
                    Toolkit.getDefaultToolkit().beep();
                    timer.stop();
                    progressBar.setValue(progressBar.getMaximum());
                }
            }
        });


	// Listen to events from Generatorbuttun.
	generateParticle.addActionListener(this);
	runParticle.addActionListener(this);

    }



    // Implementation of ActionListener interface.
    public void actionPerformed(ActionEvent event) {
	//progressBar.setIndeterminate(true);

	if("GENERATE".equals(event.getActionCommand())){ // Generate a particle
	    flavor = particleFlavor.getSelectedIndex();
	    doublet = particleDoublet.getSelectedIndex();
	    logE = Double.parseDouble(inputEnergy.getText());
	    energy = Math.pow(10.0,logE);
	    if(Particle.isValidFlavor(flavor) && Particle.isValidDoublet(doublet)){
		if(logE>=6.0){
		    julietTask.generate();
		    timer.start();
		}else{
		    messageField.setText("Energy Must be Grater than 1.0e6 GeV");
		}
	    }else{
		particleName.setText("Invalid flavor/doublet");
	    }
	}

	if("RUN".equals(event.getActionCommand())){ // Generate a particle
	    nadirAngle_ice3 = Double.parseDouble(inputNadirAngle.getText());
	    azimuthAngle_ice3 = Double.parseDouble(inputAzimuthAngle.getText());
	    x_ice3 = Double.parseDouble(inputX.getText());
	    y_ice3 = Double.parseDouble(inputY.getText());
	    z_ice3 = Double.parseDouble(inputZ.getText());
	    julietTask.propagate();
	    timer.start();
	}
    }
	   

    // main method
    public static void main(String[] args) {
        javax.swing.SwingUtilities.invokeLater(new Runnable() {
            public void run() {
		SwingRun swingRun = new SwingRun();
		JFrame julietFrame = new JFrame("JULIeT");

		// Set the look and feel.
		try {
		    UIManager.setLookAndFeel(
				     UIManager.getCrossPlatformLookAndFeelClassName());
		} catch(Exception e) {}

		julietFrame.setContentPane(swingRun.mainPanel);
		julietFrame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		// Show the generator.
		julietFrame.pack();
		julietFrame.setVisible(true);
	    }
        });
   }
}

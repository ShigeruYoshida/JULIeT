package iceCube.uhe.event;

import iceCube.uhe.particles.*;

import java.io.*;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import java.net.URL;


/** Run the JulietEventGenerator  with the Swing User Interface */
public class SwingRunIceTray extends SwingRun {

    JFrame julietFrame;

    //constructor
    public SwingRunIceTray( ) {
        super();
	julietFrame = new JFrame("JULIeT");
    }

    public void configure(){
        // Set the look and feel.
        try {
            UIManager.setLookAndFeel(
                UIManager.getCrossPlatformLookAndFeelClassName());
        } catch(Exception e) {}
        
        julietFrame.setContentPane(mainPanel);
    }

    public void makeItVisible(){
        // Exit when the window is closed.
        julietFrame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE); 
        
        // Show the converter.
        julietFrame.pack();
        julietFrame.setVisible(true);
    }

    // main method
    public static void main(String[] args) {
	SwingRunIceTray swingRunIceTray = new SwingRunIceTray();

	// Set the look and feel.
        swingRunIceTray.configure();
        swingRunIceTray.makeItVisible();
   }
}

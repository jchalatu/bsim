package BSimCrossFeeding;

import bsim.BSim;  
import bsim.BSimChemicalField;
import bsim.BSimTicker;
import bsim.capsule.BSimCapsuleBacterium;
import bsim.capsule.Mover;
import bsim.capsule.RelaxationMoverGrid;

import java.util.ArrayList;
import java.util.Random;

public class CrossFeedingTicker extends BSimTicker {

    BSim sim;
    // logs data about time taken by ticker every LOG_INTERVAL timesteps
    final int LOG_INTERVAL;
    public static boolean WITH_GROWTH;
    static Random bacRng;
    //growth rate standard deviation
    public final double growth_stdv;
    //growth rate mean
    public final double growth_mean;
    //elongation threshold standard deviation
    public final double length_stdv;
    //elongation threshold mean
    public final double length_mean;

    final ArrayList<CrossFeedingBacterium> bacA;
    final ArrayList<CrossFeedingBacterium> bacB;
    final ArrayList<BSimCapsuleBacterium> bacteriaAll;
    final ArrayList<CrossFeedingBacterium> bac_bornA;
    final ArrayList<CrossFeedingBacterium> bac_bornB;
    final ArrayList<CrossFeedingBacterium> bac_deadA;
    final ArrayList<CrossFeedingBacterium> bac_deadB;
    
    public static BSimChemicalField amino_acid_A;
    public static BSimChemicalField amino_acid_B;

    // internal machinery - don't worry about this
    final Mover mover;

    public CrossFeedingTicker(BSim sim, ArrayList<CrossFeedingBacterium> bacA, ArrayList<CrossFeedingBacterium> bacB,
    		ArrayList<BSimCapsuleBacterium> bacteriaAll, int LOG_INTERVAL, Random bacRng, 
    		double growth_stdv, double growth_mean, double length_stdv, double length_mean,
    		BSimChemicalField amino_acid_A, BSimChemicalField amino_acid_B) {
        this.sim = sim;
        this.LOG_INTERVAL = LOG_INTERVAL;
        this.bacRng = bacRng; //random number generator
        this.growth_stdv = growth_stdv;
        this.growth_mean = growth_mean;
        this.length_stdv = length_stdv;
        this.length_mean = length_mean;
        this.bacA = bacA;
        this.bacB = bacB;
        this.bacteriaAll = bacteriaAll;
        bac_bornA = new ArrayList();
        bac_bornB = new ArrayList();
        bac_deadA = new ArrayList();
        bac_deadB = new ArrayList();
        mover = new RelaxationMoverGrid(bacteriaAll, sim);
        
        CrossFeedingTicker.amino_acid_A = amino_acid_A;
        CrossFeedingTicker.amino_acid_B = amino_acid_B;
    }
    
    /** Sets the flag for growth. **/
    public void setGrowth(boolean b) { WITH_GROWTH = b; }
    
    /** Function for bacteria growth activities. */
    public void growBacteria( ArrayList<CrossFeedingBacterium> bac, ArrayList<CrossFeedingBacterium> bacBorn ) {
    	Random bacRng = new Random(); 			// Random number generator
    	bacRng.setSeed(50); 					// Initializes random number generator
    	
        for (CrossFeedingBacterium b : bac) { 				// Loop over bac array
            b.grow();

            // Divide if grown past threshold
            if (b.L >= b.L_th) {
            	CrossFeedingBacterium daughter = b.divide();
            	
            	bacBorn.add(daughter);  		// Add daughter to newborn class, 'mother' keeps her status
            }
        }
        
        bac.addAll(bacBorn); 					// Adds all the newborn daughters to sub-population
        bacteriaAll.addAll(bacBorn); 			// Adds all the newborn daughters to total population
        
        // Allow daughter cells to grow 
        for ( CrossFeedingBacterium b : bacBorn ) {
        	
            // Assigns a division length to each bacterium according to a normal distribution
            double lengthThreshold = length_stdv*bacRng.nextGaussian()+length_mean;
            b.setElongationThreshold(lengthThreshold);
        }
        
        bacBorn.clear(); 						// Cleared for next time-step
    }
    
    /** Function to remove bacteria due to cell death or by boundary. */
    public void removeBacteria( BSim sim, ArrayList<CrossFeedingBacterium> bac, ArrayList<CrossFeedingBacterium> bac_dead ) {

        for (CrossFeedingBacterium b : bac) {
        	
            // Kick out if past any boundary
        	// Bacteria out of bounds = dead
            if(b.position.x < 0 || b.position.x > sim.getBound().x || b.position.y < 0 || b.position.y > sim.getBound().y || b.position.z < 0 || b.position.z > sim.getBound().z){
            	bac_dead.add(b);
            } 
            // Remove cell after it shrinks and becomes too small
            if ( b.L <= 1 ) {
            	bac_dead.add(b);
            }
        }
        bac.removeAll(bac_dead);
        
        // Remove from total population
        bacteriaAll.removeAll(bac_dead);
        bac_dead.clear();
    }
    
    // This one is a bit long too. Let's break it up
    // 1. Begins an "action" -> this represents one timestep
    // 2. Tells each bacterium to perform their action() function
    // 3. Updates each chemical field in the simulation
    // 4. Bacteria are then told to grow()
    // 5. bacteria which are longer than their threshold are told to divide()
    // 6. forces are applied and bacteria move around
    // 7. bacteria which are out of bounds are removed from the simulation
    @Override
    public void tick() {
        // increase lifetimes of cells
        for (CrossFeedingBacterium b : bacA) {
            b.lifetime++;
        }
        for (CrossFeedingBacterium b : bacB) {
            b.lifetime++;
        }
    	
        /********************************************** Action */
    	
        long startTimeAction = System.nanoTime(); 	// Wall-clock time, for diagnosing timing

        // Calculates and stores the midpoint of the cell.
        // Bacteria does action at each time step
        for(BSimCapsuleBacterium b : bacteriaAll) {
            b.action(); 							
        }

        // Outputs how long each step took, once every log interval.
        long endTimeAction = System.nanoTime();
        if((sim.getTimestep() % LOG_INTERVAL) == 0) {
            System.out.println("Action update for " + bacteriaAll.size() + " bacteria took " + (endTimeAction - startTimeAction)/1e6 + " ms.");
        }  

        /********************************************** Chemical fields */
        
        startTimeAction = System.nanoTime();
        
        // Update the amino acid field
        amino_acid_A.update(); 
        amino_acid_B.update();

        endTimeAction = System.nanoTime();
        if((sim.getTimestep() % LOG_INTERVAL) == 0) {
            System.out.println("Chemical field update took " + (endTimeAction - startTimeAction)/1e6 + " ms.");
        }

        /********************************************** Growth related activities if enabled. */
        if (WITH_GROWTH) {

            // ********************************************** Growth and division
            startTimeAction = System.nanoTime(); 	// Start action timer

            growBacteria( bacA, bac_bornA );		// For sub-population A
            growBacteria( bacB, bac_bornB );		// For sub-population B

            // Prints out information about bacteria when u want it to
            endTimeAction = System.nanoTime();
            if ((sim.getTimestep() % LOG_INTERVAL) == 0) {
                System.out.println("Growth and division took " + (endTimeAction - startTimeAction) / 1e6 + " ms.");
            }
            
            /********************************************** Neighbor interactions */
            
            startTimeAction = System.nanoTime();

            mover.move();

            endTimeAction = System.nanoTime();
            if ((sim.getTimestep() % LOG_INTERVAL) == 0) {
                System.out.println("Wall and neighbour interactions took " + (endTimeAction - startTimeAction) / 1e6 + " ms.");
            }

            /********************************************** Boundaries/removal */
            startTimeAction = System.nanoTime();
            
            removeBacteria( sim, bacA, bac_deadA );		// For sub-population A
            removeBacteria( sim, bacB, bac_deadB );		// For sub-population B
            
            endTimeAction = System.nanoTime();
            if ((sim.getTimestep() % LOG_INTERVAL) == 0) {
                System.out.println("Death and removal took " + (endTimeAction - startTimeAction) / 1e6 + " ms.");
            }
        }

        startTimeAction = System.nanoTime();

        endTimeAction = System.nanoTime();
        if ((sim.getTimestep() % LOG_INTERVAL) == 0) {
            System.out.println("Switch took " + (endTimeAction - startTimeAction) / 1e6 + " ms.");
        }

    }
}

package BSimCrossProtection;

import bsim.BSim;
import bsim.BSimTicker;
import bsim.capsule.BSimCapsuleBacterium;
import bsim.capsule.Mover;
import bsim.capsule.RelaxationMoverGrid;

import java.util.ArrayList;
import java.util.Random;

public class BasicTicker extends BSimTicker {

    BSim sim;
    // logs data about time taken by ticker every LOG_INTERVAL timesteps
    final int LOG_INTERVAL;
    /** Whether to enable growth in the ticker etc. or not... */
    private static boolean WITH_GROWTH = true;
    
    static Random bacRng;
    //growth rate standard deviation
    public final double growth_stdv;
    //growth rate mean
    public final double growth_mean;
    //elongation threshold standard deviation
    public final double length_stdv;
    //elongation threshold mean
    public final double length_mean;

    final ArrayList<CrossProtectionBacterium> bacA;
    final ArrayList<CrossProtectionBacterium> bacB;
    final ArrayList<BSimCapsuleBacterium> bacteriaAll;
    final ArrayList<CrossProtectionBacterium> bac_bornA;
    final ArrayList<CrossProtectionBacterium> bac_bornB;
    final ArrayList<CrossProtectionBacterium> bac_deadA;
    final ArrayList<CrossProtectionBacterium> bac_deadB;
    
    public static ChemicalField antibioticA;
    public static ChemicalField antibioticB;
    
	// Initial Conditions
	/** Used to determine the toxin distribution for the simulation. */
	private int toxin_condition;
	
	/** Toxins flow in from the left boundary. */
	private static final int FLOW_IN = 1;
	/** Simulation increases uniform distribution of toxin. */
	private static final int UNIFORM = 2;
	/** Steady state concentration level. */
	private double initial_conc = 310;
    /** How much toxin A increases per hour. */
    private double toxin_increment_A = 50;
    /** How much toxin B increases per hour. */
    private double toxin_increment_B = 50;
    /** Delay for toxin A. */
    private int delayA = 0;
    /** Delay for toxin B. */
    private int delayB = 0;
    /** Delay counter for toxin A. */
    private int delayCountA = 0;
    /** Delay counter for toxin B. */
    private int delayCountB = 0;
    
    // For a single screen
    final int MIXED_CONC = 1;
    final int CHECKER_BOARD = 2;
    int SINGLE_SCREEN = 2;
    
    /** Defines the progress of the chemical field flowing through the boundary on the x-axis. */
    int endpoint_x = 0;
    int field_box_num = 50;

    // internal machinery - don't worry about this
    final Mover mover;

    public BasicTicker(BSim sim, ArrayList<CrossProtectionBacterium> bacA, ArrayList<CrossProtectionBacterium> bacB,
    		ArrayList<BSimCapsuleBacterium> bacteriaAll, int LOG_INTERVAL, Random bacRng, 
    		double growth_stdv, double growth_mean, double length_stdv, double length_mean,
    		ChemicalField antibioticA, ChemicalField antibioticB, int toxin_condition) {
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
        
        BasicTicker.antibioticA = antibioticA;
        BasicTicker.antibioticB = antibioticB;
        this.toxin_condition = toxin_condition;
    }
    
    /** Sets the flag for growth. **/
    public void setGrowth(boolean b) { WITH_GROWTH = b; }
    
    /** Function for bacteria growth activities. */
    public void growBacteria( ArrayList<CrossProtectionBacterium> bac, ArrayList<CrossProtectionBacterium> bacBorn ) {
    	Random bacRng = new Random(); 			// Random number generator
    	bacRng.setSeed(50); 					// Initializes random number generator
    	
        for (CrossProtectionBacterium b : bac) { 				// Loop over bac array
            b.grow();

            // Divide if grown past threshold
            if (b.L >= b.L_th) {
            	CrossProtectionBacterium daughter = b.divide();
            	
            	bacBorn.add(daughter);  		// Add daughter to newborn class, 'mother' keeps her status
            }
        }
        
        bac.addAll(bacBorn); 					// Adds all the newborn daughters to sub-population
        bacteriaAll.addAll(bacBorn); 			// Adds all the newborn daughters to total population
        
        // Allow daughter cells to grow 
        for ( CrossProtectionBacterium b : bacBorn ) {
        	
            // Assigns a division length to each bacterium according to a normal distribution
            double lengthThreshold = length_stdv*bacRng.nextGaussian()+length_mean;
            b.setElongationThreshold(lengthThreshold);
        }
        
        bacBorn.clear(); 						// Cleared for next time-step
    }
    
    /** Function to remove bacteria due to cell death or by boundary. */
    public void removeBacteria( BSim sim, ArrayList<CrossProtectionBacterium> bac, ArrayList<CrossProtectionBacterium> bac_dead ) {

        for (CrossProtectionBacterium b : bac) {
        	
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
        
        /** Allow the flow of antibiotics through a boundary. */
          
        final int VER1 = 1;
        final int VER2 = 2;
        final int VER3 = 3;
        int version = 2;
        
        // State enabled where antibiotic flows in through a boundary
        if ( toxin_condition == FLOW_IN ) {
        	if ( version == VER1 ) {
        		final double conc = 1e3; 
            	if ( endpoint_x < field_box_num ) {
            		for ( int i = 0; i < field_box_num; i ++ ) {
            			antibioticA.addQuantity( endpoint_x, i, 0, conc * sim.getDt() );
            			antibioticB.addQuantity( endpoint_x, i, 0, conc * sim.getDt() );
            		}
            		endpoint_x ++;
            	}
            	else {
            		endpoint_x = 0;
            	}
        	}
        	else if ( version == VER2 ) {
        		final double conc = 50; 
            	if ( endpoint_x < field_box_num ) {
                	for ( int x = 0; x < endpoint_x; x ++ ) {
                		for ( int y = 0; y < field_box_num; y ++ ) {
                			antibioticA.addQuantity( x, y, 0, conc * sim.getDt() );
                			antibioticB.addQuantity( x, y, 0, conc * sim.getDt() );
                		}
                	}
                	endpoint_x ++;
            	}
            	else {
            		endpoint_x = 0;
            	}
        	}
        	else if ( version == VER3 ) {
        		final double conc = 2e3; 
            	for ( int i = 0; i < field_box_num; i ++ ) {
            		antibioticA.addQuantity( 0, i, 0, conc * sim.getDt() );
            		antibioticB.addQuantity( 0, i, 0, conc * sim.getDt() );
            	}
        	}
        }
        
        // Uniformly increase the toxin concentration in the simulation
        if ( toxin_condition == UNIFORM ) {
        	if ( delayCountA == delayA ) {
        		//antibioticA.fill(0, toxin_increment_A * sim.getDt(), toxin_increment_A * sim.getDt());
        		antibioticA.fillArea(toxin_increment_A, initial_conc);
        		delayCountA = 0;
        	}
        	else {
        		delayCountA ++;
        	}
        	
        	if ( delayCountB == delayB ) {
        		//antibioticB.fill(0, toxin_increment_B * sim.getDt(), toxin_increment_B * sim.getDt());
        		antibioticB.fillArea(toxin_increment_B, initial_conc);
        		delayCountB = 0;
        	}
        	else {
        		delayCountB ++;
        	}
			
        }

        // Update the antibiotic field
        antibioticA.update(); 
        antibioticB.update();

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

package BasicSimulation2D;

import bsim.BSim;
import bsim.BSimTicker;
import bsim.capsule.BSimCapsuleBacterium;
import bsim.capsule.Mover;
import bsim.capsule.RelaxationMoverGrid;
import bsim.winter2021.Bacterium;

import java.util.ArrayList;
import java.util.Random;

public class BasicTicker extends BSimTicker {

    BSim sim;
    // logs data about time taken by ticker every LOG_INTERVAL timesteps
    final int LOG_INTERVAL;
    public boolean WITH_GROWTH = true;
    static Random bacRng;
    //growth rate standard deviation
    public final double growth_stdv;
    //growth rate mean
    public final double growth_mean;
    //elongation threshold standard deviation
    public final double length_stdv;
    //elongation threshold mean
    public final double length_mean;

    final ArrayList<Bacterium> bac;
    final ArrayList<BSimCapsuleBacterium> bacteriaAll;
    final ArrayList<Bacterium> bac_born;
    final ArrayList<Bacterium> bac_dead;

    // internal machinery - don't worry about this
    final Mover mover;

    public BasicTicker(BSim sim, ArrayList<Bacterium> bac, ArrayList<BSimCapsuleBacterium> bacteriaAll, int LOG_INTERVAL,
                       Random bacRng, double growth_stdv, double growth_mean, double length_stdv, double length_mean) {
        this.sim = sim;
        this.LOG_INTERVAL = LOG_INTERVAL;
        this.bacRng = bacRng; //random number generator
        this.growth_stdv = growth_stdv;
        this.growth_mean = growth_mean;
        this.length_stdv = length_stdv;
        this.length_mean = length_mean;
        this.bac = bac;
        this.bacteriaAll = bacteriaAll;
        bac_born = new ArrayList();
        bac_dead = new ArrayList();
        mover = new RelaxationMoverGrid(bacteriaAll, sim);
    }

    /** Sets the flag for growth. **/
    public void setGrowth(boolean b) { WITH_GROWTH = b; }

    // 1. update lifetimes
    // 2. Begins an "action" -> this represents one timestep
    // 3. Tells each bacterium to perform their action() function
    // 4. Updates each chemical field in the simulation
    // 5. Bacteria are then told to grow()
    // 6. bacteria which are longer than their threshold are told to divide()
    // 7. forces are applied and bacteria move around
    // 8. bacteria which are out of bounds are removed from the simulation
    @Override
    public void tick() {
        // increase lifetimes of cells
        for (Bacterium b : bac) {
            b.lifetime++;
        }

        // ********************************************** Action
        long startTimeAction = System.nanoTime(); //wall-clock time, for diagnosing timing

        for (BSimCapsuleBacterium b : bacteriaAll) {
            b.action(); //bacteria do action at each time step
            // calculates and stores the midpoint of the cell.
        }

        long endTimeAction = System.nanoTime();
        if ((sim.getTimestep() % LOG_INTERVAL) == 0) {
            System.out.println("Action update for " + bacteriaAll.size() + " bacteria took " + (endTimeAction - startTimeAction) / 1e6 + " ms.");
        }  //outputing how long each step took, once every log interval.

        // ********************************************** Chemical fields
        startTimeAction = System.nanoTime();

        //field activity would go here.

        endTimeAction = System.nanoTime();
        if ((sim.getTimestep() % LOG_INTERVAL) == 0) {
            System.out.println("Chemical field update took " + (endTimeAction - startTimeAction) / 1e6 + " ms.");
        }

        // ********************************************** Growth related activities if enabled.
        if (WITH_GROWTH) {

            // ********************************************** Growth and division
            startTimeAction = System.nanoTime(); //start action timer

            for (Bacterium b : bac) { //loop over bac array
                b.grow();

                // Divide if grown past threshold
                if (b.L >= b.L_th) {
                    bac_born.add(b.divide());  // add daughter to newborn class,  'mother' keeps her status
                }
            }
            bac.addAll(bac_born); //adds all the newborn daughters
            bacteriaAll.addAll(bac_born); //adds all the newborn daughters
            // just added - Sohaib Nadeem
            for (Bacterium b : bac_born) {
                // assigns a growth rate and a division length to each new bacterium according to a normal distribution
                double growthRate = growth_stdv * bacRng.nextGaussian() + growth_mean;
                b.setK_growth(growthRate);

                double lengthThreshold = length_stdv * bacRng.nextGaussian() + length_mean;
                b.setElongationThreshold(lengthThreshold);
            }
            bac_born.clear(); // cleared for next time-step

            endTimeAction = System.nanoTime();
            if ((sim.getTimestep() % LOG_INTERVAL) == 0) {
                System.out.println("Growth and division took " + (endTimeAction - startTimeAction) / 1e6 + " ms.");
            }
            //above: prints out information abt bacteria when u want it to
            // ********************************************** Neighbour interactions
            startTimeAction = System.nanoTime();

            mover.move();

            endTimeAction = System.nanoTime();
            if ((sim.getTimestep() % LOG_INTERVAL) == 0) {
                System.out.println("Wall and neighbour interactions took " + (endTimeAction - startTimeAction) / 1e6 + " ms.");
            }

            // ********************************************** Boundaries/removal
            startTimeAction = System.nanoTime();
            // Removal
            for (Bacterium b : bac) {
//                         Kick out if past the top or bottom boundaries
//                        if ((b.x1.y < 0) && (b.x2.y < 0)) {
//                            act_dead.add(b);
//                        }
//                        if ((b.x1.y > sim.getBound().y) && (b.x2.y > sim.getBound().y)) {
//                            act_dead.add(b);
//                        }
                // kick out if past any boundary
                if (b.position.x < 0 || b.position.x > sim.getBound().x || b.position.y < 0 || b.position.y > sim.getBound().y || b.position.z < 0 || b.position.z > sim.getBound().z) {
                    bac_dead.add(b);
                } //bacteria out of bounds = dead
            }
            bac.removeAll(bac_dead);
            bacteriaAll.removeAll(bac_dead);
            bac_dead.clear();

            endTimeAction = System.nanoTime();
            if ((sim.getTimestep() % LOG_INTERVAL) == 0) {
                System.out.println("Death and removal took " + (endTimeAction - startTimeAction) / 1e6 + " ms.");
            }
        }

        startTimeAction = System.nanoTime();

        endTimeAction = System.nanoTime();
        if ((sim.getTimestep() % LOG_INTERVAL) == 0) {
            // System.out.println("Switch took " + (endTimeAction - startTimeAction) / 1e6 + " ms.");
        }

    }
}

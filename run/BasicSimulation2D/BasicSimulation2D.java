package BasicSimulation2D;

import bsim.BSim;
import bsim.BSimUtils;
import bsim.capsule.BSimCapsuleBacterium;
import bsim.export.BSimPngExporter;
import bsim.export.BSimMovExporter;
import bsim.winter2021.*;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;

import javax.vecmath.Vector3d;
import java.io.*;
import java.util.*;


// Winter 2021 Project
// Growth experiments from Winter 2021
public class BasicSimulation2D {
    // this is the very first function that runs in the simulation
    // this runs the simulation
    public static void main(String[] args) {
        /*
        // starts up JCommander which allows you to read options from the command line more easily
        // command line stuff
         */
        Parameters bsim_parameters = new Parameters();
        new JCommander(bsim_parameters, args);

        // begins simulation
        run(bsim_parameters);
    }

    // a multithreaded version of main; this runs 4 simulation concurrently
    // TODO: output to console and saved files are not being handled correctly
    //  (saved files by the loggers are only kept for a single thread)
    public static void multithreading(String[] args) {
        /*
        ArrayList<Parameters> par = new ArrayList();
        for(int i = 0; i < 4; i++) {
            // starts up JCommander which allows you to read options from the command line more easily
            // command line stuff
            Parameters bsim_parameters = new Parameters();
            new JCommander(bsim_parameters, args);
            par.add(bsim_parameters);
        }

        ArrayList<Thread> threads = new ArrayList();
        for(int i = 0; i < 4; i++) {
            threads.add( new Thread(() -> run(par.get(i))) );
        }
        */

        ArrayList<Thread> threads = new ArrayList();
        Parameters p1 = new Parameters();
        new JCommander(p1, args);
        threads.add( new Thread(() -> run(p1)) );
        Parameters p2 = new Parameters();
        new JCommander(p2, args);
        threads.add( new Thread(() -> run(p2)) );
        Parameters p3 = new Parameters();
        new JCommander(p3, args);
        threads.add( new Thread(() -> run(p3)) );
        Parameters p4 = new Parameters();
        new JCommander(p4, args);
        threads.add( new Thread(() -> run(p4)) );

        for(int i = 0; i < 4; i++) {
            threads.get(i).start();
        }

        for(int i = 0; i < 4; i++) {
            try {
                threads.get(i).join();
            } catch (InterruptedException e) {
                e.printStackTrace(); // if there is an error, this will just print out the message
            }
        }
    }

    public static Bacterium createBacterium(BSim sim, Vector3d x1, Vector3d x2, double growthRate, double lengthThreshold) {
        // creates a new bacterium object whose endpoints correspond to the above data
        //Bacterium bacterium = new Bacterium(sim, x1, x2, origin, -1);
        Bacterium bacterium = new Bacterium(sim, x1, x2);

        // determine the vector between the endpoints
        // if the data suggests that a bacterium is larger than L_max, that bacterium will instead
        // be initialized with length 1. Otherwise, we set the length in the following if statement
        // Earlier they were all initialised to length 1.
        Vector3d dispx1x2 = new Vector3d();
        dispx1x2.sub(x2, x1); // sub is subtract
        double length = dispx1x2.length(); // determined.
        if (length < bacterium.L_max) {
            bacterium.initialise(length, x1, x2); // redundant to record length, but ok.
        }

        // assign the growth rate and division length to the bacterium
        bacterium.setK_growth(growthRate);
        bacterium.setElongationThreshold(lengthThreshold);

        return bacterium;
    }

    // This function does a lot. let's break it down
    // 1. Initializes simulation (how long it runs, etc)
    // 2. Reads initial position data for bacteria and creates new bacterium objects accordingly
    // 3. Runs simulation loop until the simulation time is up
    // |----> Updates bacteria positions according to the forces acting on them
    // |----> Logs all data from the simulation into an excel sheet
    // |----> Saves images of simulation
    public static void run(Parameters parameters) {
        boolean export = parameters.export;
        String export_path = parameters.export_path;
        String input_data = parameters.input_data;
        List<Double> simDimensions = parameters.simDimensions;
        boolean fixedBounds = parameters.fixedBounds;
        int initialPopulation = parameters.initialPopulation;
        double populationRatio = parameters.populationRatio;
        double el_stdv = parameters.el_stdv;
        double el_mean = parameters.el_mean;
        double div_stdv = parameters.div_stdv;
        double div_mean = parameters.div_mean;
        double sim_time = parameters.sim_time;;
        double sim_dt = parameters.sim_dt;
        double export_time = parameters.export_time;
        double k_int = parameters.k_int;
        double k_cell = parameters.k_cell;
        double k_sticking = parameters.k_sticking;
        double range_sticking = parameters.range_sticking;
        double twist = parameters.twist;
        double push = parameters.push;
        double L_asym = parameters.L_asym;
        double asymmetry = parameters.asymmetry;
        double sym_growth = parameters.sym_growth;

        // Assigns the specified forces, range, and impulses
        BSimCapsuleBacterium.setIntForce(k_int);
        BSimCapsuleBacterium.setCellForce(k_cell);
        BSimCapsuleBacterium.setStickForce(k_sticking);
        BSimCapsuleBacterium.setStickingRange(range_sticking);

        Bacterium.setTwist(twist);
        Bacterium.setPush(push);
        Bacterium.setLAsym(L_asym);
        Bacterium.setAsym(asymmetry);
        Bacterium.setSym(sym_growth);

        final double pixel_to_um_ratio = parameters.pixel_to_um_ratio;

        // initializes simulation domain size
        final double simX = simDimensions.get(0);
        final double simY = simDimensions.get(1);
        final double simZ = simDimensions.get(2);

        // saves the exact time when the simulation started, real wall clock time
        long simulationStartTime = System.nanoTime();

        // create the simulation object
        final BSim sim = new BSim();
        sim.setDt(sim_dt);
        sim.setSimulationTime(sim_time);      // specified in time units, could also specify a termination condition elsewhere
        sim.setTimeFormat("0.00");            // Time Format for display on images
        sim.setBound(simX, simY, simZ);       // Simulation domain Boundaries


        // Boundaries periodicity: true means walls, false means periodic
        sim.setSolid(true, true, true);

        // Leaky -> chemicals can escape from sides of six faces of the box, only usable if fixedBounds is false
        // leaky means completely open (only the boundaries in the x and y dimensions)
        // leakyRate is for small molecules, not quite sure what this does
        if (!fixedBounds) {
            sim.setLeaky(true, true, true, true, false, false);
            sim.setLeakyRate(0.1 / 60.0, 0.1 / 60.0, 0.1 / 60.0, 0.1 / 60.0, 0, 0);
        }

        /*********************************************************
         * Create the bacteria
         */
        // Separate lists of bacteria in case we want to manipulate the species individually
        // if multiple subpopulations, they'd be initialized separately, they'd be kept in different
        // need an array for each subpopulation, members can be repeated.
        final ArrayList<Bacterium> bac = new ArrayList();
        // Track all of the bacteria in the simulation, for use of common methods etc
        // A general class, no sub-population specifics
        final ArrayList<BSimCapsuleBacterium> bacteriaAll = new ArrayList();

        Random bacRng = new Random(); //random number generator
        bacRng.setSeed(50); // initializes random number generator

        // gets the location of the file that is currently running
        //specify output file path

        // TODO: FIX ALL OF THE HARDCODED PATHS!!
        String systemPath = new File("").getAbsolutePath();

        String initial_data_path = input_data;
        RawReader reader = new RawReader(pixel_to_um_ratio);

        ArrayList<double[]> cell_endpoints = reader.readcsv(initial_data_path);
        for(int i = 0; i < cell_endpoints.size(); i++) {//double[] cell : cell_endpoints) {
            double[] cell = cell_endpoints.get(i);
            // initializes the endpoints of each bacterium from the array of endpoints; z-dimension is 0.5
            // note that the endpoint positions are scaled by pixel_to_um_ratio,
            Vector3d x1 = new Vector3d(cell[0], cell[1], 0.5);
            Vector3d x2 = new Vector3d(cell[2], cell[3], 0.5);
            // assigns a growth rate and a division length to bacterium according to a normal distribution
            // assigns a growth rate and a division length to bacterium according to a normal distribution
            double growthRate = el_stdv * bacRng.nextGaussian() + el_mean;
            double divThreshold = div_stdv * bacRng.nextGaussian() + div_mean;
            Bacterium bac0 = createBacterium(sim, x1, x2, growthRate, divThreshold);
            // adds the newly created bacterium to our lists for tracking purposes
            bac.add(bac0); // for separate subpopulations
            bacteriaAll.add(bac0);  // for all cells
        }

        /*********************************************************
         * Set up the ticker
         */
        final int LOG_INTERVAL = 100; // logs data every 100 timesteps
        BasicTicker ticker = new BasicTicker(sim, bac, bacteriaAll, LOG_INTERVAL, bacRng, el_stdv, el_mean,
                div_stdv, div_mean);
        ticker.setGrowth(true); // whether to enable growth in the ticker or not; if false nothing can grow
        sim.setTicker(ticker);

        // the rest of the code is the drawer (makes simulation images) and data logger (makes csv files)
        /*********************************************************
         * Set up the drawer
         */
        BasicDrawer drawer = new BasicDrawer(sim, 800, 600, pixel_to_um_ratio, bac);
        sim.setDrawer(drawer);

        /*********************************************************
         * Set up the exporters if in export mode and run the simulation
         */
        if (export) {
            String simParameters = "" + BSimUtils.timeStamp() + "__dim_" + simX + "_" + simY + "_" + simZ
                    + "__ip_" + initialPopulation
                    + "__pr_" + populationRatio;

            if (fixedBounds) {
                simParameters += "__fixedBounds";
            } else {
                simParameters += "__leakyBounds";
            }

            String filePath;
            if(export_path == "default") {
                filePath = BSimUtils.generateDirectoryPath(systemPath +"/" + simParameters + "/");
                //String filePath = BSimUtils.generateDirectoryPath("/home/am6465/tmp-results/" + simParameters + "/");
            } else {
                filePath = BSimUtils.generateDirectoryPath(export_path + "/");
            }


            /** Export a csv file that matches CellProfiler's output */
            CellProfilerLogger cp_logger = new CellProfilerLogger(sim, filePath + "BSim_Simulation.csv", bac, pixel_to_um_ratio);
            // Set export time step, should be the same as sim.dt for the TrackObjects_fields to be correct
            // This is because division events are identified by a lifetime of 0, and lifetime increments are based on on sim.dt
            // (Lifetime could be converted to be in terms of cp_logger.dt if we use a multiple of sim.dt)
            cp_logger.setDt(sim.getDt());
            sim.addExporter(cp_logger);

            /*
            MetaLogger metaLogger = new MetaLogger(sim, filePath + "simInfo.txt", parameters, bac.size());
            metaLogger.setDt(export_time);			// Set export time step
            sim.addExporter(metaLogger);

            PositionLogger posLogger = new PositionLogger(sim, filePath + "position.csv", bacteriaAll);
            posLogger.setDt(export_time);			// set export time step for csv file
            sim.addExporter(posLogger);

            SummaryLogger sumLogger = new SummaryLogger(sim, filePath + "summary.csv", bacteriaAll);
            sumLogger.setDt(export_time);			// Set export time step
            sim.addExporter(sumLogger);
             */

            /** Export a rendered image file */
            BSimPngExporter imageExporter = new BSimPngExporter(sim, drawer, filePath);
            imageExporter.setDt(export_time); //this is how often it will output a frame
            // separate time-resolution for images vs excel file
            sim.addExporter(imageExporter);

            /** Export a video of the simulation */
            BSimMovExporter videoExporter = new BSimMovExporter(sim, drawer, filePath + "video.mp4" );
            videoExporter.setSpeed(1); // the number of simulation time units played in one second
            videoExporter.setDt(0.01); // this is how often (in simulation time) it will output a frame for the video
            sim.addExporter(videoExporter);

            sim.export();

            /**
             * Drawing a java plot once we're done?
             * See TwoCellsSplitGRNTest
             */
        } else {
            sim.preview();
        }

        long simulationEndTime = System.nanoTime();
        System.out.println("Total simulation time: " + (simulationEndTime - simulationStartTime) / 1e9 + " sec.");
    }
}

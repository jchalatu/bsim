package BasicSimulation2D;

import bsim.BSim;
import bsim.BSimTicker;
import bsim.BSimUtils;
import bsim.capsule.BSimCapsuleBacterium;
import bsim.capsule.Mover;
import bsim.capsule.RelaxationMoverGrid;
import bsim.draw.BSimDrawer;
import bsim.draw.BSimP3DDrawer;
import bsim.export.BSimExporter;
import bsim.export.BSimLogger;
import bsim.export.BSimMovExporter;
import bsim.export.BSimPngExporter;
import bsim.winter2021.Bacterium;
import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import processing.core.PConstants;
import processing.core.PGraphics3D;

import javax.vecmath.Vector3d;
import java.awt.*;
import java.io.*;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.*;
import java.util.List;


// Winter 2021 Project
// Growth experiments from Winter 2021

public class BasicSimulation2D {
    /*
    // Sohaib Nadeem
    // pixel to um scaling: the images are a bit more than 2000 pixels wide, while the simulation is rougly 200 micrometers
    // so the conversion factor ends up being 13.89
    static final double pixel_to_um_ratio = 13.89;
    final int width_pixels = 400;
    final int height_pixels = 400;
    final double width_um = width_pixels / pixel_to_um_ratio; // should be kept as constants for cell prof.
    final double height_um = height_pixels / pixel_to_um_ratio; // need to confirm if these values are always used

    // @parameter means an optional user-specified value in the command line
    // export mode means output appears
    @Parameter(names = "-export", description = "Enable export mode.")
    private boolean export = true;

    //set dimnesions in um
    @Parameter(names = "-dim", arity = 3, description = "The dimensions (x, y, z) of simulation environment (um).")
    public List<Double> simDimensions = new ArrayList<>(Arrays.asList(new Double[]{width_um, height_um, 1.}));

    // Boundaries
    //Boolean flag: specifies whether any walls are needed
    @Parameter(names = "-fixedbounds", description = "Enable fixed boundaries. (If not, one boundary will be leaky as real uf chamber).")
    private boolean fixedBounds = true;

    // Grid ->
    // 52x42 -> 546
    // 100x86 -> 2150
    // Random:
    // 50x42 -> 250
    // 100x85 -> 1000
    // Density (cell number)
    //optional call to a default initial set of cells
    @Parameter(names = "-pop", arity = 1, description = "Initial seed population (n_total).")
    public int initialPopulation = 10;

    // A:R ratio
    // for default set of cells, set ratio of two subpopulations
    @Parameter(names = "-ratio", arity = 1, description = "Ratio of initial populations (proportion of activators).")
    public double populationRatio = 0.0;

    //growth rate standard deviation
    @Parameter(names = "-gr_stdv", arity = 1, description = "growth rate standard deviation")
    public static double growth_stdv = 0.2;
    //growth rate mean
    @Parameter(names = "-gr_mean", arity = 1, description = "growth rate mean")
    public static double growth_mean = 2.1;

    //elongation threshold standard deviation
    @Parameter(names = "-len_stdv", arity = 1, description = "elongation threshold standard deviation")
    public static double length_stdv = 0.1;
    //elongation threshold mean
    @Parameter(names = "-len_mean", arity = 1, description = "elongation threshold mean")
    public static double length_mean = 7.0;
    */


    /**
     * Whether to enable growth in the ticker etc. or not...
     */
    // if false nothing can grow
    private static final boolean WITH_GROWTH = true;

    static Random bacRng;

    /*
    public void set_parameters(Parameters par) {
        export = par.export;
        imageDimensions = par.imageDimensions;
        fixedBounds = par.fixedBounds;
        initialPopulation = par.initialPopulation;
        populationRatio = par.export;
        growth_stdv = par.export;
        growth_mean = par.export;
        length_stdv = par.export;
        length_mean = par.export;
        pixel_to_um_ratio = par.export;
        width_pixels = par.export;
        height_pixels = par.export;
        width_um = width_pixels / pixel_to_um_ratio; // should be kept as constants for cell prof.
        height_um = height_pixels / pixel_to_um_ratio; // need to confirm if these values are always used

    }

     */

    /*
    // this is the very first function that runs in the simulation
    // this runs the simulation
    public static void main(String[] args) {

        // creates new simulation data object
        // intializing storage unit for simualation
        BasicSimulation2D bsim_ex = new BasicSimulation2D();

        // starts up JCommander which allows you to read options from the command line more easily
        // command line stuff
        new JCommander(bsim_ex, args);

        // begins simulation
        bsim_ex.run();
    }
     */

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
        int[] imageDimensions = parameters.imageDimensions;
        boolean fixedBounds = parameters.fixedBounds;
        int initialPopulation = parameters.initialPopulation;
        double populationRatio = parameters.populationRatio;
        double growth_stdv = parameters.growth_stdv;
        double growth_mean = parameters.growth_mean;
        double length_stdv = parameters.length_stdv;
        double length_mean = parameters.length_mean;

        final double pixel_to_um_ratio = parameters.pixel_to_um_ratio;
        final int width_pixels = parameters.imageDimensions[0];
        final int height_pixels = parameters.imageDimensions[1];
        final double width_um = width_pixels / pixel_to_um_ratio; // should be kept as constants for cell prof.
        final double height_um = height_pixels / pixel_to_um_ratio; // need to confirm if these values are always used

        // initializes simulation domain size
        final double simX = width_um;
        final double simY = height_um;
        final double simZ = 1.0;

        // saves the exact time when the simulation started, real wall clock time
        long simulationStartTime = System.nanoTime();

        // create the simulation object
        final BSim sim = new BSim();
        sim.setDt(0.05);                    // set Simulation Timestep in time units
        sim.setSimulationTime(5.0);       // specified in time units, could also specify a termination condition elsewhere
        sim.setTimeFormat("0.00");            // Time Format for display on images
        sim.setBound(simX, simY, simZ);        // Simulation domain Boundaries


        // Boundaries periodicity: true means walls, false means periodic
        sim.setSolid(true, true, true);

        // Leaky -> bacteria can escape from sides of six faces of the box. only useable if fixedbounds allows it
        // leaky means completely open
        // leakyrate is for small molecules, not quite sure what this does
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

        bacRng = new Random(); //random number generator
        bacRng.setSeed(50); // initializes random number generator

        // gets the location of the file that is currently running
        //specify output file path
        String systemPath = new File("").getAbsolutePath() + "\\SingleCellSims";

        String initial_data_path = "C:\\Users\\sohai\\IdeaProjects\\bsim\\examples\\PhysModBsim\\twocellssidebyside2-400by400.csv";
        RawReader reader = new RawReader(initial_data_path, pixel_to_um_ratio);
        //String initial_data_path = "C:\\Users\\sohai\\IdeaProjects\\bsim\\examples\\PhysModBsim\\MyExpt_IdentifyPrimaryObjects.csv";
        //CellProfilerReader reader = new CellProfilerReader(initial_data_path, pixel_to_um_ratio, 1);
        ArrayList<double[]> cell_endpoints = reader.readcsv();
        for(int i = 0; i < cell_endpoints.size(); i++) {//double[] cell : cell_endpoints) {
            double[] cell = cell_endpoints.get(i);
            // initializes the endpoints of each bacterium from the array of endpoints; z-dimension is 0.5
            // note that the endpoint positions are scaled by pixel_to_um_ratio,
            Vector3d x1 = new Vector3d(cell[0], cell[1], 0.5);
            Vector3d x2 = new Vector3d(cell[2], cell[3], 0.5);
            // assigns a growth rate and a division length to bacterium according to a normal distribution
            double growthRate = growth_stdv * bacRng.nextGaussian() + growth_mean;
            double lengthThreshold = length_stdv * bacRng.nextGaussian() + length_mean;
            Bacterium bac0 = createBacterium(sim, x1, x2, growthRate, lengthThreshold);
            // adds the newly created bacterium to our lists for tracking purposes
            bac.add(bac0); // for separate subpopulations
            bacteriaAll.add(bac0);  // for all cells
        }

        /*********************************************************
         * Set up the ticker
         */
        final int LOG_INTERVAL = 100; // logs data every 100 timesteps
        BasicTicker ticker = new BasicTicker(sim, bac, bacteriaAll, LOG_INTERVAL, bacRng, growth_stdv, growth_mean,
                length_stdv, length_mean);
        sim.setTicker(ticker);

        // the rest of the code is the drawer (makes simulation images) and data logger (makes csv files)
        /*********************************************************
         * Set up the drawer
         */
        BasicDrawer drawer = new BasicDrawer(sim, width_pixels, height_pixels, pixel_to_um_ratio, bac);
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

            String filePath = BSimUtils.generateDirectoryPath(systemPath + "/" + simParameters + "/");
//            String filePath = BSimUtils.generateDirectoryPath("/home/am6465/tmp-results/" + simParameters + "/");

            double export_time = 0.5; // was 10, and simulation time was 100 - Sohaib Nadeem

            /**
             * Export a rendered image file
             */
            //BSimLogger metaLogger = new BSimLogger(sim, filePath + "simInfo.txt");
            //metaLogger.setDt(export_time);//3600);			// Set export time step
            //sim.addExporter(metaLogger);
            //BSimLogger posLogger = new BSimLogger(sim, filePath + "position.csv");
            //posLogger.setDt(export_time);			// set export time step for csv file
            //sim.addExporter(posLogger);
            //BSimLogger sumLogger = new BSimLogger(sim, filePath + "summary.csv");
            //sumLogger.setDt(export_time);			// Set export time step
            //sim.addExporter(sumLogger);

            BSimPngExporter imageExporter = new BSimPngExporter(sim, drawer, filePath);
            imageExporter.setDt(export_time); //this is how often it will output a frame
            // separate time-resolution for images vs excel file
            sim.addExporter(imageExporter);


            // Export a csv file that matches CellProfiler's output
            CellProfilerLogger cp_logger = new CellProfilerLogger(sim, filePath + "MyExpt_EditedObjects8_simulation.csv", bac, pixel_to_um_ratio);
            // Set export time step, should be the same as sim.dt for the TrackObjects_fields to be correct
            // This is because division events are identified by a lifetime of 0, and lifetime increments are based on on sim.dt
            // (Lifetime could be converted to be in terms of cp_logger.dt if we use a multiple of sim.dt)
            cp_logger.setDt(sim.getDt());
            sim.addExporter(cp_logger);

            //Export a video of the simulation
            BSimMovExporter videoExporter = new BSimMovExporter(sim, drawer, filePath + "video.mp4" );
            videoExporter.setSpeed(1); // the number of simulation time units played in one second
            videoExporter.setDt(0.05); // this is how often (in simulation time) it will output a frame for the video
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

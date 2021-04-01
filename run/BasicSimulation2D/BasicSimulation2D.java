package BasicSimulation2D;

import PhysModBsim.Bacterium;
import PhysModBsim.CellProfilerLogger;

import bsim.BSim;
import bsim.BSimTicker;
import bsim.BSimUtils;
import bsim.capsule.BSimCapsuleBacterium;
import bsim.capsule.Mover;
import bsim.capsule.RelaxationMoverGrid;
import bsim.draw.BSimDrawer;
import bsim.draw.BSimP3DDrawer;
import bsim.export.BSimLogger;
import bsim.export.BSimPngExporter;
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

public class BasicSimulation2D {
    // Sohaib Nadeem
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

    /**
     * Whether to enable growth in the ticker etc. or not...
     */
    // if false nothing can grow
    private static final boolean WITH_GROWTH = true;

    static Random bacRng;


    // this is the very first function that runs in the simulation
    // this runs the simulation
    public static void main(String[] args) {

        // creates new simulation data object
        // intializing storage unit for simualation
        PhysModBsim.NeighbourInteractions bsim_ex = new PhysModBsim.NeighbourInteractions();

        // starts up JCommander which allows you to read options from the command line more easily
        // command line stuff
        new JCommander(bsim_ex, args);

        // begins simulation
        bsim_ex.run();
    }

    public static Bacterium createBacterium(BSim sim, Vector3d x1, Vector3d x2) {
        // creates a new bacterium object whose endpoints correspond to the above data
        Bacterium bacterium = new Bacterium(sim, x1, x2);

        // determine the vector between the endpoints
        // if the data suggests that a bacterium is larger than L_max, that bacterium will instead
        // be initialized with length 1. Otherwise, we set the length in the following if statement
        // Earlier they were all initiatlied to lengthe 1.
        Vector3d dispx1x2 = new Vector3d();
        dispx1x2.sub(x2, x1); // sub is subtract
        double length = dispx1x2.length(); // determined.
        if (length < bacterium.L_max) {
            bacterium.initialise(length, x1, x2); // redundant to record length, but ok.
        }

        // assigns a growth rate and a division length to bacterium according to a normal distribution
        double growthRate = growth_stdv * bacRng.nextGaussian() + growth_mean;
        bacterium.setK_growth(growthRate);
        double lengthThreshold = length_stdv * bacRng.nextGaussian() + length_mean;
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
    public void run() {
        // initializes simulation domain size
        final double simX = simDimensions.get(0);
        final double simY = simDimensions.get(1);
        final double simZ = simDimensions.get(2);

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
        // if multiple subpopulations, they'd be initiallized separately, they'd be kept in different
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

        CellProfilerReader reader = new CellProfilerReader("C:\\Users\\sohai\\IdeaProjects\\bsim\\examples\\PhysModBsim\\twocellssidebyside2-400by400.csv",
                pixel_to_um_ratio, 1);
        double [] cell_endpoints = reader.readcsv();
        while(cell_endpoints != null) {
            // initializes the endpoints of each bacterium from the array of endpoints
            // z-dimension is a small number, randomly generated, not sure why.
            Vector3d x1 = new Vector3d(cell_endpoints[0], cell_endpoints[1], 0.5);
            Vector3d x2 = new Vector3d(cell_endpoints[2], cell_endpoints[3], 0.5);
            // note that the endpoint positions are scaled by 13.89, since the images are a bit more than 2000 pixels wide
            // while the simulation is rougly 200 micrometers. the conversion factor ends up being 13.89
            //pixel to um scaling

            Bacterium bac0 = createBacterium(sim, x1, x2);

            // adds the newly created bacterium to our lists for tracking purposes
            bac.add(bac0); //for separate subpopulations
            bacteriaAll.add(bac0);  // for all cells

            cell_endpoints = reader.readcsv();
        }

        /*********************************************************
         * Set up the ticker
         */
        final int LOG_INTERVAL = 100; // logs data every 100 timesteps
        BasicTicker ticker = new BasicTicker(sim, LOG_INTERVAL, bacRng);
        sim.setTicker(ticker);

        // the rest of the code is the drawer (makes simulation images) and data logger (makes csv files)
        /*********************************************************
         * Set up the drawer
         */
        BasicDrawer drawer = new BasicDrawer(sim, width_pixels, height_pixels);
        sim.setDrawer(drawer);

        export = false;
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
            BSimPngExporter imageExporter = new BSimPngExporter(sim, drawer, filePath);
            imageExporter.setDt(export_time); //this is how often it will output a frame
            // separate time-resolution for images vs excel file
            sim.addExporter(imageExporter);


            // Export a csv file that matches CellProfiler's output
            CellProfilerLogger cp_logger = new CellProfilerLogger(sim, filePath + "MyExpt_EditedObjects8_simulation.csv", bacteriaAll);
            cp_logger.setDt(export_time);            // Set export time step
            sim.addExporter(cp_logger);

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

// Fall 2020 Project
// Growth experiments form Fall 2020

package BSimPhageField;

import bsim.winter2021.CellProfilerLogger;
import bsim.BSim;

import bsim.BSimChemicalField;
import bsim.BSimUtils;
import bsim.capsule.BSimCapsuleBacterium;
import bsim.capsule.Mover;
import bsim.capsule.RelaxationMoverGrid;
import bsim.export.BSimLogger;
import bsim.export.BSimPngExporter;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;

import javax.vecmath.*;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.*;
import java.util.List;
import java.io.File;

/**
 *
 * This class represents an environment with phage.
 * A number of BSimBacterium bacteria are set up to swim around in a phage field.
 * If the density of phage around a bacteria is above a certain threshold, the bacteria
 * will become infected.
 * Infected bacteria will produce and release phages into the surrounding medium.
 *
 */
public class BSimPhageField {

    static final double pixel_to_um_ratio = 13.89;
    final int width_pixels = 800;
    final int height_pixels = 600;
    final double width_um = width_pixels / pixel_to_um_ratio; 	// should be kept as constants for cell prof.
    final double height_um = height_pixels / pixel_to_um_ratio; // need to confirm if these values are always used

    // Boundaries
    // Boolean flag: specifies whether any walls are needed
    @Parameter(names = "-fixedbounds", description = "Enable fixed boundaries. (If not, one boundary will be leaky as real uf chamber).")
    private boolean fixedBounds = false;

    // @parameter means an optional user-specified value in the command line
    // export mode means output appears
    @Parameter(names = "-export", description = "Enable export mode.")
    private boolean export = true;

    @Parameter(names = "-export_path", description = "export location")
    private String export_path = "default";

    @Parameter(names = "-grow", description = "Enable bacteria growth.")
    private boolean WITH_GROWTH = true;
    @Parameter(names = "-flow_in", description = "Enable phage from boundary.")
    private boolean FLOW_IN = false;
    @Parameter(names = "-phage", description = "Enable initial phage distribution.")
    private boolean ADD_PHAGE = false;
    @Parameter(names = "-phage_quant", description = "Initial quantity of phage distribution.")
    private double initial_quantity = 1e5;
    @Parameter(names = "-phage_pos", description = "Initial position of phage distribution.")
    private Vector3d initial_pos = new Vector3d(width_um/2, height_um/2, 1);

    // Simulation setup parameters. Set dimensions in um
    @Parameter(names = "-dim", arity = 3, description = "The dimensions (x, y, z) of simulation environment (um).")
    //public List<Double> simDimensions = new ArrayList<>(Arrays.asList(new Double[] {75.0, 50.0, 1.0} ));
    public List<Double> simDimensions = new ArrayList<>(Arrays.asList(new Double[]{width_um, height_um, 1.}));

    // Grid ->
    // 52x42 -> 546
    // 100x86 -> 2150
    // Random:
    // 50x42 -> 250
    // 100x85 -> 1000
    // Density (cell number)
    // optional call to a default initial set of cells
    @Parameter(names = "-pop", arity = 1, description = "Initial seed population (n_total).")
    public int initialPopulation = 1;

    // A:R ratio
    // for default set of cells, set ratio of two subpopulations
    @Parameter(names = "-ratio", arity = 1, description = "Ratio of initial populations (proportion of activators).")
    public double populationRatio = 0.0;

    //growth rate standard deviation
    @Parameter(names="-el_stdv",arity=1,description = "elongation rate standard deviation")
    public static double el_stdv = 0.277;
    @Parameter(names="-el_mean",arity=1,description = "elongation rate mean")
    public static double el_mean = 1.23;

    //elongation threshold standard deviation
    @Parameter(names="-div_stdv",arity=1,description = "elongation threshold standard deviation")
    public static double div_stdv = 0.1;
    //elongation threshold mean
    @Parameter(names="-div_mean",arity=1,description = "elongation threshold mean")
    public static double div_mean = 7.0;

    // Simulation Time
    @Parameter(names="-simt",arity=1,description = "simulation time")
    public static double sim_time = 6.5;
    @Parameter(names="-simdt",arity=1,description = "simulation time step")
    public static double sim_dt = 0.05;
    @Parameter(names="-export_time",arity=1,description = "export time")
    public static double export_time = 0.5;// Previously was 10, and simulation time was 100

    // internal force
    @Parameter(names="-k_int",arity=1,description = "internal force")
    public static double k_int = 50.0;
    // cell-cell collision force
    @Parameter(names="-k_cell",arity=1,description = "cell-cell collision force")
    public static double k_cell = 500.0;
    // sticking force
    @Parameter(names="-k_stick",arity=1,description = "side-to-side attraction")
    public static double k_sticking = 0.01;

    // sticking range
    @Parameter(names="-rng_stick",arity=1,description = "max range side-to-side attraction")
    public static double range_sticking = 5.0;

    // twist
    @Parameter(names="-twist",arity=1,description = "twist")
    public static double twist = 0.1;
    // push
    @Parameter(names="-push",arity=1,description = "push")
    public static double push = 0.05;

    // asymmetric growth threshold
    @Parameter(names="-l_asym",arity=1,description = "asymmetric growth threshold")
    public static double L_asym = 3.75;
    // value of asymmetry
    @Parameter(names="-asym",arity=1,description = "asymmetry")
    public static double asymmetry = 0.1;
    // symmetric growth
    @Parameter(names="-sym",arity=1,description = "symmetric growth")
    public static double sym_growth = 0.05;

    /** Main Function.
     * This is the very first function that runs in the simulation.
     * This runs the simulation. */
    public static void main(String[] args) {

        // Creates new simulation data object
        // Initializing storage unit for simulation
    	BSimPhageField bsim_ex = new BSimPhageField();

        // Starts up JCommander which allows you to read options from the command line more easily
        // Command line stuff
        new JCommander(bsim_ex, args);

        // Begins the simulation
        bsim_ex.run();
    }

    /** Creates a new Bacterium object. */
    public static PhageFieldBacterium createBacterium(BSim sim, BSimChemicalField field ) {
    	Random bacRng = new Random(); 		// Random number generator
        bacRng.setSeed(50); 				// Initializes random number generator

        // Random initial positions
        Vector3d pos1 = new Vector3d(Math.random()*sim.getBound().x,
				Math.random()*sim.getBound().y,
				Math.random()*sim.getBound().z);

        double r = div_mean * Math.sqrt(Math.random());
        double theta = Math.random() * 2 * Math.PI;
        Vector3d pos2 = new Vector3d(pos1.x + r * Math.cos(theta),
				pos1.y + r * Math.sin(theta),
				pos1.z);

        // Check if the random coordinates are within bounds
        while(pos2.x >= sim.getBound().x || pos2.x <= 0 || pos2.y >= sim.getBound().y || pos2.y <= 0) {
        	pos1 = new Vector3d(Math.random()*sim.getBound().x,
    				Math.random()*sim.getBound().y,
    				Math.random()*sim.getBound().z);
        	pos2 = new Vector3d(pos1.x + r * Math.cos(theta),
    				pos1.y + r * Math.sin(theta),
    				pos1.z);
        }

        // Creates a new bacterium object whose endpoints correspond to the above data
        PhageFieldBacterium bacterium = new PhageFieldBacterium(sim, field, pos1, pos2);

        // Determine the vector between the endpoints
        // If the data suggests that a bacterium is larger than L_max, that bacterium will instead
        // be initialized with length 1. Otherwise, we set the length in the following if statement
        // Earlier they were all initiatlied to length 1.
        Vector3d dispx1x2 = new Vector3d();
        dispx1x2.sub(pos2, pos1); 						// Sub is subtract
        double length = dispx1x2.length(); 				// Determined.

        if (length < bacterium.L_max) {
            bacterium.initialise(length, pos1, pos2); 	// Redundant to record length, but ok.
        }

        // Assigns a growth rate and a division length to bacterium according to a normal distribution
        double growthRate = el_stdv * bacRng.nextGaussian() + el_mean;
        bacterium.setK_growth(growthRate);

        double lengthThreshold = div_stdv * bacRng.nextGaussian() + div_mean;
        bacterium.setElongationThreshold(lengthThreshold);

        // Assigns the specified forces, range, and impulses
        bacterium.setIntForce(k_int);
        bacterium.setCellForce(k_cell);
        bacterium.setStickForce(k_sticking);
        bacterium.setStickingRange(range_sticking);
        bacterium.setTwist(twist);
        bacterium.setPush(push);

        bacterium.setAsym(asymmetry);

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

		/*********************************************************
		 * Define simulation domain size and start time
		 */
        final double simX = simDimensions.get(0);
        final double simY = simDimensions.get(1);
        final double simZ = simDimensions.get(2);

        // Saves the exact time when the simulation started, real wall clock time
        long simulationStartTime = System.nanoTime();

		/*********************************************************
		 * Create a new simulation object and set up simulation settings
		 */
        final BSim sim = new BSim();
        sim.setDt(sim_dt);					// Set simulation timestep in time units
        									// Let the time units be in hours
        sim.setSimulationTime(sim_time);    // Specified in time units, could also specify a termination condition elsewhere
        sim.setTimeFormat("0.00");		    // Time Format for display on images
        sim.setBound(simX, simY, simZ);		// Simulation domain Boundaries

		// NOTE - solid = false sets a periodic boundary condition.
		// This overrides leakiness!
		sim.setSolid(true, true, true);		// Solid (true) or wrapping (false) boundaries

        // Leaky -> bacteria can escape from sides of six faces of the box.
		// Only usable if fixedbounds allows it

		// Set a closed or open boundary for phage diffusion
		if ( fixedBounds ) {
			sim.setLeaky(false, false, false, false, false, false);
			sim.setLeakyRate(0, 0, 0, 0, 0, 0);
		}
		else {
			sim.setLeaky(true, true, true, true, false, false);
			sim.setLeakyRate(1, 1, 1, 1, 0, 0);
		}

		/*********************************************************
		 * Set up the phage field
		 */
		final double c = 1e3;        		// Molecules; Decrease this to see lower concentration
		final double decayRate = 0.00308;	// Decay rate of phages (M13: 0.074/day), 0.5
											// 0.074/24hrs = 0.00308/hr
		final double diffusivity = 3;	 	// (Microns)^2/sec

		final int field_box_num = 50;		// Number of boxes to represent the chemical field
		final BSimChemicalField field = new BSimChemicalField(sim, new int[]{field_box_num, field_box_num, 1}, diffusivity, decayRate);

        /*********************************************************
         * Create the bacteria
         */

        // Separate lists of bacteria in case we want to manipulate the species individually
        // If multiple subpopulations, they'd be initialized separately, they'd be kept in different
        // need an array for each subpopulation, members can be repeated.
        final ArrayList<PhageFieldBacterium> bac = new ArrayList();

        /** Track all of the bacteria in the simulation, for use of common methods etc.
        A general class, no sub-population specifics. */
        final ArrayList<BSimCapsuleBacterium> bacteriaAll = new ArrayList();

        Random bacRng = new Random(); 		// Random number generator
        bacRng.setSeed(50); 				// Initializes random number generator

        // Gets the location of the file that is currently running
        // Specify output file path
        String systemPath = new File("").getAbsolutePath()+"\\PhageFieldSims";

		// Create new phage sensing bacteria objects randomly in space
		while( bacteriaAll.size() < initialPopulation ) {

            // Creates a new bacterium object whose endpoints correspond to the above data
			PhageFieldBacterium bacterium = createBacterium( sim, field );

			// Adds the newly created bacterium to our lists for tracking purposes
			bac.add(bacterium); 			// For separate subpopulations
			bacteriaAll.add(bacterium);		// For all cells
		}

        // Set up stuff for growth. Placeholders for the recently born and dead
        final ArrayList<PhageFieldBacterium> bac_born = new ArrayList();
        final ArrayList<PhageFieldBacterium> bac_dead = new ArrayList();

        // Internal machinery - dont worry about this line
        // Some kind of initialize of mover
        final Mover mover;
        mover = new RelaxationMoverGrid(bacteriaAll, sim);

        /*********************************************************
         * Add initial distribution of phage
         */
        if (ADD_PHAGE) {
        	field.addQuantity(initial_pos, initial_quantity);
        }

        /*********************************************************
         * Set up the ticker
         */
        final int LOG_INTERVAL = 100; 			// logs data every 100 timesteps
        PhageFieldTicker ticker = new PhageFieldTicker(sim, bac, bacteriaAll, LOG_INTERVAL, bacRng, el_stdv, el_mean,
                div_stdv, div_mean, field, field_box_num);
        ticker.setGrowth(WITH_GROWTH);			// enables bacteria growth
        ticker.setFlow(FLOW_IN);				// enables flow in of phage from boundary
        sim.setTicker(ticker);

        /*********************************************************
         * Set up the drawer
         */
        PhageFieldDrawer drawer = new PhageFieldDrawer(sim, width_pixels, height_pixels, pixel_to_um_ratio, bac,
        		field, c);
        sim.setDrawer(drawer);

        /*********************************************************
         * Export data
         */
        if(export) {
            String simParameters = "" + BSimUtils.timeStamp() + "__dim_" + simX + "_" + simY + "_" + simZ
                    + "__ip_" + initialPopulation
                    + "__pr_" + populationRatio;

            if(fixedBounds){
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

            /*********************************************************
             * Various properties of the simulation, for future reference.
             */
            BSimLogger metaLogger = new BSimLogger(sim, filePath + "simInfo.txt") {
                @Override
                public void before() {
                    super.before();
                    write("Simulation metadata.");
                    write("Catie Terrey Fall 2020."); //change name when new person :)
                    write("Simulation dimensions: (" + simX + ", " + simY + ", " + simZ + ")");
                    write("Initial population: "+ initialPopulation);
                    write("Ratio " + populationRatio);

                    if(fixedBounds){
                        write("Boundaries: fixed");
                    } else {
                        write("Boundaries: leaky");
                    }
                }

                @Override
                public void during() {

                }
            };
            metaLogger.setDt(export_time);			// Set export time step
            sim.addExporter(metaLogger);

            BSimLogger posLogger = new BSimLogger(sim, filePath + "position.csv") {
                DecimalFormat formatter = new DecimalFormat("###.##", DecimalFormatSymbols.getInstance( Locale.ENGLISH ));

                @Override
                public void before() {
                    super.before();
                    write("per Act; per Rep; id, p1x, p1y, p1z, p2x, p2y, p2z, growth_rate,directions");
                }

                @Override
                public void during() {
                    String buffer = new String();

                    buffer += sim.getFormattedTime() + "\n";
                    write(buffer);

                    write("acts");

                    buffer = "";
                    for(BSimCapsuleBacterium b : bacteriaAll) {
                        buffer += b.id + "," + formatter.format(b.x1.x)
                                + "," + formatter.format(b.x1.y)
                                + "," + formatter.format(b.x1.z)
                                + "," + formatter.format(b.x2.x)
                                + "," + formatter.format(b.x2.y)
                                + "," + formatter.format(b.x2.z)
                                + "," + formatter.format(b.getK_growth())
                                + "\n";
                    }

                    write(buffer);

                }
            };
            posLogger.setDt(export_time);			// set export time step for csv file
            sim.addExporter(posLogger);

            BSimLogger sumLogger = new BSimLogger(sim, filePath + "summary.csv") {

                @Override
                public void before() {
                    super.before();
                    write("time,id, p1x, p1y, p1z, p2x, p2y, p2z, px, py, pz, growth_rate, directions");
                }

                @Override
                public void during() {
                    String buffer = new String();
                    buffer = "";
                    for(BSimCapsuleBacterium b : bacteriaAll) {
                        buffer += sim.getFormattedTime()+","+b.id
                                + "," + b.x1.x
                                + "," + b.x1.y
                                + "," + b.x1.z
                                + "," + b.x2.x
                                + "," + b.x2.y
                                + "," + b.x2.z
                                + "," + b.position.x
                                + "," + b.position.y
                                + "," + b.position.z
                                + "," + b.getK_growth()
                                + "," + b.direction()
                                + "\n";
                    }

                    write(buffer);
                }
            };
            sumLogger.setDt(export_time);			// Set export time step
            sim.addExporter(sumLogger);

            /**
             * Export a rendered image file
             */
            BSimPngExporter imageExporter = new BSimPngExporter(sim, drawer, filePath );
            imageExporter.setDt(export_time); //this is how often it will output a frame
            // separate time-resolution for images vs excel file
            sim.addExporter(imageExporter);

            /**
             * Export a csv file to save information about infection
             */
            //PhageFieldLogger logger = new PhageFieldLogger(sim, filePath + "BSim_Simulation.csv", bac, pixel_to_um_ratio);
            //logger.setDt(sim.getDt());			// Set export time step
            //sim.addExporter(logger);

            /** Export a csv file that matches CellProfiler's output */
            CellProfilerLogger cp_logger = new CellProfilerLogger(sim, filePath + "BSim_Simulation.csv", bac, pixel_to_um_ratio);
            // Set export time step, should be the same as sim.dt for the TrackObjects_fields to be correct
            // This is because division events are identified by a lifetime of 0, and lifetime increments are based on on sim.dt
            // (Lifetime could be converted to be in terms of cp_logger.dt if we use a multiple of sim.dt)
            cp_logger.setDt(sim.getDt());
            sim.addExporter(cp_logger);

            sim.export();

            /**
             * Drawing a java plot once we're done?
             * See TwoCellsSplitGRNTest
             */

        }

        // Run the simulation
        else {
            sim.preview();
        }

        long simulationEndTime = System.nanoTime();

        System.out.println("Total simulation time: " + (simulationEndTime - simulationStartTime)/1e9 + " sec.");

    }


}

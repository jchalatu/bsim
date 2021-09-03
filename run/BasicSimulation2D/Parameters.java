package BasicSimulation2D;

import com.beust.jcommander.Parameter;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class Parameters {
    // Simulation Time
    @Parameter(names="-simt",arity=1,description = "simulation time")
    public static double sim_time = 2.0;
    @Parameter(names="-simdt",arity=1,description = "simulation time step")
    public static double sim_dt = 0.01;

    // @parameter means an optional user-specified value in the command line
    // export mode means output appears
    @Parameter(names = "-export", description = "Enable export mode.")
    public boolean export = false;

    @Parameter(names = "-export_path", description = "export location")
    public String export_path = "default";

    @Parameter(names = "-input_data", description = "path to input data")
    public String input_data = "default";

    @Parameter(names="-export_time",arity=1,description = "export time")
    public static double export_time = 0.5;// Previously was 10, and simulation time was 100

    //set dimensions in um
    @Parameter(names = "-dim", arity = 3, description = "The dimensions (x, y, z) of simulation environment (um).")
    //public List<Double> simDimensions = new ArrayList<>(Arrays.asList(new Double[]{198.0, 159.0, 1.0}));
    //public List<Double> simDimensions = new ArrayList<>(Arrays.asList(new Double[]{1870 / 13.89, 2208 / 13.89, 1.0}));
    public List<Double> simDimensions = new ArrayList<>(Arrays.asList(new Double[]{1800 / 13.89, 1800 / 13.89, 1.0}));

    // pixel to um scaling: the images are a bit more than 2000 pixels wide, while the simulation is rougly 200 micrometers
    // so the conversion factor ends up being 13.89
    @Parameter(names = "-pixel_to_um_ratio", description = "The image scaling in pixels to um (microns).")
    public double pixel_to_um_ratio = 13.89;

    // Boundaries
    //Boolean flag: specifies whether any walls are needed
    @Parameter(names = "-fixedbounds", description = "Enable fixed boundaries. (If not, one boundary will be leaky as real uf chamber).")
    public boolean fixedBounds = true;

    // Grid ->
    // 52x42 -> 546
    // 100x86 -> 2150
    // Random:
    // 50x42 -> 250
    // 100x85 -> 1000
    // Density (cell number)
    //optional call to a default initial set of cells
    @Parameter(names = "-pop", arity = 1, description = "Initial seed population (n_total).")
    public int initialPopulation = 1;

    // A:R ratio
    // for default set of cells, set ratio of two subpopulations
    @Parameter(names = "-ratio", arity = 1, description = "Ratio of initial populations (proportion of activators).")
    public double populationRatio = 0.0;


    //growth rate standard deviation
    @Parameter(names="-el_stdv",arity=1,description = "elongation rate standard deviation")
    public static double el_stdv = 0.2;//0.277;
    @Parameter(names="-el_mean",arity=1,description = "elongation rate mean")
    public static double el_mean = 2.0;//1.23;

    //elongation threshold standard deviation
    @Parameter(names="-div_stdv",arity=1,description = "elongation threshold standard deviation")
    public static double div_stdv = 0.1;
    //elongation threshold mean
    @Parameter(names="-div_mean",arity=1,description = "elongation threshold mean")
    public static double div_mean = 7.0;

    // internal force
    @Parameter(names="-k_int",arity=1,description = "internal force")
    public static double k_int = 50.0;
    // cell-cell collision force
    @Parameter(names="-k_overlap",arity=1,description = "cell-cell collision force")
    public static double k_overlap = 500.0;
    // sticking force
    @Parameter(names="-k_stick",arity=1,description = "side-to-side attraction")
    public static double k_sticking = 0.01;


    // contact threshold (for sticking force)
    @Parameter(names="-contact_damping",arity=1,description = "contact threshold for sticking interaction")
    public static double contact_damping = 2.0;
    // contact range extension (pilus length)
    @Parameter(names="-contact_rng",arity=1,description = "contact range for pilus interaction")
    public static double contact_range = 0.25;

    // twist
    @Parameter(names="-birth_twist",arity=1, description = "twist")
    public static double birth_twist = 0.25;
    // push
    @Parameter(names="-push",arity=1,description = "push")
    public static double push = 0.05;
    // initial degree of asymmetrical growth
    @Parameter(names="-init_growth_asym",arity=1,description = "initital asymmetry")
    public static double init_growth_asym = 0.5;
    // length scale at which asymmetrical growth stops
    @Parameter(names="-asymmetry_scale",arity=1,description = "asymmetry scale")
    public static double asymmetry_scale = 0.75;
}

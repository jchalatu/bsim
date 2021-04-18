package BasicSimulation2D;

import com.beust.jcommander.Parameter;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class Parameters {
    // @parameter means an optional user-specified value in the command line
    // export mode means output appears
    @Parameter(names = "-export", description = "Enable export mode.")
    public boolean export = false;

    // image dimension in pixels
    @Parameter(names = "-dim", arity = 2, description = "The dimensions (x, y) of image in pixels.")
    //public List<Integer> imageDimensions = new ArrayList<>(Arrays.asList(new Integer[]{400, 400}));
    public int[] imageDimensions = {800, 800};

    /*
    //set dimensions in um
    @Parameter(names = "-dim", arity = 3, description = "The dimensions (x, y, z) of simulation environment (um).")
    public List<Double> simDimensions = new ArrayList<>(Arrays.asList(new Double[]{198.0, 159.0, 1.}));
     */

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
}

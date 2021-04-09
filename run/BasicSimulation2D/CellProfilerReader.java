package BasicSimulation2D;

import PhysModBsim.Bacterium;

import javax.vecmath.Vector3d;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;

public class CellProfilerReader {
    BufferedReader csvReader;
    final double pixel_to_um_ratio;
    int image_number = 1;

    // TODO: add exception when column names don't match
    public CellProfilerReader(String filepath, double pixel_to_um_ratio, int image_number) {
        this.pixel_to_um_ratio = pixel_to_um_ratio;
        this.image_number = image_number;

        try {
            // try reading the initial position file
            csvReader = new BufferedReader(new FileReader(filepath));
        } catch (FileNotFoundException e) {
            // if that doesn't work, print out an error
            e.printStackTrace();
        }

        // check if fields match those of the required csv file from CellProfiler
        try {
            String row = csvReader.readLine();
            String columns_names[] = row.split(",");
            if(columns_names[0] != "ImageNumber" || columns_names[8] != "AreaShape_Center_X" ||
                    columns_names[9] != "AreaShape_Center_Y" || columns_names[16] != "AreaShape_MajorAxisLength" ||
                    columns_names[22] != "AreaShape_MinorAxisLength") {
                // throw an exception
            }
        } catch(IOException e) {
            e.printStackTrace(); // if there is an error, this will just print out the message
        }
    }

    public ArrayList<double[]> readcsv() {
        ArrayList<double[]> initEndpoints = new ArrayList();

        try {
            String row;
            while ( (row = csvReader.readLine()) != null ) {
                double cell_info[] = Arrays.stream(row.split(",")).mapToDouble(Double::parseDouble).toArray();

                if ((int) cell_info[0] == image_number) {
                    double cell_length = (cell_info[16] - cell_info[22]) / pixel_to_um_ratio;
                    double cell_orientation = -1 * (cell_info[23] + (Math.PI / 2));
                    double cell_center_x = cell_info[8] / pixel_to_um_ratio;
                    double cell_center_y = cell_info[9] / pixel_to_um_ratio;

                    double axis_x = cell_length * Math.cos(cell_orientation); // use trig. identity to simplify?
                    double axis_y = cell_length * Math.sin(cell_orientation); // use trig. identity to simplify?

                    //Vector3d x1 = new Vector3d(cell_center_x + (axis_x / 2), cell_center_y + (axis_y / 2), 0.5/*bacRng.nextDouble()*0.1*(simZ - 0.1)/2.0*/);
                    //Vector3d x2 = new Vector3d(cell_center_x - (axis_x / 2), cell_center_y - (axis_y / 2), 0.5/*bacRng.nextDouble()*0.1*(simZ - 0.1)/2.0*/);

                    initEndpoints.add(new double[] {cell_center_x + (axis_x / 2), cell_center_y + (axis_y / 2),
                            cell_center_x - (axis_x / 2), cell_center_y - (axis_y / 2)});
                }
            }
        } catch (IOException e) {
            e.printStackTrace(); // if there is an error, this will just print out the message
        }

        return initEndpoints;
    }
}

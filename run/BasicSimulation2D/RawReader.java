package BasicSimulation2D;

import PhysModBsim.Bacterium;

import javax.vecmath.Vector3d;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;

public class RawReader {
    BufferedReader csvReader;
    double[][] initEndpoints = new double[4][];
    int rows;
    int row_number;
    final double pixel_to_um_ratio;

    public RawReader(String filepath, double pixel_to_um_ratio) {
        this.pixel_to_um_ratio = pixel_to_um_ratio;

        try {
            // try reading the initial position file
            csvReader = new BufferedReader(new FileReader(filepath));
        } catch (FileNotFoundException e) {
            // if that doesn't work, print out an error
            e.printStackTrace();
        }

        try {
            // goes through each row of the excel sheet and pulls out the initial positions
            String row = csvReader.readLine();
            int i=0;
            while (row != null) {
                // row.split takes a single line of the excel sheet and chops it up into the columns
                // maptodouble then takes the values in those columns and converts them to Java double data format
                // toarray then finally converts the data into an array
                initEndpoints[i] = Arrays.stream(row.split(",")).mapToDouble(Double::parseDouble).toArray();

                row = csvReader.readLine();
                i++;
            }
            csvReader.close(); //finally, close the file once all data is extracted
        } catch(IOException e) {
            e.printStackTrace(); // if there is an error, this will just print out the message
        }

        rows = initEndpoints[0].length;
        row_number = 0;
    }

    // TODO: improve
    public double[] readcsv() {
        if (row_number < rows) {
            double[] endpoints = {initEndpoints[0][row_number] / pixel_to_um_ratio,
                    initEndpoints[1][row_number] / pixel_to_um_ratio,
                    initEndpoints[2][row_number] / pixel_to_um_ratio,
                    initEndpoints[3][row_number] / pixel_to_um_ratio};
            row_number++;
            return endpoints;
        }
        else {
            return null;
        }
    }
}

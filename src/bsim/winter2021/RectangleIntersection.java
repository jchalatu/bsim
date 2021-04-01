package bsim.winter2021;

import javax.vecmath.Vector2d;
import java.util.ArrayList;
import java.util.Arrays;

public class RectangleIntersection {
    public static Vector2d line_intersection(double a1, double b1, double c1, double a2, double b2, double c2) {
        // See e.g. https://en.wikipedia.org/wiki/Line%E2%80%93line_intersection#Using_homogeneous_coordinates
        double w = a1 * b2 - b1 * a2;
        return new Vector2d((b1 * c2 - c1 * b2) / w, (c1 * a2 - a1 * c2) / w);
    }

    public static double distance(Vector2d v1, Vector2d v2) {
        return Math.sqrt(Math.pow((v2.x - v1.x), 2) + Math.pow((v2.y - v1.y), 2));
    }

    public static double rectangle_intersection_area(Vector2d[] r_1, Vector2d[] r_2) {
        // Use the vertices of the first rectangle as the starting vertices of the intersection polygon.
        ArrayList<Vector2d> intersection = new ArrayList<>(Arrays.asList(r_1));

        // Loop over the edges of the second rectangle
        // Determine the intersection of the two rectangle r_1 and r_2
        for (int i = 0; i < r_2.length; i++) {
            if (intersection.size() <= 2) {
                break; // No intersection
            }

            // define a line with the points r_2[i] and r_2[(i + 1) % r_2.length] (adjacent vertices of r_2)
            // a, b, c are constants that determine the equation of the line: ax + by + c = 0
            Vector2d v1 = r_2[i];
            Vector2d v2 = r_2[(i + 1) % r_2.length];
            double a = v2.y - v1.y;
            double b = v1.x - v2.x;
            double c = v2.x * v1.y - v2.y * v1.x;
            // Any point p with line(p) <= 0 is on the "inside" (or on the boundary),
            // any point p with line(p) > 0 is on the "outside".

            // Loop over the edges of the intersection polygon, and determine which part
            // is inside and which is outside
            ArrayList<Vector2d> new_intersection = new ArrayList<Vector2d>();
            double[] line_values = new double[intersection.size()];
            for (int j = 0; j < intersection.size(); j++) {
                line_values[j] = a * intersection.get(j).x + b * intersection.get(j).y + c;
            }
            for(int j = 0; j < intersection.size(); j++) {
                double s_value = line_values[j];
                double t_value = line_values[(j + 1) % line_values.length];

                if(s_value <= 0) {
                    // TODO: change this to add a newly created vector that is a copy of intersection.get(j) ?
                    new_intersection.add(intersection.get(j));
                }
                if (s_value * t_value < 0) {
                    // Points are on opposite sides. Add the intersection of the lines to new_intersection.
                    Vector2d s = intersection.get(j);
                    Vector2d t = intersection.get((j + 1) % intersection.size());
                    double a2 = t.y - s.y;
                    double b2 = s.x - t.x;
                    double c2 = t.x * s.y - t.y * s.x;
                    Vector2d intersection_point = line_intersection(a, b, c, a2, b2, c2);
                    new_intersection.add(intersection_point);
                }
            }

            // set intersection to new intersection
            intersection = new_intersection;
        }

        // Calculate area
        int intersection_length = intersection.size();
        if (intersection_length <= 2) {
            return 0.0;
        }
        double area = intersection.get(intersection_length - 1).x * intersection.get(0).y -
                intersection.get(intersection_length - 1).y * intersection.get(0).x;
        for (int i = 0; i < intersection_length - 1; i++) {
            area += intersection.get(i).x * intersection.get(i + 1).y -
                    intersection.get(i).y * intersection.get(i + 1).x;
        }
        area *= 0.5;
        return area;
    }

    public static double rectangle_intersection_perimeter(Vector2d[] r_1, Vector2d[] r_2) {
        // Use the vertices of the first rectangle as the starting vertices of the intersection polygon.
        ArrayList<Vector2d> intersection = new ArrayList<>(Arrays.asList(r_1));

        // Loop over the edges of the second rectangle
        // Determine the intersection of the two rectangle r_1 and r_2
        for (int i = 0; i < r_2.length; i++) {
            if (intersection.size() <= 2) {
                break; // No intersection
            }

            // define a line with the points r_2[i] and r_2[(i + 1) % r_2.length] (adjacent vertices of r_2)
            // a, b, c are constants that determine the equation of the line: ax + by + c = 0
            Vector2d v1 = r_2[i];
            Vector2d v2 = r_2[(i + 1) % r_2.length];
            double a = v2.y - v1.y;
            double b = v1.x - v2.x;
            double c = v2.x * v1.y - v2.y * v1.x;
            // Any point p with line(p) <= 0 is on the "inside" (or on the boundary),
            // any point p with line(p) > 0 is on the "outside".

            // Loop over the edges of the intersection polygon, and determine which part
            // is inside and which is outside
            ArrayList<Vector2d> new_intersection = new ArrayList<Vector2d>();
            double[] line_values = new double[intersection.size()];
            for (int j = 0; j < intersection.size(); j++) {
                line_values[j] = a * intersection.get(j).x + b * intersection.get(j).y + c;
            }
            for(int j = 0; j < intersection.size(); j++) {
                double s_value = line_values[j];
                double t_value = line_values[(j + 1) % line_values.length];

                if(s_value <= 0) {
                    // TODO: change this to add a newly created vector that is a copy of intersection.get(j) ?
                    new_intersection.add(intersection.get(j));
                }
                if (s_value * t_value < 0) {
                    // Points are on opposite sides. Add the intersection of the lines to new_intersection.
                    Vector2d s = intersection.get(j);
                    Vector2d t = intersection.get((j + 1) % intersection.size());
                    double a2 = t.y - s.y;
                    double b2 = s.x - t.x;
                    double c2 = t.x * s.y - t.y * s.x;
                    Vector2d intersection_point = line_intersection(a, b, c, a2, b2, c2);
                    new_intersection.add(intersection_point);
                }
            }

            // set intersection to new intersection
            intersection = new_intersection;
        }

        // Calculate perimeter
        int intersection_length = intersection.size();
        if (intersection_length <= 2) {
            return 0.0;
        }
        double perimeter = distance(intersection.get(intersection_length - 1), intersection.get(0));
        for (int i = 0; i < intersection_length - 1; i++) {
            perimeter += distance(intersection.get(i), intersection.get(i + 1));
        }
        return perimeter;
    }

    public static Vector2d[] rectangle_vertices(double cx, double cy, double w, double h, double angle) {
        double dx = w / 2;
        double dy = h / 2;
        double dxcos = dx * Math.cos(angle);
        double dxsin = dx * Math.sin(angle);
        double dycos = dy * Math.cos(angle);
        double dysin = dy * Math.sin(angle);
        Vector2d[] rectangle = {new Vector2d(cx + -dxcos - -dysin, cy + -dxsin + -dycos),
                new Vector2d(cx + dxcos - -dysin, cy + dxsin + -dycos),
                new Vector2d(cx + dxcos - dysin, cy + dxsin + dycos),
                new Vector2d(cx + -dxcos - dysin, cy + -dxsin + dycos)};
        return rectangle;
    }
}

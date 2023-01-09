import java.util.Arrays;

/**
 * Class to solve cubic equations using the cubic formula
 * 
 * @author Akhilesh Khambekar
 *
 */
public class CubicSolver {
    
    /**
     * 
     * @param a
     * @param b
     * @param c
     * @param d
     *          coefficients of the cubic best fit of the stress-strain curve of each modular element
     *          representing the equation as ax^3 + bx^2 + cx + d
     *          
     * @return solutions to the equation, in ascending order
     */
    public double[] solve(double a, double b, double c, double d) {
        
        double bNorm = b/a;
        double cNorm = c/a;
        double dNorm = d/a;
        
        double Q = ((3.0 * cNorm) - (bNorm * bNorm)) / 9.0;
        double R = ((9.0 * bNorm * cNorm) - (27.0 * dNorm) - (2.0 * Math.pow(bNorm, 3.0))) / 54.0;
        
        // discriminant that determines number of roots
        double D = Math.pow(Q, 3.0) + Math.pow(R, 2.0);
        
        if (D < 0.0) {
            
            double theta = Math.acos(R / (Math.sqrt(0.0 - Math.pow(Q, 3))));
            double[] res = new double[3];
            
            res[0] = 2 * Math.sqrt(0.0 - Q)*Math.cos(theta / 3.0) - bNorm / 3.0;
            res[1] = 2 * Math.sqrt(0.0 - Q)*Math.cos((theta + 2.0 * Math.PI) / 3.0) - bNorm / 3.0;
            res[2] = 2 * Math.sqrt(0.0 - Q)*Math.cos((theta + 4.0 * Math.PI) / 3.0) - bNorm / 3.0;
            
            Arrays.sort(res);
            return res;
            
        } else if (D > 0.0) {
            
            double S = Math.cbrt(R + Math.sqrt(D));
            double T = Math.cbrt(R - Math.sqrt(D));
            double[] res = {S + T - (bNorm / 3.0)};
            return res;
            
        } else {
            double[] res = new double[3];
            
            res[0] = 2 * Math.cbrt(R) - (bNorm / 3.0);
            res[1] = Math.cbrt(R) - (bNorm / 3.0);
            res[2] = Math.cbrt(R) - (bNorm / 3.0);
            
            Arrays.sort(res);
            return res;
        }
    }
}
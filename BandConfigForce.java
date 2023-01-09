import java.util.Arrays;
/**
 * 
 * This class is a representation of one individual design or modular element/rubber-band configuration
 * This class stores information about the configuration itself and also calculates
 * mechanical properties such as estimated spring constant, energy absorbed and extension under applied force
 * 
 * @author Akhilesh Khambekar
 *
 */
public class BandConfigForce {
    
    // 2d array to store the configuration, each row represents a group of parallel elements
    public int[][] config;
    
    // estimated spring constant of each band for estimates in N/m. Should be manually assigned  
    double bandK;

    // minimum parallel elements per layer Should be manually assigned
    double minPerLayer;
    
    // maximum parallel elements per layer Should be manually assigned
    double maxPerLayer;
    
    // coefficients of the cubic best fit of the stress-strain curve of each modular element
    // representing the equation as ax^3 + bx^2 + cx + d
    double a = 280.35;
    double b = -147.69;
    double c = 40.545;
    double d = 0.4301;
    
    // cubic equation solver for 
    CubicSolver cubic;
    
    /**
     * 
     * Constructor for a new BandConfigForce
     * 
     * @param init the 2d matrix that represents the configuration
     *             0 represents no element and 1 means element exists in the corresponding cell
     *        
     * @param a
     * @param b
     * @param c
     * @param d
     * 
     *             coefficients of the cubic best fit of the stress-strain curve of each modular element
     *             representing the equation as ax^3 + bx^2 + cx + d
     */
    public BandConfigForce (int[][] init, double a, double b, double c, double d) {
        
        config = init;
        
        //parameters to manually tune
        bandK = 22.5;
        minPerLayer = 4.0;
        maxPerLayer = Math.max(init[0].length, 100);
        
        // initialize equation solver
        cubic = new CubicSolver();
    }
    
    /**
     * 
     * Estimates the Spring Constant based on treating each element as a linear spring
     * Uses the manually entered spring constant for each element (assumes all elements are identical)
     * 
     * @return returns the Spring constant in N/m, returns 0 for invalid designs which have at least a layer with 0 elements
     */
    public double calculateK() {
        
        // series spring add in reciprocal
        double kValInv = 0.0;
        
        // calculate spring constant for each layer
        for(int i = 0; i< config.length; i++) {
            double layerK = 0.0;
            int layerBands = 0;
            for (int j = 0; j < config[i].length; j++) {
                if (config[i][j] == 1) {
                    layerK += bandK;
                    layerBands++;
                }
                                
            }
            
            // ensures number of bands per layer is within limits
            if (layerBands < minPerLayer) {
                return 0.0;
            }
            if (layerBands > maxPerLayer) {
                return 0.0;
            }
            
            // calculating the series result
            kValInv += 1.0/layerK;
        }    
        double kVal = 1.0/kValInv;
        return kVal;  
    }
    
    /**
     * Integrates the cubic stress-strain curve line of best fit that describes each rubber band or other modular element
     * from  an extension 0 to y to calculate energy absorbed
     * 
     * @param y the upper limit of integration (0 is the lower limit)
     * 
     * @return the value of the integral
     */
    public double integral(double y) {
        
        double coeff4 = (1.0 / 4.0) * a;
        double coeff3 = (1.0 / 3.0) * b;
        double coeff2 = (1.0 / 2.0) * c;
        double coeff1 = (1.0) * d;
        
        return (coeff4 * Math.pow(y, 4)) + (coeff3 * Math.pow(y, 3)) + (coeff2 * Math.pow(y, 2)) + (coeff1 * Math.pow(y, 1));
        
    }
    
    /**
     * 
     * Calculates the extension of the rubber band or modular element based on the stress-strain curve
     * The stress-strain curve is represented as its cubic line of best fit
     * 
     * @param force the force applied to the element
     * @return the extension after applying the aforementioned force
     */
    public double solveCubic(double force) {
        
        double[] result = cubic.solve(a, b, c, d - force);
        
        // return the smallest positive solution
        return result[0];
    }
    
    /**
     * 
     * Calculates the energy absorbed by the Configuration depending on the applied force
     * 
     * @param force The applied force
     * @return The absorbed energy in Joules, returns 0 for invalid designs
     */
    public double calculateEnergy(double force) {
        
        double energy = 0.0;
        
        
        for (int i = 0;  i < config.length; i++) {
            
            double layerBands = 0;
            
            // calculate energy absorbed by each layer
            for (int j = 0; j < config[i].length; j++) {
                
                if (config[i][j] == 1) {
                    
                    layerBands = layerBands + 1;
                    
                }               
            }
            
            // return 0 for invalid designs with at least one layer with 0 elements
            if (layerBands < minPerLayer) {
                return 0.0;
             }
            
            double forcePerBand = force/((double)layerBands);

            double layerExten = solveCubic(forcePerBand);

            if (layerExten > 0.5) {
                return 0.0;
            }
            energy += integral(layerExten)*((double)layerBands);            
        }
        return energy;
    }
    
    /**
     * 
     * Calculates the amount the configuration will extend when a force is applied
     * 
     * @param force The applied force
     * @return The resulting extension in meters. returns 0 for invalid designs
     */
    public double calculateExten(double force) {
        double exten = 0.0;
        
        for (int i = 0;  i < config.length; i++) {
            int layerBands = 0;
            
            // calculate extension for each layer
            for (int j = 0; j < config[i].length; j++) {
                if (config[i][j] == 1) {
                    layerBands++;
                }
                
            }
            
            // return 0 for invalid designs
            if (layerBands < minPerLayer) {
                return 0.0;
            }
            double forcePerBand = force/((double)layerBands);
            double layerExten = solveCubic(forcePerBand);
            if (layerExten > 0.5) {
                
            }
            exten += layerExten;
        }
        return exten;
    }
    
    
    /**
     * 
     * Prints the configuration layer by layer as a matrix of 0s and 1s
     * A 0 represents no element in the corresponding position and 1 represents an element
     * 
     */
    public void printConfig() {
        for (int i = 0; i < config.length; i++) {
            for (int j = 0; j < config[i].length; j++) {
                System.out.print(config[i][j] + " ");
            }
            System.out.println();    
        }
    }
    
    /**
     * 
     *  Prints the configuration layer by layer as a matrix of 0s and 1s
     *  A 0 represents no element in the corresponding position and 1 represents an element
     *  Also enumerates the number of parallel elements in each layer
     */
    public void printConfigBetter() {
        for (int i = 0; i < config.length; i++) {
            Arrays.sort(config[i]);
            int countBands = 0;
            for (int j = config[i].length - 1; j > -1; j--) {
                if (config[i][j] == 1) {
                    countBands++;
                }
                System.out.print(config[i][j] + " ");
            }
            System.out.println(" ---> Layer "+(i+1)+": " + countBands + " bands parallel" );    
        }
    }  
}
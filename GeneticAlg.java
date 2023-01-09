
/**
 * The goal of this project is to design the optimal bungee cord out of modular elements
 * that are rubber bands or springs, given constraints of Drop height, Max G load and Jumper mass.
 * 
 *  This class implements a genetic/evolutionary algorithm to determine the optimal configuration
 *  of the modular elements to satisfy the aforementioned constraints
 * 
 * @author Akhilesh Khambekar
 *
 */
public class GeneticAlg {
   // the number of individuals per generation
   int populationSize;
   
   // the number of the best performing individuals selected for the next generation
   int selected;
   
   // stores the individuals in the current generation
   // Please read the BandConfigForce constructor for parameters that need to be manually tuned
   BandConfigForce[] population;
   
   // stores error for each individual
   double[] Edevs;
   
   // acceleration due to gravity (m/s^2)
   final double g = 9.8;
   
   // energy to be absorbed by the bungee cord. Equal to mass * g * drop height 
   double targetE;
   
   // force exerted on the cord at max G load. Equal to mass* G limit
   double targetForce;
   
   // max number of elements in series
   int configLen;
   
   // max number of element in parallel
   int configWid;
   
   //error percent for a design to be considered 'good'
   //should be a decimal. 5% is 0.05
   double goodCrit;
   
   // coefficients of the cubic best fit of the stress-strain curve of each modular element
   // representing the equation as ax^3 + bx^2 + cx + d
   double a;
   double b;
   double c;
   double d;
   
   /**
    * 
    * Constructor for EvoltutionaryAlg class
    * Please read the BandConfigForce constructor for parameters that need to be manually tuned
    * 
    * @param popSize        initialize the number of individuals per generation
    * @param numSel         initialize the number of best performing individuals selected for the next generation
    * @param maxLen         initialize the max number of elements in series
    * @param maxWid         initialize the max number of element in parallel
    * @param mass           initialize the mass of the jumper (in kilograms)
    * @param height         initialize the drop height (in meters)
    * @param Gmax           initialize the max load (in Gs, 49 m/s^2 is 5 Gs)
    * @param err Allowable  initialize the allowable error percent (should be a decimal. 5% is 0.05)
    * 
    * 
    * @param a
    * @param b
    * @param c
    * @param d
    *  
    *              coefficients of the cubic best fit of the stress-strain curve of each modular element
    *              representing the equation as ax^3 + bx^2 + cx + d
    */
   public GeneticAlg(int popSize, int numSel,  int maxLen, int maxWid, double mass, double height, double Gmax, 
           double errAllowable, double a, double b, double c, double d) {
       
       populationSize = popSize;
       selected = numSel;
       
       population = new BandConfigForce[populationSize];
       Edevs = new double[populationSize];
       
       configLen = maxLen;
       configWid = maxWid;
       goodCrit = errAllowable;
       
       this.a = a;
       this.b = b;
       this.c = c;
       this.d = d;
       
       targetE = mass * g * height;
       targetForce = mass * g * Gmax;      
   }
   
   /**
    * This method generates the individuals for the first generation randomly
    */
   public void generate() {
       
       // performance metrics: counting number of valid and invalid individuals
       // an invalid design is one which is not continuous, one where there is at least one layer with 0 elements 
       int invalid = 0;
       int numValid = 0;
       
       // initialize total error
       double totalError = 0.0;
       
       // initialize the min and max energy absorbed by the bungee cord
       double min = Double.MAX_VALUE;
       double max = Double.MIN_VALUE;
       
       // count the well performing designs
       int countGood = 0;
       
       // count the designs which match the target parameters perfectly
       int countPerfect = 0;
       
       // generating the individuals, invalid ones are dropped
       while (numValid < populationSize) {
           
           // Initializing the individual
           int[][] init = new int[configLen][configWid];
           
           // each row
           for (int i = 0; i < configLen; i++) {
               
               // each element in row
               for (int j = 0; j < configWid; j++) {
                   
                   // 1 means cell contains element, 0 means not
                   // randomly generates 1 or 0, based on comparing a random number between 0 and 1 with 0.5
                   double r = Math.random();
                   if (r < 0.5) {  
                       
                       init[i][j] = 1;
                       
                   } else {
                       
                       init[i][j] = 0;
                       
                   }
               }
           }
       
       // creating the instance
       BandConfigForce bNew = new BandConfigForce(init, a, b, c, d);
       
       // calculating the energy absorbed by the current instance
       double currEnergy = bNew.calculateEnergy(targetForce);
           // only add it to the population if valid
           if (currEnergy != 0) {     
               
               population[numValid] = bNew;
               
               // calculate error for current individual
               Edevs[numValid] = Math.abs(currEnergy - targetE);  
               
               // update max, min, number of good and perfect designs
               if (Math.abs(currEnergy - targetE) < targetE * goodCrit) {
                   
                   countGood++;
                   
               }
               if (currEnergy == targetE) {
                   
                   countPerfect++;
                   
               }
               if (min > currEnergy) {
                   
                   min = currEnergy;
                   
               }
               if (max < currEnergy) {
                   
                   max = currEnergy;
                   
               }
               
               // updating total err for average err calculation
               totalError += Math.abs(currEnergy - targetE);
               numValid++;
           } else {
               
               invalid++;
               
           }               
       }
       
       // Sort the individuals by performance
       quickSort(Edevs, population, 0, populationSize - 1);
       
       // printing the performance metrics
 
       System.out.println("Smallest Error: " + Edevs[0]);
       System.out.println("Best design energy absorbed (Joules): " + population[0].calculateEnergy(targetForce));
       System.out.println("Best design stretch amount (Meters): " + population[0].calculateExten(targetForce));
       System.out.println("Best design: ");
       
       population[0].printConfigBetter();
       System.out.println("\nGeneration info:");
       System.out.println("average error (Joules): " + (totalError / populationSize));
       System.out.println("Max energy absorbed (J): " + max);
       System.out.println("Min energy absorbed (J): " + min);
       System.out.println("Good designs count: " + countGood);
       System.out.println("Perfect designs count: " + countPerfect);
       System.out.println("Valid Designs: "+ numValid);
       System.out.println("Invalid Designs: "+ invalid + "\n");
   }
   
   /**
    * 
    * This method 'breeds' 2 individuals in order to produce individuals for the next generation
    * It combines features from parent 1 and 2 by choosing the corresponding cell from parent 1 or from parent 2
    * with equal probability
    * 
    * @param b1 parent 1
    * @param b2 parent 2
    * @return the 'child' from the 2 parents
    * 
    * @throws IllegalArgumentException if sizes do not match
    */
   public BandConfigForce breed(BandConfigForce b1, BandConfigForce b2) {
       
       // throws IllegalArgumentException if parents have incorrect size
       int len1 = b1.config.length;
       int wid1 = b1.config[0].length;
       
       int len2 = b2.config.length;
       int wid2 = b2.config[0].length;
       
       if (len1 != configLen || len2 != configLen || wid1 != configWid || wid2 != configWid) {
           
           throw new IllegalArgumentException();
       }
       
       // initialize the child
       int[][] init = new int[configLen][configWid];
       
       // breeding loop
       for (int i = 0; i < configLen; i++) {
           
           for (int j = 0; j < configWid; j++) {
               
               // if random num is less than 0.5, choose corresponding cell from parent 1, else from parent 2 
               double r = Math.random();
               if (r < 0.5) {
                   
                   init[i][j] = b1.config[i][j];
                   
               } else {
                   
                   init[i][j] = b2.config[i][j];
                   
               }
           }
       }
       return new BandConfigForce(init, a, b, c, d);
   }
   
   /**
    *  produces the next generation by eliminating the worst performing individuals from the previous
    *  generation and replacing them with the children of the best performing ones
    */
   public void evolve() {
       
       // count number of valid and invalid designs
       int numValid = selected;
       int invalid = 0;
       
       // initialize total error, max, min for energy absorbed in this generation
       double totalError = 0.0;
       double min = Double.MAX_VALUE;
       double max = Double.MIN_VALUE;
       
       // count the well performing designs
       int countGood = 0;
       
       // count the designs which match the target parameters perfectly
       int countPerfect = 0;
       
       // replace all except the best performing individuals
       for (int i = selected; i < populationSize; i++) {
           
           // initialize the child
           BandConfigForce b = null;
           
           // only add valid individuals to the population
           boolean validConf = false;
           while (!validConf) {
               
               // choose parents randomly
               int parent1 = (int)(Math.random() * selected);
               int parent2 = (int)(Math.random() * selected);
               
               // create the child
               b = breed(population[parent1], population[parent2]);
               
               // calculate energy absorbed by current child
               double currEnergy = b.calculateEnergy(targetForce);
               
               // only add valid individuals to the populations
               if (currEnergy != 0) {
                   
                   validConf = true;
                   
                   // add child to the population
                   population[i] = b;
                   Edevs[i] = Math.abs(currEnergy - targetE);
                   
                   // update max, min, good, perfect counts 
                   if (Math.abs(currEnergy - targetE) < targetE * goodCrit) {
                       
                       countGood++;
                       
                   }
                   if (currEnergy == targetE) {
                       
                       countPerfect++;
                       
                   }
                   if (min > currEnergy) {
                       
                       min = currEnergy;
                       
                   }
                   if (max < currEnergy) {
                       
                       max = currEnergy;
                       
                   }
                   
                   // for average error
                   totalError += Math.abs(currEnergy - targetE);
                   
                   numValid++;
                   
               } else {
                   
                   invalid++;
                   
               }
           }
           population[i] = b;
       }
       
       // sort by performance
       quickSort(Edevs, population, 0, populationSize-1);

       // printing the performance metrics
       
       System.out.println("Smallest Error: " + Edevs[0]);
       System.out.println("Best design energy absorbed (Joules): " + population[0].calculateEnergy(targetForce));
       System.out.println("Best design stretch amount (Meters): " + population[0].calculateExten(targetForce));
       System.out.println("Best design: ");
       
       population[0].printConfigBetter();
       System.out.println("\nGeneration info:");
       System.out.println("average error (Joules): " + (totalError / populationSize));
       System.out.println("Max energy absorbed (J): " + max);
       System.out.println("Min energy absorbed (J): " + min);
       System.out.println("Good designs count: " + countGood);
       System.out.println("Perfect designs count: " + countPerfect);
       System.out.println("Valid Designs: "+ numValid);
       System.out.println("Invalid Designs: "+ invalid + "\n");
   }
   
   /**
    *
    * Main method to run the evolutionary algorithm
    * 
    */
   public static void main(String[] args) {
       
       //store the parameters to be used for the algorithm run
       int populationSize = 1000;
       int numSelected = 100;
       int maxSeries = 8;
       int maxParallel = 7;
       double jumperMassKg = 0.75;
       double dropHeightMeters = 6.8;
       double maxG = 5.0;
       double allowedErr = 0.05;
       double a = 280.35;
       double b = -147.69;
       double c = 40.545;
       double d = 0.4301;
       
       //max number of generations
       int maxGens = 20;
       
       // instantiate the driver class
       GeneticAlg model = new GeneticAlg(populationSize, numSelected, maxSeries,
               maxParallel, jumperMassKg, dropHeightMeters, maxG, allowedErr, a, b, c, d);
       
       System.out.println("generation 1");
       model.generate();
       
       //evolution loop
       for (int i = 2; i <= maxGens; i++) {
           
           System.out.println("Generation " + i + ":");
           model.evolve();
           
       }              
   }
   
   /**
    * Swaps 2 elements in 2 arrays in parallel, for QuickSorting
    * the individuals by performance
    * 
    * @param devArr performance metric array
    * @param configArr the array of individuals to be sorted
    * 
    *        devArr and configArr must be aligned such that the performance metric for individual i
    *        is located at devArr[i]
    * 
    * @param i index of first element
    * @param j index of other element
    */
    void swap(double[] devArr, BandConfigForce[] configArr, int i, int j) {
 
       double temp = devArr[i];
       BandConfigForce temp2 = configArr[i];
       devArr[i] = devArr[j];
       configArr[i] = configArr[j];
       devArr[j] = temp;
       configArr[j] = temp2;
    }
   
    /**
     * partitions the arrays around a pivot
     * 
     * @param devArr performance metric array
     * @param configArr the array of individuals to be sorted
     * 
     *        devArr and configArr must be aligned such that the performance metric for individual i
     *        is located at devArr[i]
     * 
     * @param i index of first element
     * @param j index of other element
     */
    int partition(double[] devArr, BandConfigForce[] configArr, int low, int high) {
 
        // initializing the pivot
        double pivot = devArr[high];

        int i = (low - 1);
 
        for (int j = low; j <= high - 1; j++) {
            
            if (devArr[j] < pivot) {
                
                i++;
                swap(devArr, configArr, i, j);
                
            }
        }
        swap(devArr, configArr, i + 1, high);
        return (i + 1);
    } 
    
    /**
     * 
     * Implements Quicksort to sort the individuals in a generation by performance
     * 
     *@param devArr performance metric array
     * @param configArr the array of individuals to be sorted
     * 
     *        devArr and configArr must be aligned such that the performance metric for individual i
     *        is located at devArr[i]
     *        
     * @param low start index for sorting
     * @param high end index for sorting
     */
    
     void quickSort(double[] devArr, BandConfigForce[] configArr, int low, int high) {
  
        if (low < high) {
 
            int pi = partition(devArr, configArr, low, high);
            
            // recursively sort each half
            quickSort(devArr, configArr, low, pi - 1);
            quickSort(devArr, configArr, pi + 1, high);
            
        }
    }
}
package com.company;
import java.util.Random;

public class Personnel {

    private static int maxC = 20000; // Max jobs
    private static int maxRun = 10; // Max number of runs
    private static int  maxS = 10; // Max number of servers per station
    private static int  maxAS = 5; // Max number of arrival sources
    private static int maxNrStations = 10; // max number of stations
    private static int maxNrJobTypes = 10; // max number of job types

    /* COUNTER */
    private double j1, j2, j3, l1; // define float numbers
    private int i1, i2, i6, run, i3;  // define counters (integers)
    private double K, s0, L;
    private double[] avg = new double[30];
    private char[] naam = new char[300];
    private char[] sProblem = new char[10];
    private double leftVarTriangular, right_var_triangular; // inputs for triangular distribution

    /* INPUT DATA RELATED TO RADIOLOGY DPT */
    private int nrStations;   //Number of workstations
    private int[] nrServersPerStation = new int[maxNrStations]; //Number of servers per workstation

    /* INPUT DATA RELATED TO SYSTEM JOBS */
    private int nrJobTypes; // Number of job types
    private int[] nrWorkstationsPerJobType = new int[maxNrJobTypes]; // Number of workstations per job type
    private int[][] route = new int[maxNrJobTypes][maxNrStations]; // Route to follow for each job type
    private int[] currentStation = new int[maxC];  // Matrix that denotes the current station of a scan (sequence number)

    /* GENERAL DISCRETE EVENT SIMULATION PARAMETERS */
    private double t; // Simulation time
    private int N; // number of scans (Stop criterion)

    /* VARIABLES RELATED TO system SCANS */
    private int n; // Number of scans in the system
    private int[] nWS = new int[maxNrStations]; // Number of scans at a particular workstation

    private double[] meanCustomersSystem = new double[maxRun];
    private double[] totN = new double[maxRun]; // Number of customers in the system over time
    private double[][] totNWS = new double[maxRun][maxNrStations]; // Number of customers in a workstation over time

    /* PARAMETERS RELATED TO ARRIVAL OF SCANS */
    private int nrArrivalSources; // Number of arrival sources;
    private double[] lambda = new double[maxAS]; // Arrival rate
    private double[][] cumDistrScans = new double[maxAS][maxNrJobTypes]; // Cumulative(!) Distribution of job types
    // per source
    private int nA; // Number of scans arrived to the system
    private int[] nAWS = new int[maxNrStations]; // Number of scans arrived to a particular workstation
    private double[] tA = new double[maxAS];                        // Time of next arrival for each source
    private double firstTA;                     // First arrival time over all sources;
    private int indexArr;                       // Source of next arrival;
    private double tLambda;
    private double[][] totLambda = new double[maxRun][maxAS];  // Total arrival rate and arrival rate per source and
    // run
    private int[] scanType = new int[maxC];     // Type of scan arriving
    private double[][] timeArrival = new double[maxRun][maxC]; // Time of arrival of the scan to the ancillary services
    private double[][][] timeArrivalWS = new double[maxRun][maxNrStations][maxC];        // Time of arrival of a
    // particular scan
    // to a
    // particular
    // workstation
    private double[] meanInterarrivalTime = new double[maxRun];      // Calculated average interarrival time per source

    /* PARAMETERS RELATED TO PROCESSING OF SCANS */
    private double[][] mu = new double[maxNrStations][maxNrJobTypes];     // Processing rate per station and job type
    private double[][] var = new double[maxNrStations][maxNrJobTypes];     // Variance per station and job type
    private double[][] sigma = new double[maxNrStations][maxNrJobTypes];     // Standard deviation per station and job type
    private int nD;        // Number of scans handled
    private int[] nDWS = new int[maxNrStations];     // Number of scans handled in a particular workstation
    private double[][] tD = new double[maxNrStations][maxS];         // Time of next departure for each server in each workstation
    private double firstTD;                 // First departure time over all sources
    private int indexDepStation;                // Station with first departure
    private int indexDepServer;                 // Server with first departure
    private double[] meanServiceTime = new double [maxNrStations];   // Calculated average service time
    private double tMu;                         // Generated service time
    private double[] totMu = new double[maxRun];                     // Total service time generated
    private double[][][] timeService = new double[maxRun][maxNrStations][maxC]; // Service time per customer and workstation
    private int[][] currentCust = new int [maxNrStations][maxS];       // Customer handles by a particular workstation and server
    private int[][] listScan = new int [maxNrStations][maxC];    // List of customers processed at a particular workstation on a particular point in time

    /* PARAMETERS RELATED TO WAITING OF SCANS */
    private double[] meanWaitingTime = new double[maxRun];           // Calculated average waiting time per source
    private double[] waitingTime = new double[maxRun];               // Waiting time for a particular customer
    private double[][][] waitingTimeJobWS = new double[maxRun][maxNrStations][maxC];      // Waiting time for a job on a particular workstation
    private double[] meanCustomersQueue = new double[maxRun];        // Calculated average number of scans in queue per source
    private double[] totNQueue = new double[maxRun];                 // Total number of scans in queue over time
    private double[][] totNQueueWS = new double[maxRun][maxNrStations];             // Total number of scans in queue at workstation over time

    /* VARIABLES RELATED TO PROCESSED SCANS */
    private double[][] timeDeparture = new double[maxRun][maxC];           // Time of departure of the scan from the ancillary services
    private double[][][] timeDepartureWS= new double[maxRun][maxNrStations][maxC];       // Time of departure of a particular scan from a particular workstation
    private double[][] timeSystem= new double[maxRun][maxC];              // Time in system for a particular customer
    private double[][][] timeSystemJobWS= new double[maxRun][maxNrStations][maxC];       // Time in system for a particular job on a particular workstation
    private int[] orderOut= new int[maxC];                     // Order in which jobs are coming out of the system
    private double[] meanSystemTime= new double[maxRun];            // Calculated average time in system per source

    /* OTHER PARAMETERS */
    private int infinity;                       // Value for infinity
    private double[][][] idle = new double[maxRun][maxNrStations][maxS];                  // Idle time for server s at workstation w
    private double[][] rhoWSS = new double [maxNrJobTypes][maxS];                  // Utilization of server s at workstation w
    private double[] rhoWS= new double[maxNrStations];                     // Utilization of workstation w
    private double rho;                         // Overall utilization

    /* VARIABLES RELATED TO CLOCK TIME */
    private double elapsedTime, timeSubproblem;
    private long startTime, interTime, projectStartTime; // Time measurements to compute run times


    public double exponentialDistribution(double lambda) {

        Random rand = new Random();
        j1 = rand.nextDouble(); // Random number
        if (j1 == 0) { // Random number should lie in interval ]0,1[
            j1 += 0.0001;
        }
        j2 = -Math.log(j1)/lambda; // Inversion method
        return j2;
    }

    public int poissonDistribution(double lambda) {
        Random rand = new Random();
        double k, L;
        int p;
        j1 = rand.nextDouble(); // Random number
        k = 0;
        L = Math.exp(-lambda);
        j3 = 0;
        do { // Inversion method
            j2 = L * Math.pow(lambda, k);
            p = 1;
            for (i6 = 0; i6 <= k; i6++) {
                if (i6 == 0) {
                    p = 1;
                } else {
                    p *= i6;
                }
            }
            j2 /= p;
            j3 += j2;
            k++;
        } while (j1 >= j3);
        return (int) (k-1);
    }

    public int normalDistribution(double mean, double stdev) {
        Random rand = new Random();
        double v1, v2, t;
        int x;
        do {
            v1 = rand.nextDouble() * 2;
            v1 -= 1;
            v2 = rand.nextDouble() * 2;
            v2 -= 1;
            t = v1 * v1 + v2 * v2;
        } while (t >= 1 || t == 0);
        double multiplier = Math.sqrt(-2 * Math.log(t) / t);
        x = (int) (v1 * multiplier * stdev + mean);
        return x;
    }

    public int bernouilliDistribution(double prob) {
        Random rand = new Random();
        j1 = rand.nextDouble(); // random number
        if (j1 < prob) { // Inversion method
            return 0;
        } else {
            return 1;
        }
    }

    public int uniformDistribution(double a, double b) {
        Random rand = new Random();
        int x;
        j1 = rand.nextDouble(); // random number
        x = (int) (a + (b-a) * j1); // Inversion method
        return x;
    }

    public int triangularDistribution(int a, int b, int c) {
        double mean, stdev;
        double x, L;
        mean = (a + b + c) / 3;
        stdev = (Math.pow(a, 2) + Math.pow(b, 2) + Math.pow(c, 2) - a * b - a * c - b * c) / 18;
        stdev = Math.sqrt(stdev);
        Random rand = new Random();
        j1 = rand.nextDouble();
        x = a;
        do {
            if (x <= b) {
                L = Math.pow((x-a), 2) / ((c-a) * (b-a));
            } else {
                L = 1 - (Math.pow(c-x, 2) / ((c-a) * (c-b)));
            }
            x++;
        } while (j1 >= L);
        return (int) x-1;
    }


    public static void main(String[] args) {
	// write your code here
    }
}

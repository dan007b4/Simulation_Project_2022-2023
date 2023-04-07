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

    public void initialize_functions(){

        /* PUT ALL VARIABLES TO ZERO */

        /* INPUT DATA RELATED TO RADIOLOGY DPT */
        nrStations = 0;

        for (i1 = 0; i1 < maxNrStations; i1++)
            nrServersPerStation[i1] = 0;

        /* INPUT DATA RELATED TO SYSTEM JOBS */
        nrJobTypes = 0;


        for (i2 = 0; i2 < maxNrJobTypes; i2++)
        {   nrWorkstationsPerJobType[i2] = 0;
            for (i1 = 0; i1 < maxNrStations; i1++)
                route[i2][i1] = 0;
        }
        for (i1 = 0; i1 < maxC; i1++)
            currentStation[i1] = 0;

        /* GENERAL DISCRETE EVENT SIMULATION PARAMETERS */

        t = 0;
        N = 0;

        /* VARIABLES RELATED TO system SCANS */
        n = 0;
        for (i1 = 0; i1 < maxNrStations; i1++)
            nWS[i1] = 0;

        for(i2 = 0; i2 < maxRun; i2++)
        {   meanCustomersSystem[i2] = 0;
            totN[i2] = 0;
            for (i1 = 0; i1 < maxNrStations; i1++)
                totNWS[i2][i1] = 0;
        }

        /* PARAMETERS RELATED TO ARRIVAL OF SCANS */
        nrArrivalSources = 0;
        nA = 0;
        firstTA = 0;
        indexArr = 0;
        tLambda = 0;

        for (i1 = 0; i1 < maxNrStations; i1++)
            nAWS[i1] = 0;

        for(i2 = 0; i2 < maxRun; i2++)
        {   meanInterarrivalTime[i2] = 0;
            for (i6 = 0; i6 < maxAS; i6++)
                totLambda[i2][i6] = 0;


            for (i3 = 0; i3 < maxC; i3++)
            {    timeArrival[i2][i3] = 0;
                for (i1 = 0; i1 < maxNrStations; i1++)
                    timeArrivalWS[i2][i1][i3] = 0;

            }
        }
        for (i3 = 0; i3 < maxC; i3++)
        {   scanType[i3] = 0;
        }


        for (i6 = 0; i6 < maxAS; i6++)
        {   tA[i6] = 0;
            lambda[i6] = 0;
            for (i3 = 0; i3 < maxNrJobTypes; i3++)
                cumDistrScans[i6][i3] = 0;

        }

        /* PARAMETERS RELATED TO Processing OF SCANS */

        for (i1 = 0; i1 < maxNrStations; i1++)
        {   nDWS[i1] = 0;
            for (i3 = 0; i3 < maxNrJobTypes; i3++)
            {   mu[i1][i3] = 0;
                var[i1][i3] = 0;
                sigma[i1][i3] = 0;
            }
            for (i6 = 0; i6 < maxS; i6++)
            {   tD[i1][i6] = 0;
                currentCust[i1][i6] = 0;
            }
            for (i6 = 0; i6 < maxC; i6++)
                listScan[i1][i6] = -1;

        }
        nD = 0;
        firstTD = 0;
        indexDepStation = 0;
        indexDepServer = 0;
        tMu = 0;

        for(i2 = 0; i2 < maxRun; i2++)
        {   meanServiceTime[i2] = 0;
            totMu[i2] = 0;
            for (i1 = 0; i1 < maxNrStations; i1++)
            {   for (i3 = 0; i3 < maxC; i3++)
                timeService[i2][i1][i3] = 0;


            }
        }
        /* PARAMETERS RELATED TO waiting OF SCANS */

        for(i2 = 0; i2 < maxRun; i2++)
        {
            meanWaitingTime[i2] = 0;
            waitingTime[i2] = 0;
            meanCustomersQueue[i2] = 0;
            totNQueue[i2] = 0;
            for (i1 = 0; i1 < maxNrStations; i1++)
            {   totNQueueWS[i2][i1] = 0;
                for (i3 = 0; i3 < maxC; i3++)
                    waitingTimeJobWS[i2][i1][i3] = 0;
            }
        }


        /* VARIABLES RELATED TO Processed SCANS */

        for(i2 = 0; i2 < maxRun; i2++)
        {   meanSystemTime[i2] = 0;
            for (i3 = 0; i3 < maxC; i3++)
            {   timeDeparture[i2][i3] = 0;
                timeSystem[i2][i3] = 0;

                for (i1 = 0; i1 < maxNrStations; i1++)
                {
                    timeDepartureWS[i2][i1][i3] = 0;
                    timeSystemJobWS[i2][i1][i3] = 0;
                }
            }
        }

        for (i3 = 0; i3 < maxC; i3++)
        {   orderOut[i3] = 0;

        }

        /* OTHER PARAMETERS */
        infinity = 0;

        for(i2 = 0; i2 < maxRun; i2++)
        {
            for (i3 = 0; i3 < maxNrStations; i3++)
            {

                for (i1 = 0; i1 < maxS; i1++)
                {
                    idle[i2][i3][i1] = 0;
                }
            }
        }
        rho = 0;
        for (i3 = 0; i3 < maxNrStations; i3++)
        {   rhoWS[i3] = 0;

            for (i1 = 0; i1 < maxS; i1++)
            {
                rhoWSS[i3][i1] = 0;
            }
        }
    }

    public Personnel() {
    }

    public void procedure(){
        L = 1;
        for (i3 = 0; i3 < L; i3++)                          // Count number of runs
        {   K = 1;
            for (run = 0; run < K; run++)                  // Count number of replications per run
            {   init();
                radiology_system();
                output();
            }
        }

        //TODO hier moet nog die outputstream komen

    }

    public void init(){
        /* PUT ALL VARIABLES TO ZERO */

        initialize_functions();


        /* SET INPUT VALUES */

        Random rand = new Random((long) ((i3+1)*K-run));               // Ensure you each time use a different seed to get IID replications

        /* INPUT RADIOLOGY DPT */

        nrStations = 5;                        // Number of workstations

        nrServersPerStation[0] = 3;                     // Input number of servers per workstation
        nrServersPerStation[1] = 2;
        nrServersPerStation[2] = 4;
        nrServersPerStation[3] = 3;
        nrServersPerStation[4] = 1;

        /* INPUT JOB TYPES */
        nrJobTypes = 4;                       // Number of scans types
        nrWorkstationsPerJobType[0] = 4;             // Number of workstations per job type
        nrWorkstationsPerJobType[1] = 3;
        nrWorkstationsPerJobType[2] = 5;
        nrWorkstationsPerJobType[3] = 3;

        route[0][0] = 2;                        // Route to follow for each job type (JOB = 1)
        route[0][1] = 0;                        // Note: Workstation i in assignment corresponds to workstation i-1 in code as here we start counting from 0
        route[0][2] = 1;
        route[0][3] = 4;

        route[1][0] = 3;                        // Route to follow for each job type (JOB = 2)
        route[1][1] = 0;
        route[1][2] = 2;

        route[2][0] = 1;                        // Route to follow for each job type (JOB = 3)
        route[2][1] = 4;
        route[2][2] = 0;
        route[2][3] = 3;
        route[2][4] = 2;

        route[3][0] = 1;                        // Route to follow for each job type (JOB = 4)
        route[3][1] = 3;
        route[3][2] = 4;

        /* INPUT ARRIVAL PROCESS */
        nrArrivalSources = 2;                 // Number of arrival sources
        // Arrival from radiology department
        lambda[0] = 1/0.25;                     // Input arrival rate = 1/mean interarrival time
        cumDistrScans[0][0] = 0.2;                   // Distribution scans (SOURCE = 1) - Cumulative distribution
        cumDistrScans[0][1] = 0.4;
        cumDistrScans[0][2] = 0.5;
        cumDistrScans[0][3] = 1;

        // Arrival from other services
        lambda[1] = 1/1;                        // Input arrival rate = 1/mean interarrival time
        cumDistrScans[1][0] = 0;                   // Distribution scans (SOURCE = 2) - Cumulative distribution
        cumDistrScans[1][1] = 0.4;
        cumDistrScans[1][2] = 0.4;
        cumDistrScans[1][3] = 1;


        /* INPUT SERVICE PROCESS */

        mu[0][0] = 12;                               //Processing time per ws and job type (WS1, J1)
        mu[0][1] = 15;
        mu[0][2] = 15;
        mu[0][3] = 0;
        mu[1][0] = 20;                               //Processing time per ws and job type (WS2, J1)
        mu[1][1] = 0;
        mu[1][2] = 21;
        mu[1][3] = 18;
        mu[2][0] = 16;                               //Processing time per ws and job type (WS3, J1)
        mu[2][1] = 14;
        mu[2][2] = 10;
        mu[2][3] = 0;
        mu[3][0] = 0;                               //Processing time per ws and job type (WS4, J1)
        mu[3][1] = 20;
        mu[3][2] = 24;
        mu[3][3] = 13;
        mu[4][0] = 25;                               //Processing time per ws and job type (WS5, J1)
        mu[4][1] = 0;
        mu[4][2] = 20;
        mu[4][3] = 25;
        var[0][0] = 2;                               //Processing variance per ws and job type (WS1, J1)
        var[0][1] = 2;
        var[0][2] = 3;
        var[0][3] = 0;
        var[1][0] = 4;                               //Processing variance per ws and job type (WS2, J1)
        var[1][1] = 0;
        var[1][2] = 3;
        var[1][3] = 3;
        var[2][0] = 4;                               //Processing variance per ws and job type (WS3, J1)
        var[2][1] = 2;
        var[2][2] = 1;
        var[2][3] = 0;
        var[3][0] = 0;                               //Processing variance per ws and job type (WS4, J1)
        var[3][1] = 3;
        var[3][2] = 4;
        var[3][3] = 2;
        var[4][0] = 5;                               //Processing variance per ws and job type (WS5, J1)
        var[4][1] = 0;
        var[4][2] = 3;
        var[4][3] = 5;
        for (int i = 0; i < nrStations; i++) {
            for (int j = 0; j < nrJobTypes; j++) {          //vond Broos zijn methode te groot + omslachtig
                sigma[i][j]=Math.sqrt(var[i][j]);
            }
        }
//        sigma[0][0] = Math.sqrt(var[0][0]);               //Processing stdev per ws and job type (WS1, J1)
//        sigma[0][1] = Math.sqrt(var[0][1]);
//        sigma[0][2] = Math.sqrt(var[0][2]);
//        sigma[0][3] = Math.sqrt(var[0][3]);
//        sigma[1][0] = Math.sqrt(var[1][0]);               //Processing stdev per ws and job type (WS2, J1)
//        sigma[1][1] = Math.sqrt(var[1][1]);
//        sigma[1][2] = Math.sqrt(var[1][2]);
//        sigma[1][3] = Math.sqrt(var[1][3]);
//        sigma[2][0] = Math.sqrt(var[2][0]);               //Processing stdev per ws and job type (WS3, J1)
//        sigma[2][1] = Math.sqrt(var[2][1]);
//        sigma[2][2] = Math.sqrt(var[2][2]);
//        sigma[2][3] = Math.sqrt(var[2][3]);
//        sigma[3][0] = Math.sqrt(var[3][0]);               //Processing stdev per ws and job type (WS4, J1)
//        sigma[3][1] = Math.sqrt(var[3][1]);
//        sigma[3][2] = Math.sqrt(var[3][2]);
//        sigma[3][3] = Math.sqrt(var[3][3]);
//        sigma[4][0] = Math.sqrt(var[4][0]);               //Processing stdev per ws and job type (WS5, J1)
//        sigma[4][1] = Math.sqrt(var[4][1]);
//        sigma[4][2] = Math.sqrt(var[4][2]);
//        sigma[4][3] = Math.sqrt(var[4][3]);


        /* STOP CRITERION */
        N = 1000;                                // Number of scans

        infinity = 999999999;


        /* 3. INITIALISE SYSTEM */
        /************************/

        /* DETERMINE FIRST ARRIVAL + FIRST DEPARTURE */
        for (i2 = 0; i2 < maxNrStations; i2++)
        {   for (i1 = 0; i1 < maxS; i1++)
            tD[i2][i1] =infinity;          // Put all departure times for all servers to +infty (system is idle and no departures have been scheduled yet
        }

        for (i1 = 0; i1 < nrArrivalSources; i1++)
            tA[i1] = exponentialDistribution(lambda[i1]);                     // Generate first arrival for all sources
        indexArr = 0;                                                                  // Initialise arrival source indicator
        firstTA = infinity;
        for(i1 = 0; i1 < nrArrivalSources; i1++)                             // Get next arrival = Smallest arrival time
        {   //printf("%lf\t", tA[i1]);
            if (firstTA > tA[i1])
            {   firstTA = tA[i1];
                indexArr = i1;
            }

        }
        //printf("\n");


        totLambda[run][indexArr] = firstTA;                                  // Add interarrival time to the counter for calculating the average interarrival time


    }

    public void radiology_system(){
        while (nD < N)                                                         // Perform simulation until prespecified number of customers have departed
        {
            firstTD = infinity;                                                //Identify next departure event = The customer that leaves the soonest (=minimum departure time)
            for (i2 = 0; i2 < nrStations; i2++)
            {   for(i1 = 0; i1 < nrServersPerStation[i2]; i1++)      // Loop over all servers and workstations
            {
                if (firstTD > tD[i2][i1])
                {   firstTD = tD[i2][i1];
                    indexDepServer = i1;
                    indexDepStation = i2;

                }
            }
            }

            //printf("First departure time %lf\n", firstTD);                // You may want to print the first departure time

            indexArr = 0;                                                      // Identify next arrival event
            firstTA = infinity;
            for(i1 = 0; i1 < nrArrivalSources; i1++)
            {   if (firstTA > tA[i1])
            {   firstTA = tA[i1];
                indexArr = i1;
            }
            }



            if (tA[indexArr] < tD[indexDepStation][indexDepServer])                 // Is first event an arrival or a departure?
            {   arrivalEvent();                                                                                // Arrival event
            }
            else                                                                                                        // Departure event
            {
                /* UPDATE STATISTICS FOR PERIOD [t,tD]*/
                totN[run] += (tD[indexDepStation][indexDepServer] - t) * n;                                     // Update statistics: continuous time statistics to count the average number of jobs in process for specific workstation and server
                for (i2 = 0; i2 < nrStations; i2++)
                    totNWS[run][i2] += (tD[indexDepStation][indexDepServer] - t) * nWS[i2];         // Update statistics: continuous time statistics to count the average number of jobs in process for specific workstation


                for (i2 = 0; i2 < nrStations; i2++)
                {   if (nWS[i2] >= nrServersPerStation[i2])
                {    totNQueueWS[run][i2] += (tD[indexDepStation][indexDepServer] - t) * (nWS[i2] - nrServersPerStation[i2]); totNQueue[run] += (tD[indexDepStation][indexDepServer] - t) * (nWS[i2] - nrServersPerStation[i2]);
                }
                }// Update statistics: continuous time statistics to count the average number of jobs in queue for specific workstation and server

                // Calculate idle time servers
                for (i2 = 0; i2 < nrStations; i2++)
                {   for (i1 = 0; i1 < nrServersPerStation[i2]; i1++)
                {   if (tD[i2][i1] == infinity)
                {   idle[run][i2][i1] += (tD[indexDepStation][indexDepServer] - t);
                }
                }
                }   // If no departure is planned, the workstation/server is idle (continuous time statistic)

                /* Increment simulation time t = tD*/
                t = tD[indexDepStation][indexDepServer];

                timeDepartureWS[run][indexDepStation][currentCust[indexDepStation][indexDepServer]] = t;                    // Store departure time customer nD in source center
                timeSystemJobWS[run][indexDepStation][currentCust[indexDepStation][indexDepServer]] += (t-timeArrivalWS[run][indexDepStation][currentCust[indexDepStation][indexDepServer]]);             // Calculate system time in source center of served customer

                nDWS[indexDepStation]++;    // Count number of departures at a particular station
                nWS[indexDepStation]--;      // Count number of scans at the station

                currentStation[currentCust[indexDepStation][indexDepServer]]++;                                               //Selected scan should go to next station


                for (i3 = 0; i3 < nWS[indexDepStation]+1; i3++)                                                                  // Selected customer is identified from scan list current station
                {   if (listScan[indexDepStation][i3] == currentCust[indexDepStation][indexDepServer])
                    break;
                }


                listScan[indexDepStation][i3] = -1;      // Selected customer is removed from scan list current station

                for (i2 = i3; i2 < nWS[indexDepStation]; i2++)                                                                   // Scan list should be updated
                    listScan[indexDepStation][i2] = listScan[indexDepStation][i2+1];  // Move all customers one position up
                listScan[indexDepStation][i2] = -1;
                //TODO fix deze lijn: printf("Departure event\tJob %d\tWorkstation %d\tServer %d\tTime %lf\tCurrent station %d of %d\tReal %d\n",currentCust[indexDepStation][indexDepServer], indexDepStation, indexDepServer, t, currentStation[currentCust[indexDepStation][indexDepServer]],nrWorkstationsPerJobType[scanType[currentCust[indexDepStation][indexDepServer]]], route[scanType[currentCust[indexDepStation][indexDepServer]]][currentStation[currentCust[indexDepStation][indexDepServer]]-1]);



                // Does the finished job has more tasks to be done on if route?
                // Evaluate number of tasks done versus required number

                if (currentStation[currentCust[indexDepStation][indexDepServer]] < nrWorkstationsPerJobType[scanType[currentCust[indexDepStation][indexDepServer]]])
                { // For this customer there are still jobs to do
                    nWS[route[scanType[currentCust[indexDepStation][indexDepServer]]][currentStation[currentCust[indexDepStation][indexDepServer]]]]++;   // Count number of scans at a particular workstation, defined by the route (scan type and current                                                            station)
                    nAWS[route[scanType[currentCust[indexDepStation][indexDepServer]]][currentStation[currentCust[indexDepStation][indexDepServer]]]]++;// Count the number of arrivals at a particular workstation
                    timeArrivalWS[run][route[scanType[currentCust[indexDepStation][indexDepServer]]][currentStation[currentCust[indexDepStation][indexDepServer]]]][currentCust[indexDepStation][indexDepServer]]=t;// Set the time of arrival at a particular workstation of a particular job



                    if (nWS[route[scanType[currentCust[indexDepStation][indexDepServer]]][currentStation[currentCust[indexDepStation][indexDepServer]]]] <= nrServersPerStation[route[scanType[currentCust[indexDepStation][indexDepServer]]][currentStation[currentCust[indexDepStation][indexDepServer]]]])  // Arrival to a system where one of the servers is idle?
                    {   for (i1 = 0; i1 < nrServersPerStation[route[scanType[currentCust[indexDepStation][indexDepServer]]][currentStation[currentCust[indexDepStation][indexDepServer]]]]; i1++)                      // Identify the system that was idle
                    {   if (tD[route[scanType[currentCust[indexDepStation][indexDepServer]]][currentStation[currentCust[indexDepStation][indexDepServer]]]][i1] == infinity)
                        break;
                    }

                        for (i2 = 0; i2 < nrServersPerStation[route[scanType[currentCust[indexDepStation][indexDepServer]]][currentStation[currentCust[indexDepStation][indexDepServer]]]]; i2++)
                        {   if (listScan[route[scanType[currentCust[indexDepStation][indexDepServer]]][currentStation[currentCust[indexDepStation][indexDepServer]]]][i2] == -1)
                            listScan[route[scanType[currentCust[indexDepStation][indexDepServer]]][currentStation[currentCust[indexDepStation][indexDepServer]]]][i2] = currentCust[indexDepStation][indexDepServer];
                        }// change scan list



                        //listScan[route[scanType[currentCust[indexDepStation][indexDepServer]]][currentStation[currentCust[indexDepStation][indexDepServer]]]][i1] = currentCust[indexDepStation][indexDepServer];// CHANGE

                        tMu = normalDistribution(mu[route[scanType[currentCust[indexDepStation][indexDepServer]]][currentStation[currentCust[indexDepStation][indexDepServer]]]][scanType[currentCust[indexDepStation][indexDepServer]]], sigma[route[scanType[currentCust[indexDepStation][indexDepServer]]][currentStation[currentCust[indexDepStation][indexDepServer]]]][scanType[currentCust[indexDepStation][indexDepServer]]]);                        // Derive service time of the system that was idle determine the departure time for the newly arrived scan
                        tMu = tMu/60;             // Change service time to an hourly basis
                        timeService[run][route[scanType[currentCust[indexDepStation][indexDepServer]]][currentStation[currentCust[indexDepStation][indexDepServer]]]][currentCust[indexDepStation][indexDepServer]] = tMu;       // Store service time

                        waitingTimeJobWS[run][route[scanType[currentCust[indexDepStation][indexDepServer]]][currentStation[currentCust[indexDepStation][indexDepServer]]]][currentCust[indexDepStation][indexDepServer]] = 0;       /* Tally a delay of 0 for this job*/
                        /* Make a machine for this workstation busy */
                        tD[route[scanType[currentCust[indexDepStation][indexDepServer]]][currentStation[currentCust[indexDepStation][indexDepServer]]]][i1] = t + tMu; totMu[run] += tMu;// Calculate departure time
                        currentCust[route[scanType[currentCust[indexDepStation][indexDepServer]]][currentStation[currentCust[indexDepStation][indexDepServer]]]][i1] = currentCust[indexDepStation][indexDepServer];                                     // Track the current customer being served
                        tD[indexDepStation][indexDepServer] = infinity;                                                // Set departure time of current station to infty
                    }
                    else
                    {   listScan[route[scanType[currentCust[indexDepStation][indexDepServer]]][currentStation[currentCust[indexDepStation][indexDepServer]]]][nWS[route[scanType[currentCust[indexDepStation][indexDepServer]]][currentStation[currentCust[indexDepStation][indexDepServer]]]]-1] = currentCust[indexDepStation][indexDepServer];// Add new scan to the last place n the scan list

                    }




                }
                else
                {
                    orderOut[nD] = currentCust[indexDepStation][indexDepServer];                             // Store the order in which scans leave the system

                    nD++;                                                                                          // Increment the number of orders that have left the system

                    n--;                                                                                            // Decrease the number of items in the system
                    timeDeparture[run][currentCust[indexDepStation][indexDepServer]] = t;                     // Set time of departure
                    timeSystem[run][currentCust[indexDepStation][indexDepServer]] = t-timeArrival[run][currentCust[indexDepStation][indexDepServer]];// Calculate system time
                    //TODO zelfde als hierboven, idk waarom deze printf hier staat: printf("Job %d has left\tTime %lf\n", currentCust[indexDepStation][indexDepServer], t);
                }




                // Determine next scan in departure station or idle status
                if (nWS[indexDepStation] >= nrServersPerStation[indexDepStation])
                {   // Compute next job by calculating the waiting time: Select the job with the largest waiting time
                    // Job with the largest waiting time = Next job on the scan list of workstation
                    currentCust[indexDepStation][indexDepServer] = listScan[indexDepStation][nrServersPerStation[indexDepStation]-1];

                    // Make idle server busy with next job and generate next departure event

                    tMu = normalDistribution(mu[route[scanType[currentCust[indexDepStation][indexDepServer]]][currentStation[currentCust[indexDepStation][indexDepServer]]]][scanType[currentCust[indexDepStation][indexDepServer]]], sigma[route[scanType[currentCust[indexDepStation][indexDepServer]]][currentStation[currentCust[indexDepStation][indexDepServer]]]][scanType[currentCust[indexDepStation][indexDepServer]]]);                        // For the system that was idle determine the departure time for the newly arrived scan
                    tMu = tMu/60;
                    timeService[run][route[scanType[currentCust[indexDepStation][indexDepServer]]][currentStation[currentCust[indexDepStation][indexDepServer]]]][currentCust[indexDepStation][indexDepServer]] = tMu;           // Store service time
                    /* Calculate the delay for this job and workstation*/
                    waitingTimeJobWS[run][route[scanType[currentCust[indexDepStation][indexDepServer]]][currentStation[currentCust[indexDepStation][indexDepServer]]]][currentCust[indexDepStation][indexDepServer]] = t - timeArrivalWS[run][indexDepStation][currentCust[indexDepStation][indexDepServer]];



                    /* Make a machine for this workstation busy */
                    tD[indexDepStation][indexDepServer] = t + tMu;
                    totMu[run] += tMu;

                }
                else
                {   // Make server in this station idle
                    tD[indexDepStation][indexDepServer] = infinity;

                }



            }
        }

    }
    public void arrivalEvent(){
        /* UPDATE STATISTICS FOR PERIOD [t,tA]*/
        totN[run] += (tA[indexArr] - t) * n;                                                     // Number of scans in system (Continuous time statistic)
        for (i2 = 0; i2 < nrStations; i2++)
            totNWS[run][i2] += (tA[indexArr] - t) * nWS[i2];                           // Number of scans in workstation (Continuous time statistic)

        for (i2 = 0; i2 < nrStations; i2++)
        {   if (nWS[i2] >= nrServersPerStation[i2])                                                         // Number of scans in queues
        {    totNQueueWS[run][i2] += (tA[indexArr] - t) * (nWS[i2] - nrServersPerStation[i2]);
            totNQueue[run] += (tA[indexArr] - t) * (nWS[i2] - nrServersPerStation[i2]);
        }
        }

        // Calculate idle time servers
        for (i2 = 0; i2 < nrStations; i2++)
        {   for (i1 = 0; i1 < nrServersPerStation[i2]; i1++)
        {   if (tD[i2][i1] == infinity)
        {   idle[run][i2][i1] += (tA[indexArr] - t);
        }
        }
        }


        /* Increment simulation time */
        t = tA[indexArr];


        /* Generate the job type*/
        Random rand = new Random();
        double j1 = rand.nextInt(1000); // Use a random number to determine the scan type
        j1 = (j1)/1000;
        scanType[nA] = 0;
        while (j1 > cumDistrScans[indexArr][scanType[nA]])     // Inversion method for discrete distribution
            scanType[nA]++;                                                           // Determine scan type of arrival

        /* Set task = 1 for this job */
        /* Determine the station for this job*/
        currentStation[nA] = 0;                                               // New arrival needs to go to first station


        for (i1 = 0; i1 < nWS[route[scanType[nA]][currentStation[nA]]]; i1++)
        {   if (listScan[route[scanType[nA]][currentStation[nA]]][i1] == -1)
            break;

        }
        listScan[route[scanType[nA]][currentStation[nA]]][i1] = nA;//add new arrival to the scan list
        timeArrivalWS[run][route[scanType[nA]][currentStation[nA]]][nA]=t;       // Set the time of arrival at a particular workstation of a particular job
        if (nA <= N)
        {   timeArrival[run][nA] = t;                                                 // store time of arrival
        }


        nAWS[route[scanType[nA]][currentStation[nA]]]++;                          // Count the number of arrivals at a particular workstation

        nWS[route[scanType[nA]][currentStation[nA]]]++;                        // Count number of scans at a particular workstation, defined by the route (scan type and current station)

        /* Are all machine busy on this station?*/
        if (nWS[route[scanType[nA]][currentStation[nA]]] <= nrServersPerStation[route[scanType[nA]][currentStation[nA]]])  // Arrival to a system where one of the servers was idle
        {   for (i1 = 0; i1 < nrServersPerStation[route[scanType[nA]][currentStation[nA]]]; i1++)                      // Identify the server that was idle
        {   if (tD[route[scanType[nA]][currentStation[nA]]][i1] == infinity)
            break;
        }
            tMu = normalDistribution(mu[route[scanType[nA]][currentStation[nA]]][scanType[nA]], sigma[route[scanType[nA]][currentStation[nA]]][scanType[nA]]);                                                       // For the system that was idle determine the departure time for the newly arrived scan
            tMu = tMu/60;
            timeService[run][route[scanType[nA]][currentStation[nA]]][nA] = tMu;
            /* Tally a delay of 0 for this job*/
            waitingTimeJobWS[run][route[scanType[nA]][currentStation[nA]]][nA] = 0;
            /* Make a machine for this workstation busy */
            tD[route[scanType[nA]][currentStation[nA]]][i1] = t + tMu; totMu[run] += tMu;
            currentCust[route[scanType[nA]][currentStation[nA]]][i1] = nA;        // Track the current customer being served
        }
        //TODO nog een printf lijn printf("Arrival event\tJob %d\tSource %d\tTime %lf\tScan Type %d\n", nA, indexArr, t, scanType[nA]);

        n++;                    // Increase number of jobs in system
        nA++;                  // Increase the number of arrivals

        /* Schedule next arrival*/
        tLambda = exponentialDistribution(lambda[indexArr]);
        tA[indexArr] = t + tLambda;
        totLambda[run][indexArr] += tLambda;


    }

    public void departureEvent(){

    }

    public void output(){
        //TODO Outputstream
    }

    public static void main(String[] args) {
        // write your code here
    }
}

package com.company;
import java.util.Random;

public class Personnel {

    private double j1, j2, j3;
    private int i6;

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

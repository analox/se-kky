package mytest;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;

import mytest.evaluation.real.*;
import optimization.search.CoopMA;
import optimization.search.Search;
import optimization.searchspace.Chromosome;
import optimization.searchspace.Individual;
import optimization.searchspace.Population;
import optimization.tools.Matrix;
import optimization.tools.Utils;
import runtime.ConfigContainer;
/**
 * Real-coded adaptive optimization test-bed
 * @author Le Minh Nghia, NTU-Singapore
 *
 */
public class AdaptiveEARun extends EARunTemplate{

	public static void main(String[] args) {
		// Get runtime parameters in Evolutionary Strategy configuration
		//ConfigContainer profile = new ConfigContainer("esconf.xml");
		ConfigContainer profile = new ConfigContainer("conf.xml");
		profile.parseParams();
		// Initialize runtime parameters
		int numOfRuns = profile.runtime.numRuns;
		double [] bestFounds = new double[numOfRuns];
		double [] numEvals = new double[numOfRuns];
		int [] numGens = new int[numOfRuns];
		
		OSS2Params evalFunc = new OSS2Params(0, args);
		int nDim = evalFunc.getTaskFitPotential().nparam; // OSS2 problem
		System.out.println("Number of activated variable: " + nDim);
		
		Population f;
		int [] chromDim;
		byte [] chromType;
		// Initialize to ES configuration
		chromDim = new int[2]; 
		chromDim[0] = nDim;
		chromDim[1] = nDim;
		chromType= new byte[2]; 
		chromType[0] = Chromosome.FLOAT;
		chromType[1] = Chromosome.FLOAT;

		/* Setting input/output stream */
		boolean success = (new File(profile.runtime.savedDir)).mkdir();
	    if (success) {
	        System.out.println("Directory: " + profile.runtime.savedDir + " created");
	    } 
	    /* Record Message and Error to files */
	    PrintStream outputStream= null;
		PrintStream errStream= null;
	    if (profile.runtime.writeIO2File) {
	    	String outputFilename = profile.runtime.savedDir + "/" + 
	    							profile.runtime.savedDir + ".out";
    		String errFilename = profile.runtime.savedDir + "/" +
    							 profile.runtime.savedDir + ".err";
    		try {
    			outputStream = new PrintStream(new FileOutputStream(outputFilename));
    			errStream = new PrintStream(new FileOutputStream(errFilename));
    		}
    		catch (IOException e) {
    			System.out.println("OOps");
    		}
    	
    		/* Redirect Output message to file */
    		//System.setOut(outputStream);
    		/* Redirect Error message to file */
	    	System.setErr(errStream);
	    }
	    /* Start run the search */
    	Search optimizer = null;
    	
		f = new Population(chromDim, chromType, profile.method.populationDim);
		// Get fitness function
		if (profile.problem.problemID < 0)
		{
			System.out.println("Missing problem description...");
			System.exit(0);
		}
		/*
		 * true: rotated function
		 * false: standard function
		 */
		
		if (f != null)
		{
			/*
			 * Two ways to use use fitness functions: 
			 * 1. Use default search on [0,1]^n space ->
			 * 		Create fitness function. Then normalize it to the actual range
			 * 2. Use search on specified range in f.initPopulation()
			 */
			// Boolean variable: (true) for Method 1, (false) for Method 2
			boolean normalize = false;
			/* Method 1 */
			// Normalize if the range is different from [0,1]
			//int nDim = profile.method.chromosomeDim;
			
			if (normalize) {
				double lowBound [] = new double[nDim];
				double highBound [] = new double[nDim];
				for (int i = 0; i <nDim; i++) {
					//lowBound[i] = -RuntimeConfiguration.maxPopRange;
					//highBound[i] = RuntimeConfiguration.maxPopRange;
					lowBound[i] = -(i+1)*10;
					highBound[i] = (i+1)*10;
				}
				evalFunc.normalize(lowBound, highBound);
			}
			/* If not normalize, (optional) randomly shift fitness function
			 * for debugging purpose */
			else {
				double [] s = new double[nDim];
				double shiftRange = (evalFunc.getUpBound() - evalFunc.getLowBound())/1000;
				for (int i = 0; i < nDim; i++)
				{
					s[i] = Math.random();
					if (s[i] > shiftRange)
						s[i] = shiftRange * s[i];
					
					s[i] = 0.0;
					
					System.out.print(s[i]+ " ");
				}
				System.out.println();
				evalFunc.setShiftVector(s);
//				System.out.println("Evaluated: " + 
//						evalFunc.getFitnessFunc(s));
			}
			
			// Start optimizing
			for (int i = 0; i < numOfRuns; i++)
			{
				// Reset evaluation count for different runs
				evalFunc.resetEvalCount();
				evalFunc.setWriteOption(profile.runtime.writeEvaluation2File);
				evalFunc.setOutputFile(i, profile.runtime.savedDir,
						profile.runtime.evalFilename);
				// Get optimizer
				optimizer = new CoopMA(profile, evalFunc); 
				// initialize function
				if (normalize)
					f.initPopulation();
				else {
					/* Method 2
					double lowerBound [] = {-5, -5, -0.5};
					double upperBound [] = {5, 5, 0.5};
					*/
					double lowerBound [] = new double[nDim];
					double upperBound [] = new double[nDim];
					for (int j = 0; j <nDim; j++) {
						// Using setting ranges
						if (profile.method.maxPopRange > 0) 
						{
							lowerBound[j] = -profile.method.maxPopRange;
							upperBound[j] = profile.method.maxPopRange;
						}
						// Using default function ranges
						else {
							lowerBound[j] = evalFunc.getDefaultRanges()[0];
							upperBound[j] = evalFunc.getDefaultRanges()[1];
						}
					}
					f.initPopulation(lowerBound, upperBound);
				}
				
				/* Add initial parameters in the OSS2 model */ 
				double [] initParam = evalFunc.getTaskFitPotential().param;
				double initRMS = evalFunc.getTaskFitPotential().Evaluate(initParam, 0);
				System.out.println("Initial RMS in SMA: " + initRMS);
				double [] lb = evalFunc.getTaskFitPotential().lb;
				double [] hb = evalFunc.getTaskFitPotential().ub;
				
				for (int c = 0; c < initParam.length; c++) {
					double val = 0;
					if (hb[c] > lb[c])
						val = (initParam[c] - lb[c])/(hb[c] - lb[c]);
					else val = lb[c];
					System.out.println("Index -" + c + " : " + val);
					f.getIndividualAtIndex(0).getChromAtIndex(0).setGene(c, val);
				}
				evalFunc.evaluate(f.getIndividualAtIndex(0));
				System.out.println("Initial fitness: " + f.getIndividualAtIndex(0).getFitnessValue());
				
				
				//f.initPopulation(lowerBound, upperBound, "p0.txt");
			   
		        try
		        {	// initialize the strategy parameters
		            for (int k = 0; k < profile.method.populationDim; k++)
		            {
		            	for (int j = 0; j < nDim; j++)
		            	{
		            		Individual indi = f.getIndividualAtIndex(k);
		            		Chromosome variance = indi.getChromAtIndex(1);
		            		double value = (variance.uBound[j] - variance.lBound[j])/2;
		            		variance.setGene(j, value);
		            	}
		            }
		            // do search
		        	optimizer.search(f);
		        	if (f.measure != null) {
		        		f.measure.evaluation = evalFunc.getEvalCount();
		        		if (profile.runtime.writeResults2File)
		        			f.measure.toFile(profile.runtime.savedDir + "/" + 
		        					profile.runtime.evalFilename + "-" + i + ".csv");
		        	}
		        	// store best found solution to array
		        	bestFounds[i] = f.measure.getBestFound();
		        	numGens[i] = f.measure.maxGens;
		        	numEvals[i] = evalFunc.getEvalCount();
		        	System.out.println("Number of eval: "
		        			+ evalFunc.getEvalCount());
		        }
		        catch (Exception gae)
		        {
		            System.out.println(gae.getMessage() + "\n" + gae.toString()
		            		+ "\n" + gae.getStackTrace().toString());
		        }
		        evalFunc.close();
			}
			/* Print out final results */
			System.out.println("/********* RESULTS *********/");
			System.out.println(optimizer.toString());
			System.out.println("Test function: " + evalFunc.toString());
			
			//f.toFile("finalPop.txt"); 
			System.out.println(f.getPopInfo());
			
			// Calculate performance statistics: mean & variance 
			System.out.println("Best solutions found: ");
			double [][] a = new double[3][numOfRuns];
			for (int i = 0; i < numOfRuns; i++) 
			{
//				System.out.println(bestFounds[i]);
				a[0][i] = bestFounds[i];
			}
			for (int i = 0; i < numOfRuns; i++) 
			{
//				System.out.println(numEvals[i]);
				a[1][i] = numEvals[i];
			}
			for (int i = 0; i < numOfRuns; i++) 
			{
//				System.out.println(numGens[i]);
				a[2][i] = numGens[i];
			}
			// Output best solutions found/function evaluations/generations to file  
			Matrix matA = new Matrix(a);
			Matrix Trans_MatA = Matrix.transpose(matA);
			System.out.println(Trans_MatA.toString());

			Matrix.DMat2File(profile.runtime.savedDir + "/" + 
					evalFunc.getName() + "_" + nDim + "D.csv", Trans_MatA);
			// Best solutions found			
			System.out.println("Average fitness: " + Utils.average(a[0]) 
					+ " +/- " + Utils.deviance(a[0]) 
					+ " | Median: " + Utils.median(a[0]));
			// Function evaluations
			System.out.println("Average number of evaluations: " + Utils.average(a[1]) 
					+ " +/- " + Utils.deviance(a[1]) 
					+ " | Median: " + Utils.median(a[1]));
			// Number of generations
			System.out.println("Average number of generations: " + Utils.average(a[2]) 
									+ " +/- " + Utils.deviance(a[2]) 
									+ " | Median: " + Utils.median(a[2]));
		}
		// Close stream files
		if (profile.runtime.writeIO2File) {
			outputStream.close();
			errStream.close();
		}
	}
}

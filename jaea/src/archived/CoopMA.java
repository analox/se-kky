package archived;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Date;
import java.util.LinkedList;
import java.util.Vector;

import mytest.evaluation.FitnessFunction;
import optimization.operator.individual.*;
import optimization.operator.population.Crossover;
import optimization.operator.population.ESMutation;
import optimization.operator.population.Evaluation;
import optimization.operator.population.Merging;
import optimization.operator.population.Mutation;
import optimization.operator.population.OperatorTemplate;
import optimization.operator.population.Performance;
import optimization.operator.population.Selection;
import optimization.operator.population.SmartLocalLearning;
import optimization.operator.population.SwarmMove;
import optimization.search.Search;
import optimization.searchspace.Chromosome;
import optimization.searchspace.Individual;
import optimization.searchspace.Population;
import optimization.tools.Matrix;
import optimization.tools.Utils;
import optimization.tools.database.DBEntry;
import optimization.tools.database.ILDatabase;
import optimization.tools.distribution.*;
import runtime.ConfigContainer;
import runtime.operator.CrossoverOption;
import runtime.operator.MergeOption;
import runtime.operator.MutationOption;
import runtime.operator.Operator;
import runtime.operator.SelectionOption;
import runtime.operator.SwarmOption;
/**
 * Search strategy - default designed for MINIMIZATION problem
 * @author Le Minh Nghia, NTU-Singapore
 *
 */
public class CoopMA extends Search {

	public static final double P_EXP = 1/10.0;
	/**
	 * Parametric & runtime configuration of the search 
	 */
	public ConfigContainer profile = null;
	/**
	 * TRUE to print the search trace. FALSE if otherwise
	 */
	private boolean debug = true;
	/**
	 * Variable to control recording: best_found, current_found, evaluations
	 * per generation to memory
	 */
	private boolean keepGenRecord = false;
	/**
	 * Variable to control recording the percentage of the population
	 * undergoing each profile of operator
	 */
	private boolean keepProfileDecision = true;
	/**
	 * Setting to debug mode
	 * @param val
	 */
	public void setDebugMode(boolean val) {
		debug = val;
	}
	
	/**
	 * List of operators used in the composite algorithm
	 */
	int nRep = 0, nIL = 0;
	ILDatabase database = null;
	OperatorTemplate [] reproduce = null;
	SmartLocalLearning [] teacher = null;
	
	Crossover crossover = new Crossover();
	Mutation mutation = new Mutation();	
	Evaluation evaluation = new Evaluation();
	Selection selector = new Selection();	
	Merging merger = new Merging();
	ESMutation mutationES;
	SwarmMove swarm = new SwarmMove();
	/**
	 * Constructors
	 * @param config
	 * @param fc
	 */
	public CoopMA(ConfigContainer config, FitnessFunction fc) 
	{
		super(fc);
		profile = config;
	}

	/**
	 * Start one complete EA search
	 * @param pop Initial population (haven't evaluated fitness)
	 */
	@SuppressWarnings("unchecked")
	public void search(Population pop) 
	{
		int iGen;
		int nImprove = 0;
		int nNonImprove = 0;
		int MAXCOUNT = 4;
		int LEVEL = 1; 
		int poolsizeF = 1;
        double best_so_far = Double.MAX_VALUE;
        Individual bestSol = null;
        
        System.out.println("Search by " + profile.method.searchMethod + "...");
        this.initOperators();
        // init Performance Record
        pop.measure = new Performance(profile.runtime.maxGenerations, keepGenRecord);
        
        System.out.println("EA start time: " + new Date().toString());

        iGen = -1;
        boolean stopCond = false;
        double currBest = Double.MAX_VALUE;
        double prevBest = Double.MAX_VALUE; 

        
        evalFunc.evaluate(pop);
        pop.computeRankings(); // since chromosomes are just copied, the fitness are available
        
        /* USE if DELTA_CLOSE stopping condition is used 
        double initBest = pop.getFittestIndividualsFitness();
        */
        int nDim = pop.getIndividualAtIndex(0).getChromAtIndex(0).n_Dim;
        double [][] cov = new double[nDim][nDim];
		for (int j = 0; j < nDim; j++)
			cov[j][j] = 1;
		Matrix covMat = new Matrix(cov);
    	GaussianDistFunc gaussian = new GaussianDistFunc(null, covMat);
        /* Start iterations - EDIT from here*/
    	/* get record file opened */
    	String filename = "decisions.csv";
	    PrintStream outputStream= null;
	    if (keepProfileDecision) {
	    	try {
				outputStream = new PrintStream(new FileOutputStream(filename));
			}
			catch (IOException e) {
				System.out.println("Recording file can not be created");
			}
	    }
    	// running
        while (!stopCond)
        {	    
        	iGen++;
        	if (profile.runtime.writePop2File)
        		pop.toFile(profile.runtime.savedDir + "/p"+iGen+".txt");
        	
        	prevBest = currBest;
        	currBest = pop.getFittestIndividualsFitness();
        	/*   */
        	if (this.debug) {
        		System.out.println(iGen+ " generation: " + evalFunc.getEvalCount() 
    								+ "\t" + currBest);
        	}
    		
        	if (currBest <  best_so_far) {
        		best_so_far = currBest;
        		bestSol = new Individual(pop.getFittestIndividual());
        	}
        	/* Check if progressively good  */
        	if (Utils.l(currBest,prevBest)) {
        		nImprove++;
        		nNonImprove = 0;
        		if ((nImprove >= MAXCOUNT) &&
        				(poolsizeF > LEVEL))
        		{
    				poolsizeF = poolsizeF / 2;
    				nImprove = 0;
    				System.out.println("Reduce pool size F = " + poolsizeF);
        		}
        	}
        	else {
        		nNonImprove++;
        		nImprove = 0;
        		if (nNonImprove >= MAXCOUNT) {
    				poolsizeF = poolsizeF * 2;
    				nNonImprove = 0;
    				System.out.println("Increase pool size F = " + poolsizeF);
        		}
        	}
        	
        	// record
        	pop.measure.setBestFound(best_so_far);
        	pop.measure.maxGens++;
        	if (keepGenRecord)
        	{
	        	pop.measure.bestFound[iGen] = best_so_far;
	        	pop.measure.currentBestFound[iGen] = currBest;
	        	pop.measure.evals[iGen] = evalFunc.getEvalCount();
        	}
        	
        	/*
        	 * Restart if poolsize factor = 16 (too large)
        	 */
        	if (poolsizeF >= 16) {
        		pop.reInitPopulation();
        		evalFunc.evaluate(pop);
                pop.computeRankings(); // since chromosomes are just copied, the fitness are available
                poolsizeF = 1;
                database.resetDB(); // reset database
                System.out.println("Premature convergence - Restart the search. Pool size F = " + poolsizeF);
                /* Record the distribution to file */
                if (keepProfileDecision) {
	                outputStream.print(currBest + ", ");
		        	for (int u = 0; u < nRep; u++)
	        			for (int v = 0; v < nIL; v++) {
	        				double ratio = 0;
	        				if (u < (nRep-1) || v < (nIL-1))
	        					outputStream.print(ratio + ", ");
	        				else
	        					outputStream.println(ratio);
	        			}
                }
        	}
        	else 
        	{
	        	Vector<Individual> temp = pop.individuals;
	        	int nIndivs = temp.size();
	        	
	        	/* 
	        	 * Run the pipeline/workflow 
	        	 * */
	        	// Selection
	        	selector.setPoolSize(pop.populationDim);
	        	temp = selector.doProcess(temp);
	        	
	        	/*
	        	 * Prepare database
	        	 */
	        	database.prepareListArray();
	        	/* 
	        	 * Select individual stream for each profile
	        	 */
	        	nIndivs = temp.size();
	        	/*
	        	 * Index 1: individual index
	        	 * Index 2: reproduction index
	        	 * Index 3: learning method index
	        	 * Index 4: FI / cost
	        	 */
	        	double [][][][] profileMetrices = new double[nIndivs][nRep][nIL][2];
	        	
	        	Vector<Individual> [][]list = new Vector[nRep][nIL];
	        	
	        	for (int k = 0; k < nRep; k++)
	        		for (int v = 0; v < nIL; v++)
	        			list[k][v] = new Vector<Individual>();
	        	
	        	// Calculate some global fitness improvement (efficient)
	        	// distribution function for crossover (POP_TYPE)
	        	// find the hyper-rectangular for reproduction pool
	        	double [] uBound = new double[nDim];
	        	double [] lBound = new double[nDim];
	        	for (int i = 0; i < nDim; i++) 
	        	{
	        		uBound[i] = -Double.MAX_VALUE;
	        		lBound[i] = Double.MAX_VALUE;
	        	}
	        
	        	
	        	for (int i = 0; i < nIndivs; i++) {
	        		Individual tempIndiv = temp.elementAt(i);
	        		for (int j = 0; j < nDim; j++) {
	        			double geneVal = tempIndiv.getChromAtIndex(0).getGeneAsDouble(j);
	        			if (Utils.le(uBound[j], geneVal)) {
	        				uBound[j] = geneVal;
	        			}
	        			if (Utils.ge(lBound[j], geneVal)) {
	        				lBound[j] = geneVal;
	        			}
	        		}
	        	}
	        	// form distribution
	        	UniformDistFunc uni = new UniformDistFunc(lBound, uBound);
	        	
	//        	System.out.println("Crossover Estimation:");
	//        	System.out.println(uni.toString());
	        	
	        	// measuring the likelihood
//	        	System.out.println("Crossover Estimation:");
	        	double [][] metricG = this.metricMeasure(uni, database);
	        	for (int i = 0; i < nIndivs; i++)
	        		for (int k =0; k < metricG.length; k++)
	        		{
	        			// fitness
	        			profileMetrices[i][0][k][0] = metricG[k][0];
	        			// cost
	        			profileMetrices[i][0][k][1] = metricG[k][1];
	        		}
	        	
	        	// Calculate local FI & compare
	//        	System.out.println("Mutation Estimation:");
	        	for (int i= 0; i < nIndivs; i++)
	        	{
	        		// Generate distribution
	        		// distribution function for mutation (INDIVIDUAL TYPE)
	        		double [] mean = temp.elementAt(i).getChromAtIndex(0).getDoubleArray();
	            	// customize distribution to individual x
	        		gaussian.setMean(mean);
	        		
	//            	System.out.println(gaussian.toString());
//	        		System.out.println("Mutation Estimation:");
	        		double [][] metricL = this.metricMeasure(gaussian, database);
	        		
	            	for (int k =0; k < metricL.length; k++) {
	            		// fitness
	            		profileMetrices[i][1][k][0] = metricL[k][0];
	            		// cost
	            		profileMetrices[i][1][k][1] = metricL[k][1];
	            	}
	        	}
	        	// Assign profile
	        	for (int i= 0; i < nIndivs; i++)
	        	{
	        		int iRep=-1, iIL=-1;
	        		double best = -Double.MAX_VALUE;
	        		boolean randomMode = true;
	        		// find indices with best profile's metric
	        		for (int u = 0; u < nRep; u++)
	        			for (int v = 0; v < nIL; v++)
	        			{
	        				// heuristic to choose the best profile
//OLDCODE	        				double fi = profileMetrices[i][u][v][0] - temp.elementAt(i).getFitnessValue();
	        				double fi = temp.elementAt(i).getFitnessValue() - profileMetrices[i][u][v][0];
	        				
//	        				System.out.println("(" + u + "," + v + ") " + "FI: " + fi + 
//	        						", Cost: " + profileMetrices[i][u][v][1]);
	        				
	        				// if positive improvement and relevant samples are used to estimate
	        				if ((fi >= 0) && (profileMetrices[i][u][v][1] > 0) ) {
//	        				if (profileMetrices[i][u][v][1] >0) {
	        					randomMode = false;	
	        					// if FI is positive, then calculate fitness/cost metric
	        					double test = fi/profileMetrices[i][u][v][1];
	//        					double test = fi;
		        				if (test >= best) {
		        					iRep = u;
		        					iIL = v;
		        					best = test;
		        				}
	        				}
	        			}
	        		// check if randomize is required
	        		if (randomMode) {
	        			// Choose nRep
	        			double threshold = 1.0/ nRep;
	        			double rand = Utils.getRandom(1.0) ;
	        			for (int v = 0;  v < nRep; v++)
	        				if ((v+1) * threshold > rand)  
	        				{
	        					iRep = v;
	        					break;
	        				}
	    				
	        			// Choose iIL
	        			threshold = 1.0/ nIL;
	        			rand = Utils.getRandom(1.0) ;
	        			for (int v = 0;  v < nIL; v++)
	        				if ((v+1) * threshold > rand) 
	        				{
	        					iIL = v;
	        					break;
	        				}
	        			
//	        			System.out.println("Random profile selection: iRep =" 
//	        					+ iRep + ", iIL = " + iIL);
	        		}
	        		else {
//	        			System.out.println("Best profile selection: iRep =" 
//	        					+ iRep + ", iIL = " + iIL);
	        		}
	        		if ((iRep >= 0) && (iIL >= 0))
	        			list[iRep][iIL].addElement(temp.elementAt(i));
	        	}
	        	/* Record the distribution to file */
	        	if (keepProfileDecision) {
		        	outputStream.print(currBest + ", ");
		        	for (int u = 0; u < nRep; u++)
	        			for (int v = 0; v < nIL; v++) {
	        				double ratio =((double) list[u][v].size())/nIndivs;
	        				if (u < (nRep-1) || v < (nIL-1))
	        					outputStream.print(ratio + ", ");
	        				else
	        					outputStream.println(ratio);
	        			}
	        	}
	        	/*
	        	 * Execute reproduction & Individual Learning
	        	 */
	        	// Crossover + Teacher 1/2
	        	for (int k = 0; k < nIL; k++) {
		        	if (list[0][k].size() > 1) {
		        		// Crossover
		        		crossover.setPoolSize(poolsizeF * list[0][k].size()); 
		        		list[0][k] = crossover.doProcess(temp);	
		        		
			        	// Evaluation
			        	list[0][k] = evaluation.doProcess(list[0][k]);
			        	// Individual Learning
			        	int size = (int)Math.round(list[0][k].size()*P_EXP);
//			        	size = (size >= 1) ? size : 1;
			        	
			        	teacher[k].setPoolSize(size);
			        	list[0][k] = teacher[k].doProcess(list[0][k]);
		        	}
	        	}
	           	// Mutation + Teacher 1/2
	        	for (int k = 0; k < nIL; k++) {
		        	if (list[1][k].size() > 0) {
		        		// Crossover
		        		mutation.setFactor(poolsizeF);
			        	list[1][k] = mutation.doProcess(list[1][k]);
			        	// Evaluation
			        	list[1][k] = evaluation.doProcess(list[1][k]);
			        	// Individual Learning
			        	int size = (int)Math.round(list[1][k].size()*P_EXP);
//			        	size = (size >= 1) ? size : 1;
			        	
			        	teacher[k].setPoolSize(size);
			        	list[1][k] = teacher[k].doProcess(list[1][k]);
		        	}
	        	}
	        	/*
	        	 * Copy back to main stream
	        	 */
	        	temp.clear();
	        	for (int k = 0; k < nRep; k++)
	        		for (int v = 0; v < nIL; v++)
	        			for (int s = 0; s < list[k][v].size(); s++)
	        			temp.add(list[k][v].elementAt(s));
	        	
	        	// Merging
	    		merger.setParents(pop.individuals);
	    		temp = merger.doProcess(temp);
	    
	    		// Swarm Move
	    		// swarm.nBest = new Individual(pop.getFittestIndividual());
				
	            pop.individuals = temp;
	            pop.computeRankings();
        	}
            /* Check stop condition */
            /* Stopping condition by limited number of fitness evaluations */
            if ((profile.runtime.maxEvaluation > 0) &&
            		(evalFunc.getEvalCount() > profile.runtime.maxEvaluation))
            {
            	stopCond = true;
            	System.out.println("Max evaluation reached at " + iGen + " generation...");
            	System.out.println("Stop running!");
            	break;
            }
            /* Stopping condition by fitness threshold (MINIMIZATION problem)*/
            if (best_so_far <= profile.runtime.fitThreshold)
            {
            	stopCond = true;
            	System.out.println("Below fitness threshold. Solution found...");
            	System.out.println("Stop running!");
            	break;
            }
            /* Stopping condition by limited number of generations/iterations */
        	if (iGen >= (profile.runtime.maxGenerations -1)) 
			{
    			System.out.println("Use up all iterations at " + iGen + " generation...");
				System.out.println("Stop running!");
    			stopCond = true;	
			}
            
        }
        /* close recording file */
        if (keepProfileDecision) outputStream.close();
        // computeFitnessRankings();
        if (profile.runtime.writePop2File)
    		pop.toFile(profile.runtime.savedDir + "/p"+(iGen+1)+".txt");
        Individual curBestIndividual = pop.individuals.elementAt(pop.bestFitnessChromIndex);
        
		if (curBestIndividual.getFitnessValue() < best_so_far) {
			best_so_far = curBestIndividual.getFitnessValue();
			bestSol.copyIndividual(curBestIndividual);
		}
		pop.bestpool.add(bestSol);
    	System.out.println("Best-found solutions: ");
    	for (int i = 0; i < pop.bestpool.size(); i++)
    	{
    		System.out.println(i+ " solution :" + pop.bestpool.elementAt(i).toString());
			System.out.println("Fittest value: " 
					+ pop.bestpool.elementAt(i).getFitnessValue());
    	}
        
    	// record last item
    	pop.measure.setBestFound(best_so_far);
        pop.measure.maxGens++;
    	if (this.keepGenRecord)
    	{
	        pop.measure.bestFound[iGen+1] = best_so_far;
	        pop.measure.currentBestFound[iGen+1] = curBestIndividual.getFitnessValue();
	        pop.measure.evals[iGen+1] = evalFunc.getEvalCount();
    	}
        System.out.println("MA end time: " + new Date().toString());
//        System.out.println(database.toString());
	}
	/**
	 * Metric measurement for global distribution.
	 * Return double [][] for each stream of individual learning
	 * [IL Index][0]: expected fitness can be obtained
	 * [IL Index][1]: expected learning cost
	 */
	public double [][] metricMeasure(DistributionFunc func, ILDatabase db)
	{
		double [][] val = new double[db.nMethods][2];
		
		// for each individual learning method
		for (int i = 0; i <  db.nMethods; i++) {
			LinkedList<DBEntry> list = db.getListForMethod(i);
			// check if list empty
			if (list.isEmpty()) {
				val[i][0] = Double.MAX_VALUE;
//				val[i][1] = 1;
				val[i][1] = 0;
			}
			else {
				double[][] array = db.getListArrayForMethod(i);
				// calculate weights for samples
				double [] w = WeightCal.weightCal(func, array);
				
				double sum = 0;
//				System.out.print("Weight ("+ w.length + "): ");
				for (int k = 0; k < w.length; k++) {
//					System.out.print(w[k] + ", ");
					sum += w[k];
				}
//				System.out.println(" | " + sum);
				
				// measure expect fitness obtained and cost
				// if no relevant data obtained
				if (sum <= 1E-9) {
//					System.out.println("No relevance samples");
					val[i][0] = Double.MAX_VALUE;
					val[i][1] = 0;
				}
				else { // measure expected fitness
//					System.out.println("Have relevance samples");
					val[i][0] = 0;
					val[i][1] = 0;
					for (int j= 0; j < w.length; j++)
					{
						// fitness
						val[i][0] += w[j] * list.get(j).getEndPoint().getFitnessValue();
						// cost
						val[i][1] += w[j] * list.get(j).getCost();
					}
				}
			}
		}
//		System.out.println();
		return val;
	}
	
	/**
	 * Initiate the operators in the work-flow list of the search
	 */
	@SuppressWarnings("unchecked")
	public void initOperators()
	{
		ArrayList template = (ArrayList) profile.method.operators;
		for (int i = 0; i < template.size(); i++)
		{
			Operator opt = (Operator)template.get(i);
			System.out.println(opt.Name);
			if (opt.Name.compareTo("Crossover")==0)
			{
				CrossoverOption option = (CrossoverOption) opt;
				crossover.setXType(option.crossoverType);
				crossover.setXProb(option.crossoverProb);
				crossover.setPoolSize(profile.method.populationDim);
			}
			else if (opt.Name.compareTo("Mutation")== 0)
			{
				MutationOption option = (MutationOption) opt;
				mutation.setMutProb(option.mutationProb);
				mutation.setMutType(option.mutationType);
				mutation.setMutRadius(option.mutationRadius);
			}
			else if (opt.Name.compareTo("Evaluation")==0)
			{
				evaluation.setFitnessFunc(evalFunc);
			}
			else if (opt.Name.compareTo("Selection")==0)
			{			
				SelectionOption option = (SelectionOption) opt;
				selector.setType(option.selectionType);
				if (option.poolsize < 0)
					// default setting
					selector.setPoolSize(profile.method.populationDim);
				else
					selector.setPoolSize(option.poolsize);
			}	
			else if (opt.Name.compareTo("Merging")==0)
			{
				MergeOption option = (MergeOption) opt;
				merger.setType(option.selectionType);
				if (option.poolsize < 0)
					// default setting
					merger.setPoolSize(profile.method.populationDim);
				else
					merger.setPoolSize(option.poolsize);
			}
			else if (opt.Name.compareTo("ES Mutation")==0)
			{
				mutationES = new ESMutation(profile.method.chromosomeDim);
			}
			else if (opt.Name.compareTo("SwarmMove")==0)
			{
				
				SwarmOption option = (SwarmOption) opt;
				if (option.inertia > 0)
					swarm.setInertia(option.inertia);
				if (option.neighborLearnRate > 0)
					swarm.setNeighborLearnRate(option.neighborLearnRate);
				if (option.localLearnRate > 0)
					swarm.setLocalLearnRate(option.localLearnRate);
				
				// Set up prelimChrom used by PSO
				swarm.pBest.clear();
				swarm.velocity.clear();
				for (int k = 0; k < profile.method.populationDim; k++) {
		    		Vector<Chromosome> temp = new Vector<Chromosome>();
		    		temp.addElement(new Chromosome(Chromosome.FLOAT, profile.method.chromosomeDim, false));
		    		
		    		Individual tempInd = new Individual(temp);
		    		tempInd.setFitnessValue(Double.MAX_VALUE);
		    		swarm.pBest.addElement(tempInd);
		    		Individual velocityVector = new Individual(temp);
		    		swarm.velocity.addElement(velocityVector);
		    	}
			}
		}
		/*
		 * 2 genetic variations
		 */
		nRep = 2;
		
		/* 3 individual learnings */		
//		nIL = 3;
		 
		/* 3 individual learnings */ 
		nIL = 2;
		
		database = new ILDatabase(nIL);
		reproduce = new OperatorTemplate[nRep];
		teacher = new SmartLocalLearning[nIL];
		/* Initialize Reproduction (TWO)*/
		this.reproduce[0] = this.crossover;
		this.reproduce[1] = this.mutation;
	
		/* Initialize Local Learning */
		int learnPoolsize = 5;
		// normal run int intensity = 300;
		int intensity = 150; // for Step function
		double accuracy = 1E-8;
		/* Set step size */
		double stepRange = (evalFunc.getUpBound() - evalFunc.getLowBound())/10;
		double stepsize = ((stepRange > 1.0) ? 1.0 : stepRange);
		System.out.println("Step size: " + stepsize);
		// set local learning
		IndivSearch [] iSearch = new IndivSearch[nIL];
		/* 3 cases */
//		iSearch[0] = new DFP(evalFunc,intensity);
//		iSearch[1] = new DSCG(evalFunc, intensity);
//		iSearch[2] = new ESIL(evalFunc, intensity);
		/* 2 cases */
		iSearch[0] = new DFP(evalFunc,intensity);
		iSearch[1] = new ESIL(evalFunc, intensity);
		
		for (int i = 0; i < nIL; i++) {
			teacher[i] = new SmartLocalLearning(database,i);
			// set poolsize
			teacher[i].setPoolSize(learnPoolsize);
			// set selection
			teacher[i].setSelectionType(Selection.KBest);
			iSearch[i].ACC = accuracy;
			iSearch[i].StepSize = stepsize;
			teacher[i].setLearningOpt(iSearch[i]);
		}
	}
	
	public String toString() 
	{
		String res = "";
		res = profile.method.toString();
		return res;
	}
}

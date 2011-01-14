package optimization.operator.individual;

import java.util.Random;

import mytest.EARunTemplate;
import mytest.evaluation.FitnessFunction;
import optimization.searchspace.Individual;
import optimization.tools.Matrix;

/**
 * Individual Learning in ES style
 * @author  Le Minh Nghia, NTU-Singapore
 *
 */
public class ESIL extends IndivSearch {
	Random rd = null;
	/**
	 * Enable debug = true to print out search trace
	 */
	public boolean debug = false;
	
	public ESIL(FitnessFunction fc, int eval) {
		super(fc, eval);
		rd = new Random(System.currentTimeMillis());
		// TODO Auto-generated constructor stub
	}

	public ESIL(FitnessFunction fc) {
		super(fc);
		rd = new Random(System.currentTimeMillis());
		// TODO Auto-generated constructor stub
	}

	@Override
	public void search(Individual indiv) {
		// TODO Auto-generated method stub
		double [] xInit = indiv.getChromAtIndex(0).getDoubleArray();
		double yInit = indiv.getFitnessValue();
		double fit = search(xInit, yInit);
		for (int i = 0; i < xInit.length; i++)
			indiv.getChromAtIndex(0).setGene(i, xInit[i]);
		indiv.setFitnessValue(fit); 
	}
	/**
	 * (1,nDim) ES with 1/5 rules
	 */
	public double search(double[] xInit, double yInit)
	{
		double [] fret = new double[1];
		double [] step = new double[1];
		if (!this.StopIfConverge) 
		{
			if (debug) {
				System.out.println("ESIL at budget: " + this.evalMax + " evaluations");
			}
			step[0] = this.StepSize;
			fret[0] = this.move(xInit, yInit, this.evalMax, step);
		}
		else 
		{
			if (debug) {
				System.out.println("ESIL runs until convergence");
			}
			
			/* Check converge by phenotype - OLD convergence */
			double [] xTemp = new double[xInit.length];
			int count = 0;
			step[0] = this.StepSize;
			
			do {
				// Copy xInit to xTemp
				for (int i = 0; i < xInit.length; i++)
					xTemp[i] = xInit[i];
				// Search from xInit
				fret[0] = move(xInit, yInit, xInit.length, step);
				// After search, xInit is updated. Now update yInit
				yInit = fret[0];
				// Calculate the length of this move
				double distance = distance(xInit, xTemp);
				
				if (distance >= this.CONV_ACC)
				{
					if (debug)
						System.out.println("Improved distance: " + distance);
					count = 0;
				}
				else
				{
					if (debug)
						System.out.println("Converged distance: " + distance);
					count++;
				}
			}
			while (count < this.MAX_COUNT);
			
		}
		return fret[0];
	}
	/**
	 * Function to calculate Euclidean distance
	 * @param x 
	 * @param y
	 * @return
	 */
	private double distance(double [] x, double [] y)
	{
		 double res = 0;
		 if (x.length != y.length) return Double.MIN_VALUE;
		 for (int i = 0;  i < x.length; i++)
		 {
			 double temp = x[i] - y[i]; 
			 res += temp * temp;
		 }
		 return Math.sqrt(res);
	}
	public double move(double[] xInit, double yInit, int maxEval, double [] step)
	{
		double fret = yInit;
		long startEval = evalFunc.getEvalCount();
		int nDim = xInit.length;
		while((evalFunc.getEvalCount() < startEval + maxEval) &&
				(step[0] >= this.ACC))
		{
			// copy x
			double [] x_ = new double[nDim];
			for (int i = 0; i < nDim; i++) x_[i] = xInit[i];
			
			int nSuccess = 0;
			// generate nDim x offspring and check
			for (int i = 0; i < nDim; i++)
			{
				// create offspring by Gaussian distribution
				double [] offspring = new double[nDim];
				for (int j = 0; j < nDim; j++)
					offspring[j] = x_[i] + step[0]*rd.nextGaussian();
				// evaluate
				double offVal = evalFunc.getFitnessFunc(offspring);
				// if got improve offspring
				if (offVal < yInit) {
					nSuccess++;
					// compare with current best
					if (offVal < fret) {
						// copy new value to xInit
						for (int j = 0; j < nDim; j++)
							xInit[j] = offspring[j];
						fret = offVal;
					}
				}
			}
			yInit = fret;
			double ratio = 0.5;
			if (nSuccess/((double)nDim) >= 0.2)
				step[0] = step[0]/ratio;
			else
				step[0] = step[0]*ratio;
			if (debug) {
				System.out.println("Success rate = " + nSuccess/((double)nDim) 
						+ " - New step = " + step[0]);
			}
		}
		return fret;
	}
	static public void main(String [] args)
	{
//		int ndim = 30;
//		double [] x = new double[ndim];
		
		double [] x = {1, 1,0.9, 0.3};
		FitnessFunction f = 
			EARunTemplate.getFitnessFunc(FitnessFunction.F_Sphere, x.length, false);
//		double a = f.getUpBound();
		double a = 2;
		
		double [] v = {-a, a};
		f.setDefaultRanges(v);
		
//		int [] index = new int[ndim];
//		for (index[0] = 0; index[0] < 8; index[0]++)
//			for (index[1] = 0; index[1] < 8; index[1]++)
//				for (index[2] = 0; index[2] < 4; index[2]++)
//					for (index[3] = 0; index[3] < 4; index[3]++)
//						for (index[4] = 0; index[4] < 4; index[4]++)
			{
//				for (int j =0; j < x.length; j++)
//						x[j] = -a + 2*a* Math.random();
					double y = f.getFitnessFunc(x);
					System.out.print(y);
					System.out.println("Initial: \n" + (new Matrix(x)).toString() + "->" + y);
		
					ESIL directSearch = new ESIL(f, 300);
					directSearch.ACC = 1E-4;
					directSearch.StepSize = 200.0;
					directSearch.debug = true;
					
					y = directSearch.search(x, y);
					System.out.println("Final: \n" + (new Matrix(x)).toString() + "->" + y + "| " + f.getFitnessFunc(x));
//					System.out.println(", " + y);
					System.out.println("Eval: " + (f.getEvalCount()-1));
			}
	}
}

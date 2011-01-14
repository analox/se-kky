package optimization.tools.distribution;

import optimization.tools.Matrix;

/**
 * Provide weights to samples according to given distribution
 * @author lemi0005
 *
 */
public class WeightCal {
	static public final double ESP = Double.MIN_VALUE;
	static public final double STEP = 1E-4;
	/**
	 * 
	 * @param p
	 * @param array
	 * @return Vector of weighted RELEVANT samples, otherwise vector of 0s for irrelevant samples
	 */
	static public double [] weightCal(DistributionFunc p, double[][] array)
	{
		if (array == null) return null;
		int nSamples = array.length;
		int nDim = array[0].length;
		double [] w = new double [nSamples];
		double sum = 0;
		
		// for each samples
		for (int i = 0; i < nSamples; i++)
		{
			// initialize 
			w[i] = 0;
			
			// copy values to temp
			double [] temp = new double [nDim];
			for (int k = 0; k < nDim; k++)
				temp[k] = array[i][k];
			
			// for each dimension
			for (int j = 0; j < nDim; j++)
			{		
//				double oldVal;
//				oldVal = temp[j];
				// move right
//				temp[j] = oldVal + STEP;
				w[i] += p.getDensityVal(temp);
				// move left
//				temp[j] = oldVal - STEP;
//				w[i] += p.getDensityVal(temp);
				
//				temp[j] = oldVal;
			}
			// average the weight
//			w[i] = w[i]/(2*nDim);
			// for normalization
			sum += w[i];
		}
		if (sum < nDim*ESP) {
			for (int i = 0; i < nSamples; i++) w[i] = 0;
		}
		else {
			for (int i = 0; i < nSamples; i++) w[i] = w[i]/sum;
		}
		return w;
	}
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// distribution code
		double [] mean = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
		double [] [] sigma = {{3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
							  {0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
							  {0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0},
							  {0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0},
							  {0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0},
							  {0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0},
							  {0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0},
							  {0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0},
							  {0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0},
							  {0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0},
							  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0},
							  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3},
							};
		Matrix sigMat = new Matrix(sigma);
		GaussianDistFunc gauss = new GaussianDistFunc(mean, sigMat);
		// uniform 
		double [] l = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};
		double [] u = {3, 4, 3, 3, 4, 3, 3, 4, 3, 3, 4, 3};
		UniformDistFunc uni = new UniformDistFunc(l, u);
		// input
		double [][] array = {
				{1.01, 1.02, 1.01, 1.02, 1.01, 1.02, 1.01, 1.02, 1.01, 1.02, 1.01, 1.02},
				{2, 3, 2, 3, 2, 3, 2, 3, 2, 3, 2, 3}, 
				{1, 3, 1, 3, 1, 3, 1, 3, 1, 3, 1, 3}, 
				{1.5, 2.5, 1.5, 2.5, 1.5, 2.5, 1.5, 2.5, 1.5, 2.5, 1.5, 2.5}
				};
		double [] fitness = {3, 1, 1, 2};
		
		// weight
		double [] weight;
		double sum;
		// calculate weight
		System.out.println("Gaussian distribution: ");
		weight = WeightCal.weightCal(gauss, array);
		sum = 0;
		for (int i = 0; i < weight.length; i++) {
			System.out.print(weight[i] + ", ");
			sum += weight[i] * fitness[i];
		}
		System.out.println();
		System.out.println("Average fitness = " + sum);
		
		System.out.println("Uniform distribution: ");
		weight = WeightCal.weightCal(uni, array);
		sum = 0;
		for (int i = 0; i < weight.length; i++) {
			System.out.print(weight[i] + ", ");
			sum += weight[i] * fitness[i];
		}
		System.out.println();
		System.out.println("Average fitness = " + sum);
		
		
	}

}

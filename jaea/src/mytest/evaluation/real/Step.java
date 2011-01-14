package mytest.evaluation.real;
import mytest.evaluation.*;
/**
 * Implementation of Step benchmark function. <p></p>
 * Step - Global opt : 0.0 at x[] = 0
 * @author Le Minh Nghia, NTU-Singapore
 *
 */
public class Step extends FitnessFunction{
	public String getName() {return "Step";}
	public Step(long initCount) {
		super(initCount);
		this.defaultRanges[0] = -5.12;
		this.defaultRanges[1] = 5.12;
	}
	
	public double calculate(double [] inputs) 
	{
		double res = 0;

		for(int i = 0; i < inputs.length; i++ ) {
			res += Math.floor(inputs[i]);
		}

		return res + 6*inputs.length;
	}
	/**
	 * Override function - Using analytical gradient function instead of approximaton
	 */
	public void getFitnessGrad(double [] inputs, double f, double [] g)
	{
		/* Exact gradient formula */
		for (int i = 0; i < inputs.length; i++)
			g[i] = 0.0;
		/* Increase number of evaluations */
		NUM_EVAL += inputs.length; 
	}
	
}

package mytest.evaluation.real;
import java.util.Random;

import mytest.evaluation.FitnessFunction;
/**
 * Implementation of Schwefel benchmark function. <p></p>
 * Schwefel's problem 1.2 - Unimodal
 * Global optimum: f = 0 at x[] = 0
 * @author Le Minh Nghia, NTU-Singapore
 */
public class SchwefelNoise extends FitnessFunction {
	public String getName() {return "SchwefelNoise";}
	public SchwefelNoise(long initCount) {
		super(initCount);
		this.defaultRanges[0] = -100;
		this.defaultRanges[1] = 100;
	}

	public double calculate(double[] x) {

		double prev_sum, curr_sum, outer_sum;

		curr_sum = x[0];
		outer_sum = (curr_sum * curr_sum);

		for (int i = 1 ; i < x.length ; i ++) {
			prev_sum = curr_sum;
			curr_sum = prev_sum + x[i];
			outer_sum += (curr_sum * curr_sum);
		}
		
		/* Noise generator */
		Random gene = new Random();
		double rand = gene.nextGaussian();
		outer_sum = outer_sum * ( 1 + 0.4* Math.abs(rand));

		return (outer_sum);
	}
}

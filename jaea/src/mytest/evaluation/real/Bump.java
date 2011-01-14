package mytest.evaluation.real;
import mytest.evaluation.*;
public class Bump extends FitnessFunction {
	public String getName() {return "Bump";}
	public Bump(long initCount) {
		super(initCount);
		this.defaultRanges[0] = 0;
		this.defaultRanges[1] = 10;
	}
	
	public double calculate(double [] inputs) 
	{
		double res = 0;
		double sum_cosX, prod_cosX, sum_X_2, prod_X, sum_X;
		sum_cosX = 0;
		sum_X_2 = 0; 
		sum_X = 0;
		prod_X = 1;
		prod_cosX = 1;
		
		for(int i = 0; i < inputs.length; i++ ) {
			//calculate each term
			double temp = Math.pow(Math.cos(inputs[i]), 2);
			
			sum_cosX = sum_cosX + temp * temp;
			prod_cosX = prod_cosX * temp;
			sum_X_2 = sum_X_2 + (i+1) * inputs[i] * inputs[i];
			prod_X = prod_X * inputs[i];
			sum_X = sum_X + inputs[i];
		}
		
		if ((prod_X > 0.75) && (sum_X < 7.5*inputs.length))
			res = - Math.abs(sum_cosX - 2*prod_cosX)/Math.sqrt(sum_X_2); 
		else
			res = 0;
		return res;
	}
	public static void main(String [] args) {
	    double [] inputs = {6.283185309654639, 3.106113581516227, 3.0711887359410572, 3.0265337603665823, 3.0010159426160987, 2.9491757724802428, 2.9118234898106965, 1.5172059382470031, 0.6796812756864266, 0.33480708967672684, 0.36265212997447266, 0.6501123355299889, 0.4465067413985092, 0.6237423478117559, 0.41925103084249377, 0.33586065269823523, 0.5640270589002103, 0.3228923562046965, 0.5700749201058536, 0.4840600720652139};
	    Bump a = new Bump(0);
	    System.out.println(a.calculate(inputs));
	}
}

package optimization.tools.distribution;
/**
 * 
 * @author Le Minh Nghia, NTU-Singapore
 *
 */
public class UniformDistFunc extends DistributionFunc{

	public double [] lBound = null;
	public double [] uBound = null;
	private int nDim = 0;
	private double hyperVol = 1;
	private double dValue = 1;
	
	/**
	 * 
	 * @param l Lower bound
	 * @param u Upper bound
	 */
	public UniformDistFunc(double [] l, double [] u)
	{
		setBounds(l,u);
	}
	
	public void setBounds(double [] l, double [] u)
	{
		if (l.length != u.length)
		{
			System.err.println("Dimension mismatched...");
			return;
		}
		this.lBound = l;
		this.uBound = u;
		this.nDim = lBound.length;
		hyperVol = 1;
		for (int i = 0; i < nDim; i++)
		{
			double temp = uBound[i] - lBound[i];
			if (temp > 0) hyperVol *= temp;
			if (temp < 0)
			{
				System.err.println("Dimension mismatched...");
				return;
			}	
		}
		if (hyperVol != 0) dValue = 1/hyperVol;
	}
	@Override
	public double getDensityVal(double[] x) {
		// TODO Auto-generated method stub
		boolean inRange = true;
		for (int i = 0; i < nDim; i++)
			if ((x[i] < lBound[i]) || (x[i] > uBound[i]))
				inRange = false;
		if (!inRange)
			return 0;
		return dValue;
	}
	public String toString() {
		String res = "Uniform distribution - \n";
		for (int i = 0; i < lBound.length; i++) {
			res += "[" + lBound[i] +"," + uBound[i] + "]x";
		}
		return res;
	}
	static public void main(String [] args)
	{
		double [] l = {1, 1};
		double [] u = {1, 3};
		double [] x = {1, 1.02};
		UniformDistFunc func = new UniformDistFunc(l, u);
		System.out.println("PDF: " + func.getDensityVal(x));
	}
}

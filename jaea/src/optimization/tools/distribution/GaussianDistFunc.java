package optimization.tools.distribution;
import optimization.tools.*;
/**
 * Gaussian distribution
 * @author Le Minh Nghia, NTU-Singapore
 *
 */
public class GaussianDistFunc extends DistributionFunc {

	public double [] mean;
	public Matrix covMat = null;
	public Matrix inverseCovMat = null;
	public double rsDet = 0;
	static public final double rsPi = Math.sqrt(2*Math.PI);
	/**
	 * 
	 * @param m Mean
	 * @param cov Covariance matrix
	 */
	public GaussianDistFunc(double [] m, Matrix cov)
	{
		this.mean = m;
		this.covMat = cov;
		this.inverseCovMat = Matrix.inverse(covMat);
		this.rsDet = Math.sqrt(covMat.getDeterminant());
	}
	public void setMean(double [] a)
	{
		this.mean = a;
	}
	public double [] getMean()
	{
		return this.mean;
	}
	public Matrix getCovMat() {return this.covMat; }
	public void setCovMat(Matrix mat)
	{
		this.covMat = mat;
		this.inverseCovMat = Matrix.inverse(covMat);
		this.rsDet = Math.sqrt(covMat.getDeterminant());
	}
	@Override
	public double getDensityVal(double[] x) {
		// TODO Auto-generated method stub
		double res = 0;
		if (x.length != mean.length) return Double.NaN;
		double [] dist2Mean = new double [x.length];
		for (int i = 0; i < x.length; i++)
			dist2Mean[i] = x[i] - mean[i];
		Matrix d2Mat = new Matrix(dist2Mean);
		Matrix d2MatTran = Matrix.transpose(d2Mat);
		Matrix temp = Matrix.multiply(d2Mat, 
				Matrix.multiply(inverseCovMat, d2MatTran));
		
		double kernelVal = temp.mat[0][0];
		double upper = Math.exp(-kernelVal/2);
		res = upper;
		for (int i = 0; i < x.length; i++)
			res = res/ rsPi;
		res = res/rsDet;
		return res;
	}
	public String toString() {
		String res = "Gassian distribution - \n";
		res += "[Mean]: " ;
		for (int i = 0; i < mean.length; i++) res+= mean[i] + ",";
		res += "\n[Cov]: \n";
		res += this.covMat.toString();	
		return res;
	}
	static public void main(String [] args)
	{
		double [] mean = {1, 1};
		double [] [] sigma = {{0.25 ,0.3},{0.3 , 1}};
		Matrix sigMat = new Matrix(sigma);
		double [] x = {1.01, 1.02};
		GaussianDistFunc func = new GaussianDistFunc(mean, sigMat);
		System.out.println("PDF: " + func.getDensityVal(x));
	}
}

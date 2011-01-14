package mytest;

import mytest.evaluation.FitnessFunction;
import mytest.evaluation.real.*;
import optimization.search.*;
import optimization.tools.Matrix;
import runtime.ConfigContainer;

/**
 * Abstract class contains basic functions - to be extended by MAIN class/optimizer 
 * @author Le Minh Nghia, NTU-Singapore
 *
 */
public abstract class EARunTemplate {
	
	/**
	 * Get optimization algorithm instance
	 * @param c Container of runtime and algorithm parameters
	 * @param evalFunc Fitness function to be optimized
	 * @return search algorithm in Search template
	 */
	static public Search getOptimizer(ConfigContainer c, FitnessFunction evalFunc) {
		Search optimizer = null;
		optimizer = new SimpleEA(c, evalFunc);	
		return optimizer;
	}
	
	/**
	 * Get fitness function instance 
	 * @param fType Function code (FitnessFunction class)
	 * @return fitness Function in FitnessFunction type
	 */
	static public FitnessFunction getFitnessFunc(byte fType) {
		FitnessFunction res = null;
		switch (fType) {
		case FitnessFunction.F_Ackley:
			res = new Ackley(0);
			break;
		case FitnessFunction.F_Elliptic:
			res = new Elliptic(0);
			break;
		case FitnessFunction.F_Equal:
			res = new Equality(0);
			break;
		case FitnessFunction.F_Griewank:
			res = new Griewank(0);
			break;
		case FitnessFunction.F_MultiCos:
			res = new MultiCosine(0);
			break;
		case FitnessFunction.F_Rastrigin:
			res = new Rastrigin(0);
			break;
		case FitnessFunction.F_Rastrigin_NonCon:
			res = new RastriginNonCont(0);
			break;
		case FitnessFunction.F_Rosenbrock:
			res = new Rosenbrock(0);
			break;
		case FitnessFunction.F_Sphere_Noise:
			res = new SphereNoise(0);
			break;
		case FitnessFunction.F_Sphere:
			res = new Sphere(0);
			break;
		case FitnessFunction.F_Weierstrass:
			res = new Weierstrass(0);
			break;
		case FitnessFunction.F_Schwefel_102:
			res = new Schwefel(0);
			break;
		case FitnessFunction.F_ExpandedScaffer:
			res = new ExpandedScaffer(0);
			break;
		case FitnessFunction.F_ScaledSphere:
			res = new ScaledSphere(0);
			break;
		case FitnessFunction.F_DeceptiveCore:
			res = new DeceptiveCore(0);
			break;
		case FitnessFunction.F_DeceptiveRastrigin:
			res = new DeceptiveRastrigin(0);
			break;
		case FitnessFunction.F_SchwefelNoise:
			res = new SchwefelNoise(0);
			break;
		case FitnessFunction.F_Step:
			res = new Step(0);
			break;
		case FitnessFunction.F_Bump:
			res = new Bump(0);
			break;
		}
		return res;
	}
	/**
	 * Get fitness function instance 
	 * @param fType Function code (FitnessFunction class)
	 * @return fitness Function in FitnessFunction type
	 */
	static public FitnessFunction getFitnessFunc(byte fType, int nDim, boolean fromFile) {
		FitnessFunction res = null;
		res = getFitnessFunc(fType);
		if (fromFile) {
			switch (fType) {
			case FitnessFunction.F_Ackley:
				if (nDim == 2) {
					Matrix rotMat = Matrix.File2DMat("input_data/ackley_M_D2.txt", "\\s+");
					res.setRotationMatrix(rotMat.mat);
				}
				else if (nDim == 10) {
					Matrix rotMat = Matrix.File2DMat("input_data/ackley_M_D10.txt", "\\s+");
					res.setRotationMatrix(rotMat.mat);
				}
				else if (nDim == 30) {
					Matrix rotMat = Matrix.File2DMat("input_data/ackley_M_D30.txt", "\\s+");
					res.setRotationMatrix(rotMat.mat);
				}
				else if (nDim == 50) {
					Matrix rotMat = Matrix.File2DMat("input_data/ackley_M_D50.txt", "\\s+");
					res.setRotationMatrix(rotMat.mat);
				}
				break;
			case FitnessFunction.F_Elliptic:
				if (nDim == 2) {
					Matrix rotMat = Matrix.File2DMat("input_data/elliptic_M_D2.txt", "\\s+");
					res.setRotationMatrix(rotMat.mat);
				}
				else if (nDim == 10) {
					Matrix rotMat = Matrix.File2DMat("input_data/elliptic_M_D10.txt", "\\s+");
					res.setRotationMatrix(rotMat.mat);
				}
				else if (nDim == 30) {
					Matrix rotMat = Matrix.File2DMat("input_data/elliptic_M_D30.txt", "\\s+");
					res.setRotationMatrix(rotMat.mat);
				}
				else if (nDim == 50) {
					Matrix rotMat = Matrix.File2DMat("input_data/elliptic_M_D50.txt", "\\s+");
					res.setRotationMatrix(rotMat.mat);
				}
				break;
			case FitnessFunction.F_Equal:
				
				break;
			case FitnessFunction.F_Griewank:
				if (nDim == 2) {
					Matrix rotMat = Matrix.File2DMat("input_data/griewank_M_D2.txt", "\\s+");
					res.setRotationMatrix(rotMat.mat);
				}
				else if (nDim == 10) {
					Matrix rotMat = Matrix.File2DMat("input_data/griewank_M_D10.txt", "\\s+");
					res.setRotationMatrix(rotMat.mat);
				}
				else if (nDim == 30) {
					Matrix rotMat = Matrix.File2DMat("input_data/griewank_M_D30.txt", "\\s+");
					res.setRotationMatrix(rotMat.mat);
				}
				else if (nDim == 50) {
					Matrix rotMat = Matrix.File2DMat("input_data/griewank_M_D50.txt", "\\s+");
					res.setRotationMatrix(rotMat.mat);
				}
				break;
			case FitnessFunction.F_MultiCos:
				
				break;
			case FitnessFunction.F_Rastrigin:
				if (nDim == 2) {
					Matrix rotMat = Matrix.File2DMat("input_data/rastrigin_M_D2.txt", "\\s+");
					res.setRotationMatrix(rotMat.mat);
				}
				else if (nDim == 10) {
					Matrix rotMat = Matrix.File2DMat("input_data/rastrigin_M_D10.txt", "\\s+");
					res.setRotationMatrix(rotMat.mat);
				}
				else if (nDim == 30) {
					Matrix rotMat = Matrix.File2DMat("input_data/rastrigin_M_D30.txt", "\\s+");
					res.setRotationMatrix(rotMat.mat);
				}
				else if (nDim == 50) {
					Matrix rotMat = Matrix.File2DMat("input_data/rastrigin_M_D50.txt", "\\s+");
					res.setRotationMatrix(rotMat.mat);
				}
				break;
			case FitnessFunction.F_Rastrigin_NonCon:
				
				break;
			case FitnessFunction.F_Rosenbrock:
				
				break;
			case FitnessFunction.F_Sphere_Noise:
				
				break;
			case FitnessFunction.F_Sphere:
				
				break;
			case FitnessFunction.F_Weierstrass:
				if (nDim == 2) {
					Matrix rotMat = Matrix.File2DMat("input_data/weierstrass_M_D2.txt", "\\s+");
					res.setRotationMatrix(rotMat.mat);
				}
				else if (nDim == 10) {
					Matrix rotMat = Matrix.File2DMat("input_data/weierstrass_M_D10.txt", "\\s+");
					res.setRotationMatrix(rotMat.mat);
				}
				else if (nDim == 30) {
					Matrix rotMat = Matrix.File2DMat("input_data/weierstrass_M_D30.txt", "\\s+");
					res.setRotationMatrix(rotMat.mat);
				}
				else if (nDim == 50) {
					Matrix rotMat = Matrix.File2DMat("input_data/weierstrass_M_D50.txt", "\\s+");
					res.setRotationMatrix(rotMat.mat);
				}
				break;
			case FitnessFunction.F_Schwefel_102:
				
				break;
			case FitnessFunction.F_ExpandedScaffer:
				if (nDim == 2) {
					Matrix rotMat = Matrix.File2DMat("input_data/E_ScafferF6_M_D2.txt", "\\s+");
					res.setRotationMatrix(rotMat.mat);
				}
				else if (nDim == 10) {
					Matrix rotMat = Matrix.File2DMat("input_data/E_ScafferF6_M_D10.txt", "\\s+");
					res.setRotationMatrix(rotMat.mat);
				}
				else if (nDim == 30) {
					Matrix rotMat = Matrix.File2DMat("input_data/E_ScafferF6_M_D30.txt", "\\s+");
					res.setRotationMatrix(rotMat.mat);
				}
				else if (nDim == 50) {
					Matrix rotMat = Matrix.File2DMat("input_data/E_ScafferF6_M_D50.txt", "\\s+");
					res.setRotationMatrix(rotMat.mat);
				}
				break;
			case FitnessFunction.F_ScaledSphere:
				
				break;
			case FitnessFunction.F_DeceptiveCore:
				
				break;
			case FitnessFunction.F_DeceptiveRastrigin:
				
				break;
			case FitnessFunction.F_SchwefelNoise:
				
				break;
				
			case FitnessFunction.F_Step:	
				break;
				
			case FitnessFunction.F_Bump:	
				break;
			}
		}
		return res;
	}
}

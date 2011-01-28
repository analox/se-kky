package mytest.evaluation.real;
import java.util.Vector;

import optimization.tools.Matrix;

import mytest.evaluation.*;
import mytest.evaluation.real.oss2.TaskFitPotential;

/**
 * Implementation of Sphere benchmark function. <p></p>
 * Unimodal - Global opt : 0.0 at x[] = 0
 * @author Le Minh Nghia, NTU-Singapore
 *
 */
public class OSS2Params extends FitnessFunction{
	
	
	//[, ]
	//[1, 1]
	TaskFitPotential calc= null;
	public String getName() {return "OSS2Params";}
	public TaskFitPotential getTaskFitPotential() {return this.calc;}
	public OSS2Params(long initCount, String [] filenames) {
		super(initCount);
		
		this.defaultRanges[0] = 0;
		this.defaultRanges[1] = 1;
		
		Vector<String> INdatafile = new Vector<String>();
		Vector<Double> INdataWeight = new Vector<Double>();
		String INsFileParam = "kkydata/KKY-water.param";
		String INsFileRef = "kkydata/W.ref";
		String INsFileBound = "kkydata/KKY-water.bound";
		String INsFileTargetLM = "kkydata/W1-4opt-mp2.xyz";
		
		for (int i = 0; i < filenames.length; i++) {
			INdatafile.addElement(filenames[i]); INdataWeight.addElement(1.0);
		}
		
		calc = new TaskFitPotential("KKY",INsFileParam, INsFileRef, INsFileBound, INsFileTargetLM, INdatafile , INdataWeight );
		calc.Initialize();
	}
	
	public OSS2Params(long initCount, String [] filenames, String INsFileParam, 
			String INsFileRef, String INsFileBound) {
		super(initCount);
		
		this.defaultRanges[0] = 0;
		this.defaultRanges[1] = 1;
		
		Vector<String> INdatafile = new Vector<String>();
		Vector<Double> INdataWeight = new Vector<Double>();
		
		for (int i = 0; i < filenames.length; i++) {
			INdatafile.addElement(filenames[i]); INdataWeight.addElement(1.0);
		}
		
		calc = new TaskFitPotential("KKY",INsFileParam, INsFileRef, INsFileBound,  "",INdatafile , INdataWeight );
		calc.Initialize();
	}
	
	public double calculate(double [] inputs) 
	{
		double res = 0;
		
		double [] act_inputs = new double[inputs.length];
		for (int i = 0; i < inputs.length; i++)
			act_inputs[i] = calc.lb[i] + inputs[i] * (calc.ub[i] - calc.lb[i]);
		/*
		for (int i = 0; i < inputs.length; i++)
			System.out.print(act_inputs[i] + " ");
		System.out.println("\n");
		*/
		res = calc.Evaluate(act_inputs, 0);
		return res;
	}
	
	public static void main(String []args)
	{
		
		//OSS2Params f = new OSS2Params(0, new String[]{"kkydata/all.xyz"});
//		if (args.length < 5) {
//			System.out.println("Input: INsFileParam INsFileRef INsFileBound DataXYZ input");
//			System.exit(0);
//		}
		
		//OSS2Params f = new OSS2Params(0, new String[]{"kkydata/all-filtered.xyz"}, args[0], args[1], args[2]);
		
		OSS2Params f = new OSS2Params(0, new String[]{"kkydata/W1-l1.xyz","kkydata/W1-l2.xyz","kkydata/W1-l3.xyz",
													  "kkydata/W2-l1.xyz","kkydata/W2-l2.xyz","kkydata/W2-l3.xyz",
													  "kkydata/W3-l1.xyz","kkydata/W3-l2.xyz"});
		
		//Matrix inputMat = Matrix.File2DMat(args[4]);
		Matrix inputMat = Matrix.File2DMat("kkydata/input.txt");
		double [] inputs = inputMat.getFlatRow(0);
		/*
		double [] inputs = {4.967757955065471E-12, 0.16725654728240277, 6.876796967155116E-5, 
				9.266956697840489E-5, 4.912908936969936E-5, 0.9999095722031548, 5.295258043886205E-6, 
				0.9989805913241133, 0.999943079473415, 0.022925816409318792, 9.169520111001922E-6, 
				2.0204595886893104E-5, 0.9999945403112305, 0.9999570040183053, 0.5742331967444042, 
				0.9999461567784033, 0.1615109554497059, 0.37749922467342034, 0.999796728871809};
		*/		
		double [] act_inputs = new double[inputs.length];
		for (int i = 0; i < inputs.length; i++)
			act_inputs[i] = f.getTaskFitPotential().lb[i] + 
					inputs[i] * (f.getTaskFitPotential().ub[i] - f.getTaskFitPotential().lb[i]);
		
		double res = f.getTaskFitPotential().Evaluate(act_inputs, 1);
		//f.getTaskFitPotential().WriteParam(act_inputs, "H2O-fDSCG2-27D.param");
		System.out.println("RMS = " + res);
	}
}

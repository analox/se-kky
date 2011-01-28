/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package mytest.evaluation.real.oss2;

import java.io.*;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author nqc
 */
public class TaskFitPotential  {	
	
    //@Option(name="-b",usage="Bound file. Each line should consist of at least three columns including [upper] [lower] [fixed].[]",metaVar="STRING")
	String sFileBound="";
	
	boolean bWeighted=true;
	//@Option(name="-a",usage="Parameter file for potential if applicant.",metaVar="FILE")
    String sFileParam="";
    
  //@Option(name="-ref",usage="to calculate binding energies. It MUST be in F_XYZ format and in the following order: neutral, protonated, deprotonated. ",metaVar="FILE")
	String sFileRef="";
	
	String sFileTargetLM="";
	
//	Vector<String> datafile=new Vector<String>();
//	Vector<Double> dataWeight=new Vector<Double>();

	//@Option(name = "-i", usage = "Input file name. [xyz]", metaVar = "FILE")
	String sFileIn = "";
	String sFormatIn = "xyz";
	//@Option(name = "-o", usage = "Output file name. [xyz]", metaVar = "FILE")
	String sFileOut = "";
	String sFormatOut = "xyz";
	
	
	//@Option(name = "-v", usage = "Verbose level. [0]", metaVar = "INTEGER")
	int verbose = 0;
	
	//@Option(name = "-h", usage = "Print out the help")
	boolean isHelp = false;
	//CmdLineParser parser = null;
	BufferedWriter stdwriter = new BufferedWriter(new OutputStreamWriter(System.out));
	XmlWriter xmllog = new XmlWriter(stdwriter);

	Scanner fileIn=null;
	FileWriter fileOut=null;
	
	
	//@Option(name="-ff",usage="Potential model. [LJ}",metaVar="POTENTIAL")
    String sPotential="OSS2";

    
    protected static Logger logger=Logger.getLogger(Potential.class.getName());

	

	//@Option(name="-u",usage="Unit of energy if applicant: Hartree, kcal/mol, kJ/mol, eV, au.",metaVar="STRING")
    String sUnit="Hartree";

	//@Option(name="-med",usage="Optimization methods if applicant. METHOD must be one of DFPMIN(quasi-newton), CG (conjugate gradient)",metaVar="METHOD")
	String sMethod="DFPMIN";

	
	protected Potential pot;

	Cluster mol = new Cluster();

	public String getName() {
		return "FitPotential";
	}


//    @Override
//    public void Execute(String[] args) {
//         for(int i=0;i<args.length;i++){
//            if(args[i].contentEquals("-mpi")) bMPI=true;
//        }
//
//        if(!bMPI) super.Execute(args);
//        else{
//            String[] newargs=MPI.Init(args);
//            rank = MPI.COMM_WORLD.Rank();
//            ncpu = MPI.COMM_WORLD.Size();
//            ParseArguments(newargs);
//
//            try {
//                String sFileLog="main-t"+rank+".xml";
//                stdwriter=new BufferedWriter(new FileWriter(new File(sFileLog)));
//                xmllog = new XmlWriter(stdwriter);
//                System.out.println("Hi from <>"+rank);
//
//                xmllog.writeEntity(getName());
//                Initialize();
//                if(rank==MASTER) {//if u r MASTER, proceed
//                    Process();
//                }else{
//                    while(true){//if u r slaves, just waiting for the orders
//                        mpi.Status status=MPI.COMM_WORLD.Probe(MPI.ANY_SOURCE,MPI.ANY_TAG);
//                        if(status.tag==-1)                         break;
//                        if(status.tag==1){
//                            //MPI.COMM_WORLD.Recv(&a[0],mess.size, newtype,mess.destfrom,mess.tag);
//                            //return true;
//                        }
//                    }
//                }
//                Finalize();
//                xmllog.endEntity().close();
//                stdwriter.close();
//            } catch (IOException ex) {
//                logger.severe(ex.getMessage());
//            }
//            MPI.Finalize();
//        }
//    }


	
	public void Initialize(){
        try {
            xmllog.writeAttribute("InputFile", sFileIn).writeAttribute("FormatInput", sFormatIn);
            xmllog.writeAttribute("OutputFile", sFileOut).writeAttribute("FileBound", sFileBound);
            //xmllog.writeAttribute("RandomParam", Boolean.toString(isRandom));
            xmllog.writeAttribute("Verbose", Integer.toString(verbose));

            xmllog.writeEntity("Initialize");
                pot = MolExtra.SetupPotential(sPotential, sFileParam, sUnit,sMethod);
                
                xmllog.writeNormalText(pot.XMLInfo(1));
                xmllog.writeEntity("ReadBound");
                    ReadBound();
                xmllog.endEntity().flush();

                xmllog.writeEntity("ReadReference");
                    ReadReference();
                xmllog.endEntity().flush();

                xmllog.writeEntity("ReadData");
                    ReadData();
                    xmllog.writeEntity("TotalSize");
                    xmllog.writeNormalText(Integer.toString(data.size()));
                    xmllog.endEntity().flush();
                xmllog.endEntity().flush();

            xmllog.endEntity().flush();

            } catch (IOException ex) {
                //logger.severe(ex.getMessage());
            }
    }


//	protected void Process(){
//			xmllog.writeEntity("InitialEvaluation");
//				Evaluate(param,4);
//			xmllog.endEntity().flush();
//
//			xmllog.writeEntity("Optimization").writeAttribute("Method",sMethod );
//				long duration=System.currentTimeMillis();
//
//                nEvaluations=0;
//                if(nOpts>0)
//                    if(sMethod.contentEquals("GA"))
//                        FitUsingGeneticAlgorithm();
//                    else
//                        FitUsingApacheLevenbergMarquardt();
//
//				duration=(System.currentTimeMillis()-duration);
//
//                 xmllog.writeEntity("ConvergenceInfo");
//                    xmllog.writeAttribute("TotalRMS", Double.toString(finalRMS));
//                    //xmllog.writeAttribute("MaxInterations", Integer.toString(fit.getMaxIterations()));
//                    //xmllog.writeAttribute("MaxEvaluations", Integer.toString(fit.getMaxEvaluations()));
//                    xmllog.writeAttribute("Iterations", Integer.toString(nIterations));
//                    xmllog.writeAttribute("Evaluations", Integer.toString(nEvaluations));
//                    xmllog.writeAttribute("Duration", Double.toString(duration/1000.0));
//                xmllog.endEntity().flush();
//			xmllog.endEntity().flush();
//
//
//            if(finalParam!=null){
//                xmllog.writeEntity("FinalEvaluation");
//                    Evaluate(finalParam, 4);
//                xmllog.endEntity().flush();
//
//                xmllog.writeEntity("FinalParameters");
//                    WriteParam(finalParam);
//                xmllog.endEntity().flush();
//            }
//	}
    

    public int nparam; //!< number of non-fixed parameters
	public double[] param; //!< non-fixed parameters
	public double[] ub;//!< non-fixed upper bounds, all fixed parameters are removed
	public double[] lb;//!< non-fixed lower bounds, all fixed parameters are removed
	String[] label; //!< label of non-fixed parameters

	//int param_inp.length; //!< number of parameters (in total)
	double[] param_inp;//=new double[cParamMax]; //!< starting value of parameters
	public double[] ub_inp;//=new double[cParamMax];//!< input upper and lower bounds
	public double[] lb_inp;//=new double[cParamMax];
	int[] fixed;//=new int[cParamMax]; //!< determine wheter parameter are fixed
	String[] label_inp; //!< label of parameter

	Vector<Cluster> molRef= new Vector<Cluster>();//!< reference clusters
	Vector<Cluster> data=new Vector<Cluster>();//!< data/fitting clusters
	Vector<Cluster> targetLM=new Vector<Cluster>();//!<
	double[] weight;//!< weighted numbers of each clusters

	double[] target; //!< contain target data, in that case is binding energy +rms

	Vector<String> datafile=new Vector<String>();
	Vector<Double> dataWeight=new Vector<Double>();


    int nIterations=0;
    int nEvaluations=0;
    double[] finalParam;
    double finalRMS=0;



    private double CalcBindingEnergy(Cluster m,double ener,double[] ener_shift){
        double be=0;
        int nIons=m.getNonHydrogenNum();
        switch(m.getTotalCharge()){
            case  0:    be = ener - ener_shift[0]*(nIons); break;
            case  1:    be = ener - ener_shift[0]*(nIons-1) -     ener_shift[1]; break;
            case -1:    be = ener - ener_shift[0]*(nIons-1) -     ener_shift[2]; break;
            case  2:    be = ener - ener_shift[0]*(nIons-2) - 2.0*ener_shift[1]; break;
        }
        return be;
    }

    protected void ReadReference() throws FileNotFoundException{
        Scanner scanRef = new Scanner(new File(sFileRef));
			molRef.clear();
			while (scanRef.hasNext()) {
				Cluster tmpMol = new Cluster();
				tmpMol.Read(scanRef, "xyz");

				tmpMol.CorrectOrder();
				if(tmpMol.getNAtoms()==0) break;
				molRef.add(tmpMol);

				xmllog.writeEntity("Cluster");
				xmllog.writeAttribute("id", Integer.toString(molRef.size()-1));
				xmllog.writeAttribute("nAtoms", Integer.toString(tmpMol.getNAtoms()));
				xmllog.writeAttribute("Energy", Double.toString(tmpMol.getEnergy()));
				xmllog.endEntity().flush();
			}
		scanRef.close();
    }

    protected void ReadData() throws FileNotFoundException{
        ArrayList<ArrayList<Double> > arrayData=new ArrayList<ArrayList<Double>>();
        //reading the data
        for (int i = 0; i < datafile.size();i++) {
            ArrayList<Double> tmpData= new ArrayList<Double>();
			Scanner scanData = new Scanner(new File(datafile.get(i)));
			int count=0;

			while (scanData.hasNext()) {
				Cluster tmpMol = new Cluster();
				tmpMol.Read(scanData, "xyz");
				if(tmpMol.getNAtoms()==0) break;
				tmpMol.CorrectOrder();

				data.add(tmpMol);

                tmpData.add(tmpMol.getEnergy());
				count++;
			}
            arrayData.add(tmpData);

			xmllog.writeEntity("Data");
			xmllog.writeAttribute("File", datafile.get(i));
			xmllog.writeAttribute("Size", Integer.toString(count));
			xmllog.endEntity().flush();
			scanData.close();
		}

        //calculate target energies
        double[] ener_ref={molRef.get(0).getEnergy(),molRef.get(1).getEnergy(),molRef.get(2).getEnergy()};
		target=new double[data.size()];
		for (int i=0; i< data.size(); i++){
			Cluster m=data.get(i);
			target[i]= CalcBindingEnergy(m,m.getEnergy(),ener_ref);
		}


        //calculate weighted numbers
        weight = new double[data.size()];
        if(bWeighted){
            int i=0;
            for (ArrayList<Double> d : arrayData) {
                //finding the smallest ones
                double minE=Double.MAX_VALUE;
                double maxE=-Double.MAX_VALUE;
                for(Double k : d){
                    minE=Math.min(minE,k);
                    maxE=Math.max(maxE,k);
                }

                //System.out.println(" Min = " + minE + " max = "+maxE);
                double c=maxE+0.1*(maxE-minE)/d.size();
                double normalizationFactor=0;

                for(Double k : d){
                   normalizationFactor+=(c-k);
                }


                //scale weighted numbers
                for(Double k : d){
                    weight[i]=(c-k)/normalizationFactor*d.size()/weight.length;
                    i++;
                }

            }
        }else //if not weighted
            for(int i=0;i<weight.length;i++){
                weight[i]=1.0/weight.length;
            }
        
        
        //read target LM
        Scanner scanData = new Scanner(new File(sFileTargetLM));
        targetLM.clear();
        while (scanData.hasNext()) {
			Cluster tmpMol = new Cluster();
			tmpMol.Read(scanData, "xyz");
			if(tmpMol.getNAtoms()==0) break;
			tmpMol.CorrectOrder();
			tmpMol.CalcUSRsignature();

			targetLM.add(tmpMol);
		}
        scanData.close();

    }

    protected void ReadBound() throws FileNotFoundException, IOException{		
        param_inp=pot.getParam().clone();
        ub_inp=new double[param_inp.length];//!< input upper and lower bounds
        lb_inp=new double[param_inp.length];
        fixed=new int[param_inp.length]; //!< determine wheter parameter are fixed
        label_inp=new String[param_inp.length];

		int nfixed = 0;

        File file=new File(sFileBound);
        
        if(file.canRead()){//if bound file exists, read it
            Scanner fin = new Scanner(file);
            for(int i=0;i<param_inp.length;i++){
                String line = fin.nextLine();
                if (line.isEmpty()) break;

                StringTokenizer tokenizer;
                tokenizer = new StringTokenizer(line, "\t ,;"); //read no of atoms
                String info;
                //info = tokenizer.nextToken();
                //param_inp[i] = Double.parseDouble(info);
                info = tokenizer.nextToken();			lb_inp[i] = Double.parseDouble(info);
                info = tokenizer.nextToken();			ub_inp[i] = Double.parseDouble(info);
                info = tokenizer.nextToken();			fixed[i] = Integer.parseInt(info);
                label_inp[i] = "";
                while (tokenizer.hasMoreTokens()) {
                    label_inp[i] = tokenizer.nextToken();
                }
                //sscanf(line.c_str(),"%lf %lf %lf %d %line", &param_inp[param_inp.length], &lb_inp[param_inp.length], &ub_inp[param_inp.length], &fixed[param_inp.length]);

//                if((isRandom)&&(fixed[i]==0)){
//                    param_inp[i]=lb_inp[i]+Math.random()*(ub_inp[i]-lb_inp[i]);
//                }

                nfixed += fixed[i];
                xmllog.writeEntity("Param").writeAttribute("id", Integer.toString(i));
                xmllog.writeAttribute("Value", Double.toString(param_inp[i]));
                xmllog.writeAttribute("Lower", Double.toString(lb_inp[i]));
                xmllog.writeAttribute("Upper", Double.toString(ub_inp[i]));
                xmllog.writeAttribute("Fixed", Integer.toString(fixed[i]));
                xmllog.writeAttribute("Label", label_inp[i]);
                xmllog.endEntity().flush();
            }
            fin.close();
        }else{ //if not, predict from input value
            double fluc=0.2;
            if(!sFileBound.isEmpty()) fluc=Double.parseDouble(sFileBound);
            for(int i=0;i<param_inp.length;i++){
                lb_inp[i]=param_inp[i]-fluc*Math.abs(param_inp[i]);
                ub_inp[i]=param_inp[i]+fluc*Math.abs(param_inp[i]);
                fixed[i]=0;
                label_inp[i] = "";
            }

        }

		nparam = param_inp.length - nfixed;		

		xmllog.writeEntity("ActualParam");
		xmllog.writeAttribute("NParam", Integer.toString(nparam));
		xmllog.writeAttribute("NFixed", Integer.toString(nfixed));
		param = new double[nparam];
		ub = new double[nparam];
		lb = new double[nparam];
		label = new String[nparam];
		int t = 0;
		for (int i = 0; i < param_inp.length; i++) {
			if (fixed[i] == 0) {
				lb[t] = lb_inp[i];
				ub[t] = ub_inp[i];
				param[t] = param_inp[i];
				label[t] = label_inp[i];
				//cout<<i << ":  " << param_inp[i] << " " << lb_inp[i] << " " << ub_inp[i] << " " <<  fixed[i] << " "<<label_inp[i]<<endl;
				assert (lb[t] <= ub[t]);
				t++;
			}
		}
		xmllog.endEntity().flush();
    }

    protected void WriteParam(double[] param){      
            double[] result = new double[param_inp.length];
            ConvertParam(param, result);

            if (!sFileOut.isEmpty()) {
                  try {
                        fileOut = new FileWriter(new File(sFileOut));
                        pot.setParam(result);
                        pot.writeParam(fileOut);
                        fileOut.close();
                   } catch (IOException ex) {
                       logger.severe("Cannot write to "+sFileOut);
                    }
            }

            for(int i=0;i<result.length;i++){
					xmllog.writeEntity("Param").writeAttribute("id", Integer.toString(i));
					xmllog.writeAttribute("Value", Double.toString(result[i]));
					xmllog.writeAttribute("Lower", Double.toString(lb_inp[i]));
					xmllog.writeAttribute("Upper", Double.toString(ub_inp[i]));
					xmllog.writeAttribute("Fixed", Integer.toString(fixed[i]));
					xmllog.writeAttribute("Label", label_inp[i]);
					xmllog.endEntity().flush();
			}       
    }


	public TaskFitPotential(String INsPotential, String INsFileParam, String INsFileRef, String INsFileBound, String INsFileTargetLM,
			Vector<String> INdatafile, Vector<Double> INdataWeight )
	{
		//@Option(name="-a",usage="Parameter file for potential if applicant.",metaVar="FILE")
	    sFileParam= INsFileParam;
	    
	    sFileBound = INsFileBound;
	    
	  //@Option(name="-ref",usage="to calculate binding energies. It MUST be in F_XYZ format and in the following order: neutral, protonated, deprotonated. ",metaVar="FILE")
		sFileRef= INsFileRef;
		
		datafile= INdatafile;
		dataWeight= INdataWeight;
		
		sFileTargetLM=INsFileTargetLM;
		
		sPotential=INsPotential;
		sUnit="Hartree";
	}


    public double Evaluate(double[] p,int verbose){
		double[] y=new double[data.size()];
		return Evaluate(p,y,verbose);
	}

    

	public double Evaluate(double[] p,double[] y,int verbose){
        double[] alpha=new double[param_inp.length];
        ConvertParam(p,alpha);
        pot.setParam(alpha);

        if(verbose>=2){
                xmllog.writeEntity("Alpha");
                String s = "";
                for (int i = 0; i < alpha.length; i++) {
                    s += String.format("%f ", alpha[i]);
                }
                xmllog.writeNormalText(s);
                xmllog.endEntity().flush();
        }
        //calculate
        double[] ener_shift=new double[molRef.size()];
        double rms=0;

        
        //calculate new references
        for(int i=0; i<molRef.size(); i++){
            Cluster m=molRef.get(i);
            pot.setCluster(m);
            ener_shift[i] = pot.getEnergy(m.getCoords());
            
        }
        
        if (verbose>=3){
            xmllog.writeEntity("Reference");
            for (int i = 0; i < ener_shift.length; i++) {
                Cluster m = molRef.get(i);
                xmllog.writeEntity("Cluster");
                xmllog.writeAttribute("id", Integer.toString(i));
                xmllog.writeAttribute("nAtoms", Integer.toString(m.getNAtoms()));
                xmllog.writeAttribute("Energy", Double.toString(ener_shift[i]));
                xmllog.endEntity().flush();
            }
            xmllog.endEntity().flush();
        }
        
        //optimized the targeted local minima
        
        double[] disimilarity=new double[targetLM.size()];       
                
        if (verbose>=3){   
            xmllog.writeEntity("TargetLM");
        }
        
        double avgDis=0;
        for(int i=0;i<targetLM.size();i++){
        	Cluster m=(Cluster)targetLM.get(i).clone();
        	pot.Optimize(m);
        	m.CalcUSRsignature();
        	
        	disimilarity[i]= 1 - m.CalcUSRSimilarity(targetLM.get(i));
        	avgDis+=disimilarity[i]/targetLM.size();
        	 if (verbose>=3){   
                     xmllog.writeEntity("Cluster");
                     xmllog.writeAttribute("id", Integer.toString(i));
                     xmllog.writeAttribute("nAtoms", Integer.toString(m.getNAtoms()));
                     xmllog.writeAttribute("Energy", Double.toString(m.getEnergy()));
                     xmllog.writeAttribute("Dissimilarity", Double.toString(disimilarity[i]));
                     xmllog.writeAttribute("RMSGrad", Double.toString(pot.getRMSGrad()));
                     xmllog.endEntity().flush();
                 }
        	 

            if(Double.isNaN(m.getEnergy())) return 999999;
        }
        
        if (verbose>=3){
            xmllog.endEntity().flush();
        }

        //omp for private(m) reduction(+:rms)
        for (int i=0; i<data.size(); i++){
            Cluster m=data.get(i);
            //calculate binding energy
            pot.setCluster(m);
            double ener = pot.getEnergy(m.getCoords());

            y[i] = CalcBindingEnergy(m,ener,ener_shift);
            //y[i]*=1000;

            double dE = Math.abs(y[i] - target[i]);
            rms+=dE*dE* weight[i];

            if(verbose>=4){
                xmllog.writeEntity("Eval");
                xmllog.writeAttribute("id", Integer.toString(i));
                xmllog.writeAttribute("nAtoms", Integer.toString(m.getNAtoms()));
                xmllog.writeAttribute("Weight", Double.toString(weight[i]));
                xmllog.writeAttribute("ObservedBE", Double.toString(target[i]));
                xmllog.writeAttribute("CalcE", Double.toString(ener));
                xmllog.writeAttribute("CalcBE", Double.toString(y[i]));
                xmllog.writeAttribute("DeltaE", Double.toString(dE));
                xmllog.endEntity().flush();
            }
        }

        rms=Math.sqrt(rms);
        
        double score=rms+avgDis/50;

        if (verbose>=1){
             xmllog.writeEntity("Step").writeAttribute("id", Integer.toString(nEvaluations));
					xmllog.writeAttribute("RMS", Double.toString(rms));
					xmllog.writeAttribute("AvgDis", Double.toString(avgDis));
					xmllog.writeAttribute("Score", Double.toString(score));
		     xmllog.endEntity().flush();
//                xmllog.writeEntity("TotalRMS");
//                xmllog.writeNormalText(Double.toString(rms));
//                xmllog.endEntity().flush();
        }

        nEvaluations++;
        
        if(Double.isNaN(score)) score = 999999;
        return score;
    }

	private void ConvertParam(double[] actual,double[] p){
		int c=0;
		for(int i=0; i<param_inp.length; i++)
			if (fixed[i]==0){
				p[i] = actual[c];
				c++;
			}else p[i] = param_inp[i];
	}

	

   
}


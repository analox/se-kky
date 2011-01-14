package optimization.operator.population;

import java.util.Vector;

import optimization.searchspace.Individual;
import optimization.tools.database.DBEntry;
import optimization.tools.database.ILDatabase;
import optimization.operator.individual.*;

/**
 * Operator to perform individual learning on set/subset of individuals.<p></p>
 * (Population-based Operator)
 * @author Le Minh Nghia, NTU-Singapore
 *
 */
public class SmartLocalLearning extends OperatorTemplate{
	Selection selector = new Selection();
	IndivSearch learningOperator;
	ILDatabase db;
	int methodID;
	/**
	 * Set learning procedure used
	 * @param val
	 */
	public void setLearningOpt(IndivSearch val)
	{
		learningOperator = val;
	}
	/**
	 * Number of individuals undergo learning
	 */
	int poolsize;
	public void setPoolSize(int val)
	{
		poolsize = val;
	}
	/**
	 * Scheme to selection individuals that undergo learning
	 */
	int selectionType;
	public void setSelectionType(int val)
	{
		this.selectionType = val;
	}
	/**
	 * Constructor
	 */
	public SmartLocalLearning(ILDatabase database, int id)
	{
		super();
		this.db = database;
		this.methodID = id;
	}
	
    public Vector<Individual> doProcess(Vector<Individual> indivs){
    	Vector<Individual> res = null;
    	/* Local learning */
        double [] fits = new double[indivs.size()];
        for (int i = 0; i < indivs.size(); i++)
        {
    		fits[i] = indivs.elementAt(i).getFitnessValue();
        }
        // Select individuals for learning
        int [] selectedIndivs;
        switch (selectionType) {
        case Selection.RouletteWheel:
        	Scaling.scaling4selection(fits, Scaling.MINIMIZE);
        	selectedIndivs = selector.RWSelection(fits, poolsize);
        	break;
        case Selection.SUS:
        	Scaling.scaling4selection(fits, Scaling.MINIMIZE);
        	selectedIndivs = selector.SUSelection(fits, poolsize);
        	break;
        case Selection.Tournament:
        	//Scaling.scaling4selection(fits, Scaling.MAXIMIZE);
        	selectedIndivs = selector.TournamentSelection(fits, poolsize);
        	break;
        case Selection.KBest:
        	//Scaling.scaling4selection(fits, Scaling.MAXIMIZE);
        	selectedIndivs = selector.KBestSelection(fits, poolsize);
        	break;
        default:
        	// Copy whole population
        	selectedIndivs = new int[fits.length];
        	for (int i = 0; i < fits.length; i++)
        		selectedIndivs[i] = i;
        	break;
        }
        // Start learning 
        for (int i = 0; i < selectedIndivs.length; i++) {
        	long start= learningOperator.evalFunc.getEvalCount();
        	Individual startP = new Individual(indivs.elementAt(selectedIndivs[i]));
        	
        	learningOperator.search(indivs.elementAt(selectedIndivs[i]));
        	
        	long end = learningOperator.evalFunc.getEvalCount();
        	Individual endP = new Individual(indivs.elementAt(selectedIndivs[i]));
        	int cost = (int) (end - start);
        	double fi = startP.getFitnessValue() - endP.getFitnessValue();
        	DBEntry entri = new DBEntry(0, startP, endP, cost, fi);
        	db.add(entri, methodID);
        }
        res = indivs;
    	return res;
    }
}

package optimization.tools.database;

import java.util.LinkedList;

/**
 * Class to manage list of (starting, end point) for different 
 * individual learning
 * @author Le Minh Nghia, NTU-Singapore
 *
 */

public class ILDatabase {
	public static final int MAX_ENTRIES_PER_LIST = 500;
	public int nMethods;
	public LinkedList<DBEntry> [] DB = null;
	public double [][][] listArray = null;
	@SuppressWarnings("unchecked")
	public ILDatabase(int n)
	{
		if (n > 0) {
			nMethods = n;
			DB = new LinkedList [nMethods];
			listArray = new double[nMethods][][];
			for (int i = 0; i < nMethods; i++) 
			{
				DB[i] = new LinkedList<DBEntry>();
				listArray[i] = null;
			}
		}
		else
			System.out.println("Invalid input...");
	}
	/**
	 * Get list of entries for methodID
	 * @param methodID
	 * @return
	 */
	public LinkedList<DBEntry> getListForMethod(int methodID)
	{
		if (methodID >= nMethods) return null;
		return DB[methodID];
	}
	/**
	 * Add a record of some methodID to list
	 * @param entri
	 * @param methodID
	 */
	public void add(DBEntry entri, int methodID)
	{
		LinkedList<DBEntry> aList = this.getListForMethod(methodID);
		if (aList.size() >= MAX_ENTRIES_PER_LIST)
			aList.removeFirst();
		aList.addLast(entri);
	}
	public void prepareListArray() {
		for (int i = 0; i <  this.nMethods; i++) {
			LinkedList<DBEntry> list = this.getListForMethod(i);
			// check if list empty
			if (list.isEmpty()) {
				this.listArray[i] = null;
			}
			else {
				int nDim = list.getFirst().getStartPoint().getChromAtIndex(0).n_Dim;
				this.listArray[i] = new double[list.size()][nDim];
				// get array of start points
				for (int k = 0; k < list.size(); k++) {
					listArray[i][k] = list.get(k).getStartPoint().getChromAtIndex(0).getDoubleArray();
				}
			}
		}
	}
	/**
	 * Get list of start points for methodID
	 * @param methodID
	 * @return
	 */
	public double [][] getListArrayForMethod(int methodID) {
		if (methodID >= nMethods) return null;
		return this.listArray[methodID];
	}
	/**
	 * Reset the database
	 */
	@SuppressWarnings("unchecked")
	public void resetDB() {
		DB = new LinkedList [nMethods];
		listArray = new double[nMethods][][];
		for (int i = 0; i < nMethods; i++) 
		{
			DB[i] = new LinkedList<DBEntry>();
			listArray[i] = null;
		}
	}
	public String toString() {
		String res ="";
		for (int i = 0; i < nMethods; i++)
		{
			res += "\t List ID = " + i + ", " + DB[i].size() + " entries \n";
			for (int j = 0; j < DB[i].size(); j++)
			{
				res += ((DBEntry)DB[i].get(j)).toString();
			}
		}
		return res;
	}
}

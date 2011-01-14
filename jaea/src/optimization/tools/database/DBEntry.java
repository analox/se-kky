package optimization.tools.database;

import optimization.searchspace.Individual;

public class DBEntry {
	int t;
	Individual startPoint, endPoint;
	long cost;
	double fi;
	public DBEntry(int t, Individual s, Individual e, long c, double f) {
		this.t = t;
		this.startPoint = new Individual(s);
		this.endPoint = new Individual(e);
		this.cost = c;
		this.fi= f;
	}
	public int getGenNum() {return t;}
	public long getCost() {return cost;}
	public Individual getStartPoint() {return startPoint;}
	public Individual getEndPoint() {return endPoint;}
	public String toString() {
		String res = "=========== ";
		res += t + "\n";
		res += startPoint.toString() + "\n";
		res += endPoint.toString() + "\n";
		res += "FI: " + fi + ", Cost: " + cost + "\n";
		return res;
	}
}

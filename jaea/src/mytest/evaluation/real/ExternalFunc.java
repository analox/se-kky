package mytest.evaluation.real;

import java.io.*;
import java.net.*;
import mytest.evaluation.*;


/**
 * Implementation of calling external evaluation function. <p></p>
 * @author Le Minh Nghia, NTU-Singapore
 */
public class ExternalFunc extends FitnessFunction {
	String dirName = "";
	String hostname = "localhost";
	public String getName() {return "External Evaluation Function";}
// 	public ExternalFunc(long initCount, String dir) {
// 		super(initCount);
// 		this.defaultRanges[0] = -100;
// 		this.defaultRanges[1] = 100;
// 		dirName = dir;
// 	}
	public ExternalFunc(long initCount, String host) {
		super(initCount);
		this.defaultRanges[0] = -100;
		this.defaultRanges[1] = 100;
		this.hostname = host;
	}
/*
	public double calculate(double[] inputs) {
		// increase number of evaluations
//		String fileName = "output" + NUM_EVAL+ ".txt";
		String fileName = "output.txt";
		double res = 0;
		// output input vector to file
		outputAParam(inputs);
		
//		System.out.println("Input ready. Waiting for result...");
		
		// wait for evaluation result
		File dir = new File(dirName);
		File file;
		
		do 
		{
			file = new File(dir, "otok.txt");
		} 
		while(!file.exists());
//		System.out.println("Result ready. Prepare tor read...");

		String s = null;	
		try {
			FileReader fr = new FileReader(dirName + fileName);
			BufferedReader input = new BufferedReader(fr);
			
			s = input.readLine();
			fr.close();
			if (s != null) {
					res = Double.parseDouble(s);
			}
			else {
				System.out.println("NUL NUL NUL");
				System.exit(0);
			}
			// read the return evaluation 
//			System.out.println("Evaluated " + NUM_EVAL + ": " + s + ", " +  res);
		} catch(IOException e) {
			System.err.println(e.getMessage());
		}
		
		while (!file.delete()); 
//		System.out.println("Delete output flag...");	
		return res;
	}
*/
	/** 
	Using TPC/IP Network protocol 
	*/
	public double calculate(double[] inputs) {
		double res = 0;
		int DATAPACK = 15;
		Socket evalSocket = null;
		DataOutputStream os = null;
		DataInputStream is = null;

		// evaluation port 30000
		try {
			evalSocket = new Socket(hostname,30000);
			os = new DataOutputStream(evalSocket.getOutputStream());
			is = new DataInputStream(evalSocket.getInputStream());
		}
		catch(UnknownHostException e) {
			System.err.println("Don't know hostname");
			System.exit(0);
		}
		catch(IOException e) {
			System.err.println("Couldn't get IO for connection");
			System.exit(0);
		}
		
		if (evalSocket != null && os !=null && is != null) {
			try {
				// Send data
				int nDim = inputs.length;
				String message = "";
				message += nDim;
				// Send number of inputs
				os.write(message.getBytes(), 0, message.getBytes().length);
				os.flush();

				String confirm = "?";
				System.out.println("...Waiting for confirm");
				while ((confirm = is.readLine()) != null) {
					System.out.println("Server: " + confirm);
					if (confirm.indexOf("1") != -1) break;
				}
				// Send pack of 15 variables
				message = "";
				for (int i = 0; i < nDim; i++) 
				{
    					if (i == nDim - 1) {
    						message += inputs[i];
						message = nDim + " " + message;
						System.out.println("Send: " + message + "\n at " + message.length() + "-" + message.getBytes().length);
						os.write(message.getBytes(), 0, message.getBytes().length);
						os.flush();

						System.out.println("...Waiting for confirm");
						while ((confirm = is.readLine()) != null) {
							System.out.println("Server: " + confirm);
							if (confirm.indexOf("1") != -1) break;
						}
					}
    					else 
					{
						if (i%DATAPACK ==0 && i > 0)
						{
							message = i + " " + message;
							System.out.println("Send: " + message + "\n at " + message.length() + "-" + message.getBytes().length);
							os.write(message.getBytes(), 0, message.getBytes().length);
							os.flush();
					
							System.out.println("...Waiting for confirm");
							while ((confirm = is.readLine()) != null) {
								System.out.println("Server: " + confirm);
								if (confirm.indexOf("1") != -1) break;
							}
							// reset message
							message = ""+ inputs[i] + " ";
						}
						else
    							message += inputs[i] + " ";
						
					}
				}
				//os.writeBytes(message);
				//os.writeBytes("4 1.1 2.2 3");
				//message = "4 1.1 2.2 3";
				


				// Receive fitness
				String responseLine;
				while ((responseLine = is.readLine()) != null) {
					System.out.println("Server: " + responseLine);
					res = Double.parseDouble(responseLine);
				}
				System.out.println("Done communication!");
				os.close();
				is.close();
				evalSocket.close();
			}
			catch(UnknownHostException e) {
				System.err.println("Don't know hostname");
				System.exit(0);
			}
			catch(IOException e) {
				System.err.println("Couldn't get IO for connection: " + e.getMessage());
				System.exit(0);
			}
		}
		return res;
	}
	public void outputAParam(double[] inputs) {
//		long count = NUM_EVAL;
		/* prepare log file */
    	PrintStream MyOutput = null;
    	try {
//    	       MyOutput = new PrintStream(new FileOutputStream(dirName + 
//    	    		   "input" + (count-1) +".txt"));
    		MyOutput = new PrintStream(new FileOutputStream(dirName + 
 	    		   "input.txt"));
    	} catch(Exception e) {}
    	int nDim = inputs.length;
    	MyOutput.println(nDim);
    	for (int i = 0; i < nDim; i++)
    		if (i == nDim - 1)
//    			MyOutput.print(inputs[i] + "\n");
    			MyOutput.print(inputs[i]);
    		else
    			MyOutput.print(inputs[i] + " ");
    	MyOutput.close();
    	/* Rename to input%counter%.txt 
    	File oldFile = new File(dirName + "input" + (count-1) +".txt");
    	File newFile = new File(dirName + "input" + count+".txt");
    	boolean success = oldFile.renameTo(newFile);
    	if (!success) {
    		System.err.println("IO Error - Can not rename the file");
    		System.exit(0);
    	}
    	*/
    	try {
  	       MyOutput = new PrintStream(new FileOutputStream(dirName + 
  	    		   "itok.txt"));
  	       MyOutput.close();
 		} catch(Exception e) {}
    	/* Rename to to%counter%.txt 
    	if (count == 1) {
    		try {
     	       MyOutput = new PrintStream(new FileOutputStream(dirName + 
     	    		   "itok" + count +".txt"));
     	       MyOutput.close();
    		} catch(Exception e) {}
    	}
    	else 
    	{
    		oldFile = new File(dirName + "itok" + (count-1) + ".txt");
    		newFile = new File(dirName + "itok" + count + ".txt");
    		success = oldFile.renameTo(newFile);
    		if (!success) {
    			System.err.println("IO Error - Can not rename the file");
    			System.exit(0);
    		}
    	}
    	*/
	}
}

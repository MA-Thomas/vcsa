package selfassembly;


public class returnRun {
	double[][] One;
	double[][] Two;
	int i_vec;
	int i_sim;
  
	public returnRun(double[][] vec_matrix, double[][] sim_matrix, int numVec, int numSim){
		One = vec_matrix;
		Two = sim_matrix;
		i_vec = numVec;
		i_sim = numSim;
	}
		
	public returnRun() {
		One = null;
		Two = null;
		i_vec = 0;
		i_sim = 0;
	}
}

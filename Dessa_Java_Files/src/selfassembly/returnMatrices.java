package selfassembly;


public class returnMatrices {
	double[][] One;
	double[][] Two;
  
	public returnMatrices(double[][] vec_matrix, double[][] sim_matrix) {
		One = vec_matrix;
		Two = sim_matrix;
	}
		
	public returnMatrices() {
		One = null;
		Two = null;
	}
}

package igtools.analyses.clustering;

public class SimilarityInfDiagonalTable {

	
	private double[][] table;
	
	
	public SimilarityInfDiagonalTable(int nofElements, Double defaultValue){
		
		table = new double[nofElements][];
		
		for(int i=0; i<nofElements; i++){
			
			table[i] = new double[i+1];
			
			for(int j=0; j<table[i].length; j++)
				table[i][j] = defaultValue;
		}
	}
	
	
	public void set(int i, int j, double value){
		if(i>=j){
			table[i][j] = value;
		}
		else{
			table[j][i] = value;
		}
	}
	
	public double get(int i, int j){
		if(i>=j){
			return table[i][j];
		}
		else{
			return table[j][i];
		}
	}
}

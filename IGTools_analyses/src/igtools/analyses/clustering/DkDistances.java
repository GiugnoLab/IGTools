package igtools.analyses.clustering;


/**
 * Distances between weighted vectors.
 * 
 * @author vbonnici
 *
 */
public abstract class DkDistances {

	
	public enum DISTANCE{JACCARD_SIMILARITY, GENERALIZED_JACCARD_DISTANCE, GENERALIZED_JACCARD_SIMILARITY};
	
	public abstract double distance(double[] v1, double[] v2, double[] weights);
		
	
	
	public static DkDistances factory(DISTANCE d){
		switch(d){
			case JACCARD_SIMILARITY:
				return new DkDistances(){
					@Override
					public double distance(double[] v1, double[] v2, double[] weights) {
						double intersection = 0.0;
						double union = 0.0;
						
						for(int i=0; i<v1.length; i++){
							if(v1[i]!=0 || v2[i]!=0){
								union++;
							}
							if(v1[i]!=0 && v2[i]!=0){
								intersection++;
							}
						}
						return intersection/union;
					}
				};
				
				
			case GENERALIZED_JACCARD_DISTANCE:
				return new DkDistances(){
					@Override
					public double distance(double[] v1, double[] v2, double[] weights) {
						double min = 0.0;
						double max = 0.0;
						
						for(int i=0; i<v1.length; i++){
							if(v1[i]==0 && v2[i]==0){
							}
							else if(v1[i]<v2[i]){
								min += v1[i] * weights[i];
								max += v2[i] * weights[i];
							}
							else{
								min += v2[i] * weights[i];
								max += v1[i] * weights[i];
							}
						}
						return min/max;
					}
				};
				
				
			case GENERALIZED_JACCARD_SIMILARITY:
				return new DkDistances(){
					@Override
					public double distance(double[] v1, double[] v2, double[] weights) {
						double dist = 0.0;
						
						for(int i=0; i<v1.length; i++){
							if(v1[i]==0 && v2[i]==0){
							}
							else if(v1[i]<v2[i]){
								dist += v1[i] / v2[i];
							}
							else{
								dist += v2[i] / v1[i];
							}
						}
						return dist;
					}
				};
				
				
			default:
				return null;
		}
	}
}

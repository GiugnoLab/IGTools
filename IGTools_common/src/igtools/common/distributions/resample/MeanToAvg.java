package igtools.common.distributions.resample;

import java.util.Map;

public class MeanToAvg extends DistributionResampler{

	@Override
	public void resample(double[][] distr) {
		double avg = 0.0;
		for(int i=0; i<distr.length; i++){
			avg += distr[i][1];
		}
		avg /= (double)distr.length;
		
		for(int i=0; i<distr.length; i++){
			distr[i][1] = (distr[i][1] + avg) / 2.0;
		}
	}

	@Override
	public void resample(Map<Double, Double> distr) {
		double avg = 0.0;
		for(Map.Entry<Double, Double> entry : distr.entrySet()){
			avg += entry.getValue();
		}
		avg /= (double)distr.size();
		for(Map.Entry<Double, Double> entry : distr.entrySet()){
			entry.setValue((entry.getValue() + avg)/2.0 );
		}
	}

}

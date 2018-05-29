package igtools.common.distributions.resample;

import java.util.Map;

public abstract class DistributionResampler {
	public abstract void resample(double[][] distr);
	public abstract void resample(Map<Double,Double> distr);
}

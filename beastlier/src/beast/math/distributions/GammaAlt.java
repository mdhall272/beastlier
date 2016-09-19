package beast.math.distributions;


import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import org.apache.commons.math.distribution.ContinuousDistribution;
import org.apache.commons.math.distribution.GammaDistributionImpl;
/*
If a _very_ informative gamma prior is put on latent periods, the run may have trouble starting up because it needs
short latent periods to find a starting state and the pdf of those periods will round to zero. This minor modification
fixes this by using a gamma distribution whose log pdf is calculated directly rather than as the logarithm of a
calculated pdf.
*/

@Description("Gamma distribution. for x>0  g(x;alpha,beta) = 1/Gamma(alpha) beta^alpha} x^{alpha - 1} e^{-\frac{x}{beta}}" +
        "If the input x is a multidimensional parameter, each of the dimensions is considered as a " +
        "separate independent component.")
public class GammaAlt extends ParametricDistribution {
    final public Input<RealParameter> alphaInput = new Input<>("alpha", "shape parameter, defaults to 2");
    final public Input<RealParameter> betaInput = new Input<>("beta", "scale parameter, defaults to 2");

    org.apache.commons.math.distribution.GammaDistribution m_dist = new GammaDistributionImplAlt(1, 1);

    @Override
    public void initAndValidate() {
        refresh();
    }

    /**
     * make sure internal state is up to date *
     */
    @SuppressWarnings("deprecation")
    void refresh() {
        double alpha;
        double beta;
        if (alphaInput.get() == null) {
            alpha = 2;
        } else {
            alpha = alphaInput.get().getValue();
        }
        if (betaInput.get() == null) {
            beta = 2;
        } else {
            beta = betaInput.get().getValue();
        }
        m_dist.setAlpha(alpha);
        m_dist.setBeta(beta);
    }

    @Override
    public ContinuousDistribution getDistribution() {
        refresh();
        return m_dist;
    }

    @Override
    public double getMean() {
        return offsetInput.get() + m_dist.getAlpha() * m_dist.getBeta();
    }
    // class Gamma

    public class GammaDistributionImplAlt extends GammaDistributionImpl {

        public GammaDistributionImplAlt(double alpha, double beta) {
            super(alpha, beta);
        }

        public GammaDistributionImplAlt(double alpha, double beta, double inverseCumAccuracy) {
            super(alpha, beta, inverseCumAccuracy);
        }

        public double logDensity(double x) {

            double shape = getAlpha();
            double scale = getBeta();

            return ((shape - 1.0) * (Math.log(x) - Math.log(scale)) - x / scale -
                    org.apache.commons.math.special.Gamma.logGamma(shape) - Math.log(scale));
        }

    }

}
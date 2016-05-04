package beast.core.util;

import beast.core.*;
import beast.core.Input.Validate;
import beast.core.parameter.BooleanParameter;
import beast.core.parameter.IntegerParameter;

import java.io.PrintStream;


@Description("calculates the exponential function")
public class ExponentialFunction extends CalculationNode implements Function, Loggable {
    final public Input<Function> functionInput = new Input<>("arg", "argument",
            Validate.REQUIRED);


    boolean needsRecompute = true;
    double value = 1;
    double storedValue = 1;

    @Override
    public void initAndValidate() {
        Function valuable = functionInput.get();
        if(valuable.getDimension()!=1){
            throw new IllegalArgumentException("Negative takes only one-dimensional variables as input");
        }
    }

    @Override
    public int getDimension() {
        return 1;
    }

    @Override
    public double getArrayValue() {
        if (needsRecompute) {
            compute();
        }
        return value;
    }

    /**
     * do the actual work, and reset flag *
     */
    void compute() {
        final double v = functionInput.get().getArrayValue(0);
        value = Math.exp(v);
        needsRecompute = false;
    }

    @Override
    public double getArrayValue(int dim) {
        if (dim == 0) {
            return getArrayValue();
        }
        return Double.NaN;
    }

    /**
     * CalculationNode methods *
     */
    @Override
    public void store() {
        storedValue = value;
        super.store();
    }

    @Override
    public void restore() {
        value = storedValue;
        super.restore();
    }

    @Override
    public boolean requiresRecalculation() {
        needsRecompute = true;
        return true;
    }

    /**
     * Loggable interface implementation follows
     */
    @Override
    public void init(PrintStream out) {
        out.print("product(" + ((BEASTObject) functionInput.get()).getID() + ")\t");
    }

    @Override
    public void log(int sampleNr, PrintStream out) {
        double value = Math.exp(functionInput.get().getArrayValue(0));
        out.print(value + "\t");
    }

    @Override
    public void close(PrintStream out) {
        // nothing to do
    }

} // class Sum

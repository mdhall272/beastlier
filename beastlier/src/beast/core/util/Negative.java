package beast.core.util;

import beast.core.*;
import beast.core.Input.Validate;
import beast.core.parameter.BooleanParameter;
import beast.core.parameter.IntegerParameter;

import java.io.PrintStream;


@Description("multiplies a one-dimensional valuable by -1")
public class Negative extends CalculationNode implements Function, Loggable {
    final public Input<Function> functionInput = new Input<>("arg", "argument to be multiplied by -1",
            Validate.REQUIRED);

    enum Mode {integer_mode, double_mode}

    Mode mode;

    boolean needsRecompute = true;
    double value = -1;
    double storedValue = -1;

    @Override
    public void initAndValidate() {
        Function valuable = functionInput.get();
        if(valuable.getDimension()!=1){
            throw new IllegalArgumentException("Negative takes only one-dimensional variables as input");
        }
        if (valuable instanceof IntegerParameter || valuable instanceof BooleanParameter) {
            mode = Mode.integer_mode;
        } else {
            mode = Mode.double_mode;
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
        value = -1*v;
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
        double value = -1*functionInput.get().getArrayValue(0);
        if (mode == Mode.integer_mode) {
            out.print((int) value + "\t");
        } else {
            out.print(value + "\t");
        }
    }

    @Override
    public void close(PrintStream out) {
        // nothing to do
    }

} // class Sum

package beast.core.util;

import beast.core.*;
import beast.core.Input.Validate;
import beast.core.parameter.BooleanParameter;
import beast.core.parameter.IntegerParameter;

import java.io.PrintStream;


@Description("calculates product of a valuable")
public class Product extends CalculationNode implements Function, Loggable {
    final public Input<Function> functionInput = new Input<>("arg", "argument to be multiplied", Validate.REQUIRED);

    enum Mode {integer_mode, double_mode}

    Mode mode;

    boolean needsRecompute = true;
    double product = 1;
    double storedProduct = 1;

    @Override
    public void initAndValidate() {
        Function valuable = functionInput.get();
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
        return product;
    }

    /**
     * do the actual work, and reset flag *
     */
    void compute() {
        product = 0;
        final Function v = functionInput.get();
        for (int i = 0; i < v.getDimension(); i++) {
            product *= v.getArrayValue(i);
        }
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
        storedProduct = product;
        super.store();
    }

    @Override
    public void restore() {
        product = storedProduct;
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
        Function valuable = functionInput.get();
        final int dimension = valuable.getDimension();
        double product = 1;
        for (int i = 0; i < dimension; i++) {
            product *= valuable.getArrayValue(i);
        }
        if (mode == Mode.integer_mode) {
            out.print((int) product + "\t");
        } else {
            out.print(product + "\t");
        }
    }

    @Override
    public void close(PrintStream out) {
        // nothing to do
    }

} // class Sum

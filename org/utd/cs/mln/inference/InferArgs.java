package org.utd.cs.mln.inference;

import com.beust.jcommander.*;


/**
 * This class is used for parsing command line arguments when running {@link InferTest}.
 * For this, it uses @see <a href="http://jcommander.org">JCommander</a> library.
 * Data structures for storing the command line arguments are also declared here, and initialized with their
 * default values.
 *
 *
 * For each argument, create an "@Parameter" annotation. It takes several arguments.
 * names : argument name i.e. flag
 * description : description of that flag
 * required : whether that flag is mandatory or not
 * order : The order in which arguments are shown in usage manual (default order is alphabetical)
 * Corresponding to each "@Parameter", create a data structure to store value for that flag.
 * @author Happy
 * @since 31/03/18
 * @see InferTest
 */

public class InferArgs {
    @Parameter(names="-i", description = "Input MLN file", required = true, order = 0)
    public String mlnFile;
    @Parameter(names="-o", description = "Output file", required = true, order = 1)
    public String outFile;
    @Parameter(names="-g", description = "Gold file", required = true, order = 2)
    public String goldFile;
    @Parameter(names="-e", description = "Evidence file", order = 3)
    public String evidFile;
    @Parameter(names="-priorSoftEvidence", description = "If specified, include priorSoftEvidence", order = 6)
    public boolean priorSoftEvidence;
    @Parameter(names="-se", description = "Softevidence file (-priorSoftEvidence flag must be set)", order = 7)
    public String softEvidenceFile;
    @Parameter(names="-sePred", description = "Soft Evidence predicate name (-priorSoftEvidence flag must be set)", order = 8)
    public String sePred;
    @Parameter(names="-seLambda", description = "Lambda for softevidence (-priorSoftEvidence flag must be set)", order = 9)
    public double seLambda = 1.0;
    @Parameter(names="-queryEvidence", description = "If specified, make all query ground predicates not specified in training file as evidence)", order = 10)
    public boolean queryEvidence;@Parameter(names="-agg", description = "If specified, use aggregator method", order = 2)
    public boolean agg;

    /**
     * @see InferArgs.GibbsParamConverter
     */
    @Parameter(names="-gibbs", description = "Inference arguments provided between double quotes", converter = InferArgs.GibbsParamConverter.class)
    public GibbsParams gibbsParam = new GibbsParams();

    /**
     * To print final values of parameters after parsing.
     * @return Final values of parameters as string
     */
    @Override
    public String toString() {
        return "Parameters given : \n" +
                "-i = " + mlnFile + "\n" +
                "-o = " + outFile + "\n"+
                "-g = " + goldFile + "\n"+
                "-e = " + evidFile + "\n"+
                "-priorSoftEvidence = " + priorSoftEvidence + "\n"+
                "-se = " + softEvidenceFile + "\n"+
                "-sePred = " + sePred  + "\n"+
                "-seLambda = " + seLambda + "\n"+
                "-queryEvidence = " + queryEvidence + "\n"+
                "-agg = " + agg + "\n"
                ;
    }

    /**
     * Tester function
     * To use JCommander, create its object and then call its parse method.
     * @param args command line arguments
     */
    public static void main(String args[])
    {
        InferArgs ia = new InferArgs();
        // Create JCommander object with object of InferArgs class passed. We pass the object of that class in which
        // @Parameter annotations were given.
        JCommander jc = JCommander.newBuilder().addObject(ia).build();
        // Screen column width i.e number of characters printed in one line. Default is 80.
        jc.setColumnSize(200);
        // When showing usage manual, show the program name as LearnTest.
        jc.setProgramName("InferTest");
        try
        {
            // parse the arguments and automatically assign the values to corresponding fields. If any error occurs,
            // it throws ParameterException
            jc.parse(args);
            // Validate some of the things in parameters.
            ia.validate();
            // Print the fields' set value after parsing
            System.out.println(ia);
        }
        catch(ParameterException p)
        {
            System.out.println(p.getMessage());
            // Print usage manual
            jc.usage();
        }

    }

    /**
     * Validates parameters.
     * Parameters related to softevidence must be given if priorSoftEvidence flag is set.
     */
    void validate()
    {
        if(priorSoftEvidence)
        {
            if(softEvidenceFile.isEmpty() || sePred.isEmpty())
                throw new ParameterException("Error !!! Need to specify softEvidence files and sePred");
        }
    }

    /**
     * This class is created for parsing the argument after -infer flag.
     * After -infer flag, user must give arguments for inference in double quotes.
     * This class reads those arguments and create a single gibbsParam object storing those arguments.
     */

    public static class GibbsParamConverter implements IStringConverter<GibbsParams>{
        @Parameter(names="-numChains", description = "Number of chains in gibbs sampling")
        public int numChains = 10;
        @Parameter(names="-maxSteps", description = "Max Number of samples per pred per chain")
        public int maxSteps = 1000;
        @Parameter(names="-samplesPerTest", description = "Number of samples after which convergence to be checked")
        public int samplesPerTest = 100;
        @Parameter(names="-burnMaxSteps", description = "Max Number of samples per pred per chain for burning")
        public int burnMaxSteps = 10;
        @Parameter(names="-testConvergence", description = "If specified, test for convergence", arity=1)
        public boolean testConvergence = true;


        /**
         * We created jc as static member here because we need to use it in above class LearnArgs when printing usage manual.
         */

        public static JCommander jc = JCommander.newBuilder().addObject(new InferArgs.GibbsParamConverter()).build();

        /**
         * This function will automatically gets called by above class's jc.parse() method when it tries to parse
         * arguments after -infer. It receives a single string (in double quotes).
         * @param s String in double quotes in which all inference parameters are given
         * @return GibbsParams object with fields set according to parameters given
         */
        @Override
        public GibbsParams convert(String s) throws ParameterException{
            // Remove double quotes
            s = s.replaceAll("\"","");

            // Create args array by splitinng arguments
            String args[] = s.split("\\s+");
            jc.setProgramName("-gibbs");

            jc.parse(args);

            GibbsParams gibbsparams = new GibbsParams();
            gibbsparams.numChains = ((GibbsParamConverter)jc.getObjects().get(0)).numChains;
            gibbsparams.maxSteps = ((GibbsParamConverter)jc.getObjects().get(0)).maxSteps;
            gibbsparams.samplesPerTest = ((GibbsParamConverter)jc.getObjects().get(0)).samplesPerTest;
            gibbsparams.burnMaxSteps = ((GibbsParamConverter)jc.getObjects().get(0)).burnMaxSteps;
            gibbsparams.testConvergence = ((GibbsParamConverter)jc.getObjects().get(0)).testConvergence;
            return gibbsparams;
        }
    }

}

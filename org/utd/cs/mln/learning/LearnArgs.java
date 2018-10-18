package org.utd.cs.mln.learning;

import com.beust.jcommander.*;
import com.beust.jcommander.converters.FileConverter;
import org.utd.cs.mln.inference.GibbsParams;

import java.util.List;

/**
 * This class is used for parsing command line arguments when running {@link LearnTest}.
 * For this, it uses @see <a href="http://jcommander.org">JCommander</a> library.
 * Data structures for storing the command line arguments are also declared here, and initialized with their
 * default values.
 * <p>
 * For each argument, create an "@Parameter" annotation. It takes several arguments.
 * names : argument name i.e. flag
 * description : description of that flag
 * required : whether that flag is mandatory or not
 * order : The order in which arguments are shown in usage manual (default order is alphabetical)
 * Corresponding to each "@Parameter", create a data structure to store value for that flag.
 * </p>
 * @author Happy
 * @since 3/31/18
 */

public class LearnArgs {
    @Parameter(names="-i", description = "Input MLN file", required = true, order = 0)
    public String mlnFile;
    @Parameter(names="-o", description = "Output file", required = true, order = 1)
    public String outFile;
    @Parameter(names="-t", description = "Comma separated training files", required = true, order = 2)
    public List<String> truthFiles;
    @Parameter(names="-e", description = "Comma separated evidence files", order = 3)
    public List<String> evidenceFiles;
    @Parameter(names="-priorSoftEvidence", description = "If specified, include priorSoftEvidence", order = 8)
    public boolean priorSoftEvidence;
    @Parameter(names="-se", description = "Comma separated softevidence files (-priorSoftEvidence flag must be set)", order = 9)
    public List<String> softEvidenceFiles;
    @Parameter(names="-sePred", description = "Soft Evidence predicate name (-priorSoftEvidence flag must be set)", order = 10)
    public String sePred;
    @Parameter(names="-seLambda", description = "Lambda for softevidence (-priorSoftEvidence flag must be set)", order = 11)
    public double seLambda = 1.0;
    @Parameter(names="-withEm", description = "If specified, learn with EM", order = 12)
    public boolean withEM;
    @Parameter(names="-numEMSamples", description = "Number of samples in E step of EM (-withEM flag must be set)", order = 14)
    public int numEMSamples = 5;
    @Parameter(names="-queryEvidence", description = "If specified, make all query ground predicates not specified in training file as evidence)", order = 15)
    public boolean queryEvidence;
    @Parameter(names="-useMlnWts", description = "If specified, initialize weights with MLN weights", order = 16)
    public boolean useMlnWts;
    @Parameter(names="-method", description = "Method Name (cg/lbfgs)", order = 19)
    public String method = "cg";
    @Parameter(names="-pll", description = "If specified, use neg pseudo log likelihood as loss function", order = 20)
    public boolean pll = false;
    @Parameter(names="-usePrior", description = "If specified, use prior gaussian distribution", order = 21)
    public boolean usePrior = false;
    @Parameter(names="-debug", description = "If specified, prints gradient and counts at every iteration", order = 19)
    public boolean debug;
    @Parameter(names="-agg", description = "If specified, use aggregator method", order = 2)
    public boolean agg;
    /**
     * @see GibbsParamConverter
     */
    @Parameter(names="-infer", description = "Inference arguments provided between double quotes", converter = GibbsParamConverter.class)
    public GibbsParams gibbsParam = new GibbsParams();

    /**
     * @see CGParamsConverter
     */
    @Parameter(names="-cg", description = "cg params provided between double quotes", converter = CGParamsConverter.class)
    public CGParams cgParam = new CGParams();

    /**
     * To print final values of parameters after parsing.
     * @return Final values of parameters as string
     */
    @Override
    public String toString() {
        return "Parameters given : \n" +
                "-i = " + mlnFile + "\n" +
                "-o = " + outFile + "\n"+
                "-t = " + truthFiles + "\n"+
                "-e = " + evidenceFiles + "\n"+
                "-priorSoftEvidence = " + priorSoftEvidence + "\n"+
                "-se = " + softEvidenceFiles + "\n"+
                "-sePred = " + sePred  + "\n"+
                "-seLambda = " + seLambda + "\n"+
                "-withEM = " + withEM + "\n"+
                "-numEMSamples = " + numEMSamples + "\n"+
                "-queryEvidence = " + queryEvidence + "\n"+
                "-useMlnWts = " + useMlnWts + "\n"+
                "-method = " + method + "\n" +
                "-pll = " + pll + "\n" +
                "-debug = " + debug + "\n" +
                "-agg = " + agg + "\n" +
                "-usePrior = " + usePrior + "\n"
                ;
    }

    /**
     * Tester function
     * To use JCommander, create its object and then call its parse method.
     * @param args command line arguments
     */
    public static void main(String args[])
    {
        LearnArgs la = new LearnArgs();
        // Create JCommander object with object of LearnArgs class passed. We pass the object of that class in which
        // @Parameter annotations were given.
        JCommander jc = JCommander.newBuilder().addObject(la).build();
        // Screen column width i..e number of characters printed in one line. Default is 80.
        jc.setColumnSize(200);
        // When showing usage manual, show the program name as LearnTest.
        jc.setProgramName("LearnTest");
        try
        {
            // parse the arguments and automatically assign the values to corresponding fields. If any error occurs,
            // it throws ParameterException
            jc.parse(args);
            // Validate some of the things in parameters.
            la.validate();
            // Print the fields' set value after parsing
            System.out.println(la);
        }
        catch(ParameterException p)
        {
            System.out.println(p.getMessage());
            // Print usage manual
            jc.usage();
            // Print usage manual of -infer also
            GibbsParamConverter.jc.usage();
        }

    }

    /**
     * Validates parameters
     * <ol>
     *     <li>Both genLearn and discLearn can't be on at the same time.</li>
     *     <li>Parameters related to softevidence must be given if priorSoftEvidence flag is set.</li>
     *     <li>EM related parameters must be given if withEM is set.</li>
     * </ol>
     */
    void validate()
    {
        if(priorSoftEvidence)
        {
            if(softEvidenceFiles == null || sePred.isEmpty())
                throw new ParameterException("Error !!! Need to specify softEvidence files and sePred");
        }
        if(!method.equals("cg") && !method.equals("lbfgs"))
        {
            throw new ParameterException("Error !!! Invalid method " + method);
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
        public int maxSteps = 100;
        @Parameter(names="-samplesPerTest", description = "Number of samples after which check for convergence")
        public int samplesPerTest = 10;
        @Parameter(names="-burnMaxSteps", description = "Max Number of samples per pred per chain for burning")
        public int burnMaxSteps = 10;
        @Parameter(names="-testConvergence", description = "If specified, test for convergence")
        public boolean testConvergence = true;

        /**
         * We created jc as static member here because we need to use it in above class LearnArgs when printing usage manual.
         */

        public static JCommander jc = JCommander.newBuilder().addObject(new GibbsParamConverter()).build();

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
            jc.setProgramName("-infer");

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

    /**
     * This class is created for parsing the argument after -cg flag.
     * After -cg flag, user must give arguments for CG in double quotes.
     * This class reads those arguments and create a single cgParam object storing those arguments.
     */

    public static class CGParamsConverter implements IStringConverter<CGParams>{
        @Parameter(names="-lambda", description = "Initial Lambda for CG")
        public double lambda = 100.0;
        @Parameter(names="-maxLambda", description = "Max Lambda for CG")
        public double maxLambda = Double.MAX_VALUE;
        @Parameter(names="-maxBackTracks", description = "Max Lambda for CG")
        public int maxBackTracks = 1000;
        @Parameter(names="-preCondition", description = "If specified do preconditioning")
        public boolean preCondition = true;
        @Parameter(names="-minllChange", description = "convergence criteria for CG")
        public double minllChange = 0.001;
        @Parameter(names="-numIter", description = "max number of iterations")
        public int numIter = 100;

        /**
         * We created jc as static member here because we need to use it in above class LearnArgs when printing usage manual.
         */

        public static JCommander jc = JCommander.newBuilder().addObject(new CGParamsConverter()).build();

        /**
         * This function will automatically gets called by above class's jc.parse() method when it tries to parse
         * arguments after -cg. It receives a single string (in double quotes).
         * @param s String in double quotes in which all cg parameters are given
         * @return CGParams object with fields set according to parameters given
         */
        @Override
        public CGParams convert(String s) throws ParameterException{
            // Remove double quotes
            s = s.replaceAll("\"","");

            // Create args array by splitinng arguments
            String args[] = s.split("\\s+");
            jc.setProgramName("-cg");

            jc.parse(args);
            CGParams cgParams = new CGParams();
            cgParams.cg_lambda = ((CGParamsConverter)jc.getObjects().get(0)).lambda;
            cgParams.cg_max_lambda = ((CGParamsConverter)jc.getObjects().get(0)).maxLambda;
            cgParams.maxBackTracks = ((CGParamsConverter)jc.getObjects().get(0)).maxBackTracks;
            cgParams.preConditionCG = ((CGParamsConverter)jc.getObjects().get(0)).preCondition;
            cgParams.numIter = ((CGParamsConverter)jc.getObjects().get(0)).numIter;
            cgParams.min_ll_change = ((CGParamsConverter)jc.getObjects().get(0)).minllChange;
            return cgParams;
        }
    }
}

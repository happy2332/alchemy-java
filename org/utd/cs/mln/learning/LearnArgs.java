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
    @Parameter(names="-g", description = "If specified, do generative learning", order = 4)
    public boolean genLearn = false;
    @Parameter(names="-d", description = "If specified, do discriminative learning", order = 5)
    public boolean discLearn = true;
    @Parameter(names="-ep", description = "Comma separated evidence predicates (default : all except query and hidden predicates", order = 6)
    public List<String> evidPreds;
    @Parameter(names="-qp", description = "Comma separated query predicates", order = 7)
    public List<String> queryPreds;
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
    @Parameter(names="-hp", description = "Comma separated hidden predicates (-withEM flag must be set)", order = 13)
    public List<String> hiddenPreds;
    @Parameter(names="-numEMSamples", description = "Number of samples in E step of EM (-withEM flag must be set)", order = 14)
    public int numEMSamples = 10;
    @Parameter(names="-queryEvidence", description = "If specified, make all query ground predicates not specified in training file as evidence)", order = 15)
    public boolean queryEvidence;
    @Parameter(names="-usePrior", description = "If specified, initialize weights with MLN weights", order = 16)
    public boolean usePrior;
    @Parameter(names="-numIter", description = "Number of learning iterations in dicriminative learning", order = 17)
    public int numIter = 100;
    @Parameter(names="-minllChange", description = "convergence criteria for discriminative learning", order = 17)
    public double minllChange = 0.001;
    /**
     * @see GibbsParamConverter
     */
    @Parameter(names="-infer", description = "Inference arguments provided between double quotes", converter = GibbsParamConverter.class)
    public GibbsParams gibbsParam;

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
                "-g = " + genLearn + "\n"+
                "-d = " + discLearn + "\n"+
                "-ep = " + evidPreds + "\n"+
                "-qp = " + queryPreds + "\n"+
                "-priorSoftEvidence = " + priorSoftEvidence + "\n"+
                "-se = " + softEvidenceFiles + "\n"+
                "-sePred = " + sePred  + "\n"+
                "-seLambda = " + seLambda + "\n"+
                "-withEM = " + withEM + "\n"+
                "-hp = " + hiddenPreds + "\n"+
                "-numEMSamples = " + numEMSamples + "\n"+
                "-queryEvidence = " + queryEvidence + "\n"+
                "-usePrior = " + usePrior + "\n"+
                "-numIter = " + numIter + "\n"+
                "-minllChange = " + minllChange + "\n"
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
        if(genLearn) {
            discLearn = false;
            if(queryPreds != null || evidPreds != null)
                throw new ParameterException("Error !!! Can't provide query or evidence predicates in generative learning");
        }
        if(discLearn)
        {
            if(queryPreds == null)
                throw new ParameterException("Error !!! Necessary to provide query predicates in discriminative learning (use -qp)");
        }
        if(priorSoftEvidence)
        {
            if(softEvidenceFiles == null || sePred.isEmpty())
                throw new ParameterException("Error !!! Need to specify softEvidence files and sePred");
        }
        if(withEM)
        {
            if(hiddenPreds == null)
                throw new ParameterException("Error !!! Need to specify hidden preds");
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
        @Parameter(names="-numSamples", description = "Number of samples per pred per chain")
        public int numSamples = 100;

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
            gibbsparams.numChains = numChains;
            gibbsparams.samplesPerTest = numSamples;
            return gibbsparams;
        }
    }
}

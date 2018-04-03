package org.utd.cs.mln.inference;

import com.beust.jcommander.*;


import java.util.List;

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
    @Parameter(names="-ep", description = "Comma separated evidence predicates (default : all except query predicates", order = 4)
    public List<String> evidPreds;
    @Parameter(names="-qp", description = "Comma separated query predicates", order = 5)
    public List<String> queryPreds;
    @Parameter(names="-priorSoftEvidence", description = "If specified, include priorSoftEvidence", order = 6)
    public boolean priorSoftEvidence;
    @Parameter(names="-se", description = "Softevidence file (-priorSoftEvidence flag must be set)", order = 7)
    public String softEvidenceFile;
    @Parameter(names="-sePred", description = "Soft Evidence predicate name (-priorSoftEvidence flag must be set)", order = 8)
    public String sePred;
    @Parameter(names="-seLambda", description = "Lambda for softevidence (-priorSoftEvidence flag must be set)", order = 9)
    public double seLambda = 1.0;
    @Parameter(names="-queryEvidence", description = "If specified, make all query ground predicates not specified in training file as evidence)", order = 10)
    public boolean queryEvidence;
    @Parameter(names="-numChains", description = "Number of chains in gibbs sampling", order = 11)
    public int numChains = 10;
    @Parameter(names="-numSamples", description = "Number of samples per pred per chain", order = 12)
    public int numSamples = 100;

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
                "-ep = " + evidPreds + "\n"+
                "-qp = " + queryPreds + "\n"+
                "-priorSoftEvidence = " + priorSoftEvidence + "\n"+
                "-se = " + softEvidenceFile + "\n"+
                "-sePred = " + sePred  + "\n"+
                "-seLambda = " + seLambda + "\n"+
                "-queryEvidence = " + queryEvidence + "\n"+
                "-numChains = " + numChains + "\n"+
                "-numSamples = " + numSamples + "\n"
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
}

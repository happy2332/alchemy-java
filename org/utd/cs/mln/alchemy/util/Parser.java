package org.utd.cs.mln.alchemy.util;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.InputStreamReader;
import java.util.*;

//import com.sun.tools.javac.comp.Infer;
import org.utd.cs.gm.core.LogDouble;
import org.utd.cs.gm.utility.Pair;
import org.utd.cs.mln.alchemy.core.*;
import org.utd.cs.mln.inference.InferTest;

public class Parser {
	public static final String DOMAINSTART = "#domains";
	public static final String VALUESSTART = "#values";
	public static final String PREDICATESTART = "#predicates";
	public static final String FORMULASTART = "#formulas";
	public static final String LEFTPRNTH = "(";
	public static final String RIGHTPRNTH = ")";
	public static final String NOTOPERATOR = "!";
	public static final String ANDOPERATOR = "^";
	public static final String OROPERATOR = "|";
	public static final String LEFTFLOWER = "{";
	public static final String RIGHTFLOWER = "}";
	public static final String WEIGHTSEPARATOR = "::";
	public static final String COMMASEPARATOR = ",";
	public static final String EQUALSTO = "=";
	public static final String ELLIPSIS = "...";
    public static final String COMMENT = "//";

	private static final String REGEX_ESCAPE_CHAR = "\\";
    private List<String> truthEvidFiles;
    public Set<String> queryPreds = new HashSet<>(), evidPreds = new HashSet<>(), hiddenPreds = new HashSet();


    public void collectDomain(List<String> files) throws FileNotFoundException {

        for(String file : files)
        {
            if(file != null)
            {
                Scanner scanner = new Scanner(new BufferedReader(new InputStreamReader(new FileInputStream(file))));
                while (scanner.hasNextLine()) {
                    String line = scanner.nextLine().replaceAll("\\s", "");

                    if (line.isEmpty()) {
                        continue;
                    }
                    String[] predArr = line.split(REGEX_ESCAPE_CHAR + LEFTPRNTH);
                    String symbolName = predArr[0];
                    String[] predArr2 = predArr[1].split(EQUALSTO);
                    String[] termNames = predArr2[0].replace(RIGHTPRNTH, "").split(COMMASEPARATOR);
                    for(PredicateSymbol symbol : mln.symbols)
                    {
                        if(symbol.symbol.equals((symbolName)))
                        {
                            for(int i = 0 ; i < termNames.length ; i++)
                            {
                                String var_type = symbol.variable_types.get(i);
                                mln.varTypeToDomainMap.get(var_type).add(Integer.parseInt(termNames[i]));
                            }
                            break;
                        }
                    }

                }
                scanner.close();
            }
        }

    }

    public void setTruthEvidFiles(List<String> truthEvidFiles) {
        this.truthEvidFiles = truthEvidFiles;
    }


    private enum ParserState {
		Domain,
		Values, // Added by Happy
		Predicate,
		Formula;
	};

	private MLN mln;

	private int predicateId;

	//set of domains with values in String format
	//private List<Domain> domainList = new ArrayList<Domain>();
    //private Map<String, Set<Integer>> varTypeToDomainMap = new HashMap<>();

	//List of Values
	private List<Values> valuesList = new ArrayList<Values>();

	//map key:predicate Id, value:type of each term
	//private Map<Integer,List<String>> predicateDomainMap = new HashMap<>();


	public Parser(MLN mln_) {
		mln = mln_;
		predicateId = 0;
	}


	boolean isTermConstant(String term)
	{
		//if its a number, it is taken as a constant
        try
        {
            int num = Integer.parseInt(term);
            return true;
        }
        catch(Exception e)
        {
            return false;
        }
	}

	Term create_new_term(String var_type)
	{
		//map domain to an integer domain for ease of manipulation
		List<Integer> iDomain = new ArrayList<Integer>(mln.varTypeToDomainMap.get(var_type));
		//create a new term
		Term term = new Term(var_type,iDomain);
		return term;
	}

	WClause create_new_clause(List<Integer> predicateSymbolIndex,List<Boolean> sign,
			List<List<Term> > iTermsList, List<Integer> valTrueList)
	{
		int numAtoms = predicateSymbolIndex.size();
		WClause clause = new WClause();
		clause.atoms = new ArrayList<Atom>(numAtoms);
		clause.satisfied = false;
		clause.sign = sign;
		clause.valTrue = valTrueList;
		Set<Term> termsSet = new HashSet<>();
		for(int i=0;i<numAtoms;i++)
		{
			PredicateSymbol symbol = MLN.create_new_symbol(mln.symbols.get(predicateSymbolIndex.get(i)));
			List<Term> terms = iTermsList.get(i);
			for(int j = 0 ; j < terms.size() ; j++)
            {
                terms.get(j).type = new String(symbol.variable_types.get(j));
            }
			Atom atom = new Atom(symbol,terms);
			termsSet.addAll(terms);
			clause.atoms.add(atom);
		}
		clause.terms = new ArrayList<>(termsSet);
		return clause;
	}

    private Formula create_new_formula(Double weight, List<Integer> predicateSymbolIndex, List<Boolean> sign,
                                       List<Integer> valTrueList, List<List<Term>> iTermsList,
                                       List<Integer> clausePartitionIndex, boolean tautology)
    {
        List<WClause> CNF = new ArrayList<WClause>();
        Set<Term> termsSet = new HashSet<>();
        for(int i = 0 ; i < clausePartitionIndex.size()-1 ; i++)
        {
            int from = clausePartitionIndex.get(i);
            int to = clausePartitionIndex.get(i+1);
            WClause new_clause = create_new_clause(predicateSymbolIndex.subList(from,to), sign.subList(from,to), iTermsList.subList(from, to), valTrueList.subList(from,to));
            new_clause.weight = new LogDouble(weight, true);
            termsSet.addAll(new_clause.terms);
            CNF.add(new_clause);
        }
        Formula formula = new Formula(CNF, new LogDouble(weight, true), tautology);
        formula.terms = new ArrayList<>(termsSet);
        return formula;
    }

	// Parse formula string. String ex : !S(x)=1 ^ C(x)=3 | F(x,y)=1::5
    // This will be 2 clauses : 1. !S(x)=1, 2. C(x)=3 V F(x,y)=1, each with weight 5
    // This formula is true if (i) S(x) is not 1, and (ii) Either C(x) is 3 or F(x,y)=1

	void parseCNFString(String line)
	{
	    boolean tautology;
        String[] formulaArr = line.split(WEIGHTSEPARATOR);
        Double weight;
        try {
            tautology = false;
            weight = Double.parseDouble(formulaArr[1]);
        } catch(Exception e){
            tautology = true;
            weight = Double.valueOf(1000.0);
        }
        List<Integer> clausePartitionIndex = new ArrayList<>(); // Stores list of indexes of atoms at which clause partitions. For example, if formula a | b ^ c ^ d|e, then this list will be [0,2,3,5] i.e. partition at atom b, then atom c, and then atom e.
        clausePartitionIndex.add(0);
        String formulaString = formulaArr[0];

        //If a formula starts with parenthesis, remove it
        if(formulaString.startsWith(LEFTPRNTH)) {

            if(!formulaString.endsWith(RIGHTPRNTH)) {
                System.out.println("Missing right parenthesis in clause " + formulaString);
                System.exit(-1);
            }

            formulaString = formulaString.substring(1, formulaString.length() - 1);
        }

        List<Boolean> sign = new ArrayList<Boolean>();
        List<Integer> valTrueList = new ArrayList<>();
        List<Integer> predicateSymbolIndex = new ArrayList<Integer>();
        List<List<String>> sTermsList = new ArrayList<List<String>>();

        String[] clauseStrings = formulaString.split(REGEX_ESCAPE_CHAR + ANDOPERATOR);
        ArrayList<String> atomStrings = new ArrayList<>();
        for(String clauseString : clauseStrings) {
            String[] clauseAtomStrings = clauseString.split(REGEX_ESCAPE_CHAR + OROPERATOR);
            for (int i = 0; i < clauseAtomStrings.length; i++) {
                atomStrings.add(clauseAtomStrings[i]);
            }
            clausePartitionIndex.add(atomStrings.size());
        }

        for (int i = 0; i < atomStrings.size(); i++) {
            sign.add(false);
            predicateSymbolIndex.add(null);
        }
        //System.out.println("line = " + line);
        for(int i=0; i<atomStrings.size(); i++) {
            String[] atomStrings1 = atomStrings.get(i).split(REGEX_ESCAPE_CHAR + EQUALSTO);
            String atomString = atomStrings1[0];
            int valTrue = Integer.parseInt(atomStrings1[1]);
            valTrueList.add(valTrue);
            //find opening and closing braces
            int startpos = atomString.indexOf(LEFTPRNTH);
            String predicateName = atomString.substring(0, startpos);

            if(predicateName.startsWith(NOTOPERATOR))
            {
                //negation
                predicateName = predicateName.substring(1,predicateName.length());
                sign.set(i, true);
            }

            for (int k = 0; k < mln.symbols.size(); k++) {
                //found the predicate
                if (mln.symbols.get(k).symbol.equals(predicateName)) {
                    predicateSymbolIndex.set(i, k);
                    break;
                }
            }

            int endpos = atomString.indexOf(RIGHTPRNTH);
            String termsString = atomString.substring(startpos + 1, endpos);
            //System.out.println("termsString = "+termsString);
            String[] terms = termsString.split(COMMASEPARATOR);
            sTermsList.add(new ArrayList<String>(Arrays.asList(terms)));

            //check if the number of terms is equal to the declared predicate
            if (terms.length != mln.symbols.get(predicateSymbolIndex.get(i)).variable_types.size()) {
                System.out.println("Error! Number/domain of terms in the predicate delcaration does not match in formula. " + predicateName);
                System.exit(-1);
            }
        }

        //create required terms
        List<List<Term> > iTermsList = new ArrayList<List<Term>>();
        for (int i = 0; i < atomStrings.size(); i++) {
            iTermsList.add(null);
        }
        for(int j=0;j<atomStrings.size();j++)
        {
            //for each term of atom j, check if it has already appeared in previous atoms of formula
            List<Term> iTerms = new ArrayList<Term>();
            for (int i = 0; i < sTermsList.get(j).size(); i++) {
                iTerms.add(null);
            }

            for(int k=0; k<sTermsList.get(j).size(); k++)
            {
                int pIndex = predicateSymbolIndex.get(j);
                PredicateSymbol ps = mln.symbols.get(pIndex);
                String var_type = ps.variable_types.get(k);
                //int domainIndex = predicateDomainMap.get(pIndex).get(k);

                //if term is a constant must be a unique term
                if(isTermConstant(sTermsList.get(j).get(k)))
                {
                    //find the id of the term
                    int id = Integer.parseInt((sTermsList.get(j).get(k)));
                    mln.varTypeToDomainMap.get(var_type).add(id);
//                    for(int m=0;m<domainList.get(domainIndex).values.size();m++)
//                    {
//                        if(domainList.get(domainIndex).values.get(m).equals(sTermsList.get(j).get(k)))
//                        {
//                            id=m;
//                            break;
//                        }
//                    }
//                    if(id==-1)
//                    {
//                        System.out.println("Constant does not match predicate's domain. "  + domainList.get(domainIndex).name );
//                        System.exit(-1);
//                    }
                    iTerms.set(k, new Term(var_type,id));
                }
                else
                {
                    boolean isExistingTerm = false;
                    boolean sameTerminPred = false;
                    int atomIndex=-1;
                    int termIndex=-1;
                    //check in term lists for atoms 0 to j;
                    for(int m=0;m<=j;m++)
                    {
                        for(int n=0; n < (m != j ? sTermsList.get(m).size() : k) ; n++)
                        {
                            if(sTermsList.get(m).get(n).equals(sTermsList.get(j).get(k)))
                            {
                                //check if the domains of the matched variables are the same
                                int atomSymbolIndex1 = predicateSymbolIndex.get(m);
                                int atomId1 = mln.symbols.get(atomSymbolIndex1).id;
                                String varType1 = mln.predicateDomainMap.get(atomId1).get(n);

                                int atomSymbolIndex2 = predicateSymbolIndex.get(j);
                                int atomId2 = mln.symbols.get(atomSymbolIndex2).id;
                                String varType2 = mln.predicateDomainMap.get(atomId2).get(k);
                                if(!varType1.equals(varType2))
                                {
                                    System.out.println("Error! variables do not match type ." + atomStrings.get(j) + "(" + varType1 + ", " + varType2 + ")" );
                                    System.exit(-1);
                                }
                                //variable is repeated, use the term created for atom m, term n
                                isExistingTerm = true;
                                atomIndex = m;
                                termIndex = n;
                                if(m == j)
                                {
                                    sameTerminPred = true;
                                }
                                break;
                            }
                        }
                        if(isExistingTerm)
                            break;
                    }

                    if(sameTerminPred)
                    {
                        iTerms.set(k, iTerms.get(termIndex));
                    }
                    else if(isExistingTerm)
                    {
                        //use the terms created for previous atoms
                        iTerms.set(k, iTermsList.get(atomIndex).get(termIndex));
                    }
                    else
                    {
                        //create a new Term
                        iTerms.set(k, create_new_term(var_type));
                    }
                }
            }
            iTermsList.set(j, iTerms);
        }//j atoms

        Formula newFormula = create_new_formula(weight, predicateSymbolIndex, sign, valTrueList, iTermsList, clausePartitionIndex, tautology);
        newFormula.formulaId = mln.formulas.size();
        mln.formulas.add(newFormula);
	}


    void parseDomainString(String line)
	{
		String[] domainArr = line.split(EQUALSTO);
		String domainName = domainArr[0];

		String[] domValArr = domainArr[1].replace(LEFTFLOWER, "").replace(RIGHTFLOWER, "").split(COMMASEPARATOR);
		List<String> domainValues = new ArrayList<String>();
		for (int i = 0; i < domValArr.length; i++) {
			if(domValArr[i].equals(ELLIPSIS)) {
				Integer startNumber = Integer.parseInt(domValArr[i-1]);
				Integer endNumber   = Integer.parseInt(domValArr[i+1]);
				for (int j = startNumber+1; j < endNumber; j++) {
					domainValues.add(Integer.toString(j));
				}
			} else {
				domainValues.add(domValArr[i]);
			}
		}

		//Domain domain = new Domain(domainName,domainValues);
        mln.varTypeToDomainMap.put(domainName, new HashSet<Integer>());
		//domainList.add(domain);
	}

	// Added by Happy
    private void parseValuesString(String line) {
        String[] valuesArr = line.split(EQUALSTO);
        String valuesName = valuesArr[0];

        String[] valuesValArr = valuesArr[1].replace(LEFTFLOWER, "").replace(RIGHTFLOWER, "").split(COMMASEPARATOR);
        ArrayList<Integer> valuesVals = new ArrayList<Integer>();
        for (int i = 0; i < valuesValArr.length; i++) {
            if(valuesValArr[i].equals(ELLIPSIS)) {
                Integer startNumber = Integer.parseInt(valuesValArr[i-1]);
                Integer endNumber   = Integer.parseInt(valuesValArr[i+1]);
                for (int j = startNumber+1; j < endNumber; j++) {
                    valuesVals.add(j);
                }
            } else {
                valuesVals.add(Integer.parseInt(valuesValArr[i]));
            }
        }

        Values values = new Values(valuesName, valuesVals);
        valuesList.add(values);
    }

	void parsePredicateString(String line)
	{
		String[] predArr = line.split(REGEX_ESCAPE_CHAR + LEFTPRNTH);
		String symbolName = predArr[0];
		boolean isQuery = false, isEvidence = false, isHidden = false;
		if(symbolName.charAt(0) == '*')
        {
            isQuery = true;
            if(symbolName.charAt(1) == '@')
            {
                isEvidence = true;
            }
        }
        else if(symbolName.charAt(0) == '@')
        {
            isEvidence = true;
            if(symbolName.charAt(1) == '*')
            {
                isQuery = true;
            }
        }
        else
        {
            isHidden = true;
        }
        symbolName = symbolName.replaceAll("[*@]","");
        if(isQuery)
        {
            queryPreds.add(symbolName);
        }
        if(isEvidence)
        {
            evidPreds.add(symbolName);
        }
        if(isHidden)
        {
            hiddenPreds.add(symbolName);
        }
		String[] predArr2 = predArr[1].split(EQUALSTO);
		String valuesName = predArr2[1];
		String[] termNames = predArr2[0].replace(RIGHTPRNTH, "").split(COMMASEPARATOR);

		List<String> var_types = new ArrayList<>();
		for(int m=0; m < termNames.length; m++) {
		    if(!mln.varTypeToDomainMap.containsKey(termNames[m])) {
                System.out.println("Error !!! " + termNames[m] + " type doesn't exist");
                System.exit(-1);
            }
            var_types.add(termNames[m]);
		}

		int matchingIndex = -1;
		for(int i = 0 ; i < valuesList.size() ; i++)
        {
            if(valuesName.equals(valuesList.get(i).name))
            {
                matchingIndex = i;
                break;
            }
        }
        if(matchingIndex == -1)
        {
            System.out.println("Error! Value name does not exist for predicate. " + symbolName );
            System.exit(-1);
        }
		//create a new predicate symbol
		PredicateSymbol p = new PredicateSymbol(predicateId,symbolName,var_types, valuesList.get(matchingIndex), LogDouble.ONE,LogDouble.ONE);
		//predicateList.push_back(p);
        mln.symbols.add(p);

		//Build the map for this predicate;
		//For predicateid, generate a List of domainIds that index Domains
		List<Integer> domainIndex = new ArrayList<Integer>();
		for(int i=0; i<termNames.length; i++)
		{
            if(!mln.varTypeToDomainMap.containsKey(termNames[i]))
			{
				System.out.println("Error! Domain does not exist for predicate. " + symbolName );
				System.exit(-1);
			}
			domainIndex.add(matchingIndex);
		}
		mln.predicateDomainMap.put(predicateId, var_types);
		//increment predicateid
		predicateId++;
	}


    /**
     *
     * @param filename Input MLN file name
     * @throws FileNotFoundException
     */
	public void parseInputMLNFile(String filename) throws FileNotFoundException
	{
		Scanner scanner = new Scanner(new BufferedReader(new InputStreamReader(new FileInputStream(filename))));
		ParserState state = null;

		while(scanner.hasNextLine()) {
			String line = scanner.nextLine().replaceAll("\\s","");

			if(line.isEmpty() || line.contains(COMMENT)) {
				continue;
			}

			if(line.contains(DOMAINSTART)) {
				state = ParserState.Domain;
				continue;
			} else if(line.contains(VALUESSTART)) {
                state = ParserState.Values;
                continue;
            }
			else if (line.contains(PREDICATESTART)) {
				state = ParserState.Predicate;
				continue;
			} else if (line.contains(FORMULASTART)) {
				state = ParserState.Formula;
				continue;
			}

			switch (state) {
			case Domain:
				parseDomainString(line);
				break;

            case Values:
                parseValuesString(line);
                break;

			case Predicate:
				parsePredicateString(line);
				mln.max_predicate_id = predicateId-1;
				break;

			case Formula:
			    // Fill in varTypeToDomain
                collectDomain(truthEvidFiles);
				parseCNFString(line);
				break;

			default:
				break;
			}
		}

		scanner.close();
	}


    public ArrayList<FirstEvidence> parseInputEvidenceFile(String filename) throws FileNotFoundException, PredicateNotFound
    {
        Scanner scanner = new Scanner(new BufferedReader(new InputStreamReader(new FileInputStream(filename))));
        ArrayList<FirstEvidence> evidList = new ArrayList<FirstEvidence>();
        while(scanner.hasNextLine()) {
            String line = scanner.nextLine().replaceAll("\\s","");

            if(line.isEmpty() || line.startsWith(COMMENT)) {
                continue;
            }

            evidList.add(parseEvidenceString(line));
        }
        scanner.close();
        return evidList;
    }

    private FirstEvidence parseEvidenceString(String line) throws PredicateNotFound{
        String[] predArr = line.split(REGEX_ESCAPE_CHAR + LEFTPRNTH);
        String symbolName = predArr[0];
        if(symbolName.contains(NOTOPERATOR))
            throw new PredicateNotFound("! oprator not allowed in evidence...");
        String[] predArr2 = predArr[1].split(EQUALSTO);
        String value = predArr2[1];
        int truthVal = Integer.parseInt(value);
        String[] termNames = predArr2[0].replace(RIGHTPRNTH, "").split(COMMASEPARATOR);

        ArrayList<Integer> values = new ArrayList<Integer>();
        for(String term : termNames){
            values.add(Integer.parseInt(term));
        }
        for(PredicateSymbol symbol : mln.symbols){
            if(symbolName.equals(symbol.symbol)){
                return new FirstEvidence(MLN.create_new_symbol(symbol),values,truthVal);
            }
        }
        System.out.println("line : "+line);
        throw new PredicateNotFound("wrong predicate in evidence");
    }

    public Evidence parseEvidence(GroundMLN groundMln, String evidence_file, Set<String> predsToParse) throws FileNotFoundException, PredicateNotFound {
        Evidence evidence = new Evidence();
        if(evidence_file == null)
            return evidence;
        Scanner scanner = new Scanner(new BufferedReader(new InputStreamReader(new FileInputStream(evidence_file))));
        while (scanner.hasNextLine()) {
            String line = scanner.nextLine().replaceAll("\\s", "");

            if (line.isEmpty()) {
                continue;
            }
            String[] predArr = line.split(REGEX_ESCAPE_CHAR + LEFTPRNTH);
            String symbolName = predArr[0];
            if(symbolName.contains(NOTOPERATOR))
                throw new PredicateNotFound("! oprator not allowed in evidence...");
            if(!predsToParse.contains(symbolName))
                continue;
            String[] predArr2 = predArr[1].split(EQUALSTO);
            String value = predArr2[1];
            String[] termNames = predArr2[0].replace(RIGHTPRNTH, "").split(COMMASEPARATOR);

            GroundPredicate gp = new GroundPredicate();
            // Find groundPredicate symbol with string=symbolname
            for(PredicateSymbol sym : mln.symbols)
            {
                if(sym.symbol.equals(symbolName))
                {
                    gp.symbol = new GroundPredicateSymbol(sym.id,symbolName, sym.values, sym.variable_types);
                    break;
                }
            }
            for(String term : termNames)
            {
                gp.terms.add(Integer.parseInt(term));
            }

//            MyAssert.assume(groundMln.groundPredToIntegerMap.containsKey(gp));
            if(!groundMln.groundPredToIntegerMap.containsKey(gp))
                continue;
            int gpIndex = groundMln.groundPredToIntegerMap.get(gp);
            if(gpIndex == -1)
                System.out.println("error in Parser::parseEvidence()");
            evidence.predIdVal.put(gpIndex, Integer.parseInt(value));
        }
        return evidence;
    }

}

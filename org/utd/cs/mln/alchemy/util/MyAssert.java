package org.utd.cs.mln.alchemy.util;

/**
 * Created by Happy on 4/16/18.
 */
public class MyAssert {
    public static void assume(boolean cond)
    {
        try{
            if(cond == false)
            {
                throw new AssertionError();
            }
        }
        catch(Exception e)
        {
            e.printStackTrace();
            System.exit(-1);
        }

    }
}

package org.utd.cs.mln.learning;

import org.utd.cs.mln.inference.GibbsSampler_v2;

/**
 * Created by piyush on 13/9/17.
 */
public class InferPerMLN implements Runnable {
    private GibbsSampler_v2 gs;
    private boolean burningIn, isInit;
    public InferPerMLN(GibbsSampler_v2 gs, boolean burningIn, boolean isInit){
        this.gs = gs;
        this.burningIn = burningIn;
        this.isInit = isInit;
    }

    @Override
    public void run() {
        gs.infer(burningIn, isInit);
    }
}

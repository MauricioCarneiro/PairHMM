package org.broadinstitute;

import edu.syr.pcpratts.rootbeer.runtime.Rootbeer;
import edu.syr.pcpratts.rootbeer.runtime.Kernel;

import java.util.List;
import java.util.ArrayList;
import java.util.Iterator;
import org.apache.log4j.Logger;

public class RootbeerEvaluate {

    private final double PRECISION;
    private static Logger logger = Logger.getLogger("Main");

    public RootbeerEvaluate(final double PRECISION){
        this.PRECISION = PRECISION;
    }

    public void runTests(final Iterator<Evaluate.TestRow> testCache) {
        final long startTime = System.currentTimeMillis();

        boolean testsPassed = true;

        List<Kernel> kernels = new ArrayList<Kernel>();
        int thread_counter = 0;
        int block_counter = 0;

        Rootbeer rootbeer = new Rootbeer();

        while (testCache.hasNext()) {
            Evaluate.TestRow currentTest = testCache.next();

            HmmKernel kernel = new HmmKernel(currentTest.getHaplotypeBases(), 
              currentTest.getReadBases(), currentTest.getReadQuals(), 
              currentTest.getReadInsQuals(), currentTest.getReadDelQuals(), 
              currentTest.getOverallGCP(), currentTest.getHaplotypeStart(), 
              currentTest.getReachedReadValye(), currentTest.getLikelihood());

            kernels.add(kernel);

            thread_counter++;
            if(thread_counter == 32){
                thread_counter = 0;
                block_counter++;
                if(block_counter == 14){
                    block_counter = 0;

                    rootbeer.setThreadConfig(32, 14);
                    rootbeer.runAll(kernels);

                    for(Kernel result_kernel : kernels){
                        HmmKernel result_hmm_kernel = (HmmKernel) result_kernel;
                        double result = result_hmm_kernel.getResult();
                        double likelihood = result_hmm_kernel.getLikelihood();

                        logger.debug(String.format(" Result:%4.3f",result));
                        logger.debug("==========================================================%n");
                        if (Math.abs(likelihood - result) > PRECISION) {
                            logger.error("Wrong result. Expected " + likelihood + " , actual: " + result);
                            testsPassed = false;
                        }
                    }
                    kernels.clear();
                }
            }
        }
        if (testsPassed) {
            logger.info(String.format("All tests PASSED in %.3f secs", (double) (System.currentTimeMillis() - startTime)/1000));
        }
    }
}

package org.broadinstitute;

import edu.syr.pcpratts.rootbeer.runtime.Rootbeer;
import edu.syr.pcpratts.rootbeer.runtime.Kernel;

import java.util.List;
import java.util.ArrayList;
import java.util.Iterator;

import java.io.FileWriter;
import java.io.IOException;

public class RootbeerEvaluate {

    private static Logger logger = new Logger();

    public RootbeerEvaluate(){
    }

    public void runTests(final Iterator<Evaluate.TestRow> testCache, String hmmName, String testSet, String filename) throws IOException {
        long totalTime = 0L;
        final long startTime = System.currentTimeMillis();

        final String runName = hmmName + "." + testSet;
        final FileWriter out = new FileWriter(filename);

        List<Kernel> kernels = new ArrayList<Kernel>();
        int thread_counter = 0;
        int block_counter = 0;

        Rootbeer rootbeer = new Rootbeer();

        while (testCache.hasNext()) {
            Evaluate.TestRow currentTest = testCache.next();

            HmmKernel kernel = new HmmKernel(currentTest.haplotypeBases, 
              currentTest.readBases, currentTest.readQuals, 
              currentTest.readInsQuals, currentTest.readDelQuals, 
              currentTest.overallGCP, currentTest.haplotypeStart, 
              currentTest.reachedReadValue);

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
                        double likelihood = result_hmm_kernel.getLikelihood();
                        out.write("" + likelihood + "\n");
                    }
                    kernels.clear();
                }
            }
        }

        long stopTime = System.currentTimeMillis();
        totalTime = stopTime - startTime;

        logger.info(String.format("%s test completed in %.3f secs, results written to: %s", hmmName, totalTime/1000000.0, filename));
        out.close();
    }
}

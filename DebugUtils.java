/**
 * @author carneiro
 * @since 8/16/12
 */
public class DebugUtils {
    public static void printMatrix(double[][] m) {
        for (double[] x : m) {
            for (double v : x) {
                System.out.print(v + " ");
            }
            System.out.println();
        }
        System.out.println();
    }
}

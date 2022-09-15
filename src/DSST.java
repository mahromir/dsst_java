import org.hipparchus.ode.nonstiff.AdaptiveStepsizeIntegrator;
import org.hipparchus.ode.nonstiff.DormandPrince853Integrator;
import org.orekit.data.*;


import org.orekit.errors.OrekitException;

import org.orekit.orbits.Orbit;
import org.orekit.orbits.OrbitType;
import org.orekit.propagation.SpacecraftState;
import org.orekit.propagation.analytical.tle.TLE;
import org.orekit.propagation.analytical.tle.TLEPropagator;
import org.orekit.propagation.numerical.NumericalPropagator;
import org.orekit.propagation.semianalytical.dsst.DSSTPropagator;
import org.orekit.time.AbsoluteDate;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Locale;


public class DSST {

    public static void main(String[] args) throws Exception {

        final File home = new File(System.getProperty("user.home"));
        final File orekitData = new File(home, "orekit-data");
        if (!orekitData.exists()) {
            System.err.format(Locale.US, "Failed to find %s folder%n", orekitData.getAbsolutePath());
            System.exit(1);
        }

        final DataProvidersManager manager = DataContext.getDefault().getDataProvidersManager();
        manager.addProvider(new DirectoryCrawler(orekitData));

        TLE myTle = new DSST().getTleFromInput();
        AbsoluteDate extrapDate = myTle.getDate();

        TLEPropagator tle_prop = TLEPropagator.selectExtrapolator(myTle);
        SpacecraftState scState = tle_prop.propagate(extrapDate);
        Orbit orbit = scState.getOrbit();

        NumericalPropagator myNM = DSST.createNumProp(orbit, 100.0);
        System.out.println(myNM);

    }

    private TLE getTleFromInput() throws IOException, OrekitException{
        String[] tle = {"", ""};

        File file = new File("C:/tles/iss/iss.txt");
        BufferedReader fr = new BufferedReader(new FileReader(file));
        String st;

        while ((st = fr.readLine()) != null){
            String[] splitedLine = st.split("");
            if (splitedLine[0].equals("1")) {
//                System.out.println(st);
                tle[0] = st;
            }
            if (splitedLine[0].equals("2")){
                tle[1] = st;
            }
        }
        return new TLE(tle[0], tle[1]);
    }

    private static NumericalPropagator createNumProp(final Orbit orbit, final double mass) {
        final double[][] tol = NumericalPropagator.tolerances(1.0, orbit, orbit.getType());
        System.out.println(tol.length);
        final double minStep = 1.e-3;
        final double maxStep = 1.e+3;
        final AdaptiveStepsizeIntegrator integrator = new DormandPrince853Integrator(minStep, maxStep, tol[0], tol[1]);
        integrator.setInitialStepSize(100.);

        final NumericalPropagator numProp = new NumericalPropagator(integrator);
        numProp.setInitialState(new SpacecraftState(orbit, mass));

        return numProp;
    }


}



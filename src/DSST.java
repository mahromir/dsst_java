import org.hipparchus.ode.nonstiff.AdaptiveStepsizeIntegrator;
import org.hipparchus.ode.nonstiff.DormandPrince853Integrator;
import org.orekit.bodies.CelestialBody;
import org.orekit.bodies.CelestialBodyFactory;
import org.orekit.bodies.OneAxisEllipsoid;
import org.orekit.data.*;


import org.orekit.errors.OrekitException;

import org.orekit.forces.gravity.potential.GravityFields;
import org.orekit.forces.gravity.potential.UnnormalizedSphericalHarmonicsProvider;
import org.orekit.frames.FramesFactory;
import org.orekit.frames.ITRFVersion;
import org.orekit.frames.VersionedITRF;
import org.orekit.models.earth.atmosphere.NRLMSISE00;
import org.orekit.models.earth.atmosphere.data.CssiSpaceWeatherData;
import org.orekit.orbits.Orbit;
import org.orekit.orbits.OrbitType;
import org.orekit.propagation.SpacecraftState;
import org.orekit.propagation.analytical.tle.TLE;
import org.orekit.propagation.analytical.tle.TLEPropagator;
import org.orekit.propagation.numerical.NumericalPropagator;
import org.orekit.propagation.semianalytical.dsst.DSSTPropagator;
import org.orekit.propagation.semianalytical.dsst.forces.DSSTAtmosphericDrag;
import org.orekit.propagation.semianalytical.dsst.forces.DSSTSolarRadiationPressure;
import org.orekit.propagation.semianalytical.dsst.forces.DSSTTesseral;
import org.orekit.propagation.semianalytical.dsst.forces.DSSTThirdBody;

import org.orekit.forces.gravity.potential.UnnormalizedSphericalHarmonicsProvider;
import org.orekit.forces.gravity.potential.GravityFieldFactory;

import org.orekit.time.AbsoluteDate;
import org.orekit.utils.Constants;
import org.orekit.utils.IERSConventions;

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

        AdaptiveStepsizeIntegrator myIntegrator = DSST.createInegrator(orbit);


        // Declaration of Bodies that are involved in simulation
        OneAxisEllipsoid Earth = new OneAxisEllipsoid(  Constants.WGS84_EARTH_EQUATORIAL_RADIUS,
                Constants.WGS84_EARTH_ANGULAR_VELOCITY,
                getItrf()
        );
        CelestialBody SunBody = CelestialBodyFactory.getSun();
        CelestialBody MoonBody = CelestialBodyFactory.getMoon();

        // Spaceweather with NRLMSISE atmoshphere model
        CssiSpaceWeatherData myWeatherData = new CssiSpaceWeatherData("SpaceWeather-All-v1.2.txt");
        NRLMSISE00 myMSISE = new NRLMSISE00(myWeatherData,SunBody, Earth);

        int degree = 20; int order = 20;
        final UnnormalizedSphericalHarmonicsProvider unnormalized_20x20 = GravityFieldFactory.getConstantUnnormalizedProvider(degree, order);

        final double cr = 1.4;
        final double cd = 2.2;
        final double cf = 1.2;
        final double meanCrossSection = 20.843234;
        final double mass = 1000.0;

        // Declaration of forcemodels
        DSSTTesseral DsstTessFM                 = new DSSTTesseral(getItrf(), Constants.WGS84_EARTH_ANGULAR_VELOCITY, unnormalized_20x20);
        DSSTSolarRadiationPressure DsstSR_FM    = new DSSTSolarRadiationPressure(cr, meanCrossSection, SunBody, Constants.WGS84_EARTH_EQUATORIAL_RADIUS, Constants.WGS84_EARTH_MU);
        DSSTThirdBody DsstMoonFM                = new DSSTThirdBody(MoonBody, Constants.JPL_SSD_MOON_GM);
        DSSTThirdBody DsstSunFM                 = new DSSTThirdBody(SunBody, Constants.JPL_SSD_SUN_GM);
        DSSTAtmosphericDrag DsstAtmFm           = new DSSTAtmosphericDrag(myMSISE, cd, meanCrossSection, Constants.WGS84_EARTH_MU);




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

    private static AdaptiveStepsizeIntegrator createInegrator(final Orbit orbit) {
        final double[][] tol = NumericalPropagator.tolerances(1.0, orbit, orbit.getType());
        System.out.println(tol.length);
        final double minStep = 1.e-3;
        final double maxStep = 1.e+3;
        final AdaptiveStepsizeIntegrator integrator = new DormandPrince853Integrator(minStep, maxStep, tol[0], tol[1]);
        integrator.setInitialStepSize(100.);

//        final NumericalPropagator numProp = new NumericalPropagator(integrator);
//        numProp.setInitialState(new SpacecraftState(orbit, mass));

        return integrator;
    }

    private static VersionedITRF getItrf(){
        return FramesFactory.getITRF(ITRFVersion.ITRF_2020, IERSConventions.IERS_2010, false);
    }

    //
    private static void createForceModels(DSSTTesseral DsstTessFM,
                                          DSSTThirdBody Moon,
                                          DSSTThirdBody Sun,
                                          DSSTAtmosphericDrag Atm,
                                          DSSTSolarRadiationPressure Srp){

        // Declaration of Bodies that are involved in simulation
        OneAxisEllipsoid Earth = new OneAxisEllipsoid(  Constants.WGS84_EARTH_EQUATORIAL_RADIUS,
                                                        Constants.WGS84_EARTH_ANGULAR_VELOCITY,
                                                        getItrf()
        );
        CelestialBody SunBody = CelestialBodyFactory.getSun();
        CelestialBody MoonBody = CelestialBodyFactory.getMoon();

        // Spaceweather with NRLMSISE atmoshphere model
        CssiSpaceWeatherData myWeatherData = new CssiSpaceWeatherData("SpaceWeather-All-v1.2.txt");
        NRLMSISE00 myMSISE = new NRLMSISE00(myWeatherData,SunBody, Earth);

        int degree = 20; int order = 20;
        final UnnormalizedSphericalHarmonicsProvider unnormalized_20x20 = GravityFieldFactory.getConstantUnnormalizedProvider(degree, order);

        // Declaration of forcemodels
//        DSSTTesseral DsstTessFM = new DSSTTesseral(getItrf(), Constants.WGS84_EARTH_ANGULAR_VELOCITY, unnormalized_20x20);







    }
}



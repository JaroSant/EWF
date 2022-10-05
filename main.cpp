#include "WrightFisher.h"
#include <boost/math/constants/constants.hpp>
#include <boost/math/distributions/beta.hpp>
#include <boost/math/distributions/binomial.hpp>
#include <boost/math/special_functions/beta.hpp> // for beta function
#include <boost/math/special_functions/binomial.hpp> // for function binomial_coefficient
#include <boost/math/special_functions/erf.hpp> // for erf function
#include <boost/timer/timer.hpp>
#include <boost/unordered/unordered_map.hpp>
#include <fstream>
#include <iostream>
#include <libconfig.h++>
#include <random>
#include <sstream>
#include <string>

using namespace std;

int main(int argc, char **argv) {
  if ((argc == 2) && (static_cast<string>(argv[1]) == "horses")) {
    cout << "================================" << endl;
    cout << "Demo invoked - allele frequency paths for horse coat coloration "
            "(ASIP locus) will be generated, using data from Ludwig et al. 2009"
         << endl;
    cout << "================================" << endl;

    vector<double> thetaP_in;
    bool non_neut_in, Absorption;
    double100 sigma_in, genGap;
    int SelPolyDeg_in, SelSetup_in, Ne_in, nEndpoints, nSim, nInterTimes;
    double dom_in;
    vector<double> selCoefs_in;

    libconfig::Config cfg;
    try {
      cfg.readFile("configHorseCoat.cfg");
    } catch (libconfig::FileIOException &e) {
      cerr << "FileIOException occurred. Could not read configHorseCoat.cfg!"
           << endl;
      return (EXIT_FAILURE);
    } catch (libconfig::ParseException &e) {
      cerr << "Parse error at " << e.getFile() << ":" << e.getLine() << "-"
           << e.getError() << endl;
      return (EXIT_FAILURE);
    }

    if ((cfg.lookupValue("g_entry", genGap)) && (genGap > 0.0)) {
      cout << "Generation age set to " << genGap << " years." << endl;
    } else {
      cerr
          << "Error in reading generation age - please check input is correct!";
      return (EXIT_FAILURE);
    }

    if (cfg.lookupValue("nonneutral_entry", non_neut_in)) {
      if (non_neut_in) {
        cout << "You have requested non-neutral paths." << endl;
        if (cfg.lookupValue("s_entry", sigma_in)) {
          cout << "Prelimiting selection coefficient set to " << sigma_in << "."
               << endl;
        } else {
          cerr << "Error in reading prelimiting selection coefficient - please "
                  "check input is correct!";
          return (EXIT_FAILURE);
        }

        if (cfg.lookupValue("selSetup_entry", SelSetup_in) &&
            (SelSetup_in == 0 || SelSetup_in == 1 || SelSetup_in == 2)) {
          if (SelSetup_in == 0) {
            cout << "You have chosen to have genic selection." << endl;
          } else if (SelSetup_in == 1) {
            cout << "you have chosen to have diploid selection. " << endl;
          } else {
            cout << "You have chosen to have polynomial selection." << endl;
          }
        } else {
          cerr << "Error in reading selection setup indicator function - "
                  "please check input is correct!";
          return (EXIT_FAILURE);
        }

        if (SelSetup_in == 1) {
          if ((cfg.lookupValue("dominance_entry", dom_in)) && (dom_in >= 0.0) &&
              (dom_in <= 1.0)) {
            cout << "You have set the dominance parameter to be " << dom_in
                 << "." << endl;
          } else {
            cerr << "Error in reading dominance parameter - please check input "
                    "is correct!";
            return (EXIT_FAILURE);
          }
        }

        if (SelSetup_in == 2) {
          if (cfg.lookupValue("polyDeg_entry", SelPolyDeg_in) &&
              (SelPolyDeg_in > 0)) {
            cout << "You have chosen a degree " << SelPolyDeg_in
                 << " polynomial as selection function." << endl;
          } else {
            cerr << "Error in reading degree of selection polynomial - please "
                    "check input is correct!";
            return (EXIT_FAILURE);
          }
        }
      } else {
        cout << "You have requested neutral paths." << endl;
      }
    } else {
      cerr << "Error in reading non-neutral indicator function - please check "
              "input is correct!";
      return (EXIT_FAILURE);
    }

    if (cfg.lookupValue("Ne_entry", Ne_in) && (Ne_in > 0)) {
      cout << "You have chosen an effective population size of " << Ne_in
           << endl;
    } else {
      cerr << "Error in reading effective population size - please check input "
              "is correct!";
      return (EXIT_FAILURE);
    }

    if (cfg.lookupValue("Absorption_entry", Absorption)) {
      if (Absorption) {
        cout << "You have chosen to not condition on non-absorption." << endl;
      } else {
        cout << "You have chosen to condition on non-absorption." << endl;
      }
    } else {
      cerr << "Error in reading absorption indicator function - please check "
              "input is correct!";
      return (EXIT_FAILURE);
    }

    if (cfg.lookupValue("nEndpoints", nEndpoints) && (nEndpoints > 0)) {
      cout << "You have chosen to have " << nEndpoints << " endpoints." << endl;
    } else {
      cerr << "Error in number of endpoints - please check input is correct!";
      return (EXIT_FAILURE);
    }

    if (cfg.lookupValue("nSim_entry", nSim) && (nSim > 0)) {
      cout << "You have chosen to run " << nSim << " simulations." << endl;
    } else {
      cerr << "Error in reading number of simulations - please check input is "
              "correct!";
      return (EXIT_FAILURE);
    }

    if (cfg.lookupValue("nInterTimes", nInterTimes) && (nInterTimes > 0)) {
      cout << "You have chosen to simulate " << nInterTimes
           << " intermediate time points." << endl;
    } else {
      cerr << "Error in reading number of intermediate times to sample - "
              "please check input is correct!";
      return (EXIT_FAILURE);
    }

    const libconfig::Setting &root = cfg.getRoot();
    if (root["mu_entry"].getLength() > 2) {
      cerr << "Mutation vector should contain only two entries!";
    }
    for (int i = 0; i < root["mu_entry"].getLength(); i++) {
      thetaP_in.push_back(2.0 * static_cast<double>(Ne_in) *
                          ((double)root["mu_entry"][i]));
      if (i == 0 && (thetaP_in.back() > 0.0)) {
        cout << "theta_1 = " << thetaP_in.back() << "." << endl;
      } else if (thetaP_in.back() > 0.0) {
        cout << "theta_2 = " << thetaP_in.back() << "." << endl;
      }
    }
    if (!(thetaP_in[0] > 0.0) && !(thetaP_in[1] > 0.0)) {
      thetaP_in.clear();
      cout << "(theta_1,theta_2) = (0,0)." << endl;
    }

    if (SelSetup_in == 2) {
      if (root["polyCoeffs_entries"].getLength() != SelPolyDeg_in) {
        cerr << "Mismatch in given selection coefficient vector and specified "
                "degree of polynomial!";
        return (EXIT_FAILURE);
      }
      cout << "The inputted selection polynomial has: " << endl;
      for (int i = 0; i < root["polyCoeffs_entries"].getLength(); i++) {
        selCoefs_in.push_back((double)root["polyCoeffs_entries"][i]);
        cout << " a coefficient " << selCoefs_in.back() << " for the "
             << root["polyCoeffs_entries"].getLength() - 1 - i
             << "-th power of x." << endl;
      }
    }

    if (non_neut_in) {
      sigma_in = 2.0 * static_cast<double>(Ne_in) * sigma_in;
      cout << "With the above input, the population rescaled selection "
              "coefficient sigma is "
           << sigma_in << endl;
    }
    cout << "================================" << endl;

    WrightFisher Horses(thetaP_in, non_neut_in, sigma_in, SelSetup_in, dom_in,
                        SelPolyDeg_in, selCoefs_in);

    Options o;
    boost::random::mt19937 gen;

    if (!(Absorption)) {
      Horses.ThetaResetter();
    }

    if ((root["observationTimes_entry"].getLength() !=
         root["observationSamples_entry"].getLength()) ||
        (root["observationTimes_entry"].getLength() !=
         root["observationCount_entry"].getLength()) ||
        (root["observationSamples_entry"].getLength() !=
         root["observationCount_entry"].getLength())) {
      cerr << "Mismatch in the inputted observation times, samples and counts!"
           << endl;
      return (EXIT_FAILURE);
    }

    vector<double100> obsTimes;
    vector<int> obsSamples, obsCount;

    for (int i = 0; i < nEndpoints; i++) {
      obsTimes.push_back(root["observationTimes_entry"][i]);
      obsSamples.push_back(root["observationSamples_entry"][i]);
      obsCount.push_back(root["observationCount_entry"][i]);
    }

    vector<double100> obsTimesDiff, obsFrequencies;
    vector<int>::iterator oSi = obsSamples.begin(), oCi = obsCount.begin();
    for (vector<double100>::iterator yearIter = obsTimes.begin();
         yearIter != obsTimes.end(); yearIter++, oSi++, oCi++) {
      obsTimesDiff.push_back(
          static_cast<double100>(abs(*yearIter - obsTimes.front())) /
          (2.0 * static_cast<double100>(Ne_in) * genGap));
      obsFrequencies.push_back(static_cast<double100>(*oCi) /
                               static_cast<double100>(*oSi));
    }

    ofstream originalTrajectoryFile, originalTimesFile;
    string originalTrajectoryFilename = "OGHT.txt",
           originalTimesFilename = "OGT.txt";
    originalTrajectoryFile.open(originalTrajectoryFilename);
    originalTimesFile.open(originalTimesFilename);

    for (vector<double100>::iterator tit = obsTimesDiff.begin(),
                                     xit = obsFrequencies.begin();
         tit != obsTimesDiff.end(); tit++, xit++) {
      originalTrajectoryFile << *xit << " ";
      originalTimesFile << *tit << " ";
      if (xit == obsFrequencies.end() - 1) {
        originalTrajectoryFile << "\n";
        originalTimesFile << "\n";
      }
    }

    originalTrajectoryFile.close();
    originalTimesFile.close();

    vector<double100> skeletonPoints, skeletonTimes, samplingTimes,
        samplingValues;

    for (int j = 0; j != nInterTimes; j++) {
      samplingTimes.push_back(
          obsTimesDiff.back() *
          (static_cast<double100>(j) / static_cast<double100>(nInterTimes)));
    }

    ofstream imputedTrajectoryFile, imputedTimesFile;
    string imputedTrajectoryFilename = "ImpHT.txt",
           imputedTimesFilename = "ImpT.txt";
    imputedTrajectoryFile.open(imputedTrajectoryFilename);
    imputedTimesFile.open(imputedTimesFilename);

    for (int i = 0; i != nSim; i++) {
      vector<double100>::iterator leftTime = samplingTimes.begin(), rightTime;
      for (vector<double100>::iterator tIt = obsTimesDiff.begin(),
                                       xIt = obsFrequencies.begin();
           tIt != (obsTimesDiff.end() - 1); tIt++, xIt++) {
        vector<vector<double100>> Skeleton = Horses.NonNeutralDrawBridge(
            *xIt, *tIt, *(tIt + 1), *(xIt + 1), Absorption, o, gen);
        vector<double100> skeletonPoints(Skeleton[0]),
            skeletonTimes(Skeleton[1]);
        skeletonPoints.push_back(*(xIt + 1));
        skeletonTimes.push_back(*(tIt + 1));

        vector<double100>::iterator timeRunner = leftTime;
        if (tIt != obsTimesDiff.end() - 2) {
          while (*timeRunner < *(tIt + 1)) {
            timeRunner++;
          }

          rightTime = timeRunner - 1;
        } else {
          rightTime = samplingTimes.end() - 1;
        }

        vector<double100>::iterator skeltimeIt = skeletonTimes.begin(),
                                    skelvalIt = skeletonPoints.begin();
        for (vector<double100>::iterator timeIt = leftTime;
             timeIt != rightTime - 1; timeIt++) {
          if (timeIt == leftTime) {
            samplingValues.push_back(*xIt);
            imputedTrajectoryFile << samplingValues.back() << " ";
            imputedTimesFile << *tIt << " ";
          } else {
            if (*(timeIt + 1) > *(skeltimeIt)) {
              skeltimeIt++;
              skelvalIt++;
            }
            samplingValues.push_back(
                Horses
                    .DrawBridgepoint(samplingValues.back(), *skelvalIt, *timeIt,
                                     *skeltimeIt, *(timeIt + 1), o, gen)
                    .first);
            imputedTrajectoryFile << samplingValues.back() << " ";
            imputedTimesFile << *timeIt << " ";
          }
        }

        leftTime = rightTime;
      }

      imputedTrajectoryFile << obsFrequencies.back() << "\n";
      imputedTimesFile << samplingTimes.back() << "\n";

      cout << "path " << i << " ready!" << endl;
    }

    imputedTrajectoryFile.close();
    imputedTimesFile.close();
    cout << "================================" << endl;
    cout << "Simulations complete - the generated paths can be found in "
            "'ImpHT.txt', whilst 'ImpT.txt' contains the corresponding time "
            "stamps."
         << endl;
    cout << "The files 'OGT.txt' and 'OGHT.txt' contain the observation times "
            "and observed frequencies respectively."
         << endl;
    cout << "================================" << endl;

  } else {
    vector<double> thetaP_in;
    bool non_neut_in;
    double100 sigma_in;
    int SelPolyDeg_in, SelSetup_in;
    double dom_in;
    vector<double> selCoefs_in;

    libconfig::Config cfg;
    try {
      cfg.readFile("config.cfg");
    } catch (libconfig::FileIOException &e) {
      cerr << "FileIOException occurred. Could not read config.cfg!" << endl;
      return (EXIT_FAILURE);
    } catch (libconfig::ParseException &e) {
      cerr << "Parse error at " << e.getFile() << ":" << e.getLine() << "-"
           << e.getError() << endl;
      return (EXIT_FAILURE);
    }

    if (cfg.lookupValue("nonneutral_entry", non_neut_in)) {
      if (non_neut_in) {
        cout << "You have requested non-neutral paths." << endl;
        if ((cfg.lookupValue("sigma_entry", sigma_in))) {
          cout << "You have set the population rescaled selection parameter to "
                  "be "
               << sigma_in << "." << endl;
        } else {
          cerr << "Error in reading population rescaled selection parameter - "
                  "please check input is correct!";
          return (EXIT_FAILURE);
        }

        if (cfg.lookupValue("selSetup_entry", SelSetup_in) &&
            (SelSetup_in == 0 || SelSetup_in == 1 || SelSetup_in == 2)) {
          if (SelSetup_in == 0) {
            cout << "You have chosen to have genic selection." << endl;
          } else if (SelSetup_in == 1) {
            cout << "you have chosen to have diploid selection. " << endl;
          } else {
            cout << "You have chosen to have polynomial selection." << endl;
          }
        } else {
          cerr << "Error in reading selection setup indicator function - "
                  "please check input is correct!";
          return (EXIT_FAILURE);
        }

        if (SelSetup_in == 1) {
          if ((cfg.lookupValue("dominance_entry", dom_in)) && (dom_in >= 0.0) &&
              (dom_in <= 1.0)) {
            cout << "You have set the dominance parameter to be " << dom_in
                 << "." << endl;
          } else {
            cerr << "Error in reading dominance parameter - please check input "
                    "is correct!";
            return (EXIT_FAILURE);
          }
        }

        if (SelSetup_in == 2) {
          if (cfg.lookupValue("polyDeg_entry", SelPolyDeg_in) &&
              (SelPolyDeg_in > 0)) {
            cout << "You have chosen a degree " << SelPolyDeg_in
                 << " polynomial as selection function." << endl;
          } else {
            cerr << "Error in reading degree of selection polynomial - please "
                    "check input is correct!";
            return (EXIT_FAILURE);
          }
        }
      } else {
        cout << "You have requested neutral paths." << endl;
      }
    } else {
      cerr << "Error in reading non-neutral indicator function - please check "
              "input is correct!";
      return (EXIT_FAILURE);
    }

    const libconfig::Setting &root = cfg.getRoot();
    if (root["theta_entries"].getLength() > 2) {
      cerr << "Mutation vector should contain only two entries!";
    }
    for (int i = 0; i < root["theta_entries"].getLength(); i++) {
      thetaP_in.push_back((double)root["theta_entries"][i]);
      if (i == 0 && (thetaP_in.back() > 0.0)) {
        cout << "theta_1 = " << thetaP_in.back() << "." << endl;
      } else if (thetaP_in.back() > 0.0) {
        cout << "theta_2 = " << thetaP_in.back() << "." << endl;
      }
    }
    if (!(thetaP_in[0] > 0.0) && !(thetaP_in[1] > 0.0)) {
      thetaP_in.clear();
      cout << "(theta_1,theta_2) = (0,0)." << endl;
    }

    if (SelSetup_in == 2) {
      if (root["polyCoeffs_entries"].getLength() != SelPolyDeg_in) {
        cerr << "Mismatch in given selection coefficient vector and specified "
                "degree of polynomial!";
        return (EXIT_FAILURE);
      }
      cout << "The inputted selection polynomial has: " << endl;
      for (int i = 0; i < root["polyCoeffs_entries"].getLength(); i++) {
        selCoefs_in.push_back((double)root["polyCoeffs_entries"][i]);
        cout << " a coefficient " << selCoefs_in.back() << " for the "
             << root["polyCoeffs_entries"].getLength() - 1 - i
             << "-th power of x." << endl;
      }
    }

    WrightFisher test(thetaP_in, non_neut_in, sigma_in, SelSetup_in, dom_in,
                      SelPolyDeg_in, selCoefs_in);

    Options o;
    boost::random::mt19937 gen;

    cout << "================================" << endl;
    cout << "Please choose whether you wish to simulate draws from the "
            "Wright-Fisher diffusion, or from the diffusion bridge."
         << endl;
    cout << "Enter: 1 for diffusion, 2 for diffusion bridge" << endl;
    cout << "================================" << endl;
    int DiffOrBridge;
    cin >> DiffOrBridge;
    bool stop =
        (cin.fail()
             ? false
             : ((DiffOrBridge != 1 && DiffOrBridge != 2) ? false : true));

    while (!stop) {
      cout << "Please enter either 1 for diffusion simulation or 2 for "
              "diffusion bridge simulation!"
           << endl;
      cin.clear();
      cin.ignore(256, '\n');
      cin >> DiffOrBridge;
      stop = (cin.fail()
                  ? false
                  : ((DiffOrBridge != 1 && DiffOrBridge != 2) ? false : true));
    }

    if (DiffOrBridge == 1) {
      cout << "================================" << endl;
      cout << "You have chosen to simulate a diffusion!" << endl;
      cout << "================================" << endl;
      cout << "Would you like to additionally compute the truncated transition "
              "density?"
           << endl;
      cout << "Please enter 1 for yes, 2 for no." << endl;
      cout << "================================" << endl;

      int Density;
      cin >> Density;
      bool stopyn =
          (cin.fail() ? false
                      : ((Density != 1 && Density != 2) ? false : true));

      while (!stopyn) {
        cout << "Please enter either 1 for yes or 2 for no!" << endl;
        cin.clear();
        cin.ignore(256, '\n');
        cin >> Density;
        stopyn = (cin.fail() ? false
                             : ((Density != 1 && Density != 2) ? false : true));
      }

      vector<double100> startPoints, startTimes, sampleTimes;
      vector<int> nSim, meshSize;

      libconfig::Config cfgDiff;
      try {
        cfgDiff.readFile("configDiffusion.cfg");
      } catch (libconfig::FileIOException &e) {
        cerr << "FileIOException occurred. Could not read configDiffusion.cfg!"
             << endl;
        return (EXIT_FAILURE);
      } catch (libconfig::ParseException &e) {
        cerr << "Parse error at " << e.getFile() << ":" << e.getLine() << "-"
             << e.getError() << endl;
        return (EXIT_FAILURE);
      }

      bool Absorption;
      if (cfgDiff.lookupValue("Absorption_entry", Absorption)) {
        if (Absorption) {
          cout << "You have chosen to not condition on non-absorption." << endl;
        } else {
          cout << "You have chosen to condition on non-absorption." << endl;
        }
      } else {
        cerr << "Error in reading absorption indicator function - please check "
                "input is correct!";
        return (EXIT_FAILURE);
      }

      const libconfig::Setting &root = cfgDiff.getRoot();
      int xlen = root["startDiff_entry"].getLength(),
          tlen = root["startDiffTime_entry"].getLength(),
          slen = root["sampleDiffTime_entry"].getLength();
      int nlen = root["nSim_entry"].getLength(),
          mlen = root["meshSize_entry"].getLength();

      if ((xlen == tlen) && (tlen == slen) && (slen == nlen)) {
        if (Density == 1) {
          if (mlen != xlen) {
            cout << "Mismatch in the configuration file input!" << endl;
            if (mlen < xlen) {
              cout << "There are too few mesh size entries!" << endl;
            } else {
              cout << "There are too many mesh size entries!" << endl;
            }
            cerr << "Simulation aborted - please fix configDiffusion.cfg as "
                    "per the above suggestions."
                 << endl;
            return (EXIT_FAILURE);
          }
        }
        cout << "Reading in diffusion simulation configuration from "
                "configDiffusion.cfg"
             << endl;
        cout << "================================" << endl;
      } else {
        cout << "There is a mismatch in the configuration file input!" << endl;
        if (xlen < tlen) {
          cout << "There are more diffusion start points than start times!"
               << endl;
        }
        if (xlen > tlen) {
          cout << "There are more diffusion start times than start points!"
               << endl;
        }
        if (slen < tlen) {
          cout << "There are more diffusion start times than sample times!"
               << endl;
        }
        if (slen > tlen) {
          cout << "There are more diffusion sample times than start times!"
               << endl;
        }
        if (xlen < slen) {
          cout << "There are more diffusion start points than sample times!"
               << endl;
        }
        if (xlen > slen) {
          cout << "There are more diffusion sample times than start points!"
               << endl;
        }
        if (xlen < nlen) {
          cout
              << "There are more number of samples than diffusion start points!"
              << endl;
        }
        if (xlen > nlen) {
          cout
              << "There are more diffusion start points than number of samples!"
              << endl;
        }
        if (slen < nlen) {
          cout
              << "There are more number of samples than diffusion sample times!"
              << endl;
        }
        if (slen > nlen) {
          cout
              << "There are more diffusion sample times than number of samples!"
              << endl;
        }
        if (tlen < nlen) {
          cout << "There are more number of samples than diffusion start times!"
               << endl;
        }
        if (tlen > nlen) {
          cout << "There are more diffusion start times than number of samples!"
               << endl;
        }
        cout << "Simulation aborted - please fix configDiffusion.cfg as per "
                "the above suggestions."
             << endl;
        return (EXIT_FAILURE);
      }

      cout << "================================" << endl;
      cout << "Diffusion input read in without errors." << endl;
      cout << "================================" << endl;

      for (int i = 0; i < root["startDiff_entry"].getLength(); i++) {
        startPoints.push_back(root["startDiff_entry"][i]);
        startTimes.push_back(root["startDiffTime_entry"][i]);
        sampleTimes.push_back(root["sampleDiffTime_entry"][i]);
        nSim.push_back(root["nSim_entry"][i]);
        if (Density == 1) {
          meshSize.push_back(root["meshSize_entry"][i]);
        }
      }

      cout << "You have chosen to perform the following simulations:" << endl;
      int i = 1;
      vector<int>::iterator nit = nSim.begin();
      for (vector<double100>::iterator xit = startPoints.begin(),
                                       tit = startTimes.begin(),
                                       sit = sampleTimes.begin();
           xit != startPoints.end(); xit++, tit++, sit++, nit++) {
        cout << "Diffusion " << i << " : x = " << *xit << ", t_0 = " << *tit
             << ", t = " << *sit << ", nSim = " << *nit << endl;
        i++;
      }

      cout << "================================" << endl;
      nit = nSim.begin();
      vector<int>::iterator mit = meshSize.begin();
      for (vector<double100>::iterator xit = startPoints.begin(),
                                       tit = startTimes.begin(),
                                       sit = sampleTimes.begin();
           xit != startPoints.end(); xit++, tit++, sit++, nit++, mit++) {
        cout << " Simulating diffusion draws started at x " << *xit
             << " at time " << *tit << ", sampled at time " << *sit << endl;
        time_t t = time(0); // get time now
        struct tm *now = localtime(&t);

        char bufferSaveFile[80];
        strftime(bufferSaveFile, 80, "%Y-%m-%d-%H-%M", now);
        string absnoabs;
        if (Absorption) {
          absnoabs = "Unconditioned";
        } else {
          absnoabs = "Conditioned";
        }
        string saveFilename =
            static_cast<string>(bufferSaveFile) + absnoabs + "Diffusion" +
            "SamplesX" +
            boost::lexical_cast<std::string>(static_cast<int>(100.0 * (*xit))) +
            "T" +
            boost::lexical_cast<std::string>(static_cast<int>(100.0 * (*tit))) +
            "S" +
            boost::lexical_cast<std::string>(static_cast<int>(100.0 * (*sit))) +
            ".txt";
        cout << "Output will be saved in the file " << saveFilename << endl;

        test.DiffusionRunner(*nit, *xit, *tit, *sit, Absorption, saveFilename,
                             o, gen);

        if (Density == 1) {
          string densityFilename = static_cast<string>(bufferSaveFile) +
                                   absnoabs + "Diffusion" + "DensityX" +
                                   boost::lexical_cast<std::string>(
                                       static_cast<int>(100.0 * (*xit))) +
                                   "T" +
                                   boost::lexical_cast<std::string>(
                                       static_cast<int>(100.0 * (*tit))) +
                                   "S" +
                                   boost::lexical_cast<std::string>(
                                       static_cast<int>(100.0 * (*sit))) +
                                   ".txt";
          cout << "Truncated density will be saved in the file "
               << densityFilename << endl;
          test.DiffusionDensityCalculator(*mit, *xit, *tit, *sit, Absorption,
                                          densityFilename, o);
        }
        cout << "================================" << endl;
      }

    } else {
      cout << "================================" << endl;
      cout << "You have chosen to simulate a diffusion bridge!" << endl;
      cout << "Would you like to additionally compute the truncated transition "
              "density?"
           << endl;
      cout << "Please enter 1 for yes, 2 for no." << endl;
      cout << "================================" << endl;

      int Density;
      cin >> Density;
      bool stopyn =
          (cin.fail() ? false
                      : ((Density != 1 && Density != 2) ? false : true));

      while (!stopyn) {
        cout << "Please enter either 1 for yes or 2 for no!" << endl;
        cin.clear();
        cin.ignore(256, '\n');
        cin >> Density;
        stopyn = (cin.fail() ? false
                             : ((Density != 1 && Density != 2) ? false : true));
      }

      libconfig::Config cfgBridge;
      try {
        cfgBridge.readFile("configBridge.cfg");
      } catch (libconfig::FileIOException &e) {
        cerr << "FileIOException occurred. Could not read configBridge.cfg!"
             << endl;
        return (EXIT_FAILURE);
      } catch (libconfig::ParseException &e) {
        cerr << "Parse error at " << e.getFile() << ":" << e.getLine() << "-"
             << e.getError() << endl;
        return (EXIT_FAILURE);
      }

      bool Absorption;
      if (cfgBridge.lookupValue("Absorption_entry", Absorption)) {
        if (Absorption) {
          cout << "You have chosen to not condition on non-absorption." << endl;
        } else {
          cout << "You have chosen to condition on non-absorption." << endl;
        }
      } else {
        cerr << "Error in reading absorption indicator function - please check "
                "input is correct!";
        return (EXIT_FAILURE);
      }

      const libconfig::Setting &root = cfgBridge.getRoot();
      int nBridges = root["nEndpoints"].getLength(),
          xlen = root["bridgePoints_entry"].getLength(),
          tlen = root["bridgeTimes_entry"].getLength();
      int nSamples = root["nSampleTimes_entry"].getLength(),
          slen = root["sampleTimes_entry"].getLength(),
          nlen = root["nSim_entry"].getLength(),
          mlen = root["meshSize_entry"].getLength();
      int nBridgeChecker = -1;

      if ((xlen == tlen) && (tlen == slen) && (slen == nlen)) {
        if (Density == 1) {
          if (mlen != xlen) {
            cout << "Mismatch in the configuration file input!" << endl;
            if (mlen < xlen) {
              cout << "There are too few mesh size entries!" << endl;
            } else {
              cout << "There are too many mesh size entries!" << endl;
            }
            cerr << "Simulation aborted - please fix configBridge.cfg as per "
                    "the above suggestions."
                 << endl;
            return (EXIT_FAILURE);
          }
        }
        cout << "Reading in bridge simulation configuration from "
                "configBridge.cfg"
             << endl;
        cout << "================================" << endl;
      } else {
        cout << "There is a mismatch in the configuration file input!" << endl;

        if (xlen < tlen) {
          cout << "There are more diffusion start points than start times!"
               << endl;
        }
        if (xlen > tlen) {
          cout << "There are more diffusion start times than start points!"
               << endl;
        }
        if (slen < tlen) {
          cout << "There are more diffusion start times than sample times!"
               << endl;
        }
        if (slen > tlen) {
          cout << "There are more diffusion sample times than start times!"
               << endl;
        }
        if (xlen < slen) {
          cout << "There are more diffusion start points than sample times!"
               << endl;
        }
        if (xlen > slen) {
          cout << "There are more diffusion sample times than start points!"
               << endl;
        }
        if (xlen < nlen) {
          cout
              << "There are more number of samples than diffusion start points!"
              << endl;
        }
        if (xlen > nlen) {
          cout
              << "There are more diffusion start points than number of samples!"
              << endl;
        }
        if (slen < nlen) {
          cout
              << "There are more number of samples than diffusion sample times!"
              << endl;
        }
        if (slen > nlen) {
          cout
              << "There are more diffusion sample times than number of samples!"
              << endl;
        }
        if (tlen < nlen) {
          cout << "There are more number of samples than diffusion start times!"
               << endl;
        }
        if (tlen > nlen) {
          cout << "There are more diffusion start times than number of samples!"
               << endl;
        }
        cout << "Simulation aborted - please fix configDiffusion.cfg as per "
                "the above suggestions."
             << endl;
        return (EXIT_FAILURE);
      }

      if ((nBridges < 1) || (nSamples < 1)) {
        cerr << "No bridge or sampling information provided! Please amend "
                "configBridge.cfg appropriately!";
        return (EXIT_FAILURE);
      } else {
        nBridgeChecker = 0;
        for (int i = 0; i < nBridges; i++) {
          nBridgeChecker += (int)root["nEndpoints"][i];
        }

        if (nBridgeChecker != nBridges + nSamples) {
          cerr << "Mismatch in nEndpoints and nSampleTimes_entry input in "
                  "configBridge.cfg! Please consult config file for info on "
                  "how to set these two quantities up!";
          return (EXIT_FAILURE);
        }
      }

      cout << "================================" << endl;
      cout << "Bridge input read in without errors." << endl;
      cout << "================================" << endl;

      vector<double100> bridgePoints, bridgeTimes, sampleTimes;
      vector<int> nSim, meshSize, nEndpoints, nSampleTimes;

      for (int i = 0; i < nBridges; i++) {
        nEndpoints.push_back(root["nEndpoints"][i]);
      }

      for (int i = 0; i < nSamples; i++) {
        nSampleTimes.push_back(root["nSampleTimes_entry"][i]);
      }

      for (int i = 0; i < xlen; i++) {
        bridgePoints.push_back(root["bridgePoints_entry"][i]);
        bridgeTimes.push_back(root["bridgeTimes_entry"][i]);
      }
      for (int i = 0; i < slen; i++) {
        sampleTimes.push_back(root["sampleTimes_entry"][i]);
        nSim.push_back(root["nSim_entry"][i]);
        meshSize.push_back(root["meshSize_entry"][i]);
      }

      cout << "You have chosen to perform the following simulations:" << endl;
      vector<double100>::iterator indexbP = bridgePoints.begin();
      vector<double100>::iterator indexbT = bridgeTimes.begin();
      vector<int>::iterator indexnS = nSampleTimes.begin();
      vector<double100>::iterator indexsT = sampleTimes.begin();
      int counter = 1;
      for (vector<int>::iterator nEp = nEndpoints.begin();
           nEp != nEndpoints.end(); nEp++) {
        cout << "Bridge " << counter << " : " << endl;
        vector<double100> brPt(indexbP, indexbP + (*nEp)),
            brTs(indexbT, indexbT + (*nEp));
        vector<int> nsT(indexnS, indexnS + (*nEp) - 1);
        indexbP += *nEp;
        indexbT += *nEp;
        indexnS += (*nEp) - 1;
        vector<int>::iterator nsti = nsT.begin();
        for (vector<double100>::iterator bP = brPt.begin(), bT = brTs.begin();
             bP != brPt.end() - 1; bP++, bT++) {
          vector<double100> nsTs(indexsT, indexsT + *nsti);
          indexsT += *nsti;
          for (vector<double100>::iterator sT = nsTs.begin(); sT != nsTs.end();
               sT++) {
            cout << "x = " << *bP << ", z = " << *(bP + 1) << ", t1 = " << *bT
                 << ", t2 = " << *(bT + 1) << ", s = " << *sT << endl;
          }
          nsti++;
        }
        counter++;
      }

      indexbP = bridgePoints.begin();
      indexbT = bridgeTimes.begin();
      indexnS = nSampleTimes.begin();
      indexsT = sampleTimes.begin();
      vector<int>::iterator nS = nSim.begin(), nM = meshSize.begin();
      for (vector<int>::iterator nEp = nEndpoints.begin();
           nEp != nEndpoints.end(); nEp++) {
        vector<double100> brPt(indexbP, indexbP + (*nEp)),
            brTs(indexbT, indexbT + (*nEp));
        vector<int> nsT(indexnS, indexnS + (*nEp) - 1);
        indexbP += *nEp;
        indexbT += *nEp;
        indexnS += (*nEp) - 1;
        vector<int>::iterator nsti = nsT.begin();
        for (vector<double100>::iterator bP = brPt.begin(), bT = brTs.begin();
             bP != brPt.end() - 1; bP++, bT++) {
          vector<double100> nsTs(indexsT, indexsT + *nsti);
          indexsT += *nsti;
          for (vector<double100>::iterator sT = nsTs.begin(); sT != nsTs.end();
               sT++) {
            cout << "================================" << endl;
            cout << "Simulating diffusion bridges draws for x = " << *bP
                 << ", z = " << *(bP + 1) << ", t1 = " << *bT
                 << ", t2 = " << *(bT + 1) << ", s = " << *sT << endl;
            time_t t = time(0); // get time now
            struct tm *now = localtime(&t);

            char bufferSaveFile[80];
            strftime(bufferSaveFile, 80, "%Y-%m-%d-%H-%M", now);
            string absnoabs;
            if (Absorption) {
              absnoabs = "Unconditioned";
            } else {
              absnoabs = "Conditioned";
            }
            string saveFilename = static_cast<string>(bufferSaveFile) +
                                  absnoabs + "Bridge" + "SamplesX" +
                                  boost::lexical_cast<std::string>(
                                      static_cast<int>(100.0 * (*bP))) +
                                  "Z" +
                                  boost::lexical_cast<std::string>(
                                      static_cast<int>(100.0 * (*(bP + 1)))) +
                                  "T1" +
                                  boost::lexical_cast<std::string>(
                                      static_cast<int>(100.0 * (*bT))) +
                                  "T2" +
                                  boost::lexical_cast<std::string>(
                                      static_cast<int>(100.0 * (*(bT + 1)))) +
                                  "S" +
                                  boost::lexical_cast<std::string>(
                                      static_cast<int>(100.0 * (*sT))) +
                                  ".txt";
            cout << "Output will be saved in the file " << saveFilename << endl;

            test.BridgeDiffusionRunner(*nS, *bP, *(bP + 1), *bT, *(bT + 1), *sT,
                                       Absorption, saveFilename, o, gen);
            nS++;

            if (Density == 1) {
              string densityFilename =
                  static_cast<string>(bufferSaveFile) + absnoabs + "Bridge" +
                  "DensityX" +
                  boost::lexical_cast<std::string>(
                      static_cast<int>(100.0 * (*bP))) +
                  "Z" +
                  boost::lexical_cast<std::string>(
                      static_cast<int>(100.0 * (*(bP + 1)))) +
                  "T1" +
                  boost::lexical_cast<std::string>(
                      static_cast<int>(100.0 * (*bT))) +
                  "T2" +
                  boost::lexical_cast<std::string>(
                      static_cast<int>(100.0 * (*(bT + 1)))) +
                  "S" +
                  boost::lexical_cast<std::string>(
                      static_cast<int>(100.0 * (*sT))) +
                  ".txt";
              cout << "Truncated density will be saved in the file "
                   << densityFilename << endl;
              test.BridgeDiffusionDensityCalculator(*nM, *bP, *(bP + 1), *bT,
                                                    *(bT + 1), *sT, Absorption,
                                                    densityFilename, o);
            }
            cout << "================================" << endl;
          }
          nsti++;
        }
      }
    }
  }

  return 0;
}

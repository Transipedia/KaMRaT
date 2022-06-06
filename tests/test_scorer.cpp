#include <string>
#include <vector>
#include <iostream>

#include "lest.hpp"
#include "scorer.hpp"

using namespace std;


#define VECT_SIZE 200

const lest::test module[] =
{
    CASE( "Comparison test scores after refactoring" ) {
        srand(time(NULL));
        vector<float> v;

        // Randomly fill the test vector
        for (uint i=0 ; i<VECT_SIZE ; i++) {
            float x = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX/9997.0));
            v.push_back(x);
        }

        // Scorer binary conditions
        vector<string> headers;
        for (uint i=0 ; i<VECT_SIZE ; i++) {
            headers.push_back(to_string(i % 2));
        }

        // Ttest
        SETUP( "Ttest" ) {
            cout << "Ttest" << endl;
            // Create the scorer
            Scorer scorer("ttest.padj", 0, headers);

            double res_old = scorer.EstimateScore_old(v);
            double res_new = scorer.EstimateScore(v);
            EXPECT( abs(res_old - res_new) < 1.0/pow(10, 10) );
            cout << "   ok" << endl;
        }

        // Ttest
        SETUP( "Ttest pi value" ) {
            cout << "Ttest pi value" << endl;
            // Create the scorer
            Scorer scorer("ttest.pi", 0, headers);

            double res_old = scorer.EstimateScore_old(v);
            double res_new = scorer.EstimateScore(v);
            EXPECT( abs(res_old - res_new) < 1.0/pow(10, 10) );
            cout << "   ok" << endl;
        }

        // SNR
        SETUP( "Signal to noise" ) {
            cout << "Signal to noise" << endl;
            // Create the scorer
            Scorer scorer("snr", 0, headers);

            double res_old = scorer.EstimateScore_old(v);
            double res_new = scorer.EstimateScore(v);
            EXPECT( abs(res_old - res_new) < 1.0/pow(10, 10) );
            cout << "   ok" << endl;
        }

        // Dids score
        SETUP( "DIDS score" ) {
            cout << "DIDS score" << endl;
            // Create the scorer
            Scorer scorer("dids", 0, headers);

            double res_old = scorer.EstimateScore_old(v);
            double res_new = scorer.EstimateScore(v);
            EXPECT( abs(res_old - res_new) < 1.0/pow(10, 10) );
            cout << "   ok" << endl;
        }

        // Standard dev
        SETUP( "Standard deviation" ) {
            cout << "Standard deviation" << endl;
            // Create the scorer
            Scorer scorer("sd", 0, headers);

            double res_old = scorer.EstimateScore_old(v);
            double res_new = scorer.EstimateScore(v);
            EXPECT( abs(res_old - res_new) < 1.0/pow(10, 10) );
            cout << "   ok" << endl;
        }

        // Mean stddev
        SETUP( "Mean stddev" ) {
            cout << "Mean stddev" << endl;
            // Create the scorer
            Scorer scorer("rsd1", 0, headers);

            double res_old = scorer.EstimateScore_old(v);
            double res_new = scorer.EstimateScore(v);
            EXPECT( abs(res_old - res_new) < 1.0/pow(10, 10) );
            cout << "   ok" << endl;
        }

        // Min stddev
        SETUP( "Min stddev" ) {
            cout << "Min stddev" << endl;
            // Create the scorer
            Scorer scorer("rsd2", 0, headers);

            double res_old = scorer.EstimateScore_old(v);
            double res_new = scorer.EstimateScore(v);
            EXPECT( abs(res_old - res_new) < 1.0/pow(10, 10) );
            cout << "   ok" << endl;
        }

        // Median stddev
        SETUP( "Median stddev" ) {
            cout << "Median stddev" << endl;
            // Create the scorer
            Scorer scorer("rsd3", 0, headers);

            double res_old = scorer.EstimateScore_old(v);
            double res_new = scorer.EstimateScore(v);
            EXPECT( abs(res_old - res_new) < 1.0/pow(10, 10) );
            cout << "   ok" << endl;
        }

        // Entropie
        SETUP( "Entropie" ) {
            cout << "Entropie" << endl;
            // Create the scorer
            Scorer scorer("entropy", 0, headers);

            double res_old = scorer.EstimateScore_old(v);
            double res_new = scorer.EstimateScore(v);
            // Error propagation is bad due to log computations
            EXPECT( abs(res_old - res_new) < 1.0/pow(10, 2) );
            cout << "   ok" << endl;
        }
    }
};

extern lest::tests & specification();

MODULE( specification(), module )

#include <iostream>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <cmath>

#include "lest.hpp"
#include "vect_opera.hpp"

using namespace std;



const lest::test specification[] =
{
    CASE( "test CalcPearsonCorr (new vs old version)" )
    {
        cout << "Pearson correlation verification" << endl;
        srand(time(NULL));
        vector<float> x, y;
        for (uint test_idx=0 ; test_idx<100 ; test_idx++) {
            for (uint i=0 ; i<2000 ; i++) {
                float xi = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX/9997.0));
                x.push_back(xi);
                float yi = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX/9997.0));
                y.push_back(yi);
            }

            double res_old = CalcPearsonCorr_old(x, y);
            double res_new = CalcPearsonCorr(x, y);

            // cout << res_old << " " << res_new << " " << (res_new - res_old) << endl;
            EXPECT( res_old - res_new < 1.0/pow(10, 10) );
        }
        cout << "   ok" << endl;
    }
};


int main( int argc, char * argv[] )
{
    return lest::run( specification, argc, argv );
}
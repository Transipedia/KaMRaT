#include <iostream>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <cmath>

#include "lest.hpp"
#include "vect_opera.hpp"

using namespace std;



const lest::test module[] =
{
    CASE( "test CalcPearsonCorr (new vs old version)" )
    {
        cout << "Pearson correlation verification" << endl;
        srand(time(NULL));
        vector<float> x, y;
        for (uint test_idx=0 ; test_idx<20 ; test_idx++) {
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
    },


    // --- Spearman requirements ---

    CASE( "Test rank generation for pearson test" ) {
        cout << "float order test" << endl;

        vector<float> lst_test {1, 1, 3, 0, 2, 1, 2, 1};
        vector<uint> real_order{3, 0, 1, 5, 7, 4, 6, 2};
        vector<uint> order;

        getOrder(lst_test, order);

        for (uint i=0 ; i<lst_test.size() ; i++) {
            EXPECT(order[i] == real_order[i]);
        }
        cout << "   ok" << endl;

        cout << "rank from order test" << endl;

        vector<float> ranks;
        vector<float> real_ranks{2.5, 2.5, 7, 0, 5.5, 2.5, 5.5, 2.5};
        orderToRank(lst_test, order, ranks);

        for (uint i=0 ; i<lst_test.size() ; i++) {
            EXPECT(ranks[i] == real_ranks[i]);
        }

        cout << "   ok" << endl;
    },


    CASE( "test CalcSpearmanCorr (new vs old version)" )
    {
        cout << "Spearman correlation verification" << endl;
        srand(time(NULL));
        vector<float> x, y;
        for (uint test_idx=0 ; test_idx<20 ; test_idx++) {
            for (uint i=0 ; i<200 ; i++) {
                float xi = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX/9997.0));
                x.push_back(xi);
                float yi = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX/9997.0));
                y.push_back(yi);
            }

            double res_old = CalcSpearmanCorr_old(x, y);
            double res_new = CalcSpearmanCorr(x, y);
            
            EXPECT( abs(res_old - res_new) < 1.0/pow(10, 10) );
        }
        cout << "   ok" << endl;
    }
};


extern lest::tests & specification();

MODULE( specification(), module )
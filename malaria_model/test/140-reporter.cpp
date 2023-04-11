#include <iostream>
#include <map>

#include "third_party/catch2/catch.hpp"

#include "common/common.h"
#include "util/util.h"

#include "village/village.h"
#include "village/mosquito_dispersion.h"
#include "human/human.h"
#include "human/blood_system.h"
#include "simulation/reporter.h"

TEST_CASE("140.1: Reporter, drug failure rate", "[rp][rp_drug_failure]") {
    
    std::cout << "Reporter: testing daily_report_drug_fail_rate_calc\n";

    std::vector<int> intake_list {
        1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 11, 12,
        0, 0, 0
    };
    std::vector<int> fail_list {
        0, 0, 0,
        0, 0, 1, 2, 1, 2, 4, 3, 7, 0, 10, 1, 
    };
    std::vector<float> exp_rate_list_len1 {
        0, 0, 0,
        0, 0, 1.0/3.0, 0.5, 0.2, 2.0/6.0, 4.0/7.0,
        3.0/8.0, 7.0/9.0, 0, 10.0/11.0, 1.0/12
    };

    std::vector<float> exp_rate_list_len5 {
        0, 0, 0,
        0, 0, 0, 0,
        (0.0 + 0 + 1 + 2 + 1)/(1.0 + 2 + 3 + 4 + 5),
        (0.0 + 1 + 2 + 1 + 2)/(2.0 + 3 + 4 + 5 + 6),
        (1.0 + 2 + 1 + 2 + 4)/(3.0 + 4 + 5 + 6 + 7),
        (2.0 + 1 + 2 + 4 + 3)/(4.0 + 5 + 6 + 7 + 8),
        (1.0 + 2 + 4 + 3 + 7)/(5.0 + 6 + 7 + 8 + 9),
        (2.0 + 4 + 3 + 7 + 0)/(6.0 + 7 + 8 + 9 + 0),
        (4.0 + 3 + 7 + 0 + 10)/(7.0 + 8 + 9 + 0 + 11),
        (3.0 + 7 + 0 + 10 + 1)/(8.0 + 9 + 0 + 11 + 12)
    };

    int drug_failure_if_clinical_on_day = 4;
    // size_t measurement_length = 1;

    std::cout << "\tmeasurement length of 1";
    for (size_t ii=0; ii < fail_list.size(); ii++) {
        // std::cout << "ii=" << ii << "\n";
        float rate = simulation::ReportManager::daily_report_drug_fail_rate_calc(
            intake_list[ii],
            fail_list[ii],
            drug_failure_if_clinical_on_day,
            1,
            (ii == 0)
        );
        REQUIRE(rate == exp_rate_list_len1[ii]);
    }
    std::cout << " ... pass\n";

    std::cout << "\tmeasurement length of 5";
    for (size_t ii=0; ii < fail_list.size(); ii++) {
        // std::cout << "ii=" << ii << "\n";
        float rate = simulation::ReportManager::daily_report_drug_fail_rate_calc(
            intake_list[ii],
            fail_list[ii],
            drug_failure_if_clinical_on_day,
            5,
            (ii == 0)
        );
        REQUIRE(rate == exp_rate_list_len5[ii]);
    }
    std::cout << " ... pass\n";

}
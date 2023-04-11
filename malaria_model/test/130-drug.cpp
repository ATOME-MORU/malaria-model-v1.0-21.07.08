#include <iostream>

#include "third_party/catch2/catch.hpp"

#include "common/common.h"


TEST_CASE("130.4: Drug", "[dg][dg_print]") {
    common::DrugName dg_last = common::DrugName::Last;
    common::DrugName dg_first = common::DrugName::First;
    std::vector<common::DrugName> dg_good_list{
        dg_first,
        static_cast<common::DrugName>(static_cast<int>(dg_first)-1),
        dg_last
    };
    std::cout << "Printing defined values:\n";
    for (const auto& dg : dg_good_list) {
        REQUIRE_NOTHROW(
            std::cout << "\t" << static_cast<int>(dg)
            << ", " << dg << "\n");
    }

    std::vector<common::DrugName> dg_bad_list{
        static_cast<common::DrugName>(static_cast<int>(dg_first)-2),
        static_cast<common::DrugName>(static_cast<int>(dg_last)+1)
    };
    std::cout << "Printing undefined values:\n";
    for (const auto& dg : dg_bad_list) {
        std::cout << "\t";
        REQUIRE_THROWS(std::cout << dg);
        std::cout << "\n";
    }
}

TEST_CASE("130.1: Drug", "[dg][dg_init]") {
    // common::DrugManager dm("./data/drug/tact_efficacy_table.csv");
    common::DrugManager dm(
        "./data/test/drug_efficacy_table.csv",
        "./data/test/drug_loss_array.csv",
        true
        );
    
    // common::DrugManager dm("./data/regions_datajuly2.txt");

}

TEST_CASE("130.2: drug_effect", "[dg][dg_effect]") {

    common::DrugManager dm(
        "./data/test/drug_efficacy_table.csv",
        "./data/test/drug_loss_array.csv"
        );

    std::cout << "DM: testing drug_effect function\n";
    int cleared = 0;

    common::genotype_t gt = 1;
    common::DrugName dg = common::DrugName::kDP;
    common::StageName stg = common::StageName::kInfectious;

    const int kNumTests = 100000;
    for (int tt = 0; tt < kNumTests; tt++) {
        if (dm.drug_effect(gt, dg, common::StageName::kLiver)) {
            cleared++;
        }
    }
    std::cout << "\tNo effect on liver stage parasites";
    REQUIRE(cleared == 0);
    std::cout << " ... pass\n";
    for (int tt = 0; tt < kNumTests; tt++) {
        if (dm.drug_effect(gt, dg, stg)) {
            cleared++;
        }
    }

    float expected_efficacy = dm.get_drug_effect_rate(gt, dg);
    float actual_efficacy = static_cast<float>(cleared)/kNumTests;
    std::cout << "\tGenotype: " << static_cast<int>(gt) << "\n";
    std::cout << "\tDrug: " << dg << "\n";
    std::cout << "\tExpected efficacy: " << expected_efficacy << "\n";
    std::cout << "\tActural #cleared: " << cleared << " of "
                << kNumTests <<" ("<< 100.0*actual_efficacy << "%)";

    REQUIRE( (actual_efficacy - expected_efficacy)/expected_efficacy < 0.01 );
    std::cout << " ... pass\n";

}

TEST_CASE("130.3: drug_loss", "[dg][dg_loss]") {
    common::DrugManager dm(
        "./data/test/drug_efficacy_table.csv",
        "./data/test/drug_loss_array.csv"
    );

    std::cout << "DM: testing drug_loss function\n";

    //////////////////////////////
    // Part 1 one to one

    std::vector<common::DrugName> vDgFrom {
        common::DrugName::kNoDrug,
        common::DrugName::kAS,
        common::DrugName::kLM,
        common::DrugName::kAQ,
        common::DrugName::kPPQ,
        common::DrugName::kMQ,
        common::DrugName::kSP,
        common::DrugName::kCQ,

        common::DrugName::kAL,
        common::DrugName::kASAQ,
        common::DrugName::kDP,
        common::DrugName::kASMQ,

        common::DrugName::kALb,
        common::DrugName::kASAQb,
        common::DrugName::kDPb,
        common::DrugName::kASMQb,

        common::DrugName::kASMQPPQ,
        common::DrugName::kALAQ

    };
    std::vector<common::DrugName> vDgTo {
        common::DrugName::kNoDrug,
        common::DrugName::kNoDrug,
        common::DrugName::kNoDrug,
        common::DrugName::kNoDrug,
        common::DrugName::kNoDrug,
        common::DrugName::kNoDrug,
        common::DrugName::kNoDrug,
        common::DrugName::kNoDrug,
        
        common::DrugName::kALb,
        common::DrugName::kASAQb,
        common::DrugName::kDPb,
        common::DrugName::kASMQb,

        common::DrugName::kNoDrug,
        common::DrugName::kNoDrug,
        common::DrugName::kNoDrug,
        common::DrugName::kNoDrug,

        common::DrugName::kASMQbPPQ,
        common::DrugName::kALbAQ

    };
    for (size_t idd = 0; idd < vDgFrom.size(); idd++) {
        // const common::DrugName kDgFrom = common::DrugName::kDP;
        // const common::DrugName kDgTo = common::DrugName::kDPb;
        const common::DrugName kDgFrom = vDgFrom[idd];
        const common::DrugName kDgTo = vDgTo[idd];
        const int kNumTests = 1000000;
        int num_drug_loss = 0;
        for (int tt = 0; tt < kNumTests; tt++) {
            common::DrugName dd = kDgFrom;
            if (dm.drug_loss(dd)) {
                REQUIRE(dd == kDgTo);
                num_drug_loss++;
            } else {
                REQUIRE(dd == kDgFrom);
            }
        }
        if (kDgFrom == common::DrugName::kNoDrug) {
            REQUIRE(num_drug_loss == 0);
            std::cout << "> no loss registered for " << kDgFrom <<" ... pass\n";
        } else {
            float actual_rate = static_cast<float>(num_drug_loss) / kNumTests;
            float expected_rate = dm.get_drug_loss_rate(kDgFrom);
            std::cout << "> simple one to one, from " << kDgFrom << " to " << kDgTo << ":\n";
            std::cout << "\t exp:\t"  << expected_rate << "\n"
                        << "\t act:\t" << actual_rate;

            REQUIRE((actual_rate - expected_rate)/expected_rate < 0.01);
            std::cout << " ... pass\n";
        }
    }

    //////////////////////////////
    // Part 2 kASMQbPPQ, kALbAQ

    const std::vector< std::vector<common::DrugName> > kListOutcomes{
        {   common::DrugName::kASMQbPPQ,
            common::DrugName::kASMQb,
            common::DrugName::kPPQ,
            common::DrugName::kNoDrug,
        },
        {   common::DrugName::kALbAQ,
            common::DrugName::kALb,
            common::DrugName::kAQ,
            common::DrugName::kNoDrug,
        }
    };

    for (auto& outcomes : kListOutcomes) {

        const common::DrugName kDgFrom = outcomes[0];
        const common::DrugName kDgOutcomeA = outcomes[1];
        const common::DrugName kDgOutcomeB = outcomes[2];
        const common::DrugName kDgOutcome0 = outcomes[3];

        std::cout << "> one to many, from " << kDgFrom << ":\n";
        const std::vector<float> exp_outcome{
                        // ab
                        (1.0f-dm.get_drug_loss_rate(kDgOutcomeA))
                        *(1.0f-dm.get_drug_loss_rate(kDgOutcomeB)),
                        // a
                        (1.0f-dm.get_drug_loss_rate(kDgOutcomeA))
                        *dm.get_drug_loss_rate(kDgOutcomeB),
                        // b
                        dm.get_drug_loss_rate(kDgOutcomeA)
                        *(1.0f-dm.get_drug_loss_rate(kDgOutcomeB)),
                        // 0
                        dm.get_drug_loss_rate(kDgOutcomeA)
                        *dm.get_drug_loss_rate(kDgOutcomeB)
        };
        float exp_sum = 0.0;
        for (auto& ee : exp_outcome) {
            exp_sum += ee;
        }
        std::cout << "\t sum of all possible outcomes' probabilities is 1 (actual:"
                    << exp_sum << ")\n";
        REQUIRE(exp_sum - 1.0 < 0.00000001);

        const int kNumTests = 5000000;
        std::vector<int> act_outcome_cout(4,0);

        for (int tt = 0; tt < kNumTests; tt++) {
            common::DrugName dd = kDgFrom;
            if (!dm.drug_loss(dd)) {
                REQUIRE(dd == kDgFrom);
                act_outcome_cout[0]++;
            } else {
                if (dd == kDgOutcomeA) {
                    act_outcome_cout[1]++;
                } else if (dd == kDgOutcomeB) {
                    act_outcome_cout[2]++;
                } else if (dd == kDgOutcome0) {
                    act_outcome_cout[3]++;
                } else {
                    std::cout << "fail: invalid loss outcome ("
                                << dd << ") for " << kDgFrom << "\n";
                    REQUIRE(false);
                }
            }
        }
        for (size_t ioo = 0; ioo < exp_outcome.size(); ioo++) {
            float act_outcome = static_cast<float>(act_outcome_cout[ioo])/kNumTests;
            std::cout << "\t outcome: " << outcomes[ioo] << "\n"
                        << "\t exp: " << exp_outcome[ioo]
                        << "\n\t act: " << act_outcome;
            REQUIRE((act_outcome - exp_outcome[ioo])/exp_outcome[ioo] < 0.01);
            std::cout << " ... pass\n";
        }

    }
}
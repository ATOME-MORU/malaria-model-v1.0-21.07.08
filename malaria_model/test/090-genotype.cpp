#include <iostream>
#include <bitset>
// #include <iomanip>

#include <string>
#include <cmath> //pow

#include <algorithm> // std::fill_n
#include <vector>

#include "third_party/catch2/catch.hpp"
#include "util/randomness.h"
#include "util/util.h"
#include "util/statistics.h"

#include "common/genotype.h"


TEST_CASE("90.8", "[gt][gt_sum_vector_on_sites]") {
    std::cout << "\nTEST_CASE 90.8\nTesting sum_vector_on_sites\n";
    std::vector<float> v(common::kParasiteType_size, 0);

    std::cout << "common::kParasiteType_size: " << common::kParasiteType_size << "\n";

    std::vector<common::ParasiteType> pt_list{
        common::ParasiteType::kKNY__Y1, // one site mutated, 580Y
        common::ParasiteType::kKNY__Y2 // two sites mutated
        // common::ParasiteType::kKNYNYY2 // three sites mutated
    };
    std::vector<float> pt_prob{
        0.5,    // half have mutation
        // 0.75   // 3/4 have one of the two mutations, for sum_vector_on_sites_rm1
        0.25   // 1/4 have both of the two mutations, for sum_vector_on_sites
        // 7.0/8.0 // 7/8 have one of the three mutations
    };
    std::vector<size_t> mod_4_min{
        1,
        // 0 // sum_vector_on_sites_rm1
        2 // sum_vector_on_sites
    };

    for (size_t pp = 0; pp < pt_list.size(); pp++) {

        std::cout << "ParasiteType: " << pt_list[pp] << "\n";

        const int kNum_tests = 100000;
        float sum_total = 0.0;
        for (int tt = 0; tt < kNum_tests; tt++) {

            float gw_sum = 0.0;
            for (float& vv : v) {
                vv = util::get_rand_uniform();
                gw_sum += vv;
            }
            float sum_func = common::GenotypeManager::sum_vector_on_sites(v, pt_list[pp]);
            float sum_target = 0.0;
            for (size_t ii = 0; ii < common::kParasiteType_size; ii++) {
                if (ii % 4 > mod_4_min[pp]) {
                    sum_target += v[ii];
                }
            }

            if (tt < 5) {
                std::cout << "\tsum_target: " << sum_target;
                std::cout << ",\tsum_func: " << sum_func;
            }
            
            Approx sum_target_appx = Approx(sum_target).epsilon(0.01);
            REQUIRE(sum_func == sum_target_appx);
                
            // sum_vector_on_sites should sum all genotype weights when given wild type
            REQUIRE(gw_sum == common::GenotypeManager::sum_vector_on_sites(v, common::ParasiteType::kWild));

            if (tt < 5) {
                std::cout << " ... pass\n";
            }
            sum_total += sum_func;
        }

        const float kSum_total_target = (common::kParasiteType_size-1.0) // exclude kNoParasite at 64
                                            *0.5 // uniform 0.0 ~ 1.0, avg. 0.5
                                            *pt_prob[pp] // prob to have one of the mutations
                                            *kNum_tests;
        std::cout << "\tAvg. over "<< kNum_tests << " tests:\n\t"
                    <<"\tsum_total: " << sum_total
                    << ",\n\t\tsum_total_target: " << kSum_total_target;
        Approx sum_total_target_appx = Approx(kSum_total_target).epsilon(0.01);
        REQUIRE(sum_total == sum_total_target_appx);
        std::cout << " ... pass\n";

    }


}


TEST_CASE("90.7", "[gt][gt_map_str_to_pt]") {
    std::cout << "\nTEST_CASE 90.7\nTesting map_str_to_pt function\n";
    std::cout << " - with good values";
    REQUIRE(common::GenotypeManager::map_str_to_pt("KNY--C1") == common::ParasiteType::kKNY__C1);
    REQUIRE(common::GenotypeManager::map_str_to_pt("wild") == common::ParasiteType::kKNY__C1);
    REQUIRE(common::GenotypeManager::map_str_to_pt("TYFYFY2") == common::ParasiteType::kTYFYFY2);
    std::cout << " ... pass\n";
    std::cout << " - with bad values\n";
    REQUIRE_THROWS(common::GenotypeManager::map_str_to_pt("KNYC1"));
    std::cout << " ... pass\n";
}


TEST_CASE( "90.6", "[gt][gt_fun_mutation_combined]" ) {
    
    std::cout << "\nTEST_CASE 90.6\nTesting fun_mutation_combined function\n";

    const std::vector<float> kMutation_probability{
        0.0001,
        0.0002,
        0.0003,
        0.0004,
        0.0005,
        0.0006
    };
    const int kNum_bits = common::GenotypeManager::kNum_bits_in_genotype;
    const int kNum_tests = 10000000;

    std::vector<bool> back_mutation_options{true, false};
    for (bool bm_option: back_mutation_options) {

        std::cout << "Back mutation : " << ((bm_option)? "on":"off") << "\n";

        common::GenotypeManager gm(
            kMutation_probability,
            static_cast<int>(kMutation_probability.size()),
            bm_option
        );

        // genotype_t gg = 0;
        common::genotype_t gg = static_cast<common::genotype_t>(common::ParasiteType::kKNY__C1);

        int num_print = 10;
        for (size_t tt = 0; tt < kNum_tests; tt++) {
            common::genotype_t gg_orig = gg;
            int rr = gm.fun_mutation_combined(gg);
            if (rr == 0) {
                REQUIRE(gg_orig == gg);
            } else {
                if (num_print-- > 0 || rr > 1){
                    std::cout << "\ttt = " << tt << " from "
                                << std::bitset<kNum_bits>(gg_orig)
                                << " (" << std::bitset<kNum_bits>(gg_orig).count() << " ones)"
                                << " to "
                                << std::bitset<kNum_bits>(gg)
                                << " (" << std::bitset<kNum_bits>(gg).count() << " ones) returns "
                                << rr << "\n";
                }
                REQUIRE(std::bitset<kNum_bits>(gg^gg_orig).count() == rr);
            }

        }
    }



}


TEST_CASE( "90.5", "[gt][gt_pt_conversion]" ) {
    // common::genotype_t gt = 60;
    std::vector<common::genotype_t> gt_list_good{0,32,33,63,64}; // error from 256
    std::vector<common::genotype_t> gt_list_bad{65,100,255}; // compile errors from 256, or less than 0
    std::cout << "Converting defined values:\n";
    for (const auto& gt : gt_list_good) {
        std::cout << "\tConvert from gt = " << static_cast<int>(gt) << ": ";
        common::ParasiteType pt = static_cast<common::ParasiteType>(gt);
        std::cout << "\t(int)pt: " << static_cast<int>(pt);
        REQUIRE_NOTHROW(std::cout << "\tpt: " << pt << "\n");
    }
    std::cout << "Converting undefined values:\n";
    for (const auto& gt : gt_list_bad) {
        std::cout << "\tConvert from gt = " << static_cast<int>(gt) << ": ";
        common::ParasiteType pt = static_cast<common::ParasiteType>(gt);
        std::cout << "\t(int)pt: " << static_cast<int>(pt);
        REQUIRE_THROWS(std::cout << "\tpt: " << pt << "\n");
        std::cout << "\n";
        // std::cout << "\tpt: " << pt << "\n";
    }    
}

TEST_CASE( "90.1: Genotype: bitwise operations for mutation definition", "[gt][gt_basics]" ) {

    const int kNum_tests = 5000;

    const int kNum_bits_in_byte = 8;
    const int kNum_positions = sizeof(common::genotype_t) * kNum_bits_in_byte;

    // const float kMutation_probability = 0.1;
    const int kNum_probabilities = 4; // each probability is applied to all positions when tested
    const float kMutation_probability[kNum_probabilities] = {0.0, 0.2, 0.7, 1.0};

    int num_mutations = 0;

    common::genotype_t gt = 0;
    common::genotype_t gt_temp = 0;
    common::genotype_t mask = 0;

    float rnd_nmb = 0.0;

    int num_printed_mutations = 3;

    for (int pp = 0; pp < kNum_probabilities; pp++){
        num_mutations = 0;

        std::cout << "Test 90.1/" << pp
                <<":\nFor a genotype with " << kNum_positions
                << " positions, each with " << kMutation_probability[pp] * 100 << "%"
                << " chance to mutate:\n";
        
        for (int tt = 0; tt < kNum_tests; tt++) {
            rnd_nmb = 0.0;
            mask = 0;

            if (tt < num_printed_mutations) {
                std::cout << "\nTime step " << tt;
                std::cout << " Genotype (before): " << std::bitset<kNum_positions>(gt) << ". ";
                std::cout << "Mask construction:\n";
            }

            for (uint bb = 0; bb < kNum_positions; bb++) {
                util::set_random_numbers_uniform(&rnd_nmb, 1);
                if (tt < num_printed_mutations){
                    std::cout << "Mutation Mask: from " << std::bitset<kNum_positions>(mask);
                }
                mask <<= 1;
                mask += (rnd_nmb < kMutation_probability[pp]) ? 1 : 0;
                if (tt < num_printed_mutations){
                    std::cout << " to " << std::bitset<kNum_positions>(mask)
                                << ", rnd_nmb: " << rnd_nmb << "(" << kMutation_probability[pp] << ")\n";
                }
            }
            gt_temp = gt;
            gt ^= mask; //allows back mutation
            // gt |= mask;

            gt_temp ^= gt; // the number of ones in gt_temp is the number of bits that were flipped by mask
            num_mutations += std::bitset< kNum_positions >(gt_temp).count();

            if (tt < num_printed_mutations){
                std::cout << "Genotype (after): " << std::bitset<kNum_positions>(gt)
                        << " with " << std::bitset< kNum_positions >(gt_temp).count() << " mutations. \n";
            }

        }


        int target_num_mutations = kNum_tests*kNum_positions*kMutation_probability[pp];
        float tolerance = 0.02; // 1%

        std::cout << "\nAfter " << kNum_tests << " time steps, "
                  << num_mutations << " mutations happend, expecting "
                  << target_num_mutations << " +/-" << tolerance * 100 << "% ... ";
        REQUIRE(static_cast<float>(std::abs(num_mutations - target_num_mutations)) <= target_num_mutations * tolerance);
        std::cout << "PASS\n\n";
        
    }

    std::cout << "\nExample: \nGenotype (before mutation): " << std::bitset<kNum_positions>(gt)
                << " (" << std::bitset<kNum_positions>(gt).count() << " ones)" << "\n";

    std::bitset<kNum_positions> mask_b(mask);
    std::cout << "Mutation Mask: " << mask_b
                << " (" << mask_b.count() << " ones)" << "\n";
    gt ^= mask;
    std::cout << "Genotype (after mutation): "
                << std::bitset<kNum_positions>(gt)
                << " (" << std::bitset<kNum_positions>(gt).count() << " ones)\n";

}

TEST_CASE( "90.1.1: Genotype: Initilisation", "[gt][gt_init]" ) {

    const int kMax_num_sites_in_combined_mutation = 2;
    const std::vector<float> kMutation_probability{0.1, 0.2, 0.3, 0.4};
    
    common::GenotypeManager gm(kMutation_probability, kMax_num_sites_in_combined_mutation);


    gm.print_all_combined_probabilities();

    std::vector<float> all_probabilities = gm.get_all_combined_probabilities();


    REQUIRE(
        all_probabilities.size() ==
            util::from_n_choose_up_to_k_get_num_combinations(
                kMutation_probability.size(), kMax_num_sites_in_combined_mutation
            )
    );

    float non_mutation_probability = 1.0;
    for (uint ss = 0; ss < kMutation_probability.size(); ss++) {
        non_mutation_probability *= (1-kMutation_probability.at(ss));
    }

    REQUIRE( all_probabilities.at(0) == non_mutation_probability );
}

TEST_CASE( "90.2: Genotype: step_mutation_independent_kernel function", "[gt][gt_step_mutation_independent_kernel]" ) {

    const int kNum_positions = 8;
    const float kMutation_probability_max = 1.0;
    float mutation_probability_array[kNum_positions] = {};
    int position_mask_array[kNum_positions] = {};
    for (int ii = 0; ii < kNum_positions; ii++) {
        mutation_probability_array[ii] = kMutation_probability_max 
                                            / static_cast<float>(kNum_positions)
                                            * static_cast<float>(ii);
        position_mask_array[ii] = 1<<(kNum_positions-ii-1);
    }
    mutation_probability_array[kNum_positions] = kMutation_probability_max;


    const int kNum_infections = 100;
    common::genotype_t gt_array[kNum_infections] = {};
    std::fill(gt_array, gt_array + kNum_infections, 0);


    const int kNum_time_steps = 3650;

    int num_mutations[kNum_positions];
    std::fill(num_mutations, num_mutations + kNum_positions, 0);
    common::genotype_t gt_array_temp[kNum_infections];

    for (int tt = 0; tt < kNum_time_steps; tt++) {

        std::copy(gt_array, gt_array + kNum_infections, gt_array_temp);

        common::GenotypeManager::step_mutation_independent_kernel(
            gt_array,
            kNum_infections,

            mutation_probability_array,
            kNum_positions
        );

        for (int gg = 0; gg < kNum_infections; gg++) {
            gt_array_temp[gg] ^= gt_array[gg];
            for (int ii = 0; ii < kNum_positions; ii++) {
                if ((gt_array_temp[gg] & position_mask_array[ii]) == position_mask_array[ii]) {
                    num_mutations[ii]++;
                }
            }
        }

    }

    float tolerance = 0.01; // 1%

    std::cout << "After " << kNum_time_steps << " steps of " << kNum_infections << " infections:\n";
    for (int ii = 0; ii < kNum_positions; ii++) {
        int num_expected_mutations = mutation_probability_array[ii] * kNum_time_steps * kNum_infections;
        std::cout << "Position " << ii 
                    << " with mutation probability " << mutation_probability_array[ii]
                    << ", number of actual mutations:" << num_mutations[ii]
                    << ", expecting " << num_expected_mutations
                    << " +/- " << tolerance*100 << "% ...";
        REQUIRE(std::abs(num_expected_mutations - num_mutations[ii]) <= num_expected_mutations * tolerance);
        std::cout << " PASS\n";
    }

}


TEST_CASE("90.3: Genotype: killing function performance test", "[gt][gt_kill_perf]") {
    const int kNum_tests = 10;

    const int kNum_genotypes = 128;
    const int kNum_majority_types = 1;
    const int kNum_minority_types = kNum_genotypes - kNum_majority_types;
    // const int kMajority_type_count = 2000000*10;
    const int64_t kMajority_type_count = 10000000000;
    const int64_t kMinority_type_count = 200;
    const int64_t kTotal_count = kMajority_type_count*kNum_majority_types
                            + kMinority_type_count*kNum_minority_types;
    // const int64_t kKill_count = kTotal_count * 0.9;
    const int64_t kKill_count_left = 1000000;
    const int64_t kKill_count = kTotal_count - kKill_count_left;

    int64_t parasite_count[kNum_genotypes] = {};
    int64_t parasite_count_survives[kNum_genotypes] = {};

    std::fill(parasite_count, parasite_count+kNum_majority_types, kMajority_type_count);
    std::fill(parasite_count+kNum_majority_types, parasite_count+kNum_genotypes, kMinority_type_count);

    std::fill(parasite_count_survives, parasite_count_survives+kNum_genotypes, 0);



    double wall_time_begin = 0.0;
    double wall_time_end = 0.0;


    wall_time_begin = util::get_wall_time();
    
    // float uniform_random_array[kKill_count] = {};
    float rnd_nmb= 0.0;
    for (int tt = 0; tt < kNum_tests; tt++ ){

        std::fill(parasite_count, parasite_count+kNum_majority_types, kMajority_type_count);
        std::fill(parasite_count+kNum_majority_types, parasite_count+kNum_genotypes, kMinority_type_count);

        std::fill(parasite_count_survives, parasite_count_survives+kNum_genotypes, 0);



        // util::set_random_numbers_uniform(uniform_random_array, kKill_count);

        int64_t total_cout = kTotal_count;
        // for (int64_t kk = 0; kk < kKill_count; kk++) {
        for (int64_t kk = 0; kk < kKill_count_left; kk++) {

            // int kill_at = uniform_random_array[kk] * total_cout;
            util::set_random_numbers_uniform(&rnd_nmb, 1);
            // int kill_at = rnd_nmb * total_cout;

            int survive_at = rnd_nmb * kTotal_count;

            int accumulator = 0;

            for (int tt = 0; tt < kNum_genotypes; tt++) {
                // if (parasite_count[tt] + accumulator >= kill_at) {
                if (parasite_count[tt] + accumulator >= survive_at) {
                    // parasite_count[tt]--;
                    // total_cout--;

                    parasite_count_survives[tt]++;
                    break;
                } else {
                    accumulator += parasite_count[tt];
                }
            }
        }



        int64_t remain_count = 0;
        
        for (int tt = 0; tt < kNum_genotypes; tt++) {
            // remain_count += parasite_count[tt];
            remain_count += parasite_count_survives[tt];
        }

        // REQUIRE(remain_count == (kTotal_count - kKill_count));
        REQUIRE(remain_count == kKill_count_left);

    }
    wall_time_end = util::get_wall_time();

    std::cout << "Per-kill sample with cascaded loop locator: " << std::setprecision(5) << wall_time_end - wall_time_begin << " seconds.\n";




    // wall_time_begin = util::get_wall_time();
    
    // // float uniform_random_array[kKill_count] = {};
    // // std::vector<int> all_parasites(kTotal_count, 0);
    // // int kill_array[kKill_count/200] = {};
    // int* kill_array = new int[kKill_count];
    // std::vector<int> parasite_histogram(kNum_genotypes, 0);

    // for (int tt = 0; tt < kNum_tests; tt++ ){

    //     std::fill(parasite_count, parasite_count+kNum_majority_types, kMajority_type_count);
    //     std::fill(parasite_count+kNum_majority_types, parasite_count+kNum_genotypes, kMinority_type_count);


    //     for (int tt = 0; tt < kNum_genotypes; tt++) {
    //         parasite_histogram[tt] = parasite_count[tt];
    //     }

    //     util::set_random_numbers_discrete<int>(kill_array, kKill_count, parasite_histogram);

    //     for (int kk = 0; kk < kKill_count; kk++) {
    //         parasite_count[kill_array[kk]]--;
    //     }

    //     int remain_count = 0;
    //     for (int tt = 0; tt < kNum_genotypes; tt++) {
    //         remain_count += parasite_count[tt];
    //     }
    //     REQUIRE(remain_count == (kTotal_count - kKill_count));

    // }
    // delete[] kill_array;
    // wall_time_end = util::get_wall_time();

    // std::cout << "Per-kill sample from discrete distribution: " << std::setprecision(5) << wall_time_end - wall_time_begin << " seconds.\n";







    // wall_time_begin = util::get_wall_time();
    
    // std::vector<int> all_parasites(kTotal_count, 0);
    // for (int tt = 0; tt < kNum_tests; tt++ ){

    //     std::fill(parasite_count, parasite_count+kNum_majority_types, kMajority_type_count);
    //     std::fill(parasite_count+kNum_majority_types, parasite_count+kNum_genotypes, kMinority_type_count);

    //     int accumulator = 0;
    //     for (int tt = 0; tt < kNum_genotypes; tt++) {
    //         std::fill(  all_parasites.begin() + accumulator, 
    //                     all_parasites.begin() + accumulator + parasite_count[tt],
    //                     tt
    //         );
    //         accumulator += parasite_count[tt];
    //     }
    //     std::random_shuffle(all_parasites.begin(), all_parasites.end());

    //     std::fill(parasite_count, parasite_count + kNum_genotypes, 0);
    //     for (int pp = kKill_count; pp < kTotal_count; pp++) {
    //         parasite_count[all_parasites.at(pp)]++;
    //     }

    //     int remain_count = 0;
    //     for (int tt = 0; tt < kNum_genotypes; tt++) {
    //         remain_count += parasite_count[tt];
    //     }
    //     REQUIRE(remain_count == (kTotal_count - kKill_count));

    // }
    // wall_time_end = util::get_wall_time();

    // std::cout << "Lineup with shuffle and truncate: " << std::setprecision(5) << wall_time_end - wall_time_begin << " seconds.\n";





    // wall_time_begin = util::get_wall_time();

    // std::vector<bool> all_parasites_kill_flag(kTotal_count, false);
    // std::fill(all_parasites_kill_flag.begin(), all_parasites_kill_flag.end(), false);

    
    // for (int tt = 0; tt < kNum_tests; tt++ ){

    //     std::fill(parasite_count, parasite_count+kNum_majority_types, kMajority_type_count);
    //     std::fill(parasite_count+kNum_majority_types, parasite_count+kNum_genotypes, kMinority_type_count);


    //     // util::set_random_numbers_uniform(uniform_random_array, kKill_count);

    //     int accumulator = 0;
    //     for (int tt = 0; tt < kNum_genotypes; tt++) {
    //         std::fill(  all_parasites.begin() + accumulator, 
    //                     all_parasites.begin() + accumulator + parasite_count[tt],
    //                     tt
    //         );
    //         accumulator += parasite_count[tt];
    //     }
    //     std::fill(all_parasites_kill_flag.begin(), all_parasites_kill_flag.end(), false);

    //     for (int kk = 0; kk < kKill_count; kk++) {

    //         // int kill_at = uniform_random_array[kk] * total_cout;
    //         float rnd_nmb = 0.0;
    //         util::set_random_numbers_uniform(&rnd_nmb, 1);
    //         int kill_at = rnd_nmb * kTotal_count;

    //         while(all_parasites_kill_flag[kill_at]) {
    //             kill_at++;
    //             if (kill_at >= kTotal_count) {
    //                 kill_at = 0;
    //             }
    //         }

    //         all_parasites_kill_flag[kill_at] = true;
    //         parasite_count[all_parasites.at(kill_at)]--;

    //     }

    //     int remain_count = 0;
    //     for (int tt = 0; tt < kNum_genotypes; tt++) {
    //         remain_count += parasite_count[tt];
    //     }
    //     REQUIRE(remain_count == (kTotal_count - kKill_count));

    // }
    // wall_time_end = util::get_wall_time();

    // std::cout << "Per-kill sample from lineup with flag for removal: " << std::setprecision(5) << wall_time_end - wall_time_begin << " seconds.\n";




    // Not a good idea (erase slow) >>

    // wall_time_begin = util::get_wall_time();
    
    // for (int tt = 0; tt < kNum_tests; tt++ ){

    //     std::fill(parasite_count, parasite_count+kNum_majority_types, kMajority_type_count);
    //     std::fill(parasite_count+kNum_majority_types, parasite_count+kNum_genotypes, kMinority_type_count);

    //     util::set_random_numbers_uniform(uniform_random_array, kKill_count);

    //     int accumulator = 0;
    //     all_parasites.resize(kTotal_count);
    //     for (int tt = 0; tt < kNum_genotypes; tt++) {
    //         std::fill(  all_parasites.begin() + accumulator, 
    //                     all_parasites.begin() + accumulator + parasite_count[tt],
    //                     tt
    //         );
    //         accumulator += parasite_count[tt];
    //     }

    //     for (int kk = 0; kk < kKill_count; kk++) {

    //         int kill_at = uniform_random_array[kk] * (kTotal_count-kk);

    //         parasite_count[all_parasites.at(kill_at)]--;
    //         all_parasites.erase(all_parasites.begin() + kill_at);

    //         std::cout << kk << "\n";

    //     }

    //     int remain_count = 0;
    //     for (int tt = 0; tt < kNum_genotypes; tt++) {
    //         remain_count += parasite_count[tt];
    //     }
    //     REQUIRE(remain_count == (kTotal_count - kKill_count));

    // }
    // wall_time_end = util::get_wall_time();

    // std::cout << "Per-kill sample from lineup with per-kill lineup reduction: " << std::setprecision(5) << wall_time_end - wall_time_begin << " seconds.\n";

}

TEST_CASE("90.4: Genotype: combination mutation", "[gt][gt_func_mut]") {

    // const std::vector<float> kMutation_probability{0.1, 0.2, 0.3, 0.4};
    const std::vector<float> kMutation_probability{
        0.0001,
        0.0002,
        0.0003,
        0.0004,
        0.0005,
        0.0006
    };
    common::GenotypeManager gm(kMutation_probability);
    const int kNum_bits = common::GenotypeManager::kNum_bits_in_genotype;

    const int kNum_hosts = 10000000;
    const float kTolerance = 0.1;
    std::vector<int> num_expected_mutations(kMutation_probability.size(),0);
    for (size_t pp=0; pp < kMutation_probability.size(); pp++) {
        num_expected_mutations[pp] = kMutation_probability[pp] * kNum_hosts;
    }

    gm.print_all_combined_probabilities();

    std::vector<std::string> mut_methods {
        "fun_mutation_combined"//,
        // "fun_mutation_individual" //this method does not give the desired behaviour
    };

    for (size_t mm = 0; mm < mut_methods.size(); mm++) {
    
        std::vector<common::genotype_t> gt_list(kNum_hosts, 0);

        double wall_time_begin = 0.0;
        double wall_time_end = 0.0;

        int tot_mutations = 0;

        wall_time_begin = util::get_wall_time();
        for (auto& gg : gt_list) {
            if (mm == 0) {
                tot_mutations += gm.fun_mutation_combined(gg);
            } else {
                tot_mutations += gm.fun_mutation_individual(gg);
            }
        }
        wall_time_end = util::get_wall_time();
        std::cout << "Method: " << mut_methods[mm] << "\n";
        std::cout << "Time taken for " << kNum_hosts
                    << " hosts to complete mutation test per day:"
                    << std::setprecision(5) << wall_time_end - wall_time_begin << " seconds.\n";

        // Check per-position mutation
        std::vector<int> num_mutations_at_position(kMutation_probability.size(), 0);
        int tot_mutations_check = 0;
        for (const auto& gg : gt_list) {
            for (size_t pp = 0; pp < num_mutations_at_position.size(); pp++){
                if (std::bitset<kNum_bits>(gg)[pp]) {
                    num_mutations_at_position[pp]++;
                    tot_mutations_check++;
                }
            }
        }

        std::cout << "position" << " : number of mutations (expected +/- tolerance)\n";
        for (size_t pp = 0; pp < num_mutations_at_position.size(); pp++){
            std::cout << pp << " : "
                << num_mutations_at_position[pp]
                << " (" << num_expected_mutations[pp] <<" +/- "<< kTolerance*100 <<"%)";
            REQUIRE(num_mutations_at_position[pp] - num_expected_mutations[pp]
                < num_expected_mutations[pp] * kTolerance);
            std::cout << " ... pass\n";
        }

        std::cout << "tot_mutations : " << tot_mutations;
        REQUIRE(tot_mutations == tot_mutations_check);
        std::cout << " ... match\n";

        // Check genotype distribution
        std::vector<int> num_hosts(pow(2,int(kMutation_probability.size())), 0);
        // std::cout << "int(kMutation_probability.size():" << int(kMutation_probability.size()) << "\n";
        // std::cout << num_hosts.size() << "\n";
        for (const common::genotype_t& gg : gt_list) {
            num_hosts[unsigned(gg)]++;
        }
        for (size_t gg = 0; gg < num_hosts.size(); gg++) {
            std::cout << "gt" << gg << "-" << num_hosts[gg] << ", ";
        }
        std::cout << "\n";

    }
    
}

TEST_CASE("90.5: Genotype: gt_has_mutations", "[gt][gt_has_mutations]") {

    std::vector<common::ParasiteType> true_pair_gt {
        common::ParasiteType::kKNY__C1,
        common::ParasiteType::kKNY__C2,
        common::ParasiteType::kKYY__C1,
        common::ParasiteType::kKNY__Y2,

        common::ParasiteType::kKNY__Y2,
        common::ParasiteType::kKNY__Y2,
        common::ParasiteType::kKNY__Y2,

        common::ParasiteType::kTNY__Y1,
        common::ParasiteType::kTNY__Y1,
        common::ParasiteType::kTNY__Y1

    };
    std::vector<common::ParasiteType> true_pair_mt {
        common::ParasiteType::kKNY__C1,
        common::ParasiteType::kKNY__C2,
        common::ParasiteType::kKYY__C1,
        common::ParasiteType::kKNY__Y2,
        
        common::ParasiteType::kKNY__Y1,
        common::ParasiteType::kKNY__C2,
        common::ParasiteType::kKNY__C1,
        
        common::ParasiteType::kTNY__Y1,
        common::ParasiteType::kKNY__Y1,
        common::ParasiteType::kTNY__C1
    };

    for (size_t ii = 0; ii < true_pair_gt.size(); ii++) {
        REQUIRE(common::GenotypeManager::gt_has_mutations(
            static_cast<common::genotype_t>(true_pair_gt[ii]),
            static_cast<common::genotype_t>(true_pair_mt[ii])
        ));
    }

    std::vector<common::ParasiteType> false_pair_gt {
        common::ParasiteType::kKNY__C1,
        common::ParasiteType::kKNY__C2,
        common::ParasiteType::kKYY__C2,
        common::ParasiteType::kKNY__Y2,

        common::ParasiteType::kKNY__Y2,
        common::ParasiteType::kKNY__Y2,
        common::ParasiteType::kKNY__Y2,

        common::ParasiteType::kTNY__Y1,
        common::ParasiteType::kTNY__Y1,
        common::ParasiteType::kTNY__Y1

    };
    std::vector<common::ParasiteType> false_pair_mt {
        common::ParasiteType::kKNY__C2,
        common::ParasiteType::kKNY__Y1,
        common::ParasiteType::kKYYYYC1,
        common::ParasiteType::kKYY__Y2,
        
        common::ParasiteType::kTNY__Y2,
        common::ParasiteType::kTNY__C2,
        common::ParasiteType::kKNYNYC1,
        
        common::ParasiteType::kTNYNYY1,
        common::ParasiteType::kKNY__Y2,
        common::ParasiteType::kTNY__C2
    };

    for (size_t ii = 0; ii < true_pair_gt.size(); ii++) {
        REQUIRE(!common::GenotypeManager::gt_has_mutations(
            static_cast<common::genotype_t>(false_pair_gt[ii]),
            static_cast<common::genotype_t>(false_pair_mt[ii])
        ));
    }



}

#include "third_party/catch2/catch.hpp"
#include <iostream>

#include <algorithm>

#include "common/common.h" 

#include "village/village.h"

#include "human/human.h"


TEST_CASE( "110.1: Forest Goers Init test: init_get_villagers", "[village:init_get_villagers]" ) {

    /* 
     This is to test whether new VillageManager 'init_get_villagers' function set the numbers of forest_goers
     stochastically and correctly within the accepted range specified by the forest_goer_percentage_of each village. 
     This is done by counting the numbers forest_goers (through 'get_is_forest_goer' method) for each village and then compare
     their percentage with the exact (no-stochastic) numbers.
     The test requires the difference between the model and exact calculation to be lower than a tolerated (acceptable) rate (e.g. 0.09)
    */

    village::VillageManager vll_mgr(
            true,    // const bool if_batch_mode,
            "./data/test/regions_forest01.txt", //const std::string input_file_name,
            '\t', //const char input_file_delimiter,
            20, //1500, const int rows_to_read,
            0, //const int rows_to_skip,
            true, //const bool init_overwrite_ra_rate_if,
            true, //const float init_overwrite_ra_rate_with,
            
            "mp", //const std::string func_treatment_general_version,

            3, //const float G_timecomp,
            0.9, //const float G_fullcourse,
            0.7, //const float G_covab,
            0.1, //const float G_nomp,
            0.01, //const float asymptomatic_to_clinical_ratio,
            false, //const float establish_malaria_post_during_mda,
            false, //const float establish_malaria_post_during_survey,
            0, //const int establish_malaria_post_on_day,

            true, //const bool malaria_post_rate_reduction_enabled,
            0.5, //const float malaria_post_rate_reduction_size,
            10, //const int malaria_post_rate_reduction_start_on_day,
            30, //const int malaria_post_rate_reduction_slope_days,
            335, //const int malaria_post_rate_reduction_length_days,

            0.5333, //const float treatment_rate_mda,
            1800, //const int mda_max_age,
            0, //const int mda_min_age,

            0.6, //const float treatment_coverage_mmc_wp2,
            3, //const int treatment_drug_mmc_wp2_id,

            true, //const bool init_with_uniform_prevalence_if,
            0.1, //const float init_with_uniform_prevalence_with,
            5.1, //const float init_transmission_coefficient_scaling_factor_a,
            0.75, //const float init_transmission_coefficient_scaling_factor_b,
            1, //const float init_infected_mosquito_to_infected_human_ratio,
            1, //const float init_infectious_mosquito_to_infected_human_ratio,
            70, //const float transmission_coefficient_to_beta_scaling_factor,

            true, //const bool mosquito_seasonality_switched_on,
            false, //const bool mosquito_seasonality_init_if_from_file,
            "./data/data/seasonality_100years4_avg_0.7.txt", //const std::string mosquito_seasonality_init_file_name,
            40, //const int mosquito_seasonality_init_file_horizontal_shift_days,
            0.0, //const float seasonality_base_percent,
            1.0, //const float seasonality_amplitude,

            0.5, //const float seasonality_amplitude_multiplier,
            2.0, //const float cos_pi_multiplier,
            90.0, //const float cos_days_offset,

            2, //const int func_infection_m2h_version,
            1.0, //const float infection_susceptibility_multiplier,
            1.0, //const float infection_infectiousness_multiplier,

            0.5, //const float mosquito_biting_rate,
            14.0, //const float mosquito_incubation_days,
            20.0, //const float mosquito_infected_days,
            7.0, //const float mosquito_infectious_days,

            0, //const float mosquito_min_filter_threshold,
            365000, //const int mosquito_min_filter_start_on_day,

            // Mobility parameters
            false, //const bool mobility_enabled, 
            0.0685, //const float static_population_move_out_probability,
            0.3333, //const float static_population_return_home_probability,
            0.0685, //const float non_static_population_move_out_probability,
            0.3333, //const float non_static_population_return_home_probability

            // Mobility forest 
            true, //const bool mobility_forest_enabled,
            0.5, //const float forest_precedence_over_static_mobility,
            0.5, //const float forest_set_off_probability,
            0.5, //const float forest_return_home_probability
            true,  // const bool stay_at_entry_points 
            true, //return_stay_at_entry_points
            1,    // const int func_mobility_forest_version
            true, //const bool return_early,
            2,   // const int min_num_way_points 
            false //const bool follow_routes
        ); 
    
    human::HumanManager hmn_mgr(vll_mgr.sum_total_population);

    int num_of_villages = vll_mgr.sum_num_villages;

    vll_mgr.init_get_villagers(&hmn_mgr);

    int* should_be_forest_goers_nb_per_village = new int[num_of_villages];

    for (int vv = 0; vv < num_of_villages; vv++){
        int pop = vll_mgr.get_home_population_of_village(vv);
        should_be_forest_goers_nb_per_village[vv] = pop*vll_mgr.forest_goer_percentage_of[vv];
    }

    int* calculated_forest_goers_nb_per_village = new int[num_of_villages];
    float* calculated_forest_goers_pc_per_village = new float[num_of_villages];
    village::VillageManager::Register_t reg = vll_mgr.at_village_register;

    for (int vv = 0; vv < num_of_villages; vv++){
        village::VillageManager::Register_row_t hum_reg_for_this_village= reg.at(vv);
        int pop_for_this_village = vll_mgr.get_home_population_of_village(vv);

        int total_forest_goers_for_this_village =0;
        for ( int ii = 0; ii < pop_for_this_village; ii++){
            if (hmn_mgr.get_is_forest_goer(hum_reg_for_this_village[ii])) {
                total_forest_goers_for_this_village++;
            }
        }
        calculated_forest_goers_nb_per_village[vv]=total_forest_goers_for_this_village;
    
        calculated_forest_goers_pc_per_village[vv]=((float)total_forest_goers_for_this_village/ (float)pop_for_this_village);
        if (pop_for_this_village ==0) {
            calculated_forest_goers_pc_per_village[vv]=0;
        }
    }

    float rate_tolerance = 0.09;
    
    for (int vv = 0; vv < num_of_villages; vv++){
        std::cout.precision(5);
        std::cout << vv << ", forest goers #, should be: "<< should_be_forest_goers_nb_per_village[vv] <<" vs model calculated: "<< calculated_forest_goers_nb_per_village[vv] << ", diff %: "<< calculated_forest_goers_pc_per_village[vv]<<"\n";
        float rate_diff = std::abs(vll_mgr.forest_goer_percentage_of[vv] - calculated_forest_goers_pc_per_village[vv]);
        //std::cout << "rate_diff= "  << rate_diff <<"\n";
        REQUIRE(rate_diff < rate_tolerance);
    }
    
}

TEST_CASE( "110.2: Forest Goers set test: set_human_is_forest_goer", "[human_manager:set_human_is_forest_goer]" ) {
    /* 
     This is to test new HumanManager 'set_human_is_forest_goer' function by creating ad-hoc village_register 
     (along with 3 villages and their desirable forest_goer numbers percentages).
     The number of these forest-goer agents are then calcuated again for each village by going through villlage_register 
     and are required to be the same as those calcuated by multiplying village population to their specified 'forest_goer_percentage_of' each village
    
    */
    
    const int num_villages = 3;
    int village_population[] = {100,0,10};
    int num_humans = 0;
    for (int vv = 0; vv < num_villages; vv++){
        num_humans += village_population[vv];
    }

    float forest_goer_percentage_of[] = {0.10, 0, 0.20};

    human::HumanManager hmn_mgr(num_humans);

    int* should_be_forest_goers_nb_per_village = new int[num_villages];
    int* actual_forest_goers_nb_per_village = new int[num_villages];

    village::VillageManager::Register_t at_village_register;

    int human_index = 0;
    for (int vv = 0; vv < num_villages; vv++) {

        at_village_register.resize(num_villages);
        int max_forest_goers = village_population[vv]*forest_goer_percentage_of[vv];
        int total_forest_goers_for_village=0;
        
        for (int hh = 0; hh < village_population[vv]; hh++) {
            at_village_register[vv].push_back(human_index);

            hmn_mgr.set_human_at_village(human_index, vv);
            hmn_mgr.set_human_home_village(human_index, vv);

            if (total_forest_goers_for_village < max_forest_goers) {
                hmn_mgr.set_human_is_forest_goer(human_index, true);
                total_forest_goers_for_village++;
            } else {
                hmn_mgr.set_human_is_forest_goer(human_index, false);
            }

            human_index++;
        }
    }
    REQUIRE(human_index == hmn_mgr.sum_num_humans);

    for (int vv = 0; vv < num_villages; vv++){
        village::VillageManager::Register_row_t hum_reg_for_this_village= at_village_register.at(vv);
        int pop_for_this_village = village_population[vv];

        should_be_forest_goers_nb_per_village[vv] = pop_for_this_village*forest_goer_percentage_of[vv];

        int total_forest_goers_for_this_village =0;
        for ( int ii = 0; ii < pop_for_this_village; ii++){
            if (hmn_mgr.get_is_forest_goer(hum_reg_for_this_village[ii])) {
                total_forest_goers_for_this_village++;
            }
        }
        actual_forest_goers_nb_per_village[vv] = total_forest_goers_for_this_village;
        std::cout << vv << ", forest goers #, should be: "<< should_be_forest_goers_nb_per_village[vv] <<" vs actual # "<< actual_forest_goers_nb_per_village[vv] <<"\n";
        REQUIRE(should_be_forest_goers_nb_per_village[vv] == actual_forest_goers_nb_per_village[vv]);
    }
}
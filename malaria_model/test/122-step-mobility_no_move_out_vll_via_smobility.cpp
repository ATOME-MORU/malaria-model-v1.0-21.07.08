#include "third_party/catch2/catch.hpp"
#include <iostream>
#include <limits> //std::numeric_limits

#include <algorithm>

#include <vector>

#include <cmath>

#include "common/common.h" 

#include "village/village.h"
#include "common/common.h"

#include "util/randomness.h"

#include "human/human.h"

#include <map>

TEST_CASE("122.1: Step_mobility Functions: step_mobility", "[step_mobility:no_move_outside_village_via_smobility]" ) {

    /*
    This test is desinged to verify that no agent will visit outside village via static mobility.
    Function 'step_population_movement_static()' is used for this test. 
    The test can be done using different number of simulated days (e.g. 1 month/ 1 year, etc)
    and using differnt 'region' file (including list of different villages, forest entries/ POIs)
    */

    int nb_of_simulated_days = 30;
    //int nb_of_simulated_days = 360;

    const std::string villages_forests_file =  "./data/test/regions_forest01.txt";
    const int rows_to_read = 20; //this must be selected for regions_forest01.txt 

    //const std::string villages_forests_file =  "./data/test/regions_forest_small_test.txt";
    //const int rows_to_read = 5; //this must be selected for regions_forest_small_test.txt 


   //===========================

    village::VillageManager vll_mgr(
            true,    // const bool if_batch_mode,
            villages_forests_file, //"./data/test/regions_forest_small_test.txt", //const std::string input_file_name,
            '\t', //const char input_file_delimiter,
            rows_to_read, //20, //1500, const int rows_to_read,
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

    vll_mgr.print_summary();
    //vll_mgr.print_all();

    human::HumanManager hmn_mgr(vll_mgr.sum_total_population);
    vll_mgr.init_get_villagers(&hmn_mgr);

    std::cout << "print village distances:" << "\n";
    vll_mgr.print_village_distances();

    bool an_agent_moved_outside_village_via_static_mobility = false;
    
    std::cout << "\n ========== BEFORE call to step_population_movement_static/ simulated days start/:" << "\n";

    for (int vv = 0; vv < vll_mgr.sum_num_villages; vv++){      
        std::cout << " village: " << vv<< ", type: " <<  vll_mgr.get_location_type(vv);
        std::cout <<", village_reg[i].size: " <<  vll_mgr.at_village_register.at(vv).size();
        std::cout << " , village pop (home): " <<  vll_mgr.get_home_population_of_village(vv)<< ", (current): "<< vll_mgr.get_current_population_of_village(vv)<<"\n";
        std::cout <<"  -Agents(*forest-goer): ";
        for (unsigned int ii = 0; ii < vll_mgr.at_village_register.at(vv).size(); ii++){
            std::cout <<  vll_mgr.at_village_register.at(vv)[ii];
            if (hmn_mgr.get_is_forest_goer(vll_mgr.at_village_register.at(vv)[ii])) {
                std::cout << "*";
            }
            std::cout << ", ";
        }
        std::cout <<"\n";
    }

    for (int dd = 0; dd < nb_of_simulated_days; dd++){     

        vll_mgr.step_population_movement_static();

        for (int vv = 0; vv < vll_mgr.sum_num_villages; vv++){      
            for (unsigned int ii = 0; ii < vll_mgr.at_village_register.at(vv).size(); ii++){
                int agentID = vll_mgr.at_village_register.at(vv)[ii];
                if ((vll_mgr.get_location_type(vv) == common::LocationType::kForestPOI) || 
                   ((vll_mgr.get_location_type(vv) == common::LocationType::kForestEntry))){    
                    std::cout << "forestEntry/ forestPOI ("<< vv<< "): an agent visited, id: " << agentID << "\n";   
                    an_agent_moved_outside_village_via_static_mobility= true;                    
                }
            }
        }
    }


    REQUIRE(an_agent_moved_outside_village_via_static_mobility == false);
}
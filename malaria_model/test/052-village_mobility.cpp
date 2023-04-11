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

TEST_CASE( "52.1: Mobility Functions: step_population_movement_static_kernel", "[mobility:step_movement_static]" ) {
    
    village::VillageManager vll_mgr(
            true,    // const bool if_batch_mode,
            "./data/test/regions_forest01.txt", //const std::string input_file_name,
            '\t', //const char input_file_delimiter,
            20, //20, //1500, const int rows_to_read,
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
            2,    // const int func_mobility_forest_version
            true, //const bool return_early,
            2,   // const int min_num_way_points 
            false //const bool follow_routes

        ); 

    std::cout << "total villages: " <<  vll_mgr.sum_num_villages<<"\n";
    std::cout << "total population: " <<  vll_mgr.sum_total_population<<"\n";

    for (int vv = 0; vv < vll_mgr.sum_num_villages; vv++){
        std::cout << vv<< " , village pop (home): " <<  vll_mgr.get_home_population_of_village(vv)<< ", current: "<< vll_mgr.get_current_population_of_village(vv)<<"\n";
    }

    human::HumanManager hmn_mgr(vll_mgr.sum_total_population);

    vll_mgr.init_get_villagers(&hmn_mgr);

    village::VillageManager::Register_t reg = vll_mgr.at_village_register;

    for (int vv = 0; vv < vll_mgr.sum_num_villages; vv++){
        village::VillageManager::Register_row_t hum_reg_for_this_village= reg.at(vv);
        std::cout <<vv << ", registered population: " <<  hum_reg_for_this_village.size()<<"\n";
    }

    std::vector<std::vector<int>> villager_reg;

    villager_reg.resize(vll_mgr.sum_num_villages);

    int* human_home_village_array = new int[vll_mgr.sum_total_population];

    common::HumanMobilityType* human_mobility_type_array = new common::HumanMobilityType[vll_mgr.sum_total_population];

    int human_index = 0;

    for (int vv = 0; vv < vll_mgr.sum_num_villages; vv++){       
        for (int hh = 0; hh < vll_mgr.get_home_population_of_village(vv); hh++) {
            villager_reg[vv].push_back(human_index);
            human_home_village_array[human_index] = vv;
            human_mobility_type_array[human_index] = common::HumanMobilityType::kStatic;
            human_index++;
        }
    }

    std::vector<std::vector<int>> mobility_network_static;
    mobility_network_static.resize(vll_mgr.sum_num_villages);
    
    mobility_network_static[0].push_back(0);
    mobility_network_static[1].push_back(0);
    mobility_network_static[2].push_back(0);
    mobility_network_static[3].push_back(0);
    mobility_network_static[4].push_back(0);
    mobility_network_static[5].push_back(0);
    mobility_network_static[6].push_back(0);
    mobility_network_static[7].push_back(0);
    mobility_network_static[8].push_back(0);
    mobility_network_static[9].push_back(0);
    mobility_network_static[10].push_back(0);
    mobility_network_static[11].push_back(0);
    mobility_network_static[12].push_back(0);
    mobility_network_static[13].push_back(0);
    mobility_network_static[14].push_back(0);
    mobility_network_static[15].push_back(0);
    mobility_network_static[16].push_back(0);
    mobility_network_static[17].push_back(0);
    mobility_network_static[18].push_back(0);
    mobility_network_static[19].push_back(0); 

    float static_population_move_out_rate = 1;
    float static_population_return_home_rate = 1;
    float non_static_population_move_out_rate = 0;
    float non_static_population_return_home_rate = 0;

    std::cout << "------- Before call to step_population_movement_static_kernel:" << "\n";

    for (int vv = 0; vv < vll_mgr.sum_num_villages; vv++){
        std::cout << vv<<", village_reg[i].size: " <<  villager_reg.at(vv).size() <<  "\n";
    }

   vll_mgr.step_population_movement_static_kernel(
        vll_mgr.sum_num_villages, //num_villages,
        villager_reg,  //vll_mgr.at_village_register, 
        mobility_network_static,
        vll_mgr.sum_total_population, //num_humans,
        human_mobility_type_array,
        human_home_village_array,
        static_population_move_out_rate,    // G_mobility_static
        static_population_return_home_rate, // G_static_return
        non_static_population_move_out_rate,   // G_mobility_mobile
        non_static_population_return_home_rate // G_mobile_return
    ); 


    std::cout << "\n ++++++++ After call to step_population_movement_static_kernel:" << "\n";

    for (int vv = 0; vv < vll_mgr.sum_num_villages; vv++){
        //std::cout << vv<<", village_reg[i].size: " <<  vll_mgr.at_village_register.at(vv).size() << ", village_population[i] (home): " << vll_mgr.get_home_population_of_village(vv) << ", (current): " << vll_mgr.get_current_population_of_village(vv) << "\n";
        std::cout << vv<<", village_reg[i].size: " <<  villager_reg.at(vv).size() << "\n";
    }

    REQUIRE(villager_reg[0].size()==vll_mgr.sum_total_population);

    mobility_network_static.clear();
    mobility_network_static.resize(vll_mgr.sum_num_villages);

    mobility_network_static[0].push_back(1);  //pop. movement of village 0 to village 1
    mobility_network_static[1].push_back(0);
    mobility_network_static[2].push_back(0);
    mobility_network_static[3].push_back(0);
    mobility_network_static[4].push_back(0);
    mobility_network_static[5].push_back(0);
    mobility_network_static[6].push_back(0);
    mobility_network_static[7].push_back(0);
    mobility_network_static[8].push_back(0);
    mobility_network_static[9].push_back(0);
    mobility_network_static[10].push_back(0);
    mobility_network_static[11].push_back(0);
    mobility_network_static[12].push_back(0);
    mobility_network_static[13].push_back(0);
    mobility_network_static[14].push_back(0);
    mobility_network_static[15].push_back(0);
    mobility_network_static[16].push_back(0);
    mobility_network_static[17].push_back(0);
    mobility_network_static[18].push_back(0);
    mobility_network_static[19].push_back(0); 

    std::cout << "\n \n S2 ------- Before call to step_population_movement_static_kernel:" << "\n";

    for (int vv = 0; vv < vll_mgr.sum_num_villages; vv++){
        std::cout << vv<<", village_reg[i].size: " <<  vll_mgr.at_village_register.at(vv).size() << ", village_population[i] (home): " << vll_mgr.get_home_population_of_village(vv) << ", (current): " << vll_mgr.get_current_population_of_village(vv) << "\n";
    }

   vll_mgr.step_population_movement_static_kernel(
        vll_mgr.sum_num_villages, //num_villages,
        vll_mgr.at_village_register, //villager_reg,
        mobility_network_static,
        vll_mgr.sum_total_population, //num_humans,
        human_mobility_type_array,
        human_home_village_array,
        static_population_move_out_rate,    // G_mobility_static
        static_population_return_home_rate, // G_static_return
        non_static_population_move_out_rate,   // G_mobility_mobile
        non_static_population_return_home_rate // G_mobile_return
    ); 

    std::cout << "\n S2 ++++++++ After call to step_population_movement_static_kernel:" << "\n";

    for (int vv = 0; vv < vll_mgr.sum_num_villages; vv++){
        std::cout << vv<<", village_reg[i].size: " <<  vll_mgr.at_village_register.at(vv).size() << ", village_population[i] (home): " << vll_mgr.get_home_population_of_village(vv) << ", (current): " << vll_mgr.get_current_population_of_village(vv) << "\n";
    }

    REQUIRE(vll_mgr.at_village_register.at(1).size()==vll_mgr.get_home_population_of_village(0)+ vll_mgr.get_home_population_of_village(1));

    delete[] human_home_village_array;
    delete[] human_mobility_type_array;
}
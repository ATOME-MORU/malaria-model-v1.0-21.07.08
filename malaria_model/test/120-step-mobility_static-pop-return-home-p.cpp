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

TEST_CASE("120.1: Step_mobility Functions: step_mobility", "[step_mobility:static_pop_return_home_p]" ) {

    /*
    This test aims to help evaluating the effects of value change for 'static_population_return_home_probability' parameter.
    It generates a simple simulation by calling new VillageManager 'step_mobility' function at desirable and settable 
    number of days (e.g. 30 or 360 days, etc.)
    It currenlty works using a simplified region file version caled 'regions_forest_small_test.txt', so that the movements only
    exist between the first village and the only single forest in the file. The number of population in village
    are also small (e.g. 30 agents).
    The test keeps track of number of days an agent spend in the forest and provides some statitics at the 
    end of the simulation period, including the last longest days spent in forest by a given agent. 
    */

    //const float static_population_return_home_probability = 0.50; // 1/2
    const float static_population_return_home_probability = 0.3333; // 1/3
    //const float static_population_return_home_probability = 0.25; // 1/4
    //const float static_population_return_home_probability = 0.2; // 1/5
    //const float static_population_return_home_probability = 0.1111; // 1/9

    int num_of_days = 30;
    //int num_of_days = 360;

    village::VillageManager vll_mgr(
            true,    // const bool if_batch_mode,
            "./data/test/regions_forest_small_test.txt", //const std::string input_file_name,
            '\t', //const char input_file_delimiter,
            5, //20, //1500, const int rows_to_read,
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
            static_population_return_home_probability, //0.3333, //const float static_population_return_home_probability,
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

    //std::cout << "total villages: " <<  vll_mgr.sum_num_villages<<"\n";
    //std::cout << "total population: " <<  vll_mgr.sum_total_population<<"\n";

    vll_mgr.print_summary();
    //vll_mgr.print_all();

    human::HumanManager hmn_mgr(vll_mgr.sum_total_population);
    vll_mgr.init_get_villagers(&hmn_mgr);

    /*for (int vv = 0; vv < vll_mgr.sum_num_villages; vv++){
        std::cout << vv<< " , village pop (home): " <<  vll_mgr.get_home_population_of_village(vv)<< ", current: "<< vll_mgr.get_current_population_of_village(vv)<<"\n";
    } */

    /*std::cout << "Villages info: " <<"\n";
    for (int vv = 0; vv < vll_mgr.sum_num_villages; vv++){
        village::VillageManager::Register_row_t hum_reg_for_this_village= vll_mgr.at_village_register.at(vv);
        std::cout <<vv << ", registered population: " <<  hum_reg_for_this_village.size()<<"\n";
        //int vll_pop = hum_reg_for_this_village.size();
        //int vll_pop = vll_mgr.get_home_population_of_village(vv);
        int vll_pop = vll_mgr.get_current_population_of_village(vv);
        for ( int ii = 0; ii < vll_pop; ii++){
            std::cout <<  hum_reg_for_this_village[ii]<< ", ";
        }
        std::cout <<"\n";
    } */

    //std::cout << "Print forest entry network:" << "\n";
    //vll_mgr.print_forest_entry_network();

    //std::cout << "Print forest poi network:" << "\n";
    //vll_mgr.print_forest_poi_network();

    std::cout << "print village distances:" << "\n";
    vll_mgr.print_village_distances();

    //std::cout << "update summary:" << "\n";
    //vll_mgr.update_summary();
    
    std::cout << "\n ====================== BEFORE call to step_mobility:" << "\n";

    for (int vv = 0; vv < vll_mgr.sum_num_villages; vv++){      
        std::cout << " village: " << vv<< ", type: " <<  vll_mgr.get_location_type(vv);
        std::cout <<", village_reg[i].size: " <<  vll_mgr.at_village_register.at(vv).size();
        std::cout << " , village pop (home): " <<  vll_mgr.get_home_population_of_village(vv)<< ", (current): "<< vll_mgr.get_current_population_of_village(vv)<<"\n";
        std::cout <<"  -Agents(*forest-goer): ";
        for (unsigned int ii = 0; ii < vll_mgr.at_village_register.at(vv).size(); ii++){
            std::cout <<  vll_mgr.at_village_register.at(vv)[ii];
            if (hmn_mgr.get_is_forest_goer(vll_mgr.at_village_register.at(vv)[ii]))
                std::cout << "*";
            std::cout << ", ";
        }
        std::cout <<"\n";
    }
     
    //==================

    std::map<int, int> f_map;
    std::vector<int> in_forest_list;
    int total_nb_of_returned_agents=0;
    int total_nb_of_days_in_forest_for_returned_agents=0;
    int longest_stay_in_forest = -1;
    int shortest_stay_in_forest = 99999;
    int day_longest_stay_in_forest_occured = -1;
    int agent_id_of_longest_stay_in_forest = -1;
    
    for (int dd = 0; dd < num_of_days; dd++){     

        vll_mgr.step_mobility();

        std::cout << "\n Day: "<< dd << "  ======== AFTER call to step_mobility:" << "\n";
        for (int vv = 0; vv < vll_mgr.sum_num_villages; vv++){           
            std::cout << " village: " << vv<< ", type: " <<  vll_mgr.get_location_type(vv);
            std::cout <<", village_reg[i].size: " <<  vll_mgr.at_village_register.at(vv).size();            
            std::cout << " , village pop (home): " <<  vll_mgr.get_home_population_of_village(vv)<< ", (current): "<< vll_mgr.get_current_population_of_village(vv)<<"\n";
            std::cout <<"  -Agents(*forest-goer): ";
            for (unsigned int ii = 0; ii < vll_mgr.at_village_register.at(vv).size(); ii++){
                std::cout <<  vll_mgr.at_village_register.at(vv)[ii];
                if (hmn_mgr.get_is_forest_goer(vll_mgr.at_village_register.at(vv)[ii]))
                    std::cout << "*";
                std::cout << ", ";
        }
            std::cout <<"\n";
        }
        std::cout <<"-------\n";
        for (int vv = 0; vv < vll_mgr.sum_num_villages; vv++){           

            if (vll_mgr.get_location_type(vv) == common::LocationType::kForestPOI) {

                in_forest_list.clear();
                for (unsigned int jj = 0; jj < vll_mgr.at_village_register.at(vv).size(); jj++){
                    int agentID = vll_mgr.at_village_register.at(vv)[jj];
                    //std::cout <<"+ ID added to list: " <<agentID <<"\n";
                    in_forest_list.push_back(agentID);
                }
                for (auto it = begin(in_forest_list); it != end(in_forest_list); ++it) {
                    int agentID = *it;
                    if (!f_map.count(agentID)) { //map does not contain ID
                        f_map.insert(std::pair<int,int>(agentID,0));
                        //std::cout << "added ID to map with 0 day: " <<agentID <<"\n";
                    }
                }                                        

                for (std::map<int,int>::iterator it = f_map.begin(); it != f_map.end(); ) {
                    int agentID = it->first;
                    int nb_days_in_forest = it->second;
                    //std::cout <<" map: does map contains key?: " <<agentID;

                    if (std::find(in_forest_list.begin(), in_forest_list.end(), agentID) != in_forest_list.end()){
                        it->second++;
                        nb_days_in_forest++;
                        //std::cout << ", ID: " <<it->first<<" is present in map for days: "<< it->second<<"\n";
                        std::cout << "+ AgentID: " << agentID <<" is present in map for days: "<< nb_days_in_forest<<"\n";
                        it++;
                    }
                    else { //does not contain
                        //std::cout << "Deleted agentID: " <<it->first<<" WAS present in map for days: "<< it->second<<"\n";
                        std::cout << "- Deleted agentID: " <<agentID <<" WAS present in map for days: "<< nb_days_in_forest <<"\n";

                        f_map.erase(it++);

                        if (nb_days_in_forest>longest_stay_in_forest) {
                            longest_stay_in_forest = nb_days_in_forest;
                            day_longest_stay_in_forest_occured= dd;
                            agent_id_of_longest_stay_in_forest=agentID;
                        }
                        if (nb_days_in_forest<shortest_stay_in_forest) {
                            shortest_stay_in_forest = nb_days_in_forest;
                        }
                        total_nb_of_returned_agents++;
                        total_nb_of_days_in_forest_for_returned_agents += nb_days_in_forest;
                    }                 
                }
            } 
        }
    }
    std::cout << "\n ========> Summary report: during " <<num_of_days<< " days  <=========:\n";
    std::cout << "[Expected average of each trip to forest to last about: " <<  1/static_population_return_home_probability << " days]\n";
    std::cout << "Total nb of returned agents: " <<total_nb_of_returned_agents<<"\n";
    std::cout <<"Total nb of days returned agents spent in forest: "<< total_nb_of_days_in_forest_for_returned_agents<<"\n";
    float avg_days_in_forest_for_returned_agents = (float)total_nb_of_days_in_forest_for_returned_agents/(float)total_nb_of_returned_agents;
    std::cout << "Avg days agents spent in forest: " <<avg_days_in_forest_for_returned_agents<< "\n";
    std::cout << "Last shortest stay: " <<shortest_stay_in_forest<<"\n";
    std::cout << "Last longest stay: " <<longest_stay_in_forest;
    std::cout << " (by agent ID: " <<agent_id_of_longest_stay_in_forest;
    std::cout << ", who was returned home on day: " <<day_longest_stay_in_forest_occured<<")\n";

    //REQUIRE(avg_days_in_forest_for_returned_agents <= 1/static_population_return_home_probability);
    //REQUIRE(longest_stay_in_forest <= 1/static_population_return_home_probability);

}
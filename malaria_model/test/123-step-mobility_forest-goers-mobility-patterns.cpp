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

TEST_CASE("123.1: Step_mobility Functions: step_mobility", "[step_mobility:forest_goers_mobility_patterns]" ) {

    /*
    This test aims to track and provide a short statitic on forest_goers mobility pattern (e.g. nb of days stayed at each location)
    It generates a simple simulation by calling new VillageManager 'step_mobility' function at desirable and settable 
    number of days (e.g. 30 or 360 days, etc.)
    This can be done using the orignal region file (('region_forest_01.txt') or the simplified/smaller version ('region_forest_02.txt').
    The stats also include the the last longest days spent in forestPOI and no forest location. 
    */

    int nb_of_simulated_days = 10;
    //int nb_of_simulated_days = 360;

    const std::string villages_forests_file =  "./data/test/regions_forest01.txt";
    const int rows_to_read = 20; //this must be selected for regions_forest01.txt 

    //const std::string villages_forests_file =  "./data/test/regions_forest_small_test.txt";
    //const int rows_to_read = 5; //this must be selected for regions_forest_small_test.txt

    //const float static_population_return_home_probability = 0.50; // 1/2
    const float static_population_return_home_probability = 0.3333; // 1/3
    //const float static_population_return_home_probability = 0.25; // 1/4
    //const float static_population_return_home_probability = 0.2; // 1/5
    //const float static_population_return_home_probability = 0.1111; // 1/9

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

    vll_mgr.print_summary();
    //vll_mgr.print_all();

    //std::cout << "total villages: " <<  vll_mgr.sum_num_villages<<"\n";
    //std::cout << "total population: " <<  vll_mgr.sum_total_population<<"\n";

    vll_mgr.print_summary();
    //vll_mgr.print_all();

    human::HumanManager hmn_mgr(vll_mgr.sum_total_population);
    vll_mgr.init_get_villagers(&hmn_mgr);

    std::cout << "print village distances:" << "\n";
    vll_mgr.print_village_distances();

    //std::cout << "update summary:" << "\n";
    //vll_mgr.update_summary();

    std::vector<int> agents;
    std::vector<std::vector<int>> agents_itinerary(vll_mgr.sum_total_population); // [agentID][itinerary]

    int h_ind;
    //hmn_mgr.human_reg

    std::cout << "\n ====================== BEFORE call to step_mobility:" << "\n";
    

    for (int vv = 0; vv < vll_mgr.sum_num_villages; vv++){      
        std::cout << " village: " << vv<< ", type: " <<  vll_mgr.get_location_type(vv);
        std::cout <<", village_reg[i].size: " <<  vll_mgr.at_village_register.at(vv).size();
        std::cout << " , village pop (home): " <<  vll_mgr.get_home_population_of_village(vv)<< ", (current): "<< vll_mgr.get_current_population_of_village(vv)<<"\n";
        std::cout <<"  -Agents(*forest-goer): ";

        for (unsigned int ii = 0; ii < vll_mgr.at_village_register.at(vv).size(); ii++){        
            std::cout <<  vll_mgr.at_village_register.at(vv)[ii];
            int agentID = vll_mgr.at_village_register.at(vv)[ii];
            agents_itinerary.at(agentID).push_back(hmn_mgr.get_at_village(agentID));
            //agents_itinerary.at(agentID).push_back( hmn_mgr.home_village_of[agentID]);
            if (hmn_mgr.get_is_forest_goer(vll_mgr.at_village_register.at(vv)[ii]))
                std::cout << "*";
            std::cout << ", ";
        }
        std::cout <<"\n";
    }

    std::cout << "\n ++++Agent itenaries:" << "\n";
    for (unsigned int ii = 0; ii < agents_itinerary.size(); ii++){   
        std::cout << "itinerary for agent: " <<ii <<" is: ";
        for (unsigned int jj = 0; jj < agents_itinerary.at(ii).size(); jj++){       
            std::cout << agents_itinerary.at(ii).at(jj) << ", ";
        }
        std::cout <<"\n";
    }

    // Simulation: 
    for (int dd = 0; dd < nb_of_simulated_days; dd++){     

        vll_mgr.step_mobility();

        for (int vv = 0; vv < vll_mgr.sum_num_villages; vv++){      
            for (unsigned int ii = 0; ii < vll_mgr.at_village_register.at(vv).size(); ii++){        
                int agentID = vll_mgr.at_village_register.at(vv)[ii];
                agents_itinerary.at(agentID).push_back(hmn_mgr.get_at_village(agentID));
            }
        }
    }

    std::cout << "\n ++++AFTER Simulation Agent itenaries :" << "\n";
    for (unsigned int ii = 0; ii < agents_itinerary.size(); ii++){   
        std::cout << "itinerary for agent: " <<ii <<" is: ";
        for (unsigned int jj = 0; jj < agents_itinerary.at(ii).size(); jj++){       
            std::cout << agents_itinerary.at(ii).at(jj) << ", ";
        }
        std::cout <<"\n";
    }
    std::cout <<"\n";
    
    struct stat {
        int loc;
        int stayed_day=1;
        stat(int l) : loc(l) {}
    };

    std::vector<std::vector<stat>> agents_stats(vll_mgr.sum_total_population); // [agentID][stat]

    for (unsigned int ii = 0; ii < agents_itinerary.size(); ii++){   

        unsigned int jj=0;
        std::cout << "making stat for agent " <<ii << ": ";
        int itenary_index =0;

        while(jj < agents_itinerary.at(ii).size()) {  
            bool moved=false; 
            int from_loc = agents_itinerary.at(ii).at(jj);
            agents_stats.at(ii).push_back(stat(from_loc));
            int to_loc = -1;
            
            std::cout <<from_loc << "(";

            while(moved == false)  {
                if (jj < agents_itinerary.at(ii).size()-1) {
                    jj++;
                    to_loc=agents_itinerary.at(ii).at(jj);
                    if (from_loc==to_loc) {
                        agents_stats.at(ii).at(itenary_index).stayed_day++;
                    }
                    else{
                        std::cout <<agents_stats.at(ii).at(itenary_index).stayed_day << "), ";
                        itenary_index++;
                        moved=true;
                    }
                }
                else {
                    std::cout <<agents_stats.at(ii).at(itenary_index).stayed_day << "), ";
                    std::cout <<"\n";
                    itenary_index++;
                    jj++;
                    moved = true;
                } 
            }
        }
    }

    std::cout << "\n ================= Final Summary ===============: \n";

    int nb_of_longest_no_forest_stay_by_agent=-1;
    int agent_id_of_longest_no_forest_stayer=-1;
    int location_of_longest_stay=-1;

    int nb_of_longest_stay_by_moving_agent_at_forest=-1;
    int moving_agent_id_of_longest_stayer_at_forest=-1;
    int location_of_longest_stay_at_forest=-1;

    std::cout << "\n All agents' stats at the end of " << nb_of_simulated_days <<" days simulation (inlduing the original location night):" << "\n";
    std::cout << "  [*: forest-goer, v: village, f: forestPOI, e:forestEntry, (): #of day/night at location]" << "\n";
    std::cout << "  [Expected average of each trip to forest to last about: " <<  1/static_population_return_home_probability << " days]\n \n";

    for (unsigned int ii = 0; ii < agents_stats.size(); ii++){   
        if (hmn_mgr.get_is_forest_goer(ii)) {
            std::cout << "Itinerary for Agent: " <<ii;
            if (hmn_mgr.get_is_forest_goer(ii))
                    std::cout << "*";
            std::cout << " -> ";
            for (unsigned int jj = 0; jj < agents_stats.at(ii).size(); jj++){    
                int loc =   agents_stats.at(ii).at(jj).loc;
                int nb_days =  agents_stats.at(ii).at(jj).stayed_day;
                if (agents_stats.at(ii).size()>1) {
                    if (vll_mgr.get_location_type(loc) == common::LocationType::kForestPOI) {
                        if (nb_days >= nb_of_longest_stay_by_moving_agent_at_forest) {
                            nb_of_longest_stay_by_moving_agent_at_forest=nb_days;
                            moving_agent_id_of_longest_stayer_at_forest=ii;
                            location_of_longest_stay_at_forest=loc;
                        }
                    }
                    else {
                        if (nb_days >= nb_of_longest_no_forest_stay_by_agent) {
                            nb_of_longest_no_forest_stay_by_agent = nb_days;
                            agent_id_of_longest_no_forest_stayer=ii;
                            location_of_longest_stay=loc;
                        }
                    }
                }
                std::cout << vll_mgr.get_location_type(loc) << loc<< "(" <<nb_days<<"), ";
            }
            if (agents_stats.at(ii).size()==1) {
                std::cout << "[immobile agent]";
            }
            else {
                std::cout << "[moving agent]";
            }
            std::cout <<"\n";
        }
    }

    std::cout << "\n --- Last observed longest stays: \n";
    std::cout << "In No-forest location, by agentID: " <<agent_id_of_longest_no_forest_stayer;
    if (hmn_mgr.get_is_forest_goer(agent_id_of_longest_no_forest_stayer)) std::cout << "*";
    std::cout << " [home v: " <<  hmn_mgr.get_home_village(agent_id_of_longest_no_forest_stayer)<<"]";
    std::cout << ", at location: " <<vll_mgr.get_location_type(location_of_longest_stay)<<location_of_longest_stay;
    std::cout << ", for " <<nb_of_longest_no_forest_stay_by_agent <<" days/nights \n";

    std::cout << "In Forest location, by moving agentID: " <<moving_agent_id_of_longest_stayer_at_forest;
    if (hmn_mgr.get_is_forest_goer(moving_agent_id_of_longest_stayer_at_forest)) std::cout << "*";
    std::cout << ", at location: " <<vll_mgr.get_location_type(location_of_longest_stay_at_forest)<<location_of_longest_stay_at_forest;
    std::cout << ", for " <<nb_of_longest_stay_by_moving_agent_at_forest <<" days/nights \n";

}
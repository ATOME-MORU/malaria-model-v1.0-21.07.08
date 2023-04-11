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

TEST_CASE("125.1: Step_mobility Functions: step_mobility", "[step_mobility:forest-precedence-over-static-mobility]" ) {

    /*
    
    This test aims to help to observe the effects of the 'forest_precedence_over_static_mobility' parameter 
    on the dynamic of model/population moves. 

    This can be done by setting different values to this parameter and observe certain statistics related to
    forest and static moves.
    
    The difference between these outputs can help to see the effect of the parameter change. 
           
    The test makes a call to  VillageManager 'step_mobility' function at a desirable and settable 
    number of days (e.g. 30 or 360 days, etc.)
    
    This can be done using any region file.

    The test can be also done by choosing whether the agents 'stay_at_entry_points' or not
    as well as specifying the min_num_way_points'. 
    
    */

    //int nb_of_simulated_days =10;
    int nb_of_simulated_days = 360;
    //int nb_of_simulated_days = 1800;

        //float forest_precedence_over_static_mobility=0.1;
    float forest_precedence_over_static_mobility=0.5;
    //float forest_precedence_over_static_mobility=0.9;

    const bool stay_at_entry_points = true;
    //const bool stay_at_entry_points = false;

    //const int min_num_way_points = 1;
    const int min_num_way_points = 2;
    //const int min_num_way_points = 3;
    //const int min_num_way_points = 4;

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
            true, //const bool mobility_enabled, 
            0.0685, //const float static_population_move_out_probability,
            static_population_return_home_probability, //0.3333, //const float static_population_return_home_probability,
            0.0685, //const float non_static_population_move_out_probability,
            0.3333, //const float non_static_population_return_home_probability

            // Mobility forest 
            true, //const bool mobility_forest_enabled,
            forest_precedence_over_static_mobility, //0.5, //const float forest_precedence_over_static_mobility,
            0.5, //const float forest_set_off_probability,
            0.5, //const float forest_return_home_probability
            stay_at_entry_points, //true,  // const bool stay_at_entry_points 
            true, //return_stay_at_entry_points
            1,    // const int func_mobility_forest_version
            true, //const bool return_early,
            min_num_way_points,  //2 // const int min_num_way_points 
            false //const bool follow_routes
        ); 

    std::cout << "print summary:" << "\n";
    vll_mgr.print_summary();
    //vll_mgr.print_all();

    //std::cout << "total villages: " <<  vll_mgr.sum_num_villages<<"\n";
    //std::cout << "total population: " <<  vll_mgr.sum_total_population<<"\n";

    //vll_mgr.print_summary();
    //vll_mgr.print_all();

    human::HumanManager hmn_mgr(vll_mgr.sum_total_population);
    vll_mgr.init_get_villagers(&hmn_mgr);

    std::cout << "print village distances:" << "\n";
    vll_mgr.print_village_distances();

    std::cout << "print forest entry network:" << "\n";
    vll_mgr.print_forest_entry_network();

    std::vector<std::vector<int>> forest_route_network_lists;
    for (unsigned int ee=0; ee<vll_mgr.forest_route_network.size(); ee++) {
        //std::cout << "\t";
        for (unsigned int rr=0; rr<vll_mgr.forest_route_network.at(ee).size(); rr++) {
            //std::cout << "\n \t";
            std::vector<int> vec;
            for (unsigned int vv=0; vv<vll_mgr.forest_route_network.at(ee).at(rr).size(); vv++) {
                //std::cout << vll_mgr.forest_route_network.at(ee).at(rr).at(vv) << ", ";
                vec.push_back(vll_mgr.forest_route_network.at(ee).at(rr).at(vv));
            }
            forest_route_network_lists.push_back(vec);
        }
        //std::cout << "\n";
    }

    //std::cout << "update summary:" << "\n";
    //vll_mgr.update_summary();

    int total_agents = vll_mgr.sum_total_population;
    int total_villages =0;
    int total_forest_poi=0;
    int total_forest_goers=0;

    std::vector<std::vector<int>> agents_itinerary(vll_mgr.sum_total_population); // [agentID][itinerary]

    std::cout << "\n ====================== BEFORE call to step_mobility:" << "\n";

    for (int vv = 0; vv < vll_mgr.sum_num_villages; vv++){      
        std::cout << " village: " << vv<< ", type: " <<  vll_mgr.get_location_type(vv);
        std::cout <<", village_reg[i].size: " <<  vll_mgr.at_village_register.at(vv).size();
        std::cout << " , village pop (home): " <<  vll_mgr.get_home_population_of_village(vv)<< ", (current): "<< vll_mgr.get_current_population_of_village(vv)<<"\n";
        std::cout <<"  -Agents(*forest-goer): ";

        if (vll_mgr.get_location_type(vv) == common::LocationType::kVillage) 
            total_villages++;
        
        if (vll_mgr.get_location_type(vv) == common::LocationType::kForestPOI) 
            total_forest_poi++;

        for (unsigned int ii = 0; ii < vll_mgr.at_village_register.at(vv).size(); ii++){        
            std::cout <<  vll_mgr.at_village_register.at(vv)[ii];
            int agentID = vll_mgr.at_village_register.at(vv)[ii];
            agents_itinerary.at(agentID).push_back(hmn_mgr.get_at_village(agentID));
            if (hmn_mgr.get_is_forest_goer(vll_mgr.at_village_register.at(vv)[ii])) {
                std::cout << "*";
                total_forest_goers++;
            }
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
    }  // end of simulation 



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

    //vll_mgr.build_forest_route_network();

    std::cout << "\n ================= Final Summary ===============: \n";

    //bool stay_home_between_mobilities = true;

    int nb_of_longest_stay_at_home_vilage=-1; //longest stay by a moving agenet at own home village (after move)
    int agent_id_of_longest_stay_at_home_vilage=-1;
    int location_of_longest_stay_at_home=-1;  

    int nb_of_longest_stay_at_forest=-1;
    int agent_id_of_longest_stay_at_forest=-1;
    int location_of_longest_stay_at_forest=-1;

    bool static_move_took_place = false;  //an agent moved  to another village

    int nb_of_longest_stay_at_another_village=-1;
    int agent_id_of_longest_stay_at_another_village=-1; //agent made a static move (moved to another village)
    int location_of_longest_stay_at_another_village=-1;

    int total_nb_of_forest_moves =0;
    int total_days_spent_outside_home_in_forest_moves=0;

    int total_nb_of_forest_moves_including_in_trans_routes =0;

    int total_nb_of_static_moves =0;
    int total_nb_of_static_moves_by_forest_goers =0;
    int total_days_spent_outside_home_in_static_moves=0;


    std::cout << "\n All agents' stats at the end of " << nb_of_simulated_days <<" days simulation (inlduing the original location night):" << "\n";
    std::cout << "  [*: forest-goer, ^: static move, v: village, f: forestPOI, e:forestEntry, (): #of day/night at location]" << "\n";
    std::cout << "  [Expected average of each trip to forest to last about: " <<  1/static_population_return_home_probability << " days]\n \n";
    std::cout << "  [the first day when the agent is at home is also counted in the calculation of days] \n";

    for (unsigned int hh = 0; hh < agents_stats.size(); hh++){   
        std::cout << "Itinerary for Agent: " <<hh;
        if (hmn_mgr.get_is_forest_goer(hh))
                std::cout << "*";
        std::cout << " -> ";

        std::vector<int> out_of_home_route;
        std::vector<int> out_of_home_route_days;
        bool record_out_of_home_route =false;
 
        for (unsigned int jj = 0; jj < agents_stats.at(hh).size(); jj++){    
            int loc =   agents_stats.at(hh).at(jj).loc;
            int nb_days =  agents_stats.at(hh).at(jj).stayed_day;

            if (agents_stats.at(hh).size()>1) {
                int home_village = hmn_mgr.get_home_village(hh);

                if ((home_village != loc) && (vll_mgr.get_location_type(loc) == common::LocationType::kVillage)) { 
                    std::cout << "^" ; 
                    static_move_took_place=true;

                    if (nb_days >= nb_of_longest_stay_at_another_village) {
                        nb_of_longest_stay_at_another_village=nb_days;
                        agent_id_of_longest_stay_at_another_village=hh;
                        location_of_longest_stay_at_another_village=loc;
                    }
                }

                if (vll_mgr.get_location_type(loc) == common::LocationType::kForestPOI) {
                    if (nb_days >= nb_of_longest_stay_at_forest) {
                        nb_of_longest_stay_at_forest=nb_days;
                        agent_id_of_longest_stay_at_forest=hh;
                        location_of_longest_stay_at_forest=loc;
                    }
                }
                else if (home_village == loc){
                    if (nb_days >= nb_of_longest_stay_at_home_vilage) {
                        nb_of_longest_stay_at_home_vilage = nb_days;
                        agent_id_of_longest_stay_at_home_vilage=hh;
                        location_of_longest_stay_at_home=loc;
                    }
                }

                //-------------------------

                if (home_village != loc)   
                    record_out_of_home_route = true;
                else record_out_of_home_route = false;
                               
                if (record_out_of_home_route) {
                    out_of_home_route.push_back(loc);
                    out_of_home_route_days.push_back(nb_days);
                    //std::cout << " Added: "<<loc<<", ";  
                } 
                //process journey route: 
                else if (out_of_home_route.size() >= 1) { //end of recording and agent moved/ not immobile

                    bool move_static = vll_mgr.get_location_type(out_of_home_route.at(0)) == common::LocationType::kVillage; 

                    if (move_static) {

                        if (out_of_home_route.size() > 1) {
                            //stay_home_between_mobilities=false;
                            std::cout << ">>>>>ERROR: agent must return home after static move (visiting another village) <<<<< \n ";
                        } 

                        /*std::cout << "\n    static move(s: " <<out_of_home_route.size() <<"): ";  
                        for (unsigned int ii=0; ii<out_of_home_route.size(); ii++) 
                            std::cout << out_of_home_route.at(ii) << ", ";
                        std::cout << "\n "; */

                        total_nb_of_static_moves++;
                        total_days_spent_outside_home_in_static_moves += out_of_home_route_days.at(0); //number of days spent at ther village  visisted
                        
                        if (hmn_mgr.get_is_forest_goer(hh))
                            total_nb_of_static_moves_by_forest_goers++;

                        out_of_home_route.clear();
                        out_of_home_route_days.clear();
                        
                    }
                    else {  //forest/ no static move case
                        bool route_exist = false;
                        bool stayed_only_one_day_in_transient_forest_routes = true;
                        
                        if (stay_at_entry_points)
                            out_of_home_route.erase(out_of_home_route.begin()); //remove the entry point if stay_at_entry_points is true
                        
                        for (unsigned int ii=0; ii<forest_route_network_lists.size(); ii++) {
                            if (forest_route_network_lists.at(ii) == out_of_home_route ) 
                                route_exist= true;        
                        }
                        if (!route_exist) {
                            //stay_home_between_mobilities=false;
                            std::cout << ">>>>>ERROR: route does not exist / staying home not observed<<<<< \n ";
                            for (unsigned int ii=0; ii<out_of_home_route.size(); ii++) {
                                    std::cout << out_of_home_route.at(ii) << ", ";
                            }
                            std::cout << "\n ";                                
                        }
                        /*std::cout << "route exist: ";  
                        for (unsigned int ii=0; ii<out_of_home_route.size(); ii++) 
                            std::cout << out_of_home_route.at(ii) << ", ";
                        std::cout << "\n "; */

                        for (unsigned int ii=0; ii<out_of_home_route_days.size()-1; ii++) {
                            if (out_of_home_route_days.at(ii) != 1) // when stay_at_entry_points is true, we still check that agents can stay there only 1 day/night
                                stayed_only_one_day_in_transient_forest_routes = false;        
                        }
                        if (!stayed_only_one_day_in_transient_forest_routes) {
                            //stay_home_between_mobilities=false;
                            std::cout << ">>>>>ERROR: stayed more than one day in one of the transient forest routes <<<<< \n ";
                            for (unsigned int ii=0; ii<out_of_home_route_days.size(); ii++) {
                                    std::cout << out_of_home_route_days.at(ii) << ", ";
                            }
                            std::cout << "\n ";                                
                        }

                        /*std::cout << " trans stay days: ";  
                        for (unsigned int ii=0; ii<out_of_home_route_days.size()-1; ii++) {
                            std::cout << out_of_home_route_days.at(ii) << ", ";
                        }
                        std::cout << "\n "; */

                        total_nb_of_forest_moves++;
                        total_nb_of_forest_moves_including_in_trans_routes += out_of_home_route.size();

                        total_days_spent_outside_home_in_forest_moves += (out_of_home_route_days.size()-1); // for those one nights before last forest stop to home
                        total_days_spent_outside_home_in_forest_moves += out_of_home_route_days.at(out_of_home_route_days.size()-1); // for last.main stop in forest before get home

                        out_of_home_route.clear();
                        out_of_home_route_days.clear();

               
                    }
                }
            } 
            std::cout << vll_mgr.get_location_type(loc) << loc<< "(" <<nb_days<<"), ";
        }
        if (agents_stats.at(hh).size()==1) {
            std::cout << "[immobile agent]";
        }
        else {
            std::cout << "[moving agent]";
        }
        std::cout <<"\n";
    }

    std::cout << "\n --- Last observed longest stay incidences: \n";
    std::cout << "In No-forest/ own home village location, by agentID: " <<agent_id_of_longest_stay_at_home_vilage;
    if (hmn_mgr.get_is_forest_goer(agent_id_of_longest_stay_at_home_vilage)) std::cout << "*";
    std::cout << " [home v: " <<  hmn_mgr.get_home_village(agent_id_of_longest_stay_at_home_vilage)<<"]";
    std::cout << ", at location: " <<vll_mgr.get_location_type(location_of_longest_stay_at_home)<<location_of_longest_stay_at_home;
    std::cout << ", for " <<nb_of_longest_stay_at_home_vilage <<" days \n";

    std::cout << "In Forest location, by moving agentID: " <<agent_id_of_longest_stay_at_forest;
    if (hmn_mgr.get_is_forest_goer(agent_id_of_longest_stay_at_forest)) std::cout << "*";
    std::cout << ", at location: " <<vll_mgr.get_location_type(location_of_longest_stay_at_forest)<<location_of_longest_stay_at_forest;
    std::cout << ", for " <<nb_of_longest_stay_at_forest <<" days \n";

    if (static_move_took_place) {
        std::cout << "In another village (static move took place), ";
        std::cout << "by moving agentID: " <<agent_id_of_longest_stay_at_another_village;
        if (hmn_mgr.get_is_forest_goer(agent_id_of_longest_stay_at_another_village)) std::cout << "*";
        std::cout << " [home v: " <<  hmn_mgr.get_home_village(agent_id_of_longest_stay_at_another_village)<<"]";
        std::cout << ", at location: " <<vll_mgr.get_location_type(location_of_longest_stay_at_another_village)<<location_of_longest_stay_at_another_village;
        std::cout << ", for " <<nb_of_longest_stay_at_another_village <<" days \n";
    }
    else 
        std::cout << "  No static move took place during simulation! \n";
    
    
    //REQUIRE(stay_home_between_mobilities == true);

    std::cout << "\n ==============================================================================================  \n";

    std::cout << "\n >>>>>>>>> Statistics related to effects of 'forest_precedence_over_static_mobility' parameter, set to: " << forest_precedence_over_static_mobility << "\n";
    std::cout << "The entire route to forest until return home move is counted as 1, numbers of days there :   \n";

    std::cout << "\nTotal number of villages: "  << total_villages <<"\n";
    std::cout << "Total number of forest POI: "  << total_forest_poi <<"\n";

    std::cout << "\nNumber of simulation days: "  << nb_of_simulated_days <<"\n";
    std::cout << "\nTotal number of agents (at t=0): "  << total_agents <<"\n";
    std::cout << "Total number of forest_goer agents (at t=0): "  << total_forest_goers <<"\n";
    std::cout << "\n";

    std::cout << "Total number of forest moves: "  << total_nb_of_forest_moves << "\n";
    std::cout << "                   , per day: "   << (total_nb_of_forest_moves/nb_of_simulated_days) << "\n";
    //std::cout << "           , per forest_goer: "   << (total_nb_of_forest_moves/total_forest_goers)  <<  ", per forest_goer/day: " << ((float)total_nb_of_forest_moves/(float)total_forest_goers)/(float)nb_of_simulated_days << "\n";
    //std::cout << "                , per agents: "   << ((float)total_nb_of_forest_moves/(float)total_agents)  <<  ", per agents/day: " << ((float)total_nb_of_forest_moves/(float)total_agents)/(float)nb_of_simulated_days<< "\n";
    
    std::cout << "\n";
    std::cout << "Total number of forest moves (including in trans routes): "  << total_nb_of_forest_moves_including_in_trans_routes << "\n";
    std::cout << "                                               , per day: "   << (total_nb_of_forest_moves_including_in_trans_routes/nb_of_simulated_days) << "\n";


    std::cout << "\n";
    std::cout << "Total number of stay days in forest: "  << total_days_spent_outside_home_in_forest_moves<< "\n";
    std::cout << "                          , per day: "  << (total_days_spent_outside_home_in_forest_moves/nb_of_simulated_days) << "\n";
    //std::cout << "                 , per forest_goers: "  << (total_days_spent_outside_home_in_forest_moves/total_forest_goers)  << ", per forest_goers/day: "   << ((float)total_days_spent_outside_home_in_forest_moves/(float)total_forest_goers)/(float)nb_of_simulated_days<< "\n";
    //std::cout << "                       , per agents: "  << ((float)total_days_spent_outside_home_in_forest_moves/(float)total_agents)  << ", per agents/day: "   << ((float)total_days_spent_outside_home_in_forest_moves/(float)total_agents)/(float)nb_of_simulated_days<< "\n";

    std::cout << "\n";
    std::cout << "Total number of static moves: "  << total_nb_of_static_moves << "\n";
    std::cout << "                   , per day: "   << (total_nb_of_static_moves/nb_of_simulated_days) << "\n";
    //std::cout << "          , per forest_goers: " << (total_nb_of_static_moves/total_forest_goers)  <<  ", per forest_goers/day: " << ((float)total_nb_of_static_moves/(float)total_forest_goers)/(float)nb_of_simulated_days<< "\n";
    //std::cout << "                , per agents: "   << ((float)total_nb_of_static_moves/(float)total_agents)  <<  ", per agents/day: " << ((float)total_nb_of_static_moves/(float)total_agents)/(float)nb_of_simulated_days<< "\n";

    std::cout << "\n";
    std::cout << "Total number of static moves by forest_goers: "  << total_nb_of_static_moves_by_forest_goers << "\n";
    std::cout << "                                   , per day: "   << (total_nb_of_static_moves_by_forest_goers/nb_of_simulated_days) << "\n";

    std::cout << "\n";
    std::cout << "Total number of stay days in other villages: "  << total_days_spent_outside_home_in_static_moves << "\n";
    std::cout << "                                  , per day: "   << (total_days_spent_outside_home_in_static_moves/nb_of_simulated_days) << "\n";
    //std::cout << "                         , per forest_goers: " << (total_days_spent_outside_home_in_static_moves/total_forest_goers)  << ", per forest_goers/day: " << ((float)total_days_spent_outside_home_in_static_moves/(float)total_forest_goers)/(float)nb_of_simulated_days<< "\n";
    //std::cout << "                               , per agents: "  << ((float)total_days_spent_outside_home_in_static_moves/(float)total_agents)  << ", per agents/day: "   << ((float)total_days_spent_outside_home_in_static_moves/(float)total_agents)/(float)nb_of_simulated_days<< "\n";


    int both_moves = total_nb_of_forest_moves + total_nb_of_static_moves;
    std::cout << "\n";
    std::cout << "Total number of moves (forest+static): "  << both_moves << "\n";
    std::cout << "                            , per day: "   << (both_moves/nb_of_simulated_days) << "\n";
    //std::cout << "          , per forest_goers: " << (both_moves/total_forest_goers)  <<  ", per forest_goers/day: " << ((float)both_moves/(float)total_forest_goers)/(float)nb_of_simulated_days<< "\n";
    //std::cout << "                , per agents: "   << ((float)both_moves/(float)total_agents)  <<  ", per agents/day: " << ((float)both_moves/(float)total_agents)/(float)nb_of_simulated_days<< "\n";

    int total_days_spent_outside_home_in_both_moves =total_days_spent_outside_home_in_forest_moves +total_days_spent_outside_home_in_static_moves;

    std::cout << "\n";
    std::cout << "Total number of stay days (forest+other villages): "  << total_days_spent_outside_home_in_both_moves << "\n";
    std::cout << "                                        , per day: "   << (total_days_spent_outside_home_in_both_moves/nb_of_simulated_days) << "\n";
    //std::cout << "                               , per forest_goers: " << (total_days_spent_outside_home_in_both_moves/total_forest_goers)  << ", per forest_goers/day: " << ((float)total_days_spent_outside_home_in_both_moves/(float)total_forest_goers)/(float)nb_of_simulated_days<< "\n";
    //std::cout << "                                     , per agents: "  << ((float)total_days_spent_outside_home_in_both_moves/(float)total_agents)  << ", per agents/day: "   << ((float)total_days_spent_outside_home_in_both_moves/(float)total_agents)/(float)nb_of_simulated_days<< "\n";

    std::cout << "\nTotal forest/static moves ratio: "  << (float)total_nb_of_forest_moves/(float)total_nb_of_static_moves <<"\n";
    //std::cout << "Total static/forest moves ratio: "  << (float)total_nb_of_static_moves/(float)total_nb_of_forest_moves <<"\n";

    std::cout << "Total forest [incl. trans.]/static moves ratio: "  << (float)total_nb_of_forest_moves_including_in_trans_routes/(float)total_nb_of_static_moves <<"\n";

}
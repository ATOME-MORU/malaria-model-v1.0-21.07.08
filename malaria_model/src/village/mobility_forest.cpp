#include <iostream>
#include <vector>
#include <cassert>
#include <algorithm>
#include <numeric>   // std::iota

#include <iomanip> //std::setw

#include "common/common.h"
#include "village/village.h"
#include "human/human.h"
#include "util/randomness.h"
#include "util/statistics.h"

#include "util/graph.h"

namespace village {

void VillageManager::build_forest_networks() {
    int idx_ent = 0;
    int idx_poi = 0;

    this->indices_of_ids.insert(this->indices_of_ids.end(), this->sum_num_villages, -1);
    
    for (int ll = 0; ll < this->sum_num_villages; ll++) {

        std::cout << ll << "," << this->location_type_of[ll] << "\n" ;
        if (this->location_type_of[ll] == common::LocationType::kForestEntry) {
            this->ids_of_forest_entry_points.push_back(ll);
            this->indices_of_ids[ll] = idx_ent;
            idx_ent++;
        }
        if (this->location_type_of[ll] == common::LocationType::kForestPOI) {
            this->ids_of_forest_pois.push_back(ll);
            this->indices_of_ids[ll] = idx_poi;
            idx_poi++;
        }

    }

    this->build_forest_entry_network();
    this->print_forest_entry_network();

    this->build_forest_poi_network();
    this->print_forest_poi_network();

    this->fg = this->build_forest_graph(
        this->distances,

        this->ids_of_forest_entry_points,
        this->ids_of_forest_pois,
        true
    );

    this->build_forest_route_network();

    this->build_forest_walk_network();
}

void VillageManager::build_forest_entry_network() {
    this->forest_entry_network = this->build_forest_entry_network_kernel_v1 (
        this->sum_num_villages,
        this->distances,
        this->ids_of_forest_entry_points
    );
    assert(static_cast<int>(this->forest_entry_network.size()) == this->sum_num_villages);
}

std::vector<std::vector<float>> VillageManager::build_forest_entry_network_kernel_v1(
        int num_locations,
        float** location_distances,
        const std::vector<int>& ids_of_forest_entry_points
    ) {

    std::vector<std::vector<float>> forest_entry_network;

    for (int ll = 0; ll < num_locations; ll++) {
        std::vector<float> ll_distance_to_entry_points;
        for (const int& ee : ids_of_forest_entry_points) {
            // ll_distance_to_entry_points.push_back(1.0/location_distances[ll][ee]);
            ll_distance_to_entry_points.push_back(
                VillageManager::func_dist_to_weight(location_distances[ll][ee])
            );
        }
        forest_entry_network.push_back(ll_distance_to_entry_points);
    }
    return forest_entry_network;
}

void VillageManager::print_forest_entry_network() const {
    std::cout << "VM::print_forest_entry_network(): ("
                << static_cast<int>(this->forest_entry_network.size())
                << " x "
                << static_cast<int>(this->forest_entry_network[0].size())
                << ")\n";
    for (const auto& rr:this->forest_entry_network) {
        std::cout << "\t";
        for (const auto& cc: rr) {
            std::cout << cc << ", ";
        }
        std::cout << "\n";
    }
}

void VillageManager::build_forest_poi_network() {
    this->forest_poi_network = this->build_forest_poi_network_kernel_v1(
        this->ids_of_forest_entry_points,
        this->ids_of_forest_pois,
        this->forest_entry_network
    );
    assert(this->forest_poi_network.size() == this->ids_of_forest_entry_points.size());
    assert(this->forest_poi_network[0].size() == this->ids_of_forest_pois.size());
}
std::vector<std::vector<float>> VillageManager::build_forest_poi_network_kernel_v1(
        const std::vector<int>& ids_of_forest_entry_points,
        const std::vector<int>& ids_of_forest_pois,
        const std::vector<std::vector<float>> forest_entry_network
    ) {

    std::vector<std::vector<float>> forest_poi_network;
    forest_poi_network.resize(ids_of_forest_entry_points.size());

    for (size_t iee = 0; iee < ids_of_forest_entry_points.size(); iee++) {
        for (const int& pp : ids_of_forest_pois) {
            forest_poi_network[iee].push_back(forest_entry_network[pp][iee]);
        }
    }

    return forest_poi_network;
}
void VillageManager::print_forest_poi_network() const {
    std::cout << "VM::print_forest_poi_network(): ("
                << static_cast<int>(this->forest_poi_network.size())
                << " x "
                << static_cast<int>(this->forest_poi_network[0].size())
                << ")\n";
    for (const auto& rr: this->forest_poi_network) {
        std::cout << "\t";
        for (const auto& cc: rr) {
            std::cout << cc << ", ";
        }
        std::cout << "\n";
    }
}

util::VillageGraph* VillageManager::build_forest_graph(
        const float* const* distances,

        const std::vector<int>& ids_of_forest_entry_points,
        const std::vector<int>& ids_of_forest_pois,
        bool verbose
    ) {

    std::vector<int> ids_of_locations;
    std::vector<common::LocationType> types_of_locations;
    for (int id: ids_of_forest_entry_points) {
        ids_of_locations.push_back(id);
        types_of_locations.push_back(common::LocationType::kForestEntry);
    }
    for (int id: ids_of_forest_pois) {
        ids_of_locations.push_back(id);
        types_of_locations.push_back(common::LocationType::kForestPOI);
    }

    size_t num_forest_locations = ids_of_forest_entry_points.size() + ids_of_forest_pois.size();
    assert(num_forest_locations == ids_of_locations.size());
    assert(num_forest_locations == types_of_locations.size());

    float ** forest_dist_matrix = new float *[num_forest_locations];
    for (size_t ll1 = 0; ll1 < num_forest_locations; ll1++) {
        forest_dist_matrix[ll1] = new float[num_forest_locations];
        std::fill_n(
            forest_dist_matrix[ll1],
            num_forest_locations,
            std::numeric_limits<float>::infinity()
        );

        for (size_t ll2 = 0; ll2 < num_forest_locations; ll2++) {
            if (types_of_locations[ll1] != common::LocationType::kForestEntry
                || types_of_locations[ll2] != common::LocationType::kForestEntry
                ) {
                forest_dist_matrix[ll1][ll2] = distances[ids_of_locations[ll1]][ids_of_locations[ll2]];
            }
        }
    }

    util::VillageGraph* fg = new util::VillageGraph(num_forest_locations, forest_dist_matrix, ids_of_locations);

    for (size_t ll = 0; ll < num_forest_locations; ll++) {
        fg->find_shortest_paths_route(0, ll, true);
    }

    if (verbose) {
        std::cout << "reducing graph to min spanning tree\n";
    }

    fg->mark_graph_with_min_spanning_tree();
    
    if (verbose) {
        fg->print_graph_to_file_dot("temp/temp_forest_routes");
    }
    
    fg->reduce_graph_to_min_spanning_tree();

    for (size_t ll = 0; ll < num_forest_locations; ll++) {
        fg->find_shortest_paths_route(0, ll, true);
    }

    for (size_t ll1 = 0; ll1 < num_forest_locations; ll1++) {
        if (forest_dist_matrix[ll1]) {
            delete[] forest_dist_matrix[ll1];
        }
    }
    if (forest_dist_matrix) {
        delete[] forest_dist_matrix;
    }

    return fg;
}

void VillageManager::build_forest_route_network() {
    // this->build_forest_route_network_tester();

    this->build_forest_route_network_kernel_v1(
        this->forest_route_network,
        this->forest_route_network_degree,
        this->forest_route_network_dist,

        // this->distances,
        this->fg,

        this->ids_of_forest_entry_points,
        this->ids_of_forest_pois,
        true
    );

    this->sort_routes_by_dist(
        true,
        this->forest_route_network,
        this->forest_route_network_degree,
        this->forest_route_network_dist
    );

    if (this->kMin_num_way_points > 1) {
        int rm = this->remove_short_routes (
                        this->kMin_num_way_points,
                        this->sum_num_villages,
                        this->forest_route_network,
                        this->forest_route_network_degree,
                        this->forest_route_network_dist
                    );
        std::cout << rm << " routes removed for including less than "
                    << this->kMin_num_way_points
                    << " way points (excl. entry points).\n";
    }

    this->build_forest_route_network_weight_v1(
        this->forest_route_network_dist,
        this->forest_route_network_weight
    );
    assert(this->forest_route_network_dist.size() == this->forest_route_network_weight.size());

    assert(this->forest_route_network.size() == this->ids_of_forest_entry_points.size());

    for (size_t ee = 0; ee < this->ids_of_forest_entry_points.size(); ee++) {
        std::cout << "Available routes at entry point" << ids_of_forest_entry_points[ee] << ":\n";
        assert(this->forest_route_network[ee].size() > 0);
        for (auto rr : this->forest_route_network[ee]) {
            for (auto ll : rr) {
                std::cout << ll << ", ";
            }
            std::cout << "\n";
        }
    }
    // std::cin.ignore();
}
void VillageManager::build_forest_route_network_kernel_v1(
        std::vector<std::vector<std::vector<int>>>& forest_route_network,
        std::vector<std::vector<std::vector<int>>>& forest_route_network_degree,
        std::vector<std::vector<float>>& forest_route_network_dist,

        util::VillageGraph* fg,

        const std::vector<int>& ids_of_forest_entry_points,
        const std::vector<int>& ids_of_forest_pois,
        bool verbose
    ) {

    std::vector<common::LocationType> types_of_locations(
                                            ids_of_forest_entry_points.size(),
                                            common::LocationType::kForestEntry
                                        );
    types_of_locations.insert(
        types_of_locations.end(),
        ids_of_forest_pois.size(),
        common::LocationType::kForestPOI
    );

    size_t num_forest_locations = types_of_locations.size();

    if (verbose) {
        std::cout << "building forest route network\n";
    }

    for (size_t ll1 = 0; ll1 < num_forest_locations; ll1++) {
        if (types_of_locations[ll1] == common::LocationType::kForestEntry) {
            
            std::vector<std::vector<int>> ll1_routes;
            std::vector<std::vector<int>> ll1_routes_degree;
            std::vector<float> ll1_dists;

            for (size_t ll2 = 0; ll2 < num_forest_locations; ll2++) {
                if (types_of_locations[ll2] == common::LocationType::kForestPOI) {

                    std::vector<int> route_ll1_to_ll2;
                    std::vector<int> route_degree_ll1_to_ll2;
                    
                    float dist_ll1_to_ll2 = fg->find_shortest_paths_route(
                                                ll1, ll2,
                                                route_ll1_to_ll2,
                                                route_degree_ll1_to_ll2,
                                                verbose
                                            );
                    
                    bool route_include_entry_point = false;
                    for (int ll_id : route_ll1_to_ll2) {
                        if (find(
                                ids_of_forest_entry_points.begin(),
                                ids_of_forest_entry_points.end(),
                                ll_id
                            ) != ids_of_forest_entry_points.end()) {
                            route_include_entry_point = true;
                        }
                    }
                    if (!route_include_entry_point) {
                        ll1_routes.push_back(route_ll1_to_ll2);
                        ll1_routes_degree.push_back(route_degree_ll1_to_ll2);
                        ll1_dists.push_back(dist_ll1_to_ll2);
                    }

                }
            }

            forest_route_network.push_back(ll1_routes);
            forest_route_network_degree.push_back(ll1_routes_degree);
            forest_route_network_dist.push_back(ll1_dists);
        }
    }
}
void VillageManager::build_forest_route_network_tester() {

    const int num_locations = 6;
    const float dist[num_locations][num_locations] = {
        {0.0, 1.0, 9.0, 1.0, 9.0, 1.0 },
        {0.0, 0.0, 2.0, 7.0, 1.0, 3.0 },
        {0.0, 0.0, 0.0, 3.0, 1.0, 9.0 },
        {0.0, 0.0, 0.0, 0.0, 2.0, 9.0 },
        {0.0, 0.0, 0.0, 0.0, 0.0, 1.0 },
        {0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
    };

    const std::vector<common::LocationType> types_of_locations{
        common::LocationType::kForestEntry,
        common::LocationType::kForestEntry,
        common::LocationType::kForestEntry,
        common::LocationType::kForestPOI,
        common::LocationType::kForestPOI,
        common::LocationType::kForestPOI
    };

    float** distances;

    distances = new float *[num_locations];
    for (int ll = 0; ll < num_locations; ll++) {
        distances[ll] = new float[num_locations];
        std::fill_n(distances[ll], num_locations, std::numeric_limits<float>::infinity());
    }


    std::cout << "Village distances:\n";
    for (int ll1 = 0; ll1 < num_locations; ll1++) {
        for (int ll2 = 0; ll2 < num_locations; ll2++) {
            distances[ll1][ll2] = dist[ll1][ll2];
            std::cout << "[" << std::setw(1) << ll1
                      << "][" << std::setw(1) << ll2
                      << "]" << std::setw(2) << std::setprecision(1)
                      << distances[ll1][ll2] << " ";
        }
        std::cout << std::endl;
    }

    std::cout << "Set distance between entry points to infinity:\n";
    for (int ll1 = 0; ll1 < num_locations; ll1++) {
        for (int ll2 = 0; ll2 < num_locations; ll2++) {
            if (types_of_locations[ll1] == common::LocationType::kForestEntry
                && types_of_locations[ll2] == common::LocationType::kForestEntry
                ) {
                distances[ll1][ll2] = std::numeric_limits<float>::infinity();
            }
            std::cout << "[" << std::setw(1) << ll1
                      << "][" << std::setw(1) << ll2
                      << "]" << std::setw(2) << std::setprecision(1)
                      << distances[ll1][ll2] << " ";
        }
        std::cout << std::endl;
    }


    std::vector<int> ids_of_vll(num_locations);
    std::iota(std::begin(ids_of_vll), std::end(ids_of_vll), 10); // 10, 11, 12, ..., 10+num_vll-1 

    std::vector<int> ids_of_entry_points;
    for (int ll = 0; ll < num_locations; ll++ ) {
        if (types_of_locations[ll] == common::LocationType::kForestEntry) {
            ids_of_entry_points.push_back(ids_of_vll[ll]);
        }
    }

    util::VillageGraph vg(num_locations, distances, ids_of_vll);

    for (int ll = 0; ll < num_locations; ll++) {
        vg.find_shortest_paths_route(2, ll, true);
    }


    std::cout << "reducing graph to min spanning tree\n";
    vg.mark_graph_with_min_spanning_tree();
    vg.print_graph_to_file_dot("temp/temp_forest_routes");
    vg.reduce_graph_to_min_spanning_tree();
    

    for (int ll = 0; ll < num_locations; ll++) {
        vg.find_shortest_paths_route(2, ll, true);
    }

    std::cout << "building forest route network\n";
    std::vector<std::vector<std::vector<int>>> forest_route_network;

    for (int ll1 = 0; ll1 < num_locations; ll1++) {
        if (types_of_locations[ll1] == common::LocationType::kForestEntry) {

            std::vector<std::vector<int>> ll1_routes;

            for (int ll2 = 0; ll2 < num_locations; ll2++) {
                if (types_of_locations[ll2] == common::LocationType::kForestPOI) {
                    std::vector<int> route_ll1_to_ll2 = vg.find_shortest_paths_route(ll1, ll2, true);

                    bool route_include_entry_point = false;
                    for (int ll_id : route_ll1_to_ll2) {
                        if (find(
                                ids_of_entry_points.begin(),
                                ids_of_entry_points.end(),
                                ll_id)
                            != ids_of_entry_points.end()) {
                            route_include_entry_point = true;
                        }
                    }
                    if (!route_include_entry_point) {
                        ll1_routes.push_back(route_ll1_to_ll2);
                    }
                }
            }
            forest_route_network.push_back(ll1_routes);
        }
    }

    for (size_t ll = 0; ll < forest_route_network.size(); ll++ ) {
        std::cout << "Available routes at location " << ids_of_entry_points[ll] << ":\n";
        for (auto rr : forest_route_network[ll]) {
            for (auto ll : rr) {
                std::cout << ll << ", ";
            }
            std::cout << "\n";
        }
    }

    std::cin.ignore();

}

void VillageManager::sort_routes_by_dist(
        bool ascend,
        std::vector<std::vector<std::vector<int>>>& forest_route_network,
        std::vector<std::vector<std::vector<int>>>& forest_route_network_degree,
        std::vector<std::vector<float>>& forest_route_network_dist
    ) {
    for (size_t ee = 0; ee < forest_route_network.size(); ee++) {
        std::vector<size_t> new_order;
        if (ascend) {
            new_order = util::sort_by_value_and_return_indices_ascend(forest_route_network_dist[ee]);
        } else {
            new_order = util::sort_by_value_and_return_indices_descend(forest_route_network_dist[ee]);
        }
        std::vector<std::vector<int>> ee_new_route_list;
        std::vector<std::vector<int>> ee_new_degree_list;
        std::vector<float> ee_new_dist_list;
        for (size_t& ii : new_order) {
            ee_new_route_list.push_back(forest_route_network[ee][ii]);
            ee_new_degree_list.push_back(forest_route_network_degree[ee][ii]);
            ee_new_dist_list.push_back(forest_route_network_dist[ee][ii]);
        }
        forest_route_network[ee] = ee_new_route_list;
        forest_route_network_degree[ee] = ee_new_degree_list;
        forest_route_network_dist[ee] = ee_new_dist_list;
    }
}

int VillageManager::remove_short_routes(
        int min_num_way_points,
        int num_locations,
        std::vector<std::vector<std::vector<int>>>& forest_route_network,
        std::vector<std::vector<std::vector<int>>>& forest_route_network_degree,
        std::vector<std::vector<float>>& forest_route_network_dist
    ) {
    int total_removed = 0;
    std::vector<int> poi_app_count(num_locations,0); // poi appearance count
    for (const auto& ee_routes_list : forest_route_network) {
        for (const auto& route : ee_routes_list) {
            for (const int& id: route) {
                assert(id < num_locations);
                poi_app_count[id]++;
            }
        }
    }
    for (size_t ee = 0; ee < forest_route_network.size(); ee++) {
        for (size_t rr = 0; rr < forest_route_network[ee].size(); rr++) {
            if (forest_route_network[ee].size() <= 1) {break;}
            if (static_cast<int>(forest_route_network[ee][rr].size()) < min_num_way_points) {
                bool all_way_points_has_alternative_route = true;
                // if route can be removed without excluding any poi from the network
                for (const int& ll : forest_route_network[ee][rr]) {
                    if (poi_app_count[ll] <= 1) {
                        all_way_points_has_alternative_route = false;
                    }
                }
                if (all_way_points_has_alternative_route) {
                    std::cout << "Removing route: ";
                    for (const int& ll : forest_route_network[ee][rr]) {
                        std::cout << ll << ", ";
                        poi_app_count[ll]--;
                    }
                    std::cout << "\n";

                    forest_route_network[ee].erase(forest_route_network[ee].begin()+rr);
                    forest_route_network_degree[ee].erase(forest_route_network_degree[ee].begin()+rr);
                    forest_route_network_dist[ee].erase(forest_route_network_dist[ee].begin()+rr);

                    rr--;
                    total_removed++;
                }
            }
        }
    }
    return total_removed;
}

void VillageManager::build_forest_route_network_weight_v1(
        const std::vector<std::vector<float>>& forest_route_network_dist,
        std::vector<std::vector<float>>& forest_route_network_weight
    ) {
    forest_route_network_weight.clear();
    for (auto ee_routes_dist: forest_route_network_dist) {
        std::vector<float> ee_weight;
        for (float rr_dist: ee_routes_dist) {
            // ee_weight.push_back(1.0/rr_dist);
            ee_weight.push_back(
                VillageManager::func_dist_to_weight(rr_dist)
            );
        }
        forest_route_network_weight.push_back(ee_weight);
    }
}

void VillageManager::build_forest_walk_network() {

    bool kVerbose = true;
    
    this->forest_walk_network.clear();
    this->forest_walk_network_degree.clear();

    if (this->kFollow_routes) {
        this->build_forest_walk_network_kernel_v0(
            this->forest_walk_network,
            this->forest_walk_network_degree,
            this->fg,
            this->ids_of_forest_entry_points,
            this->ids_of_forest_pois,
            kVerbose
        );
    } else {
        this->build_forest_walk_network_kernel_v1(
            this->forest_walk_network,
            this->forest_walk_network_degree,
            this->distances,
            this->ids_of_forest_pois,
            kVerbose
        );
    }
}

void VillageManager::build_forest_walk_network_kernel_v0 (
        std::vector<std::vector<float>>& forest_walk_network,
        std::vector<float>& forest_walk_network_degree,

        util::VillageGraph* fg,

        const std::vector<int>& ids_of_forest_entry_points,
        const std::vector<int>& ids_of_forest_pois,
        bool verbose
    ) {

    if (verbose) {
        std::cout << "building forest walk network (follow routes)\n";
    }

    for (size_t ii = 0; ii < ids_of_forest_pois.size(); ii++) {
        std::vector<float> ii_row;
        size_t ll1 = ii + ids_of_forest_entry_points.size();

        int dg = fg->get_degree_of_village(ll1);
        // int dg = 0;
        for (size_t jj = 0; jj < ids_of_forest_pois.size(); jj++) {
            size_t ll2 = jj + ids_of_forest_entry_points.size();

            float d = fg->get_edge_weight(ll1, ll2);
            if (ll1 == ll2) {d = 0.0;}
            if (d > 0.0) {
                // ii_row.push_back(1.0/d);
                ii_row.push_back(
                    VillageManager::func_dist_to_weight(d)
                );
                // dg++;
                if (verbose) {
                    std::cout << ids_of_forest_pois[ii] << " to "
                              << ids_of_forest_pois[jj] << ":"
                              << d << "\n";
                }
            } else {
                ii_row.push_back(0.0);
            }
        }
        if (verbose) {
            std::cout << ids_of_forest_pois[ii] << " degree: "
                      << dg << "\n";
        }
        forest_walk_network.push_back(ii_row);
        forest_walk_network_degree.push_back(dg);
    }

}

void VillageManager::build_forest_walk_network_kernel_v1 (
        std::vector<std::vector<float>>& forest_walk_network,
        std::vector<float>& forest_walk_network_degree,
        float** location_distances,
        const std::vector<int>& ids_of_forest_pois,
        bool verbose
    ) {

    if (verbose) {
        std::cout << "building forest walk network (straightline distance)\n";
    }

    for (size_t ii = 0; ii < ids_of_forest_pois.size(); ii++) {
        std::vector<float> ii_row;
        for (size_t jj = 0; jj < ids_of_forest_pois.size(); jj++) {
            ii_row.push_back( VillageManager::func_dist_to_weight(
                            location_distances[ids_of_forest_pois[ii]][ids_of_forest_pois[jj]]
                        ));
            if (verbose) {
                std::cout << ids_of_forest_pois[ii] << " to "
                          << ids_of_forest_pois[jj] << ":"
                          << ii_row.back() << "\n";
            }
        }
        forest_walk_network.push_back(ii_row);
        forest_walk_network_degree.push_back(ids_of_forest_pois.size()-1);
    }

}

int VillageManager::step_population_movement_forest() {
    if (this->kMobility_forest_enabled) {
        // return this->step_population_movement_forest_kernel_v1();
        switch(this->kFunc_mobility_forest_version) {
            case 0 : 
                return this->step_population_movement_forest_kernel_v0();
            case 1 :
                return this->step_population_movement_forest_kernel_v1();
            case 2 :
                return this->step_population_movement_forest_kernel_v2();
            default :
                std::cout << "VM: Unknown kFunc_mobility_forest_version = "
                    << this->kFunc_mobility_forest_version << "\n";
                assert(false);
        }
    }
    return 0;
}
int VillageManager::proc_forest_entry_point(
        int ll,
        int value_for_erase
    ) {

    int num_moves = 0;

    if (this->kStay_at_entry_points || this->kReturn_stay_at_entry_points) {
        for (int& hh: this->at_village_register[ll]) {
            if (!this->human_already_moved_in_this_step[hh]) {
                
                if (this->v_hmn_mgr->journey_is_empty(hh)) { // exiting forest
                    
                    assert(this->kReturn_stay_at_entry_points);
                    assert(this->last_entry_point_visited[hh] == ll);

                    this->proc_return_human_to_home(
                            hh,
                            value_for_erase
                    );
                    num_moves++;
                
                } else { // entering forest

                    assert(this->kStay_at_entry_points);

                    this->last_entry_point_visited[hh] = ll;
                    if (this->forest_entry_turn_back) {
                        this->proc_return_human_to_home(
                                hh,
                                value_for_erase
                        );
                    } else {
                        this->proc_move_human_to_location(
                            hh,
                            this->v_hmn_mgr->journey_pop_next_location(hh),
                            value_for_erase
                        );
                    }
                    num_moves++;

                }

            }
        }
    } else {
        if (this->at_village_register[ll].size() != 0) {
            std::cout << "VM::proc_forest_entry_point:\n"
                        << "\tLocation " << ll <<" (entry point) is occupied by "
                        << this->at_village_register[ll].size()
                        << " when staying at entry points is not allowed.\n";
            assert(false);
        }
    }
    return num_moves;
}
void VillageManager::proc_return_human_to_home_forest(
        int& human_id,
        int value_for_erase
    ) {
    if (this->kReturn_stay_at_entry_points) {
        this->proc_move_human_to_location(
            human_id,
            this->last_entry_point_visited[human_id],
            value_for_erase
        );
    } else {
        this->proc_return_human_to_home(
            human_id,
            value_for_erase
        );
    }
}
int VillageManager::step_population_movement_forest_kernel_v0() {
    const int kValue_for_erase = -1;
    int num_moves = 0;
    for (int vv = 0; vv < this->sum_num_villages; vv++) {
        int vv_num_moves = 0;
        switch(this->get_location_type(vv)){

            case common::LocationType::kVillage:
                for (int& hh: this->at_village_register[vv]) {

                    if (this->forest_banned
                        && util::get_rand_uniform() < this->kForest_ban_compliance
                    ) {
                        continue;
                    }

                    if (this->v_hmn_mgr->get_is_forest_goer(hh)
                        && this->v_hmn_mgr->get_if_at_home(hh)
                        && !this->human_already_moved_in_this_step[hh]
                        && util::get_rand_uniform() < this->kForest_set_off_probability 
                    ) {
                        // select entry point
                        int entry_point_idx = util::get_rand_discrete(this->forest_entry_network[vv]);
                        // select poi
                        int poi_idx = util::get_rand_discrete(this->forest_poi_network[entry_point_idx]);

                        int move_to = 0;
                        
                        if (this->kStay_at_entry_points) {
                            move_to = this->ids_of_forest_entry_points[entry_point_idx];
                            this->v_hmn_mgr->journey_add_location(
                                                hh,
                                                this->ids_of_forest_pois[poi_idx]
                                            );
                        } else {
                            move_to = this->ids_of_forest_pois[poi_idx];
                        }

                        // int poi = this->ids_of_forest_pois[poi_idx];
                        // this->at_village_register[poi].push_back(hh);
                        // this->v_hmn_mgr->set_human_at_village(hh, poi);
                        // this->human_already_moved_in_this_step[hh] = true;
                        // hh = kValue_for_erase;
                        this->proc_move_human_to_location(
                            hh,
                            move_to,
                            kValue_for_erase
                        );
                        vv_num_moves++;
                    }
                }
                break;

            case common::LocationType::kForestPOI:
                for (int& hh: this->at_village_register[vv]) {
                    if (!this->human_already_moved_in_this_step[hh]
                        && util::get_rand_uniform() < this->kForest_return_home_probability
                    ) {
                        // this->at_village_register[this->v_hmn_mgr->get_home_village(hh)].push_back(hh);
                        // this->v_hmn_mgr->set_human_return_home_village(hh);
                        // this->human_already_moved_in_this_step[hh] = true;
                        // hh = kValue_for_erase;
                        this->proc_return_human_to_home_forest(
                            hh, kValue_for_erase
                        );
                        vv_num_moves++;
                    }
                }
                break;

            case common::LocationType::kForestEntry:

                vv_num_moves += this->proc_forest_entry_point(vv, kValue_for_erase);
                break;

            default:
                std::cout << "VM:step_population_movement_forest_kernel_v0():\n"
                            << "\tUnkonwn village location type (vv=" << vv
                            << ")\n";
                assert(false);
        }
        if (vv_num_moves > 0) {
            this->at_village_register[vv].erase(
                std::remove(
                    this->at_village_register[vv].begin(),
                    this->at_village_register[vv].end(),
                    kValue_for_erase
                ),
                this->at_village_register[vv].end()
            );
        }
        num_moves += vv_num_moves;
    }
    return num_moves;
}
int VillageManager::step_population_movement_forest_kernel_v1() {
    const int kValue_for_erase = -1;
    
    int num_journey_started = 0;
    int num_moves = 0;
    for (int ll = 0; ll < this->sum_num_villages; ll++) {

        int ll_num_moves = 0;
        switch(this->get_location_type(ll)) {

            case common::LocationType::kVillage:
                for (int& hh : this->at_village_register[ll]) {

                    if (this->forest_banned
                        && util::get_rand_uniform() < this->kForest_ban_compliance
                    ) {
                        continue;
                    }

                    if ( this->v_hmn_mgr->get_is_forest_goer(hh)
                        && this->v_hmn_mgr->get_if_at_home(hh)
                        && !this->human_already_moved_in_this_step[hh]
                        && util::get_rand_uniform() < this->kForest_set_off_probability 
                        ) {

                        // Start a new journey

                        assert(this->v_hmn_mgr->journey_is_empty(hh));

                        // 1. select entry point
                        int entry_point_idx = util::get_rand_discrete(this->forest_entry_network[ll]);
                        // 2. select route
                        int route_idx = util::get_rand_discrete(this->forest_route_network_weight[entry_point_idx]);
                        std::vector<int> route = this->forest_route_network[entry_point_idx][route_idx];

                        // 3. register journey route
                        for (auto ill = route.rbegin(); ill != route.rend(); ++ill) {
                            this->v_hmn_mgr->journey_add_location(hh, *ill);
                        }
                        if (this->kStay_at_entry_points) {
                            this->v_hmn_mgr->journey_add_location(hh, this->ids_of_forest_entry_points[entry_point_idx]);
                        }
                        num_journey_started++;

                        // 4. move to first location
                        // int first_location = this->v_hmn_mgr->journey_pop_next_location();
                        // this->at_village_register[first_location].push_back(hh);
                        // this->v_hmn_mgr->set_human_at_village(hh, first_location);
                        // this->human_already_moved_in_this_step[hh] = true;
                        
                        // hh = kValue_for_erase;
                        this->proc_move_human_to_location(
                            hh,
                            this->v_hmn_mgr->journey_pop_next_location(hh),
                            kValue_for_erase
                        );

                        ll_num_moves++;

                    }
                }
                break;

            case common::LocationType::kForestPOI:
                for (int& hh : this->at_village_register[ll]) {
                    if (!this->human_already_moved_in_this_step[hh]) {

                        if (this->kReturn_early) { // allows returning before end of route
                            if (util::get_rand_uniform() < this->kForest_return_home_probability) {
                                this->v_hmn_mgr->journey_clear(hh);
                                this->proc_return_human_to_home_forest(
                                    hh, kValue_for_erase
                                );
                                ll_num_moves++;
                            } else {
                                if (!this->v_hmn_mgr->journey_is_empty(hh)) {
                                    this->proc_move_human_to_location(
                                        hh,
                                        this->v_hmn_mgr->journey_pop_next_location(hh),
                                        kValue_for_erase
                                    );
                                    ll_num_moves++;
                                }
                            }
                        } else {
                            if (this->v_hmn_mgr->journey_is_empty(hh)) {
                                if (util::get_rand_uniform() < this->kForest_return_home_probability) {
                                    this->proc_return_human_to_home_forest(
                                        hh, kValue_for_erase
                                    );
                                    ll_num_moves++;
                                }
                            } else {
                                this->proc_move_human_to_location(
                                    hh,
                                    this->v_hmn_mgr->journey_pop_next_location(hh),
                                    kValue_for_erase
                                );
                                ll_num_moves++;
                            }
                        }

                    }
                }
                break;

            case common::LocationType::kForestEntry:

                ll_num_moves += this->proc_forest_entry_point(ll, kValue_for_erase);
                break;

            default:
                std::cout << "VM:step_population_movement_forest_kernel_v1():\n"
                            << "\tUnkonwn location type (ll=" << ll
                            << ")\n";
                assert(false);

        }

        if (ll_num_moves > 0) {
            this->at_village_register[ll].erase(
                std::remove(
                    this->at_village_register[ll].begin(),
                    this->at_village_register[ll].end(),
                    kValue_for_erase
                ),
                this->at_village_register[ll].end()
            );
        }

        num_moves += ll_num_moves;

    }
    return num_moves;
}

int VillageManager::step_population_movement_forest_kernel_v2() {
    const int kValue_for_erase = -1;

    int num_moves = 0;
    for (int ll = 0; ll < this->sum_num_villages; ll++) {

        int ll_num_moves = 0;
        switch (this->get_location_type(ll)) {

            case common::LocationType::kVillage:
            // TODO: possibility to combine the 3 types into one function

                for (int& hh: this->at_village_register[ll]) {

                    if (this->forest_banned
                        && util::get_rand_uniform() < this->kForest_ban_compliance
                    ) {
                        continue;
                    }

                    if (this->v_hmn_mgr->get_is_forest_goer(hh)
                        && this->v_hmn_mgr->get_if_at_home(hh)
                        && !this->human_already_moved_in_this_step[hh]
                        && util::get_rand_uniform() < this->kForest_set_off_probability 
                    ) {

                        assert(this->v_hmn_mgr->journey_is_empty(hh));

                        // select entry point
                        int entry_point_idx = util::get_rand_discrete(this->forest_entry_network[ll]);
                        // select poi
                        int poi_idx = util::get_rand_discrete(this->forest_poi_network[entry_point_idx]);

                        int move_to = 0;
                        
                        if (this->kStay_at_entry_points) {
                            move_to = this->ids_of_forest_entry_points[entry_point_idx];
                            this->v_hmn_mgr->journey_add_location(
                                                hh,
                                                this->ids_of_forest_pois[poi_idx]
                                            );
                        } else {
                            move_to = this->ids_of_forest_pois[poi_idx];
                        }

                        // int poi = this->ids_of_forest_pois[poi_idx];
                        // this->at_village_register[poi].push_back(hh);
                        // this->v_hmn_mgr->set_human_at_village(hh, poi);
                        // this->human_already_moved_in_this_step[hh] = true;
                        // hh = kValue_for_erase;
                        this->proc_move_human_to_location(
                            hh,
                            move_to,
                            kValue_for_erase
                        );
                        ll_num_moves++;
                    }
                }

                break;

            case common::LocationType::kForestPOI:
                for (int& hh: this->at_village_register[ll]) {
                    if (!this->human_already_moved_in_this_step[hh]) {

                        if (util::get_rand_uniform() < this->kForest_return_home_probability) {
                            this->proc_return_human_to_home_forest(
                                hh, kValue_for_erase
                            );
                            ll_num_moves++;
                        } else {

                            int idx_ll = this->indices_of_ids[ll];
                            if (this->forest_walk_network_degree[idx_ll] > 1) {
                                int poi_idx = util::get_rand_discrete(this->forest_walk_network[idx_ll]);
                                int move_to = this->ids_of_forest_pois[poi_idx];
                                assert(move_to >= 0);
                                assert(move_to < this->sum_num_villages);
                                this->proc_move_human_to_location(
                                    hh,
                                    move_to,
                                    kValue_for_erase
                                );
                                ll_num_moves++;
                            } else {
                                // no move at a leaf, degree = 0 or 1
                            }

                        }

                    }

                }

                break;

            case common::LocationType::kForestEntry:

                ll_num_moves += this->proc_forest_entry_point(ll, kValue_for_erase);
                break;

            default:
                std::cout << "VM:step_population_movement_forest_kernel_v2:\n"
                            << "\tUnkonwn location type (ll=" << ll
                            << ")\n";
                assert(false);

        }

        if (ll_num_moves > 0) {
            this->at_village_register[ll].erase(
                std::remove(
                    this->at_village_register[ll].begin(),
                    this->at_village_register[ll].end(),
                    kValue_for_erase
                ),
                this->at_village_register[ll].end()
            );
        }

        num_moves += ll_num_moves;

    }

    return num_moves;

}

}
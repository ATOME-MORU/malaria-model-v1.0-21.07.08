#include <vector>
#include <iostream>
#include <cassert>
#include <cmath>
#include <algorithm> // std::min

#include "common/common.h"
#include "util/randomness.h"


#include "village/mosquito_dispersion.h"
#include "village/village.h"

namespace village {

MosquitoDP::MosquitoDP(
        bool enabled,
        float radius,
        float strength,
        bool individual_based,

        float retaining_strength_human_presence_each,
        float retaining_strength_human_presence_max,
        float retaining_strength_spawn_point_each,
        float retaining_strength_spawn_point_max,

        bool pit_enabled,
        float pit_lost_prob,
        float pit_found_prob,

        const float* const* distances,
        const float* mosquito_spawning_site_size_of,
        const std::vector<int>& ids_of_dispersion_locations
    ) : 
        kEnabled(enabled),
        kRadius(radius),
        kStrength(strength),
        kIndividualBased(individual_based),

        kRetaining_strength_human_presence_each(retaining_strength_human_presence_each),
        kRetaining_strength_human_presence_max(retaining_strength_human_presence_max),
        kRetaining_strength_spawn_point_each(retaining_strength_spawn_point_each),
        kRetaining_strength_spawn_point_max(retaining_strength_spawn_point_max),

        kPit_enabled(pit_enabled),
        kPit_lost_prob(pit_lost_prob),
        kPit_found_prob(pit_found_prob),

        ids_of_swamps(ids_of_dispersion_locations)
    {

    assert(kPit_lost_prob <= 1.0 && kPit_lost_prob >= 0.0);
    assert(kPit_found_prob <= 1.0 && kPit_found_prob >= 0.0);

    for (const int& ss1 : this->ids_of_swamps) {

        // Dispersion Network
        std::vector<float> ss1_weight_to_swamps;
        float weight_sum = 0.0;
        for (const int& ss2 : this->ids_of_swamps) {
            float ww = 0.0;
            if (distances[ss1][ss2] <= kRadius) {
                if (ss1 == ss2 || distances[ss1][ss2] == 0) {
                    ww = 0.0;
                } else {
                    // ww = 1.0/distances[ss1][ss2];
                    ww = VillageManager::func_dist_to_weight(distances[ss1][ss2]);
                }
            }
            ss1_weight_to_swamps.push_back(ww);
            weight_sum += ww;
        }
        if (weight_sum != 0.0) {
            for (float& ww : ss1_weight_to_swamps) {
                ww = ww/weight_sum;
            }
        }
        this->dispersion_network.push_back(ss1_weight_to_swamps);



        // Retaining Strength by Spawn Point
        float rs_sp = std::min(
                this->kRetaining_strength_spawn_point_max,
                mosquito_spawning_site_size_of[ss1] * this->kRetaining_strength_spawn_point_each
              );
        assert(rs_sp >= 0.0);
        this->retaining_strength_spawn_point.push_back(rs_sp);

    }
    
}

void MosquitoDP::simple_dispersion(
        int** mosquito_composition_of,
        const int* population_of,
        bool if_infectious
    ) {
    // const size_t kNumSwamps = this->ids_of_swamps.size();
    
    // std::vector<int> dispersion_vector(kNumSwamps, 0);

    if (!this->kEnabled) { return; }

    int pt_type = static_cast<int>(common::ParasiteType::kR0);

    int lost_to_pit = 0;
    
    for (size_t ss1_idx = 0; ss1_idx < this->ids_of_swamps.size(); ss1_idx++) {
        int vv1 = this->ids_of_swamps[ss1_idx];


        float dispersion_strength = this->kStrength;
        float retaining_strength_human_presence = std::min(
                    this->kRetaining_strength_human_presence_max,
                    population_of[vv1] * this->kRetaining_strength_human_presence_each
                );
        dispersion_strength -= retaining_strength_human_presence;
        dispersion_strength -= this->retaining_strength_spawn_point[ss1_idx];


        if (dispersion_strength <= 0.0) {
            continue;
        }


        int num_mosq_dp = std::floor(mosquito_composition_of[vv1][pt_type]*this->kStrength);

        if (this->kPit_enabled) {
            int lost_to_pit_vv = std::floor(num_mosq_dp*this->kPit_lost_prob);
            lost_to_pit += lost_to_pit_vv;
            num_mosq_dp -= lost_to_pit_vv;

            mosquito_composition_of[vv1][pt_type] -= lost_to_pit_vv;

        }

        if (num_mosq_dp > 0) {

            int total_dispersed = 0;

            if (this->kIndividualBased) {
                for (int mm = 0; mm < num_mosq_dp; mm++) {
                    int target_idx = util::get_rand_discrete(this->dispersion_network[ss1_idx]);
                    int target_vv = this->ids_of_swamps[target_idx];
                    // dispersion_vector[target_idx]++;
                    mosquito_composition_of[target_vv][pt_type]++;
                    total_dispersed++;
                }
            } else {
                for (size_t target_idx = 0; target_idx < this->dispersion_network[ss1_idx].size(); target_idx++) {
                    int dd = std::floor(num_mosq_dp * this->dispersion_network[ss1_idx][target_idx]);
                    int target_vv = this->ids_of_swamps[target_idx];
                    // dispersion_vector[target_idx] += dd;
                    mosquito_composition_of[target_vv][pt_type] += dd;
                    total_dispersed += dd;
                }
            }

            assert(total_dispersed <= num_mosq_dp);
            mosquito_composition_of[vv1][pt_type] -= total_dispersed;

        }

    }

    if (this->kPit_enabled) {

        int found_from_pit = 0;
        if (if_infectious) {
             found_from_pit = std::floor(this->pit_num_infectious * this->kPit_found_prob);
        } else {
             found_from_pit = std::floor(this->pit_num_infected * this->kPit_found_prob);
        }
        for (int mm = 0; mm < found_from_pit; mm++) {
            size_t target_idx = util::get_rand_uniform() * this->ids_of_swamps.size();
            int target_vv = this->ids_of_swamps[target_idx];
            mosquito_composition_of[target_vv][pt_type] += 1;
        }

        // Remove & Add to pit
        if (if_infectious) {
            this->pit_num_infectious += lost_to_pit;
            this->pit_num_infectious -= found_from_pit;
        } else {
            this->pit_num_infected += lost_to_pit;
            this->pit_num_infected -= found_from_pit;
        }

    }
}

void MosquitoDP::mosquito_incubation_step(
        float incubation_rate
    ) {

    if (!this->kEnabled || !this->kPit_enabled) { return; }

    int incubated = 0;
    if (this->pit_num_infected > 0) {
        for (int mm = 0; mm < this->pit_num_infected; mm++) {
            if (util::get_rand_uniform() < incubation_rate ) {
                incubated++;
            }
        }
        this->pit_num_infected -= incubated;
        this->pit_num_infectious += incubated;
    }
    this->log_step_mosq_pit_incubation = incubated;
}

void MosquitoDP::mosquito_survival_step(
        float death_rate_infected,
        float death_rate_infectious
    ) {
    
    if (!this->kEnabled || !this->kPit_enabled) { return; }

    this->log_step_mosq_pit_death_infected = this->mosquito_survival_step_kernel(
                                                death_rate_infected,
                                                this->pit_num_infected
                                            );
    this->log_step_mosq_pit_death_infectious = this->mosquito_survival_step_kernel(
                                                death_rate_infectious,
                                                this->pit_num_infectious
                                            );
}

int MosquitoDP::mosquito_survival_step_kernel(
        float death_rate,
        int& num_mosq
    ) {

    int killed = 0;
    if (num_mosq > 0) {
        // if (num_mosq * death_rate >= 1.0) {
        //     killed = std::round(num_mosq * death_rate);
        // } else {
            for (int mm = 0; mm < num_mosq; mm++) {
                if (util::get_rand_uniform() < death_rate) {
                    killed++;
                }
            }
        // }
        num_mosq -= killed;
    }

    return killed;
}

}
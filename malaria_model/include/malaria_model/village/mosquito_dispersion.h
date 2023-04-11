#ifndef MOSQUITO_DP_H
#define MOSQUITO_DP_H

namespace village {

class MosquitoDP {

    const bool kEnabled;
    const float kRadius;
    const float kStrength;

    const bool kIndividualBased;

    const float kRetaining_strength_human_presence_each;
    const float kRetaining_strength_human_presence_max;
    const float kRetaining_strength_spawn_point_each;
    const float kRetaining_strength_spawn_point_max;

    // Pit
    const bool kPit_enabled;
    const float kPit_lost_prob;
    const float kPit_found_prob;

    int pit_num_infected = 0;
    int pit_num_infectious = 0;

    // const int kVersion;

    // const int kNumLocations;
    const std::vector<int> ids_of_swamps;

    std::vector<float> retaining_strength_spawn_point;
    std::vector<std::vector<float>> dispersion_network;

public:

    int log_step_mosq_pit_death_infected = 0;
    int log_step_mosq_pit_death_infectious = 0;
    int log_step_mosq_pit_incubation = 0;

    MosquitoDP(
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
    );

    // MosquitoDP(
    //     const std::vector<std::vector<float>> distances,
    //     const std::vector<int>& ids_of_dispersion_locations
    // );

    void simple_dispersion(
        int** mosquito_composition_of,
        const int* population_of,
        bool if_infectious
    );

    void mosquito_incubation_step(
        float incubation_rate
    );

    void mosquito_survival_step(
        float death_rate_infected,
        float death_rate_infectious
    );

    ///////////////////////////////////////////
    // 
    inline int get_num_infected() const {
        return this->pit_num_infected;
    }
    inline int get_num_infectious() const {
        return this->pit_num_infectious;
    }

private:
    int mosquito_survival_step_kernel(
        float death_rate,
        int& num_mosq
    );

};


}

#endif
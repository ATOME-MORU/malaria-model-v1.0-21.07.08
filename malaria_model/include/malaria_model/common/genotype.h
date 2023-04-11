#ifndef GENOTYPE_H
#define GENOTYPE_H

#include <unordered_map>

namespace common {

typedef uint8_t genotype_t;
enum class ParasiteType: genotype_t { 
    kKNY__C1 = 0,
    kKNY__C2 = 1,
    kKNY__Y1 = 2,
    kKNY__Y2 = 3,
    kKNYNYC1 = 4,
    kKNYNYC2 = 5,
    kKNYNYY1 = 6,
    kKNYNYY2 = 7,
    kKNF__C1 = 8,
    kKNF__C2 = 9,
    kKNF__Y1 = 10,
    kKNF__Y2 = 11,
    kKNFNFC1 = 12,
    kKNFNFC2 = 13,
    kKNFNFY1 = 14,
    kKNFNFY2 = 15,
    kKYY__C1 = 16,
    kKYY__C2 = 17,
    kKYY__Y1 = 18,
    kKYY__Y2 = 19,
    kKYYYYC1 = 20,
    kKYYYYC2 = 21,
    kKYYYYY1 = 22,
    kKYYYYY2 = 23,
    kKYF__C1 = 24,
    kKYF__C2 = 25,
    kKYF__Y1 = 26,
    kKYF__Y2 = 27,
    kKYFYFC1 = 28,
    kKYFYFC2 = 29,
    kKYFYFY1 = 30,
    kKYFYFY2 = 31,
    kTNY__C1 = 32,
    kTNY__C2 = 33,
    kTNY__Y1 = 34,
    kTNY__Y2 = 35,
    kTNYNYC1 = 36,
    kTNYNYC2 = 37,
    kTNYNYY1 = 38,
    kTNYNYY2 = 39,
    kTNF__C1 = 40,
    kTNF__C2 = 41,
    kTNF__Y1 = 42,
    kTNF__Y2 = 43,
    kTNFNFC1 = 44,
    kTNFNFC2 = 45,
    kTNFNFY1 = 46,
    kTNFNFY2 = 47,
    kTYY__C1 = 48,
    kTYY__C2 = 49,
    kTYY__Y1 = 50,
    kTYY__Y2 = 51,
    kTYYYYC1 = 52,
    kTYYYYC2 = 53,
    kTYYYYY1 = 54,
    kTYYYYY2 = 55,
    kTYF__C1 = 56,
    kTYF__C2 = 57,
    kTYF__Y1 = 58,
    kTYF__Y2 = 59,
    kTYFYFC1 = 60,
    kTYFYFC2 = 61,
    kTYFYFY1 = 62,
    kTYFYFY2 = 63,

    kWild = kKNY__C1,

    // old definition
    kNoParasite = 64,
    kRa = kKNY__Y1,
    kRb = kKNY__C2,
    kR0 = kWild,
    kRab = kKNY__Y2,
    First = kKNY__C1,
    Last = kNoParasite

};
std::ostream& operator<<(std::ostream& os, ParasiteType pt);
const size_t kParasiteType_size = static_cast<size_t>(
                static_cast<int>(ParasiteType::Last)
                - static_cast<int>(ParasiteType::First) + 1 
             );

const std::unordered_map<std::string, ParasiteType> kStr_to_pt_map = {
    {"KNY--C1", ParasiteType::kKNY__C1}, {"wild", ParasiteType::kKNY__C1}, 
    {"KNY--C2", ParasiteType::kKNY__C2},
    {"KNY--Y1", ParasiteType::kKNY__Y1}, {"580Y", ParasiteType::kKNY__Y1},
    {"KNY--Y2", ParasiteType::kKNY__Y2},
    {"KNYNYC1", ParasiteType::kKNYNYC1},
    {"KNYNYC2", ParasiteType::kKNYNYC2},
    {"KNYNYY1", ParasiteType::kKNYNYY1},
    {"KNYNYY2", ParasiteType::kKNYNYY2},
    {"KNF--C1", ParasiteType::kKNF__C1},
    {"KNF--C2", ParasiteType::kKNF__C2},
    {"KNF--Y1", ParasiteType::kKNF__Y1},
    {"KNF--Y2", ParasiteType::kKNF__Y2},
    {"KNFNFC1", ParasiteType::kKNFNFC1},
    {"KNFNFC2", ParasiteType::kKNFNFC2},
    {"KNFNFY1", ParasiteType::kKNFNFY1},
    {"KNFNFY2", ParasiteType::kKNFNFY2},
    {"KYY--C1", ParasiteType::kKYY__C1},
    {"KYY--C2", ParasiteType::kKYY__C2},
    {"KYY--Y1", ParasiteType::kKYY__Y1},
    {"KYY--Y2", ParasiteType::kKYY__Y2},
    {"KYYYYC1", ParasiteType::kKYYYYC1},
    {"KYYYYC2", ParasiteType::kKYYYYC2},
    {"KYYYYY1", ParasiteType::kKYYYYY1},
    {"KYYYYY2", ParasiteType::kKYYYYY2},
    {"KYF--C1", ParasiteType::kKYF__C1},
    {"KYF--C2", ParasiteType::kKYF__C2},
    {"KYF--Y1", ParasiteType::kKYF__Y1},
    {"KYF--Y2", ParasiteType::kKYF__Y2},
    {"KYFYFC1", ParasiteType::kKYFYFC1},
    {"KYFYFC2", ParasiteType::kKYFYFC2},
    {"KYFYFY1", ParasiteType::kKYFYFY1},
    {"KYFYFY2", ParasiteType::kKYFYFY2},
    {"TNY--C1", ParasiteType::kTNY__C1},
    {"TNY--C2", ParasiteType::kTNY__C2},
    {"TNY--Y1", ParasiteType::kTNY__Y1},
    {"TNY--Y2", ParasiteType::kTNY__Y2},
    {"TNYNYC1", ParasiteType::kTNYNYC1},
    {"TNYNYC2", ParasiteType::kTNYNYC2},
    {"TNYNYY1", ParasiteType::kTNYNYY1},
    {"TNYNYY2", ParasiteType::kTNYNYY2},
    {"TNF--C1", ParasiteType::kTNF__C1},
    {"TNF--C2", ParasiteType::kTNF__C2},
    {"TNF--Y1", ParasiteType::kTNF__Y1},
    {"TNF--Y2", ParasiteType::kTNF__Y2},
    {"TNFNFC1", ParasiteType::kTNFNFC1},
    {"TNFNFC2", ParasiteType::kTNFNFC2},
    {"TNFNFY1", ParasiteType::kTNFNFY1},
    {"TNFNFY2", ParasiteType::kTNFNFY2},
    {"TYY--C1", ParasiteType::kTYY__C1},
    {"TYY--C2", ParasiteType::kTYY__C2},
    {"TYY--Y1", ParasiteType::kTYY__Y1},
    {"TYY--Y2", ParasiteType::kTYY__Y2},
    {"TYYYYC1", ParasiteType::kTYYYYC1},
    {"TYYYYC2", ParasiteType::kTYYYYC2},
    {"TYYYYY1", ParasiteType::kTYYYYY1},
    {"TYYYYY2", ParasiteType::kTYYYYY2},
    {"TYF--C1", ParasiteType::kTYF__C1},
    {"TYF--C2", ParasiteType::kTYF__C2},
    {"TYF--Y1", ParasiteType::kTYF__Y1},
    {"TYF--Y2", ParasiteType::kTYF__Y2},
    {"TYFYFC1", ParasiteType::kTYFYFC1},
    {"TYFYFC2", ParasiteType::kTYFYFC2},
    {"TYFYFY1", ParasiteType::kTYFYFY1},
    {"TYFYFY2", ParasiteType::kTYFYFY2}

};

class GenotypeManager {

    const int kNum_sites_in_genotype_code;
    const int kNum_genotypes;
    const bool kAllow_back_mutation;

    // Probability of mutation at one site regardless of other sites
    // float* kMutation_independent_probability;
    // float kMutation_independent_probability_sum = 0;
    std::vector<float> kMutation_independent_probability;
    std::vector<float> kMutation_independent_probability_positive;

    // Probability of combination of mutation events at all sites
    const int kMax_num_sites_in_combined_mutation;
    const int kNum_mutation_combinations;
    // float* kMutation_combined_probability;
    // genotype_t* kMutation_combined_masks;
    std::vector<float> kMutation_combined_probability;
    std::vector<float> kMutation_combined_probability_positives;
    std::vector<genotype_t> kMutation_combined_masks;
    // std::vector<int> kMutation_counts;


public:

    static constexpr int kNum_bits_in_a_byte = 8;
    static constexpr int kNum_bits_in_genotype = sizeof(genotype_t) * kNum_bits_in_a_byte;

    GenotypeManager(
        const std::vector<float>& mutation_probability_vector,
        int max_num_sites_in_combined_mutation,
        bool allow_back_mutation
    );
    GenotypeManager(
        const std::vector<float>& mutation_probability_vector,
        int max_num_sites_in_combined_mutation
    );
    GenotypeManager(
        const std::vector<float>& mutation_probability_vector
    );
    ~GenotypeManager();

    void print_all_combined_probabilities() const;
    std::vector<float> get_all_combined_probabilities() const;

    void step_mutation(
        genotype_t* genotype_array,
        const int genotype_array_size
    );
    static void step_mutation_independent_kernel(
        genotype_t* genotype_array,
        const int genotype_array_size,

        const float* mutation_probability_independent_array,
        const int mutation_probability_array_size
    );

    int fun_mutation_combined (
        ParasiteType& pt
    ) const;
    int fun_mutation_combined (
        genotype_t& g
    ) const;
    int fun_mutation_individual(
        genotype_t& g
    ) const;


    inline int get_num_genotypes() const {
        return this->kNum_genotypes;
    }

    static inline bool gt_has_mutations(
        genotype_t gt,
        genotype_t mutations
        ) {
        return ((gt & mutations) == mutations);
    }

    static ParasiteType map_str_to_pt(std::string);

    template <typename T>
    static T sum_vector_on_sites(const std::vector<T>& vec, ParasiteType pt);

    template <typename T>
    static std::vector<T> sum_vector_on_sites(
                            const std::vector<T>& vec,
                            const std::vector<ParasiteType>& pt_list
                        );

};



}
#endif
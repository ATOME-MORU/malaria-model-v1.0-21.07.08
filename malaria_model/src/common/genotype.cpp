#include <iostream>
#include <vector>
#include <algorithm>
#include <cstdint> //uint8_t

#include <numeric> //accumulate
#include <cmath> //pow
#include <cassert>
#include <stdexcept>

#include <bitset>

#include "util/randomness.h"
#include "util/statistics.h"

#include "common/genotype.h"

namespace common {

std::ostream& operator<<(std::ostream& os, ParasiteType pt) {
    switch(pt) {
        case  ParasiteType::kKNY__C1 : os << "KNY--C1(0w)"; break;
        case  ParasiteType::kKNY__C2 : os << "KNY--C2(1)"; break;
        case  ParasiteType::kKNY__Y1 : os << "KNY--Y1(2)"; break;
        case  ParasiteType::kKNY__Y2 : os << "KNY--Y2(3)"; break;
        case  ParasiteType::kKNYNYC1 : os << "KNYNYC1(4)"; break;
        case  ParasiteType::kKNYNYC2 : os << "KNYNYC2(5)"; break;
        case  ParasiteType::kKNYNYY1 : os << "KNYNYY1(6)"; break;
        case  ParasiteType::kKNYNYY2 : os << "KNYNYY2(7)"; break;
        case  ParasiteType::kKNF__C1 : os << "KNF--C1(8)"; break;
        case  ParasiteType::kKNF__C2 : os << "KNF--C2(9)"; break;
        case  ParasiteType::kKNF__Y1 : os << "KNF--Y1(10)"; break;
        case  ParasiteType::kKNF__Y2 : os << "KNF--Y2(11)"; break;
        case  ParasiteType::kKNFNFC1 : os << "KNFNFC1(12)"; break;
        case  ParasiteType::kKNFNFC2 : os << "KNFNFC2(13)"; break;
        case  ParasiteType::kKNFNFY1 : os << "KNFNFY1(14)"; break;
        case  ParasiteType::kKNFNFY2 : os << "KNFNFY2(15)"; break;
        case  ParasiteType::kKYY__C1 : os << "KYY--C1(16)"; break;
        case  ParasiteType::kKYY__C2 : os << "KYY--C2(17)"; break;
        case  ParasiteType::kKYY__Y1 : os << "KYY--Y1(18)"; break;
        case  ParasiteType::kKYY__Y2 : os << "KYY--Y2(19)"; break;
        case  ParasiteType::kKYYYYC1 : os << "KYYYYC1(20)"; break;
        case  ParasiteType::kKYYYYC2 : os << "KYYYYC2(21)"; break;
        case  ParasiteType::kKYYYYY1 : os << "KYYYYY1(22)"; break;
        case  ParasiteType::kKYYYYY2 : os << "KYYYYY2(23)"; break;
        case  ParasiteType::kKYF__C1 : os << "KYF--C1(24)"; break;
        case  ParasiteType::kKYF__C2 : os << "KYF--C2(25)"; break;
        case  ParasiteType::kKYF__Y1 : os << "KYF--Y1(26)"; break;
        case  ParasiteType::kKYF__Y2 : os << "KYF--Y2(27)"; break;
        case  ParasiteType::kKYFYFC1 : os << "KYFYFC1(28)"; break;
        case  ParasiteType::kKYFYFC2 : os << "KYFYFC2(29)"; break;
        case  ParasiteType::kKYFYFY1 : os << "KYFYFY1(30)"; break;
        case  ParasiteType::kKYFYFY2 : os << "KYFYFY2(31)"; break;
        case  ParasiteType::kTNY__C1 : os << "TNY--C1(32)"; break;
        case  ParasiteType::kTNY__C2 : os << "TNY--C2(33)"; break;
        case  ParasiteType::kTNY__Y1 : os << "TNY--Y1(34)"; break;
        case  ParasiteType::kTNY__Y2 : os << "TNY--Y2(35)"; break;
        case  ParasiteType::kTNYNYC1 : os << "TNYNYC1(36)"; break;
        case  ParasiteType::kTNYNYC2 : os << "TNYNYC2(37)"; break;
        case  ParasiteType::kTNYNYY1 : os << "TNYNYY1(38)"; break;
        case  ParasiteType::kTNYNYY2 : os << "TNYNYY2(39)"; break;
        case  ParasiteType::kTNF__C1 : os << "TNF--C1(40)"; break;
        case  ParasiteType::kTNF__C2 : os << "TNF--C2(41)"; break;
        case  ParasiteType::kTNF__Y1 : os << "TNF--Y1(42)"; break;
        case  ParasiteType::kTNF__Y2 : os << "TNF--Y2(43)"; break;
        case  ParasiteType::kTNFNFC1 : os << "TNFNFC1(44)"; break;
        case  ParasiteType::kTNFNFC2 : os << "TNFNFC2(45)"; break;
        case  ParasiteType::kTNFNFY1 : os << "TNFNFY1(46)"; break;
        case  ParasiteType::kTNFNFY2 : os << "TNFNFY2(47)"; break;
        case  ParasiteType::kTYY__C1 : os << "TYY--C1(48)"; break;
        case  ParasiteType::kTYY__C2 : os << "TYY--C2(49)"; break;
        case  ParasiteType::kTYY__Y1 : os << "TYY--Y1(50)"; break;
        case  ParasiteType::kTYY__Y2 : os << "TYY--Y2(51)"; break;
        case  ParasiteType::kTYYYYC1 : os << "TYYYYC1(52)"; break;
        case  ParasiteType::kTYYYYC2 : os << "TYYYYC2(53)"; break;
        case  ParasiteType::kTYYYYY1 : os << "TYYYYY1(54)"; break;
        case  ParasiteType::kTYYYYY2 : os << "TYYYYY2(55)"; break;
        case  ParasiteType::kTYF__C1 : os << "TYF--C1(56)"; break;
        case  ParasiteType::kTYF__C2 : os << "TYF--C2(57)"; break;
        case  ParasiteType::kTYF__Y1 : os << "TYF--Y1(58)"; break;
        case  ParasiteType::kTYF__Y2 : os << "TYF--Y2(59)"; break;
        case  ParasiteType::kTYFYFC1 : os << "TYFYFC1(60)"; break;
        case  ParasiteType::kTYFYFC2 : os << "TYFYFC2(61)"; break;
        case  ParasiteType::kTYFYFY1 : os << "TYFYFY1(62)"; break;
        case  ParasiteType::kTYFYFY2 : os << "TYFYFY2(63)"; break;

        // case ParasiteType::kWild    : os << "Wild"; break; //KNY--C1

        case ParasiteType::kNoParasite  : os << "NP(64)"; break;
        // case ParasiteType::kRa          : os << "RA"; break; //KNY--Y1
        // case ParasiteType::kRb          : os << "RB"; break; //KNY--C2
        // case ParasiteType::kR0          : os << "R0"; break; //KNY--C1
        // case ParasiteType::kRab         : os << "RAB"; break;//KNY--Y2
        default :
            // if (pt < ParasiteType::Last && pt > ParasiteType::First) {
            //     os << "unnamed";
            // } else {
                os << "Exception: Unknown ParasiteType:" << static_cast<int>(pt);
                throw std::out_of_range("Unknown ParasiteType");
            // }
    }
    return os;
}

GenotypeManager::GenotypeManager(

        const std::vector<float>& mutation_probability_vector,
        int max_num_sites_in_combined_mutation,
        bool allow_back_mutation
    
    ) : kNum_sites_in_genotype_code(static_cast<int>(mutation_probability_vector.size())),
        kNum_genotypes(pow(2,kNum_sites_in_genotype_code)),
        kAllow_back_mutation(allow_back_mutation),

        kMutation_independent_probability(mutation_probability_vector),
        kMutation_independent_probability_positive(std::vector<float>(kNum_sites_in_genotype_code, 0.0)),
        kMax_num_sites_in_combined_mutation(std::min(
                max_num_sites_in_combined_mutation,
                kNum_sites_in_genotype_code
        )),
        kNum_mutation_combinations(
            util::from_n_choose_up_to_k_get_num_combinations(
                kNum_sites_in_genotype_code,
                kMax_num_sites_in_combined_mutation
            )
        ),
        kMutation_combined_probability(std::vector<float>(kNum_mutation_combinations, 1.0)),
        kMutation_combined_probability_positives(std::vector<float>(kNum_mutation_combinations-1, 1.0)),
        kMutation_combined_masks(std::vector<genotype_t>(kNum_mutation_combinations, 0))//,
        // kMutation_counts(std::vector<int>(kNum_mutation_combinations,0))
        {

    // Independent Mutation Probability
    float kMutation_independent_probability_sum = std::accumulate(
        this->kMutation_independent_probability.begin(),
        this->kMutation_independent_probability.end(),
        0.0
    );
    for( size_t pp = 0; pp < this->kMutation_independent_probability_positive.size(); pp++) {
        this->kMutation_independent_probability_positive[pp] =
            this->kMutation_independent_probability[pp] / kMutation_independent_probability_sum;
        // std::cout << "pp = " << pp << "\n";
    }

    // Combined Mutation Probability
    // this->kMutation_combined_probability = new float[this->kNum_mutation_combinations];
    // this->kMutation_combined_masks = new genotype_t[this->kNum_mutation_combinations];
    // std::fill(this->kMutation_combined_probability,
    //     this->kMutation_combined_probability +kNum_mutation_combinations,
    //     1
    // );
    // std::fill(this->kMutation_combined_masks,
    //     this->kMutation_combined_masks +kNum_mutation_combinations,
    //     0
    // );

    int num_combinations_accumulator = 0;
    // no mutation at any site
    this->kMutation_combined_masks[0] = 0;
    for (int ss = 0; ss < kNum_sites_in_genotype_code; ss++) {
        this->kMutation_combined_probability[0]
            *= (1 - this->kMutation_independent_probability[ss]);
    }
    // this->kMutation_counts[0] = 0;
    num_combinations_accumulator++;

    // mutation at least one site
    for(int kk = 1; kk <= kMax_num_sites_in_combined_mutation; kk++) {
        std::vector<std::vector<int>> choose_kk_combination_set = util::from_n_choose_k(kNum_sites_in_genotype_code, kk);

        for (auto& each_combination : choose_kk_combination_set) {

            // for (auto cc : each_combination) {
            //     std::cout << cc << ",";
            // }
            // std::cout << "\n";
            
            // for (int ss = 0; ss < kNum_sites_in_genotype_code; ss++) {
            for (int ss = kNum_sites_in_genotype_code-1; ss >=0 ; ss--) {
                if (std::find(each_combination.begin(), each_combination.end(), ss) != each_combination.end()) {
                    this->kMutation_combined_masks[num_combinations_accumulator] <<= 1;
                    this->kMutation_combined_masks[num_combinations_accumulator] += 1;
                    // this->kMutation_counts[num_combinations_accumulator] += 1;

                    this->kMutation_combined_probability[num_combinations_accumulator]
                        *= this->kMutation_independent_probability[ss];

                } else {

                    this->kMutation_combined_masks[num_combinations_accumulator] <<= 1;

                    this->kMutation_combined_probability[num_combinations_accumulator]
                        *= (1 - this->kMutation_independent_probability[ss]);

                }
            }
            this->kMutation_combined_probability_positives[num_combinations_accumulator-1] =
                this->kMutation_combined_probability[num_combinations_accumulator];
            num_combinations_accumulator++;
        }
    }

}
GenotypeManager::GenotypeManager (
        const std::vector<float>& mutation_probability_vector,
        int max_num_sites_in_combined_mutation
    ) : GenotypeManager(
        mutation_probability_vector,
        max_num_sites_in_combined_mutation,
        true
    ) {

}
GenotypeManager::GenotypeManager (
        const std::vector<float>& mutation_probability_vector
    ) : GenotypeManager(
        mutation_probability_vector,
        static_cast<int>(mutation_probability_vector.size()),
        true // kAllow_back_mutation default
    ) {

}

GenotypeManager::~GenotypeManager() {
    // if (kMutation_independent_probability) delete[] kMutation_independent_probability;
    // if (kMutation_combined_probability) delete[] kMutation_combined_probability;
    // if (kMutation_combined_masks) delete[] kMutation_combined_masks;
}

void GenotypeManager::print_all_combined_probabilities() const {

    std::cout << "Given that the *independent* mutation probabilities at each of the "
            << this->kNum_sites_in_genotype_code << " sites are as [";
    for (int ss = 0; ss < this->kNum_sites_in_genotype_code; ss++) {
        std::cout << " " << this->kMutation_independent_probability[ss];
    }
    std::cout << " ]\n";

    std::cout << "Followings are the probabilities of possible combinations of mutations:\n";
    std::cout << "(Max number of concurrent mutation events: " << kMax_num_sites_in_combined_mutation << ")\n";

    std::cout << "[combined mutation event: probability]\n";
    for (int pp = 0; pp < this->kNum_mutation_combinations; pp++) {
        genotype_t mm = this->kMutation_combined_masks[pp];
        std::cout << std::bitset<this->kNum_bits_in_genotype>(mm)
                    << "(" << unsigned(mm) << ")"
                    // << "(" << std::bitset<this->kNum_bits_in_genotype>(mm)[0] << ")"
                    << ": " << this->kMutation_combined_probability[pp] << "\n";
    }
    std::cout << "with [" << std::bitset<this->kNum_bits_in_genotype>(0)
        << "] being the situation where no mutation occurs at any sites.\n";
}

std::vector<float> GenotypeManager::get_all_combined_probabilities() const {
    std::vector<float> all_probabilities;
    for (int pp = 0; pp < this->kNum_mutation_combinations; pp++) {
        all_probabilities.push_back(this->kMutation_combined_probability[pp]);
    }
    return all_probabilities;
}



void GenotypeManager::step_mutation(
        genotype_t* genotype_array,
        const int genotype_array_size
    ){

    for (int gg = 0; gg < genotype_array_size; gg++) {
        this->fun_mutation_combined(genotype_array[gg]);
        // this->fun_mutation_individual(genotype_array[gg]);
    }

    // this->step_mutation_independent_kernel(
    //     genotype_array,
    //     genotype_array_size,

    //     this->kMutation_independent_probability,
    //     this->kNum_sites_in_genotype_code
    // );

}

void GenotypeManager::step_mutation_independent_kernel(
        genotype_t* genotype_array,
        const int genotype_array_size,

        const float* mutation_probability_independent_array,
        const int mutation_probability_array_size
    ) {

    genotype_t mask = 0;
    float* rnd_nmb_array = new float[mutation_probability_array_size];

    for (int gg = 0; gg < genotype_array_size; gg++) {

        mask = 0;
        util::set_random_numbers_uniform(rnd_nmb_array, mutation_probability_array_size);
        
        for (int bb = 0; bb < mutation_probability_array_size; bb++) {
            mask <<= 1;
            mask += (rnd_nmb_array[bb] < mutation_probability_independent_array[bb]) ? 1 : 0;
        }

        genotype_array[gg] ^= mask; // back mutation

    }

    delete[] rnd_nmb_array;

}

int GenotypeManager::fun_mutation_combined(
        ParasiteType& pt
    ) const {

    genotype_t g = static_cast<genotype_t>(pt);
    int c = this->fun_mutation_combined(g);
    pt = static_cast<ParasiteType>(g);
    return c;
}

int GenotypeManager::fun_mutation_combined(
        genotype_t& g
    ) const {

    if (util::get_rand_uniform() < this->kMutation_combined_probability[0]) {
        return 0;
    } else {
        genotype_t g_orig = g;

        int mut_mask_idx = 1 + util::get_rand_discrete(this->kMutation_combined_probability_positives);
        genotype_t mut_mask = this->kMutation_combined_masks[mut_mask_idx];

        if (this->kAllow_back_mutation) {
            g ^= mut_mask;
        } else {
            g |= mut_mask;
        }

        // returns the number of sites that mutated in this mutation event
        // not actual number when back mutation is not allowed
        // return this->kMutation_counts[mut_mask_idx];

        return std::bitset<this->kNum_bits_in_genotype>(g^g_orig).count();

    }

}

int GenotypeManager::fun_mutation_individual(
        genotype_t& g
    ) const {

    assert(false && "fun_mutation_individual should not be used");

    if (util::get_rand_uniform() < this->kMutation_combined_probability[0]) {
        return 0;
    } else {
        genotype_t mask = 0;
        int num_mut = 0;
        
        for (int bb = 5; bb >= 0; bb--) {
            mask <<= 1;
            if (util::get_rand_uniform()
                < this->kMutation_independent_probability_positive[bb]
            ) {
                mask += 1;
                num_mut += 1;
            }
        }

        g |= mask;

        return num_mut;

    }

}

ParasiteType GenotypeManager::map_str_to_pt(
        std::string str
    ) {
    auto it_res = kStr_to_pt_map.find(str);
    if (it_res != kStr_to_pt_map.end()) {
        return it_res->second;
    }
    std::cout << "invalid string value: " << str << "\n";
    throw std::invalid_argument("Can't map str to parasite type.");
}

template <typename T>
T GenotypeManager::sum_vector_on_sites(const std::vector<T>& vec, ParasiteType pt) {
    assert(vec.size() == kParasiteType_size);
    T sum = 0.0;
    for (genotype_t ii = 0; ii < kParasiteType_size; ii++) {
        // if ( (ii & static_cast<genotype_t>(pt)) > 0 ){ // sum_vector_on_sites_rm1
        if ( gt_has_mutations(ii, static_cast<genotype_t>(pt)) ) {
            sum += vec[ii];
            // if (sum >0.0) {
            //     std::cout << "ii=" << (int)ii << ", vec[ii]=" << vec[ii]
            //                 << ", sum=" << sum << "\n";
            // }
        }
    }
    return sum;
}
template int GenotypeManager::sum_vector_on_sites(const std::vector<int>& vec, ParasiteType pt);
template float GenotypeManager::sum_vector_on_sites(const std::vector<float>& vec, ParasiteType pt);

template <typename T>
std::vector<T> GenotypeManager::sum_vector_on_sites(
        const std::vector<T>& vec,
        const std::vector<ParasiteType>& pt_list
    ) {
    assert(vec.size() == kParasiteType_size);
    std::vector<T> sum_list(pt_list.size(), 0.0);
    for (size_t pp = 0; pp < pt_list.size(); pp++) {
        sum_list[pp] = sum_vector_on_sites(vec, pt_list[pp]);
    }

    // for (genotype_t ii = 0; ii < kParasiteType_size; ii++) {
    //     for (size_t pp = 0; pp < pt_list.size(); pp++) {
    //         if ( (ii & static_cast<genotype_t>(pt_list[pp])) > 0 ){
    //             sum_list[pp] += vec[ii];
    //         }
    //     }
    // }
    return sum_list;
}
template std::vector<int> GenotypeManager::sum_vector_on_sites(
    const std::vector<int>& vec, const std::vector<ParasiteType>& pt
);
template std::vector<float> GenotypeManager::sum_vector_on_sites(
    const std::vector<float>& vec, const std::vector<ParasiteType>& pt
);


}
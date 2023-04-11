#ifndef DRUG_H
#define DRUG_H

#include <array>

namespace common {

// enum DrugName {
// enum class DrugName: int {
enum class DrugName: char {
    kNoDrug = 0,
    kAS  = 1,
    kLM  = 2,
    kAQ  = 3,
    kPPQ = 4,
    kMQ  = 5,
    kSP  = 6,
    kCQ  = 7,
    kAL  = 8,
    kALb = 9,
    kASAQ  = 10,
    kASAQb = 11,
    kDP  = 12,
    kDPb = 13,
    kASMQ  = 14,
    kASMQb = 15,
    kASMQPPQ  = 16,
    kASMQbPPQ = 17,
    kALAQ  = 18,
    kALbAQ = 19,

    /*1*/    kArtesunate = kAS,
    /*2*/    kPiperaquine = kDPb,
    /*3*/    kPiparte = kDP,
    // /*4*/    kAQ,
    // /*5*/    kASAQ,
    /*6*/    kLumefantrine = kLM,
    // /*7*/    kAL,
    // /*8*/    kSP,
    /*9*/    kChloroquine = kCQ,

    First = kAS,
    Last = kALbAQ
};
std::ostream& operator<<(std::ostream& os, DrugName drg);


class DrugManager {

    static const int kNumDrugs = static_cast<int>(DrugName::Last) - static_cast<int>(DrugName::First) + 1;
    static const int kNumGenotypesMax = 64;

    std::array<std::array<float, kNumDrugs+1>, kNumGenotypesMax> kGd; // [Genotype][DrugName] = daily clearance prob
    std::array<float, kNumDrugs+1> kDl; //[DrugName] = daily drug loss rate

    std::array<DrugName, kNumDrugs+1> kDlOutcome {{
            /*00*/   DrugName::kNoDrug,
            /*01*/   DrugName::kNoDrug,
            /*02*/   DrugName::kNoDrug,
            /*03*/   DrugName::kNoDrug,
            /*04*/   DrugName::kNoDrug,
            /*05*/   DrugName::kNoDrug,
            /*06*/   DrugName::kNoDrug,
            /*07*/   DrugName::kNoDrug,
            /*08*/   DrugName::kALb,
            /*09*/   DrugName::kNoDrug,
            /*10*/   DrugName::kASAQb,
            /*11*/   DrugName::kNoDrug,
            /*12*/   DrugName::kDPb,
            /*13*/   DrugName::kNoDrug,
            /*14*/   DrugName::kASMQb,
            /*15*/   DrugName::kNoDrug,
            /*16*/   DrugName::kASMQbPPQ,
            /*17*/   DrugName::kPPQ, // or kASMQb
            /*18*/   DrugName::kALbAQ,
            /*19*/   DrugName::kAQ  // or kALb
        }};

    // std

public:

    DrugManager(
        std::string gd_matrix_file,
        std::string dl_array_file,
        bool verbose
    );
    DrugManager(
        std::string gd_matrix_file,
        std::string dl_array_file
    );

    bool drug_effect(
        genotype_t gt,
        DrugName dg,
        StageName stg
    ) const;
    int compare_drug_effect_rate(
        genotype_t gt_orig, // original genotype
        genotype_t gt_prop, // proposed genotype
        DrugName dg        // drug presence
    ) const;

    bool drug_loss(DrugName& dg) const;

    /////////////////////////
    // Access
    inline float get_drug_effect_rate(genotype_t gt, DrugName dg) const {
        return this->kGd[static_cast<int>(gt)][static_cast<int>(dg)];
    }
    inline float get_drug_loss_rate(DrugName dg) const {
        return this->kDl[static_cast<int>(dg)];
    }
    inline DrugName get_drug_loss_result(DrugName dg) const {
        return this->kDlOutcome[static_cast<int>(dg)];
    }

};

}

#endif
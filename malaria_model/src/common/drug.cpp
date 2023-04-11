#include <iostream>
#include <vector>
#include <array>
#include <cassert>

#include <fstream> //std::ifstream
#include <sstream> //std::stringstream

#include <string>

#include <stdexcept> //std::invalid_argument

#include "util/randomness.h"
#include "common/common.h"

namespace common {

std::ostream& operator<<(std::ostream& os, DrugName dg) {
    switch(dg) {
        case DrugName::kNoDrug  : os << "NDG"; break;
        case DrugName::kAS    : os << "AS";   break;
        case DrugName::kLM    : os << "LM";   break;
        case DrugName::kAQ    : os << "AQ";   break;
        case DrugName::kPPQ   : os << "PPQ";  break;
        case DrugName::kMQ    : os << "MQ";   break;
        case DrugName::kSP    : os << "SP";   break;
        case DrugName::kCQ    : os << "CQ";   break;
        case DrugName::kAL    : os << "AL";   break;
        case DrugName::kALb   : os << "ALb";  break;
        case DrugName::kASAQ  : os << "ASAQ";     break;
        case DrugName::kASAQb : os << "ASAQb";    break;
        case DrugName::kDP    : os << "DP";   break;
        case DrugName::kDPb   : os << "DPb";  break;
        case DrugName::kASMQ  : os << "ASMQ";     break;
        case DrugName::kASMQb : os << "ASMQb";    break;
        case DrugName::kASMQPPQ   : os << "ASMQPPQ";  break;
        case DrugName::kASMQbPPQ  : os << "ASMQbPPQ";  break;
        case DrugName::kALAQ  : os << "ALAQ";     break;
        case DrugName::kALbAQ : os << "ALbAQ";    break;

        // case DrugName::First : os << "NDG";    break;
        // case DrugName::Last : os << "ALbAQ";    break;
        
        default :
            os << "Exception: Unknown DrugName " << static_cast<int>(dg);
            throw std::out_of_range("Unknown DrugName");
    }
    return os;
}


DrugManager::DrugManager(
        std::string gd_matrix_file,
        std::string dl_array_file,
        bool verbose
    ) {

    if (verbose) {
        std::cout << "DM: initilisation, ";
        std::cout << "#drugs: " << this->kNumDrugs << "\n";
        std::cout << "DM: reading gd matrix from " + gd_matrix_file + "\n";
    }

    std::ifstream input_file_stream(gd_matrix_file);

    if (!input_file_stream) {

        throw std::invalid_argument("Can not open file: " + gd_matrix_file);

    } else {


        std::string str_row;

        int r = 0;
        while (std::getline(input_file_stream, str_row)) {
            // std::cout << str_row << "\n";
            if (r >= this->kNumGenotypesMax) {
                throw "DM: input Gd matrix longer than expected. (Expecting "
                        + std::to_string(this->kNumGenotypesMax)
                        + ", got more than that)\n";
            }

            std::string str_cell;
            std::stringstream ss(str_row);

            char delimiter = ',';

            this->kGd[r][0] = 0.0;
            int c = 1; //
            while (std::getline(ss, str_cell, delimiter)) {
                
                if (c > this->kNumDrugs) {
                    throw "DM: input Gd matrix wider than expected at row"
                            + std::to_string(r) + ": expecting "
                            + std::to_string(this->kNumDrugs)
                            + " columns, got at least" + std::to_string(c) + ")\n";
                }
                try {
                    this->kGd[r][c++] = std::stof(str_cell);
                } catch (const std::invalid_argument& e ) {
                    throw "DM: invalid value " + str_cell
                        + " at r-" + std::to_string(r)
                        + ", c-"  + std::to_string(c)
                        + "\n";
                }
                
            }
            // std::cout << "\n";
            if (c != this->kNumDrugs + 1) {
                throw  "DM: input Gd matrix narrower than expected at row "
                        + std::to_string(r) + ": expecting "
                        + std::to_string(this->kNumDrugs)
                        + " columns, got " + std::to_string(c-1) + "\n";
            }

            r++;
        }
        if (r != this->kNumGenotypesMax) {
            throw "DM: input Gd matrix shorter than expeted : expecting "
                    + std::to_string(this->kNumGenotypesMax)
                    + ", got " + std::to_string(r) + "\n";
        }

        if (verbose) {
            std::cout << "\tkGd[0][DrugName::First]: "
                    << this->kGd[0][static_cast<int>(DrugName::First)]
                    << "\n";
            std::cout << "\tkGd[kNumGenotypesMax-1][DrugName::Last]: "
                    << this->kGd[kNumGenotypesMax-1][static_cast<int>(DrugName::Last)]
                    << "\n";
        
            // DrugName dg = DrugName::kDP;
            std::vector<DrugName> vdg{DrugName::kDP, DrugName::kALbAQ};

            for (const auto& dg : vdg) {
                std::cout << "\tkGd[" << kNumGenotypesMax-1 << "][" << dg << "]: "
                        << this->kGd[kNumGenotypesMax-1][static_cast<int>(dg)] << "\n";
            }
        
            std::cout << "DM: read " << r << " lines from " << gd_matrix_file << "\n";
        }
    }

    input_file_stream.close();

    if (verbose) {
        std::cout << "DM: reading dl array from " + dl_array_file + "\n";
    }

    input_file_stream.open(dl_array_file);
    if (!input_file_stream) {
        throw std::invalid_argument("can not open file:" + dl_array_file + "\n");
    } else {
        std::string str_row;
        this->kDl[0] = 0.0;
        int d = 1;
        while (std::getline(input_file_stream, str_row)) {
            std::string str_cell;
            std::istringstream iss(str_row);
            char delimiter = ',';
            while (std::getline(iss, str_cell, delimiter)) {
                if (d > this->kNumDrugs) {
                    throw "DM: Input Dl array has more elements than expected. (Expecting "
                            + std::to_string(this->kNumDrugs)
                            + ", got at least " + std::to_string(d) + ")\n";
                }
                try {
                    float loss_num_days = std::stof(str_cell);
                    if (loss_num_days <= 0.0) {
                        throw "DM: Drug #" + std::to_string(d)
                                + " has invalid loss_num_days of "
                                + std::to_string(loss_num_days)
                                + ".\n";
                    }
                    this->kDl[d++] = 1.0/loss_num_days;
                } catch (const std::invalid_argument& e) {
                    throw "DM: Invalid value '" + str_cell
                        + "' for drug #" + std::to_string(d)
                        + "\n";
                }
            }
        }
        d--;
        if (d < this->kNumDrugs) {
            throw "DM: Input Dl array shorter than expected. (Expecting"
                    + std::to_string(this->kNumDrugs)
                    + ", got " + std::to_string(d-1)
                    + "\n";
        }
        if (verbose) {
            std::vector<DrugName> vdg{DrugName::kDP, DrugName::kALAQ, DrugName::kALbAQ};
            for (const auto& dg : vdg) {
                std::cout << "\tkDl[" << dg << "]:" << this->kDl[static_cast<int>(dg)] << "\n";
            }
            std::cout << "DM: read " << d << " elements from " << gd_matrix_file << "\n";
        }
    }

    if (verbose) {
        std::cout << "DM: init complete. \n";
    }

}
DrugManager::DrugManager(
    std::string gd_matrix_file,
    std::string dl_array_file
    ) : DrugManager(
            gd_matrix_file,
            dl_array_file,
            false // verbose
    ) {

}

bool DrugManager::drug_effect(
        genotype_t gt,
        DrugName dg,
        StageName stg
    ) const {
    if (stg == StageName::kLiver) {return false;}
    assert(stg != StageName::kNotInSystem);
    return (util::get_rand_uniform() < this->kGd[static_cast<int>(gt)][static_cast<int>(dg)]);
}

int DrugManager::compare_drug_effect_rate(
        genotype_t gt_orig,
        genotype_t gt_prop,
        DrugName dg
    ) const {
    if (dg == DrugName::kNoDrug) {
        return 0;
    }
    float rate_orig = this->get_drug_effect_rate(gt_orig, dg);
    float rate_prop = this->get_drug_effect_rate(gt_prop, dg);
    if (rate_prop > rate_orig) {
        return 1;
    }
    if (rate_prop < rate_orig) {
        return -1;
    }
    return 0;
}

bool DrugManager::drug_loss(DrugName& dg) const {
    if (dg == DrugName::kNoDrug) {
        return false;
    } else {
        if (dg == DrugName::kASMQbPPQ) {
            bool l1 = util::get_rand_uniform() < this->kDl[static_cast<int>(DrugName::kASMQb)];
            bool l2 = util::get_rand_uniform() < this->kDl[static_cast<int>(DrugName::kPPQ)];
            if (l1 && l2) {
                dg = DrugName::kNoDrug;
                return true;
            }
            if (!l1 && l2) {
                dg = DrugName::kASMQb;
                return true;
            }
            if (l1 && !l2) {
                dg = DrugName::kPPQ;
                return true;
            }
            return false;
        } else if (dg == DrugName::kALbAQ) {
            bool l1 = util::get_rand_uniform() < this->kDl[static_cast<int>(DrugName::kALb)];
            bool l2 = util::get_rand_uniform() < this->kDl[static_cast<int>(DrugName::kAQ)];
            if (l1 && l2) {
                dg = DrugName::kNoDrug;
                return true;
            }
            if (!l1 && l2) {
                dg = DrugName::kALb;
                return true;
            }
            if (l1 && !l2) {
                dg = DrugName::kAQ;
                return true;
            }
            return false;
        } else { // only one drug to lose
            if (util::get_rand_uniform() < this->kDl[static_cast<int>(dg)]) {
                dg = this->kDlOutcome[static_cast<int>(dg)];
                return true;
            }
            return false;
        }


    }
}

}
#include <iostream>
#include <feature_elem.hpp>
#include "FeatureStreamer.hpp"


using namespace std;


/** Feature streamer constructed from the count matrix and a position file.
 * @param pos_path File containing the feature positions in the matrix
 * @param mat_path Count matrix per feature
 * @param k_len k-mer size
 **/
FeatureStreamer::FeatureStreamer(const string pos_path, const string mat_path, const size_t k_len, const size_t nb_smp) : pos_file(pos_path), mat_file(mat_path), k_len(k_len), nb_smp(nb_smp), already_loaded(false), from_file(false), merged_features(false)
{
    // File opening verification
    if (! this->pos_file.is_open()) {
        throw std::invalid_argument("loading index-pos failed, KaMRaT index folder not found or may be corrupted");
    }
    if (! this->mat_file.is_open()) {
        throw std::invalid_argument("loading index-mat failed, KaMRaT index folder not found or may be corrupted");
    }

    this->feature_pos.push_back(0);
}

/** Feature streamer constructed from a feature file.
 * @param pos_path File containing the features
 **/
FeatureStreamer::FeatureStreamer(const string file) : pos_file(file), mat_file(""), k_len(0), nb_smp(0), already_loaded(false), from_file(true), merged_features(false)
{
    // File opening verification
    if (! this->pos_file.is_open()) {
        throw std::invalid_argument("loading feature list failed, file not found or may be corrupted");
    }

    if (! this->hasNext()) {
        throw std::invalid_argument("not valid file for input sequences: " + file);
    }
}

FeatureStreamer::FeatureStreamer(FeatureStreamer&& fs) : k_len(fs.k_len), nb_smp(fs.nb_smp), already_loaded(fs.already_loaded), feature_pos(std::move(fs.feature_pos)), current_feature(fs.current_feature), from_file(fs.from_file), merged_features(fs.merged_features)

{
    std::swap(this->pos_file, fs.pos_file);
    std::swap(this->mat_file, fs.mat_file);
}

FeatureStreamer::~FeatureStreamer() {
    if (this->pos_file.is_open())
        this->pos_file.close();
    if (this->mat_file.is_open())
        this->mat_file.close();
}


/** Status of the streamer
 * @return True if the stream still contains features
 **/
bool FeatureStreamer::hasNext() {
    if (this->from_file) {
        return this->hasNext_file();
    } else {
        return this->hasNext_matrix();
    }
}

bool FeatureStreamer::hasNext_file() {
    if (this->already_loaded)
        return true;

    float _; // for "place-holder" column of feature scores
    size_t nb_mem_pos(1);

    // Read the file
    if (! (this->pos_file >> this->current_feature >> _ >> nb_mem_pos))
        return false;

    // Update flags
    if (nb_mem_pos > 0)
        this->merged_features = true;

    // Read the positions
    this->feature_pos.resize(nb_mem_pos);
    this->pos_file.ignore(1);
    this->pos_file.read(reinterpret_cast<char *>(&(this->feature_pos[0])), nb_mem_pos * sizeof(size_t));

    // Conclude
    this->already_loaded = true;
    return true;
}

bool FeatureStreamer::hasNext_matrix() {
    if (this->already_loaded)
        // Already in read
        return true;
    else if (this->pos_file.read(reinterpret_cast<char *>(&(this->feature_pos[0])), sizeof(uint64_t))) {
        // Increment the file position
        if (this->k_len > 0) {
            if (! this->pos_file.read(reinterpret_cast<char *>(&(this->feature_pos[0])), sizeof(size_t)))
                return false;
        }
        this->already_loaded = true;
        return true;
    }

    // End of file
    return false;
}


/** Return the next feature in the stream.
 * @return The next feature. nullptr if hasNext() == false
 **/
feature_t FeatureStreamer::next() {
    if (this->from_file)
        return this->next_file();
    else
        return this->next_matrix();
}

feature_t FeatureStreamer::next_file() {
    // No remaining feature
    if (!this->hasNext_file())
        return unique_ptr<FeatureElem>(nullptr);

    this->already_loaded = false;
    return std::make_unique<FeatureElem>(this->current_feature, this->feature_pos);
}

feature_t FeatureStreamer::next_matrix() {
    // No remaining k-mer
    if (! this->hasNext_matrix())
        return unique_ptr<FeatureElem>(nullptr);

    // Get the k-mer from the matrix
    size_t matrix_position = this->feature_pos[0] + this->nb_smp * sizeof(float); // skip the indexed count vector
    // Go to the right matrix position
    this->mat_file.seekg(matrix_position);
    
    string feature;
    this->mat_file >> feature;
    this->already_loaded = false;
    
    return std::make_unique<FeatureElem>(feature, this->feature_pos[0]);
}


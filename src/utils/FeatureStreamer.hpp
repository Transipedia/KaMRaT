#include <string>
#include <fstream>
#include <memory>
#include <vector>

#ifndef FEAT_STREM
#define FEAT_STREM

using feature_t = std::unique_ptr<FeatureElem>;


class FeatureStreamer {
private:
	std::ifstream pos_file;
	std::ifstream mat_file;
	size_t k_len;
	size_t nb_smp;

	size_t current_pos_pos;
	size_t current_mat_pos;

	bool already_loaded;
	std::vector<size_t> feature_pos;
	std::string current_feature;

	bool from_file;

	feature_t next_matrix();
	feature_t next_file();
	bool hasNext_matrix();
	bool hasNext_file();

	
public:
	bool merged_features;
	/** Open the files and position the file pointers on the right position.
	 * @param pos_path File containing the feature positions in the matrix
	 * @param mat_path Count matrix per feature
	 * @param k_len k-mer size
	 * @param nb_smp;
	 **/
	FeatureStreamer(const std::string pos_path, const std::string mat_path, const size_t k_len, const size_t nb_smp);
	FeatureStreamer(const std::string file);
	~FeatureStreamer();
	/** Status of the streamer
	 * @return True if the stream still contains features
	 **/
	bool hasNext();
	/** Return the next feature in the stream.
	 * @return The next feature. nullptr if hasNext() == false
	 **/
	feature_t next();
};



#endif

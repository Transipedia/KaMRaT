#include <string>
#include <fstream>
#include <memory>
#include <vector>

#ifndef IRA_HPP
#define IRA_HPP


class IndexRandomAccess {
private:
	std::ifstream pos_file;
	std::ifstream mat_file;

	size_t matrix_line_size;
	size_t pos_line_size;

	bool kmers;

	uint64_t mat_length;
	uint64_t mat_position;
	uint64_t pos_length;
	uint64_t pos_position;
	
public:
	size_t k;
	bool stranded;
	uint64_t nb_smp;
	uint64_t nb_rows;
	std::vector<std::string> matrix_header;
	
	/** Open the files and position the file pointers on the right position.
	 * @param pos_path File containing the feature positions in the matrix
	 * @param mat_path Count matrix per feature
	 * @param meta_path Metadata path
	 * @param kmers True if the features are kmers
	 **/
	IndexRandomAccess(const std::string pos_path, const std::string mat_path, const std::string meta_path, bool kmers=true);

	~IndexRandomAccess();

	/** Random access to the matrix. Extract the count vector and the feature.
	 * Restriction: This function only works for kmer matricies
	 * @param idx The position (line idx) to extract from the matrix
	 * @param counts A float array already allocated with at least nb_smp elements
	 * @param feature A char * containing at least k+1 bytes
	 **/
	void load_counts_by_row (const uint64_t idx, float * counts, char * feature);

	/** Random access to the matrix. Extract the count vector and the feature.
	 * @param file_position The position (in Bytes) in the matrix
	 * @param counts A float array already allocated with at least nb_smp elements
	 * @param feature A char * containing at least k+1 bytes
	 **/
	void load_counts_by_file_position (const uint64_t file_position, float * counts, char * feature);

	/** Random access to the position file to extract a matrix position
	 * @param feature_idx The index of the feature (0 based).
	 * @return Position of the feature in the matrix (in Bytes).
	 **/
	uint64_t feature_to_position(const uint64_t feature_idx);

	/** Random access to the matrix through the feature position in the pos index.
	 * Step 1, load the position from the pos file. Step 2 load the counts + feature from mat file.
	 * @param row Row index in the pos file. This address contains the matrix index where to 
	 * find the feature counts.
	 * @param counts A float array already allocated with at least nb_smp elements
	 * @param feature A char * containing at least k+1 bytes
	 **/
	void indirect_load_counts (const uint64_t row, float * counts, char * feature);
};



#endif

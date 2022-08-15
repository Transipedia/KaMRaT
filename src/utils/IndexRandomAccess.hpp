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

	uint64_t mat_length;
	uint64_t mat_position;
	uint64_t pos_length;
	
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
	 **/
	IndexRandomAccess(const std::string pos_path, const std::string mat_path, const std::string meta_path);

	~IndexRandomAccess();

	/** Random access to the matrix. Extract the count vector and the feature.
	 * @param idx The position (line idx) to extract from the matrix
	 * @param counts A float array already allocated with at least nb_smp elements
	 * @param feature A char * containing at least k+1 bytes
	 **/
	void load_counts_by_idx (const uint64_t idx, float * counts, char * feature);
};



#endif

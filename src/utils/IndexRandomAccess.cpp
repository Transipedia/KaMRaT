#include <iostream>

#include "IndexRandomAccess.hpp"


using namespace std;


IndexRandomAccess::IndexRandomAccess(const std::string pos_path, const std::string mat_path, const std::string meta_path) : pos_file(pos_path), mat_file(mat_path), stranded(false)
{
	// Opening
	ifstream meta_file(meta_path);
	if (!meta_file.is_open())
    {
        throw std::invalid_argument("loading index-meta failed, KaMRaT index folder not found or may be corrupted");
    }

    // Variables
    size_t nb_sample;
    std::string term;

    // Read the first line
    meta_file >> nb_sample >> this->k;
    if (this->k > 0)
    {
        meta_file >> term;
        this->stranded = (term == "T" ? true : false);
    }

    // Read the header
    meta_file >> term;
    for (size_t i(0) ; i < nb_sample ; ++i) // excluding the first column name
    {
        meta_file >> term;
        this->matrix_header.emplace_back(term);
    }

	// Cleaning
	meta_file.close();

	// Define usefull variables
	this->nb_smp = this->matrix_header.size();
	this->matrix_line_size = this->matrix_header.size() * sizeof(float) + this->k + 1;
	this->pos_line_size = sizeof(size_t) + this->k == 0 ? 0 : sizeof(uint64_t);
	
	this->mat_file.seekg (0, this->mat_file.end);
    this->mat_length = this->mat_file.tellg();
	this->mat_file.seekg (0);
	this->mat_position = 0;
	this->nb_rows = this->mat_length / (uint64_t)this->matrix_line_size;

	this->pos_file.seekg (0, this->pos_file.end);
    this->pos_length = this->pos_file.tellg();
	this->pos_file.seekg (0);
	this->pos_position = 0;
}


IndexRandomAccess::~IndexRandomAccess()
{
	if (this->pos_file.is_open())
		this->pos_file.close();
	if (this->mat_file.is_open())
		this->mat_file.close();
}



// ------------------------------------ Matrix Access ------------------------------------
#include <cassert>


void IndexRandomAccess::load_counts_by_file_position (const uint64_t file_position, float * counts, char * feature)
{
	// Extract counts
	uint64_t go_position = file_position;
	if (go_position != this->mat_position) {
		this->mat_file.seekg(go_position);
		this->mat_position = go_position;
	}
	this->mat_file.read(reinterpret_cast<char *>(counts), this->nb_smp * sizeof(float));
	
	// Extract feature
	this->mat_file.read(feature, this->k);
	feature[this->k] = '\0';
	
	// Read the extra \n byte
	char c;
	this->mat_file.read(&c, 1);
	assert(c == '\n');

	this->mat_position += this->nb_smp * sizeof(float) + this->k + 1;
}


void IndexRandomAccess::load_counts_by_row (const uint64_t row, float * counts, char * feature) 
{
	// Make sure that we are not reading random feature matrices
	if (this->k == 0)
		throw runtime_error("Impossible to have random access in non-kmer matrix");

	uint64_t go_position = row * this->matrix_line_size;
	this->load_counts_by_file_position(go_position, counts, feature);
}


uint64_t IndexRandomAccess::feature_to_position(const uint64_t feature_idx) {
	uint64_t go_position = feature_idx * this->pos_line_size;

	if (this->pos_position != go_position) {
		this->pos_file.seekg(go_position);
		this->pos_position = go_position;
	}

	uint64_t kmer_code;
	if (this->k > 0) {
		this->pos_file.read(reinterpret_cast<char *>(&kmer_code), sizeof(uint64_t));
	}
	size_t matrix_index;
	this->pos_file.read(reinterpret_cast<char *>(&matrix_index), sizeof(size_t));

	this->pos_position += sizeof(size_t) + this->k == 0 ? 0 : sizeof(uint64_t);

	return matrix_index;
}


void IndexRandomAccess::indirect_load_counts (const uint64_t row, float * counts, char * feature)
{
	size_t matrix_index = this->feature_to_position(row);
	this->load_counts_by_file_position(matrix_index, counts, feature);
}
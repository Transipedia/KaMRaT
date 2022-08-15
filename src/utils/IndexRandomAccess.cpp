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
	
	this->mat_file.seekg (0, this->mat_file.end);
    this->mat_length = this->mat_file.tellg();
	this->mat_file.seekg (0);
	this->mat_position = 0;
	this->nb_rows = this->mat_length / (uint64_t)this->matrix_line_size;

	this->pos_file.seekg (0, this->pos_file.end);
    this->pos_length = this->pos_file.tellg();
	this->pos_file.seekg (0);
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

void IndexRandomAccess::load_counts_by_idx (const uint64_t idx, float * counts, char * feature) 
{
	// Make sure that we are not reading random feature matrices
	if (this->k == 0)
		throw runtime_error("Impossible to have random access in non-kmer matrix");

	// Extract counts
	uint64_t go_position = idx * this->matrix_line_size;
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
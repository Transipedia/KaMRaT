target : kamrat kamratEvaluate clean

kamrat : kamrat.cpp kamratRank.o kamratNorm.o kamratMerge.o kamratMask.o kamratFilter.o
	g++ -o ../bin/kamrat kamrat.cpp *.o -fopenmp -lmlpack -larmadillo -I/home/haoliang.xue/.conda/envs/mlpack/include -lboost_iostreams

kamratRank.o : kamratRank.cpp data_struct/scorer.cpp data_struct/count_tab_header.cpp data_struct/count_tab_by_string.cpp data_struct/seq_elem.cpp
	g++ -c kamratRank.cpp data_struct/scorer.cpp data_struct/count_tab_header.cpp data_struct/count_tab_by_string.cpp data_struct/seq_elem.cpp -fopenmp -lmlpack -larmadillo -I/home/haoliang.xue/.conda/envs/mlpack/include

kamratNorm.o : kamratNorm.cpp data_struct/count_tab_header.cpp data_struct/sample_info.cpp 
	g++ -c kamratNorm.cpp data_struct/count_tab_header.cpp data_struct/sample_info.cpp -lboost_iostreams

kamratMerge.o : kamratMerge.cpp data_struct/count_tab_header.cpp data_struct/count_tab_by_fields.cpp data_struct/seq_elem.cpp data_struct/contig_elem.cpp data_struct/merge_knot.cpp
	g++ -c kamratMerge.cpp data_struct/count_tab_header.cpp data_struct/count_tab_by_fields.cpp data_struct/seq_elem.cpp data_struct/contig_elem.cpp data_struct/merge_knot.cpp

kamratEvaluate : kamratEvaluate.cpp data_struct/seq_elem.cpp data_struct/count_tab_header.cpp data_struct/count_tab_by_fields.cpp
	g++ -o ../bin/kamratEvaluate kamratEvaluate.cpp data_struct/seq_elem.cpp data_struct/count_tab_header.cpp data_struct/count_tab_by_fields.cpp

kamratMask.o : kamratMask.cpp data_struct/seq_elem.cpp
	g++ -c kamratMask.cpp data_struct/seq_elem.cpp

kamratFilter.o : kamratFilter.cpp data_struct/count_tab_header.cpp
	g++ -c kamratFilter.cpp data_struct/count_tab_header.cpp

clean :
	rm -f *.o

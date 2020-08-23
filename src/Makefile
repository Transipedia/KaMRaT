target: kamrat kamratReduce kamratEvaluate

kamrat: kamrat.cpp kamratReduce kamratRank kamratNorm kamratMerge kamratEvaluate kamratMask kamratFilter
	g++ -o kamrat.cpp ../bin/*.o

# kamratQuery: kamratQuery.cpp query/contig_sample_info.cpp
# 	g++ -o ../bin/kamratQuery kamratQuery.cpp query/contig_sample_info.cpp

kamratReduce: kamratReduce.cpp utils/sample_info.cpp reduce/model_info.cpp
	g++ -o ../bin/kamratReduce kamratReduce.cpp utils/sample_info.cpp reduce/model_info.cpp -fopenmp -lmlpack -larmadillo -I/home/haoliang.xue/.conda/envs/mlpack/include

kamratRank: kamratRank.cpp data_struct/scorer.cpp data_struct/count_tab_header.cpp data_struct/count_tab_by_string.cpp data_struct/seq_elem.cpp
	g++ -c ../bin/kamratRank.o kamratRank.cpp data_struct/scorer.cpp data_struct/count_tab_header.cpp data_struct/count_tab_by_string.cpp data_struct/seq_elem.cpp -fopenmp -lmlpack -larmadillo -I/home/haoliang.xue/.conda/envs/mlpack/include

kamratNorm: kamratNorm.cpp data_struct/count_tab_header.cpp data_struct/sample_info.cpp 
	g++ -c ../bin/kamratNorm.o kamratNorm.cpp data_struct/count_tab_header.cpp data_struct/sample_info.cpp -lboost_iostreams

kamratMerge: kamratMerge.cpp data_struct/count_tab_header.cpp data_struct/count_tab_by_fields.cpp data_struct/seq_elem.cpp data_struct/contig_elem.cpp data_struct/merge_knot.cpp
	g++ -c ../bin/kamratMerge.o kamratMerge.cpp data_struct/count_tab_header.cpp data_struct/count_tab_by_fields.cpp data_struct/seq_elem.cpp data_struct/contig_elem.cpp data_struct/merge_knot.cpp

kamratEvaluate: kamratEvaluate.cpp data_struct/seq_elem.cpp data_struct/count_tab_header.cpp data_struct/count_tab_by_fields.cpp
	g++ -o ../bin/kamratEvaluate kamratEvaluate.cpp data_struct/seq_elem.cpp data_struct/count_tab_header.cpp data_struct/count_tab_by_fields.cpp

kamratMask: kamratMask.cpp data_struct/seq_elem.cpp
	g++ -c ../bin/kamratMask.o kamratMask.cpp data_struct/seq_elem.cpp

kamratFilter: kamratFilter.cpp data_struct/count_tab_header.cpp
	g++ -c ../bin/kamratFilter.o kamratFilter.cpp data_struct/count_tab_header.cpp

rm -f ../bin/*.o

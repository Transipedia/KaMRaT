target: kamratReduce kamratRank kamratNorm kamratMerge kamratEvaluate kamratMask

# kamratQuery: kamratQuery.cpp query/contig_sample_info.cpp
# 	g++ -o ../bin/kamratQuery kamratQuery.cpp query/contig_sample_info.cpp

kamratReduce: kamratReduce.cpp utils/sample_info.cpp reduce/model_info.cpp
	g++ -o ../bin/kamratReduce kamratReduce.cpp utils/sample_info.cpp reduce/model_info.cpp -fopenmp -lmlpack -larmadillo -I/home/haoliang.xue/.conda/envs/mlpack/include

kamratRank: kamratRank.cpp data_struct/scorer.cpp data_struct/count_tab_header.cpp data_struct/count_tab.cpp data_struct/seq_elem.cpp
	g++ -o ../bin/kamratRank kamratRank.cpp data_struct/scorer.cpp data_struct/count_tab_header.cpp data_struct/count_tab.cpp data_struct/seq_elem.cpp -fopenmp -lmlpack -larmadillo -I/home/haoliang.xue/.conda/envs/mlpack/include

kamratNorm: kamratNorm.cpp data_struct/count_tab_header.cpp data_struct/sample_info.cpp 
	g++ -o ../bin/kamratNorm kamratNorm.cpp data_struct/count_tab_header.cpp data_struct/sample_info.cpp -lboost_iostreams

kamratMerge: kamratMerge.cpp data_struct/count_tab_header.cpp data_struct/count_tab.cpp data_struct/seq_elem.cpp data_struct/contig_elem.cpp data_struct/merge_knot.cpp
	g++ -o ../bin/kamratMerge kamratMerge.cpp data_struct/count_tab_header.cpp data_struct/count_tab.cpp data_struct/seq_elem.cpp data_struct/contig_elem.cpp data_struct/merge_knot.cpp

kamratEvaluate: kamratEvaluate.cpp data_struct/seq_elem.cpp data_struct/count_tab_header.cpp data_struct/count_tab.cpp
	g++ -o ../bin/kamratEvaluate kamratEvaluate.cpp data_struct/seq_elem.cpp data_struct/count_tab_header.cpp data_struct/count_tab.cpp

kamratMask: kamratMask.cpp data_struct/seq_elem.cpp
	g++ -o ../bin/kamratMask kamratMask.cpp data_struct/seq_elem.cpp

target: kamratMerge kamratReduce kamratEvaluate

# kamratQuery: kamratQuery.cpp query/contig_sample_info.cpp
# 	g++ -o ../bin/kamratQuery kamratQuery.cpp query/contig_sample_info.cpp

kamratReduce: kamratReduce.cpp utils/sample_info.cpp reduce/model_info.cpp
	g++ -o ../bin/kamratReduce kamratReduce.cpp utils/sample_info.cpp reduce/model_info.cpp -fopenmp -L/home/haoliang.xue/.conda/envs/mlpack/lib -lmlpack -larmadillo -I/home/haoliang.xue/.conda/envs/mlpack/include

# kamratNorm: kamratNorm.cpp
# 	g++ -o ../bin/kamratNorm kamratNorm.cpp -lboost_iostreams

kamratMerge: kamratMerge.cpp utils/column_info.cpp merge/contig_info.cpp merge/overlap_knot.cpp
	g++ -o ../bin/kamratMerge kamratMerge.cpp utils/column_info.cpp merge/contig_info.cpp merge/overlap_knot.cpp

kamratEvaluate: kamratEvaluate.cpp utils/column_info.cpp data_struct/seq_elem.cpp
	g++ -o ../bin/kamratEvaluate kamratEvaluate.cpp utils/column_info.cpp data_struct/seq_elem.cpp

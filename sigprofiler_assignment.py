import SigProfilerAssignment as spa
from SigProfilerAssignment import Analyzer as Analyze
import pandas as pd



def runsigassignment(matrix, output, sinature_database_input, exclude_signature ):
	

	Analyze.cosmic_fit(matrix, output, input_type="matrix", context_type="96",
					collapse_to_SBS96=True, cosmic_version=3.3, exome=True,
					genome_build="GRCh38", signature_database=sinature_database_input,
					exclude_signature_subgroups=exclude_signature, export_probabilities=False,
					export_probabilities_per_mutation=False,
					sample_reconstruction_plots='pdf', verbose=False)
	
	
if __name__=="__main__":
	###target assignment
	# path_to_vcf ='/data/project/TRIUMPH/4.analysis/GATK/Mutect2/vcfs/output/SBS/TRIUMPH.SBS96.exome'
	# data = pd.read_csv(path_to_vcf, sep = '\t', index_col=0) # you can put the path to your tab delimited file containing the mutational catalog matrix/table
	# # cosmic_target = '/data/project/TRIUMPH/tcga/download/COSMIC/COSMIC_target.txt'
	# # out = '/data/project/TRIUMPH/4.analysis/signature/sigprofilerassignment/TRIUMPH.SBS96.exome.target'
	# # runsigassignment(data, out, cosmic_target, None)
	# cosmic_target = '/data/project/TRIUMPH/tcga/download/COSMIC/COSMIC_intersect.txt'
	# out = '/data/project/TRIUMPH/4.analysis/signature/sigprofilerassignment/TRIUMPH.SBS96.exome.intersect'
	# runsigassignment(data, out, cosmic_target, None)
	# cosmic_target = '/data/project/TRIUMPH/tcga/download/COSMIC/COSMIC_union.txt'
	# out = '/data/project/TRIUMPH/4.analysis/signature/sigprofilerassignment/TRIUMPH.SBS96.exome.union'
	# runsigassignment(data, out, cosmic_target, None)
 
	
	# path_to_vcf ='/data/project/TRIUMPH/4.analysis/signature/TP53_vcf/output/SBS/TRIUMPH.SBS96.exome'
	# data = pd.read_csv(path_to_vcf, sep = '\t', index_col=0) # you can put the path to your tab delimited file containing the mutational catalog matrix/table
	# cosmic_target = '/data/project/TRIUMPH/tcga/download/COSMIC/COSMIC_union_plus.txt'
	# out = '/data/project/TRIUMPH/4.analysis/signature/sigprofilerassignment/TRIUMPH.SBS96.exome.TP53.union_plus'
	# runsigassignment(data, out, cosmic_target, None)
	# path_to_vcf ='/data/project/TRIUMPH/4.analysis/signature/CDKN2A_vcf/output/SBS/TRIUMPH.SBS96.exome'
	# data = pd.read_csv(path_to_vcf, sep = '\t', index_col=0) # you can put the path to your tab delimited file containing the mutational catalog matrix/table
	# cosmic_target = '/data/project/TRIUMPH/tcga/download/COSMIC/COSMIC_union_plus.txt'
	# out = '/data/project/TRIUMPH/4.analysis/signature/sigprofilerassignment/TRIUMPH.SBS96.exome.CDKN2A.union_plus'
	# runsigassignment(data, out, cosmic_target, None)
	# cosmic_target = '/data/project/TRIUMPH/tcga/download/COSMIC/COSMIC_intersect_plus.txt'
	# out = '/data/project/TRIUMPH/4.analysis/signature/sigprofilerassignment/TRIUMPH.SBS96.exome.intersect_plus'
	# runsigassignment(data, out, cosmic_target, None)
	# path_to_vcf ='/data/project/TRIUMPH/4.analysis/GATK/Mutect2/vcfs/output/SBS/TRIUMPH.SBS96.exome'
	# data = pd.read_csv(path_to_vcf, sep = '\t', index_col=0) # you can put the path to your tab delimited file containing the mutational catalog matrix/table
	# cosmic_target = '/data/project/TRIUMPH/tcga/download/COSMIC/COSMIC_union_plus.txt'
	# out = '/data/project/TRIUMPH/4.analysis/signature/sigprofilerassignment/TRIUMPH.SBS96.exome.union_plus'
	# runsigassignment(data, out, cosmic_target, None)
	##Before exome
	# cosmic_target = '/data/project/TRIUMPH/tcga/download/COSMIC/COSMIC_intersect.txt'
	# out = '/data/project/TRIUMPH/4.analysis/signature/sigprofilerassignment/TRIUMPH.SBS96.exome.intersect'
	# runsigassignment(data, out, cosmic_target, None)
	# cosmic_target = '/data/project/TRIUMPH/tcga/download/COSMIC/COSMIC_union.txt'
	# out = '/data/project/TRIUMPH/4.analysis/signature/sigprofilerassignment/TRIUMPH.SBS96.exome.union'
	# runsigassignment(data, out, cosmic_target, None)
	##240513
	path_to_vcf ='/data/project/TRIUMPH/4.analysis/GATK/Mutect2/vcfs/output/SBS/TRIUMPH.SBS96.region'
	data = pd.read_csv(path_to_vcf, sep = '\t', index_col=0) # you can put the path to your tab delimited file containing the mutational catalog matrix/table
	cosmic_target = '/data/project/TRIUMPH/tcga/download/COSMIC/COSMIC_union_plus.txt'
	out = '/data/project/TRIUMPH/4.analysis/signature/sigprofilerassignment/TRIUMPH.SBS96.region.union_plus'
	runsigassignment(data, out, cosmic_target, None)
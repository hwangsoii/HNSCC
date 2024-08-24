import SigProfilerAssignment as spa
from SigProfilerAssignment import Analyzer as Analyze
import pandas as pd



def runsigassignment(matrix, output, sinature_database_input, exclude_signature ):
	

	Analyze.cosmic_fit(matrix, output, input_type="matrix", context_type="ID",
					collapse_to_SBS96=False, cosmic_version=3.3, exome=True,
					genome_build="GRCh37", signature_database=sinature_database_input,
					exclude_signature_subgroups=exclude_signature, export_probabilities=False,
					export_probabilities_per_mutation=False,
					sample_reconstruction_plots='pdf', verbose=False)
	
	
if __name__=="__main__":


	###ID
	path_to_vcf ='/data/project/TRIUMPH/4.analysis/GATK/Mutect2/vcfs/hg19vcfs/output/ID/TRIUMPH.ID83.exome'
	data = pd.read_csv(path_to_vcf, sep = '\t', index_col=0) # you can put the path to your tab delimited file containing the mutational catalog matrix/table
	# cosmic_target = '/data/project/TRIUMPH/tcga/download/COSMIC/COSMIC_target.txt'
	# out = '/data/project/TRIUMPH/4.analysis/signature/sigprofilerassignment/TRIUMPH.SBS96.exome.target'
	# runsigassignment(data, out, cosmic_target, None)
	cosmic_target = '/data/project/TRIUMPH/download/COSMIC/COSMIC_v3.3_ID_GRCh37.txt'
	out = '/data/project/TRIUMPH/4.analysis/signature/sigprofilerassignment/TRIUMPH.ID96.exome'
	runsigassignment(data, out, cosmic_target, None)
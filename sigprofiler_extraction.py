from SigProfilerExtractor import sigpro as sig
def main_function():    
	# path_to_example_table ='/data/project/TRIUMPH/4.analysis/GATK/Mutect2/vcfs/output/ID/TRIUMPH.ID96.exome'
	# data = path_to_example_table # you can put the path to your tab delimited file containing the mutational catalog matrix/table
	# out = '/data/project/TRIUMPH/4.analysis/signature/sigprofiler/TRIUMPH.ID96.exome'
	# sig.sigProfilerExtractor("matrix", out, data, opportunity_genome="GRCh38",exome=True, minimum_signatures=4, maximum_signatures=5,cpu=8)
	# ##SNV  ###main 은 4, 5
	# path_to_example_table ='/data/project/TRIUMPH/4.analysis/GATK/Mutect2/vcfs/output/SBS/TRIUMPH.SBS96.exome'
	# data = path_to_example_table # you can put the path to your tab delimited file containing the mutational catalog matrix/table
	# out = '/data/project/TRIUMPH/4.analysis/signature/sigprofiler/TRIUMPH.SBS96.exome'
	# sig.sigProfilerExtractor("matrix", out, data, opportunity_genome="GRCh38",exome=True, minimum_signatures=4, maximum_signatures=7,cpu=10)
	# ###TP53
	# path_to_example_table ='/data/project/TRIUMPH/4.analysis/signature/TP53nega_vcf/output/SBS/TRIUMPH.SBS96.exome'
	# data = path_to_example_table # you can put the path to your tab delimited file containing the mutational catalog matrix/table
	# out = '/data/project/TRIUMPH/4.analysis/signature/sigprofiler/TRIUMPH.TP53_nega.SBS96.exome'
	# sig.sigProfilerExtractor("matrix", out, data, opportunity_genome="GRCh38",exome=True, minimum_signatures=4, maximum_signatures=7,cpu=10)
	# ##CDKN2A
	# path_to_example_table ='/data/project/TRIUMPH/4.analysis/signature/CDKN2Anega_vcf/output/SBS/TRIUMPH.SBS96.exome'
	# data = path_to_example_table # you can put the path to your tab delimited file containing the mutational catalog matrix/table
	# out = '/data/project/TRIUMPH/4.analysis/signature/sigprofiler/TRIUMPH.CDKN2A_nega.SBS96.exome'
	# sig.sigProfilerExtractor("matrix", out, data, opportunity_genome="GRCh38",exome=True, minimum_signatures=4, maximum_signatures=7,cpu=10)
	# ###TP53
	# path_to_example_table ='/data/project/TRIUMPH/4.analysis/signature/TP53_vcf/output/SBS/TRIUMPH.SBS96.exome'
	# data = path_to_example_table # you can put the path to your tab delimited file containing the mutational catalog matrix/table
	# out = '/data/project/TRIUMPH/4.analysis/signature/sigprofiler/TRIUMPH.TP53.SBS96.exome'
	# sig.sigProfilerExtractor("matrix", out, data, opportunity_genome="GRCh38",exome=True, minimum_signatures=4, maximum_signatures=7,cpu=10)
	# ##CDKN2A
	# path_to_example_table ='/data/project/TRIUMPH/4.analysis/signature/CDKN2A_vcf/output/SBS/TRIUMPH.SBS96.exome'
	# data = path_to_example_table # you can put the path to your tab delimited file containing the mutational catalog matrix/table
	# out = '/data/project/TRIUMPH/4.analysis/signature/sigprofiler/TRIUMPH.CDKN2A.SBS96.exome'
	# sig.sigProfilerExtractor("matrix", out, data, opportunity_genome="GRCh38",exome=True, minimum_signatures=4, maximum_signatures=7,cpu=10)
	# ###HPV
	# path_to_example_table ='/data/project/TRIUMPH/4.analysis/signature/HPV_vcf/output/SBS/TRIUMPH.SBS96.exome'
	# data = path_to_example_table # you can put the path to your tab delimited file containing the mutational catalog matrix/table
	# out = '/data/project/TRIUMPH/4.analysis/signature/sigprofiler/TRIUMPH.HPV.SBS96.exome'
	# sig.sigProfilerExtractor("matrix", out, data, opportunity_genome="GRCh38",exome=True, minimum_signatures=4, maximum_signatures=7,cpu=10)
	# ##HPVnega
	# path_to_example_table ='/data/project/TRIUMPH/4.analysis/signature/HPVnega_vcf/output/SBS/TRIUMPH.SBS96.exome'
	# data = path_to_example_table # you can put the path to your tab delimited file containing the mutational catalog matrix/table
	# out = '/data/project/TRIUMPH/4.analysis/signature/sigprofiler/TRIUMPH.HPVnega.SBS96.exome'
	# sig.sigProfilerExtractor("matrix", out, data, opportunity_genome="GRCh38",exome=True, minimum_signatures=4, maximum_signatures=7,cpu=10)

	##24.05.12
	path_to_example_table ='/data/project/TRIUMPH/4.analysis/GATK/Mutect2/vcfs_target/output/SBS/TRIUMPH.SBS96.region'
	data = path_to_example_table # you can put the path to your tab delimited file containing the mutational catalog matrix/table
	out = '/data/project/TRIUMPH/4.analysis/signature/sigprofiler/TRIUMPH.SBS96.region'
	sig.sigProfilerExtractor("matrix", out, data, opportunity_genome="GRCh38",exome=False, minimum_signatures=4, maximum_signatures=7,cpu=10)	
if __name__=="__main__":
	main_function()
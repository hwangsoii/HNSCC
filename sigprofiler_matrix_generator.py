from SigProfilerMatrixGenerator import install as genInstall

# genInstall.install('GRCh37', rsync=False, bash=True)
# genInstall.install('GRCh38', rsync=False, bash=True)


from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen


INTERVAL='/data/project/TRIUMPH/bed/1612AHP-0021_KimHR_3033241_Cho_MutScape_V2_1_Regions.liftover.sorted.bed'
# matrices = matGen.SigProfilerMatrixGeneratorFunc("TRIUMPH", "GRCh38", "/data/project/TRIUMPH/4.analysis/GATK/Mutect2/vcfs",plot=True, exome=True, bed_file=INTERVAL, chrom_based=False, tsb_stat=False, seqInfo=False, cushion=100)

##24.05.12_for vcfs downsample to target
matrices = matGen.SigProfilerMatrixGeneratorFunc("TRIUMPH", "GRCh38", "/data/project/TRIUMPH/4.analysis/GATK/Mutect2/vcfs_target",plot=True, exome=False, bed_file=INTERVAL, chrom_based=False, tsb_stat=False, seqInfo=False, cushion=100)
# matrices = matGen.SigProfilerMatrixGeneratorFunc("TRIUMPH", "GRCh38", "/data/project/TRIUMPH/4.analysis/GATK/Mutect2/vcfs/filtered_vcfs",plot=True, exome=True, bed_file=INTERVAL, chrom_based=False, tsb_stat=True, seqInfo=False, cushion=100)
# matrices = matGen.SigProfilerMatrixGeneratorFunc("TRIUMPH", "GRCh37", "/data/project/TRIUMPH/4.analysis/GATK/Mutect2/vcfs/hg19vcfs",plot=True, exome=True, bed_file=INTERVAL, chrom_based=False, tsb_stat=True, seqInfo=False, cushion=100)

# matrices = matGen.SigProfilerMatrixGeneratorFunc("TRIUMPH", "GRCh38", "/data/project/TRIUMPH/4.analysis/signature/TP53_vcf",plot=True, exome=True, bed_file=INTERVAL, chrom_based=False, tsb_stat=False, seqInfo=False, cushion=100)
# matrices = matGen.SigProfilerMatrixGeneratorFunc("TRIUMPH", "GRCh38", "/data/project/TRIUMPH/4.analysis/signature/CDKN2A_vcf",plot=True, exome=True, bed_file=INTERVAL, chrom_based=False, tsb_stat=False, seqInfo=False, cushion=100)


# matrices = matGen.SigProfilerMatrixGeneratorFunc("TRIUMPH", "GRCh38", "/data/project/TRIUMPH/4.analysis/signature/TP53nega_vcf",plot=True, exome=True, bed_file=INTERVAL, chrom_based=False, tsb_stat=False, seqInfo=False, cushion=100)
# matrices = matGen.SigProfilerMatrixGeneratorFunc("TRIUMPH", "GRCh38", "/data/project/TRIUMPH/4.analysis/signature/CDKN2Anega_vcf",plot=True, exome=True, bed_file=INTERVAL, chrom_based=False, tsb_stat=False, seqInfo=False, cushion=100)
# matrices = matGen.SigProfilerMatrixGeneratorFunc("TRIUMPH", "GRCh38", "/data/project/TRIUMPH/4.analysis/signature/HPV_vcf",plot=True, exome=False, bed_file=INTERVAL, chrom_based=False, tsb_stat=False, seqInfo=False, cushion=100)
# matrices = matGen.SigProfilerMatrixGeneratorFunc("TRIUMPH", "GRCh38", "/data/project/TRIUMPH/4.analysis/signature/HPVnega_vcf",plot=True, exome=False, bed_file=INTERVAL, chrom_based=False, tsb_stat=False, seqInfo=False, cushion=100)
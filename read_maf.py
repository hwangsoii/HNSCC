import sys
import os
import numpy as np
from scipy.stats import chi2_contingency
from scipy.stats import fisher_exact
from collections import Counter

def makedirs(path): 
	try: 
		os.makedirs(path) 
	except OSError: 
	   if not os.path.isdir(path): 
		   raise

black_gene = '/data/project/TRIUMPH/download/blacklist/blacklist_gene.txt'

black_dic = {}
with open(black_gene, 'r') as fb:
	for line in fb:
		# print(line.rstrip())
		gene = line.rstrip()
		if not gene in black_dic:
			black_dic[gene] = 1

# key_list = list(black_dic.keys())
# print(key_list)
dis = '/data/project/TRIUMPH'
datapath= '/data/project/TRIUMPH/4.analysis/maf/all'
datapath= '/data/project/TRIUMPH/4.analysis/vcf2maf/germline/integ_all'
datapath= '/data/project/TRIUMPH/4.analysis/vcf2maf/germline/integ_all_hard_freq'
# filename = sys.argv[4]
subpath = '/data/project/TRIUMPH/4.analysis/maf/germline/SNV/BL'
makedirs(subpath)
# sample = sys.argv[5]
mutation_li = ["Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation", "Nonsense_Mutation", "Silent", "Splice_Site", "Translation_Start_Site", "Nonstop_Mutation", "3'UTR", "3'Flank", "5'UTR", "5'Flank", "IGR", "Intron", "RNA", "Targeted_Region", "De_novo_Start_InFrame", "De_novo_Start_OutOfFrame", "Splice_Region", "Unknown"]
file_list= os.listdir(datapath)
p53_t_pik_t =0 
p53_t_pik_f =0 
p53_f_pik_t =0 
p53_f_pik_f =0 
p53_num =0
pik3_num = 0
p53_pik3 =0
pik3_final_li = []
pik3_p53_li = []
for file_name in file_list:
	# if file_name.endswith('.filtered_snps.vcf'):
	###### vep.maf
 	# sample = file_name[:-12]
	####black.maf
	sample = file_name[:-18]
	# if file_name.endswith('.black.maf'):
	# 	print(sample)
	# if file_name.endswith('.vep.maf'):
	if file_name.endswith('.black.hard.maf'):
	# if file_name.endswith('.black.soft.maf'):
		p53 = False
		pik3 = False
		pik3_li = []
		with open('{}/{}'.format(datapath, file_name), 'r') as fr:
			for line in fr:
				if line.startswith('#'):
					# fw.write(line)
					continue
				
				fact = line.rstrip().split('\t')
				gene = fact[0]
				assess = fact[8]
				IMPACT = fact[92]
				clinsig_dbsnp= fact[85]
				POLYphen = fact[72]
				SIFT = fact[71]
				gnomad = fact[122]
				af = fact[76]
				vcf_pos = fact[-2]
				vcf = fact[131]
				dbSNP = fact[14]
				# if gene == 'TP53':
				# # if gene == 'KMT2C':
				# 	#print('TP53', assess)
				# 	if assess in mutation_li:
				# 		p53_num += 1
				# 		p53 = True
				# 		break
				# 	else:
				# 		print(assess)
				# 		# print(SIFT, POLYphen, SIFT, clinsig_dbsnp)
				if gene == 'FAT3':
					if assess in mutation_li:
						pik3 = True
						HGVSp_short = fact[36]
						pik3_li.append(HGVSp_short)
	  					# p53_num += 1
						# pik3_num +=1
					HGVSp_short = fact[36]
					print(HGVSp_short)
					# break
				# print(fact[42], fact[43], fact[44])
				
				# print('vcf_pos', vcf_pos)
				# print('vcf', vcf)
				# print(assess)
				# if not '' == clinsig_dbsnp:
				# if 'likely pathogenic' in clinsig_dbsnp or 'pathogenic' in clinsig_dbsnp:
				# 	# print('clinsig', clinsig_dbsnp)
				# if 'damaging' in POLYphen:
				# 	print('polyphen', POLYphen)
				# print(POLYphen)
				#     continue
				# if 'deleterious' in SIFT:
				# 	print('sift', SIFT)
				# if 'MODERATE' in IMPACT or 'HIGH' in IMPACT:
				# 	print(IMPACT)
				#     continue
				# elif 'LOW' in IMPACT:
				#     continue
				
				# print('polyphen', POLYphen)
				# print('sift', SIFT)
				# print('impact', IMPACT)
				# if not gene in black_dic:
				#     fw.write(line)

	
				# if IMPACT == 'high':
				# 	print(line)
				# print(assess)
				# if not gnomad == '':	
				# 	print('gnomad, {}'.format(gnomad))
		if len(pik3_li) >= 1:
			for pik in pik3_li:
				if pik =='':
					continue
				pik3_final_li.append(pik)
		# if len(pik3_li) >= 1:
		# 	for pik in pik3_li:
		# 		pik3_final_li.append(pik)
		# if p53 == True and pik3 == True:
		# 	p53_t_pik_t +=1 
		# elif p53 == True and pik3!= True:
		# 	p53_t_pik_f +=1
		# elif p53 != True and pik3== True:
		# 	p53_f_pik_t +=1
		# 	if len(pik3_li) >= 1:
		# 		for pik in pik3_li:
		# 			pik3_final_li.append(pik)
		# 	# if len(pik3_li) == 1:
		# 	# 	pik3_final_li.append(pik3_li[0])
		# 	# elif len(pik3_li) >= 2:
		# 	# 	pik3_final_li.append(pik3_li[0])
		# 	# 	print('multiple', pik3_li)
		# elif p53 != True and pik3!= True:
		# 	p53_f_pik_f +=1

# print('p53_t_pik_t', p53_t_pik_t, 'p53_t_pik_f', p53_t_pik_f, 'p53_f_pik_t', p53_f_pik_t, 'p53_f_pik_f', p53_f_pik_f )


# observed = np.array([[p53_t_pik_t, p53_t_pik_f], [p53_f_pik_t, p53_f_pik_f]])
# # chi_val, p_val, dof, expected =  chi2_contingency(observed)
# # print(chi_val, p_val, dof, expected)
# oddsr, p = fisher_exact(observed, alternative='two-sided')
# print('pval', p)
print(pik3_final_li)
print(len(pik3_final_li))
counter = Counter(pik3_final_li)
# for item in pik3_final_li:
#     counter[item] +=1
print(counter)
print(counter.most_common(n=5))
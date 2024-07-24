###############################################################################################
# Code adapted from source provided by Paras Garg [@AndyMSSMLab] 
###############################################################################################


### requires file with list of top 100 SNPs, must be 1 column (suffix.snplist)
###		variantID
### requires file with list of TR to be conditioned, must be 1 column
###		variantID 
function getTopSNPsandTRGeno {
    bcftools view  -Ou --include ID==@${suffix}.snplist $SNP_VCF |  
        bcftools query -f'[%ID\t%SAMPLE\t%GT\n]' | 
        sed -e 's/0\/0/0/' -e 's/0\/1/1/' -e 's/1\/0/1/' -e 's/1\/1/2/' -e  's/.\/./NA/' > ${suffix}.geno.txt
	
	bcftools view  -Ou --include ID==@${suffix}.selectTR $TR_VCF -S ${PREFIX}.selectSamples --force-samples| 
        bcftools query -f'[%ID\t%SAMPLE\t%DS\n]' | sed -e  's/.\/./NA/'  >> ${suffix}.geno.txt
}


function runCAVIAR {
	### Requires table of top 100 SNPs from METAL reference allele dosage genotypes and TR dosage genotype in 3 columns
	### 	variantID, SampleID, Genotype
	### Requires LD matrix. In R: cor = round(cor(G, use = "pairwise.complete.obs"),6)
	### 	101 x 101 square matrix of top 100 SNPs and TR
	### Requires table of Z-scores generated from REGENIE summary statistics. In R: Z = sign(BETA) * sqrt(CHISQ)
	###		For top 100 SNPs and TR, in 2 columns
	### 	variantID, Z
	
	CAVIAR -l ${suffix}.LD -z ${suffix}.Zscores -o ${suffix}.caviar -r 0.95 -g 0.01 -c 2 -f1
}

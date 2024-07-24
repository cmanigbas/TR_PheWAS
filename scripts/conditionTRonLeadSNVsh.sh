####################################################################################
# Code adapted from source provided by Paras Garg 
####################################################################################



### requires SNPs +/- 500 kb of TR  plink .bed .bim .fam files
### 	basic filtering for AF > 0.001 and biallelic variants,
### 	mac > 100, missing call rate > 0.05, hardy weinberg > 1e-300
### require SNPs VCF file 
### requires phenotype table 
###		FID, IID, traitA .., traitZ
### require covariate file
###		FID, IID, covarA ..., covarZ
### require SNPs VCF file 
### requires TR dosage psudeogenotype VCF file

VarIn=$1				# prefix for flanking SNPS plink files (without .fileExtension)
PhenoIn=$2				# phenotype table
CovarIn=$3				# table of covariates
Phenoname=$4				# phenotype ID
suffix=$5				# file suffix				
mode=$6					# binary or quant trait type

if [[ $mode == "Binary "]]; 
	command=" --bt --spa "
elif [[ $mode == "Quant" ]];
	command=" --qt --apply-rint "
fi

### FIND LEAD SNP +/- 500 bp of TR
function SNPregenie_step1 {
    regenie \
        --step 1 \
        --bed ${VarIn} \
        $command \
		--phenoFile ${PhenoIn} \
        --covarFile ${CovarIn} \
        --phenoCol ${Phenoname} \
        --bsize 100 \
        --lowmem \
        --lowmem-prefix tmp_rg.${suffix} 
        --loocv \
        --threads 18 \
        --catCovarList $cmdCatVar \
        --covarColList $cmdCovar \
        --out ${SUFFIX}_Reg_S1 \
}
function SNPregenie_step2 {
    regenie \
        --step 2 \
        --bed ${VarIn} \
        $command \
		--phenoFile ${PhenoIn} \
		--phenoCol ${Phenoname} \
        --covarFile ${CovarIn} \
        --bsize 100 \
        --loocv \
        --threads 18 \
		--pThresh 0.05 ${command} 
		--pred ${suffix}_Reg_S1
		--minMAC 5 \
        --out ${suffix}_Reg_S2 \
}
function runMetal {
	### requires regenie results for sub-cohorts being meta-analyzed
	SCRIPT_FILE=${suffix}.script
	METAL_OUTFILE=${suffix}.metal
	METAL_LOGFILE=${suffix}.metal.log


	rm -f $SCRIPT_FILE
	for regenie in `ls ${suffix}.*.regenie`; do 
		echo MARKER UNIQUE_ID >> $SCRIPT_FILE
		echo WEIGHT N >> $SCRIPT_FILE
		echo ALLELE ALLELE0 ALLELE1 >> $SCRIPT_FILE
		echo EFFECT BETA >> $SCRIPT_FILE
		echo PVAL P  >> $SCRIPT_FILE

		echo PROCESS ${regenie} >> $SCRIPT_FILE
		echo >> $SCRIPT_FILE
		echo >> $SCRIPT_FILE

	done 
	echo OUTFILE $METAL_OUTFILE .txt >> $SCRIPT_FILE
	echo ANALYZE >> $SCRIPT_FILE

	ml metal/2018-08-28
	metal < $SCRIPT_FILE > $METAL_LOGFILE 2> $METAL_LOGFILE
}


### Requires .txt with variant ID of lead SNV (suffix.leadSNP)
###		variantID 
### Require SNP vcf file
### $COHORT is the sequencing center sub-cohort that is most representative of all the samples 
function divideTRintoLeadSNVgenotype {
	bcftools view  -Ou --include ID==@${suffix}.LeadSNP $SNP_VCF |
        bcftools query -f'[%ID\t%SAMPLE\t%GT\n]' | 
        sed -e 's/0\/0/AA/' -e 's/0\/1/AB/' -e 's/1\/0/AB/' -e 's/1\/1/BB/' -e  's/.\/./NA/' > ${suffix}.LeadSNP.geno.txt
   
    cat ${suffix}.LeadSNP.geno.txt | awk -v cohort=$COHORT '{print $0,cohort}' OFS="\t" > ${suffix}.LeadSNP.geno.txt.temp
    mv ${suffix}.LeadSNP.geno.txt.temp ${suffix}.LeadSNP.geno.txt
    
    cat ${suffix}.LeadSNP.geno.txt | awk '$3=="AA"{print $2}' > ${suffix}.samples.AA
    cat ${suffix}.LeadSNP.geno.txt | awk '$3=="AB"{print $2}' > ${suffix}.samples.AB
    cat ${suffix}.LeadSNP.geno.txt | awk '$3=="BB"{print $2}' > ${suffix}.samples.BB
    
    for sample in `ls ${suffix}.samples.* | grep -e AA$ -e AB$ -e BB$`; do 
        echo -e "\t ... extracting $TR in $sample"
        echo $TR > ${PREFIX}.selectTR
        bcftools view  -Oz --include ID==@${PREFIX}.selectTR $TR_VCF -S $sample --force-samples > ${sample}.TR.vcf.gz
        plink2 --vcf ${sample}.TR.vcf.gz  'dosage=DS' --out ${sample}.TR  --double-id
    done
}	


### requires TR genotypes in plink .pgen .psam .pvar format, after splitting samples into lead SNV  dosages 
### plink2 --vcf ${suffix}.AA.TR.vcf.gz  'dosage=DS' --out ${sample}.TR  --double-id
###   		${suffix}.AA.TR.pgen, ${suffix}.AA.TR.pvar, ${suffix}.AA.TR.psam 
### 		${suffix}.AB.TR.pgen, ${suffix}.AB.TR.pvar, ${suffix}.AB.TR.psam 
###			${suffix}.BB.TR.pgen, ${suffix}.BB.TR.pvar, ${suffix}.BB.TR.psam 
function TRconditioned_regenie_step1 {
    regenie \
        --step 1 \
        --bed ${suffix}.AA.TR \
        $command \
		--phenoFile ${PhenoIn} \
		--phenoCol ${Phenoname} \
        --covarFile ${CovarIn} \
        --bsize 100 \
        --lowmem \
        --lowmem-prefix tmp_rg.${suffix}.AA.TR
        --loocv \
        --threads 18 \
        --catCovarList $cmdCatVar \
        --covarColList $cmdCovar \
        --out ${suffix}.AA.TR_Reg_S1 \
}
function TRconditioned_regenie_step2 {
    regenie \
        --step 2 \
        --bed ${suffix}.AA.TR \
        $command \
		--phenoFile ${PhenoIn} \
		--phenoCol ${Phenoname} \
        --covarFile ${CovarIn} \
        --bsize 100 \
        --loocv \
        --threads 18 \
		--pThresh 0.05 ${command} 
		--pred ${suffix}.AA.TR_Reg_S1
		--minMAC 5 \
        --out ${suffix}.AA.TR_Reg_S2 \
}
function runMetalTR_conditioned {
	SCRIPT_FILE=${suffix}.TRConditioned.script
	METAL_OUTFILE=${suffix}.TRConditioned.metal
	METAL_LOGFILE=${suffix}.TRConditioned.metal.log


	rm -f $SCRIPT_FILE
	for regenie in `ls ${suffix}.*.TR.regenie`; do 
		echo MARKER UNIQUE_ID >> $SCRIPT_FILE
		echo WEIGHT N >> $SCRIPT_FILE
		echo ALLELE ALLELE0 ALLELE1 >> $SCRIPT_FILE
		echo EFFECT BETA >> $SCRIPT_FILE
		echo PVAL P  >> $SCRIPT_FILE

		echo PROCESS ${regenie} >> $SCRIPT_FILE
		echo >> $SCRIPT_FILE
		echo >> $SCRIPT_FILE

	done 
	echo OUTFILE $METAL_OUTFILE .txt >> $SCRIPT_FILE
	echo ANALYZE >> $SCRIPT_FILE

	ml metal/2018-08-28
	metal < $SCRIPT_FILE > $METAL_LOGFILE 2> $METAL_LOGFILE
}










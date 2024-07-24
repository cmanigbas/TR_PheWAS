VarIn=$1
PhenoIn=$2
CovarIn=$3
Phenoname=$4
suffix=$5



regenie --step 1 \
		--pgen ${VarIn} \
		--qt \
		--apply-rint \
		--phenoFile ${PhenoIn} \
		--covarFile ${CovarIn} \
		--phenoCol ${Phenoname} \
		--bsize 100 \
		--loocv \
		--lowmem \
		--lowmem-prefix tmp_rg_${suffix} \
		--threads 18
		--out ${suffix}.step1 \		



regenie --step 2 \
		--pgen ${VarIn} \
		--qt \
		--apply-rint \
		--phenoFile ${PhenoIn} \
		--covarFile ${CovarIn} \
		--phenoCol ${Phenoname} \
		--bsize 100 \
		--loocv \
		--lowmem \
		--bsize 100 
		--lowmem \
		--lowmem-prefix tmp_rg2_${suffix} \
		--threads 18 \
		--pThresh 0.05 \
		--pred ${suffix}.step1.list
		--minMAC 5		
		
		

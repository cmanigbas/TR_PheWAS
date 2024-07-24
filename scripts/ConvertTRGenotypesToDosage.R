####################################################################################
# Code adapted from source provided by Paras Garg and Alejandro Martin Trujillo
####################################################################################



library(tidyverse)
library(scales)
library('parallel')


function convertGenotypeToDosage (region) {
	d = df %>% filter(VARID == region) %>% as.data.frame
	if (n_distinct (d$avg_repeats[!is.na(d$avg_repeats)]) > 1) { 
  		y = d %>%
			mutate(z = (avg_repeats - mean(avg_repeats, na.rm = T))/sd(avg_repeats, na.rm = T)) %>%
			mutate(r = round(rescale(z, to =c(0,1)),4)) %>%
			as.data.frame %>%
			mutate(GT = ifelse(r<0.5, "0/0", "0/1")) %>%
			mutate(GT = paste0(GT,":",r)) %>%
			mutate(end = paste0("END=",end), FORMAT = "GT:DS") %>%
			rename(CHROM = chr, POS = start,ID = VARID, INFO = end) %>%
			select(CHROM, POS, ID, INFO, FORMAT, SampleId, GT) %>%
			spread(SampleId, GT) %>%
			mutate(REF=2, ALT = 1, QUAL = ".", FILTER = ".") %>%
			select(CHROM:ID, REF, ALT, QUAL, FILTER, INFO, everything())

		rm (d)
		y[y == "NA:NA"] = "./.:."
		y = y %>% mutate(CHROM = as.numeric(sub("chr","",CHROM))) %>% arrange(CHROM, POS)
  
		colnames(y)[1] = "#CHROM" 
		y
  	}
}



### df data object is a data.frame and must contain 3 columns
### 	variant ID column (VARID), sample id column, and average allele size
### coord data object is a data.frame object containing variant coordinates
### 	chromosome, start, end, variantID
df <- reshape2::dcast (df,VARID ~ SampleId)
df <- reshape2::melt (df,id.vars = "VARID",value.name = "avg_repeats", variable.name = "SampleId") 
df <- coord %>% inner_join (df, by = "VARID")	

inputRegions = unique (df$VARID)
out = mclapply(1:length(inputRegions), convertGenotypeToDosage(i), mc.cores=10)


out = do.call(rbind, out)
cat (paste0("Dimensions of VCF \n \t rows: ", dim(out)[1], ", columns: ", dim(out)[2], "\n"))


header1 = '##fileformat=VCFv4.1'
header2='##INFO=<ID=END,Number=0,Type=Integer,Description="End">'
header3='##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">'
header4 = '##FORMAT=<ID=DS,Number=1,Type=Float,Description="dosage">'

write.table(header1, SAVENAME, row.names = F, sep = "\t", quote = F, col.names = F)
write.table(header2, SAVENAME, row.names = F, sep = "\t", quote = F, col.names = F, append = T)
write.table(header3, SAVENAME, row.names = F, sep = "\t", quote = F, col.names = F, append = T)
write.table(header4, SAVENAME, row.names = F, sep = "\t", quote = F, col.names = F, append = T)
write.table(out, SAVENAME, row.names = F, sep = "\t", quote = F, append = T)


### convert to PLINK
system(paste("plink2 --vcf", SAVENAME, ".vcf 'dosage=DS' --double-id --out", INPUT, "_plink"))

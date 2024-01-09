library(data.table)
library(openxlsx)


df1 = fread('./sexde/data/raw/GTEx_v8/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt')
df1 = df1[df1$SMAFRZE == "RNASEQ", c('SMTSD', 'SAMPID', 'SMTSISCH', 'SMRIN')]
df1$SMTSD = gsub(' - ', '_', df1$SMTSD)
df1$SMTSD = gsub('[()]', '', df1$SMTSD)
df1$SMTSD = gsub(' ', '_', df1$SMTSD)
df1$SUBJID = gsub("^(.*?)-(.*?)-.*", "\\1-\\2", df1$SAMPID)

df2 = fread('./sexde/data/raw/GTEx_v8/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt')
df2 = df2[, c('SUBJID', 'AGE', 'SEX')]
df2 = unique(df2[df2$SUBJID %in% df1$SUBJID,])
table(df2$SEX)
df2$SEX = ifelse(df2$SEX==1, 'M', 'F')
table(df2$SEX)

df1$AGE = df2$AGE[match(df1$SUBJID, df2$SUBJID)]
df1$SEX = df2$SEX[match(df1$SUBJID, df2$SUBJID)]

c = as.data.frame(table(df1$SMTSD, df1$SEX))
c1 = c[c$Var2=='M',]
c2 = c[c$Var2=='F',]
cnt = as.data.frame(table(df1$SMTSD))
colnames(cnt) = c('Tissue', 'N')
cnt$N_male = c1$Freq[match(cnt$Tissue, c1$Var1)]
cnt$N_female = c2$Freq[match(cnt$Tissue, c2$Var1)]
cnt = cnt[,c('Tissue', 'N_male', 'N_female', 'N')]
write.xlsx(cnt, './sexde/data/count_bytissue.xlsx', colNames=T, rowNames=F, quote=F, overwrite=T)

tissues = cnt$Tissue[cnt$N_male>0 & cnt$N_female>0] # without these conditions, we cannot estimate sexde for lab
write.table(tissues, './sexde/data/queue/tissue.txt', quote=F, col.names=F, row.names=F) # then delete some missing in HYFA

for (t in tissues) {
    print(t)
    df = df1[df1$SMTSD==t, c('SUBJID', 'SEX', 'SMTSISCH', 'SMRIN', 'AGE')]
    fh = paste0('./sexde/zenodo_sexde_pipeline/data/sample_covariates_updated/covs.basalcovs.',t,'.v8.txt')
    print(sum(is.na(df)))
    fwrite(df, fh, col.names=T, row.names=F, sep='\t', quote=F, na=NA)
}

imp <- aggregate(cbind(SMTSISCH, SMRIN) ~ SUBJID, data = df1, FUN = function(x) median(x, na.rm = T))
df2$SMTSISCH = imp$SMTSISCH[match(df2$SUBJID, imp$SUBJID)]
df2$SMRIN = imp$SMRIN[match(df2$SUBJID, imp$SUBJID)]
df2 = df2[,c('SUBJID', 'SEX', 'SMTSISCH', 'SMRIN', 'AGE')]
print(sum(is.na(df2)))
fwrite(df2, './sexde/zenodo_sexde_pipeline/data/sample_covariates_updated/covs.basalcovs.imp_all.v8.txt', col.names=T, row.names=F, sep='\t', quote=F, na=NA)
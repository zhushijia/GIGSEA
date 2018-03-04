tools::add_datalist("/Users/shijiazhu/Dropbox/mypaper/GIGSEA/codes/GIGSEA", force = T)


library(Matrix)
data(MSigDB.TF)
net1 = geneSet2Net( MSigDB.TF$term , MSigDB.TF$geneset , value=NULL , sep=',' )
net2 = geneSet2sparseMatrix( MSigDB.TF$term , MSigDB.TF$geneset , value=NULL , sep=',' )
object.size(net1)/1e6
object.size(net2)/1e6


net=as.matrix(net2)
net = net[,colSums(net)>=20]
imputeFC = heart
imputeFC = subset(imputeFC, !is.na(zscore) )
imputeFC = imputeFC[ order(as.character(imputeFC$gene_name)) , ]
net = net[ rownames(net) %in% as.character(imputeFC$gene_name) , ]
imputeFC = subset( imputeFC , gene_name %in% rownames(net) )

fc = imputeFC$zscore
usedFrac = imputeFC$n_snps_used / imputeFC$n_snps_in_cov
weight = imputeFC$pred_perf_r2
weights = weight * usedFrac
system.time( res <- permutationMultiLm( fc , net , weights , num=10 ) )

system.time( res <- permutationSingleLmMatrix( fc , net , weights , num=1000 ) )
sum(res$empiricalPval<=0.05,na.rm=T)







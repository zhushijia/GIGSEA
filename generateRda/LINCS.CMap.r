library("cmapR" , lib="/hpc/users/zhus02/fangg03a/sjzhu/Projects/DrugRepurposing/LINCS/cmapR")
gct = parse.gctx("/hpc/users/zhus02/fangg03a/sjzhu/Projects/DrugRepurposing/LINCS/GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx")
net = gct@mat

setwd('/hpc/users/zhus02/fangg03a/sjzhu/Projects/DrugRepurposing/LINCS/')
gene = read.table('GSE92742_Broad_LINCS_gene_info.txt',sep='\t',header=T,quote='"')
drug = read.table('GSE70138_Broad_LINCS_sig_info.txt',sep='\t',header=T,quote='"')
rangeI = match( gct@rid , as.character(gene$pr_gene_id)  )
rangeJ = match( gct@cid , as.character(drug$sig_id)  )
gene = gene[rangeI,]
drug = drug[rangeJ,]
all( gct@rid == as.character(gene$pr_gene_id) )
all( gct@cid == as.character(drug$sig_id) )
rownames(net) = as.character(gene$pr_gene_symbol)
colnames(net) = as.character(drug$pert_iname)
net = net[order(rownames(net)) , ]
group = colnames(net)

net2 = list()
for( i in 1:nrow(net) )
{
	cat(i,'\n')
	net2[[i]] = tapply(net[i,],group,mean)
}

net3 = lapply(net2,function(x) as.integer( x*1e3 ) )
net3 = do.call(rbind,net3)
colnames(net3) = names(tapply(net[1,],group,mean))
rownames(net3) = rownames(net)

annot = drug[,3:7]
LINCS.CMap.drug = list(net=net3,annot=annot)

save( LINCS.CMap.drug , file="LINCS.CMap.drug.rda" )





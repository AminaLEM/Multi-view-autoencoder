#load our resulted matrix of [patients X pathways] and clinical data
load(data_matrix)
load(clinical_data)
load(Mutation_data)

# Clustering 
library(factoextra)
library(NbClust)
fviz_nbclust(data_matrix,hcut, "silhouette")
res = NbClust(data_matrix, distance = "euclidean", min.nc=2, max.nc=8, 
              method = "kmeans", index = "all")

res$All.index
res$Best.nc
clusters = res$Best.partition
# map clinical_data from LIU and al. to awer samples
for (i in 1:length(samples))
{ if(as.character(samples[i]) %in% clinical_data[,1])
{
  clinical_dataIU[i,]=as.matrix(clinical_data[which(clinical_data[,1]== as.character(samples[i])),])
}}

#compute p values
#OS
survdiff(Surv(as.numeric(clinical_dataIU[,16]),as.numeric(clinical_dataIU[,15])) ~ group)
#DSS
survdiff(Surv(as.numeric(clinical_dataIU[,18]),as.numeric(clinical_dataIU[,17])) ~ group)
#DFI
survdiff(Surv(as.numeric(clinical_dataIU[,20]),as.numeric(clinical_dataIU[,19])) ~ group)

#chisq.test for the other categorical clinical features
for(i in 4:14)
{ print( clinical_dataIU[1,i])
  print(chisq.test(table(group,clinical_dataIU[,i]))$p.value)
}
#correction of P value of survival analysis (if correlation between age and clusters is proved)
# aply likelihood ratio test 
anova(lm(as.numeric(as.matrix(clinical_dataIU[,3])) ~ group))
null_model= coxph(Surv(as.numeric(clinical_dataIU[,16]),as.numeric(clinical_dataIU[,15])) ~ as.numeric(as.matrix(clinical_dataIU[,3])))
alternative_model= coxph(Surv(as.numeric(clinical_dataIU[,16]),as.numeric(clinical_dataIU[,15])) ~ as.numeric(as.matrix(clinical_dataIU[,3])) + group)
anova(null_model, alternative_model)

#Kaplan meier curves
fit_s=survfit(Surv(as.numeric(clinical_dataIU[,16]),as.numeric(clinical_dataIU[,15])) ~ clusters)
plot(fit_s,col=c(1:4),xlab = "time",ylab = "survival probability")
legend ("topright",c("subtype1","subtype2"), col=(1:4), lwd=0.5,cex=0.8)

#Mutation enrichment (fisher exact test)
len= dim(Mutation_data)[2] 

mutation<- vector(mode="character", length=len)
cluster1=  vector(mode="double", length=len)
cluster2 =  vector(mode="double", length=len)
test = data.frame(mutation,cluster1,cluster2)
test$mutation= ""

#create adjacency matrix
for(j in 1:len)
{
  test$mutation[j]= colnames(Mutation_data)[j]
  successsample = length(which(Mutation_data[which(clusters == 1),j]==1))
  successleftpart=  length(which(Mutation_data[,j]==1))-successsample
  failuresamples=length(which(clusters == 1)) - successsample
  failureleftpart = length(which(Mutation_data[,j]==0))-failuresamples
  test$cluster1[j] = fisher.test(matrix(c(successsample, successleftpart,failuresamples, failureleftpart), 2, 2), alternative='greater')$p.value; 
  
  
}
#BH correction
test$cluster1 = p.adjust(test$cluster1 , 'BH')


for(j in 1:len)
{
  successsample = length(which(Mutation_data[which(clusters == 2),j]==1))
  successleftpart=  length(which(Mutation_data[,j]==1))-successsample
  failuresamples=length(which(clusters == 2)) - successsample
  failureleftpart = length(which(Mutation_data[,j]==0))-failuresamples
  test$cluster2[j] = fisher.test(matrix(c(successsample, successleftpart,failuresamples, failureleftpart), 2, 2), alternative='greater')$p.value; 
  
  
}
test$cluster2 = p.adjust(test$cluster2 , 'BH')
enrichment = subset(test, mutation %in% test$mutation[which(test$cluster1 <0.05)] 
                    |        mutation %in% test$mutation[which(test$cluster2 <0.05)] )

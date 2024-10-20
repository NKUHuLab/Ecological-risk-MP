library(xlsx)
library(readxl)
library(hydroGOF)
library(randomForest)
library(iml)
library(ggplot2)
library(circlize)
library(RColorBrewer)
library(dplyr)
library(randomForestExplainer)
library(iml)
library(pdp)
library(tcltk)
library(patchwork)
library(raster)
library(ggbreak)
library(reshape2)
library(grid)
library(ggpointdensity)
library(ggsci)
library(openxlsx)
library(writexl)
library(caret)
library(igraph)
library(rfPermute)
library(plspm)
library(ggalt)
library(ggpointdensity)
library(viridis)
library(hexbin)
library(cowplot)
library(mgcv)

setwd('E:/Phd/microplastic Soil/microplastic Soil second')
land <- c('All', 'Agriculture', 'Natural ecosystems', 'Urban')
#################################model built#####################
Test_predict_wb <- createWorkbook()
Test_predict_10fold_wb <- createWorkbook()
Train_predict_10fold_wb <- createWorkbook()
for (i in 1:length(land)){
  addWorksheet(Test_predict_wb,sheetName = land[i])
}
for (i in 1:length(land)){
  addWorksheet(Test_predict_10fold_wb,sheetName = land[i])
}
for (i in 1:length(land)){
  addWorksheet(Train_predict_10fold_wb,sheetName = land[i])
}
for (i in 1:length(land)){
  data <- read.xlsx('Dataset/model data.xlsx', i)
  data_model <- data[,c(38:ncol(data),12)]
  colnames(data_model)[ncol(data_model)]<-'index'
  data_model$index<-as.numeric(data_model$index)/4824
  
  set.seed(1234);disorder <- sample(nrow(data_model),replace=F)
  fold_num <- floor(nrow(data_model)/10)
  
  n_test <- data.frame()
  n_train <- data.frame()
  predict <- data.frame()
  
  for (k in 1:10){
    o <- disorder[(fold_num*(k-1)+1):(fold_num*k)]
    rf.data <- data_model[-o,]
    rf <- randomForest(index~. , data = rf.data, 
                       ntree=1000 ,mtry=10,
                       proximity = F,
                       importance = F)
    p <- as.data.frame(predict(rf, data_model[o,]))
    a <- as.data.frame(data_model[o,ncol(data_model)])
    p_train <- as.data.frame(predict(rf, rf.data))
    a_train <- as.data.frame(rf.data[,ncol(data_model)])

    r2 <- cor(p,a)
    R2 <- R2(p,a)
    rm <- Metrics::rmse(p[,1],a[,1])
    predict <- rbind(predict,cbind(a,p))
    n_test <- rbind(n_test,cbind(r2, R2, rm))
    #train
    r2_train <- cor(p_train,a_train)
    R2_train <- R2(p_train,a_train)
    rm_train <- Metrics::rmse(p_train[,1],a_train[,1])
    n_train <- rbind(n_train,cbind(r2_train, R2_train, rm_train))

    print(k)
  }
  colnames(predict) <- c('Observation','Prediction')
  colnames(n_test) <- c('r2','R2','RMSE')
  colnames(n_train) <- c('r2','R2','RMSE')
  
  writeData(Test_predict_wb, sheet = i, predict)
  writeData(Test_predict_10fold_wb, sheet = i, n_test)
  writeData(Train_predict_10fold_wb, sheet = i, n_train)
  print(i)
  print(cor(predict))
}
saveWorkbook(Test_predict_wb, "analysis/model build/Test_predict.xlsx", overwrite = TRUE)
saveWorkbook(Test_predict_10fold_wb, "analysis/model build/Test_predict_10fold.xlsx", overwrite = TRUE)
saveWorkbook(Train_predict_10fold_wb, "analysis/model build/Train_predict_10fold.xlsx", overwrite = TRUE)
###########################importance#####################
###including MSE increase; node purity increase; SHAP
importance_wb <- createWorkbook()
for (i in 1:length(land)){
  addWorksheet(importance_wb,sheetName = land[i])
}
for (i in 1:length(land)){
  data <- read.xlsx('Dataset/model data.xlsx', i)
  data_model <- data[,c(38:ncol(data),12)]
  colnames(data_model)[ncol(data_model)]<-'index'
  data_model$index<-as.numeric(data_model$index)/4824#计算risk value
  rf <- randomForest(index~. , data =  data_model, 
                     ntree=1000 ,mtry=10,
                     proximity = F,
                     importance = T)
  imp <- as.data.frame(rf$importance)
  imp[,3] <- rownames(imp)
  colnames(imp) <- c("MSE","Node","variables")
  writeData(importance_wb, sheet = i, imp)
  print(i)
}
saveWorkbook(importance_wb, "analysis/model build/Importance_three indicator.xlsx", overwrite = TRUE)
###save rda file
rf.list<-list()
for (i in 1:length(land)){
  data <- read.xlsx('Dataset/model data.xlsx', i)
  data_model <- data[,c(38:ncol(data),12)]
  colnames(data_model)[ncol(data_model)]<-'index'
  data_model$index<-as.numeric(data_model$index)/4824
  rf.list[[land[i]]]<-local({
    randomForest(index~. , data =  data_model, 
                 ntree=1000 ,mtry=10,
                 proximity = F,
                 importance = T)})
}
save(rf.list,file='analysis/model build/Rda/rflist.rda')
#multi-importance rda file
load(file='analysis/model build/Rda/rflist.rda')
md<-list()
mi<-list()
for (i in 1:length(land)){
  data <- read.xlsx('Dataset/model data.xlsx', i)
  data_model <- data[,c(38:ncol(data),12)]
  colnames(data_model)[ncol(data_model)]<-'index'
  dataset <- data_model
  colindex <- c(colnames(dataset))
  colindexf <- factor(colindex[-length(colindex)],levels=colindex[-length(colindex)])
  
  rf<-rf.list[[i]]
  min_depth_frame<-min_depth_distribution(rf)
  md[[land[i]]]<-min_depth_frame
  im_frame<-measure_importance(rf)
  im_frame[4]<-im_frame[4]/max(im_frame[4])
  im_frame[5]<-im_frame[5]/max(im_frame[5])
  mi[[land[i]]]<-im_frame
  print(i)
}
save(md,mi,file='analysis/model build/Rda/multi-importance.rda')
###multi-way plot
load(file='analysis/model build/Rda/rflist.rda')
load(file='analysis/model build/Rda/multi-importance.rda')
mdplot<-list()
miplot<-list()
for (i in 1:length(land)){
  print(i)
  min_depth_frame<-md[[i]]
  mdplot[[land[i]]]<-local({
    min_depth_frame=min_depth_frame
    plot_min_depth_distribution(min_depth_frame,k=15)+
      theme(axis.text=element_text(colour='black',size=10),
            axis.title=element_text(colour='black',size=15))+
      theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
      theme(legend.key.size = unit(0.5,'line'),legend.title = element_text(size=rel(0.6)),
            legend.text = element_text(size=rel(0.5)))
  })
  ggsave(paste0('analysis/model build/Plot/Importance/mp_',land[i],'.pdf'),width=8,height=8)
  
  im_frame=mi[[i]]
  im_frame$p_value<-im_frame$p_value/5
  miplot[[land[i]]]<-local({
    im_frame=im_frame
    plot_multi_way_importance(im_frame, x_measure = "mse_increase",
                              y_measure = "node_purity_increase",
                              size_measure = "p_value", no_of_labels = 15)+
      theme(axis.text=element_text(colour='black',size=10),
            axis.title=element_text(colour='black',size=15))+
      theme(axis.line=element_line(color='black'),
            axis.ticks.length=unit(0.5,'line'))+
      labs(x="MSE increase", y="Node Purity Increase")+
      theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
      coord_fixed(ratio=1)+
      theme(legend.position=c(0.1,0.8))
  })
  ggsave(paste0('analysis/model build/Plot/Importance/im_',land[i],'.pdf'),width=5,height=5)
}
save(mdplot,miplot,file='analysis/model build/Rda/Importanceplot.rda')
#####Feature Interaction Calculate######
inter_list<-list()
for (i in 1:length(land)){
  print(i)
  im_frame<-mi[[i]]
  rf<-rf.list[[i]]
  vars <- important_variables(im_frame, k = 5, measures = c("mean_min_depth","no_of_trees"))
  interactions_frame <- min_depth_interactions(rf, vars)
  interactions_frame <- arrange(interactions_frame,-interactions_frame[,4])
  inter_list[[land[i]]]<-interactions_frame
}
###ps:depth or occurrence
save(inter_list,file='analysis/model build/Rda/inter.rda')
###plot
fiplot<-list()
for (i in 1:length(land)){
  interactions_frame<-inter_list[[i]]
  hlim<-ceiling(max(interactions_frame[1:25,3],interactions_frame[1:25,6]))
  fip<-plot_min_depth_interactions(interactions_frame,k=25)+
    scale_y_continuous(limits=c(0,hlim+1.5),expand=c(0,0))+
    scale_fill_gradient(low='#00d538',high='#ff5e24')+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
    theme(legend.position=c(0.3,0.8),legend.box="horizontal")
  fiplot[[land[i]]]<-fip
  ggsave(paste0('analysis/model build/Plot/inter/inter',land[i],'.pdf'),width=7,height=4)
}
save(fiplot,file='analysis/model build/Rda/inter_plot.rda')
###pdp analysis
load(file='analysis/model build/Rda/inter.rda')
pdplist<-list()
pdpplot<-list()
for (i in 1:length(land)){
  data <- read.xlsx('Dataset/model data.xlsx', i)
  data_model <- data[,c(38:ncol(data),12)]
  colnames(data_model)[ncol(data_model)]<-'index'
  data_model$index<-as.numeric(data_model$index)/4824#计算risk value
  
  rf <- randomForest(index~. , data =  data_model, 
                     ntree=1000 ,mtry=10,
                     proximity = F,
                     importance = T)
  inter_frame<-inter_list[[i]]
  j=1
  k=1
  subpdp<-list()
  subpdpplot<-list()
  while (j<=4){
    interpair<-inter_frame$interaction[k]
    v1<-strsplit(interpair,':')[[1]][1]
    v2<-strsplit(interpair,':')[[1]][2]
    k=k+1
    if (v1!=v2) {
      par<-pdp::partial(rf,pred.var = c(v1, v2), chull = TRUE, progress = "text")
      subpdp[[j]]<-par
      j<-j+1
    } else {j<-j}
  }
  print(i)
  pdplist[[land[i]]]<-subpdp
  pdpplot[[land[i]]]<-subpdpplot
}
save(pdplist,file='analysis/model build/Rda/pdplist.rda')
save(pdpplot,file='analysis/model build/Rda/pdpplot.rda')
for (i in 1:length(land)){
  subpdp<-pdplist[[i]]
  for (j in 1:4){
    par<-subpdp[[j]]
    subpdpplot[[j]]<-local({
      par=par
      ggplot(par, aes(x = par[[1L]], y = par[[2L]],
                      z = par[["yhat"]], fill = par[["yhat"]])) +
        geom_tile()+
        geom_contour(color = 'white')+
        viridis::scale_fill_viridis(name =land[i], option = 'D') +
        theme_bw()+
        xlab(colnames(par)[1])+
        ylab(colnames(par)[2])+
        theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
        theme(axis.line=element_line(color='black'),axis.text = element_text(size = rel(0.6)),
              axis.ticks.length=unit(0.3,'line'),axis.title = element_text(size = rel(0.6)))+
        theme(legend.key.size = unit(0.5,'line'),legend.title = element_text(size=rel(0.6)),
              legend.text = element_text(size=rel(0.5)),legend.position = c(0.9,0.8))
      ggsave(paste0('analysis/model build/Plot/dvpdp/',land[i],'-',
                    colnames(par)[1],'-',colnames(par)[2],
                    '.pdf'),width=5,height=5)
    })
  }
  pdpplot[[land[i]]]<-subpdpplot
}
###################################################################
#import_data_continuous <- to import a txt,csv or excel input file
#containing continuous dataset into BMDS analysis
#
#
#
###################################################################
import_data_continuous<-function(fileName,is_individual=F,formatType,has_header=T){
  
  fileFormat = which(formatType==c('table_txt','csv','xlsx'),arr.ind=FALSE)
  
  print(sprintf("the file format index is:%d",fileFormat))
  
  if(identical(fileFormat,integer(0))){
    stop("please specify one of the following file format types:
         'table_txt','csv','xlsx'")
  }else if (fileFormat==1){
    mydata=read.table(fileName,header=has_header)
    
  }else if(fileFormat==2){
    mydata=read.csv(fileName,header=has_header)
    
  }else if(fileFormat==3){
    mydata=read_excel(fileName,sheet=1,col_names=has_header)
    mydata<-data.table(mydata)
  }else 
    
    return(mydata)
  
}

###################################################################
#print_ma_continuous <- to print tables or write them into csv files 
#containing continuous MA BMDS analysis results
#
#
#
###################################################################
print_ma_continuous<-function(obj,model_name,write_file=FALSE){
  UseMethod("print_ma_continuous")
}

print_ma_continuous.summaryTable<-function(obj,write_file=FALSE){
  N_rows<-length(names(obj$estimates))
  myTable<-data.frame(matrix(NA,nrow=N_rows,ncol=7))
  myTable[,1]<-names(obj$estimates)
  myTable[,2]<-obj$dose_scale*obj$ma$modelBMDS
  myTable[,3]<-obj$dose_scale*obj$ma$modelBMDLS
  myTable[,4]<-obj$dose_scale*obj$ma$modelBMDUS
  myTable[,5]<-obj$ma$modelprobs
  colnames(myTable)<-c("Model","BMD","BMDL","BMDU","Log Posterior Probability","scaled residual at BMD","scaled residual at control")
  print(myTable)
  if(write_file){
    write.csv(myTable,sprintf("./summary_table.csv"),row.names = FALSE,quote=FALSE)
  }
}

print_ma_continuous.modelResults<-function(obj,model_name,write_file=FALSE){
  myname=which(model_name==(names(obj$estimates)))
  ##print a small table for benchmark dose for each individual model
  firstTable<-data.frame(matrix(NA,nrow=4,ncol=2))
  firstTable[,1]<-c("BMD","BMDL","BMDU","LPP")
  firstTable[,2]<-c(obj$dose_scale*obj$ma$modelBMDS[myname],obj$dose_scale*obj$ma$modelBMDLS[myname],
                    obj$dose_scale*obj$ma$modelBMDUS[myname],obj$ma$modelprobs[myname])
  
  colnames(firstTable)<-c("","Benchmark Dose")
  print(firstTable)
  ##print a small table for model parameters for each individual model
  N_rows=nrow(data.frame(obj$estimates[myname]))
  secTable<-data.frame(matrix(NA,nrow=N_rows,ncol=2))
  row_names<-character(N_rows)
  for(i in 1:N_rows){
    part1<-"Parameter"
    part2<-i
    row_names[i]<-paste(part1,part2)
  }
  
  secTable[,1]<-row_names
  secTable[,2]<-data.frame(obj$estimates[myname])
  colnames(secTable)<-c("Variable","Model Parameters estimate")
  print(secTable)
  
  if(write_file){
    sink(file = sprintf("./BMDS_MA_table_for_%s.csv",model_name))
    write.csv(firstTable,row.names = FALSE,quote=FALSE)
    cat("\n")
    write.csv(secTable,row.names = FALSE,quote=FALSE)
    sink()
  }
}

###################################################################
#plot_ma_continuous <- to plot dose-response or CDF figures  
#containing continuous MA BMDS analysis results
#
#
#
###################################################################
plot_ma_continuous<-function(obj,DATA=mydata,save_figure=FALSE){
  UseMethod("plot_ma_continuous")
}

plot_ma_continuous.doseResponse<-function(obj,DATA=mydata,model_name="all",title="Dose response for all models",save_figure=FALSE){
  myname=which(model_name==c("all",names(obj$estimates))) 
  max_x=1.2*obj$dose_scale
  newdata<-seq(0,max_x,by=(max_x-0)/300)
  colnames(mydata)<-c("dose","response")
  max_y=1.2*max(mydata$response)
  min_y=0.7*min(mydata$response)
  modelNames=names(obj$estimates)
  modelList=list()
  
  for(i in 1:length(obj$estimates)){
    get_pars<-data.frame(obj$estimates[i])
    if (identical(names(obj$estimates)[i],"Hill Normal-NCV")||identical(names(obj$estimates)[i],"Hill Normal-CV")||
        identical(names(obj$estimates)[i],"Hill Log-Normal")){
      model_function = function(x){get_pars[1,]+get_pars[2,]*x^get_pars[4,]/(get_pars[3,]^get_pars[4,]+x^get_pars[4,])} 
    }else if(identical(names(obj$estimates)[i],"Power Normal-NCV")||identical(names(obj$estimates)[i],"Power Normal-CV")||
             identical(names(obj$estimates)[i],"Power Log-Normal")){
      model_function = function(x){get_pars[1,]+get_pars[2,]*x^get_pars[3,]} 
    }else if(identical(names(obj$estimates)[i],"Exp3  Normal-NCV")||identical(names(obj$estimates)[i],"Exp3  Normal-CV")||
             identical(names(obj$estimates)[i],"Exp3  Log-Normal")){
      model_function = function(x){get_pars[1,]*exp((get_pars[2,]*x)^get_pars[3,])} 
    }else if(identical(names(obj$estimates)[i],"Exp5 Normal-NCV")||identical(names(obj$estimates)[i],"Exp5 Normal-CV")||
             identical(names(obj$estimates)[i],"Exp5 Log-Normal")){
      model_function = function(x){get_pars[1,]*(exp(get_pars[3,])-(exp(get_pars[3,])-1)*exp(-(get_pars[2,]*x)^get_pars[4,]))} 
    }else if(identical(names(obj$estimates)[i],"Poly2 Normal-NCV")||identical(names(obj$estimates)[i],"Poly2 Normal-CV")||
             identical(names(obj$estimates)[i],"Poly2 Log-Normal")){
      model_function = function(x){get_pars[1,]+get_pars[2,]*x+get_pars[3,]*x^2} 
    }else{
      print("there is a model not recognizable!")
      print(names(obj$estimates)[i])
      
    }  
    
    y_to_add<-model_function(newdata)
    modelList[[i]]<-y_to_add
  }
  
  all_models=do.call(cbind,modelList)
  colnames(all_models)<-names(obj$estimates)
  ## add the MA estimate curve
  cal_y_ma<-t(apply(all_models,1,"*",obj$ma$modelprobs))
  y_ma<-data.frame(newdata,rowSums(cal_y_ma))
  df<-data.frame(newdata,all_models)
  colnames(y_ma)<-c("new_dose","model_average")
  #print(y_ma)
  df.melted<-melt(df,id="newdata")
  #print(df.melted)
  myPlot<-ggplot(data=df.melted,aes(x=newdata,y=value,color=variable))+geom_line(linetype="dashed")
  myPlot<-myPlot+geom_line(data=y_ma,mapping=aes(x=new_dose,y=model_average),colour="black")
  myPlot<-myPlot+geom_point(data=mydata,mapping=aes(x=dose,y=response),colour="black",shape=16)+
    xlab("Dose")+ylab("Response")+xlim(0,max_x)+ylim(min_y,max_y)+theme_bw()
  myPlot<-myPlot+labs(title=title)+theme(plot.title=element_text(hjust=0.5,size=14),legend.title=element_blank())
  
  print(myPlot)
  if(save_figure==TRUE){
    save_plot("./dose_response_plot.png",myPlot)
  }
}

plot_ma_continuous.bmdCDF<-function(obj,title="MA BMD CDF",save_figure=FALSE,alpha=0.05){
  df<-data.frame(obj$ma$MA_BMD_CDF)
  df <- df[-1,]
  df <- df[df[,2] < 0.999,,drop=F]
  df[,1]<-obj$dose_scale*df[,1]
  colnames(df)<-c("BMD","Cumulative_density","no_use")
  myPlot<-ggplot2::ggplot(data=df,aes(x=BMD,y=Cumulative_density))+geom_line(size=1.3)+
                  xlab("MA BMD")+ylab("Cumulative density")+theme_bw()
  myPlot<-myPlot+labs(title=title)+theme(plot.title=element_text(hjust=0.5,size=14))+
              geom_hline(yintercept=alpha,linetype=3,size=1.1) + geom_hline(yintercept=1-alpha,linetype=3,size=1.1)
  print(myPlot)
  
  if(save_figure==TRUE){
    save_plot("./MA_BMD_CDF_plot.png",myPlot)
  }
}


get_ma_bmd<-function(obj,alpha=0.05){
  if (alpha > 1 || alpha < 0)
    stop("alpha parameter must be in the [0,1] range")
  
  BMD <- obj$ma$MA_BMD_CDF
  BMD[,1] = BMD[,1]*obj$dose_scale; 
  bmd_ma = splinefun(BMD[,2],BMD[,1],method="monoH.FC")
  return(bmd_ma(alpha))
}





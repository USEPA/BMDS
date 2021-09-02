#Copyright 2021  NIEHS <matt.wheeler@nih.gov>
#   
#
#Permission is hereby granted, free of charge, to any person obtaining a copy of this software 
#and associated documentation files (the "Software"), to deal in the Software without restriction, 
#including without limitation the rights to use, copy, modify, merge, publish, distribute, 
#sublicense, and/or sell copies of the Software, and to permit persons to whom the Software 
#is furnished to do so, subject to the following conditions:
#
#The above copyright notice and this permission notice shall be included in all copies 
#or substantial portions of the Software.

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
#INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A 
#PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT 
#HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF 
#CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
#OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

## ------------------------------------------------------------------------	
## DUNNETT'S TEST: http://www.stat.wmich.edu/wang/664/egs/Rmice.html
## ------------------------------------------------------------------------

dunnett_ntp <- function(formula, data,dose_name = "dose"){
  
  data[,c(dose_name)] = as.numeric(data[,c(dose_name)])
  temp_str = strsplit(as.character(formula)[3], " ")[[1]]
  temp_str = temp_str[temp_str != "+"]
  data = data[order(data[,c(dose_name)]),]
  
  for (ii in 1:length(temp_str)){
    data = data[order(data[,temp_str[ii]]),]
  }
  
  if (!(dose_name %in% colnames(data))){
    stop(sprintf("Dose name %s does not appear in the data frame.",dose_name))
  }
  
  jonck_data =  jonckeere_ntp( formula,dose_name = dose_name,
                               data = data,pair="Williams")
  
  
	## loop through all groups flagged as DUNNETT in jonck
	dunnett <- subset(jonck_data, mult_comp_test=='DUNNETT')
	dunnett_results <- NULL

	
	temp_colnames <- unlist(c(unlist(colnames(dunnett)),dose_name,as.character(unlist(formula[[2]]))))
	temp <-   colnames(data) %in% temp_colnames
	temp_d <- as.data.frame(data[,temp==TRUE])
	
	
	if(nrow(dunnett) > 0){
    		for(d in 1:nrow(dunnett)){
    		  
    		  temp_names <- rep(NA,ncol(dunnett))
    		  for (j in 1:length(temp_names)){
    		    temp_names[j] <- unlist(dunnett[d,j])
    		  }
    		  ##KRS - changed "phase" to "phasetype"
    		  temp_dd <- temp_d
    		  temp_dd[is.na(temp_dd)] = ''
    		  #subset the data based upon the columns that exist in the formula
    		  for (j in 1:length(dunnett)){
    		    myt <- which(colnames(temp_d) == colnames(dunnett)[j])
    		    
    		    if (length(myt)>0){
    		      temp_dd <- temp_dd[as.character(temp_dd[,myt])==as.character(temp_names[j]),]
    		    }
    		    
    		  }
    		  
    			datatemp <- temp_dd
    			temp_idx = which(colnames(datatemp) == dose_name )
    			datatemp <- datatemp[order(datatemp[,temp_idx]),]
    
    			## duplicate dose classes from old SAS conversion
    			datatemp$dose3 <- as.factor(datatemp[,temp_idx])
    	
    			if(length(unique(datatemp[,temp_idx])) > 1){
    			  temp_idx2 <- which(!(colnames(dunnett) %in% c("tau","pvalue","mult_comp_test")))
    				trts  <- unique(datatemp[,temp_idx])[-1]		## get unique treatment doses, this will be used in inner loop below
    				formula_text  = sprintf( "%s ~ dose3" ,as.character(formula[[2]]))
    				data.aov  <- aov(as.formula(formula_text), data=datatemp)
    				data.dunn <- glht(data.aov, linfct=mcp(dose3="Dunnett"))
    
    				## loop through doses within sex/endpoint combo
    				for(j in 1:length(summary(data.dunn)$test$pvalues)){
    				  line <- c(temp_names[temp_idx2],trts[j],summary(data.dunn)$test$tstat[j], summary(data.dunn)$test$pvalues[j])
    				  dunnett_results          <- rbind(dunnett_results, line)
    				  
    				}
    				colnames(dunnett_results) <- c(colnames(dunnett)[temp_idx2],dose_name,"tstat","pvalue")
    			}
    	 }
    			
    	  if(!is.null(dunn_results)){
    			rownames(dunnett_results) <- NULL
    			dunnett_results <- as.data.frame(dunnett_results)		
    			
    			dunnett_results$pvalue <- as.numeric(as.character(dunnett_results$pvalue))
    			dunnett_results$mult_comp_test <- 'DUNNETT'
    
    			## determine how many asterisks each row deserves
    			dunnett_results$mult_comp_signif <- ifelse(dunnett_results$pvalue <= .01, 2,
    								   			                  ifelse(dunnett_results$pvalue <= .05, 1, 0))
    
    
    			temp_idx3 = which(colnames(dunnett_results)== dose_name)
    			dunnett_results[,temp_idx3] =  as.numeric(dunnett_results[,temp_idx3])
    	 } 
  } 
  
	return(dunnett_results)
}
	
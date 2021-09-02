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

## -----------------------------------------------------------
## JONCKHEERE'S TEST 
## ----------------Changelog----------------------------------
## Released: 
## [1.1.0] - 08/07/2019. Jira ticket:CEBSPROV-5301
## Changed
## Fixes for Kendall test. Set Exact = FALSE
## Changed to make it based upon a general formula specified in 
## formula. To do this, we assume that data is a data frame. 
## As a default, "dose_name", is set to the column header "dose"
## -----------------------------------------------------------
jonckeere_ntp<- function(formula, data, dose_name="dose", pair = 'Williams' )
{
  if (!is.data.frame(data)){
    stop("The data do not appear to be in the data frame format.")
  }
  if ( sum (colnames(data) %in% dose_name) == 0 ){
    stop(sprintf("The specified dose column '%s' does not exist in the data frame.",dose_name))
  }
  if (!(pair %in% c('Williams','Shirley'))){
    stop("Variable 'pair' must be either 'Williams' or 'Shirley'")
  }
  
  
  jonck <- NULL
  ##KRS - added "numeric_value" on the left hand side
  jonck_list <- as.data.frame(summaryBy(formula , data=data, FUN=length))
  jonck_list[is.na(jonck_list)] = ''
  
  ## may create WARNINGS when ties are present
  analysis_var <- as.character(unlist(formula[[2]]))
  temp_colnames <- unlist(c(unlist(colnames(jonck_list)),dose_name,as.character(unlist(formula[[2]]))))
  temp <-   colnames(data) %in% temp_colnames
  temp_d <- as.data.frame(data[,temp==TRUE])
  
  for(i in 1:nrow(jonck_list))
  {
    temp_names <- rep(NA,ncol(jonck_list))
    for (j in 1:length(temp_names)){
      
      temp_names[j] <- unlist(jonck_list[i,j])
    }
    
    ##KRS - changed "phase" to "phasetype"
    temp_dd <- temp_d
    temp_dd[is.na(temp_dd)] = ''
    
    for (j in 1:length(jonck_list)){
      myt <- which(colnames(temp_d) == colnames(jonck_list)[j])
      
      if (length(myt)>0){
        temp_dd <- temp_dd[as.character(temp_dd[,myt])==as.character(temp_names[j]),]
      }
      
    }
    tempdata = temp_dd
    tidx = which(colnames(tempdata) == dose_name)
    #print(tempdata)
    ## make sure there is more than ONE record for gender i
    if(nrow(tempdata) > 1)
    {
      ## make sure variance is not zero (creates error) and control is present
      if((length(unique(tempdata$numeric_value)) > 1) & (0 %in% as.numeric(tempdata[,tidx])))
      {
        
        stat <- cor.test(as.numeric(tempdata[,tidx]), tempdata$numeric_value, method="kendall", exact=FALSE)
        jline <- c(temp_names, stat$estimate, stat$p.value)
        jonck <- rbind(jonck, jline)
        
      }
    }
  }
  
  ## check for existence
  if(!is.null(jonck))
  {
    rownames(jonck) <- NULL
    jonck <- as.data.frame(jonck)
    names(jonck) <- c(colnames(jonck_list), 'tau', 'pvalue')
    jonck$tau   <- as.numeric(as.character(jonck$tau))
    jonck$pvalue <- as.numeric(as.character(jonck$pvalue))
    
    ## need to remove nulls 
    jonck[is.na(jonck)] = ''
    
    temp_col = sprintf("%s.length",as.character(formula[[2]]))
    temp_idx = which(temp_col == colnames(jonck))
    jonck = jonck[,-temp_idx]
    
    if(pair=='Williams') { jonck$mult_comp_test <- ifelse(jonck$pvalue <= .01, 'WILLIAMS', 'DUNNETT')	}
    if(pair=='Shirley')  { jonck$mult_comp_test  <- ifelse(jonck$pvalue <= .01, 'SHIRLEY', 'DUNN')	}
  
  }
  return(jonck)
}

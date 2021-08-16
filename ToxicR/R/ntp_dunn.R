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

dunn_ntp <- function(formula,dose_name = "dose", data)
{
  
  #first do Jonckheere's Test and subset all the 'DUNN' tests
  jonck_data =  jonckeere_ntp( formula,dose_name = dose_name,
                            data = data,pair="Shirley")
	dunn <- subset(jonck_data, mult_comp_test=='DUNN')
	dunn_results <- NULL
	

  temp_colnames <- unlist(c(unlist(colnames(dunn)),dose_name,as.character(unlist(formula[[2]]))))
  temp <-   colnames(data) %in% temp_colnames
  temp_d <- as.data.frame(data[,temp==TRUE])
  
  ## loop through all groups flagged as DUNN in jonck
	if(nrow(dunn) > 0){
	  
		for(d in 1:nrow(dunn)){
		  
		  temp_names <- rep(NA,ncol(dunn))
		  for (j in 1:length(temp_names)){
		    
		    temp_names[j] <- unlist(dunn[d,j])
		  }
		  
		  ##KRS - changed "phase" to "phasetype"
		  temp_dd <- temp_d
		  temp_dd[is.na(temp_dd)] = ''
		  #subset the data based upon the columns that exist in the formula
		  for (j in 1:length(dunn)){
		    myt <- which(colnames(temp_d) == colnames(dunn)[j])
		    
		    if (length(myt)>0){
		      temp_dd <- temp_dd[as.character(temp_dd[,myt])==as.character(temp_names[j]),]
		    }
		    
		  }
		  ex = temp_dd
		  dose_idx = which(colnames(ex) == dose_name)
      nval_idx = which(colnames(ex) == as.character(unlist(formula[[2]])))
		  ex[,dose_idx] = as.numeric(ex[,dose_idx])
      
			## make sure there are multiple doses to compare, and control is present
			if((length(unique(ex[,dose_idx])) > 1) & (0 %in% unique(ex[,dose_idx])))
				{

				## rank among each sex/test_name/phase/etc.
				ties <- as.data.frame(table(table(ex[,nval_idx])))

				if(nrow(ties) > 0)
					{
					names(ties) <- c('numTies','count')
					ties$numTies <- as.integer(as.character(ties$numTies))

					## remove singletons
					ties <- subset(ties, numTies != 1)

					## tie adjustment correction
					correction <- 0		
					if(nrow(ties) > 0)
						{
						for(i in 1:nrow(ties))
							{
							correction <- correction + (ties$count[i] * (ties$numTies[i]^3 - ties$numTies[i]))
							}
						}		
				}
				## rankings ... lowest value gets rank of 1
				ranks <- as.data.frame(table(ex[,nval_idx]))
				names(ranks) <- c('value', 'count')
				ranks$value  <- as.numeric(as.character(ranks$value))
				weights <- NULL
				for(i in 1:nrow(ranks))
					{
					if(ranks$count[i]==1)
						{
						weights <- c(weights, sum(ranks[1:i,'count']))	## if a singleton, just add up previous number of values to get rank
						} else
						{
						start <- sum(ranks[1:i-1,'count'])		## get starting point (previous rank)
						ranknum <- 0
						for(j in 1:ranks$count[i])				
							{
							ranknum <- ranknum + start + j		
							}

						ranknum <- ranknum / ranks$count[i]
						weights <- c(weights, ranknum)
						}
					}

				ranks2 <- cbind(ranks, weights)
				names(ranks2)[3] <- 'rank'

				## populate ex with ranks
				## need to recast as numeric, may crash otherwise (weird!)
				ex[,nval_idx] <- as.character(ex[,nval_idx])
				ex[,nval_idx] <- as.numeric(ex[,nval_idx])

				ex2 <- merge(ex, ranks2, by.x=as.character(unlist(formula[[2]])), by.y='value', all.x=TRUE)
				dose_idx = which(names(ex2) == dose_name)
				ex2 <- ex2[order(ex2[,dose_idx]),]
				formula_text  = sprintf( "rank ~ %s" ,dose_name)
				rankMeans <- summaryBy(as.formula(formula_text) , data=ex2, FUN=mean)	

				## enumerate doseCount
				count <- 0
				for(i in 1:nrow(rankMeans))
					{
					rankMeans$doseCount[i] <- count
					count <- count + 1
					}


				## find variance
				V  <- (nrow(ex) * (nrow(ex) + 1))/12

				## get crit values ... warnings are for CONTROL group, which is fine
				prob05 <- qnorm(1 - (.05 / (2 * max(rankMeans$doseCount) )))	
				prob01 <- qnorm(1 - (.01 / (2 * max(rankMeans$doseCount) )))	

				## loop through doses, determine significance
				rankMeans$DUNSIGN <- 0
				rankMeans$num <- 0
			
				for(i in 2:nrow(rankMeans))		## skip CONTROL (first line)
				{

					## populate DUNSIGN (flag for mean of treatment group being greater than control group)
					if(rankMeans$rank.mean[i] >= rankMeans$rank.mean[1])
						{
						rankMeans$DUNSIGN[i] <- 0
						} else
						{
						rankMeans$DUNSIGN[i] <- -1
						}

					## find significance level
				  rm_idx = which(dose_name == names(rankMeans) )
				  ex_idx = which(dose_name == names(ex2))
					rankdiff <- abs(rankMeans$rank.mean[i] - rankMeans$rank.mean[1])
					num   <- nrow(ex2[ex2[,ex_idx]  == rankMeans[i,rm_idx],])
					comp2 <- V * ( (1/num) + 1/nrow(ex2[ex2[,ex_idx] == rankMeans[1,rm_idx],])) 
					comp2 <- (comp2 * (1 - correction / (nrow(ex2)^3 - nrow(ex2))))^.5
					
					sig05 <- rankdiff - (prob05 * comp2)
					
				  ptest = min(1,(1 - pnorm(rankdiff/comp2))*(2 * max(rankMeans$doseCount)))
				#  message(sprintf("%f %f %f %f %f",rankdiff,comp2,prob05,sig05,ptest))
					sig01 <- rankdiff - (prob01 * comp2)
					signific <- 0
					if(sig05 > 0)
						{
						signific <- 1
						}

					if(sig01 > 0)
						{
						signific <- 2
						}

					rankMeans$num[i] <- num
					rankMeans$pvalue[i]      <- ptest
					}

				results <- rankMeans
				
				results_len = ncol(results)
			  temp_names = colnames(dunn)[unlist(colnames(dunn)) %in% unlist(colnames(data))]
				
			  for (ii in 1:length(temp_names)){
			    idx = which(temp_names[ii] == colnames(ex2))
			    results[,ii+results_len] = ex2[1,idx]
			    
			  }
			  A = results[,results_len + 1:length(temp_names)]
			  colnames(A) <- temp_names
			  results = cbind(results[,1:results_len], A)
			  
				dunn_results <- rbind(dunn_results, results)
				}
			}
		}

	## check for existence	
	if(!is.null(dunn_results))
		{
		## coerce into same format as SHIRLEY results
		dunn_results$TEST <- 'DUNN'
	
		idx = which(colnames(dunn_results) == dose_name)
		## cut out controls
		dunn_results <- dunn_results[as.numeric(dunn_results[,idx]) != 0,]
	}
  
  names_to_drop <- colnames(dunn_results)
  temp <- sprintf("%sCount",dose_name)
  dose_idx    = which(colnames(dunn_results) == dose_name)
  p_value_idx = which(colnames(dunn_results) == "pvalue")
  test_idx    = which(colnames(dunn_results) == "TEST")
  remain_idx  = which(!(1:ncol(dunn_results) %in% c(dose_idx,p_value_idx,test_idx)))
  dunn_results = dunn_results[,c(dose_idx,remain_idx,test_idx,p_value_idx)]
  
  return(dunn_results[,-which(names_to_drop %in% c("rank.mean",temp,"DUNSIGN","num"))])

}


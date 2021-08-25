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

.compute_crit_williams <- function(william_test_data,dose_name,formulaV){
  
  
  t_idx = which(colnames(william_test_data) == dose_name)
  william_test_data[,t_idx] = as.numeric(william_test_data[,t_idx])
  william_test_data$crit05 = NA
  william_test_data$crit01 = NA
  
  for(k in 1:nrow(william_test_data)){
    ## CONTROL GROUP
    if(william_test_data[k,t_idx]==0)		## should this be datatemp or datatemp?
    { 
      william_test_data$crit05 <- FALSE
      william_test_data$crit01 <- FALSE
     
    } else if(william_test_data[k,t_idx] != 0 & (william_test_data$dof[k] %in% will005$dof)){
        col1 <- paste('w1crit', k, sep='')
        col5 <- paste('w5crit', k, sep='')
        adj1 <- paste('w1adj', k, sep='')
        adj5 <- paste('w5adj', k, sep='')
        
        w1crit <- subset(will005, dof==william_test_data$dof[k])[,c(col1)]
        w1adj  <- subset(will005, dof==william_test_data$dof[k])[,c(adj1)]
        w5crit <- subset(will025, dof==william_test_data$dof[k])[,c(col5)]
        w5adj  <- subset(will025, dof==william_test_data$dof[k])[,c(adj5)]
        
        dofactor <- ((william_test_data$dof[k] - lowdof) / (highdof - lowdof))	
        temp_name <- sprintf("%s.length",formulaV)
        t_idx2   <- which(colnames(william_test_data) == temp_name)
        
        con_num <- william_test_data[william_test_data[,t_idx] == 0,][,t_idx2]
        trt_num <- william_test_data[k,t_idx2]
        
        
        william_test_data$crit01[k] <- w1crit - (.1 * w1adj * (1 - (trt_num / con_num)))
        william_test_data$crit05[k] <- w5crit - (.1 * w5adj * (1 - (trt_num / con_num)))
      } else if(william_test_data[k,t_idx] != 0 & !(william_test_data$dof[k] %in% will005$dof))  
        {
          col1 <- paste('w1crit', k, sep='')
          col5 <- paste('w5crit', k, sep='')
          adj1 <- paste('w1adj', k, sep='')
          adj5 <- paste('w5adj', k, sep='')
          
          ## get lower bound from table
          lowdof <- max(will005$dof[william_test_data$dof[k] > will005$dof])
          
          low.w1crit <- subset(will005, dof==lowdof)[,c(col1)]
          low.w1adj  <- subset(will005, dof==lowdof)[,c(adj1)]
          low.w5crit <- subset(will025, dof==lowdof)[,c(col5)]
          low.w5adj  <- subset(will025, dof==lowdof)[,c(adj5)]
  
          ## get upper bound from table
          highdof <- min(will005$dof[william_test_data$dof[k] < will005$dof])
          
          high.w1crit <- subset(will005, dof==highdof)[,c(col1)]
          high.w1adj  <- subset(will005, dof==highdof)[,c(adj1)]
          high.w5crit <- subset(will025, dof==highdof)[,c(col5)]
          high.w5adj  <- subset(will025, dof==highdof)[,c(adj5)]
          
          dofactor <- ((william_test_data$dof[k] - lowdof) / (highdof - lowdof))	
          temp_name <- sprintf("%s.length",formulaV)
          t_idx2   <- which(colnames(william_test_data) == temp_name)
          
          con_num <- william_test_data[william_test_data[,t_idx] == 0,][,t_idx2]
          trt_num <- william_test_data[k,t_idx2]
          
          william_test_data$crit01[k] <- (low.w1crit - (dofactor * (low.w1crit - high.w1crit))) - (.01 * low.w1adj * (1 - (trt_num / con_num)))
          william_test_data$crit05[k] <- (low.w5crit - (dofactor * (low.w5crit - high.w5crit))) - (.01 * low.w5adj * (1 - (trt_num / con_num)))
        }
  }
  return(william_test_data)
}

## ----------------------
## 	WILLIAM'S TEST
## ----------------------
williams_ntp <- function(formula, data,dose_name = "dose"){
  
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
                               data = data, pair = "Williams")
  
	## loop through all groups flagged as WILLIAM in jonck
	william       <- subset(jonck_data, mult_comp_test=='WILLIAMS')
	will_results  <- NULL
	will_results2 <- NULL
	
	temp_colnames <- unlist(c(unlist(colnames(william)),dose_name,as.character(unlist(formula[[2]]))))
	temp <-   colnames(data) %in% temp_colnames
	temp_d <- as.data.frame(data[,temp==TRUE])
	
	if(nrow(william) > 0){
	  
		for(w in 1:nrow(william)){
		  
		  temp_names <- rep(NA,ncol(william))
		  for (j in 1:length(temp_names)){
		    temp_names[j] <- unlist(william[w,j])
		  }
		  ##KRS - changed "phase" to "phasetype"
		  temp_dd <- temp_d
		  temp_dd[is.na(temp_dd)] = ''
		  #subset the data based upon the columns that exist in the formula
		  for (j in 1:length(william)){
		    myt <- which(colnames(temp_d) == colnames(william)[j])
		    
		    if (length(myt)>0){
		      temp_dd <- temp_dd[as.character(temp_dd[,myt])==as.character(temp_names[j]),]
		    }
		    
		  }
		  
		  datatemp <- temp_dd
		  temp_idx = which(colnames(datatemp) == dose_name )
		  datatemp <- datatemp[order(datatemp[,temp_idx]),]
			## get william-ized dose means
			## direction of smoothing dependent on Jonckheere output
      formula_temp = sprintf("%s ~ %s + %s",as.character(formula[[2]]),as.character(dose_name),as.character(formula)[3])
			wmeans  <- summaryBy(as.formula(formula_temp), data=datatemp, FUN=c(mean,length))
			temp_idx2 <- which(colnames(wmeans) == dose_name)
			wmeans  <- wmeans[order(wmeans[,temp_idx2]),]
			smeans  <- wmeans$numeric_value.mean
			slength <- wmeans$numeric_value.length

			smeans  <- smeans[-1] 	## remove control group
			slength <- slength[-1] 	## remove control group

			trend <- ifelse(william$tau[w] < 0, william$pvalue[w] * -1, william$pvalue[w])
			## set comparison direction based on JONCK trend result
			direction <- ifelse(trend < 0, 'decreasing', 'increasing')
			
			## need more than 1 trt grp for smoothing
			if(length(smeans) > 1)
			{
				
				if(direction == 'decreasing')
					{
					for(i in 1:(length(smeans)-1) )
						{
						if(smeans[i] < smeans[i+1])		## smoothing required
							{

							if(i==1)		## FIRST ITEM
								{
								tempSmean   <- (smeans[i]*(slength[i]/(slength[i]+slength[i+1]))) + (smeans[i+1])*(slength[i+1]/(slength[i]+slength[i+1]))
								smeans[i]   <- tempSmean
								smeans[i+1] <- tempSmean
								} else

							if(i > 1)		## NOT FIRST ITEM
								{
								tempSmean   <- (smeans[i]*(slength[i]/(slength[i]+slength[i+1]))) + (smeans[i+1])*(slength[i+1]/(slength[i]+slength[i+1]))

								s_count <- 0			## initialize counter for previous consecutive smoothed means
								for(j in ((i-1):1) )	## loop back to find number of elements in smoothing
									{
									if(tempSmean > smeans[j])
										{ 
										s_count <- s_count + 1	## count once for each consecutive previous "smooth"					
										if(s_count > 0)
											{
											tempSmean <- 0
											for(k in (i-s_count:i+1))	## recalculate tempSmean
												{
												tempSmean <- tempSmean + smeans[k]*(slength[k]/sum(slength[i-s_count:i+1]))
												}
											}	
										} else
										{ break
										}
									}

								denom <- 0   					## this is denominator for weighting by sample size
								for(k in (i+1):(i - s_count))	## get count from previous smoothing
									{
									denom <- denom + slength[k]
									}


								tempSmean <- 0
								for(k in (i+1):(i - s_count))	## go back through and get smoothed mean(s)
									{
									tempSmean <- tempSmean + smeans[k]*slength[k]/denom
									}


								smeans[i]   <- tempSmean
								smeans[i+1] <- tempSmean
								if(s_count > 0)
									{
									for(k in 1:s_count)	
										{
										smeans[i-k] <- tempSmean
										}
									}			
								}
							}
						}
					} else 		## DIRECTION BREAK
					{
					for(i in 1:(length(smeans)-1) )
						{
						if(smeans[i] > smeans[i+1])		## smoothing required
							{

							if(i==1)		## FIRST ITEM
								{
								tempSmean   <- (smeans[i]*(slength[i]/(slength[i]+slength[i+1]))) + (smeans[i+1])*(slength[i+1]/(slength[i]+slength[i+1]))
								smeans[i]   <- tempSmean
								smeans[i+1] <- tempSmean
								} else

							if(i > 1)		## NOT FIRST ITEM
								{
								tempSmean   <- (smeans[i]*(slength[i]/(slength[i]+slength[i+1]))) + (smeans[i+1])*(slength[i+1]/(slength[i]+slength[i+1]))

								s_count <- 0			## initialize counter for previous consecutive smoothed means
								for(j in ((i-1):1) )	## loop back to find number of elements in smoothing
									{
									if(tempSmean < smeans[j])
										{ 
										s_count <- s_count + 1	## count once for each consecutive previous "smooth"					
										if(s_count > 0)
											{
											tempSmean <- 0
											for(k in (i-s_count:i+1))	## recalculate tempSmean
												{
												tempSmean <- tempSmean + smeans[k]*(slength[k]/sum(slength[i-s_count:i+1]))
												}
											}	
										} else
										{ break
										}
									}

								denom <- 0   					## this is denominator for weighting by sample size
								for(k in (i+1):(i - s_count))	## get count from previous smoothing
									{
									denom <- denom + slength[k]
									}


								tempSmean <- 0
								for(k in (i+1):(i - s_count))	## go back through and get smoothed mean(s)
									{
									tempSmean <- tempSmean + smeans[k]*slength[k]/denom
									}


								smeans[i]   <- tempSmean
								smeans[i+1] <- tempSmean
								if(s_count > 0)
									{
									for(k in 1:s_count)	
										{
										smeans[i-k] <- tempSmean
										}
									}							
								}
							}
						}
					}	
				}

			smeans <- c(wmeans$numeric_value.mean[1], smeans)	## pre-pend control mean to beginning of smoothed means
			wmeans <- cbind(wmeans, smeans)				## combine with means info

	    ##
			temp_value <- sprintf("%s.length",as.character(formula[[2]]))
			temp_idx4 <- which(colnames(wmeans) == temp_value)
			## simplify dof calcs ... if errors try old method above
			dof1 <- sum(wmeans[,5])
			dof2 <- nrow(wmeans)
			wmeans$dof <- (dof1 - dof2)
      formula_text<- sprintf("%s ~ %s",as.character(formula[[2]]),dose_name)
			se <- summaryBy(as.formula(formula_text), data=datatemp, FUN=c(sd,length))
			temp_value <- sprintf("%s.sd",as.character(formula[[2]]))
			temp_idx4 <- which(colnames(se) == temp_value)
			temp_value <- sprintf("%s.length",as.character(formula[[2]]))
			temp_idx5 <- which(colnames(se) == temp_value)
			se$sterr <- se[,temp_idx4]/(se[,temp_idx5]^.5)
			temp_idx4 <- which(colnames(se) == dose_name)
			temp_idx5 <- which(colnames(datatemp) == dose_name)
			mse <- 0
			for(i in 1:nrow(se)){
			  tempS = sum(datatemp[,temp_idx5] == se[i,temp_idx4])
				temp <- se$sterr[i]^2 * tempS * (tempS - 1)
				mse <- mse + temp
			}

			mse <- mse / (nrow(datatemp) - length(unique(datatemp[,temp_idx5])))

			## create WILLIAMS TEST STATISTIC
			willStat <- NULL
			
			temp_value <- sprintf("%s.length",as.character(formula[[2]]))
			temp_idx4 <- which(colnames(wmeans) == temp_value)
			
			for(i in 2:nrow(wmeans)){
				control_num  <- wmeans[1,temp_idx4]
				control_mean <- wmeans$smeans[1]
				test_num     <- wmeans[i,temp_idx4]
				test_mean    <- wmeans$smeans[i]

				willstat <- ifelse(direction=='decreasing', (control_mean - test_mean) / ((mse*((1/test_num) + (1/control_num)))^.5), (test_mean - control_mean) / ((mse*((1/test_num) + (1/control_num)))^.5))
				willStat <- c(willStat, willstat)
			}

			willStat <- c('control', willStat)
			wmeans$willStat <-  willStat
      wmeans <- .compute_crit_williams(wmeans,dose_name = dose_name,formulaV = as.character(formula[[2]])) # compute the critical values for the Williams Trend Test
			will_results <- rbind(will_results, wmeans)	
			}

	  twill = which(!(colnames(william) %in% c("tau","pvalue","mult_comp_test")))
	  idx <- c(which(colnames(will_results)%in% colnames(william)[twill] ),which(colnames(will_results) == "dose"))
	  #way too complicated loop to do something simple *sigh*
	  for (ii in length(idx):1){
	    reorder <- sort(will_results[,idx[ii]],index=TRUE)$ix
	    will_results <- will_results[reorder,]
	  }
	
		## ----------------------------------------------------------------------
		## convert williams statistic into p-value based on SAS crit levels
		## read in critical tables
		## ----------------------------------------------------------------------
		## remove control
		will_results2 <- subset(will_results, willStat != 'control')	

		will_results2$willStat <- as.numeric(as.character(will_results2$willStat)) 
		will_results2$crit05   <- as.numeric(as.character(will_results2$crit05))
		will_results2$crit01   <- as.numeric(as.character(will_results2$crit01))
   
		## determine how many asterisks each row deserves
		will_results2$mult_comp_signif <- ifelse(will_results2$willStat >= will_results2$crit01, 2,
										  ifelse(will_results2$willStat >= will_results2$crit05, 1, 0))

		will_results2$mult_comp_test <- 'WILLIAMS'
		ta1 = sprintf("%s.mean",as.character(formula[[2]]))
		ta2 = sprintf("%s.length",as.character(formula[[2]]))
		t_idx = which(!(colnames(will_results2) %in% c(ta1,ta2,"crit05","crit01","smeans","dof")))
		will_results2 = will_results2[,t_idx]

		}
	
	return(will_results2)
}
	
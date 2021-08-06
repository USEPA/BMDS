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

## ------------------
## SHIRLEY CODE
## ------------------

shirley_ntp <- function(jonck_data, shirl_data)
	{

	## loop through all groups flagged as SHIRLEY in jonck
  jonck_data =  jonckeere_ntp( formula,dose_name = dose_name,
                               data = data,pair="Shirley")
	shirley_results <- NULL

	temp_colnames <- unlist(c(unlist(colnames(dunn)),dose_name,as.character(unlist(formula[[2]]))))
	temp <-   colnames(data) %in% temp_colnames
	temp_d <- as.data.frame(data[,temp==TRUE])
	
	if(nrow(shirley) > 0)
		{
		for(s in 1:nrow(shirley))
			{
		  
		  temp_names <- rep(NA,ncol(shirley))
		  for (j in 1:length(temp_names)){
		    temp_names[j] <- unlist(shirley[d,j])
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
				exLoop <- ex
				for(g in (length(unique(ex[,dose_idx]))-1):1 )	## run once for each dose/group
					{
					## get ties for use in correction formula
					ties <- as.data.frame(table(table(exLoop[,nval_idx])))

					if(nrow(ties) > 1)
						{
						names(ties) <- c('numTies','count')
						ties$numTies <- as.integer(as.character(ties$numTies))

						## remove singletons
						ties <- subset(ties, numTies != 1)

						if(nrow(ties) > 0)
							{
							correction <- 0
							for(i in 1:nrow(ties))
								{
								correction <- correction + (ties$count[i] * (ties$numTies[i]^3 - ties$numTies[i]))
								}

							correction <- correction / (12 * (nrow(exLoop) - 1))
							} else
							{ 
							correction <- 0
							}
						} else 
						{
						correction <- 0
						}

					## rankings ... lowest value gets rank of 1
					ranks <- as.data.frame(table(exLoop[,nval_idx]))
					names(ranks) <- c('value', 'count')
					ranks$value <- as.numeric(as.character(ranks$value))

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
					names(ranks2) <- c('value', 'count', 'rank')

					## populate exLoop with ranks 
					## need to recast as numeric, may crash otherwise (weird!)
					exLoop[,nval_idx] <- as.character(exLoop[,nval_idx])
					exLoop[,nval_idx] <- as.numeric(exLoop[,nval_idx])

					exLoop[,nval_idx] <- merge(exLoop, ranks2, by.x=as.character(unlist(formula[[2]])), by.y='value', all.x=TRUE)
					exLoop2 <- exLoop2[order(exLoop2[ ,dose_idx]),]
					
					formula_text  = sprintf( "rank ~ %s" ,dose_name)
					
					rankSums <- summaryBy(rank ~ dose2 + dose, data=exLoop2, FUN = c(sum, length))		##!!
					rankMeans <- summaryBy(rank ~ dose2 + dose, data=exLoop2, FUN = c(mean, length))	## raw means, will overwrite with corrected means after next step

					mean_line <- NULL
					numer <- 0
					denom <- 0
					for(i in nrow(rankSums):1)
						{
						numer <- numer + rankSums$rank.sum[i]
						denom <- denom + rankSums$rank.length[i]

						if(i == 1)
							{
							mean_temp <- rankSums$rank.sum[i] / rankSums$rank.length[i]
							} else
							{
							mean_temp <- numer / denom
							}

						mean_line <- c(mean_line, mean_temp)
						}

					rankMeans$rank.mean <- rev(mean_line)
          ###
					###
					### STOP HERE FOR THE NIGHT 
					###
					###
					## find test statistic
					V  <- (nrow(exLoop) * (nrow(exLoop) + 1))/12 - correction
					Ri <- nrow(subset(exLoop, dose2==rankMeans$dose2[nrow(rankMeans)]))
					C  <- nrow(subset(exLoop, dose2==rankMeans$dose2[1]))

					if(shirley$tau[s] >= 0)
						{
						dosemean <- max(rankMeans$rank.mean[2:nrow(rankMeans)])
						shrl_num <- dosemean - rankMeans$rank.mean[1]
						} else
						{
						dosemean <- min(rankMeans$rank.mean[2:nrow(rankMeans)])
						shrl_num <- rankMeans$rank.mean[1] - dosemean
						}

					T <- shrl_num * (V * (1/Ri + 1/C))^-.5		## shirlstat in SAS code

					testStats <- c(testStats, T)
					doseCount <- c(doseCount, g)
					dose  <- c(dose, unique(ex$dose)[g+1])
					dose2 <- c(dose2, unique(ex$dose2)[g+1])
					num  <- c(num, nrow(subset(ex, dose2==unique(ex$dose2)[g+1])))
					sex  <- c(sex, shirley$sex[s])

					## remove latest dosage before returning to top of loop
					exLoop <- subset(exLoop, dose2 != unique(exLoop$dose2)[length(unique(exLoop$dose2))])
					}

				results <- as.data.frame(cbind(sex, dose, dose2, num, doseCount, testStats))		

				results$endpoint <- point
				results$phase_type <- phase
				results$phase_time <- time
				results$phase_start <- p_start
				results$phase_end <- p_end
				results$generation <- gen
				results$selection <- sel
				results$litter_name <- lit

				results$dose2 <- as.numeric(as.character(results$dose2))
				results$testStats <- as.numeric(as.character(results$testStats))	
				results$doseCount <- as.integer(as.character(results$doseCount))
				results$num <- as.integer(as.character(results$num))
				results$phase_time <- as.integer(as.character(results$phase_time))
				results$phase_start <- as.integer(as.character(results$phase_start))
				results$phase_end <- as.integer(as.character(results$phase_end))


				## add SAS crit values
				C01 <- c(0, 2.575, 2.607, 2.615, 2.618, 2.620, 2.621, 2.622)
				C05 <- c(0, 1.96, 2.015, 2.032, 2.040, 2.044, 2.047, 2.0485)

				B01 <- c(0, 0, 3, 4, 4, 4, 4, 4)
				B05 <- c(0, 0, 3, 4, 5, 6, 6, 6)


				## find crit values, generate number of stars
				results$mult_comp_signif <- 0
				nonsignif_flag <- 'NO'			## to match Laura's version force all doses after first non-signif dose to zero
				if(length(doseCount) <= 7)		## cannot handle more than 7 non control groups ... no crit values to compare to
					{
					for(i in 1:nrow(results))
						{
						if(nonsignif_flag == 'NO')
							{
							dosenum <- results$doseCount[i] + 1
							crit05 <- C05[dosenum] - ( (B05[dosenum] / 100) * (1 - (results$num[i] / results$num[nrow(results)]) ) )	
							crit01 <- C01[dosenum] - ( (B01[dosenum] / 100) * (1 - (results$num[i] / results$num[nrow(results)]) ) )	
							if(results$testStats[i] >= crit01)
								{
								results$mult_comp_signif[i] <- 2
								} else
								{
							if(results$testStats[i] >= crit05)
								{
								results$mult_comp_signif[i] <- 1
								} else
								{
								results$mult_comp_signif[i] <- 0
								nonsignif_flag <- 'YES'
								}
								}
							}
						}
					} else
					{
					# results$mult_comp_signif <- 'NA'
					results$mult_comp_signif <- ''
					}
				if(!is.null(shirley_results)){
				  shirley_results <- rbind(shirley_results, results)
				}
				else
          shirley_results <- results
				}
			}
		## check for existence
		if(!is.null(shirley_results))
			{

			shirley_results <- shirley_results[,c('sex', 'endpoint', 'generation', 'selection', 'litter_name', 'phase_type', 'phase_time', 'phase_start', 'phase_end', 'dose', 'dose2', 'mult_comp_signif')] 

			## coerce into same format as SHIRLEY results
			shirley_results$mult_comp_test <- 'SHIRLEY'
			shirley_results <- shirley_results[order(shirley_results$sex, shirley_results$endpoint, shirley_results$generation, shirley_results$selection, shirley_results$litter_name, shirley_results$phase_type, shirley_results$phase_time, shirley_results$phase_start, shirley_results$phase_end, shirley_results$dose2),]										  

			## remove NA's
			shirley_results[is.na(shirley_results)] <- ''
			}

		}

	return(shirley_results)
	}


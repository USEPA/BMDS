library(ToxicR)
jonck_input.shirley_dunn = read.csv("jonck_input.csv", header=TRUE);

A  = jonckeere_ntp(numeric_value ~ sex + phase_time + endpoint,data=jonck_input.shirley_dunn,pair="Shirley")
B  = dunn_ntp(numeric_value ~ sex + endpoint,data=jonck_input.shirley_dunn)
dd = shirley_ntp(numeric_value ~ sex +  endpoint, data=jonck_input.shirley_dunn)
BB = dunnett_ntp(numeric_value ~ sex + endpoint,data=jonck_input.shirley_dunn)
CC = williams_ntp(numeric_value ~ sex  + endpoint,data=jonck_input.shirley_dunn)



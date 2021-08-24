library(ToxicR)
jonck_input.shirley_dunn = read.csv("jonck_input.csv", header=TRUE);

A  = jonckeere_ntp(numeric_value ~ sex + endpoint,data=jonck_input.shirley_dunn)
B  = dunn_ntp(numeric_value ~ sex + endpoint,data=jonck_input.shirley_dunn)
AA = shirley_ntp(numeric_value ~ sex + endpoint,data=jonck_input.shirley_dunn)
BB = dunnett_ntp(numeric_value ~ sex + endpoint,data=jonck_input.shirley_dunn)
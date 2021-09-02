library(ToxicR)
jonck_input.shirley_dunn = read.csv("jonck_input.csv", header=TRUE);

A  = NTP$jonckeere(numeric_value ~ sex + phase_time + endpoint,data=jonck_input.shirley_dunn,pair="Shirley")
B  = NTP$dunn(numeric_value ~ sex + endpoint,data=jonck_input.shirley_dunn)
dd = NTP$shirley(numeric_value ~ sex +  endpoint, data=jonck_input.shirley_dunn)
BB = NTP$dunnett(numeric_value ~ sex + endpoint,data=jonck_input.shirley_dunn)
CC = NTP$williams(numeric_value ~ sex  + endpoint,data=jonck_input.shirley_dunn)



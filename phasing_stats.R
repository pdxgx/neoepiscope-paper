
summary33 <- read.table('phasing_summary_33.tsv', sep='\t', header = T)
summary33$Percent_phased_germline <- (summary33$Phased_germline/summary33$Total_somatic_mutations)*100
summary33$Percent_phased_somatic <- (summary33$Phased_somatic/summary33$Total_somatic_mutations)*100
summary33$Percent_phased_combined <- (summary33$Combined_phasing/summary33$Total_somatic_mutations)*100
summary33$Percent_FS <- (summary33$Frameshift_affects/summary33$Total_somatic_mutations)*100

summary72 <- read.table('phasing_summary_72.tsv', sep='\t', header = T)
summary72$Percent_phased_germline <- (summary72$Phased_germline/summary72$Total_somatic_mutations)*100
summary72$Percent_phased_somatic <- (summary72$Phased_somatic/summary72$Total_somatic_mutations)*100
summary72$Percent_phased_combined <- (summary72$Combined_phasing/summary72$Total_somatic_mutations)*100
summary72$Percent_FS <- (summary72$Frameshift_affects/summary72$Total_somatic_mutations)*100

# Phasing prevalence at 33 bp distance
mean(summary33$Percent_phased_germline)
max(summary33$Percent_phased_germline)
mean(summary33$Percent_phased_somatic)
max(summary33$Percent_phased_somatic)
mean(summary33$Percent_FS)
max(summary33$Percent_FS)

# Phasing prevalence at 72 bp 
mean(summary72$Percent_phased_germline)
max(summary72$Percent_phased_germline)
mean(summary72$Percent_phased_somatic)
max(summary72$Percent_phased_somatic)
mean(summary72$Percent_FS)
max(summary72$Percent_FS)

par(mar=c(4,6,0.5,0.5))
boxplot(summary72$Percent_phased_germline, 
        summary72$Percent_phased_somatic, 
        summary72$Percent_phased_combined, col = c('blue', 'red', 'dimgray'),
        names = c('Somatic-Germline', 'Somatic-Somatic', 'Combined'),
        ylab='Phased somatic variants (% per patient)', xlab='', cex.lab=2, cex.axis=1.75, boxlwd=3)


# Distance analysis

all_mutations <- sum(summary72$Total_somatic_mutations)


germline <- read.table('germline_distances.tsv', sep='\t', header = T)
germline <- germline[order(germline$Distance),]
germline_prop <- rep(0, length(germline$Distance))
for (i in 1:length(germline$Distance)){
  count <- i
  proportion <- i/all_mutations
  germline_prop[i] <- proportion
}
germline$Proportion <- germline_prop

somatic <- read.table('somatic_distances.tsv', sep='\t', header = T)
somatic <- somatic[order(somatic$Distance),]
somatic_prop <- rep(0, length(somatic$Distance))
for (i in 1:length(somatic$Distance)){
  count <- i
  proportion <- i/all_mutations
  somatic_prop[i] <- proportion
}
somatic$Proportion <- somatic_prop

combined <- read.table('combined_distances.tsv', sep='\t', header = T)
combined <- combined[order(combined$Distance),]
combined_prop <- rep(0, length(combined$Distance))
for (i in 1:length(combined$Distance)){
  count <- i
  proportion <- i/all_mutations
  combined_prop[i] <- proportion
}
combined$Proportion <- combined_prop

par(mar=c(6,6,0.5,0.5))
plot(c(-10,-15), c(-10,-11), xlim = c(0, 100), ylim = c(-0.05, 1.25),
     xlab="Minimum distance to closest variant (bp)", 
     ylab="Cumulative somatic variants (% of total)", col='black', cex.lab=1.75, cex.axis=1.5)
rect(24, -1, 33, 4.75, col='gainsboro', border=NA)
rect(36, -1, 72, 4.75, col='gainsboro', border=NA)
lines(germline$Distance[germline$Distance < 100], 
      100*germline$Proportion[germline$Distance < 100], type="l", col='blue')
lines(somatic$Distance[somatic$Distance < 100], 
      100*somatic$Proportion[somatic$Distance < 100], type="l", col='red')
lines(combined$Distance[combined$Distance < 100], 
      100*combined$Proportion[combined$Distance < 100], type="l", col='black')
rect(xleft=74, xright=77, ytop=0.23, ybottom=0.17, col="black", border="black")
rect(xleft=74, xright=77, ytop=0.13, ybottom=0.07, col="blue", border="black")
rect(xleft=74, xright=77, ytop=0.03, ybottom=-0.03, col="red", border="black")
text(84.25, 0.20, labels="Combined",cex=1.5, col="black")
text(89.5, 0.10, labels="Somatic-Germline",cex=1.5, col="black")
text(89, 0.00, labels="Somatic-Somatic",cex=1.5, col="black")
box(lty="solid") 



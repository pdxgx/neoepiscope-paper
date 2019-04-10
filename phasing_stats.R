
summary33 <- read.table('phasing_summary_33.tsv', sep='\t', header = T)
summary33$Percent_phased_germline <- (summary33$Phased_germline_tx/summary33$Total_somatic_mutations)*100
summary33$Percent_phased_somatic <- (summary33$Phased_somatic_tx/summary33$Total_somatic_mutations)*100
summary33$Percent_phased_combined <- (summary33$Combined_phasing_tx/summary33$Total_somatic_mutations)*100
summary33$Percent_FS <- (summary33$Frameshift_phasing/summary33$Total_somatic_mutations)*100
summary33$Percent_nonstop <- (summary33$Nonstop_phasing/summary33$Total_somatic_mutations)*100

summary72 <- read.table('phasing_summary_72.tsv', sep='\t', header = T)
summary72$Percent_phased_germline <- (summary72$Phased_germline_tx/summary72$Total_somatic_mutations)*100
summary72$Percent_phased_somatic <- (summary72$Phased_somatic_tx/summary72$Total_somatic_mutations)*100
summary72$Percent_phased_combined <- (summary72$Combined_phasing_tx/summary72$Total_somatic_mutations)*100
summary72$Percent_FS <- (summary72$Frameshift_phasing/summary72$Total_somatic_mutations)*100
summary72$Percent_nonstop <- (summary72$Nonstop_phasing/summary72$Total_somatic_mutations)*100

summary94 <- read.table('phasing_summary_94.tsv', sep='\t', header = T)
summary94$Percent_phased_germline <- (summary94$Phased_germline_tx/summary94$Total_somatic_mutations)*100
summary94$Percent_phased_somatic <- (summary94$Phased_somatic_tx/summary94$Total_somatic_mutations)*100
summary94$Percent_phased_combined <- (summary94$Combined_phasing_tx/summary94$Total_somatic_mutations)*100
summary94$Percent_FS <- (summary94$Frameshift_phasing/summary94$Total_somatic_mutations)*100
summary94$Percent_nonstop <- (summary94$Nonstop_phasing/summary94$Total_somatic_mutations)*100

# Phasing prevalence at 33 bp distance
mean(summary33$Percent_phased_germline)
min(summary33$Percent_phased_germline)
max(summary33$Percent_phased_germline)
mean(summary33$Percent_phased_somatic)
min(summary33$Percent_phased_somatic)
max(summary33$Percent_phased_somatic)

# Phasing prevalence at 72 bp 
mean(summary72$Percent_phased_germline)
min(summary72$Percent_phased_germline)
max(summary72$Percent_phased_germline)
mean(summary72$Percent_phased_somatic)
min(summary72$Percent_phased_somatic)
max(summary72$Percent_phased_somatic)

# Phasing prevalence at 94 bp
mean(summary94$Percent_FS)
min(summary94$Percent_FS)
max(summary94$Percent_FS)
mean(summary94$Percent_nonstop)
min(summary94$Percent_nonstop)
max(summary94$Percent_nonstop)

# Plot phasing prevalence
par(mar=c(4,6,0.5,0.5))
boxplot(summary72$Percent_phased_germline, 
        summary72$Percent_phased_somatic, 
        summary72$Percent_phased_combined, col = c('blue', 'red', 'dimgray'),
        names = c('Somatic-Germline', 'Somatic-Somatic', 'Combined'),
        ylab='Phased somatic variants (% per patient)', xlab='', cex.lab=2, cex.axis=1.75, boxlwd=3)

# Phasing prevalence by disease
adj_disease <- rep('melanoma', length(summary72$Disease))
for (i in 1:length(summary72$Disease)){
  if (immunorx_summary$Disease[i] == 'NSCLC'){
    adj_disease[i] <- 'NSCLC'
  } else if (immunorx_summary$Disease[i] %in% c('colon', 'endometrial', 'thyroid')){
    adj_disease[i] <- 'MMR_deficient'
  }
}
summary72$Adjusted_disease <- adj_disease
summary(aov(summary72$Percent_phased_combined~summary72$Adjusted_disease))

par(mar=c(4.5,5.5,0.5,0.5))
boxplot(summary72$Percent_phased_combined[summary72$Adjusted_disease == 'MMR_deficient'], 
        summary72$Percent_phased_combined[summary72$Adjusted_disease == 'NSCLC'],
        summary72$Percent_phased_combined[summary72$Adjusted_disease == 'melanoma'], 
        col = 'gray',
        names = c('MMR-deficient', 'NSCLC', 'Melanoma'),
        ylab='Phased somatic variants (% per patient)', xlab='Disease type', 
        cex.lab=2, cex.axis=1.75, boxlwd=3)

# Phasing prevalence by coverage-adjusted mutation burden
coverage_summary <- read.table('coverage_summary.tsv', sep='\t', header=T)
cov <- rep(0, length(summary72$Patient))
for (i in 1:length(summary72$Patient)){
  patient <- as.character(summary72$Patient[i])
  cov[i] <- coverage_summary$Mbp_coverage[as.character(coverage_summary$Patient) == patient]
}
summary72$Mbp_covered <- cov

x <- summary72$Total_somatic_mutations/summary72$Mbp_covered
y <- summary72$Percent_phased_combined_tx
summary(lm(y~0+x))
cor.test(summary72$Total_somatic_mutations/summary72$Mbp_covered, summary72$Percent_phased_combined
         )

par(mar=c(4.5,5,0.5,1))
plot(x, y, 
     xlab="Mutational burden (# somatic variants/Mbp genome covered)", ylab="Phased somatic variants (%)", 
     cex.lab=2, cex.axis=1.75, pch=19, cex=0.5)
abline(lm(y~0+x), col='red')
text(175, 12, labels="y = 0.088257x",cex=1.5, col="black")

# False positive/false negative epitopes when not accounting for phasing
eps <- read.table('phasing_epitope_data.tsv', header=T, sep='\t')
eps$False_pos <- eps$Unphased_only/(eps$Shared+eps$Unphased_only)
eps$False_neg <- eps$Phased_only/(eps$Shared+eps$Unphased_only)

mean(eps$False_pos)
median(eps$False_pos)
min(eps$False_pos)
max(eps$False_pos)
quantile(eps$False_pos, c(0.25, 0.75))

mean(eps$False_neg)
median(eps$False_neg)
min(eps$False_neg)
max(eps$False_neg)
quantile(eps$False_neg, c(0.25, 0.75))

# Distance analysis

all_mutations <- sum(summary72$Total_somatic_mutations)


germline <- read.table('germline_distances_tx.tsv', sep='\t', header = T)
germline <- germline[order(germline$Distance),]
germline_prop <- rep(0, length(germline$Distance))
for (i in 1:length(germline$Distance)){
  count <- i
  proportion <- i/all_mutations
  germline_prop[i] <- proportion
}
germline$Proportion <- germline_prop

somatic <- read.table('somatic_distances_tx.tsv', sep='\t', header = T)
somatic <- somatic[order(somatic$Distance),]
somatic_prop <- rep(0, length(somatic$Distance))
for (i in 1:length(somatic$Distance)){
  count <- i
  proportion <- i/all_mutations
  somatic_prop[i] <- proportion
}
somatic$Proportion <- somatic_prop

combined <- read.table('combined_distances_tx.tsv', sep='\t', header = T)
combined <- combined[order(combined$Distance),]
combined_prop <- rep(0, length(combined$Distance))
for (i in 1:length(combined$Distance)){
  count <- i
  proportion <- i/all_mutations
  combined_prop[i] <- proportion
}
combined$Proportion <- combined_prop

par(mar=c(6,6,0.5,0.5))
plot(c(-10,-15), c(-10,-11), xlim = c(0, 100), ylim = c(-0.05, 4),
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

## RNA-seq coverage

rna_support <- read.table('paired_rna_support_summary.tsv', header=T, sep='\t')
rna_support$Total_phased_pairs <- rna_support$Germline_phased_pairs + rna_support$Somatic_phased_pairs

mean(100*(rna_support$Covered_germline_pairs[rna_support$Germline_phased_pairs != 0]/rna_support$Germline_phased_pairs[rna_support$Germline_phased_pairs != 0]))
mean(100*(rna_support$Covered_somatic_pairs[rna_support$Somatic_phased_pairs != 0]/rna_support$Somatic_phased_pairs[rna_support$Somatic_phased_pairs != 0]))

mean(100*(rna_support$Supported_across_exon_pairs[rna_support$Covered_across_exon_pairs != 0]/rna_support$Covered_across_exon_pairs[rna_support$Covered_across_exon_pairs != 0]))
mean(100*(rna_support$Not_supported_across_exons[rna_support$Covered_across_exon_pairs > 0]/rna_support$Covered_across_exon_pairs[rna_support$Covered_across_exon_pairs > 0]))
mean(100*(rna_support$Novel_across_exon_pairs[rna_support$Total_phased_pairs != 0]/rna_support$Total_phased_pairs[rna_support$Total_phased_pairs != 0]))

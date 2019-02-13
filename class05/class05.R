#Class 05 Graphics and plots with R :https://bioboot.github.io/bggn213_W19/class-material/lab-5-bggn213.html 
#Section 2A: Line Plot 

read.table("bimm143_05_rstats/weight_chart.txt")
weight <- read.table("bimm143_05_rstats/weight_chart.txt", header = TRUE)

plot(weight, type = "b", pch = 1:3, lwd = 2, cex = 1.5, ylim = c(1,10), xlab = "Age (months)", ylab = "Weight (kg)", main = "Fat Babies", col = c("green", "Black"))


#Section 2B : Barplot
feat <- read.table("bimm143_05_rstats/feature_counts.txt", sep = "\t")

#barplot(feat$Count, names.arg = feat$Feature, las = 2, horiz = TRUE)


#par for parameters
#par(mar = c(5,10,4,2))
#barplot(feat$Count, names.arg = feat$Feature, las = 2, horiz = TRUE)
#Section 3 : Using color in plots
#Section 3A : Providing Color Vectors
#read.table("bimm143_05_rstats/male_female_counts.txt", header = TRUE, sep = "\t")
#read.table first, to see how it's separated!!!!
mfc <- read.table("bimm143_05_rstats/male_female_counts.txt", header = TRUE, sep = "\t")

#nrow(x) to find out # of rows, then can use in rainbow() bc assign exact amount of colors to generate
#colors() to find out color names, can set using name/#
#rainbow(), and various other coloring options
#barplot(mfc$Count, names.arg = mfc$Sample, las = 2, col = c("cornsilk", "aliceblue"))

barplot(mfc$Count, names.arg = mfc$Sample, las = 2, col = rainbow(nrow(mfc)))

#Section 3B
exp <- read.delim("bimm143_05_rstats/up_down_expression.txt")
table(exp$State)
 
plot(exp$Condition1, exp$Condition2, col = exp$State)
#will plot x as condition1, y as condition 2, color it based on State variable
plot(exp$Condition1, exp$Condition2, col = exp$State, main = "conditions?", xlab = "Condition 1", ylab = "Condition 2")


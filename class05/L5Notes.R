#Lecture 5
#rnorm(x) will generate x numbers in a normal distribution
rnorm(1000, mean = 14, sd = 1)
#functions learned
plot(x)
hist(x)
list(x)
summary(x)
# in functions(x), only x is default/required; all the rest with numrals = s;klf are optional. When you see a = sign in the (), it's optional
#command-enter will put it into the console and run it 
#file needs to be in same folder as the Project file. 
#define path by listing folder name/nameoffile
read.table("bimm143_05_rstats/weight_chart.txt")
#at this point table isn't in the environment yet, it's just read out to us
#set it to the environment
weight <- read.table("bimm143_05_rstats/weight_chart.txt")
# weight <- read.table("bimm143_05_rstats/weight_chart.txt", header = TRUE)
weight <- read.table("bimm143_05_rstats/weight_chart.txt", header = TRUE)
#setting header = TRUE in argument will make Row 1 as the header of the columns; header = false is the default
plot(weight, type = "b", pch = 15, cex = 1.5, ylim = c(1,10), xlab = "Age (months)", ylab = "Weight (kg)", main = "Fat Babies")
#to change any of the variables, go into ?plot to look at what the arguments mean. 
#pch = plot character; type = type of graph, cex = character size, lwd = line thickness
plot(weight$Age, weight$Weight, pch = 1:3, col = "red")
#by doing col = c("red", "blue") will make it alternate colors; and pch = 1:3 will alternate the various data point figures
#by using plot(weight$Age, weight$Weight) , it will pull the x values from the weight data, and the y values from the weight data; j ust remember to change the axis names; this is different from directly using plot(weight) bc you can pull x/y values from different sets of data

#2B: barplot
#need to check the type of file/data in order to call it in correctly. 
# in read.table() : arguments : header, sep, quote, 
#sep sets the field separation character, if you put sep = "\t", it will set it to tab
#if you look closely in the feature_counts file, the space betwen the columns title is actually a tab, so that is how we will separate the columns
#par() to view parameters
#To change parameters: do par( whatever parameter thing), but you have to replot the graph AFTER you call the parameter change so that your graph is changed as well with the new parameters. 

feat <- read.table("bimm143_05_rstats/feature_counts.txt", header = TRUE, sep = "\t" )
barplot(feat$Count, names.arg = feat$Feature, las = 2, horiz = TRUE)
#par for parameters

par(mar = c(5, 10, 4, 2))
barplot(feat$Count, names.arg = feat$Feature, las = 2, horiz = TRUE)


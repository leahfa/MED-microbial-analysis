##########################################################################################################################
#First run scripts data_prepration_and variable_definition, and install lmic package
# results can be plotted with heatmap_of_correlation_results in `plotting_scripts` directory
##########################################################################################################################

library(lmic)

# Make Summary  tables of 1dim (correlate many features against one feature ) or 
# 2dim (many features against many features) , using Spearmans correlations

#### prepare data ####
#verify metadata (df) adn microbiome data (dat) are perfectly matched:
identical(rownames(dat), df$ID)
#Filter rare taxa from dat - crucial to reduce spurious correlations! Preferable to use the version
#with shortened taxa names (dat.short):
dat.f<-filtCols(dat.short,10)
dim(dat.f)

#### Correlation examples ####

#### cor.1dim ####
# cor.1dim takes a matrix/dataframe of features (microbiome data for example) and correlates 
#afainst a single vector. The last argument sets number of digits to round to:

#Example1: correlate genera againt daily score, at all timepoints:
res1<-cor.1dim(dat.f, df$Daily_score,4)

#### cor.2dim ####
# cor.2dim takes 2 matrices/dataframes and correlates all features against each other. Best to 
# set the matrix which has fewer features as the first argument (this means FDR correction will be
# run across all features of teh second object, making for a more stringent correction)
# The option 'dig' sets number of digits to round to, and "cut.level" controls whether or not the taxa name shsould
# be stripped to genus level(takes "yes" or "no", default is "no)

#Example 1:
#correlate genera againt all dietary variables, at all timepoints:
res4<-cor.2dim(df[ ,diet.cols],filtCols(dat.f,10),dig=4, cut.level="no")

#Example 2:
#correlate genera againt all dietary variables and all clinical variables, at all timepoints:
res4<-cor.2dim(df[ ,c(diet.cols,clin.cols)], dat.f,4,"no")

#Example 3:
#correlate genera against specific variables of interest, at all timepoints:
res4<-cor.2dim( dat.f,df[ ,c("WBC","Fruit","age")],4,"no")


#Example 4:
#correlate genera against dietary variables, at timepoint1 only
# (after running script make_subsets_by_teimepoints):
dat.f<-filtCols(dat1.short,10)
identical(df.tp1$ID, rownames(dat.f)) #verify data is matched to metadata
res4<-cor.2dim(df.tp1[ ,diet.cols], dat.f,4,"no")
#this can also be done in a single line:
res4<-cor.2dim(df.tp1[ ,diet.cols],filtCols(dat.f,10),4,"no")


#Example5:
#Correlate data with full taxa names (use dat instead of dat.short):
res4<-cor.2dim(df[ ,diet.cols],filtCols(dat,10),4,"no")
# Run the same, but shorten taxa names down to genera level using the "cut.level" option of cor.2dim:
res4<-cor.2dim(df[ ,diet.cols], filtCols(dat,10),4,"yes")


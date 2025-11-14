
######BINF6210 - Assignment 1: Madelin Stueck ----

#Original Author: Thomas Furtado
#Edited by: Madelin Stueck
#added intermediate dataframe, investigated variables + cleaned up some plot visualization


rm(list = ls())

#### Section 1. Necessary Packages----

install.packages("knitr")
install.packages("gridExtra")
install.packages("ggridges")
install.packages("gridExtra")


library(ggridges)
library(gridExtra)
library(grid)
library(knitr)
library(tidyverse)
library(ggplot2) 
library(dplyr)
library(gridExtra)

#### Section 2. Data Import and inspection----

Tardi_data<- read_tsv("../Data/result.tsv")

summary(Tardi_data)
names(Tardi_data)


#Lets explore Realm, this seems interesting!
Tardi_data$realm

#First, I'd like to generate a graph to look at the most/least sampled Realms
Tardi_data %>%
  filter(!is.na(realm)) %>%
  count(realm) %>%
  ggplot(aes(x = reorder(realm, n), y = n)) +
  geom_col(fill = "skyblue") +
  coord_flip() +
  labs(x = "realm", y = "Number of specimens", title = "Specimens per Realm")

#Edit: We might also seek to create a new dataframe for this visualization, in case we want to go back and inspect this further. 
#Making intermediate dataframes can help collaborators follow your work and also verifies that data is being manipulated properly
#Perhaps something like:

realm_counts <- Tardi_data %>%
  filter(!is.na(realm)) %>% 
  count(realm, name = "# of Specimens") %>%
  arrange(desc(`# of Specimens`))


#I find that the realms "Palearctic" and "Nearctic" have the most samples.I will investigate some variable by these two realms.

Tardi_data$nuc_basecount
#Most other columns I see have alot of N/A values, but this column contains a lot of data, so I'll investigate this variable.

#Edit: it might be helpful to add some more information about this variable. What is nuc_basecount?
#Perhaps you could verify its identity by comparing it to the column "nuc"

#create a vector with the lengths of each value in column 'nuc'
nuc_len_vector <- nchar(gsub("[^ATGCN]", "", Tardi_data$nuc))

#compare it to column 'nuc_basecount'
identical(nuc_len_vector, Tardi_data$nuc_basecount)

##Why does nuc_basecount not = the sum of bases in the column 'nuc'? This could be interesting to explore. 


unique(Tardi_data$marker_code)
#I notice that there are many different types of markers here! 18S, D4L, COI etc. I must filter by one of them ,or any comparisons I make are meaningless! I will use COI-5P


#### Section 3. Data refinement----



#I will now filter by realms of interest, isolate COI-5P data, and remove missing nuc_basecount data. SR= sequence,realm

SR_data <- Tardi_data %>%
  filter(marker_code == "COI-5P",
         realm %in% c("Palearctic", "Nearctic"),
         !is.na(nuc_basecount)) %>%
  select(realm, nuc_basecount)

#Now, I would like to make a summary of this data, displaying desired info

seq_summary <- SR_data %>%
  group_by(realm) %>%
  summarise( count = n(), mean_length = mean(nuc_basecount), median_length = median(nuc_basecount), min_length = min(nuc_basecount), max_length = max(nuc_basecount), .groups = 'drop')

print(seq_summary)

#I notice that there is a large difference in min/max lengths, additionally the spread between the mean/median is quite large

#Upon researching I found that the approximate range of a COI seq is 200-900 bp (Guo et al., 2022). I will refine the data, so that it is more normalized and more likely to be free of errors/outliers.

#rm (seq_summary)
#rm (SR_data)

##Edit: you may not want to remove useful data summaries from your environment, especially since you changed the name in the refined version, keeping both could make for helpful comparisons later, and helps collaborators follow each stage in your investigation 



#refined version: I will use suffix (R) to denote this.

SR_dataR <- Tardi_data %>%
  filter(marker_code == "COI-5P",           
         realm %in% c("Palearctic", "Nearctic"),
         !is.na(nuc_basecount),             
         nuc_basecount >= 200,              
         nuc_basecount <= 900) %>%          
  select(realm, nuc_basecount)



seq_summaryR <- SR_dataR %>%
  group_by(realm) %>%
  summarise(Count = n(), "Mean Length" = mean(nuc_basecount), "Median Length" = median(nuc_basecount), "Min Length" = min(nuc_basecount), "Max Length" = max(nuc_basecount), .groups = 'drop') %>%
  rename( Realm = realm)

print(seq_summaryR)

#### Section 4. Data Exploration----


kable(seq_summaryR, caption = "Summary of Tardigrada COI Sequence Lengths (200–900 bp) Across Realms", digits = 1)

#This function only generates a table in text form, I'll leave it in for convenience down the road. But I will like my table to also exist in a figure format.

#grid.table( seq_summaryR, rows = NULL, cols = c("Realm", "Count", "Mean Length", "Median Length", "Min Length", "Max Length"))

#Edit: this is a fantastic table, but grid.table prints tables directly onto the existing graphics device: if you already have a plot open, the table will get drawn on top of it.
#For visualization purposes, you may want to start a new graphics page first, then draw the table.
#this outputs the figure formatted table on a blank page

grid.newpage()
grid.table(
  seq_summaryR,
  rows = NULL,
  cols = c("Realm", "Count", "Mean Length", "Median Length", "Min Length", "Max Length")
)


#Now that the summary data is shown in a table format, lets look into viewing the distribution of sequence lengths, and visually compare the two regions.

#I'd like to emphasize the mean between these two data sets, this is where the data is most interesting.

ggplot(SR_dataR, aes(x = realm, y = nuc_basecount, fill = realm)) +
  geom_violin(trim = FALSE) +
  stat_summary(fun = mean, geom = "point", color = "red", size = 4) +
  labs( title = "Distribution of COI Sequence Lengths Across Realms", x = "Realm", y = "Sequence Length (bp)") + 
  theme_minimal() +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

#yuzaR Data Science (2023, Nov. 24). Master Box-Violin Plots in {ggplot2} and Discover 10 Reasons Why They Are Useful
#[Video]. YouTube. https://www.youtube.com/watch?v=rvm94zcoKT0


#Edit: I love the violin plot! It may also be interesting to visualize the distributions as histograms, so that you can see the true count values on the y-axis. Here is a distribution histogram that allows collaborators to see the count values of sequence lengths for COI samples in the realms Nearctic and Palearctic. 


ggplot(SR_dataR, aes(x = nuc_basecount, fill = realm)) +
  geom_histogram(binwidth = 10, position = "dodge", color = "black") +
  labs(
    title = "Distribution of COI Sequence Lengths Across Realms",
    x = "Sequence Length (bp)",
    y = "Count"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

#### Section 5. Data Interpretation and Analyses----



#Now it's important to note that the width of the violin plot can be a little misleading, as it only shows the density proportionate to the group.

#I can perform an F-test to check if the variance of Nearctic is significantly smaller than Palearctic 

#Null hypothesis: the variances of the two groups are equal
#Alternative hypothesis: the variances are not equal

var.test(nuc_basecount ~ realm, data = SR_dataR)

#p-value of 0.6434, we cannot reject null hypothesis. Therefore the variance between the two is not statistically significant.

#“While the violin plot shows a high peak for Nearctic sequences, suggesting tight clustering, the F-test indicates that variance is not statistically different. Thus, we cannot claim that Nearctic sequences are truly more conserved in length than Palearctic sequences.”



#I noticed possibly distinct regions in the violin plot, 1 being around the mean obvously, but there were 2 others: 1 at lengths higher than the mean, and 1 at lengths below the mean. These slight bulges aren't so apparent with the violin structure. I'd like to view them in more detail with a ridge-plot.



ggplot(SR_dataR, aes(x = nuc_basecount, y = realm, fill = realm)) +
  geom_density_ridges(alpha = 0.7, scale = 1) +        
  geom_vline(xintercept = c(570, 660, 750), linetype = "dashed", color = "red", size = 0.8) +
  labs( title = "Ridgeline Plot of COI-5P Sequence Lengths by Realm", x = "Sequence Length (bp)", y = "Realm") +
  scale_fill_manual(values = c("Nearctic" = "skyblue", "Palearctic" = "orange")) +
  theme_minimal() +
  theme(legend.position = "none")

#The ridgeline plot shows that there are roughly 3-4 common lengths that Palearctic and Nearctic COI sequences cluter around. Nearctic peaks slightly higher on average. This is an interesting finding which could be a future avenue of investigation!

#While the distributions show multiple peaks, the F-test indicates no statistically significant difference in variance between realms (p = 0.64). 
#The minor differences in peaks can be attributed to sampling as well as different sequence processing techniques.

#The direct question being asked, as well as how the results pertain to answering it will be elaborated on in greater depth through the story board.




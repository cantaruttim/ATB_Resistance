
library(dplyr)
library(ggplot2)
library(AMR)
library(cleaner)

# install.packages(c("AMR", "cleaner"))

# Working with antimicrobial names or codes

as.ab("Amoxicillin")
# [1] AMX



as.ab("amox")


# Não precisamos importar o pacote, pois ele já vem com o AMR
View(antibiotics)
# View(microorganisms)

ab_name("J01CA04")

carbapenems(antibiotics)

# Working with antimicrobial susceptibility test results

as.rsi(as.mic(2), "E. coli", "ampicillin", guideline = "EUCAST 2021")

as.rsi(as.mic(32), "E. coli", "ampicillin", guideline = "EUCAST 2021")

# Working with defined daily doses (DDD)

atc_online_ddd("amoxicillin", administration = "O")

atc_online_groups("amoxicillin")

# Enhancing antimicrobial resistance data

#5. Analyzing antimicrobial resistance data

# 5.1. Calculation of antimicrobial resistance

example_isolates %>% summarize(
   r_gen = proportion_R(GEN), r_amx = proportion_R(AMX),
   n_gen = n_rsi(GEN), n_amx = n_rsi(AMX), n_total = n())


# EXAMPLE

library("dplyr")
library("tidyr")
library("AMR")

options(AMR_locale = "en")
data <- example_isolates_unclean
glimpse(data)


unique(data$hospital)
# [1] "A" "B" "C"

unique(data$bacteria)

#  [1] "E. coli"                  "K. pneumoniae"           
# [3] "S. aureus"                "S. pneumoniae" 
# [5] "klepne"                   "strpne"    
# [7] "esccol"                   "staaur"    
# [9] "Escherichia coli"         "Staphylococcus aureus"  
# [11] "Streptococcus pneumoniae" "Klebsiella pneumoniae"


data %>% count(bacteria)
# ?mutate

data <- data %>% mutate(bacteria = as.mo(bacteria), bacteria_name = mo_name(bacteria))
mo_uncertainties()

data %>% count(bacteria, bacteria_name)

# In a next step, we can further enrich the data with additional microbial taxonomic data based
# on the bacteria variable, such as Gram-stain and microorganism family

data <- data %>% mutate(
   gram_stain = mo_gramstain(bacteria), family = mo_family(bacteria))

data %>% count(gram_stain)

data %>% count(family)

ab_info("GEN")


# In a data set containing antimicrobial names or codes (e.g., antimicrobial prescription data),
# the as.ab() function can be used to transform all values to valid antimicrobial codes. Extra
# columns with the official name and the defined daily dose (DDD) for intravenous administration could be added using ab_name() and ab_ddd().

ab_example <- data.frame(agents = c("AMX", "Ceftriaxon", "Cipro"))
ab_example %>% mutate(
 agents = as.ab(agents), agent_names = ab_name(agents),
 ddd_iv = ab_ddd(agents, administration = "iv"),
 ddd_oral = ab_ddd(agents, administration = "oral"))



data %>%
 select(AMX:GEN) %>%
 pivot_longer(everything(), names_to = "antimicrobials",
 values_to = "interpretation") %>%
 count(interpretation)

# 1 < 0.5 S          143
# 2 I               1105
# 3 R               4607
# 4 S               6145


data <- data %>%mutate_at(vars(AMX:GEN), as.rsi)

data %>%
 select(AMX:GEN) %>%
 pivot_longer(everything(), names_to = "antimicrobials",
 values_to = "interpretation") %>%
 count(interpretation)


data <- data %>% eucast_rules()


data <- data %>% mutate(mdro = mdro(., guideline = "nl"))
data %>% count(bacteria_name, mdro)



resistance_proportion <- data %>%
 filter_first_isolate() %>%
 group_by(hospital) %>%
 proportion_df()


head(resistance_proportion)

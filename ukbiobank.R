
library(ukbtools)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(scales)
library(tidyr)

setwd("C:/Users/bossarr/OneDrive - Albany Medical Center/Biobank_R/ukbiobank_rosacea/fileset")

my_ukb_data <- ukb_df("ukb42708", path = "~/Desktop/Alex_Albany/ukbiobank_rosacea/fileset")

my_ukb_key <- ukb_df_field("ukb42708", path = "~/Desktop/Alex_Albany/ukbiobank_rosacea/fileset")
df_field <- ukb_df_field("ukb42708", path = "~/Desktop/Alex_Albany/ukbiobank_rosacea/fileset")
df <- ukb_df("ukb42708", path = "~/Desktop/Alex_Albany/ukbiobank_rosacea/fileset") #takes two arguments, the stem of your fileset and the path, and returns a dataframe with usable column names.


#Salvar e deletar objetos
save(my_ukb_data, file = "my_ukb_data.rda") # Salvar um unico objeto
save(df, df_field, my_ukb_key, x_data_sub, file = "ukb_rosacea.RData") #Salvar diversos objetos
save.image(file = "ukdata_051821.RData")


load("my_ukb_data.rda")
load("ukb_rosacea.RData")
load("ukdata_051821.RData")

remove(list=c("eye", "lps_ttr_huvec_022421", "sample", "test_data", "test_filter", "test_filter_multiple", "test_filter_multiple_keratopahty", "x_data_sub", "y_data_sub")) # limpar data no Environment

### Selecionar as 10 primeiras 

str(x_data_sub)#check the samples
y_data_sub<- sample_n(my_ukb_data,500) #random 500 rows

#Buscar ICD codes
ukb_icd_keyword("eye", icd.version = 10)

#retrieve the "meaning" of an ICD code 
ukb_icd_code_meaning(icd.code = "L718", icd.version = 10)

#calculate the prevalence of a diagnosis in the UKB study 
ukb_icd_prevalence(my_ukb_data, icd.code = "L721")

#test filters
test_filter<- filter(x_data_sub, diagnoses_main_icd10_f41202_0_0 == "L721")
test_filter_multiple<- filter(x_data_sub, diagnoses_main_icd10_f41202_0_0 %in% c("L721", "L905"))

#Selecionando os diagnóticos de uma coluna
test_filter_multiple<- filter(my_ukb_data, diagnoses_main_icd10_f41202_0_0 %in% c("L71", "L718", "L719"))

test_filter_multiple_keratopahty<- filter(my_ukb_data, diagnoses_main_icd10_f41202_0_0 %in% c("H181", "H18.1"))


load("ukb_rosacea.RData")
load("ukdata_031821.RData")

#adicionar o centro em uma nova coluna
y_data_sub<- ukb_centre(y_data_sub, centre.var = "^uk_biobank_assessment_centre.*0_0")


ukb_icd_prevalence(my_ukb_data, icd.version = 10, icd.diagnosis = "L721")

ukb_context(my_ukb_data, nonmiss.var = "my_variable_of_interest")


ukb_df_1000 <- sample_n(my_ukb_data,1000)
ukb_df_1000<- ukb_centre(ukb_df_1000, centre.var = "^uk_biobank_assessment_centre.*0_0")

# Selecionar diagnostico de rosacea em todas as colunas que começam com "diagnoses_main_icd10


## look 
test<- my_ukb_data %>% filter_at(vars(starts_with("diagnoses_main_icd10")), any_vars(.== "L71"))


test<- ukb_centre(test, centre.var = "^uk_biobank_assessment_centre.*0_0")

#L721

# Select main and second diagnosis

test_1<- my_ukb_data %>% filter_at(vars(starts_with(c("diagnoses_main_icd10", "diagnoses_secondary_icd10" ))), any_vars(.== "H03"))

p <-my_ukb_data %>%
  select(starts_with(c("diagnoses_main_icd10", "diagnoses_secondary_icd10"))) %>%
  summarise(my_ukb_data,
    count = n(),na.rm = TRUE)

p <-my_ukb_data %>%
  select(starts_with(c("diagnoses_main_icd10", "diagnoses_secondary_icd10"))) %>%                                          
                                              
                                        

#Gráfico com proporção homens e mulheres no dataset
test %>% 
  #summarise( Subjects = n() , Male =  sum(sex_f31_0_0 == 1)/Subjects , 
   #            Female = sum(sex_f31_0_0 == 0)/Subjects) %>%
    #select(-Subjects) %>% 
    #gather(sex, percentage) %>% 
  ggplot(aes(x = sex_f31_0_0,fill=sex_f31_0_0 )) +
  geom_bar(aes(y = (..count..)/sum(..count..)),  width=0.4) +
  geom_text(aes(y = ((..count..)/sum(..count..)), label = scales::percent((..count..)/sum(..count..))),
            stat = "count", vjust = -0.25)  +
  coord_flip() +
  scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
  labs(x = "", y = "", color = "", fill = "",  title = "") +
theme(title = element_text(face = "bold"), panel.grid = element_blank(),
      panel.background = element_rect(color = NULL,
                                      fill = alpha("grey", 0.10)),
      legend.key = element_blank(), axis.ticks.x = element_blank())

col_vars <- 
  tribble(
    ~col_var             , ~col_name,
    "id"                 , "eid",
    "sex"                , "sex_f31_0_0",
    "age_at_assessment"  , "age_when_attended_assessment_centre_f21003_0_0",
    "age_at_death"       , "age_at_death_f40007_2_0",
    "ethnicity"          , "ethnic_background_f21000_0_0",
    "centre"             , "uk_biobank_assessment_centre_f54_0_0",
    "Diag_icd10"         , "diagnoses_main_icd10_f41202_0_0"
  )
#Buscar ICD codes
ukb_icd_keyword("rosacea", icd.version = 9)

identify_icd <- 
  tribble(
    ~phenotype      , ~icd9         , ~icd10,
    "Rosacea"       , "6953"    ,       "L71",
    "Other Rosacea" ,  "no"     ,       "L71.8",
  )

installed.packages("pacman") 
library("pacman")
pacman::p_load(tidyverse, ukbtools, tictoc, furrr, future)

uk_rosacea <- 
  inner_join( ukb_df_1000  %>% 
                dplyr::select( col_vars$col_name,
                               "diagnoses_main_icd10_f41202_0_0"),
              icd_cancer_diagnosis %>% filter(str_detect(icd_code , "^(C82[0-9])|
                                                          ^(C83[0-9])|^(C85[0-9])|
                                                          ^(200[0-9])|^(202[0-9])")),
              by = c( "eid" = "eid"  )) 




#### Combine the threee rosacea ICD10 codes -  
# Combine all diagnosis main and secondary 
# Count the frequency of rosacea compare top 50 of other diagnosis (cardiovascular disease???)
# comorbities associated 
# Filter by eid to be unique 
# Dry eye with rosacea or not (count the patients)

## Usando pivot longer - Nao funcionou para todas as colunas
test <- p %>% 
  pivot_longer(cols = diagnoses_main_icd10_f41202_0_0:diagnoses_main_icd10_f41202_0_65, 
               names_to = "disease", 
               values_to = "icd_code")

table(test$icd_code)

test_1 <- test %>%
  group_by(icd_code) %>%
  summarise(counts = n())
test_1


ggplot(test_1, aes(x = icd_code, y = counts)) +
  geom_bar(fill = "#0073C2FF", stat = "identity") +
  geom_text(aes(label = counts), vjust = -0.3) + 
  ylim(0,10000) +
  scale_y_log10 (breaks = trans_breaks("log10", function(x) 10^x)) +
  theme_pubclean()


selection = subset(test_1, icd_code == "D375" | icd_code == "D410")

#test_2 <- selection %>%
#  group_by(icd_code) %>%
#  summarise(counts = n())
#test_2

ggplot(selection, aes(x = icd_code, y = counts)) +
  geom_bar(fill = "#0073C2FF", stat = "identity") +
  geom_text(aes(label = counts), vjust = -0.3) + 
#  ylim(0,10000) +
  scale_y_log10 (breaks = trans_breaks("log10", function(x) 10^x)) +
  theme_pubclean()

#Novo teste
# Diagnostico primario
test <- p %>% 
  select(starts_with("diagnoses_main_icd10_f41202_0_")) %>%
  gather(key, value, na.rm = TRUE) %>%
  count(value)


selection = subset(test, value == "L71" | value == "L718" | value == "L719" |  value == "H16" |  value == "H11" |  value == "H04" |  value == "H18")

selection = subset(test, value == "L71" | value == "L718" | value == "L719")


ggplot(selection, aes(x = value, y = n)) +
  geom_bar(fill = "#0073C2FF", stat = "identity") +
  geom_text(aes(label = n), vjust = -0.3) + 
  xlab("ICD 10 code - Primary diagnosis") +
  ylab("Count Patient") +
#  ylim(0,10000) +
# scale_y_log10 (breaks = trans_breaks("log10", function(x) 10^x)) +
  theme_pubclean()

top100 <- test %>% top_n(100) 


# Diagnostico secundario
test_secondary <- p %>% 
  select(starts_with("diagnoses_secondary_icd10")) %>%
  gather(key, value, na.rm = TRUE) %>%
  count(value)

selection1 = subset(test_secondary, value == "L71" | value == "L718" | value == "L719" |  value == "H16" |  value == "H11" |  value == "H04" |  value == "H18" | value == "I848")

ggplot(selection1, aes(x = value, y = n)) +
  geom_bar(fill = "#0073C2FF", stat = "identity") +
  geom_text(aes(label = n), vjust = -0.3) + 
  xlab("ICD 10 code - Secondary diagnosis") +
  ylab("Count Patient") +
  #  ylim(0,10000) +
  # scale_y_log10 (breaks = trans_breaks("log10", function(x) 10^x)) +
  theme_pubclean()




# Script that runs the analysis of complete vs incomplete cases

# Create data frames of complete and incomplete cases, tabulate means
relevant_data <- stph %>%
  select(Subject, Age.at.randomisation..calc., Bilirubin.mg.dl, Albumin, Creatinine.mg.dl, WBC, protime,
         INR, Sodium, Bilirubin.day.7, Gender, D90_surv, HE)

completedata <- relevant_data[complete.cases(relevant_data),]
table1::table1(~Bilirubin.mg.dl + Creatinine.mg.dl + Albumin + WBC + protime + 
                 INR + HE + Bilirubin.day.7 + Gender + D90_surv + Sodium + Age.at.randomisation..calc., data = completedata)

incompdata <- relevant_data[!complete.cases(relevant_data),]
table1::table1(~Bilirubin.mg.dl + Creatinine.mg.dl + Albumin + WBC + protime + 
                 INR + HE + Bilirubin.day.7 + Gender + D90_surv + Sodium + Age.at.randomisation..calc., data = incompdata)

# Calculate and save p-values for the comparison of complete and incomplete cases
p_values_means <- numeric(length = ncol(completedata))
for(i in 1:ncol(completedata)){
  stat <- t.test(completedata[,i], incompdata[,i], paired = F, alternative = "two.sided")
  p_values_means[i] <- stat$p.value
}

p_means_df <- data.frame(variable = colnames(completedata),
                         p = p_values_means)

test_gender <- prop.test(x = c(sum(completedata$Gender), sum(incompdata$Gender)), n = c(nrow(completedata), nrow(incompdata)), 
                         p = NULL, alternative = "two.sided", correct = TRUE)
test_gender$p.value

test_survival <- prop.test(x = c(sum(completedata$D90_surv), sum(incompdata$D90_surv)), n = c(nrow(completedata), nrow(incompdata)), 
                           p = NULL, alternative = "two.sided", correct = TRUE)
test_survival$p.value

# Chi-squared test for comparison of HE
HE_test <- matrix(c(table(completedata$HE), table(incompdata$HE)),
               ncol = 4,
               byrow = T)
colnames(HE_test) <- c("0","1","2", "3")
rownames(HE_test) <- c("Complete","Incomplete")
data <- as.table(HE_test)

chisq.test(HE_test)

# Investigate whether missingness relates to number of days in study by using the additional data
# on hospital discharges
library(lubridate)
# import additional data and add to stph dataframe
discharge_data <- readxl::read_xlsx("Discharge_data.xlsx")
stph_discharge <- merge(stph, discharge_data, by = "Subject")
stph_discharge <- stph_discharge[!is.na(stph_discharge$`Actual discharge date()(Discharge Visit)`),]

# Calculate the time between admission and discharge and time between admission and death
stph_discharge$Initial.Admission.date <- dmy(stph_discharge$Initial.Admission.date)
stph_discharge$discharge_date <- ymd(stph_discharge$`Actual discharge date()(Discharge Visit)`)
stph_discharge$death_date <- dmy(stph_discharge$Date_of_death)
stph_discharge$days_until_death <- difftime(stph_discharge$death_date, stph_discharge$Initial.Admission.date)
stph_discharge$days_in_study <- difftime(stph_discharge$discharge_date, stph_discharge$Initial.Admission.date, units = "days")

# Split into the complete and incomplete samples
stph_discharge_complete <- filter(stph_discharge, (Subject %in% completedata$Subject))
stph_discharge_incomplete <- filter(stph_discharge, (Subject %in% incompdata$Subject))

# Get summary statistics for the complet and incomplete samples
summary(as.numeric(stph_discharge_complete$days_in_study))
summary(as.numeric(stph_discharge_incomplete$days_in_study))

mean(as.numeric(stph_discharge_complete$days_in_study))
mean(as.numeric(stph_discharge_incomplete$days_in_study))

var(as.numeric(stph_discharge_complete$days_in_study))
var(as.numeric(stph_discharge_incomplete$days_in_study))

# Create histograms to show the distribution of days until discharge between groups
hist(as.numeric(stph_discharge_complete$days_in_study), main = " ", xlab="Days until discharge", freq=FALSE)
hist(as.numeric(stph_discharge_incomplete$days_in_study), main = " ", xlab="Days until discharge", freq=FALSE)

# Calculate p-value for comparing time to discharge between groups
discharge_test <- t.test(as.numeric(stph_discharge_complete$days_in_study), 
                         as.numeric(stph_discharge_incomplete$days_in_study), paired = F, alternative = "two.sided")
discharge_test$p.value

# Assess whether early discharge could be related to early death
length(which(as.numeric(stph_discharge_incomplete$days_until_death) <= 7))
length(which(as.numeric(stph_discharge_complete$days_until_death) <= 7))

# Count how many individuals are released before the day-7 measure
length(which(as.numeric(stph_discharge_incomplete$days_in_study) <= 7))
length(which(as.numeric(stph_discharge_complete$days_in_study) <= 7))



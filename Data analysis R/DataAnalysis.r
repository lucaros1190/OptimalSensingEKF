
# Comparing the differences between policies
#
# Created by Luca Rossini on 17 December 2025
# Last update 21 December 2025
#
# E-mail: luca.rossini@ulb.be


# Acquisition of the data - File 'AvgTrace_overMC.xlsx'

library(readxl)

filename <- file.choose()

data_ML = read_excel(filename, sheet = 'AvgTrace_overMC', col_names = T)

head(data_ML)

policy <- as.factor(data_ML$Policy)
avgTrace <- as.numeric(data_ML$Avg_Trace)

# Check the levels

levels(policy)


# Check the distribution of the dataset

hist(avgTrace)
shapiro.test(avgTrace)


# Analysis 1: Average TRACE


# GLM: 'avgTrace' variable, negative binomial as distribution.

library(MASS)

GenLin_Trace <- glm.nb(avgTrace ~ policy , data=data_ML)

summary(GenLin_Trace)


# GLM: 'avgTrace' variable, Poisson as distribution.

GenLinPois_Trace <- glm(avgTrace ~ policy, family = poisson, data=data_ML)

summary(GenLinPois_Trace)

# Check if there is overdispersion (if yes, use negative binomial)

dispersion_Trace <- sum(residuals(GenLinPois_Trace, 
                        type = 'pearson')^2) / df.residual(GenLinPois_Trace)

dispersion_Trace


# Check the distribution of the residuals

qqnorm(residuals(GenLin_Trace))
qqline(residuals(GenLin_Trace))

shapiro.test(residuals(GenLin_Trace))

# Pairwise comparison - Policy:

library(multcompView)
library(emmeans)

marginal_Trace = emmeans(GenLin_Trace, ~ policy)
pairs(marginal_Trace, adjust="bonferroni")

# Letters of significance:

library(multcomp)

lettere_Trace <- cld(marginal_Trace, alpha=0.05, Letters=letters, 
                     adjust="bonferroni")
lettere_Trace


# Plots only Times

library(ggplot2)

letters_Time <- data.frame(
  group = c('Chess', 
            'Diagonal',
            'Random',
            'Smart'),
  value = c(7.3e+42, 7.3e+42, 5.0e+42, 4.0e+42),
  letters = c('c', 'c', 'b', 'a')
)

boxPlot_Time <- ggplot(data_ML, aes(x=policy, y=Avg_Trace,
                                              fill=policy)) + 
                       geom_boxplot(width=0.5) + 
                       xlab('Policy') + 
                       ylab('Average trace') + 
                       ggtitle('Performance of the different sensing policies') +
                       ylim(0,1e43) +
                       theme_bw() + 
                       theme(plot.title = element_text(hjust=0.5), 
                       text = element_text(size=21), legend.position = 'none',
                       panel.grid = element_blank(),
                       panel.border = element_blank(),
                       axis.line = element_line(color = 'black')) +
                       geom_text(data = letters_Time, 
                       aes(x = group, y = value, label = letters),
                       vjust = -0.5, size = 8, color = 'black',
                       inherit.aes = F) +
                       scale_x_discrete(guide = guide_axis(angle = 45),
                                       labels = c('Chessboard', 
                                                  'Diagonal',
                                                  'Random',
                                                  'Smart'))


boxPlot_Time



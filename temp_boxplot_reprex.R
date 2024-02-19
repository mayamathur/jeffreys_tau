


# default Y-axis: no truncation
library(ggplot2)
ggplot( data = dp,
        
        aes(x = as.factor(k.pub),
            y = ShatMAE,
            color = method.pretty,
            fill = method.pretty) ) +
  
  geom_hline(yintercept = 0,
             lty = 2) +
  
  geom_boxplot(width=0.6,
               alpha = .5,
               outlier.shape = NA) +
  
  labs(color  = "Measure of central tendency", fill = "Measure of central tendency")


# custom Y-axis: causes truncation of the first (red) box's upper whisker
ggplot( data = dp,
        
        aes(x = as.factor(k.pub),
            y = ShatMAE,
            color = method.pretty,
            fill = method.pretty) ) +
  
  scale_y_continuous(limits = c(0, 1.5)) +
  
  geom_hline(yintercept = 0,
             lty = 2) +
  
  geom_boxplot(width=0.6,
               alpha = .5,
               outlier.shape = NA) +
  labs(color  = "Measure of central tendency", fill = "Measure of central tendency")

# for "tidy" summarize: http://stackoverflow.com/questions/27879638/use-dplyrs-summarise-each-to-return-one-row-per-function
library(tidyverse)

traits <- athaliana_traits_strong()

phen <- athaliana_phen(traits = traits)
nobs <- nrow(phen)

num_na <- phen %>% summarize_all(. %>% is.na %>% sum) %>% as.integer

t <- tibble(variable = names(phen), num_na = num_na) %>%
  arrange(desc(num_na)) %>%
  mutate( 
  variable = factor(variable, levels = unique(variable)),
  num_obs = nobs - num_na)

p <- ggplot(t, aes(variable, num_obs)) + 
  geom_point() + 
  ylim(c(0, nobs)) + 
  geom_hline(yintercept = c(nobs, nobs/2), linetype = c(3, 3)) + 
  coord_flip() + 
  theme_linedraw()
  

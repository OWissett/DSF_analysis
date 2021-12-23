# main.r
#
# MIT License
# 
# Copyright (c) 2021 Oliver Wissett
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#   
#   The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

library(ggplot2)  # Provides graphing functions
library(dplyr)    # Adds piping and data manipulation functions
library(readr)    # Allows for easy handling of CSV inputs
library(reshape2) # Adds the melt function to do some data wrangling.
library(stringr)  # Adds some nice regex stuff
library(ggthemes) # Adds nice looking ggplot themes
library(tidyr)
library(viridis)

# Load dataset
DSF_DATAPATH = '../../datasets/Tripeptides_DSF.csv'

dsf_df <- readr::read_csv(DSF_DATAPATH)

assign_peptide <- function(rowID) {
  if (rowID == 'A')
    return('None')
  
  if (rowID == 'B')
    return('KAK')
  
  if (rowID == 'C')
    return('KFK')
  
  if (rowID == 'D')
    return('KTK')
  
  if (rowID == 'E')
    return('KGK')
}

assign_condition <- function(colID) {
  if(colID < 4) {
    return('HisZ')
  } else {
    return('Control')
  }
}

calculate_gradient <- function(x, y) {
  return(diff(y)/diff(x))
}

dsf_tibble <- dsf_df %>%
  select(c(1, seq(2, ncol(.), by = 2))) %>%
  tibble() %>%
  rename(temperature = 1) %>%
  slice(-1) %>%
  melt(id.vars = 'temperature',
       value.name = 'fluorescence',
       variable.name = 'cellID') %>%
  mutate(rowID = substr(cellID,1,1),
         colID = as.numeric(str_extract(cellID, '(?<=^[A-Z])[0-9]*'))) %>%
  mutate(peptide = sapply(rowID, assign_peptide),
         condition = sapply(colID, assign_condition),
         fluorescence = as.numeric(fluorescence),
         temperature = as.numeric(temperature)) %>%
  mutate(peptide = factor(peptide, levels = c('None', 'KAK', 'KTK','KGK', 'KFK')),
         condition = factor(condition, levels = c('HisZ', 'Control'))) %>%
  dplyr::filter(colID < 7 & rowID != 'F' & rowID != 'G' & rowID != 'H')

dsf_norm <- dsf_tibble %>% 
  drop_na() %>%
  group_by(peptide, condition, cellID) %>%
  summarise(count = n(),
            fluorescence.norm = (fluorescence - min(fluorescence))/ (max(fluorescence) - min(fluorescence)),
            temperature = temperature)

dsf_summary <- dsf_norm %>%
  group_by(peptide, condition, temperature) %>%
  summarise(count = n(),
            fluorescence.mean = mean(fluorescence.norm),
            fluorescence.se = sd(fluorescence.norm) / sqrt(n()),
            temperature = temperature) %>%
  mutate(fluorescence.lower = fluorescence.mean - fluorescence.se,
         fluorescence.upper = fluorescence.mean + fluorescence.se)
  
# Plotting time...

ggplot(dsf_norm %>% filter(condition == 'HisZ'),
       aes(temperature, fluorescence.norm,
                       colour = peptide,
                       shape = condition)) +
  geom_point() +
  theme_classic()

dsf_traces <- ggplot(dsf_summary %>% filter(condition == 'HisZ'),
                     aes(temperature, fluorescence.mean,
                        colour = peptide)) +
  geom_point(size = 2) +
  theme_hc() +
  labs(y = 'Relative Fluorescence',
       x = 'Temperature (°C)',
       colour = 'Peptide') +
  theme(
    axis.title = element_text(face="bold", size = 18),
    axis.text = element_text(face="bold", size = 12),
    panel.background = element_rect(fill = 'white'),
    axis.line = element_line(colour = 'black')
  ) +
  scale_color_brewer(palette = 'Set1')

dsf_traces  

ggplot(dsf_summary %>% filter(condition == 'HisZ'), 
       aes(temperature, gradient, color = peptide)) +
  geom_point()

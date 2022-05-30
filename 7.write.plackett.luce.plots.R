######################
# This scripts produces plots out of the PlackettLuce modelling
# Authors: Kaue de Sousa
# Date: May 30th, 2022
########################

# Coerce farmers evaluation to rankings
library("PlackettLuce")
library("janitor")
library("tidyverse")
library("gosset")

load("output/PL_models.rda")

ggsave(filename = "output/pltree_location.png", 
       plot = plot(pl1), width = 30, height = 25, units = "cm", dpi = 500)

ggsave(filename = "output/pltree_gender.png", 
       plot = plot(pl2), width = 25, height = 25, units = "cm", dpi = 500)

ggsave(filename = "output/worth_map_location.png", 
       plot = worth_map(pl1) + theme_bw(), width = 25, height = 18, units = "cm", dpi = 500)

ggsave(filename = "output/worth_map_gender.png", 
       plot = worth_map(pl2) + theme_bw(), width = 25, height = 18, units = "cm", dpi = 500)

library(tidyverse)

# read stratigraphy data from C2VSim preprocessor
col_names <- c("node_id","elev","ac1","aq1","ac2","aq2","ac3","aq3","x")
d <- read_lines("F:/Box Sync/Research/Pre QE Research/C2VSim/C2VSim_CG_1921IC_R374_Model/C2VSim_CG_1921IC_R374/Preprocessor/CVstrat.dat", skip = 92)
d <- str_split_fixed(d, pattern = "[ ]{1,}", n=10)
d <- as_tibble(d) %>% select(-V1) 
colnames(d) <- col_names

# import TB element ids
tb <- read_csv("F:/Box Sync/Research/Pre QE Research/C2VSim/TB_elems.csv")
colnames(tb) <- c("index","id")

# filter for tb elems
df <- filter(d, node_id %in% tb$id) %>% select(-x)
df <- mutate_all(df, as.numeric) 
pdf <- mutate(df, 
              paq1 = aq1 / (ac1+aq1),
              paq2 = aq2 / (ac2+aq2),
              paq3 = aq3 / (ac3+aq3)) %>% 
  summarise_at(.vars = c("paq1","paq2","paq3"), .funs = mean)

bdf <- mutate(df, 
              b1 = ac1+aq1,
              b2 = ac2+aq2,
              b3 = ac3+aq3) %>% 
  summarise_at(.vars = c("b1","b2","b3"), .funs = sum) 

# regional average aquifer proportion
sum(pdf * (bdf/(sum(bdf))))





# regional average porosity



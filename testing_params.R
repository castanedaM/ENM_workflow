
# Define a function that returns the sum of two variables
sum_function <- function(x) {
  return(x[1] + x[2])
}

scale_1_0 <- function(x){
  return((x - min(x)) /(max(x) - min(x)))
}


df <- model_results %>% 
  group_by(nu, gamma) %>% 
  filter(CBP <= 0.001) %>% 
  mutate(omission_test_sc = (omission_test),
         volume_sc =  (volume)) %>% 
  # filter(set %in% c(1))
  summarise(o_max = max(omission_test_sc),
            o_mean = mean(omission_test_sc),
            o_min = min(omission_test_sc),
            v_max = max(volume_sc),
            v_mean = mean(volume_sc),
            v_min = min(volume_sc))


o_sort <- df %>% 
  arrange(o_mean) 
o_sort$o_ID <- 1:nrow(o_sort)

v_sort <- df %>% 
  arrange(v_mean) 
v_sort$v_ID <- 1:nrow(v_sort)

head(o_sort)
head(v_sort)

df_sort <- left_join(o_sort, v_sort)

# t <- ggplot(df, 
#             aes(x = omission_test_sc, y = volume_sc, color = nu, shape = gamma)) +
#   geom_point() +
#   ylab("Volume") + xlab("Omission Rate")
# t


t <- ggplot(df, 
         aes(x = o_mean, y = v_mean, color = nu, shape = gamma)) +
  geom_segment(aes(x = o_min, y = v_mean, xend = o_max, yend = v_mean)) +
  geom_segment(aes(x = o_mean, y = v_min, xend = o_mean, yend = v_max)) +
  geom_point(size = 3) +
  ylab("Volume") + xlab("Omission Rate") +
  scale_x_continuous(breaks = seq(0.05, 0.4, 0.05))

t

df_sort <- df_sort %>% 
  mutate(min_sum = v_ID + o_ID) %>% 
  arrange(min_sum)

op <- df %>% 
  filter(nu == df_sort$nu[1], 
         gamma == df_sort$gamma[1]) 

df_sort <- df_sort %>% 
  group_by(nu, gamma) %>% 
  mutate(mean_sum = mean(min_sum)) %>% 
  arrange(mean_sum)

op_mean <-  df %>% 
  filter(nu == df_sort$nu[1], 
         gamma == df_sort$gamma[1])  

t +
  geom_point(data = op, color = "red",
             pch = 21, size = 6, stroke = 1.5,
             show.legend = FALSE) +
  geom_point(data = op_mean, color = "blue",
             pch = 21, size = 8, stroke = 1.5,
             show.legend = FALSE)




# Calculate Euclidean distance for each combination
df$dis_or <- sqrt((df$omission_test_sc)^2 + (df$volume_sc)^2)
df$avg <- (df$omission_test_sc + df$volume_sc)/2

# plot(df$omission_test_sc*df$volume_sc, 
#      col = factor(df$nu))


# Find the combination with minimum distance
min_dis_combo <- df[which.min(df$dis_or), ]
mean_dis_combo <- df[which.min(df$min), ]


og_dis <- df %>% 
  filter(nu == min_dis_combo$nu, 
         gamma == min_dis_combo$gamma) 

mean_dis <- df %>% 
  filter(nu == mean_dis_combo$nu, 
         gamma == mean_dis_combo$gamma) 

a <- t +
  geom_point(data = og_dis, color = "red",
               pch = 21, size = 6, stroke = 1.5,
               show.legend = FALSE) +
  geom_point(data = mean_dis, color = "blue",
             pch = 21, size = 6, stroke = 1.5,
             show.legend = FALSE)
a


# b <-  t +
#   geom_point(data = og_dis, color = "red",
#              pch = 21, size = 6, stroke = 1.5,
#              show.legend = FALSE) +
#   geom_point(data = mean_dis, color = "blue",
#              pch = 21, size = 6, stroke = 1.5,
#              show.legend = FALSE)
# 
# b

# c <-  t +
#   geom_point(data = og_dis, color = "red",
#              pch = 21, size = 6, stroke = 1.5,
#              show.legend = FALSE) +
#   geom_point(data = mean_dis, color = "blue",
#              pch = 21, size = 6, stroke = 1.5,
#              show.legend = FALSE)
# c


a <- a + ggtitle("not scaled")
b <- b + ggtitle("z scaled")
# c <- c + ggtitle("0 to 1 scaled")

ggpubr::ggarrange(a, b, nrow = 2, common.legend = TRUE, legend = "right")

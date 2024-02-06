Default_Dataset <- read_csv("C:/Users/Ross.Whippo/Desktop/Default Dataset.csv")


Default_Dataset %>%
  ggplot(aes(x = stipe_length_cm, y = `growth_rate_mm-day`)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ splines::ns(x, 2)) +
  theme_bw()





lengths <- c(seq(0.1, 300, by = 0.1)) %>%
  as.data.frame()
  colnames(lengths) <- c("stipe_length_cm")
new_DD <- Default_Dataset %>%
  mutate(across('stipe_length_cm', \(x) round(x, digits = 1))) %>%
  right_join(lengths) %>%
  arrange(stipe_length_cm)
 

mod1 <- Default_Dataset %>%
  lm(formula = `growth_rate_mm-day` ~ splines::ns(stipe_length_cm, 2))

newcol <- as.data.frame(seq(0.1, 300, by = 0.1))
colnames(newcol) <- "stipe_length_cm"
new_predicts <- as.data.frame(predict(mod1, newdata = data.frame(stipe_length_cm = seq(0.1, 300, by = 0.1))))
colnames(new_predicts) <- c("predictions") 

all_predicts <- new_DD %>%
  mutate(predictions = new_predicts$predictions)




all_predicts %>%
  ggplot(aes(x = stipe_length_cm, y = predictions)) +
  geom_smooth(method = "lm", formula = y ~ splines::ns(x, 2)) +
  geom_point(aes(x = stipe_length_cm, y = `growth_rate_mm-day`)) +
  ylab("stipe_growth_mm_day") +
 # geom_point() +
  theme_bw() 


















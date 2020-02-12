#############################
# Report
#############################
table_1 <- matrix(NA, nrow = 4, ncol = 5)
table_1[1,] <- c(waic_chain1$waic, waic_chain2$waic, waic_chain3$waic, waic_chain4$waic, waic_chain5$waic)
table_1[2,] <- c(ML_chain1$LL, ML_chain2$LL, ML_chain3$LL, ML_chain4$LL, ML_chain5$LL)

table_1[3,] <- c(waic_chain6$waic, waic_chain7$waic, waic_chain8$waic, waic_chain9$waic, waic_chain10$waic)
table_1[4,] <- c(ML_chain6$LL, ML_chain7$LL, ML_chain8$LL, ML_chain9$LL, ML_chain10$LL)
xtable::xtable(round(table_1), digits = 0)

#############################
# Report
#############################
setwd("/home/hoanguc3m/Dropbox/WP5/")
pdf(file='img/SV1.pdf', width = 9, height = 9)
colnames(datafin) <- c("IndPro", "RPI", "OilPrice","OilProduct")
head(datafin)
i <- 1
Time <- seq(as.Date("1974/2/1"), as.Date("2016/12/1"), "months")
par(mar=c(2,5,3,1))
layout(matrix(c(1,2,3,4,5,6,7,8), 4, 2, byrow = F))
i <- 1
plot(Time, datafin[,i], type = "l", ylab = colnames(datafin)[i], xlab = "")
i <- 2
plot(Time, datafin[,i], type = "l", ylab = colnames(datafin)[i], xlab = "")
i <- 3
plot(Time, datafin[,i], type = "l", ylab = colnames(datafin)[i], xlab = "")
i <- 4
plot(Time, datafin[,i], type = "l", ylab = colnames(datafin)[i], xlab = "")

H_mat <- get_post(Chain10, element = "h")
plot_shade <- function(H_mat, i){
  max_ht <- ncol(H_mat)
  h_select <- seq(from = i, to = max_ht, by = K)
  h_mean <- apply(H_mat[, h_select], MARGIN = 2, FUN = mean)
  h_q25 <-  apply(H_mat[, h_select], MARGIN = 2, FUN = quantile, probs = c(0.25))
  h_q75 <-  apply(H_mat[, h_select], MARGIN = 2, FUN = quantile, probs = c(0.75))

  x <- tail(Time, length(h_mean))
  plot(x, h_mean, type = "l",
       ylab = colnames(datafin)[i], xlab = "", ylim = c(min(h_q25),max(h_q75)))
  xx = c(x, rev(x))
  yy = c(h_q25, rev(h_q75))
  polygon(xx, yy, border=8, col=gray(.6, alpha=.3) )
}
plot_shade(H_mat, 1)
plot_shade(H_mat, 2)
plot_shade(H_mat, 3)
plot_shade(H_mat, 4)
dev.off()

#############################
# Report
#############################
library(ggthemr)
ggthemr('light')

pdf(file='img/postNu.pdf', width = 9, height = 9)
par(mar=c(2,5,3,1))
layout(matrix(c(1,2,3,4,5,6,7,8), 4, 2, byrow = F))
Gamma_matC10 <- get_post(Chain10, element = "gamma")
Nu_matC10 <- get_post(Chain10, element = "nu")
Gamma_matC5 <- get_post(Chain5, element = "gamma")
Nu_matC5 <- get_post(Chain5, element = "nu")

i <- 1
gg_nu_mat <- function(i){
  data_nu <- data.frame(nu = c(Nu_matC5[,i], Nu_matC10[,i]), models = c(rep("NonSV", ndraws), rep("SV", ndraws)))
  ggplot(data_nu, aes(x=nu,..density.., fill = models, colour = models)) +
    geom_histogram(position = "identity", alpha = 0.8, bins = 30) +
    geom_density(position = "identity", alpha = 0.5) + xlab("nu") + ylab(colnames(datafin)[i])
}
gg_gamma_mat <- function(i){
  data_gamma <- data.frame(gamma = c(Gamma_matC5[,i], Gamma_matC10[,i]), models = c(rep("NonSV", ndraws), rep("SV", ndraws)))
  ggplot(data_gamma, aes(x=gamma,..density.., fill = models, colour = models)) +
    geom_histogram(position = "identity", alpha = 0.8, bins = 30) +
    geom_density(position = "identity", alpha = 0.5) + xlab("gamma") + ylab(colnames(datafin)[i])
}

p1 <- gg_nu_mat(1)
p2 <- gg_nu_mat(2)
p3 <- gg_nu_mat(3)
p4 <- gg_nu_mat(4)

p5 <- gg_gamma_mat(1)
p6 <- gg_gamma_mat(2)
p7 <- gg_gamma_mat(3)
p8 <- gg_gamma_mat(4)

grid.arrange(p1, p5, p2, p6, p3, p7, p4, p8, nrow = 4, ncol = 2)
dev.off()

# hist(Nu_matC10[,i],ylab = colnames(datafin)[i], xlab = "", main = "", prob = TRUE, ylim = c(0,0.1), xlim = c(0,40) )
# lines(density(Nu_matC5[,i]))
# ndraws <- nrow(Nu_matC10)

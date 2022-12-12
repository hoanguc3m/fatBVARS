# NOT RUN {
######################################################################
# Non-parsed command-line arguments
######################################################################

# args <- commandArgs(trailingOnly = TRUE)
# print(args)

library("fatBVARS")
data(MonthlyData)
View(datafin)
data_name = "MonthlyData"
Time <- as.Date(Monthly$sasdate[133:732], format = "%m/%d/%Y")


# dist = args[1]
# SV = as.logical(args[2])
# t_start = as.numeric(args[3])
# path = args[4]
# if (is.na(path)) path = "/proj/nobackup/hoanguc3m/recursive/"

dist = "Gaussian"
SV = FALSE
t_start = 500
path = ""

K = ncol(datafin)
p = 4

cat(paste(path, data_name, "_", dist, "_", SV, "_", t_start+1, ".RData", sep = ""))

# RhpcBLASctl::blas_set_num_threads(1)

out_recursive <- recursive_seperate(y = datafin, t_start = t_start,
                                    t_pred = 12,
                                    K = K,
                                    p = p,
                                    dist = dist,
                                    SV = SV,
                                    outname = paste(path, data_name, "_", dist, "_", SV, "_", t_start+1, ".RData", sep = ""),
                                    samples = 60000, burnin = 10000, thin = 5)

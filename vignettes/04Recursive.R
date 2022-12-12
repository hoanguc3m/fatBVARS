# NOT RUN {
######################################################################
# Non-parsed command-line arguments
######################################################################
# install.packages("~/Downloads/fatBVARS_0.2.0.tar.gz", repos = NULL, type = "source",
#                  lib ="/pfs/stor10/users/home/h/hoangoru/Public/R-packages-4.0")

args <- commandArgs(trailingOnly = TRUE)
print(args)

library("fatBVARS")


#data(package = "fatBVARS", "MonthlyData")
#datafin = MonthlyData
#data_name = "MonthlyData"

#length(ts(start = c(2007,1), end = c(2019,12), frequency = 12))
#length(445:588)
load("/home/hoanguc3m/Dropbox/WP5/Code/RData/MonthlyData/MonthlyData.RData")
data_name = "MonthlyData"
Time <- as.Date(Monthly$sasdate[133:732], format = "%m/%d/%Y")


dist = args[1]
SV = as.logical(args[2])
t_start = as.numeric(args[3])
path = args[4]
if (is.na(path)) path = "/proj/nobackup/hoanguc3m/recursive/"

K = ncol(datafin)
p = 4

cat(paste(path, data_name, "_", dist, "_", SV, "_", t_start+1, ".RData", sep = ""))

RhpcBLASctl::blas_set_num_threads(1)

out_recursive <- recursive_seperate(y = datafin, t_start = t_start,
                                    t_pred = 12,
                                    K = K,
                                    p = p,
                                    dist = dist,
                                    SV = SV,
                                    outname = paste(path, data_name, "_", dist, "_", SV, "_", t_start+1, ".RData", sep = ""),
                                    samples = 60000, burnin = 10000, thin = 5)

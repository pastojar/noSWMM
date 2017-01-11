# this script creates an optimized LHS and saves it to the specified file
# J. Pastorek, 3.12.2015

# SET WD TO THIS FILE'S FOLDER!

start_time=proc.time()

dimension <- 8
n <- 1000                # We keep also in mind the well-known and empirical relation N~10d

X <- DiceDesign::lhsDesign(n, dimension, seed=1)$design
plot(X[,])

Xopt <- DiceDesign::maximinSA_LHS(X,T0=10,c=0.99,it=2000)
plot(Xopt$design[,])

write.table(Xopt$design, paste(getwd(), "/LHSopt.csv", sep=""), sep=";", col.names=T, row.names=F, quote=F)

end_time=proc.time()
time_taken=end_time-start_time
time_taken




plot(Xopt$critValues,type="l")
plot(Xopt$tempValues,type="l")

## Not run:
Xopt <- maximinSA_LHS(X,T0=10,c=0.99,it=1000,profile="GEOM_MORRIS")
## End(Not run)
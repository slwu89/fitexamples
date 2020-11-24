setwd("")

# Install SSM for Linux or Mac OSX see: https://github.com/JDureau/ssm

# The model is in the file ssm.json
# /data the data files 
# /prior the definition of the parameter priors
# theta.json a random initial parameter values 
# see https://github.com/JDureau/ssm/README.md


dir_model = "./ssm_SIRS"
dir_kmcmc = "../kmcmc"
dir_pmcmc = "../pmcmc"
#dir.create(paste(dir_model, "/kmcmc",sep=""))
#dir.create(paste(dir_model, "/pmcmc",sep=""))


# Compiling the model and create the /bin
cmd <- sprintf("cd %s ; ssm", dir_model)
system(cmd)


# Running different kmcmc chaines
# see the help: 
#cmd <- sprintf("cd %s/bin; ./kmcmc --help", dir_model); system(cmd)

cmd <- sprintf("cd %s/bin; cat ../theta_ksimplex.json | ./kmcmc sde -M 50000 --traj --trace -I 1 --root %s > ../theta_kmcmc.json; cd ..",dir_model,dir_kmcmc)
system(cmd)
#cmd <- sprintf("cd %s/bin; cat ../theta_kmcmc.json |./kmcmc sde -M 500000 --traj --trace -I 2 --root %s > ../theta_kmcmc2.json; cd ..",dir_model,dir_kmcmc)
#system(cmd)
#cmd <- sprintf("cd %s/bin; cat ../theta_kmcmc2.json |./kmcmc sde -M 1000000 | ./kmcmc -M 500000--traj --trace -I 3 --root %s > ../theta_kmcmc3.json; cd ..",dir_model,dir_kmcmc)
#system(cmd)


# Running different pmcmc chaines
# Be careful, this part is computationaly intensive, and may require the use of a cluster for the last command
# see the help: 
#cmd <- sprintf("cd %s/bin; ./pmcmc --help", dir_model); system(cmd)

#cmd <- sprintf("cd %s/bin; cat ../theta_kmcmc3.json | ./pmcmc psr -J 1000 -M 10000  -N 4 --trace --traj -I 1 --root %s > ../theta_pmcmc1.json; cd ..",dir_model,dir_pmcmc)
#system(cmd)
#cmd <- sprintf("cd %s/bin; cat ../theta_pmcmc1.json | ./pmcmc psr -J 5000 -M 50000  -N 4 --trace --traj -I 2 --root %s > ../theta_pmcmc2.json; cd ..",dir_model,dir_pmcmc)
#system(cmd)

# another  possibility would be to used previous files with parameter values that have already converged in kmcmc or in pmcmc
#cmd <- sprintf("cd %s/bin; cat ../theta_kmcmc_0.json | ./pmcmc psr -J 5000 -M 50000  -N 4 --trace --traj -I 2 --root %s > ../theta_pmcmc2.json; cd ..",dir_model,dir_pmcmc)
#system(cmd)
cat ../theta_kmcmc_0.json | ./pmcmc psr -J 5000 -M  100000 -N 4 --trace --traj --hat -I 2 --root ./pmcmc > ../theta_pmcmc2.json
# or
#cmd <- sprintf("cd %s/bin; cat ../theta_pmcmc_0.json | ./pmcmc psr -J 5000 -M 50000  -N 4 --trace --traj -I 2 --root %s > ../theta_pmcmc2.json; cd ..",dir_model,dir_pmcmc)
#system(cmd)






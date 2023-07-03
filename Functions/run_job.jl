###################################################################################
## Julia script to write the runfile, execute array job, collect and recombine data
###################################################################################


#############################
# work in temporary directory
#############################

# check if previous directories exist, if not delete
#Base.run(`if "[" -d "Data_temporary" "]"";" then rm -r Data_temporary";" fi`)
#Base.run(`if "[" -d "outfiles" "]"";" then rm -r Data_temporary";" fi`)
touch("remove_directories.sh")
remove_file = open("remove_directories.sh", "w")
write(remove_file, "#!/bin/bash \n")
write(remove_file, "if [ -d 'Data_temporary' ]; then rm -r Data_temporary; fi \n")
write(remove_file, "if [ -d 'outfiles' ]; then rm -r outfiles; fi \n")
close(remove_file)
Base.run(`bash remove_directories.sh`)
rm("remove_directories.sh")

# create needed directories
mkdir("Data_temporary")
mkdir("outfiles")


################
# write job file
################

touch("array_job.jl")
job_file = open("array_job.jl", "w")
write(job_file, "using ArgParse \n")
write(job_file, "using DelimitedFiles \n")

write(job_file, "function parse_commandline() \n")
write(job_file, "s = ArgParseSettings() \n")
write(job_file, "@add_arg_table s begin \n")
write(job_file, "'--index' \n")
write(job_file, "arg_type = Int \n")
write(job_file, "default = 1 \n")
write(job_file, "end \n")
write(job_file, "return parse_args(s) \n")
write(job_file, "end \n")

write(job_file, "parsed_args = parse_commandline() \n")
write(job_file, "ind = parsed_args['index'] \n")

write(job_file, "function main(ind) \n")
write(job_file, "param1 = [1, 2, 3, 4, 5] \n")
write(job_file, "param2 = [6, 7, 8] \n")
write(job_file, "param3 = [9, 10, 11, 12] \n")
write(job_file, "param_space = collect(Iterators.product(param1, param2, param3)) \n")
write(job_file, "params = param_space[ind] \n")
write(job_file, "p1 = params[1] \n")
write(job_file, "p2 = params[2] \n")
write(job_file, "p3 = params[3] \n")
write(job_file, "val1 = p1 + p2 + p3 \n")
write(job_file, "val2 = p1*p2*p3 \n")
write(job_file, "data = [p1 p2 p3 val1 val2] \n")
write(job_file, "open('Data_temporary/data_$(ind).txt', 'w') do io \n")
write(job_file, "writedlm(io, data) \n")
write(job_file, "end \n")
write(job_file, "end \n")

write(job_file, "main(ind) \n")

close(job_file)


##############################
# write run file for array job
##############################

touch("run_array_job.sh")
run_file = open("run_array_job.jl", "w")
write(run_file, "#!/bin/bash \n")
write(run_file, "#SBATCH --time=00:15:00 \n")
write(run_file, "#SBATCH --mem=100M \n")
write(run_file, "#SBATCH --output=outfiles/array_job_%j.out \n")
write(run_file, "#SBATCH --array=1-60 \n")

write(run_file, "srun julia array_job.jl --index $SLURM_ARRAY_TASK_ID \n")

close(run_file)


##########
# run jobs
##########

Base.run(`sbatch run_array_job.sh`)


###############################
# wait until jobs have finished
###############################

# get number of jobs (pending, running) in queue
#Base.run(`num_jobs=$(squeue -u niederm1 -h -t pending,running -r | wc -l)`) 

# if there are still jobs in queue, wait until finished
touch("check_queue.sh")
check_file = open("check_queue.sh", "w")
write(check_file, "#!/bin/bash \n")
write(check_file, "num_jobs=$(squeue -u niederm1 -h -t pending,running -r | wc -l) \n")
write(check_file, "while [ $num_jobs -gt 0 ] \n")
write(check_file, "do \n")
write(check_file, "\t sleep 2 \n")
write(check_file, "done \n")
close(check_file)

# execute
Base.run(`bash check_queue.sh`)

##################################
# check for failed jobs and re-run
##################################

cd("Data_temporary")

# write file to identify the missing data
touch("find_missing_data.sh")
missing_data = open("find_missing_data.sh", "w")
write(missing_data, "#!/bin/bash \n")
write(missing_data, "for num in {1..60} \n")
write(missing_data, "do \n")
write(missing_data, "\t if [[ ! -f data_${num}.txt ]] ; then \n")
write(missing_data, "\t \t echo "${num}" \n")
write(missing_data, "\t fi \n")
write(missing_data, "done > missing_file_indices.txt \n")
close(missing_data)

# identify missing data
Base.run(`bash find_missing_data.sh`)

# feed indices back into runscript and re-run 
#
#
#

# check if data now complete, otherwise repeat until complete
#
#
#


################
# recombine data
################

# write file to regroup data
touch("regroup_data.sh")
regroup = open("regroup_data.sh", "w")
write(regroup, "#!/bin/bash \n")
write(regroup, "for num in {1..60} \n")
write(regroup, "do \n")
write(regroup, "\t cat data_${num}.txt >> final_data.txt \n")
write(regroup, "done \n")
close(regroup)

# execute and leave data in main directory
Base.run(`bash regroup_data.sh`)
mv("final_data.txt ..)



##############################
# delete temporary directories
##############################

cd(..)
#rm("Data_temporary", recursive=true)
#rm("outfiles", recursive=true)

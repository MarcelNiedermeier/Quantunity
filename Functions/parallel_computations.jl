
###################################################################
## Functions to simplify parallel computations in Triton with Slurm
###################################################################

""" The functions below can be used to simplify parallel computations
on a cluster using SLURM. """


using ArgParse


################################
# Functions to create files etc.
################################

""" Function to remove a given file if it exists. """
function remove_file(filename::String)

    # build bash file
    touch("remove_file.sh")
    remove = open("remove_file.sh", "w")

    write(remove, "#!/bin/bash \n")
    write(remove, "if [ -f "*filename*" ]; then rm "*filename*"; fi \n")
    close(remove)

    run(`bash remove_file.sh`)
    rm("remove_file.sh")
end


""" Function to remove a given directory if it exists. """
function remove_directory(directoryname::String)

    # build bash file
    touch("remove_directory.sh")
    remove = open("remove_directory.sh", "w")

    write(remove, "#!/bin/bash \n")
    write(remove, "if [ -d "*directoryname*" ]; then rm -r "*directoryname*"; fi \n")
    close(remove)

    run(`bash remove_directory.sh`)
    rm("remove_directory.sh")
end


""" Function to create the temporary directories where the data
is saved. """
function create_data_directories()

    # remove previous directories if they (still) exist
    remove_directory("Data_temporary")
    remove_directory("outfiles")

    # create needed directories
    mkdir("Data_temporary")
    mkdir("outfiles")
end


""" Function which writes the job file in Julia for a function
taking as multidimensional parameter space as an input. """
function create_job_file(func::String, file_name::String)

    # remove old file
    remove_file("job.jl")

    touch("job.jl")
    job_file = open("job.jl", "w")

    # imports
    write(job_file, "using ArgParse \n")
    write(job_file, "using DelimitedFiles \n")
    write(job_file, "include(\""*file_name*"\") \n")

    # write command line parser
    write(job_file, "function parse_commandline() \n")
    write(job_file, "s = ArgParseSettings() \n")
    write(job_file, "@add_arg_table s begin \n")
    write(job_file, """"--index" \n""")
    write(job_file, "arg_type = Int \n")
    write(job_file, "default = 1 \n")
    write(job_file, "end \n")
    write(job_file, "return parse_args(s) \n")
    write(job_file, "end \n")

    # get index for parameter space
    write(job_file, "parsed_args = parse_commandline() \n")
    write(job_file, """ind = parsed_args["index"] \n""")

    # execute function and save result in temporary directory
    write(job_file, "result = "*func*"(ind) \n")
    write(job_file, "data = reshape([result...], 1, length(result)) \n")
    write(job_file, """open("Data_temporary/data_\$(ind).txt", "w") do io \n""")
    write(job_file, "writedlm(io, data) \n")
    write(job_file, "end \n")

    close(job_file)
end


""" Function to create the bash file that submits the array
job to the SLURM cluster. Takes the configuration of the array
job as input. """
function create_array_job_file(time::String, mem::String, num_jobs::String)

    # remove old file
    remove_file("run_array_job.sh")

    touch("run_array_job.sh")
    run_file = open("run_array_job.sh", "w")

    write(run_file, "#!/bin/bash \n")
    write(run_file, "#SBATCH --time="*time*" \n")
    write(run_file, "#SBATCH --mem="*mem*" \n")
    write(run_file, "#SBATCH --output=outfiles/array_job_%j.out \n")
    write(run_file, "#SBATCH --array=1-"*num_jobs*" \n")
    write(run_file, "#SBATCH --mail-user=marcel.niedermeier@aalto.fi \n")
    write(run_file, "#SBATCH --mail-type=ALL \n")
    write(run_file, "srun julia job.jl --index \$SLURM_ARRAY_TASK_ID \n")

    close(run_file)
end


""" Function to create the bash file that submits a single
job to the SLURM cluster. Takes the configuration of the single
job as input. Meant to be used to re-run failed jobs. """
function create_single_job_file(time::String, mem::String, index::String)

    # remove old file
    remove_file("run_single_job.sh")

    touch("run_single_job.sh")
    run_file = open("run_single_job.sh", "w")

    write(run_file, "#!/bin/bash \n")
    write(run_file, "#SBATCH --time="*time*" \n")
    write(run_file, "#SBATCH --mem="*mem*" \n")
    write(run_file, "#SBATCH --output=outfiles/single_job_%j.out \n")
    write(run_file, "srun julia job.jl --index "*index*" \n")

    close(run_file)
end


""" Function to create a file the loops through the data and identifies
whether data files corresponding to the output of an array job are
missing. """
function create_missing_data_file(num_jobs::String)

    # remove old file
    remove_file("find_missing_data.sh")

    # write file to identify the missing data
    touch("find_missing_data.sh")
    missing_data = open("find_missing_data.sh", "w")

    write(missing_data, "#!/bin/bash \n")
    write(missing_data, "for num in {1.."*num_jobs*"} \n")
    write(missing_data, "do \n")
    write(missing_data, "\t if [[ ! -f data_\${num}.txt ]] ; then \n")
    write(missing_data, "\t \t echo \${num} \n")
    write(missing_data, "\t fi \n")
    write(missing_data, "done > missing_file_indices.txt \n")

    close(missing_data)
end


""" Function to create a bash file that merges the data files
of a given parallel computation into a single data file. """
function create_merge_file(num_jobs::String)

    # remove old file
    remove_file("merge_data.sh")

    touch("merge_data.sh")
    merge = open("merge_data.sh", "w")

    write(merge, "#!/bin/bash \n")
    write(merge, "for num in {1.."*num_jobs*"} \n")
    write(merge, "do \n")
    write(merge, "\t cat data_\${num}.txt >> final_data.txt \n")
    write(merge, "done \n")

    close(merge)
end


function create_incomplete_merge_file()


end


""" Function to process the data generated by a given array job.
Checks whether data files are missing and prompts the use to re-run
the corresponding jobs. If all partial data files are present,
the merging is done automatically and a single, final output
file is generated. """
function create_data_processing_script(time, mem)

    touch("data_processing.jl")
    data_processing = open("data_processing.jl", "w")

    # need functions
    write(data_processing, """include("parallel_computations.jl") \n""")

    # get missing indices and move back to work directory
    write(data_processing, """cd("Data_temporary") \n""")
    write(data_processing, """run(`bash find_missing_data.sh`) \n""")
    write(data_processing, """run(`mv missing_file_indices.txt ..`) \n""")
    write(data_processing, """cd("..") \n""")
    write(data_processing, """missing_indices = readlines("missing_file_indices.txt") \n""")

    # missing files/failed jobs detected
    write(data_processing, "index_count = length(missing_indices) \n")
    write(data_processing, "if index_count != 0 \n")
    write(data_processing, """println("There are \$(index_count) missing files. Re-run those (type y/n)?") \n""")
    write(data_processing, "decision = readline() \n")

    # re-run failed jobs
    write(data_processing, """if decision == "y" \n""")
    write(data_processing, "for ind in missing_indices \n")
    write(data_processing, """create_single_job_file("$(time)", "$(mem)", ind) \n""")
    write(data_processing, """run(`sbatch run_single_job.sh`) \n""")
    write(data_processing, """remove_file("run_single_job.sh") \n""")
    write(data_processing, "end \n")

    # abort
    write(data_processing, "else \n")
    write(data_processing, """println("Aborting merging.") \n""")
    write(data_processing, "return \n")
    write(data_processing, "end \n")

    # no missing files/failed jobs, merge
    write(data_processing, "else \n")
    write(data_processing, """cd("Data_temporary") \n""")
    write(data_processing, """run(`bash merge_data.sh`) \n""")
    write(data_processing, """run(`mv final_data.txt ..`) \n""")
    write(data_processing, """cd("..") \n""")

    # clean up
    write(data_processing, """remove_directory("Data_temporary") \n""")
    write(data_processing, """remove_directory("outfiles") \n""")
    write(data_processing, """remove_file("job.jl") \n""")
    write(data_processing, """remove_file("run_array_job.sh") \n""")
    write(data_processing, """remove_file("missing_file_indices.txt") \n""")
    write(data_processing, """rm("data_processing.jl") \n""")
    write(data_processing, "end \n")

    close(data_processing)
end


############################################
# Functions to execute parallel computations
############################################

""" Function to parallelise the computation of a given function
func over a multi-dimensional parameter space. Creates the corresponding
job and run files, as well as the data merging bash files, and submits
the array job to the queue. The user is required to check themselves
whether the computation has run. """
function parallelise_computation(func, file_name, time, mem, num_jobs)

    # directories for output, job file, run file for array job
    create_data_directories()
    create_job_file(func, file_name)
    create_array_job_file(time, mem, num_jobs)

    # files for easy data recombination
    cd("Data_temporary")
    create_missing_data_file(num_jobs)
    create_merge_file(num_jobs)
    cd("..")

    # this is the file to execute once the jobs have run
    create_data_processing_script(time, mem)

    # submit array job
    run(`sbatch run_array_job.sh`)

end

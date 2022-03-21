using CSV

working_path = "D:\\Projects\\numerical-boostrap\\oscillator-simple-prototype\\"

file = CSV.File(working_path * "nonlinear-SDP-x-power-11-standard-m-value.csv", header = false)
output_file = "nonlinear-SDP-x-power-11-standard-m-mat.jl"

len = length(file)
for i in 1 : len
    for j in 1 : len
        open(working_path * output_file, "a") do f
            println(f, "M[$i, $j] => $(file[i][j]),")
        end
    end
end
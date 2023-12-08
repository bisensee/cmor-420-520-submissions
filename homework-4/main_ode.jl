using OrdinaryDiffEq
using Plots

m = 100
n = 100
h = 1/n
nzvals = zeros((m-1)*3+4)
nzvals[1] = 1
nzvals[2] = -1
for t in 3:3:((m-1)*3+2)
    nzvals[t] = -1
    nzvals[t+1] = 2
    nzvals[t+2] = -1
end
nzvals[(m-1)*3+3] = -1
nzvals[(m-1)*3+4] = 1
num_vals = length(nzvals)

row_indices = zeros(Int, num_vals)
row_indices[1] = 1
row_indices[2] = 1
let ctr = 2
    for i in 3:3:(num_vals-2)
        row_indices[i] = ctr
        row_indices[i+1] = ctr
        row_indices[i+2] = ctr
        ctr += 1
    end
    row_indices[num_vals-1] = n+1
    row_indices[num_vals] = n+1
end
col_indices = zeros(Int, num_vals)
col_indices[1] = 1
col_indices[2] = 2
let ctr = 1
    for i in 3:3:(num_vals-2)
        col_indices[i] = ctr
        col_indices[i+1] = ctr+1
        col_indices[i+2] = ctr+2
        ctr += 1
    end
    col_indices[num_vals-1] = n
    col_indices[num_vals] = n+1
end
nzvals = nzvals/(h^2)
A_dense = zeros(m+1, m+1)
for idx in 1:num_vals
    A_dense[row_indices[idx], col_indices[idx]] = nzvals[idx]
end

b = zeros(n+1)
for i in 0:n
    b[i+1] = cos(pi*i*h)
end

p = (b_param = b, A_param = A_dense)


rhs(u, p, t) = p[1] - p[2]*u
interval = (0., 1.)
u0 = zeros(m+1)
prob = ODEProblem(rhs, u0, interval, p)
sol1 = solve(prob, Heun(), dt = 0.5*h^2)
(~,num_sol) = size(sol1)
step = floor(Int, sqrt(num_sol))
x = 0:100

@gif for i in 1:step:num_sol
    plot(x, sol1[i], ylims = [-0.11, 0.11], 
    title=string("Heun, Timestep: ", i, ", t: ", (i-1)/19975),
    label = "u solution",
    xlabel = "i", ylabel = "u_i")
end

sol2 = solve(prob, RK4(), dt = 0.5*h^2)
(~,num_sol) = size(sol2)
step = floor(Int, sqrt(num_sol))

@gif for i in 1:step:num_sol
    plot(x, sol2[i], ylims = [-0.11, 0.11], 
    title=string("RK4, Timestep: ", i, ", t: ", (i-1)/14350),
    label = "u solution",
    xlabel = "i", ylabel = "u_i")
end

sol3 = solve(prob, BS3(), dt = 0.5*h^2)
(~,num_sol) = size(sol3)
step = floor(Int, sqrt(num_sol))

@gif for i in 1:step:num_sol
    plot(x, sol3[i], ylims = [-0.11, 0.11], 
    title=string("BS3, Timestep: ", i, ", t: ", (i-1)/15903),
    label = "u solution",
    xlabel = "i", ylabel = "u_i")
end

println(sol1.stats)
println(sol2.stats)
println(sol3.stats)


### Using the rhs! method

# function rhs!(du, u, p, t) 
#     du = p[1] - p[2]*u
# end
# prob = ODEProblem(rhs!, u0, interval, p)
# sol1 = solve(prob, Heun(), dt = 0.5*h^2)
# (~,num_sol) = size(sol1)
# step = floor(Int, sqrt(num_sol))
# x = 0:100

# @gif for i in 1:step:num_sol
#     plot(x, sol1[i], ylims = [-0.11, 0.11], 
#     title=string("Timestep: ", i, ", t: ", (i-1)/19975),
#     label = "u solution",
#     xlabel = "i", ylabel = "u_i")
# end

# sol2 = solve(prob, RK4(), dt = 0.5*h^2)
# (~,num_sol) = size(sol2)
# step = floor(Int, sqrt(num_sol))

# @gif for i in 1:step:num_sol
#     plot(x, sol2[i], ylims = [-0.11, 0.11], 
#     title=string("Timestep: ", i, ", t: ", (i-1)/14350),
#     label = "u solution",
#     xlabel = "i", ylabel = "u_i")
# end

# sol3 = solve(prob, BS3(), dt = 0.5*h^2)
# (~,num_sol) = size(sol3)
# step = floor(Int, sqrt(num_sol))

# @gif for i in 1:step:num_sol
#     plot(x, sol3[i], ylims = [-0.11, 0.11], 
#     title=string("Timestep: ", i, ", t: ", (i-1)/15903),
#     label = "u solution",
#     xlabel = "i", ylabel = "u_i")
# end
struct SparseMatrixCSR{T} <: AbstractMatrix{T}
    m::Int
    n::Int
    row_ptr::Vector{Int}
    col_indices::Vector{Int}
    nzvals::Vector{T}
end

function SparseMatrixCSR_constructor(m, n, row_indices, col_indices, nzvals)
    row_ptr = Vector{Int}(undef, m+1)
    row_ptr[1] = 1
    ctr = 1
    comp = row_indices[1]
    idx = 2
    r = 2
    while r <= length(nzvals)
        if row_indices[r] == comp
            ctr += 1
        else
            row_ptr[idx] = row_ptr[idx-1] + ctr
            idx += 1
            ctr = 1
            comp = row_indices[r]
        end
        r += 1
    end
    row_ptr[m+1] = row_ptr[m] + ctr
    SparseMatrixCSR(m, n, row_ptr, col_indices, nzvals)
end

import Base: size, getindex

function size(A::SparseMatrixCSR)
    return A.m, A.n
end

function getindex(A::SparseMatrixCSR, i, j)
    start = A.row_ptr[i]
    stop = A.row_ptr[i+1]
    return_val = 0.
    for idx in start:(stop-1)
        if A.col_indices[idx] == j
            return_val = A.nzvals[idx]
            break
        else
            idx += 1
        end
    end
    return return_val
end

import LinearAlgebra: mul!

function mul!(out, A::SparseMatrixCSR, x::AbstractVector)
    m,_ = size(A)
    for i in 1:m
        row_start = A.row_ptr[i]
        row_end = A.row_ptr[i+1]
        sum = 0.
        for j in row_start:(row_end-1)
            col = A.col_indices[j]
            sum += getindex(A, i,col) * x[col]
        end
        out[i] = sum
    end
end

m = 100
n = 100
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
end
row_indices[num_vals-1] = n+1
row_indices[num_vals] = n+1

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
end
col_indices[num_vals-1] = n
col_indices[num_vals] = n+1

h = 1/n
a = 0.5*h^2
u = zeros(n+1)
b = zeros(n+1)

for i in 0:n
    b[i+1] = cos(pi*i*h)
end

nzvals = nzvals/(h^2)
A_sparse = SparseMatrixCSR_constructor(m+1, n+1, row_indices, col_indices, nzvals)
x = ones(m+1)
out = ones(m+1)

A_dense = zeros(m+1, m+1)
for idx in 1:length(nzvals)
    A_dense[row_indices[idx], col_indices[idx]] = nzvals[idx]
end

import LinearAlgebra.norm, LinearAlgebra.mul!

function iterate(A, b, u, a; tol = 1e-3)
    (m,n) = size(A)
    product = zeros(m)
    iter = 0
    r = b
    norm_r = norm(r,2)
    println("Iteration ", iter, ", norm_r = ", norm_r)
    while norm_r >= tol
        u += a*r
        mul!(product, A, u)
        r = b - product
        iter += 1
        norm_r = norm(r,2)
    end
    println("Iteration ", iter, ", norm_r = ", norm(r,2))
    return u, iter, norm_r;
end

println("Sparse Matrix")
final_u, num_iter, final_norm = iterate(A_sparse, b, u, a);
println()
println("Dense Matrix")
iterate(A_dense, b, u, a)
println()

@code_warntype mul!(out, A_sparse, x)

@code_warntype iterate(A_sparse, b, u, a)
@code_warntype iterate(A_dense, b, u, a)

println("Sparse Matrix")
@time iterate(A_sparse, b, u, a)
println()
println("Dense Matrix")
@time iterate(A_dense, b, u, a)
println()

using Plots

x = 1:101
initial_u = zeros(101)
plot(x, [initial_u, final_u], title="u solutions, Sparse matrix", label=["Initial u" "Final u"],
xlabel = "i", ylabel = "u_i")
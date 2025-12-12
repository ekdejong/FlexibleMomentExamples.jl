using Cloudy
using Cloudy.KernelFunctions
using Cloudy.KernelTensors
using CairoMakie
using Printf

function compute_kernel_from_tensor(c, x, y)
    n, m = size(c)
    output = 0.0
    for i in 1:n
        for j in 1:m
            output += c[i, j] * x^(i - 1) * y^(j - 1)
        end
    end
    return output
end

function compute_kernel_from_tensor_matrix(matrix_of_kernels, thresholds, x, y)
    num_kernels = length(matrix_of_kernels)
    output = 0.0
    for xk in 1:num_kernels
        for yk in 1:num_kernels
            # find the correct kernel index
            if x >= thresholds[xk] && x < thresholds[xk + 1]
                if y >= thresholds[yk] && y < thresholds[yk + 1]
                    c = matrix_of_kernels[xk][yk].c
                    n, m = size(c)
                    for i in 1:n
                        for j in 1:m
                            output += c[i, j] * x^(i - 1) * y^(j - 1)
                        end
                    end
                end
            end
        end
    end
    return output
end

function get_z_limits(z_kernel, z_tensor)
    min_lim = min(
        minimum(skipmissing(isnan(x) ? missing : x for x in z_kernel)),
        minimum(skipmissing(isnan(x) ? missing : x for x in z_tensor)),
    )
    max_lim = max(
        maximum(skipmissing(isnan(x) ? missing : x for x in z_kernel)),
        maximum(skipmissing(isnan(x) ? missing : x for x in z_tensor)),
    )
    return min_lim, max_lim
end

# Set up the figure with nice aesthetics
set_theme!(Theme(
    fontsize = 16,
    font = "Latin Modern Roman",
    Axis = (
        xgridvisible = false,
        ygridvisible = false,
        topspinevisible = false,
        rightspinevisible = false,
    ),
    Colorbar = (
        ticklabelsize = 14,
    )
))

# HydrodynamicKernelFunction
FT = Float64
limit = FT(1e-6)
order = 4

kernel_func = HydrodynamicKernelFunction(1e2 * π)
kernel_tensor = CoalescenceTensor(kernel_func, order, limit)

n = 100
x = range(0, limit, n)
y = range(0, limit, n)
z_kernel = zeros(n, n)
z_tensor = zeros(n, n)
for i in 1:n
    for j in 1:n
        z_kernel[i, j] = kernel_func(x[i], y[j])
        z_tensor[i, j] = compute_kernel_from_tensor(kernel_tensor.c, x[i], y[j])
    end
end

# LongKernelFunction
limit2 = FT(1e-9)
order2 = 2

kernel_func2 = LongKernelFunction(5.236e-10, 9.44e9, 5.78) # 5.236e-10 kg; 9.44e9 m^3/kg^2/s; 5.78 m^3/kg/s
matrix_of_kernels = ntuple(2) do i
    ntuple(2) do j
        if i == j == 1
            CoalescenceTensor(kernel_func2, 2, FT(5e-10))
        else
            CoalescenceTensor(kernel_func2, 2, FT(1e-6), FT(5e-10))
        end
    end
end

n2 = 100
x2 = range(0, limit2, n)
y2 = range(0, limit2, n)
z_kernel2 = zeros(n, n)
z_tensor2 = zeros(n, n)
for i in 1:n
    for j in 1:n
        z_kernel2[i, j] = kernel_func2(x2[i], y2[j])
        z_tensor2[i, j] =
            compute_kernel_from_tensor_matrix(matrix_of_kernels, (FT(0), FT(5e-10), FT(1e-6)), x2[i], y2[j])
    end
end

# Create the figure
fig = Figure(size = (800, 600), backgroundcolor = :white)

# Plot 1: Exact Hydrodynamic
ax1 = Axis(fig[1, 1],
    xlabel = L"m_1 \text{ (ng)}",
    ylabel = L"m_2 \text{ (ng)}",
    title = "Exact Hydrodynamic\nKernel Function",
    titlesize = 18,
)
min_lim1, max_lim1 = get_z_limits(z_kernel, z_tensor)
hm1 = contourf!(ax1, x * 1e9, y * 1e9, z_kernel,
    levels = 20,
    colormap = :viridis,
    colorscale = (min_lim1, max_lim1)
)
cb1 = Colorbar(fig[1, 2], 
    colorrange = (min_lim1, max_lim1),
    label = L"K \, (m^3/s)"
)

# Plot 2: Approximated Hydrodynamic
ax2 = Axis(fig[1, 3],
    xlabel = L"m_1 \text{ (ng)}",
    ylabel = L"m_2 \text{ (ng)}",
    title = "Approximated Hydrodynamic\nKernel Tensor",
    titlesize = 18,
)
hm2 = contourf!(ax2, x * 1e9, y * 1e9, z_tensor,
    levels = 20,
    colormap = :viridis,
    colorscale = (min_lim1, max_lim1)
)
cb2 = Colorbar(fig[1, 4], 
    colorrange = (min_lim1, max_lim1),
    label = L"K \, (m^3/s)",
)

# Plot 3: Exact Long 1974
ax3 = Axis(fig[2, 1],
    xlabel = L"m_1 \text{ (ng)}",
    ylabel = L"m_2 \text{ (ng)}",
    title = "Exact Long 1974\nKernel Function",
    titlesize = 18,
)
min_lim2, max_lim2 = get_z_limits(z_kernel2, z_tensor2)
hm3 = contourf!(ax3, x2 * 1e9, y2 * 1e9, z_kernel2,
    levels = 20,
    colormap = :viridis,
    colorscale = (min_lim2, max_lim2)
)
cb3 = Colorbar(fig[2, 2], 
    colorrange = (min_lim2, max_lim2),
    label = L"K \, (m^3/s)"
)

# Plot 4: Approximated Long 1974
ax4 = Axis(fig[2, 3],
    xlabel = L"m_1 \text{ (ng)}",
    ylabel = L"m_2 \text{ (ng)}",
    title = "Approximated Long 1974\nKernel Tensor",
    titlesize = 18,
)
hm4 = contourf!(ax4, x2 * 1e9, y2 * 1e9, z_tensor2,
    levels = 20,
    colormap = :viridis,
    colorscale = (min_lim2, max_lim2)
)
cb1 = Colorbar(fig[2, 4], 
    colorrange = (min_lim2, max_lim2),
    label = L"K \, (m^3/s)"
)

# Adjust layout spacing
colgap!(fig.layout, 1, 10)
colgap!(fig.layout, 2, 30)
colgap!(fig.layout, 3, 10)
rowgap!(fig.layout, 20)

# Save the figure
save(joinpath(dirname(@__DIR__), "figures/figA1/figA1_KernelFunction_Approximation.pdf"), 
    fig, pt_per_unit = 2)
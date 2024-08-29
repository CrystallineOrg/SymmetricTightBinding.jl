## Packages
using JuMP
using Gurobi

## This function will return an array of named tuples in the variables with objective value.
function find_all_long_modes(vᵀ; brs; vᴸ; verbose::Bool=true)
    ## Model
    model = Model(Gurobi.Optimizer)

    @variable(model, x[i=1:lenght(brs)-1] >= 0)
    @constraint(model, A´ * x == vᵀ + vᴸ)
    if verbose
        print(model)
    end

    ## Gurobi parameters - see https://www.gurobi.com/documentation/9.0/refman/finding_multiple_solutions.html
    JuMP.set_optimizer_attribute(model, "PoolSearchMode", 1) # set to 2 if want a more sistematic search
    JuMP.set_optimizer_attribute(model, "PoolSolutions", 100)  # 100 is an arbitrary (large enough) whole number

    if !verbose
        JuMP.set_optimizer_attribute(model, "OutputFlag", 0)
    end

    ## Optimize:
    JuMP.optimize!(model)

    ## Results
    num_results = result_count(model)
    if verbose
        println("Num of results: ", num_results)
    end

    Results = Array{NamedTuple{(:x, :y, :u, :v, :obj),NTuple{5,Int64}},1}() ## Note the conversion to Int64
    for i in 1:num_results
        sol = (x=JuMP.value(x; result=i), y=JuMP.value(y; result=i),
            u=JuMP.value(u; result=i), v=JuMP.value(v; result=i),
            obj=JuMP.objective_value(model; result=i))
        push!(Results, sol)
        if verbose
            println(Results[i])
        end
    end

    return Results
end
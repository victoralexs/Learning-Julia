# Example Pablo Guerron-Quintana - Sovereign Default

# Setting working Directory
cd("/Users/victoralexandrino/Google Drive/PhD Insper/Thesis/Paper 1 - Sovereign Default Holdings/Quant/Learning Julia")

pwd()

# Clear workspace

clearconsole()

# For doubts in Environments, always use quantecon-netbooks environment, check https://julia.quantecon.org/more_julia/tools_editors.html

using Random, Distributions
#Initialization

function tauchen(ρ, σ, Ny, P)
    #Create equally spaced pts to fill into Z
    σ_z = sqrt((σ^2)/(1-ρ^2))
    Step = 10*σ_z/(Ny-1)
    Z = -5*σ_z:Step:5*σ_z

    #Fill in entries of 1~ny, ny*(ny-1)~ny^2
    for z in 1:Ny
        P[z,1] = cdf(Normal(), (Z[1]-ρ*Z[z] + Step/2)/σ)
        P[z,Ny] = 1 - cdf(Normal(),(Z[Ny] - ρ*Z[z] - Step/2)/σ)
    end

    #Fill in the middle part
    for z in 1:Ny
        for iz in 2:(Ny-1)
            P[z,iz] = cdf(Normal(), (Z[iz]-ρ*Z[z]+Step/2)/σ) - cdf(Normal(), (Z[iz]-ρ*Z[z]-Step/2)/σ)
        end
    end
end

# Main VFI implementation for sovereign default

# Initialize functin main()
function main()

# Parameters

Ny = 7              # Number of grids for output process
Nb = 100            # Number of grids for bonds
maxIter = 500       # Number os max iterations
maxInd = Ny * Nb    # Total number of states
rstar = 0.017       # risk free interest rate
lbd = -1            # Lower bound for debt grid
ubd = 0             # Upper bound for debt grid
rrisk = 0.5         # Risk aversion
β = 0.953           # SDF
τ = 0.15            # Output cost
θ = 0.282           # Probability of re-entry
tol = 1e-10         # Tolerance for VFI
ϕ = rrisk           # Risk aversion
δ = 0.8             # Check convergence parameter
ρ = 0.9             # Output autocorrelation
σ = 0.025           # Output volatility
τ = 0.5             # Output cost

# Initial values

B = zeros(Nb)               # Initial bond guess
Y = zeros(Ny)               # Initial output guess
σ_z = sqrt((σ^2)/(1-ρ^2))   # Initial volatility of output
Step = 10*σ_z/(Ny-1)        # Initial gap output discretization
Y = -5*σ_z:Step:5*σ_z       # Initial output = equals output grid


P = zeros(Ny,Ny)                        # Initial transition matrix
V = fill(1/((1-β)*(1-rrisk)),Ny, Nb)    # Initial guess for V, equals to 1/((1-β)*(1-rrisk)) with a dimension of Ny x Nb
Price = fill(1/(1+rstar),Ny, Nb)        # Initial guess for bond price, with a dimension of Ny x Nb
Vr = zeros(Ny, Nb)                      # Initial guess for value of repayment, with a dimention of Ny x Nb
Vd = zeros(Ny)                          # Initial guess for value of default, dimension Ny
decision = ones(Ny,Nb)                  # Initial guess for prob of default, dimension Ny x Nb

U(x) = x^(1-ϕ) / (1-ϕ)                  # Utility function, x is consumption. Do not need to call function argument

#Initialize Bond grid: Discretize bond
minB = lbd                              # Lower bound bond grid
maxB = ubd                              # Upper bound bond grid
step = (maxB-minB) / (Nb-1)             # Equalized space for bond grid
B = minB:step:maxB                      # Creating bond grid

#Initialize Shock grid
tauchen(ρ, σ, Ny, P)                    # Call exogenous process. Result is P, a Ny x Ny matrix with probabilities
sumdef = 0                              # LHS of value of default

err = 2000                              # Initial error
tol = 1e-6                              # Initial tolerance
iter = 0                                # Initial tolerance

time_vd = 0                             # ?
time_vr = 0                             # ?
time_decide = 0                         # ?

# Starting VFI
#3
while (err > tol) & (iter < maxIter)
    V0 = deepcopy(V)                    # Using V as initial value
    Vd0 = deepcopy(Vd)                  # Using Vd as initial value
    Price0 = deepcopy(Price)            # Using Price as initial value
    prob = zeros(Ny,Nb)                 # Initializing probability of default
    # display(V0)

#5
    for ib in 1:Nb                      # Initialize iteration on bonds
        for iy = 1:Ny                   # Initialize iteration on output


    # Compute default and repayment
    #7

    # Value of default

            sumdef = U(exp((1-τ)*Y[iy]))                                # LHS utility of default, with output cost
            for y in 1:Ny                                               # Here, he defines a new integer "y", which corresponds to next period y
                sumdef += (β* P[iy,y]* (θ* V0[y,1] + (1-θ)* Vd0[y]))    # Interating in output the RHS of value of default
            end
            Vd[iy] = sumdef                                             # Calling total value of default as sumdef

    # Compute value of repayment
    #8

            Max = -Inf                                                  # Calling max as -∞
            for b in 1:Nb                                               # Iterating on bonds
                c = exp(Y[iy]) + B[ib] - Price0[iy,b]*B[b]              # Defining budget constraint
                if c > 0
                    sumret = 0                                          # Utility in repayment
                    for y in 1:Ny                                       # Iterating on output
                        sumret += P[iy,y]*V0[y,b]                       # Writing the total value of repayment
                    end
                    vr = U(c) + β * sumret                              # Defining value of default
                    Max = max(Max, vr)                                  # Assuring that vr > ∞
                end
            end
            Vr[iy,ib] = Max                                             # Writing vale as the max utility after iteration


            #Choose repay or default
            if (Vd[iy] < Vr[iy,ib])                                     # If value of default is lower...
                V[iy,ib] = Vr[iy,ib]                                    # It chooses to repay
                decision[iy,ib] = 0                                     # DDefault decision equals zero
            else
                V[iy,ib] = Vd[iy]                                       # If not, it chooses to default
                decision[iy,ib] = 1                                     # Default decision equals one
            end

            #calculate debt price

            for y in 1:Ny                                               # Iterating again in y, we calculate the debt price
                prob[iy,ib] += P[iy,y] * decision[y,ib]                 # Probability of default is also stochastic, so we should iterate prob() in iy and ib
            end
            Price[iy,ib] = (1-prob[iy,ib]) / (1+rstar)                  # Writing bond price


        end
    end


    err = maximum(abs.(V-V0))                                           # Convergence on continuation value
    PriceErr = maximum(abs.(Price-Price0))                              # Convergence on bond price
    VdErr = maximum(abs.(Vd-Vd0))                                       # Convergence on value of default
    Vd = δ * Vd + (1-δ) * Vd0                                           # Average value of default
    Price = δ * Price + (1-δ) * Price0                                  # Average bond price
    V = δ * V + (1-δ) * V0                                              # Average continuation value
    iter = iter + 1                                                     # Iterations
    println("Errors of round $iter: Value error: $err, price error: $PriceErr, Vd error: $VdErr") # Print number of iterations and errors in value function, bond prices and value of default

end


println("Total Round ",iter)                                            # Print number of iterations after convergence

Vd = Vd[:,:]                                                            # Saving value of default

# Print everything

println("Vr: ====================")
display(Vr)
println("Vd: ==================")
display(Vd)
println("Decision: ==================")
display(decision)
println("Price: ==================")
display(Price)

return Vr,Vd,decision,Price                                            # Return of main() function

end

@time VReturn, VDefault, Decision, Price = main()                      # Plot time

#= Storing as CSV
using Parsers
using DataFrames
using CSV
dfPrice = DataFrame(Price)
dfVr = DataFrame(VReturn)
dfVd = DataFrame(VDefault)
dfDecision = DataFrame(Decision)
CSV.write("/Users/deng/Desktop/school/ECON8873/codes/Price.csv", dfPrice)
CSV.write("/Users/deng/Desktop/school/ECON8873/codes/Vr.csv", dfVr)
CSV.write("/Users/deng/Desktop/school/ECON8873/codes/Vd.csv", dfVd)
CSV.write("/Users/deng/Desktop/school/ECON8873/codes/Decision.csv", dfDecision)
=#

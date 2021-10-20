# -----------------------------------------------------
# Simple functions to compute trends and uncertainties
# from AR1 model, so no Hector dependency for NOAA TR
# -----------------------------------------------------
# Work in progress
using Statistics
using LinearAlgebra

function EstTrend(tsteps,tseries;est_acc=false)
    # Only retain valid time steps
    acc_idx = @. isfinite(tseries)
    t = tsteps[acc_idx]
    ts = tseries[acc_idx]

    # Estimate residuals
    # est_acc ? amat = ones(sum(acc_idx),3) : amat = ones(sum(acc_idx),2) 
    # @. amat[:,2] = t - $mean(t)
    # est_acc ? (@. amat[:,3] = 0.5 * (t - $mean(t))^2) : nothing
    # β = amat\ts
    # resid = ts - amat*β

    # # First-guess covariance matrix
    # ρ = cor(resid[1:end-1],resid[2:end])
    # FillCovmat!(Σ,ρ)
    β,ϵ,covmat = GLS(amat,ts,Σ)

    sol=Dict()
    sol["bias"] = β[1]
    sol["trend"] = β[2]
    est_acc ? sol["accel"] = β[3] : nothing

    sol["bias_sigma"] = ϵ[1]
    sol["trend_sigma"] = ϵ[2]
    est_acc ? sol["accel_sigma"] = ϵ[3] : nothing
    sol["covmat"] = covmat
    return sol
end

function FillCovmat!(Σ,ρ)
    Σ[diagind(Σ)] .= 1.0
    for n ∈ 1:size(Σ,1)
        Σ[diagind(Σ,n)] .= ρ^n
        Σ[diagind(Σ,-n)] .= ρ^n
    end
    return nothing
end

function GLS(amat,ts,Σ)
    β = zeros(size(amat,2))
    resid = zeros(size(amat,1))
    covmat = zeros(size(amat,2),size(amat,2))
    Σ = zeros(sum(acc_idx),sum(acc_idx)) 
    for i ∈ 1:5
        ρ = cor(resid[1:end-1],resid[2:end])
        FillCovmat!(Σ,ρ)
        covmat[:] = inv(amat'*inv(Σ)*amat)
        β[:] = covmat*amat'*inv(Σ)*ts
        resid[:] = ts - amat*β
        ρ = cor(resid[1:end-1],resid[2:end])
        FillCovmat!(Σ,ρ)
    end
    σ2 = (resid'*resid) / (length(ts) - size(amat,2))
    ϵ = sqrt.(diag(σ2 .* covmat))
    return β,ϵ,covmat
end

function GenNoise()
end


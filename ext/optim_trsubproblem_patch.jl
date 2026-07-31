# ext/optim_trsubproblem_patch.jl
#
# TEMPORARY workaround for an Optim.jl bug in `solve_tr_subproblem!` (the trust-region
# subproblem solver used by `NewtonTrustRegion`, our default `fit` optimizer). For a
# (near-)singular Hessian, `lambda_lb = nextfloat(-min_H_ev)` is a negligible ridge, so the
# Cholesky root-finding loop never factorizes within `max_iters` and leaves the step `s`
# unwritten (NaN) — which then crashes the next `eigen!` in `fgh!`. Our moment seed for
# multi-EBR models (e.g. SG 221) can land exactly on such a singular-Hessian point, so `fit`
# fails reliably on some BLAS/LAPACK builds (observed on CI, Julia 1.12, ubuntu-latest).
#
# Fix (submitted upstream): leave `lambda_lb` untouched and, only when Cholesky fails, boost
# λ to the spectral scale of H so the next factorization succeeds (`max(2λ, √eps·‖H‖)`); the
# root-finder can still descend to a smaller optimal λ afterwards, so well-conditioned solves
# are unaffected.
#
# This shim redefines the method with that fix. It self-disables once a fixed Optim.jl is
# installed: it probes the documented failing input and only patches if the bug is present,
# so it becomes a no-op (and can be deleted) after the upstream fix is released.
#
# Fix for Optim.jl in PR #1266 (https://github.com/JuliaNLSolvers/Optim.jl/pull/1266):
# delete this once that is merged + released + compat-updated here.
#
# Included from `SymmetricTightBindingOptimExt.jl` (so every `fit`/`multistart_fit` call is
# protected, not just the benchmark harness); `benchmark/fit_benchmark.jl` relies on the same
# copy via that extension load, rather than including this file itself.
#
# Wrapped in a function called from the extension's `__init__` rather than run as top-level
# module code: `@eval`-ing into another (closed) module is disallowed during precompilation,
# and extensions — unlike plain scripts — are precompiled.
function __patch_optim_tr_subproblem!()
    _has_bug = if !isdefined(Optim, :solve_tr_subproblem!)
        false
    else
        s = fill(NaN, 2)
        try
            # H = [1 1; 1 1] is positive-semidefinite but singular (eigenvalues 0, 2); the
            # buggy solver leaves `s` at its incoming NaN for this input.
            Optim.solve_tr_subproblem!([1.0, 1.0], [1.0 1.0; 1.0 1.0], 1.0, s)
            any(isnan, s)
        catch
            false # signature changed / errors ⇒ don't attempt to patch
        end
    end

    if _has_bug
        @info "SymmetricTightBinding: applying temporary Optim `solve_tr_subproblem!` \
               patch (near-singular-Hessian NaN bug); remove once a fixed Optim.jl is \
               released (cf. JuliaNLSolvers/Optim.jl#1266)"
        @eval Optim begin
            function solve_tr_subproblem!(gr, H, delta, s; tolerance = 1e-10, max_iters = 5)
                T = eltype(gr)
                n = length(gr)
                delta_sq = delta^2

                @assert n == length(s)
                @assert (n, n) == size(H)
                @assert max_iters >= 1

                Hsym = Symmetric(H)
                if any(!isfinite, Hsym)
                    return T(Inf), false, zero(T), false, false
                end
                H_eig = eigen(Hsym)

                if !isempty(H_eig.values)
                    min_H_ev, max_H_ev = H_eig.values[1], H_eig.values[n]
                else
                    return T(Inf), false, zero(T), false, false
                end
                H_scale = max(abs(min_H_ev), abs(max_H_ev)) # ← PATCH: spectral scale of H
                H_ridged = copy(H)

                qg = H_eig.vectors' * gr

                interior = true
                hard_case = false
                reached_solution = true

                if min_H_ev >= 1e-8
                    calc_p!(zero(T), 1, n, qg, H_eig, s)
                end

                if min_H_ev >= 1e-8 && sum(abs2, s) <= delta_sq
                    interior = true
                    reached_solution = true
                    lambda = zero(T)
                else
                    interior = false

                    hard_case_candidate, min_i = check_hard_case_candidate(H_eig.values, qg)

                    lambda_lb = nextfloat(-min_H_ev)
                    lambda = lambda_lb

                    hard_case = false
                    if hard_case_candidate
                        calc_p!(lambda, min_i, n, qg, H_eig, s)
                        p_lambda2 = sum(abs2, s)
                        if p_lambda2 > delta_sq
                        else
                            hard_case = true
                            reached_solution = true
                            tau = sqrt(delta_sq - p_lambda2)
                            calc_p!(lambda, min_i, n, qg, H_eig, s)
                            LinearAlgebra.axpby!(tau, view(H_eig.vectors, :, 1), -1, s)
                        end
                    end

                    lambda = initial_safeguards(H, gr, delta, lambda)

                    if !hard_case
                        reached_solution = false
                        for iter = 1:max_iters
                            lambda_previous = lambda

                            for i in diagind(H_ridged)
                                H_ridged[i] = H[i] + lambda
                            end

                            F = cholesky(Hermitian(H_ridged), check = false)
                            if !issuccess(F)
                                # ← PATCH: a negligible λ never becomes numerically PD by
                                # doubling; jump to a ridge on the scale of H so the next
                                # factorization succeeds (root-finder may descend afterwards).
                                lambda = max(2 * lambda, sqrt(eps(T)) * H_scale)
                                continue
                            end

                            R = F.U
                            s[:] = -R \ (R' \ gr)
                            q_l = R' \ s
                            norm2_s = dot(s, s)
                            lambda_update =
                                norm2_s * (sqrt(norm2_s) - delta) / (delta * dot(q_l, q_l))
                            lambda += lambda_update

                            if lambda < lambda_lb
                                lambda = (lambda_previous + lambda_lb) / 2
                            end

                            if abs(lambda - lambda_previous) < tolerance
                                reached_solution = true
                                break
                            end
                        end
                    end
                end

                m = dot(gr, s) + dot(s, H, s) / 2

                return m, interior, lambda, hard_case, reached_solution
            end
        end
    end
end

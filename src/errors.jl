"""
    ConvergenceError(msg)

Exception raised when an iterative or numerical procedure fails to reach the
requested accuracy: an ARPACK eigensolve that does not converge, a stationary
distribution solver that exhausts its iteration budget, or a configuration-null
edge-switching chain that cannot complete the requested number of swaps under the
`swap_shortfall=:error` policy. Input validation failures keep raising
`ArgumentError` or `DimensionMismatch`.
"""
struct ConvergenceError <: Exception
    msg::String
end

Base.showerror(io::IO, error::ConvergenceError) =
    print(io, "ConvergenceError: ", error.msg)

using Aqua
using CorePeriphery
using Test

# Aqua.jl quality assurance (ambiguities, unbound type parameters, undefined
# exports, project extras, stale and compat-less dependencies, type piracy,
# persistent tasks). The full ambiguity check passes with Arpack loaded, so no
# check is disabled.
@testset "Aqua" begin
    Aqua.test_all(CorePeriphery)
end

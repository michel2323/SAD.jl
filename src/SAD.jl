module SAD
export start_context, end_context, current_context
export register, get_tangent, interpret
export Dual
mutable struct Dual{T} <: Real where {T <: Real}
    val::T
    idx::Int
    function Dual(val::T) where {T <: Real}
        return new{T}(val, 0)
    end
end

struct Op
    partials::Tuple
    args::Tuple
    output
end

mutable struct TapeContext
    opstack::Vector{Op}
    idx::Int
    TapeContext() = new(Op[], 1)
end

const ACTIVE_CONTEXT = Ref{Union{Nothing, TapeContext}}(nothing)

start_context() = (ACTIVE_CONTEXT[] = TapeContext())
end_context() = (ctx = ACTIVE_CONTEXT[]; ACTIVE_CONTEXT[] = nothing; ctx)
current_context() = ACTIVE_CONTEXT[]

import Base: +, *

using ChainRules
import ChainRules: rrule

for op in (:+, :*)
    @eval begin
        function $op(a::Dual, b::Dual)
            ctx = current_context()
            if ctx !== nothing
                y = Dual($op(a.val, b.val))
                g = rrule($op, a.val, b.val)
                partials = g[2](1.0)[2:end]
                op = Op(partials, (a.idx, b.idx), y.val)
                push!(ctx.opstack, op)
                y.idx = ctx.idx
                ctx.idx += 1
                return y
            else
                return Dual($op(a.val, b.val))
            end
        end
    end
end

function register(x::Real)
    ctx = current_context()
    dx = Dual(x)
    y = Dual(dx.val)
    y.idx = ctx.idx
    g = rrule(identity, dx.val)
    partials = g[2](1.0)[2:end]
    op = Op(partials, (dx.idx,), y.val)
    push!(ctx.opstack, op)
    ctx.idx += 1
    return y
end

function register(x::Dual)
    ctx = current_context()
    y = Dual(x.val)
    g = rrule(identity, x.val)
    partials = g[2](1.0)[2:end]
    op = Op(partials, (x.idx,), y.val)
    push!(ctx.opstack, op)
    @show ctx.idx
    y.idx = ctx.idx
    ctx.idx += 1
    return y
end

function get_tangent(x::Dual, tape::Vector{Float64})
    return tape[x.idx]
end

function interpret(ctx::TapeContext)
    nops = length(ctx.opstack)
    values_tape = zeros(ctx.idx-1)
    values_tape[end] = 1.0
    for idx in reverse(1:ctx.idx-1)
        partials = ctx.opstack[idx].partials
        args = ctx.opstack[idx].args
        for i in eachindex(args)
            if args[i] == 0
                continue
            end
            partial = partials[i]
            values_tape[args[i]] += partial*values_tape[idx]
        end
    end
    return values_tape
end

end # module SAD

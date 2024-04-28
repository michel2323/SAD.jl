module SAD
using LinearAlgebra
export start_context, end_context, current_context
export register, get_tangent, interpret
export Dual
include("functions.jl")
export SAD_FUNCTIONS
mutable struct Dual{T} <: Real where {T <: Real}
    val::T
    idx::UInt64
    function Dual(val::T) where {T <: Real}
        return new{T}(val, 0)
    end
end

Base.zero(::Type{<:Dual}) = Dual(0.0)
Base.iszero(x::Dual) = x.val == 0.0
Base.one(::Type{<:Dual}) = Dual(1.0)
Base.isfinite(x::Dual) = isfinite(x.val)
Base.:<(x::Dual, y::Dual) = x.val < y.val
Base.:>(x::Dual, y::Dual) = x.val > y.val
Dual{Float64}(x::Bool) = Dual(x ? 1.0 : 0.0)

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

for op in SAD_FUNCTIONS
    for T in (Dual,)
        @eval begin
            function $op(args::Vararg{$T,N}) where {N}
                values = getfield.(args, :val)
                ctx = current_context()
                if ctx !== nothing
                    idxs = getfield.(args, :idx)
                    y = Dual($op(values...))
                    g = rrule($op, values...)
                    partials = g[2](1.0)[2:end]
                    op = Op(partials, idxs, y.val)
                    push!(ctx.opstack, op)
                    return y
                else
                    # return Dual($op(a.val, b.val))
                    return Dual($op(values))
                end
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

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
end

mutable struct TapeContext
    opstack::Vector{Op}
    idx::UInt64
    function TapeContext(n::UInt64)
        opstack = Vector{Op}(undef, n)
        sizehint!(opstack, n)
        return new(opstack, 1)
    end
end

function Base.show(io::IO, ctx::TapeContext)
    println("Partials, Args")
    opstack = ctx.opstack
    for i in 1:ctx.idx-1
        println(io, "$(opstack[i].partials)\t, $(opstack[i].args)")
    end
end

const ACTIVE_CONTEXT = Ref{Union{Nothing, TapeContext}}(nothing)

start_context(n::UInt64) = (ACTIVE_CONTEXT[] = TapeContext(n))
end_context() = (ctx = ACTIVE_CONTEXT[]; ACTIVE_CONTEXT[] = nothing; ctx)
current_context() = ACTIVE_CONTEXT[]

import Base: +, *

using ChainRules
import ChainRules: rrule

for op in SAD_FUNCTIONS
    for T in (Dual,)
        @eval begin
            function $op(args::Vararg{$T,N}) where {N}
                @inbounds begin
                    values = getfield.(args, :val)
                    ctx = current_context()
                    if ctx !== nothing
                        idxs = getfield.(args, :idx)
                        y = Dual($op(values...))
                        g = rrule($op, values...)
                        partials = g[2](1.0)[2:end]
                        op = Op(partials, idxs)
                        ctx.opstack[ctx.idx] = op
                        y.idx = ctx.idx
                        ctx.idx += 1
                        return y
                    else
                        return Dual($op(values...))
                    end
                end
            end
        end
    end
end

function register(x::Real)
    ctx = current_context()
    dx = Dual(x)
    y = Dual(dx.val)
    g = rrule(identity, dx.val)
    partials = g[2](1.0)[2:end]
    op = Op(partials, (dx.idx,))
    ctx.opstack[ctx.idx] = op
    y.idx = ctx.idx
    ctx.idx += 1
    return y
end

function register(x::Dual)
    ctx = current_context()
    y = Dual(x.val)
    g = rrule(identity, x.val)
    partials = g[2](1.0)[2:end]
    op = Op(partials, (x.idx,))
    ctx.opstack[ctx.idx] = op
    y.idx = ctx.idx
    ctx.idx += 1
    return y
end

function get_tangent(x::Dual, tape::Vector{Float64})
    return tape[x.idx]
end

function interpret(ctx::TapeContext, dependents::Dual...)
    nops = length(ctx.opstack)
    adjoint_tape = zeros(ctx.idx-1)
    for dependent in dependents
        adjoint_tape[dependent.idx] = 1.0
    end
    for idx in reverse(1:ctx.idx-1)
        partials = ctx.opstack[idx].partials
        args = ctx.opstack[idx].args
        for i in eachindex(args)
            if args[i] == 0
                continue
            end
            partial = partials[i]
            adjoint_tape[args[i]] += partial*adjoint_tape[idx]
        end
    end
    return adjoint_tape
end

end # module SAD

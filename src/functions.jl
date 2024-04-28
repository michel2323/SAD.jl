# We want to support only operators and functions on reals

SAD_BASE_FUNCTIONS = (
    :+, :-, :/, :*, :^, :sin, :cos, :tan, :asin,
    :acos, :atan, :sinh, :cosh, :tanh, :asinh,
    :acosh, :atanh, :exp, :log, :log2, :log10,
    :sqrt, :abs, :sign, :floor, :ceil, :round,
    :trunc, :max, :min, :transpose, :adjoint,
    :factorial, :mod, :rem, :gcd, :lcm, :hypot,
    :inv,
)
SAD_LA_FUNCTIONS = (
    :dot, :cross, :norm,
    :normalize, :det,
    :pinv, :rank, :eigvals,
    :eigvecs, :svd, :svdvals,
    :qr, :lu, :ldlt,
    :cholesky,
    :hessenberg, :schur
)

for f in SAD_BASE_FUNCTIONS
    @eval begin
        import Base: $f
    end
end

for f in SAD_LA_FUNCTIONS
    @eval begin
        import LinearAlgebra: $f
    end
end

SAD_FUNCTIONS = (SAD_BASE_FUNCTIONS..., SAD_LA_FUNCTIONS...)

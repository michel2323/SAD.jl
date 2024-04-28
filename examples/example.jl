using SAD

function f(x,y)
    z = x * y
    return z
end

start_context()  # Start tracking operations
x = register(3.0)
y = register(2.0)
z = f(x,y)
dz = register(z)
ctx = end_context()  # Stop tracking and retrieve the context
tape = interpret(ctx,dz);  # Assuming the seed gradient of the output is 1
@show dx = get_tangent(x,tape)
@show dy = get_tangent(y,tape)

# Testing inv(A)
start_context()  # Start tracking operations
A = rand(2, 2)
dA = register.(A)
dinvA = inv(dA)
dy = register(dinvA[1,1])
ctx = end_context()  # Stop tracking and retrieve the context
tape = interpret(ctx, dy);  # Assuming the seed gradient of the output is 1
ad = get_tangent(dA[1,1],tape)

h = 1e-8
hA = copy(A)
hA[1,1] += h

invA = inv(A)
hinvA = inv(hA)

fd = (hinvA[1,1] - invA[1,1])/h

isapprox(ad, fd, rtol=1e-4)

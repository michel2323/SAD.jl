using SAD

function f(x,y)
    z = x * y
    return z
end

start_context()  # Start tracking operations
x = register(3.0)
y = register(2.0)
z = f(x,y)
register(z)
ctx = end_context()  # Stop tracking and retrieve the context
tape = interpret(ctx);  # Assuming the seed gradient of the output is 1
@show dx = get_tangent(x,tape)
@show dy = get_tangent(y,tape)

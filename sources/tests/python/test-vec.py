import alta

fail = False

# Init the vector using the
x = alta.vec([1.0, 1.0, 1.0])
if x[0] != 1.0:
    fail = True

# Vector arithmetic
y = x + x
z = x - x
if y[0] != 2.0 and z[0] != 0.0:
    fail = True

print "Testing 'str' method on variable 'y': " + str(y)

# Print the test result
if fail:
    print "Testing python 'vec' interface failed!"
    exit(1)
else:
    print "Testing python 'vec' interface passed!"
    exit(0)

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

# Testing str(vec)
print "Testing 'str' method on variable 'y': " + str(y)

# Testing the length of the vector
if len(y) != 3:
    fail = True

# Print the test result
if fail:
    print "Testing python 'vec' interface failed!"
    exit(1)
else:
    print "Testing python 'vec' interface passed!"
    exit(0)

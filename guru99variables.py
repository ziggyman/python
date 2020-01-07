a=100
print(a)

# Declare a variable and initialize it
f = 0
print(f)
print('f is of type ',type(f))
# re-declaring the variable works
f = 'guru99'
print(f)
print('f is of type ',type(f))

a = "Guru"
b = 99
c = a + str(b)
print(c)

# Declare a variable and initialize it
f = 101
print(f)

# Global vs. local variables in functions
def someFunction():
# global f
    f = 'I am learning Python'
    print(f)

someFunction()
print(f)

# Global vs.local variables in functions
def someFunction():
  global f
  print(f)
  f = "changing global variable"

someFunction()
print(f)

del(f)
print(f)

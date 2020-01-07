def plus(x, y):
    return x + y

print('Hello World!')
print("Hello World!")

a = 1
print(a)
print('a = ',a)
print('a has the type ',type(a))
print(type(a),' offers the following in-built functions: ',dir(a))
print('a.bit_length() = ',a.bit_length())

b = 2
c = a + b
print('c = ',c)
# a, b, c are of type integers

str = 'Hello World!'
strA = 'I am Martin'
# str is of type string
print(str)
print('str has the type ',type(str))
print(type(str),' offers the following in-built functions: ',dir(str))
print('str.upper() = ',str.upper())

d = plus(a, b)
print('d = ',d)

strB = plus(str, strA)
print(strB)

strB = plus(plus(str,' '), strA)
print(strB)

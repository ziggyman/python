var1 = "Guru99!"
var2 = "I am Martin. I am learning Python!"
print ("var1[0] = <",var1[0],">")
print ("var1[0] = <"+var1[0]+">")
print ("var2[1:5] = ",var2[1:5])

x="Guru"
print(x[0])# first element = G
print(x[1])# second element = u

y = x[1:len(x)]
print(y)

print(x[1:3])
print("u" in x)
print("l" in x)
print("l" not in x)

name = 'guru'
number = 99
print('%s %d' % (name,number))

sentense = 'Hello my name is %s and I am %d years old'
print(sentense)
print(sentense % ('Andreas',43))
print(sentense % ('Martin',17))

x="Guru"
y="99"
print(x+y)

print(x*2)

x = "Hello World!"
print(x[:6])
print(x[0:6] + "Guru99")

# replace function for variables of type String
x = 'abc abc abc'
y = x.replace('ab','d') # y = 'dc dc dc'
print('y = <'+y+'>')
z = x.replace('ab','d',2) # z = 'dc dc abc'
print('z = <'+z+'>')

oldstring = 'I like Guru99'
print(dir(oldstring))
newstring = oldstring.replace('like', 'love')
print(oldstring)
print(newstring)
print(oldstring.swapcase())

string="python at guru99"
print(string.upper())
print(string.capitalize())

string="PYTHON AT GURu99"
print('string = <'+string+'>')
print('string.lower() = <'+string.lower()+'>')

print(":"+"Python")
print(":".join("Python"))

string="12345"
print(reversed(string))# NOTE: reversed is a free function
print(''.join(reversed(string)))# NOTE: reversed is a free function
print(dir(string))

print('type(string) = ',type(string))

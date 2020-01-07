var1 = "Guru99!"
var2 = "I am Martin. I am learning Python!"
print ("var1[0]:",var1[0])
print ("var2[1:5]:",var2[1:5])

x="Guru"
print(x[1])

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

oldstring = 'I like Guru99'
newstring = oldstring.replace('like', 'love')
print(oldstring)
print(newstring)
print(oldstring.swapcase())

string="python at guru99"
print(string.upper())
print(string.capitalize())

string="PYTHON AT GURU99"
print(string.lower())
print(":".join("Python"))

string="12345"
print(''.join(reversed(string)))
print(dir(string))
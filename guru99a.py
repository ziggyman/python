def doSomething():
    print('I am a function')

def functionWithOneParameter(a):
    print('a = <',a,'>')

def functionWithTwoParameters(a,b):
    print('a + b = <',a+b,'>')

def functionWithTwoParametersAndOneOptionalParameter(a, b, c=0):
    print('a + b + c = <',a + b + c,'>\\n')

def functionWithTwoParametersAndMultipleOptionalParameters(a, b, c=0, d = 1, e = 2):
    print('a + b + c + d + e = <',a + b + c + d + e,'>')

def main():
    print("Hello World!")

if __name__== "__main__":
    main()

print("Guru99")

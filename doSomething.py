from guru99a import doSomething, functionWithOneParameter, functionWithTwoParameters
from guru99a import functionWithTwoParametersAndOneOptionalParameter,functionWithTwoParametersAndMultipleOptionalParameters

a = 1
print('a = ',a)
b = 2
print('b = ',b)

truthValue = (a == b)
print('truthValue = ',truthValue)

c = 1
truthValue = (a == c)
print('truthValue = ',truthValue)

truthValue = (a < b)
print('truthValue = ',truthValue)

truthValue = (a >= b)
print('truthValue = ',truthValue)

truthValue = (a <= c)
print('truthValue = ',truthValue)

doSomething()

functionWithOneParameter('abc')
functionWithOneParameter(7)

functionWithTwoParameters('x', 'yz')
functionWithTwoParameters(1, 2)
functionWithTwoParameters(1.3, 2.7)

functionWithTwoParametersAndOneOptionalParameter(1.3, 2.7)
functionWithTwoParametersAndOneOptionalParameter(1.3, 2.7, 1.0)
functionWithTwoParametersAndOneOptionalParameter(b=2.7, a=1.3, c=1.0)

functionWithTwoParametersAndMultipleOptionalParameters(1, 2)
functionWithTwoParametersAndMultipleOptionalParameters(1, 2, e=7)

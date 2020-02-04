with open("guru99.txt","w+") as f:
    for i in range(10):
        f.write("This is line %d\r\n" % (i+1))

with open('guru99.txt','r') as f:
    if f.mode == 'r':
        contents = f.read()
print(contents)


with open('guru99.txt','r') as f:
    if f.mode == 'r':
        fl = f.readlines()
        for x in fl:
            print(x)
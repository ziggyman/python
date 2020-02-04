def main():
    x,y = 10,8

    st = 'I am a string'
    if x < y:
        st = 'x is smaller than y'
    elif x == y:
        st = 'x is equal to y'
    else:
        st = 'x is greater than y'
    print(st)

    st = ("x is less than y" if (x < y) else "x is greater than or equal to y")
    print(st)

    total = 100
    country = "US"
    #country = "AU"
    if country == "US":
        if total <= 50:
            print("Shipping Cost is  $50")
    elif total <= 100:
        print("Shipping Cost is $25")
    elif total <= 150:
        print("Shipping Costs $5")
    else:
        print("FREE")


    if country == "AU":
        if total <= 50:
            print("Shipping Cost is  $100")
    else:
        print("FREE")

def SwitchExample(argument):
    switcher = {
        0: " This is Case Zero ",
        1: " This is Case One ",
        2: " This is Case Two ",
    }
    return switcher.get(argument, "nothing")

if __name__ == '__main__':
#    main()
    argument = 3
    print(SwitchExample(argument))

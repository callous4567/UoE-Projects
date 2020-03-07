def list_squarer(List):
    new_list = [d**2 for d in List]
    new_list.sort()
    print(new_list)

def list_getter():
    List = []
    for i in range(0,10):
        variable = float(input(("Please give your variable {0}:").format(i)))
        List.append(variable)
    print(List)
    return(List)

FirstList = list_getter()
list_squarer(FirstList)

print("Hello c:")
print("Provide coefficients for equation, please. First the x^2, then x, then the numerical addition.")
x_sq = float(input())
x = float(input())
num = float(input())
# b^2 - 4ac = the discri
# -b +- (disc)^0.5 divide by 2 a
discri = (x**2 - 4*(x_sq*num))
sol1 = (-1*x - discri)/(2*x_sq)
sol2 = (-1*x + discri)/(2*x_sq)

print(discri, sol1, sol2)

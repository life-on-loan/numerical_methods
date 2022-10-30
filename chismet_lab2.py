import math


def f(x):
    return math.log(x) + pow((x - 2), 3)


epsilon = 0.5 * pow(10, -5)

''' метод половинного деления (дихотомии) '''


def check_sign(a, b):
    return a * b > 0


def dichotomy(func, a, b):
    c = b
    c_next = a
    count = 0
    while abs(c - c_next) > epsilon:

        c_next = c
        c = (a + b) / 2.0
        if check_sign(func(a), func(c)):
            a = c
        else:
            b = c
        count += 1
    return c, count


x, count = dichotomy(f, 1, 2)

print('дихотомия:  ' + str(x) + '  n: ' + str(count))

'''__________________________________________________________________________________'''
''' метод Ньютона '''


def derivative(x):
    return 1 / x + 3 * x ** 2 - 12 * x + 12


def newtons_method(f, a, b, derivative):
    x0 = a
    xn = x0
    xn1 = xn - f(xn) / derivative(xn)
    count = 1
    while abs(xn1 - xn) > epsilon:
        xn = xn1
        xn1 = xn - f(xn) / derivative(xn)
        count += 1
    return xn1, count


x_1, count = newtons_method(f, 1, 2, derivative)

print('Ньютон:     ' + str(x_1) + '   n: ' + str(count))

'''__________________________________________________________________________________'''
''' модифицированный метод Ньютона'''


def newtons_method_modified(f, a, b, derivative):
    x0 = a
    xn = x0
    derivative_x0 = derivative(x0)
    xn1 = xn - f(xn) / derivative_x0
    count = 1
    while abs(xn1 - xn) > epsilon:
        xn = xn1
        xn1 = xn - f(xn) / derivative_x0
        count += 1
    return xn1, count


x_2, count = newtons_method_modified(f, 1, 2, derivative)

print('Ньютон мод: ' + str(x_2) + '  n: ' + str(count))

'''__________________________________________________________________________________'''
''' метод подвижных хорд '''


def chords_active_method(f, a, b):
    x0 = a
    xn = b
    xn1 = xn - (f(xn) * (xn - x0)) / (f(xn) - f(x0))
    count = 1
    while abs(xn1 - xn) > epsilon:
        xn_1 = xn
        xn = xn1
        xn1 = xn - (f(xn) * (xn - xn_1)) / (f(xn) - f(xn_1))
        count += 1
    return xn1, count


x_3, count = chords_active_method(f, 1, 2)

print('хорд под:   ' + str(x_3) + '  n: ' + str(count))

'''__________________________________________________________________________________'''
''' метод неподвижных хорд '''


def chords_no_active_method(f, a, b):
    x0 = a
    xn = b
    xn1 = xn - (f(xn) * (xn - x0)) / (f(xn) - f(x0))
    count = 1
    while abs(xn1 - xn) > epsilon:
        xn = xn1
        xn1 = xn - (f(xn) * (xn - x0)) / (f(xn) - f(x0))
        count += 1
    return xn1, count


x_4, count = chords_no_active_method(f, 1, 2)

print('хорд непод: ' + str(x_4) + '    n: ' + str(count))

'''__________________________________________________________________________________'''
''' метод простой итерации '''


def iteration_method(fi, a, b):
    x0 = a
    xn = fi(x0)
    xn1 = fi(xn)
    count = 2
    while abs(xn1 - xn) > epsilon:
        xn = xn1
        xn1 = fi(xn)
        count += 1
    return xn1, count


def fi(x):
    return (-1) * (math.log(x) + pow(x, 3) - 6 * x * x - 8) / 12


x_5, count = iteration_method(fi, 1, 2)

print('итерация:   ' + str(x_5) + '  n: ' + str(count))

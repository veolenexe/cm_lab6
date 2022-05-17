import math
import matplotlib.pyplot as plt

e = math.e
epsilon = 5e-05
y_1 = (e + (1 / e) - 2)


def f(x, y):
    return y + q(x)


def phi(x):
    return x - y_1


def q(x):
    M = 16
    alpha = 2 + 0.1 * M
    return 2 * alpha + 2 + alpha * x * (1 - x)


def getPoints(N, a=0, b=1):
    h = (b - a) / N
    points = []
    cur_point = a
    for i in range(N + 1):
        points.append(cur_point)
        cur_point += h
    return points, h


# Точное решение задачи
def get_exact_solution(x):
    return (3.6 * x * x) - (3.6 * x) + (e ** (-x)) + (
            e ** x) - 2  # y(x) = 3.6 x^2 - 3.6 x + e^(-x) + e^x - 2.


def rungeIteration(x, y, h, cur_val):
    K1 = h * f(x, y)
    K2 = h * f(x + h / 2, y + K1 / 2)
    K3 = h * f(x + h / 2, y + K2 / 2)
    K4 = h * f(x + h, y + K3)
    return cur_val + (1 / 6) * (K1 + 2 * K2 + 2 * K3 + K4)


def rungeIteration_for_y(x, cur_z, h, cur_y):
    K1 = h * cur_z
    K2 = h * rungeIteration(x + h / 2, cur_y + K1 / 2, h, cur_z)
    K3 = h * rungeIteration(x + h / 2, cur_y + K2 / 2, h, cur_z)
    K4 = h * rungeIteration(x + h, cur_y + K3, h, cur_z)
    return cur_y + (1 / 6) * (K1 + 2 * K2 + 2 * K3 + K4)


# Метод Рунге-Кутта четвертого порядка
def get_runge_kutta_solution(points, h, n):
    result_z = [n]
    cur_z = n
    result_y = [0]
    cur_y = 0

    for x in points[1:]:
        cur_y = rungeIteration_for_y(x, cur_z, h, cur_y)
        cur_z = rungeIteration(x, cur_y, h, cur_z)

        result_z.append(cur_z)
        result_y.append(cur_y)

    return result_z, result_y


# метод Адамса трехшаговый
def get_adams_solution(points, h, n):
    cur_x_list = points[:3]
    # получаем первые 3 значения методом Рунге-Кутта четвертого порядка
    result_z, result_y = get_runge_kutta_solution(cur_x_list, h, n)

    for x in points[3:]:
        cur_y = (result_y[-1] + (h / 12 *
                                 (23 * result_z[-1]
                                  - 16 * result_z[-2]
                                  + 5 * result_z[-3])))

        cur_z = (result_z[-1] + (h / 12 *
                                 (23 * f(cur_x_list[-1], result_y[-1])
                                  - 16 * f(cur_x_list[-2], result_y[-2])
                                  + 5 * f(cur_x_list[-3], result_y[-3]))))

        cur_x_list.append(x)
        result_y.append(cur_y)
        result_z.append(cur_z)

    return result_y


def main():
    get_result(10)
    get_result(20)
    get_result(1500)


def is_diagonal_dominance(h):
    return abs(2 + h ** 2) > 2


# прогонка
def get_sweep_solution(points, h, N):
    if not is_diagonal_dominance(h):
        return

    result = [0.0] + [0.0 for _ in range(1, N)] + [y_1]
    result_lambda = [0] + [0 for _ in range(N)]
    result_n = [0.0] + [0.0 for _ in range(1, N)] + [y_1]

    # Прямая прогонка
    for i in range(1, N + 1):
        cur_lambda = result_lambda[i - 1]
        cur_n = result_n[i - 1]
        A = 2 + h ** 2
        B = (h ** 2) * q(points[i - 1])
        result_lambda[i] = 1.0 / (A - cur_lambda)
        result_n[i] = (cur_n - B) / (A - cur_lambda)

    # Обратная прогонка
    for i in range(1, N):
        result[-i - 1] = result_lambda[-i] * result[-i] + result_n[-i]

    return result


def get_shoot_solution(points, h):
    n1 = -10  # наше мю1
    n2 = 10  # наше мю2
    sol1 = get_adams_solution(points, h, n1)
    sol2 = get_adams_solution(points, h, n2)
    print(sol1)
    print(sol2)

    phi_n1 = phi(sol1[-1])
    phi_n2 = phi(sol2[-1])

    print("phi_n1 = ", phi_n1)
    print("phi_n2 = ", phi_n2)
    if phi_n1 * phi_n2 >= 0:
        print('все плохо')
        return
    print("фи разных знаков, значит можно использовать метод дихотомии")

    def dichotomies_solve(a, b, phi_l, phi_r):
        l = a
        r = b
        m = (l + r) / 2
        iteration = math.ceil(math.log2((b - a) / epsilon))
        for _ in range(iteration):
            m = (l + r) / 2
            phi_m = phi(get_adams_solution(points, h, m)[-1])
            if phi_l * phi_m < 0:
                r = m
                phi_r = phi_m
            elif phi_r * phi_m < 0:
                l = m
                phi_l = phi_m
            else:
                return m

        return m

    n = dichotomies_solve(n1, n2, phi_n1, phi_n2)  # нашли мю
    print('n= ', n)
    return get_adams_solution(points, h, n)


def get_result(N):
    print()
    points, h = getPoints(N)
    print("N = " + str(N))
    print("h = ", h)
    print("точки", points)

    shoot_solution = get_shoot_solution(points, h)
    print("стрельба :", shoot_solution)

    sweep_sol = get_sweep_solution(points, h, N)
    print("прогонка :", sweep_sol)

    exact_solution = list(map(lambda x: get_exact_solution(x), points))
    print("Точное решение", exact_solution)

    plt.plot(points, shoot_solution)
    plt.plot(points, sweep_sol)
    plt.plot(points, exact_solution)
    plt.legend(["Стрельба", "Прогонка", "Точное решение"])
    plt.show()


if __name__ == '__main__':
    main()

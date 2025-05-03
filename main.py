# Добавление зависимостей для работы с FastAPI
from fastapi import Body, FastAPI, HTTPException
from fastapi.responses import FileResponse, JSONResponse
# Добавление зависимости для выполнения вычислений, связанных с диффурами
from scipy. integrate import odeint
import random
from copy import deepcopy
# Инициализация приложения FastAPI
app = FastAPI()

# Вспомогательный класс F, выражающий полином третьей степени.
# Поля:
# a - коэффициент при x^3
# b - коэффициент при x^2
# c - коэффициент при x
# d - свободный коэффициент
# Методы:
# calc(x) - вычисляет значение полинома при заданном x
class F:
    def __init__(self, a, b, c, d):
        self.a = a
        self.b = b
        self.c = c
        self.d = d

    def calc(self, x):
        return self.a* (x**3) + self.b * (x**2) + self.c * x + self.d

# Функция, выражающая систему уравнений Практической Работы 1
# Принимает:
# u - массив исследуемых параметров, представляющих собой функции xi(t), i=[1..31]
# t - массив временных точек от 0 до 1
# с - словарь констант вида {"Название": <Значение>, ...}
# f - словарь полиномов вида {"Название": <Экземпляр полинома>, ...}
# Возвращает:
# массив вычисленных выражений dxi/dt, i=[1..31]
def du_dt(u,t,c,f):

        # Извлекаем элементы, которые являются исследуемыми xi(t), из массива u
        [x1_t, x2_t, x3_t, x4_t, x5_t, x6_t, x7_t, x8_t, x9_t, x10_t, x11_t, x12_t, x13_t, x14_t, x15_t, x16_t, x17_t,
        x18_t, x19_t, x20_t, x21_t, x22_t, x23_t, x24_t, x25_t, x26_t, x27_t, x28_t, x29_t, x30_t, x31_t] = u

        # Для каждого уравнения приводим запись формулы из задания.
        # Для формулы:
        # dx1/dt = x1(t) * ( BPn / x1(t) * f1(x3(t)) * f2(x19(t)) * f6(x27(t)) + BPn / x1(t) * f3(x21(t)) * f4(x22(t)) * f5(x23(t)) * f7(x29(t)))
        # программная запись примет вид:
        dx1_dt = x1_t * (
            - c['BPn'] / x1_t *
            f['F1'].calc(x3_t) *
            f['F2'].calc(x19_t) *
            f['F6'].calc(x27_t)
            + c['BPk'] / x1_t *
            f['F3'].calc(x21_t) *
            f['F4'].calc(x22_t) *
            f['F5'].calc(x23_t) *
            f['F7'].calc(x29_t)
        )

        #Далее по аналогии со всеми уравнениями системы...
        dx2_dt = x2_t * (
            - c['MPn'] / x2_t *
            f['F8'].calc(x4_t) *
            f['F9'].calc(x20_t) *
            f['F12'].calc(x28_t)
            + c['MPk'] / x2_t *
            f['F10'].calc(x22_t) *
            f['F11'].calc(x23_t) *
            f['F13'].calc(x29_t)
        )

        dx3_dt = x3_t * (
            c['BV'] / c['BZ'] *
            f['F14'].calc(x19_t)
            - c['BC'] / c['BZ']
        )

        dx4_dt = x4_t * (
            c['MV'] / c['MZ'] *
            f['F15'].calc(x20_t)
            - c['MC'] / c['MZ']
        )

        dx5_dt = x5_t * (c['C'] / c['P'])

        dx6_dt = x6_t * (
            c['PT'] *
            f['F21'].calc(x12_t) + c['PI'] *
            f['F17'].calc(x8_t) + c['PJ'] *
            f['F19'].calc(x10_t) + c['PD'] *
            f['F22'].calc(x17_t) + c['PS']
        ) * (
            f['F20'].calc(x11_t) *
            f['F23'].calc(x25_t) *
            f['F24'].calc(x30_t)
        )

        dx7_dt = x7_t * (
            (c['G'] + c['IG'] + c['NG']) / c['F']
        ) * f['F25'].calc(x14_t)

        dx8_dt = x8_t * (-c['An'] + c['Ak']) / c['P']

        dx9_dt = x9_t * ((-c['Vn'] + c['Vk']) / c['P'] *
            f['F16'].calc(x5_t) *
            f['F18'].calc(x6_t))

        dx10_dt = x10_t * (-c['IDPn'] + c['IDPk']) / c['P']

        dx11_dt = x11_t * ((-c['DPn']) / x11_t + c['DPk'] / x11_t )

        dx12_dt = x12_t * (-c['IPn'] + c['IPk']) / c['P']

        dx13_dt = x13_t * (c['IB'] + c['IN'] + c['IL']) / c['F']

        dx14_dt = x14_t * (c['DS'] *
            f['F28'].calc(x15_t) + c['DV'] *
            f['F29'].calc(x16_t) + c['DP'])

        dx15_dt = x15_t * (-c['JDn'] + c['JDk']) / c['P']

        dx16_dt = x16_t * (-c['INn'] + c['INk']) / c['NPP']

        dx17_dt = x17_t * (-c['PPn'] + c['PPk']) / c['P']

        dx18_dt = x18_t * (c['Dn'] / c['D'] *
           f['F31'].calc(x26_t))

        dx19_dt = x19_t * (-c['PBn'] + c['PBk']) / c['PB']

        dx20_dt = x20_t * (-c['PMn'] + c['PMk']) / c['P']

        dx21_dt = x21_t * ((-c['IPn']) / x21_t + c['IPk'] / x21_t )

        dx22_dt = x22_t * ((-c['OPn']) / x22_t + c['OPk'] / x22_t )

        dx23_dt = x23_t * ((-c['JPn']) / x23_t + c['JPk'] / x23_t )

        dx24_dt = x24_t * ((-c['ISPn']) / x24_t + c['ISPk'] / x24_t *
            f['F33'].calc(x26_t) )

        dx25_dt = x25_t * (-c['MIPn'] + c['MIPk']) / c['P']

        dx26_dt = x26_t * ((-c['IIDPn']) / x26_t * f['F32'].calc(x21_t) + c['IIDPk'] / x26_t)

        dx27_dt = x27_t * ((-c['BRPn']) / x27_t + c['BRPk'] / x27_t )

        dx28_dt = x28_t * ((-c['MRPn']) / x28_t + c['MRPk'] / x28_t )

        dx29_dt = x29_t * ((-c['SRPn']) / x29_t + c['SRPk'] / x29_t )

        dx30_dt = x30_t * (-c['RPn'] + c['RPk']) / c['P']

        dx31_dt = x31_t * (c['IR'] + c['NR']) / c['F'] * f['F34'].calc(x30_t)

        # Возвращаем массив вычисленных dxi/dt, i=[1..31]
        return [dx1_dt, dx2_dt, dx3_dt, dx4_dt, dx5_dt, dx6_dt, dx7_dt, dx8_dt, dx9_dt, dx10_dt, dx11_dt, dx12_dt, dx13_dt,
        dx14_dt, dx15_dt, dx16_dt, dx17_dt, dx18_dt, dx19_dt, dx20_dt, dx21_dt, dx22_dt, dx23_dt, dx24_dt, dx25_dt, dx26_dt,
        dx27_dt, dx28_dt, dx29_dt, dx30_dt, dx31_dt]

# Указывает маршрут / до web-интерфейса, созданного для ввода значений из Практической Работы 1
@app.get("/")
async def main():
    # Возвращаем в качестве ответа HTML-страницу, расположенную по адресу /home/AlexAranara/my_fastapi/web/index.html
    return FileResponse("index.html")


# FUNCS остаётся фиксированным набором индексов полиномов
FUNCS = [
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,
    17, 18, 19, 20, 21, 22, 23, 24, 25, 28, 29, 31, 32, 33, 34
]


def _random_coefficients() -> dict[int, list[float]]:
    """Генерирует случайные коэффициенты (от 0 до 1) для всех полиномов."""
    # Сейчас установлен диапазон для всех значений [0; 0.5]
    return {k: [round(random.random() / 2, 3) for _ in range(4)] for k in FUNCS} 

def _error_score(solution: list[list[float]], thresholds: list[float]) -> int:
    """
    Считает число графиков Xi, которые
    хотя бы в одной точке ушли ниже своей нижней границы.
    """
    score = 0
    for i in range(31):                     # X1 … X31
        if any(row[i] < thresholds[i] for row in solution):
            score += 1
    return score

counter = 0
def has_outliers(solution: list[list[float]]):
    for i in range(31):                     
        if any(row[i] > 15 for row in solution): # Здесь можно контролировать вылет
            global counter
            counter += 1
            return True
    return False

def find_solution(t0: dict[str, float],
                  c: list[float],
                  uf,
                  x_threshold: list[float],
                  iters_number: int = 100):
    """
    Подбирает коэффициенты полиномов так, чтобы значения X1–X31
    нигде не опускались ниже своих порогов.
    Возвращает решение системы и найденные коэффициенты.
    """
    t_span = [i * 0.05 for i in range(21)]          # 0 … 1 с шагом 0.05

    best_score = 31                                 # максимальное возможное
    best_solution = None
    best_coeffs = None

    for _ in range(iters_number):
        coeffs = _random_coefficients()
        f = {k: F(*coeffs[int(k[1:])]) for k in uf}       # создаём экземпляры F

        try:
            usolution, d = odeint(
                du_dt,
                list(t0.values()),
                t_span,
                args=(c, f,),
                full_output=True,
                # rtol=1e-3, # Здесь можно поменять точность поиска решения
                # atol=1e-3  # Здесь можно поменять точность поиска решения
            )
        except Exception:
            continue

        # odeint отдаёт numpy‑массив; переводим в чистые list
        solution = [list(row) for row in usolution.tolist()]


        if has_outliers(solution):
            continue

        score = _error_score(solution, x_threshold)

        if score == 0:
            return solution, coeffs, score                # найдено идеальное решение

        if score < best_score:
            best_score = score
            best_solution = deepcopy(solution)
            best_coeffs = deepcopy(coeffs)

    # Если идеального решения нет, возвращаем лучшее найденное
    return best_solution, best_coeffs, best_score


# Указываем маршрут /api/count, по которому можно получить результаты вычислений для Практической работы 1
# Тело данного POST-запрос имеет вид {"x0": <Словарь начальный значений>, "c": <Словарь констант>, "f": <Словарь полиномов>}
# Записи <Словарь начальный значений> имеют вид "Имя значения": <значение>
# Записи <Словарь констант> имеют вид "Имя константы": <значение>
# Записи <Словарь полиномов> имеют вид "Имя полинома": {"a": <значение>, "b": <значение>, "c": <значение>, "d": <значение> }
@app.post("/api/count")
async def count(data  = Body()):
    # Извлекаем пришедшие начальные значения в t0
    t0 = data['x0']

    # Извлекаем пришедшие константы в с
    c = data['c']

    # Извлекаем пришедшие коэффициенты полиномов в uf
    uf = data['f']
    
    x_threashold = data['norm']

    usolution, coefficients, score  = find_solution(t0, c, uf, x_threashold, iters_number=30) # Менять кол-во итераций ЗДЕСЬ!
    print('Количество функций, которые вышли за нижнюю границу: ', score)
    print('Количество комбинаций во время перебора, при которых где-то случился вылет за верхнюю границу: ', counter)


    if usolution is None:
        print("Найти решение для таких входных данных не удалось")
        return JSONResponse(status_code=400, content={"error": "Найти решение для таких входных данных не удалось"})
    
    # Преобразуем значения для отправки их в качестве JSON-ответа
    solution = [None] * len(usolution)
    idx = 0
    for elem in usolution:
        solution[idx] = list(elem)
        idx = idx + 1
    
    return JSONResponse(content={"message": solution, "coefficients": coefficients})

# Надо подобрать такие коэфы полиномов, чтобы значения функций, которые мы строим на графике, были больше установленного порога
# Порог надо будет выставлять отдельным полем, которое можно отобразить на главном экране. 
# На фронт надо будет отдать значения коэфов
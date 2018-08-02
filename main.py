import random
import math
import copy

ARRAY_SIZE = 2
MUTATE_RATE = 0.8
TARGET_NUMBER = 4


class AIRS:

    def __init__(self, cells_number):
        self.cells_number = cells_number
        self.target_number = TARGET_NUMBER
        self.pools = {target: Pool(target, cells_number) for target in range(TARGET_NUMBER)}

    @staticmethod
    def distance(array1, array2):
        _MAX_DISTANCE = 2
        d = 0
        for i, j in zip(array1, array2):
            d += (i - j) ** 2
        return math.sqrt(d) / _MAX_DISTANCE

    def stimulate(self, pattern):

        for index in range(len(self.pools)):
            self.pools[index].stimulate(pattern=pattern)

    def get_most_stimulated_cell(self, pattern):
        self.stimulate(pattern=pattern)
        best_cell = self.pools[0].cells[0]
        for index in range(len(self.pools)):
            _p = self.pools[index].cells
            for _c in _p:
                # print(_c.target)
                if _c.affinity > best_cell.affinity:
                    best_cell = _c
        return best_cell


class Pool:
    def __init__(self, target, cells_number):
        if target is None:
            raise ValueError("Target must not be None.")
        else:
            self.target = target
            self.cells = [Cell(target=target) for _ in range(cells_number)]

    def stimulate(self, pattern):
        for cell in self.cells:
            cell.stimulate(pattern=pattern)


class Cell:
    def __init__(self, array=None, target=None):
        if array is None:
            self.array = [random.random() for _ in range(ARRAY_SIZE)]
        else:
            self.array = array
        self.target = target
        self.affinity = float('inf')
        self.stimulation = float('inf')

    def stimulate(self, pattern):
        self.affinity = 1 - AIRS.distance(pattern, self.array)

    def __str__(self):
        return str(self.array)

    def __repr__(self):
        return str(self.array)

    # todo : maybe (have to) change this
    def mutate(self):
        for idx in range(len(self.array)):
            if random.random() < self.affinity * MUTATE_RATE:
                self.array[idx] = random.random() * (max(self.array) - min(self.array)) + min(self.array)


def train_AIRS():
    airs = AIRS(cells_number=10)

    for g in range(500):
        vector = [random.random() for _ in range(2)]
        target = -1

        if vector[0] < 0.5 and vector[1] < 0.5:
            target = 0
        elif vector[0] >= 0.5 > vector[1]:
            target = 1
        elif vector[0] < 0.5 <= vector[1]:
            target = 2
        else:
            target = 3

        most = airs.get_most_stimulated_cell(vector)
        print("Affinity : " + str(most.affinity))
        print("Best : " + str(most))
        print("Best target : " + str(most.target))
        print("Vector : " + str(vector))
        print("Real target : " + str(target))
        if most.target != target:
            print("No")
            airs.pools[target].cells.append(Cell(array=vector, target=target))
        else:
            print("Yes")
            if most.affinity < 1.0:
                for c in airs.pools[target].cells:
                    c.mutate()
        print("\n")


if __name__ == '__main__':
    train_AIRS()

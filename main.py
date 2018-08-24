import math
import random

ARRAY_SIZE = 4  # Features number
MAX_ITER = 5
mapping = {"Iris-setosa": 0, "Iris-versicolor": 1, "Iris-virginica": 2}
reverseMapping = {0: "Iris-setosa", 1: "Iris-versicolor", 2: "Iris-virginica"}


class AIRS:
    def __init__(self, hyper_clonal_rate, clonal_rate, targets_number, mc_init_rate,
                 total_num_resources, affinity_threshold_scalar, k, test_size):
        self.HYPER_CLONAL_RATE = hyper_clonal_rate
        self.CLONAL_RATE = clonal_rate
        self.AFFINITY_THRESHOLD = 0
        self.TARGETS_NUMBER = targets_number
        self.MC_INIT_RATE = mc_init_rate
        self.TOTAL_NUM_RESOURCES = total_num_resources
        self.AFFINITY_THRESHOLD_SCALAR = affinity_threshold_scalar
        self.TEST_SIZE = test_size
        self.K = k
        self.MC = None
        self.AB = None

    @staticmethod
    def affinity(vector1, vector2):
        _MAX_DISTANCE = 10.91
        d = 0
        for i, j in zip(vector1, vector2):
            d += (i - j) ** 2
        return math.sqrt(d) / _MAX_DISTANCE

    def train_test_split(self):
        with open("iris.data", "r") as data:
            content = data.readlines()
            ret = [([float(x.split(",")[i]) for i in range(4)], mapping[x.split(",")[4][:-1]]) for x in content]
            min1 = 100
            min2 = 100
            min3 = 100
            max1 = 0
            max2 = 0
            max3 = 0
            for x in ret:
                if x[0][0] < min1:
                    min1 = x[0][0]
                if x[0][0] > max1:
                    max1 = x[0][0]

                if x[0][1] < min2:
                    min2 = x[0][1]
                if x[0][1] > max2:
                    max2 = x[0][1]

                if x[0][2] < min3:
                    min3 = x[0][2]
                if x[0][0] > max3:
                    max3 = x[0][2]
            random.shuffle(ret)

            print(str(min1) + " " + str(min2) + " " + str(min3) + " " + str(max1) + " " + str(max2) + " " + str(max3))
        # return [([1, 1], 1), ([0.1, 0.1], 0), ([0.9, 0.95], 1), ([0, 0], 0), ([0.11, 0.11], 0)]
        return ret[:int((1 - self.TEST_SIZE) * len(ret))], ret[int((1 - self.TEST_SIZE) * len(ret)):]

    def calculate_affinity_threshold(self):
        train_set, test_set = self.train_test_split()
        affinity_threshold = 0
        for i in range(len(train_set)):
            for j in range(i + 1, len(train_set)):
                affinity_threshold += self.affinity(train_set[i][0], train_set[j][0])

        self.AFFINITY_THRESHOLD = affinity_threshold / (len(train_set) * (len(train_set) - 1) / 2)

    def init_MC(self, train_set, MC):
        for _ in range(int(len(train_set) * self.MC_INIT_RATE)):
            seed_cell = random.choice(train_set)
            print(seed_cell)
            print(seed_cell[1])
            print(seed_cell[0])

            MC[seed_cell[1]].append(Cell(vector=seed_cell[0], target=seed_cell[1]))

    def argminARB(self, AB, target):
        minRes = 1.0
        ab = None
        abIndex = None
        for i in range(len(AB[target])):
            if AB[target][i].resources <= minRes:
                minRes = AB[target][i].resources
                ab = AB[target][i]
                abIndex = i

        return ab, abIndex

    def getMcCandidate(self, AB, target):
        maxStim = 0.0
        ab = None
        for i in range(len(AB[target])):
            if AB[target][i].stimulation >= maxStim:
                maxStim = AB[target][i].stimulation
                ab = AB[target][i]
        c = Cell(vector=ab.vector, target=ab.target)
        c.stimulation = ab.stimulation
        return c

    def train(self):
        train_set, test_set = self.train_test_split()

        self.calculate_affinity_threshold()
        print("Aff thresh " + str(self.AFFINITY_THRESHOLD))
        MC = {target: [] for target in range(self.TARGETS_NUMBER)}
        AB = {target: [] for target in range(self.TARGETS_NUMBER)}

        # MC Initialisation
        self.init_MC(train_set, MC)

        for antigene, _class in train_set:
            print("\n\n-------------------------------------- ")
            print("MC : " + str(MC))
            print("Antigene = " + str(antigene))
            print("Class = " + str(_class))
            # MC Identification
            mc_match = None
            if len(MC[_class]) == 0:
                mc_match = Cell(vector=antigene, target=_class)
                MC[_class].append(mc_match)
            else:
                best_stim = 0
                for c in MC[_class]:
                    if c.stimulation >= best_stim:
                        best_stim = c.stimulation
                        mc_match = c

            # ARB Generation
            AB[_class].append(ARB(vector=mc_match.vector, target=mc_match.target))  # add the mc_match to ARBs
            print("Mc_ match " + str(mc_match))
            stim = mc_match.stimulate(antigene)

            iterations = 0
            while True:
                iterations += 1
                print("Stim = " + str(stim))
                MAX_CLONES = int(self.HYPER_CLONAL_RATE * self.CLONAL_RATE * stim)
                print("Max Clones " + str(MAX_CLONES))
                num_clones = 0
                while num_clones < MAX_CLONES:
                    clone, mutated = mc_match.mutate()

                    if mutated:
                        AB[_class].append(clone)
                        num_clones += 1
                print("AB size : " + str(len(AB[0]) + len(AB[1])))
                print("AB : " + str(AB))

                # Competition for resources

                avgStim = sum([x.stimulate(antigene) for x in AB[_class]]) / len(AB[_class])

                print("Average Stim = " + str(avgStim))

                MIN_STIM = 1.0
                MAX_STIM = 0.0

                for c in AB.keys():
                    for ab in AB.get(c):
                        stim = ab.stimulate(antigene)
                        if stim < MIN_STIM:
                            MIN_STIM = stim
                        if stim > MAX_STIM:
                            MAX_STIM = stim
                print("Res ! ")
                for c in AB.keys():
                    for ab in AB.get(c):
                        ab.stimulation = (ab.stimulation - MIN_STIM) / (MAX_STIM - MIN_STIM)
                        ab.resources = ab.stimulation * self.CLONAL_RATE
                        print("heey " + str(ab))

                resAlloc = sum([x.resources for x in AB[_class]])
                print("Res alloc = " + str(resAlloc))
                numResAllowed = self.TOTAL_NUM_RESOURCES
                while resAlloc > numResAllowed:
                    numResRemove = resAlloc - numResAllowed
                    abRemove, abRemoveIndex = self.argminARB(AB=AB, target=_class)
                    print("Rem1 " + str(abRemove))
                    if abRemove.resources <= numResRemove:
                        AB[_class].remove(abRemove)
                        resAlloc -= abRemove.resources
                    else:
                        AB[_class][abRemoveIndex].resources -= numResRemove
                        resAlloc -= numResRemove

                if (avgStim > self.AFFINITY_THRESHOLD) or (iterations >= MAX_ITER):
                    print("MAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAX")
                    break

            print("AB[_class] = " + str(AB[_class]))
            mc_candidate = self.getMcCandidate(AB=AB, target=_class)
            print("Candidate = " + str(mc_candidate))
            print("Mc match == " + str(mc_match))

            if mc_candidate.stimulation > mc_match.stimulation:
                if AIRS.affinity(mc_candidate.vector,
                                 mc_match.vector) < self.AFFINITY_THRESHOLD * self.AFFINITY_THRESHOLD_SCALAR:
                    # The mc candidate replaces the mc match
                    MC[_class].remove(mc_match)
                # Add the mc_match to MC pool
                MC[_class].append(mc_candidate)
            print("MC == " + str(MC))
        print("MC SIZE = " + str(len(MC[0]) + len(MC[1]) + len(MC[2])))
        self.MC = MC
        self.AB = AB

        n_correct = 0
        for ag, _class in test_set:
            if self.classify(ag) == reverseMapping[_class]:
                n_correct += 1

        print("Accuracy : {}".format(n_correct / len(test_set)))

    def classify(self, antigene):
        if (self.MC is None) or (self.AB is None):
            raise Exception("AIRS must be trained first")

        vote_array = []
        for c in self.MC.keys():
            for ab in self.MC.get(c):
                ab.stimulate(antigene)
                vote_array.append(ab)

        vote_array = list(sorted(vote_array, key=lambda cell: -cell.stimulation))
        v = {0: 0, 1: 0, 2: 0}
        self.K = min(self.K, len(vote_array))
        for x in vote_array[:self.K]:
            v[x.target] += 1

        maxVote = 0
        target = 0
        for x in v.keys():
            if v[x] > maxVote:
                maxVote = v[x]
                target = x
        return reverseMapping[target]


class ARB:
    def __init__(self, vector=None, target=None):
        if vector is None:
            self.vector = [random.random() for _ in range(ARRAY_SIZE)]
        else:
            self.vector = vector
        self.target = target
        self.stimulation = float('inf')
        self.resources = 0

    def __str__(self):
        return "Vector = " + str(self.vector) + " | target = " + str(self.target) + " | stim = " + str(
            self.stimulation) + " | res = " + str(self.resources)

    def __repr__(self):
        return "Vector = " + str(self.vector) + " | target = " + str(self.target) + " | stim = " + str(
            self.stimulation) + " | res = " + str(self.resources)

    def stimulate(self, pattern):
        self.stimulation = 1 - AIRS.affinity(vector1=pattern, vector2=self.vector)
        return self.stimulation

    def mutate(self):
        _range = 1 - self.stimulation
        mutated = False
        new_vector = []

        for v in self.vector:
            print(v)
            change = random.random()
            change_to = 7 * random.random() + 0.1

            if change <= AIRS.MUTATION_RATE:
                new_vector.append(change_to)
                mutated = True
            else:
                new_vector.append(v)

        return ARB(vector=new_vector, target=self.target), mutated


class Cell:
    def __init__(self, vector=None, target=None):
        if vector is None:
            self.vector = [random.random() for _ in range(ARRAY_SIZE)]
        else:
            self.vector = vector
        self.target = target
        self.stimulation = float('inf')

    def __str__(self):
        return "Cell : Vector = " + str(self.vector) + " | target = " + str(self.target) + " | stim = " + str(
            self.stimulation)

    def __repr__(self):
        return "Cell : Vector = " + str(self.vector) + " | target = " + str(self.target) + " | stim = " + str(
            self.stimulation)

    def stimulate(self, pattern):
        self.stimulation = 1 - AIRS.affinity(vector1=pattern, vector2=self.vector)
        return self.stimulation

    def mutate(self):
        _range = 1 - self.stimulation
        mutated = False
        new_vector = []

        for v in self.vector:
            change = random.random()
            change_to = random.random()

            if change <= MUTATION_RATE:
                new_vector.append(change_to)
                mutated = True
            else:
                new_vector.append(v)

        return ARB(vector=new_vector, target=self.target), mutated


if __name__ == '__main__':
    # 6.3,2.3,4.4,1.3,Iris-versicolor
    # 5.1, 2.5, 3.0, 1.1, Iris - versicolor
    # 6.9,3.1,5.1,2.3,Iris-virginica
    MUTATION_RATE = 0.2
    airs = AIRS(hyper_clonal_rate=20, clonal_rate=0.8, targets_number=3, mc_init_rate=0.5,
                total_num_resources=10, affinity_threshold_scalar=0.8, k=6, test_size=0.4)

    print(airs.train_test_split())
    airs.train()

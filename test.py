import configs_draft
import time
import itertools
import new_attempt



pot = ((3,0),(0,2))   # representation of the decomposed potential W = x^3 + y^2


for t in range(3,15):
    print("---------------")
    for tree in configs_draft.BinTreeGen(t):
        print(configs_draft.EnumerateConfigs(tree, pot))

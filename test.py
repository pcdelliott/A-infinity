import configs_draft
import time
import itertools



pot = ((4,0,0,0,0),(0,4,0,0,0),(0,0,4,0,0),(0,0,0,4,0),(0,0,0,0,4)). #Fermat quintic


for t in range(3,8):
    print('Trees with {} leaves'.format(t))
    for tree in configs_draft.BinTreeGen(t):
        if configs_draft.EnumerateConfigs(tree, pot) > 0:
            print((tree, configs_draft.EnumerateConfigs(tree, pot)))

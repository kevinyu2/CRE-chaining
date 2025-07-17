from chaining import *

'''
Testing chaining
'''

def generate_mems_like(n=10, a_range=(0, 10), b_range=(0, 15), seed=None):
    if seed is not None:
        np.random.seed(seed)

    a_vals = np.random.randint(*a_range, size=n)
    b_vals = np.random.randint(*b_range, size=n)
    mems = list(zip(a_vals, b_vals))

    return mems

def test_np():
    for i in range(100) :
        mems = generate_mems_like(n=100, a_range=(0, 100), b_range=(0, 100), seed=42)
        mems_np = np.array(mems, dtype = np.int32)

        assert(chain_driver(mems, False) == chain_driver_np(mems_np, False))
        assert(chain_local_driver(mems, 4, -2, -1, False) == chain_local_driver_np(mems_np, 4, -2, -1, False))
        print(i)

def local_no_weight():
    mems = [(1,6), (2,5), (3, 4), (4, 3)]
    mems_np = np.array(mems, dtype=np.int32)

    print(chain_local_driver(mems, 4, -2, -1, False))
    print(chain_local_driver_np(mems_np, 4, -2, -1, False))

    mems = [(0,0), (4,5), (2,8), (5,8), (1, 6), (6,12), (4,9)]
    mems_np = np.array(mems, dtype=np.int32)

    print(chain_local_driver(mems, 4, -2, -1, False))
    print(chain_local_driver_np(mems_np, 4, -2, -1, False))

def local_weight():
    mems = [(0,0,2), (4,5,1), (2,8,3), (5,8,2), (1, 6,1.5), (6,12,3.2), (4,9,1.2)]
    print(chain_local_weighted(mems, 2, -2, -1))
    #=12.8

    mems = [(0, 5, 2), (2, 1, 10), (3, 2, 10), (5, 7, 5)]
    print(chain_local_weighted(mems, 2, -2, -1))
    #=45

    mems = [(0, 10, 5), (1, 2, 1), (5, 1, 10), (7, 8, 1)]
    print(chain_local_weighted(mems, 2, -2, -1))
    #=20

    mems = [(3, 10, 5), (1, 2, 1), (5, 1, 10), (7, 8, 1), (0, 1, 1)]
    print(chain_local_weighted(mems, 2, -2, -1))
    #=20

    mems = [(1, 10, 10), (1, 2, 1), (5, 1, 10), (7, 8, 1), (0, 1, 6)]
    print(chain_local_weighted(mems, 2, -2, -1))
    #=24

    mems = [(0, 1, 6), (3, 10, 10)]
    print(chain_local_weighted(mems, 2, -2, -1))
    #=22

def global_weight():
    mems = [(1, 1, 4), (2, 1, 5), (3, 3, 1.5), (4, 10, 9), (5, 5, 6), (6, 6, 5)]
    print(chain_driver(mems, True))
    # == 17.5

    mems = [(3, 5, 10), (4, 7, 20), (5, 6, 1), (4, 5, 1), (1, 2, 1), (2, 3, 1), (3, 4, 1), (1, 1, 4)]
    print(chain_driver(mems, True))
    # == 35

    mems = [(4, 5, 1), (2, 3, 1), (3, 4, 1), (5, 1, 7)]
    print(chain_driver(mems, True))
    # == 8

def global_no_weight():
    mems = [(1, 7), (2, 6), (3, 5), (4, 4)]
    print(chain_driver(mems, False))
    # == 4

    mems = [(1,3), (4,4), (5,3), (6,2)]
    print(chain_driver(mems, False))
    # == 3


local_weight()
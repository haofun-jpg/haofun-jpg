#coding=utf-8
import multiprocessing as mp

def job(x):
    return x*x

def multicore():
    # 使用 Pool 自動分配給 CPU 的每個一核心 (core)
    pool = mp.Pool() 
    result = pool.map(job, range(10)) # 使用 pool 還可以接到 function 的回傳值 
    print(result) # [0, 1, 4, 9, 16, 25, 36, 49, 64, 81]

     # 使用 Pool 指定分配給 CPU 的 3 核心 (core)
    pool = mp.Pool(processes=3) 
    result = pool.map(job, range(10)) # 使用 pool 還可以接到 function 的回傳值 
    print(result) # [0, 1, 4, 9, 16, 25, 36, 49, 64, 81]

    # 除了 map 功能以外，還有另一個功能： apply_async
    # apply_async 一次只能在一個核心中算一個東西
    res = pool.apply_async(job,(2,)) # 這個只能把一個值放在一個核心運算一次
    print(res.get()) # 4 

    # 如果要輸入一串的話，可以用迭代的方法搭配 apply_async，達到多核心計算的效能
    multi_res = [pool.apply_async(job, (i,)) for i in range(10)]
    print([res.get() for res in multi_res]) # [0, 1, 4, 9, 16, 25, 36, 49, 64, 81]

if __name__=='__main__':
    multicore()
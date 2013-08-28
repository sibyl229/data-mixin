import numpy as np 


def column_stack(*ndArrs):
    '''ensure columns in these arrays are named'''
    # keep the dtypes of the original columns
    dtype_arr = []
    colCount = 0
    for arr in ndArrs:
        dtype_arr += arr.dtype.descr # an array of (fieldName, type) tuples
        colCount += len(arr.dtype.descr)

    # initialize the array
    shape = (ndArrs[0].shape[0],) 
    newArr = np.zeros(shape, dtype=dtype_arr)

    # populate new array based on old value
    for arr in ndArrs:
        for fieldName in arr.dtype.names:
            # import pdb; pdb.set_trace()
            newArr[fieldName] = arr[fieldName]

    return newArr

    





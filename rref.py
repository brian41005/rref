# -*- coding: utf-8 -*-
#! python3.5.1
import numpy as np
import sympy
#------------------------------------------------------------------------
def row_interchange(matrix, rows):
    matrix[rows, :] = matrix[rows[::-1], :]
    return matrix
#------------------------------------------------------------------------
def scaled(matrix, multiple, row=None):
    if multiple == 0:
        print ('error')
        return matrix
    if row == None:
        row = range(len(matrix))
    matrix[row, :] = multiple*matrix[row, :]
    return matrix
#------------------------------------------------------------------------
def row_reduced(matrix, pivot, reduce_row, reduced_row):
    for i in reduced_row:
        matrix[i, :] += -(matrix[i, pivot]/matrix[reduce_row, pivot])*matrix[reduce_row, :]
    return matrix
#------------------------------------------------------------------------
def get_columns_sort(matrix, rows, column):
    temp = matrix.copy()
    temp = np.argmax(np.absolute(np.squeeze(np.asarray(temp[rows, column]))), axis=0)
    return temp
#------------------------------------------------------------------------
def sort_out_row(matrix):
    for i in range(len(matrix)-1):
        temp_rows = range(i, len(matrix))
        rows = [i, get_columns_sort(matrix, temp_rows, i)+i]
        matrix = row_interchange(matrix, rows)
    return matrix
#------------------------------------------------------------------------
def column_is_all_zero(matrix, column):
    return np.isclose(np.sum(matrix[:, column]), 0)
#------------------------------------------------------------------------
def get_pivot(row):
    for i, element in enumerate(np.squeeze(np.asarray(row))):
        if np.isclose(element, 0):
            continue
        else:
            return i
#------------------------------------------------------------------------
def scale_pivot(matrix):
    for i, each_row in enumerate(matrix):
        pivot = matrix[i, get_pivot(each_row)]
        matrix[i,:] *= 1.0/pivot
    return matrix
#------------------------------------------------------------------------
def forward_phase(matrix):
    matrix = sort_out_row(matrix)
    (m, n) = matrix.shape
    for i in range(m-1):
        if column_is_all_zero(matrix, i):
            continue
        else:
            matrix = row_reduced(matrix, get_pivot(matrix[i,:]), i, range(i+1, m))
        matrix = sort_out_row(matrix)
    return matrix
#------------------------------------------------------------------------
def backward_phase(matrix):
    (m, n) = matrix.shape
    for i in reversed(range(m)):
        if column_is_all_zero(matrix, i):
            continue
        else:
            matrix = row_reduced(matrix, get_pivot(matrix[i,:]), i, range(0, i))
        matrix = sort_out_row(matrix)
    return scale_pivot(matrix)
#------------------------------------------------------------------------
def row_reduced_alogrithm(matrix):
    matrix = forward_phase(matrix)
    if check_is_consistent(matrix):
        matrix = backward_phase(matrix)
    else:
        return None
    return matrix
#------------------------------------------------------------------------
def check_is_consistent(matrix):
    (m, n) = matrix.shape
    if np.isclose(np.sum(matrix[m-1, :n-1]), 0) and not np.isclose(matrix[m-1, n-1], 0):
        return False
    else:
        return True
#------------------------------------------------------------------------
def check_if_is_unique_solution(matrix):
    (m, n) = matrix.shape
    coeffiecient_matrix = matrix[:,:n-1]
    return np.sum(coeffiecient_matrix) == m
#------------------------------------------------------------------------
def get_basic_variable(row, index):
    row = np.squeeze(np.asarray(row))
    temp = 'X%d ='%(index+1)
    if not np.isclose(row[-1], 0):
        temp += ' %.3f'%(row[-1])
    row = np.delete(row, -1)
    for i, content in enumerate(row):
        if not np.isclose(content, 0) and not np.isclose(content, 1):
            temp += ' %+.3fX%d'%(-content, i+1)
    return temp
#------------------------------------------------------------------------
def get_general_solution(matrix):
    result = row_reduced_alogrithm(matrix)
    print(result)
    if result is None:
        return 'inconsistent'
    elif check_if_is_unique_solution(result):
        solution = ''
        for i, row in enumerate(result):
            solution += 'X%d = %.3f\n'%(i+1, row[:, -1])
        return solution.strip()
    else:
        (m, n) = result.shape
        solution = ''
        for i in range(n-1):
            if not np.isclose(np.sum(result[:, i]), 1.):
                solution += 'X%d is free\n'%(i+1)
            else:
                solution += get_basic_variable(result[get_pivot(result[:, i]), :], i) + '\n'
        return solution.strip()
#------------------------------------------------------------------------
if __name__ == '__main__':
    augmented_matrices = {
    'u' : np.matrix('1,1,1,4;-1,-1,1,-2;2,-1,2,2', dtype ='float64'),
    'v' : np.matrix('1,-2,1,0;0,2,-8,8;-4,5,9,-9', dtype ='float64'),
    'w' : np.matrix('0,1,-4,8;2,-3,2,1;5,-8,7,1', dtype ='float64'),
    'x' : np.matrix('4,-8,-3,2,13;3,-4,-1,-3,5;2,-4,-2,2,6', dtype ='float64'),
    'y' : np.matrix('0,3,-6,6,4,-5;3,-7,8,-5,8,9;3,-9,12,-9,6,15', dtype ='float64'),
    'z' : np.matrix('0,3,-6,6,3,-5;3,-7,8,-5,8,9;3,-9,12,-9,6,15', dtype ='float64'),
    'r' : np.matrix('1,6,0,3,0,0;0,0,1,-4,0,5;0,0,0,0,1,7', dtype ='float64')}
    for name, matrix in augmented_matrices.items():
        print('%s matrix'%(name))
        result = get_general_solution(matrix)
        print('general solution:')
        print(result)
        print('*'*50)
    # result = get_general_solution(augmented_matrices['r'])
    # print('general solution:')
    # print(result)
    # print('*'*50)

    #print("第三方lib檢查")
    #print(sympy.Matrix(u).rref())

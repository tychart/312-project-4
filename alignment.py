def align(
        seq1: str,
        seq2: str,
        match_award=-3,
        indel_penalty=5,
        sub_penalty=1,
        banded_width=-1,
        gap='-'
) -> tuple[float, str | None, str | None]:
    """
        Align seq1 against seq2 using Needleman-Wunsch
        Put seq1 on left (j) and seq2 on top (i)
        => matrix[i][j]
        :param seq1: the first sequence to align; should be on the "left" of the matrix
        :param seq2: the second sequence to align; should be on the "top" of the matrix
        :param match_award: how many points to award a match
        :param indel_penalty: how many points to award a gap in either sequence
        :param sub_penalty: how many points to award a substitution
        :param banded_width: banded_width * 2 + 1 is the width of the banded alignment; -1 indicates full alignment
        :param gap: the character to use to represent gaps in the alignment strings
        :return: alignment cost, alignment 1, alignment 2
    """

    # print(f"Aligning sequences '{seq1}' and '{seq2}'")
    # print(f"Match award: {match_award}, indel penalty: {indel_penalty}, substitution penalty: {sub_penalty}")
    # print(f"Banded width: {banded_width}")

    # matrix = TwoDimensionalListManager

    # Compress the penalty values into a single dict for easier access
    penalties = {
        'match': match_award,
        'indel': indel_penalty,
        'sub': sub_penalty
    }

    if banded_width == -1:
        matrix = edit(penalties, seq1, seq2)
    else:
        matrix = banded_edit(penalties, seq1, seq2, banded_width)

    # print("Computed Matrix:")
    # print_matrix(matrix)

    path_tuple = find_path(penalties, gap, matrix, seq1, seq2)

    # print(f"Path Tuple: {path_tuple}")

    return path_tuple



def edit(penalties: dict, x: str, y: str) -> dict:
    matrix = {}
    
    for i in range(len(x) + 1):
        matrix[(i, 0)] = i * penalties["indel"]
    for j in range(len(y) + 1):
        matrix[(0, j)] = j * penalties["indel"]
    
    for i in range(1, len(x) + 1):
        for j in range(1, len(y) + 1):

            ## Broken out for debugging
            # print(f"Diag:  {diff(penalties, x[i - 1], y[j - 1]) + matrix[(i - 1, j - 1)]}")
            # print(f"up: {penalties["indel"] + matrix[(i, j - 1)]}")
            # print(f"Left: {penalties["indel"] + matrix[(i - 1, j)]}")

            # print(f"Matrix 0,1: {matrix[(0, 1)]}")

            matrix[(i, j)] = min(
                diff(penalties, x[i - 1], y[j - 1]) + matrix[(i - 1, j - 1)],
                penalties["indel"] + matrix[(i, j - 1)],
                penalties["indel"] + matrix[(i - 1, j)]
            )
            # print("--------------------------------")
            # print_matrix(matrix)
            # print("--------------------------------")

    return matrix

def banded_edit(penalties: dict, x: str, y: str, banded_width: int) -> dict:
    matrix = {}
        
    for i in range(banded_width + 1):
        matrix[(i, 0)] = i * penalties["indel"]
    for j in range(banded_width + 1):
        matrix[(0, j)] = j * penalties["indel"]
    
    for i in range(1, len(x) + 1):

        # if i - banded_width < 1:
        #     start = 1
        # else:
        #     start = i - banded_width

        # if i + banded_width > len(y):
        #     end = len(y)
        # else:
        #     end = i + banded_width

        for j in range(get_start(i, banded_width), get_end(i, banded_width, len(y)) + 1):

            diag = calc_diag(penalties, matrix, x, y, i, j)
            up = calc_up(penalties, matrix, i, j)
            left = calc_left(penalties, matrix, i, j)

            # ## Broken out for debugging
            # print(f"Diag:  {diag}")
            # print(f"up: {up}")
            # print(f"Left: {left}")

            # print(f"Matrix 0,1: {matrix[(0, 1)]}")

            matrix[(i, j)] = min(
                diag,
                up,
                left
            )
            # print("--------------------------------")
            # print_matrix(matrix)
            # print("--------------------------------")

    return matrix

# def find_path(penalties: dict, gap: str, matrix: dict, x: str, y: str) -> tuple[int, str, str]:
    
#     outstr1 = ""
#     outstr2 = ""

#     i = len(x)
#     j = len(y)
    
#     while i > 0 and j > 0:
        
#         # print("-------------------------------------------------")

#         diag = calc_diag(penalties, matrix, x, y, i, j)
#         up = calc_up(penalties, matrix, i, j)
#         left = calc_left(penalties, matrix, i, j)

#         lowest_cost = diag
#         lowest_direction = "diag"
#         path_to_next_x = -1
#         path_to_next_y = -1

#         if up < lowest_cost:
#             lowest_cost = up
#             lowest_direction = "up"
#             path_to_next_x = 0
#             path_to_next_y = -1
        
#         if left < lowest_cost:
#             lowest_cost = left
#             lowest_direction = "left"
#             path_to_next_x = -1
#             path_to_next_y = 0
            
#         # print(f"i: {i}, j: {j}, value: {matrix[(i, j)]}, lowest_direction: {lowest_direction}")

#         i += path_to_next_y
#         j += path_to_next_x

#         # print(lowest_direction)
        

#         if lowest_direction == "diag":
#             outstr1 = x[i] + outstr1
#             outstr2 = y[j] + outstr2
#         elif lowest_direction == "up":
#             outstr1 = x[i] + outstr1
#             outstr2 = gap + outstr2
#         elif lowest_direction == "left":
#             outstr1 = gap + outstr1
#             outstr2 = y[i] + outstr2
        
#         print(f"outstr1: {outstr1}")
#         print(f"outstr2: {outstr2}")

#         print("-------------------------------------------------")
    
#     while i > 0 or j > 0:
#         if i > 0:
#             outstr1 = x[i] + outstr1
#             i -= 1
#         else:
#             outstr1 = gap + outstr1


#         if j > 0:
#             outstr2 = y[j] + outstr2
#             j -= 1
#         else:
#             outstr2 = gap + outstr2


#     return (matrix[(len(x), len(y))], outstr1, outstr2)

def find_path(penalties: dict, gap: str, matrix: dict, x: str, y: str) -> tuple[int, str, str]:
    
    # print_matrix(matrix)

    backwardsoutstr1 = ""
    backwardsoutstr2 = ""

    i = len(x)
    j = len(y)

    while i > 0 or j > 0:
        if i > 0 and j > 0:
            diag = calc_diag(penalties, matrix, x, y, i, j)
            left = calc_left(penalties, matrix, i, j)
            up = calc_up(penalties, matrix, i, j)

            """
            In case there is more than one optimal alignment, 
            break ties with the following preference order: diagonal, left, top 
            (put sequence 1 on the left of the matrix and sequence 2 on the top of the matrix).
            """

            if matrix[(i, j)] == diag:
                backwardsoutstr1 += x[i - 1]
                backwardsoutstr2 += y[j - 1]
                i -= 1
                j -= 1
            elif matrix[(i, j)] == left:
                backwardsoutstr1 += gap
                backwardsoutstr2 += y[j - 1]
                j -= 1
            elif matrix[(i, j)] == up:
                backwardsoutstr1 += x[i - 1]
                backwardsoutstr2 += gap
                i -= 1
            else:
                raise Exception("Error in backtracking") # Should not happen
        elif i > 0:
            backwardsoutstr1 += x[i - 1]
            backwardsoutstr2 += gap
            i -= 1
        elif j > 0:
            backwardsoutstr1 += gap
            backwardsoutstr2 += y[j - 1]
            j -= 1

    outstr1 = backwardsoutstr1[::-1]
    outstr2 = backwardsoutstr2[::-1]

    return (matrix[(len(x), len(y))], outstr1, outstr2)
    



def get_start(i, banded_width):
    if i - banded_width < 1:
        return 1
    return i - banded_width

def get_end(i, banded_width, len_y):
    if i + banded_width > len_y:
        return len_y
    return i + banded_width

def calc_diag(penalties: dict, matrix: dict, x, y, i, j):
    try:
        return diff(penalties, x[i - 1], y[j - 1]) + matrix[(i - 1, j - 1)]
    except KeyError:
        return float("inf")

def calc_up(penalties: dict, matrix: dict, i, j):
    try:
        return penalties["indel"] + matrix[(i - 1, j)]
        # return penalties["indel"] + matrix[(i, j - 1)]
    except KeyError:
        return float("inf")

def calc_left(penalties: dict, matrix: dict, i, j):
    try:
        return penalties["indel"] + matrix[(i, j - 1)]
        # return penalties["indel"] + matrix[(i - 1, j)]
    except KeyError:
        return float("inf")

def diff(penalties: dict, char1, char2):
    return penalties["match"] if char1 == char2 else penalties["sub"]




def print_matrix(matrix: dict):
    """
    Print the matrix in a human-readable format, pad-aligned for up to 3 digits
    (or 2 digits and a negative sign).
    :param matrix: the matrix to print
    """
    if not matrix:
        print("Matrix is empty.")
        return

    # Determine the maximum row and column indices to set the printing boundaries
    max_row = 0
    max_col = 0
    max_val_len = 0
    for (row, col), value in matrix.items():
        max_row = max(max_row, row)
        max_col = max(max_col, col)
        max_val_len = max(max_val_len, len(str(value)))

    padding = max(3, max_val_len)  # Ensure at least 3 spaces for padding

    for i in range(max_row + 1):
        for j in range(max_col + 1):
            value = matrix.get((i, j), ' ')
            print(f"{value:>{padding}}", end=' ')
        print()
    

# class TwoDimensionalListManager:
#     def __init__(self, rows, cols):
#         """
#         Initializes a new instance of the TwoDimensionalListManager class.
        
#         Args:
#             rows (int): The number of rows in the matrix.
#             cols (int): The number of columns in the matrix.
#         """
#         if rows <= 0 or cols <= 0:
#             raise ValueError("Both rows and columns must be positive.")
        
#         self.list_2d = [[0 for _ in range(cols)] for _ in range(rows)]

#     def add_row(self, row):
#         """
#         Adds a new row to the 2D list.
        
#         Args:
#             row (list): The new row to be added.
#         """
#         if len(row) != len(self.list_2d[0]):
#             raise ValueError("All rows in the 2D list must have the same number of columns.")
#         self.list_2d.append(row)

#     def get_row(self, index):
#         """
#         Gets a specific row from the 2D list.
        
#         Args:
#             index (int): The index of the row to be retrieved.
        
#         Returns:
#             list: The specified row.
#         """
#         if index < 0 or index >= len(self.list_2d):
#             raise IndexError("Index out of range.")
#         return self.list_2d[index]

#     def get_column(self, column_index):
#         """
#         Gets a specific column from the 2D list.
        
#         Args:
#             column_index (int): The index of the column to be retrieved.
        
#         Returns:
#             list: The specified column.
#         """
#         if column_index < 0 or column_index >= len(self.list_2d[0]):
#             raise IndexError("Index out of range.")
#         return [row[column_index] for row in self.list_2d]

#     def print_list(self):
#         """
#         Prints the entire 2D list.
#         """
#         for row in self.list_2d:
#             print(row)


# # Example usage:

# manager = TwoDimensionalListManager()

# # Add rows to the 2D list
# manager.add_row([1, 2, 3])
# manager.add_row([4, 5, 6])

# # Get a specific row
# print("Row at index 0:", manager.get_row(0))

# # Get a specific column
# print("Column at index 1:", manager.get_column(1))

# # Print the entire list
# manager.print_list()






# matrix = {}

# matrix[(0, 0)] = 1
# matrix[(0, 1)] = 0
# # matrix[(0, 2)] = 1
# # matrix[(0, 3)] = 0
# # matrix[(0, 4)] = 1

# # matrix[(1, 0)] = 0
# matrix[(1, 1)] = 1
# matrix[(1, 2)] = 0
# # matrix[(1, 3)] = 1
# # matrix[(1, 4)] = 0

# # matrix[(2, 0)] = 1
# # matrix[(2, 1)] = 0
# matrix[(2, 2)] = 1
# matrix[(2, 3)] = 0
# # matrix[(2, 4)] = 1

# # matrix[(3, 0)] = 0
# # matrix[(3, 1)] = 1
# # matrix[(3, 2)] = 0
# matrix[(3, 3)] = 1
# matrix[(3, 4)] = 0

# # matrix[(4, 0)] = 1
# # matrix[(4, 1)] = 0
# # matrix[(4, 2)] = 1
# # matrix[(4, 3)] = 0
# matrix[(4, 4)] = 1


# print_matrix(matrix)

# align("THARS", "OTHER")
# align("ATGCATGC", "ATGGTGC", banded_width=3)
# align("ATGCATGC", "ATGGTGC")



# def test_small_alignment(align):
#     score, aseq1, aseq2 = align('polynomial', 'exponential')
#     if score == -1 and aseq1 == 'polyn-omial' and aseq2 == 'exponential':
#         print("Test passed: Alignment is correct.")
#     else:
#         print("Test failed:")
#         if score != -1:
#             print(f"  Expected score: -1, but got: {score}")
#         if aseq1 != 'polyn-omial':
#             print(f"  Expected aseq1: 'polyn-omial', but got: {aseq1}")
#         if aseq2 != 'exponential':
#             print(f"  Expected aseq2: 'exponential', but got: {aseq2}")

# def test_small_alignment_backtrack(align):
#     score, aseq1, aseq2 = align('ATATATATAT', 'TATATATATA')
#     if score == -17 and aseq1 == 'ATATATATAT-' and aseq2 == '-TATATATATA':
#         print("Test passed: Alignment is correct.")
#     else:
#         print("Test failed:")
#         if score != -17:
#             print(f"  Expected score: -17, but got: {score}")
#         if aseq1 != 'ATATATATAT-':
#             print(f"  Expected aseq1: 'ATATATATAT-', but got: {aseq1}")
#         if aseq2 != '-TATATATATA':
#             print(f"  Expected aseq2: '-TATATATATA', but got: {aseq2}") 
#     # assert score == -17
#     # assert aseq1 == 'ATATATATAT-'
#     # assert aseq2 == '-TATATATATA'
# test_small_alignment(align)
# test_small_alignment_backtrack(align)

# def test_tiny_dna_alignment(align):
#     score, aseq1, aseq2 = align('ATGCATGC', 'ATGGTGC')
#     assert score == -12
#     assert aseq1 == 'ATGCATGC'
#     assert aseq2 == 'ATG-GTGC'

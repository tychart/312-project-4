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

    # matrix = TwoDimensionalListManager













    

class TwoDimensionalListManager:
    def __init__(self, rows, cols):
        """
        Initializes a new instance of the TwoDimensionalListManager class.
        
        Args:
            rows (int): The number of rows in the matrix.
            cols (int): The number of columns in the matrix.
        """
        if rows <= 0 or cols <= 0:
            raise ValueError("Both rows and columns must be positive.")
        
        self.list_2d = [[0 for _ in range(cols)] for _ in range(rows)]

    def add_row(self, row):
        """
        Adds a new row to the 2D list.
        
        Args:
            row (list): The new row to be added.
        """
        if len(row) != len(self.list_2d[0]):
            raise ValueError("All rows in the 2D list must have the same number of columns.")
        self.list_2d.append(row)

    def get_row(self, index):
        """
        Gets a specific row from the 2D list.
        
        Args:
            index (int): The index of the row to be retrieved.
        
        Returns:
            list: The specified row.
        """
        if index < 0 or index >= len(self.list_2d):
            raise IndexError("Index out of range.")
        return self.list_2d[index]

    def get_column(self, column_index):
        """
        Gets a specific column from the 2D list.
        
        Args:
            column_index (int): The index of the column to be retrieved.
        
        Returns:
            list: The specified column.
        """
        if column_index < 0 or column_index >= len(self.list_2d[0]):
            raise IndexError("Index out of range.")
        return [row[column_index] for row in self.list_2d]

    def print_list(self):
        """
        Prints the entire 2D list.
        """
        for row in self.list_2d:
            print(row)


# Example usage:

manager = TwoDimensionalListManager()

# Add rows to the 2D list
manager.add_row([1, 2, 3])
manager.add_row([4, 5, 6])

# Get a specific row
print("Row at index 0:", manager.get_row(0))

# Get a specific column
print("Column at index 1:", manager.get_column(1))

# Print the entire list
manager.print_list()
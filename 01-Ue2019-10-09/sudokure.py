# python3
# recursive sudoku solver

# import sudoku from file and store as grid

# format live text
def print_counter(val, msg):
    print "%s[%d] %s" % (" "*val, val, msg)

# check if the sudoku is valid
def check_sudoku(__grid)
    # something like is_unique?
            
def solve_sudoku(__grid, counter=0):

    res = check_sudoku(__grid)
    if res is None or res is False:
        return res

    grid = copy.deepcopy(__grid)

    for row in xrange(9):
        for col in xrange(9):
            if grid[row][col] == 0:
                for n in xrange(1, 10):
                    grid[row][col] = n
                    print_counter(counter,"test: (row %d, col %d) = %d" % (row,col,n))
                    new = solve_sudoku(grid, counter+1)
                    if new is not False:
                        print_counter(counter, "solve_sudoku() solved: returning")
                        return new
                # backtrack
                print_counter(counter, "backtrack")
                return False

    print_counter(counter, "**SOLVED! Returning back up to top**")
    return grid

from pprint import pprint 
solution = solve_sudoku(easy_grid)
pprint(solution)

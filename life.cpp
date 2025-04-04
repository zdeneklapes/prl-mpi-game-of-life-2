/*
 * Game of Life
 * Author: Zdenek Lapes <lapes.zdenek@gmail.com> (xlapes02)
 * Date: 2024-04-14
 * Description: Parallel implementation of the Game of Life using MPI. (Wrap-around)
 * Compile: mpic++ -o life life.cpp
 * Run: mpirun -np <number of processes> ./life <path to the game field file> <number of steps>
 *  <number of processes> accept any number of processes available on the system
 *  <path to the game field file> is the path to the file with the initial state of the game field
 *  <number of steps> is the number of steps the simulation should run
 * Output: The state of the game field after the specified number of steps (stdout)
 */

#include <mpi.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstring>


#define DEBUG 0
#define DEBUG_LITE 0
#define DEBUG_PRINT_LITE(fmt, ...) \
            do { if (DEBUG_LITE) fprintf(stderr, fmt, __VA_ARGS__); } while (0)
#define DEBUG_PRINT(fmt, ...) \
        do { if (DEBUG) fprintf(stderr, "%s:%d:%s(): " fmt, __FILE__, \
                                __LINE__, __func__, __VA_ARGS__); } while (0)
#define DEBUG_PRINT_LITE_PROCESS(process, fmt, ...) \
            do { if (DEBUG_LITE && program->process_id == process) fprintf(stderr, fmt, __VA_ARGS__); } while (0)
#define DEBUG_PRINT_PROCESS(process, fmt, ...) \
        do { if (DEBUG && program->process_id == process) fprintf(stderr, "%s:%d:%s(): " fmt, __FILE__, \
                                __LINE__, __func__, __VA_ARGS__); } while (0)

const int ALIVE = 1;
const int DEAD = 0;

void die() {
    std::cerr << "Usage: ./life <path to the game field file> <number of steps>" << std::endl;
    MPI_Abort(MPI_COMM_WORLD, 1);
    exit(1);
}

std::vector<std::vector<int>> *get_grid_from_file(char *filename) {
    MPI_File mpi_file;
    MPI_Status status;
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Open the file collectively
    if (MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &mpi_file) != MPI_SUCCESS) {
        if (rank == 0) {
            std::cerr << "File not found: " << filename << std::endl;
        }
        die();
    }

    // find file size
    MPI_Offset file_size;
    MPI_File_get_size(mpi_file, &file_size);

    // read the entire file into a buffer
    char *buffer = new char[file_size + 1];
    MPI_File_read_all(mpi_file, buffer, file_size, MPI_CHAR, &status);

    // null-terminate the buffer
    buffer[file_size] = '\0';

    // close the file
    MPI_File_close(&mpi_file);

    // Now, each process has the whole file in memory. You can distribute parsing or parse conditionally based on rank.
    auto grid = new std::vector<std::vector<int>>();
    char *line = strtok(buffer, "\n");
    while (line != nullptr) {
        std::vector<int> row;
        for (int i = 0; line[i] != '\0'; i++) {
            if (line[i] == '0') {
                row.push_back(0); // Dead cell
            } else if (line[i] == '1') {
                row.push_back(1); // Alive cell
            }
        }
        grid->push_back(row);
        line = strtok(nullptr, "\n");
    }

    // Clean up
    delete[] buffer;

    return grid;
}


/**
 * Program class
 * @tparam T
 */
class Program {
public:
    /**
     * Constructor
     * @param argc: Number of arguments
     * @param argv: Arguments
     */
    Program(int argc, char *argv[]) {
        if (argc != 3) { die(); }
        this->filename = argv[1];

        this->grid = get_grid_from_file(filename);

        // parse steps
        try {
            this->steps = std::stoi(argv[2]);
            if (this->steps <= 0) {
                throw std::invalid_argument("Number of steps must be positive");
            }
        } catch (std::invalid_argument &e) {
            std::cerr << "Invalid number of steps: " << argv[2] << std::endl;
            die();
        } catch (std::out_of_range &e) {
            std::cerr << "Number of steps is out of range: " << argv[2] << std::endl;
            die();
        }
    }

    /**
     * Error message and exit
     */

    char *filename = nullptr;
    int steps = 0;
    int process_id = -1;
    int processes = -1;
    std::vector<std::vector<int>> *grid = nullptr;
};

/**
 * Update cell state
 *
 * @param program: Program
 * @param row: Row
 * @param col: Column
 * @param next_grid: Next grid
 * @return void
 */
void update_cell_state(Program *program, int row, int col, std::vector<std::vector<int>> &next_grid) {
    int alive_neighbors = 0;
    int nrows = program->grid->size();
    int ncols = program->grid->at(0).size();

    // count neighbors
    for (int i = -1; i <= 1; ++i) {
        int ni = row + i;
        if (row == 0 && i == -1) {
            ni = nrows - 1;
        }
        if (row == nrows - 1 && i == 1) {
            ni = 0;
        }
        for (int j = -1; j <= 1; ++j) {
            int nj = (col + j + ncols) % ncols;  // wrap around the columns
            if (ni == row && nj == col) {
                continue;  // skip the cell itself
            }
            if ((*program->grid)[ni][nj] == ALIVE) {
                alive_neighbors++;
            }
        }
    }

    // apply Game of Life rules
    int new_state = program->grid->at(row).at(col);
    if ((*program->grid)[row][col] == ALIVE) {
        if (alive_neighbors < 2) {
            new_state = DEAD;
        } else if (alive_neighbors == 2 || alive_neighbors == 3) {
            new_state = ALIVE;
        } else if (alive_neighbors > 3) {
            new_state = DEAD;
        }
    } else if ((*program->grid)[row][col] == DEAD) {
        if (alive_neighbors == 3) {
            new_state = ALIVE;
        }
    } else {
        std::cerr << "Invalid cell state: " << (*program->grid)[row][col] << std::endl;
        die();
    }
//    DEBUG_PRINT_LITE("process_id: %d, row: %d, col: %d, alive_neighbors: %d, current_state: %d, new_state: %d\n", program->process_id, row, col, alive_neighbors, program->grid->at(row).at(col), new_state);
    next_grid[row][col] = new_state;
}

/**
 * Copy grid
 *
 * @param from: From grid
 * @param to: To grid
 * @return void
 */
void copy_grid(std::vector<std::vector<int>> &from, std::vector<std::vector<int>> &to) {
    for (int i = 0; i < from.size(); i++) {
        for (int j = 0; j < from[i].size(); j++) {
            to[i][j] = from[i][j];
        }
    }
}

int get_start_row(int process_id, int nrows, int process_count) {
    int rows_per_process = nrows / process_count;
    int start_row = process_id * rows_per_process;
    return start_row;
}

int get_end_row(int process_id, int nrows, int process_count) {
    int rows_per_process = nrows / process_count;
    int start_row = get_start_row(process_id, nrows, process_count);
    int end_row = (process_id == process_count - 1) ? nrows - 1 : start_row + rows_per_process - 1;
    return end_row;
}

/**
 * Print grid to stdout
 *
 * @param program: Program
 * @return void
 */
void print_grid(Program *program) {
    int nrows = program->grid->size();
    int ncols = program->grid->at(0).size();

    // get the start and end rows of each process
    int start_row = get_start_row(program->process_id, nrows, program->processes);
    int end_row = get_end_row(program->process_id, nrows, program->processes);
    int local_rows = end_row - start_row + 1;

    // prepare buffer for gathering grid data
    std::vector<int> local_grid(local_rows * ncols);

    // fill local grid data
    for (int i = 0; i < local_rows; ++i) {
        for (int j = 0; j < ncols; ++j) {
            local_grid[i * ncols + j] = program->grid->at(start_row + i).at(j);
        }
    }

    // prepare to gather grids at the root process
    std::vector<int> complete_grid;
    std::vector<int> recv_counts(program->processes);
    std::vector<int> displacements(program->processes);

    // only the root process needs to allocate space for the complete grid
    if (program->process_id == 0) {
        complete_grid.resize(nrows * ncols);
    }

    // all processes must provide the count of elements they will send
    MPI_Gather(&local_rows, 1, MPI_INT, recv_counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

    // calculate displacements for each process at the root
    if (program->process_id == 0) {
        int offset = 0;
        for (int i = 0; i < program->processes; ++i) {
            displacements[i] = offset;
            recv_counts[i] *= ncols;  // convert row count to element count
            offset += recv_counts[i];
        }
    }

    // gather all local grids to the root process using variable sizes
    MPI_Gatherv(local_grid.data(), local_rows * ncols, MPI_INT,
                complete_grid.data(), recv_counts.data(), displacements.data(), MPI_INT,
                0, MPI_COMM_WORLD);

    // only the root process prints the complete grid
    if (program->process_id == 0) {
        for (int row = 0; row < nrows; ++row) {
            // get process that owns this row
            int process_id = row / (nrows / program->processes);
            if (process_id >= program->processes) {
                process_id = program->processes - 1;
            }
            std::cout << process_id << ": ";
            for (int col = 0; col < ncols; ++col) {
                std::cout << complete_grid[row * ncols + col];
            }
            std::cout << std::endl;
        }
    }
}

/**
 * Run simulation
 *
 * @param program: Program
 * @return void
 */
void run_simulation(Program *program) {
    int nrows = program->grid->size();
    int process_count = program->processes;
    if (nrows < process_count) {
        process_count = nrows;
    }

    int ncols = program->grid->at(0).size();
    std::vector<std::vector<int>> next_grid(nrows, std::vector<int>(ncols, DEAD));

    // all from program->grid to next_grid
    copy_grid(*program->grid, next_grid);

    // determine rows per process
    int start_row = get_start_row(program->process_id, nrows, process_count);
    int end_row = get_end_row(program->process_id, nrows, process_count);

    // constants for clarity in communication
    const int UP = 0;
    const int DOWN = 1;

    bool is_first_process = program->process_id == 0;
    bool is_last_process = program->process_id == process_count - 1;
    bool is_middle_process = !is_first_process && !is_last_process;

    MPI_Request send_requests[2];
    MPI_Request recv_requests[2];

    for (int step = 0; step < program->steps; ++step) {
        MPI_Barrier(MPI_COMM_WORLD);
        if (program->process_id >= nrows) {
            continue;
        }
        for (int row = start_row; row <= end_row; ++row) {
            for (int col = 0; col < ncols; ++col) {
                update_cell_state(program, row, col, next_grid);
            }
        }

        copy_grid(next_grid, *program->grid);

        if (is_first_process) {
            // exchange data with process above
            int last_process = process_count - 1;
            MPI_Isend(&((*program->grid)[start_row][0]), ncols, MPI_INT, last_process, DOWN, MPI_COMM_WORLD,
                      &send_requests[UP]);
            MPI_Irecv(&((*program->grid)[nrows - 1][0]), ncols, MPI_INT, last_process, UP, MPI_COMM_WORLD,
                      &recv_requests[UP]);

            // exchange data with process below
            int next_process = program->process_id + 1;
            MPI_Isend(&((*program->grid)[end_row][0]), ncols, MPI_INT, next_process, UP, MPI_COMM_WORLD,
                      &send_requests[DOWN]);
            MPI_Irecv(&((*program->grid)[end_row + 1][0]), ncols, MPI_INT, next_process, DOWN, MPI_COMM_WORLD,
                      &recv_requests[DOWN]);
        } else if (is_middle_process) {
            // exchange data with process above
            int previous_process = program->process_id - 1;
            MPI_Isend(&((*program->grid)[start_row][0]), ncols, MPI_INT, previous_process, DOWN, MPI_COMM_WORLD,
                      &send_requests[UP]);
            MPI_Irecv(&((*program->grid)[start_row - 1][0]), ncols, MPI_INT, previous_process, UP, MPI_COMM_WORLD,
                      &recv_requests[UP]);

            // exchange data with process below
            int next_process = program->process_id + 1;
            MPI_Isend(&((*program->grid)[end_row][0]), ncols, MPI_INT, next_process, UP, MPI_COMM_WORLD,
                      &send_requests[DOWN]);
            MPI_Irecv(&((*program->grid)[end_row + 1][0]), ncols, MPI_INT, next_process, DOWN, MPI_COMM_WORLD,
                      &recv_requests[DOWN]);
        } else if (is_last_process) {
            // exchange data with process above
            int previous_process = program->process_id - 1;
            MPI_Isend(&((*program->grid)[start_row][0]), ncols, MPI_INT, previous_process, DOWN, MPI_COMM_WORLD,
                      &send_requests[UP]);
            MPI_Irecv(&((*program->grid)[start_row - 1][0]), ncols, MPI_INT, previous_process, UP, MPI_COMM_WORLD,
                      &recv_requests[UP]);

            // exchange data with process below
            int first_process = 0;
            MPI_Isend(&((*program->grid)[end_row][0]), ncols, MPI_INT, first_process, UP, MPI_COMM_WORLD, &send_requests[DOWN]);
            MPI_Irecv(&((*program->grid)[0][0]), ncols, MPI_INT, first_process, DOWN, MPI_COMM_WORLD, &recv_requests[DOWN]);
        }
    }

    if (program->process_id >= nrows) {
        return;
    } else {
        print_grid(program);
    }
}

/**
 * Main function
 *
 * @param argc: Number of arguments
 * @param argv: Arguments
 * @return int
 */
int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    auto program = new Program(argc, argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &program->process_id);
    MPI_Comm_size(MPI_COMM_WORLD, &program->processes);
    MPI_Barrier(MPI_COMM_WORLD);
    run_simulation(program);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    delete program;
    return 0;
}

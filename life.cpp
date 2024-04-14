#include <mpi.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>

#define INPUT_FILE "./numbers"
#define DEBUG 1
#define DEBUG_LITE 1
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

enum COORDINATES {
    ROWS, COLUMNS, DIMENSION
};

template<typename T>
class Program {
public:
    Program(int argc, char *argv[]) {
        if (argc != 3) { die(); }
        this->filename = argv[1];

        // Read file
        std::ifstream file_stream(this->filename);
        if (!file_stream.good()) {
            std::cerr << "File not found: " << this->filename << std::endl;
            die();
        } else {
            // Parse input
            std::string line;
            this->grid = new std::vector<std::vector<int>>();
            while (getline(file_stream, line)) {
                std::vector<int> row;
                for (char cell: line) {
                    if (cell == '0') {
                        row.push_back(0); // Dead cell
                    } else if (cell == '1') {
                        row.push_back(1); // Alive cell
                    }
                }
                this->grid->push_back(row);
            }
            file_stream.close();
        }
        if (this->grid->empty()) {
            std::cerr << "Empty grid" << std::endl;
            die();
        }

        // Parse steps
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

    void die() {
        std::cerr << "Usage: ./life <path to the game field file> <number of steps>" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
        exit(1);
    }

    char *filename = nullptr;
    int steps = 0;
    int process_id = -1;
    int processes = -1;
    std::vector<std::vector<int>> *grid = nullptr;
};

void update_cell_state(Program<void> *program, int row, int col, std::vector<std::vector<int>> &next_grid) {
    int alive_neighbors = 0;
    int nrows = program->grid->size();
    int ncols = program->grid->at(0).size();
    // Count neighbors with periodic boundary conditions
    for (int i = -1; i <= 1; ++i) {
        int ni = row + i;
        if (row == 0 && i == -1) ni = nrows - 1;
        if (row == nrows - 1 && i == 1) ni = 0;
        for (int j = -1; j <= 1; ++j) {
            int nj = (col + j + ncols) % ncols;  // wrap around the columns
            if (ni == row && nj == col) continue;  // Skip the cell itself
            if ((*program->grid)[ni][nj] == ALIVE) {
                alive_neighbors++;
            }
        }
    }

    // Apply Game of Life rules
    int new_state = -1;
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
        program->die();
    }
    if (new_state != -1) {
        DEBUG_PRINT_LITE("process_id: %d, row: %d, col: %d, alive_neighbors: %d, new_state: %d\n", program->process_id, row, col, alive_neighbors, new_state);
        next_grid[row][col] = new_state;
    }
}

void run_simulation(Program<void> *program) {
    int nrows = program->grid->size();
    int ncols = program->grid->at(0).size();
    std::vector<std::vector<int>> next_grid(nrows, std::vector<int>(ncols, DEAD));

    // all from program->grid to next_grid
    for (int i = 0; i < nrows; i++) {
        for (int j = 0; j < ncols; j++) {
            next_grid[i][j] = program->grid->at(i).at(j);
        }
    }

    // Determine rows per process
    int rows_per_process = nrows / program->processes;
    int start_row = program->process_id * rows_per_process;
    int end_row = (program->process_id == program->processes - 1) ? nrows - 1 : start_row + rows_per_process - 1;

    // Constants for clarity in communication
    const int UP = 0;
    const int DOWN = 1;

    MPI_Request send_requests[2];
    MPI_Request recv_requests[2];

    int num_processes = program->processes;
    int my_rank = program->process_id;

    bool is_first_process = my_rank == 0;
    bool is_middle_process = my_rank > 0 && my_rank < num_processes - 1;
    bool is_last_process = my_rank == num_processes - 1;

    // Main simulation loop
    for (int step = 0; step < program->steps; ++step) {
        for (int row = start_row; row <= end_row; ++row) {
            for (int col = 0; col < ncols; ++col) {
                update_cell_state(program, row, col, next_grid);
            }
        }

        // Update the current grid to the next grid
        std::swap(*program->grid, next_grid);

        if (is_first_process) {
            auto end_row_address_current = &((*program->grid)[end_row][0]);
            auto end_row_address_previous = &((*program->grid)[nrows - 1][0]);
            MPI_Isend(end_row_address_current, ncols, MPI_INT, (my_rank + 1) % num_processes, DOWN, MPI_COMM_WORLD,
                      &send_requests[DOWN]);
            MPI_Irecv(end_row_address_previous, ncols, MPI_INT, num_processes - 1, DOWN, MPI_COMM_WORLD,
                      &recv_requests[DOWN]);

            auto start_row_address_current = &((*program->grid)[start_row][0]);
            auto index = (my_rank + 1) * rows_per_process;
            auto start_row_address_next = &((*program->grid)[index][0]);
            MPI_Isend(start_row_address_current, ncols, MPI_INT, num_processes - 1, UP, MPI_COMM_WORLD,
                      &send_requests[UP]);
            MPI_Irecv(start_row_address_next, ncols, MPI_INT, (my_rank + 1) % num_processes, UP, MPI_COMM_WORLD,
                      &recv_requests[UP]);
        } else if (is_middle_process) {
            auto start_row_address_current = &((*program->grid)[start_row][0]);
            auto start_row_address_previous = &((*program->grid)[start_row - 1][0]);
            MPI_Isend(start_row_address_current, ncols, MPI_INT, 0, UP, MPI_COMM_WORLD, &send_requests[UP]);
            MPI_Irecv(start_row_address_previous, ncols, MPI_INT, my_rank - 1, UP, MPI_COMM_WORLD, &recv_requests[UP]);

            auto end_row_address_current = &((*program->grid)[end_row][0]);
            auto end_row_address_next = &((*program->grid)[0][0]);
            MPI_Isend(end_row_address_current, ncols, MPI_INT, 0, DOWN, MPI_COMM_WORLD, &send_requests[DOWN]);
            MPI_Irecv(end_row_address_next, ncols, MPI_INT, my_rank - 1, DOWN, MPI_COMM_WORLD, &recv_requests[DOWN]);
        } else if (is_last_process) {
            auto start_row_address_current = &((*program->grid)[start_row][0]);
            auto index = (my_rank - 1) * rows_per_process;
            auto start_row_address_previous = &((*program->grid)[index][0]);
            MPI_Isend(start_row_address_current, ncols, MPI_INT, my_rank - 1, UP, MPI_COMM_WORLD, &send_requests[UP]);
            MPI_Irecv(start_row_address_previous, ncols, MPI_INT, my_rank - 1, UP, MPI_COMM_WORLD, &recv_requests[UP]);

            auto end_row_address_current = &((*program->grid)[end_row][0]);
            auto end_row_address_next = &((*program->grid)[0][0]);
            MPI_Isend(end_row_address_current, ncols, MPI_INT, (my_rank + 1) % num_processes, DOWN, MPI_COMM_WORLD,
                      &send_requests[DOWN]);
            MPI_Irecv(end_row_address_next, ncols, MPI_INT, (my_rank + 1) % num_processes, DOWN, MPI_COMM_WORLD,
                      &recv_requests[DOWN]);
        }

        // Synchronize all processes
        MPI_Barrier(MPI_COMM_WORLD
        );
    }
}

void print_grid(Program<void> *program) {
    for (int i = 0; i < program->grid->size(); i++) {
        for (int j = 0; j < program->grid->at(i).size(); j++) {
            DEBUG_PRINT_LITE_PROCESS(0, "%d ", program->grid->at(i).at(j));
        }
        DEBUG_PRINT_LITE_PROCESS(0, "%c", '\n');
    }
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    auto program = new Program<void>(argc, argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &program->process_id);
    MPI_Comm_size(MPI_COMM_WORLD, &program->processes);

    DEBUG_PRINT_LITE_PROCESS(0, "process_id: %d\n", program->process_id);
    DEBUG_PRINT_LITE_PROCESS(0, "filename: %s\n", program->filename);
    DEBUG_PRINT_LITE_PROCESS(0, "steps: %d\n", program->steps);
    DEBUG_PRINT_LITE_PROCESS(0, "grid (rows x columns): %ld x %ld\n", program->grid->size(),
                             program->grid->at(0).size());
    // print whole grid
    print_grid(program);
    DEBUG_PRINT_LITE_PROCESS(0, "processes: %d\n", program->processes);
    DEBUG_PRINT_LITE_PROCESS(0, "====================================%c", '\n');
    MPI_Barrier(MPI_COMM_WORLD);
    run_simulation(program);
    print_grid(program);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    delete program;
    return 0;
}

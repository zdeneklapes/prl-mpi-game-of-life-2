#include <mpi.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>

#define INPUT_FILE "./numbers"
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
        if (row == 0 && i == -1) {
            ni = nrows - 1;
        }
        if (row == nrows - 1 && i == 1) {
            ni = 0;
        }
        for (int j = -1; j <= 1; ++j) {
            int nj = (col + j + ncols) % ncols;  // wrap around the columns
            if (ni == row && nj == col) {
                continue;  // Skip the cell itself
            }
            if ((*program->grid)[ni][nj] == ALIVE) {
//                DEBUG_PRINT_LITE_PROCESS(0, "process_id: %d, row: %d, col: %d, ni: %d, nj: %d\n", program->process_id, row, col, ni, nj);
                alive_neighbors++;
            }
        }
    }

    // Apply Game of Life rules
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
        program->die();
    }
    DEBUG_PRINT_LITE("process_id: %d, row: %d, col: %d, alive_neighbors: %d, current_state: %d, new_state: %d\n", program->process_id, row, col, alive_neighbors, program->grid->at(row).at(col), new_state);
    next_grid[row][col] = new_state;
}

void copy_grid(std::vector<std::vector<int>> &from, std::vector<std::vector<int>> &to) {
    for (int i = 0; i < from.size(); i++) {
        for (int j = 0; j < from[i].size(); j++) {
            to[i][j] = from[i][j];
        }
    }
}

void print_grid(Program<void> *program) {
    for (int i = 0; i < program->grid->size(); i++) {
        for (int j = 0; j < program->grid->at(i).size(); j++) {
            fprintf(stdout, "%d", program->grid->at(i).at(j));
        }
        fprintf(stdout, "%c", '\n');
    }
}

void run_simulation(Program<void> *program) {
    int nrows = program->grid->size();
    int ncols = program->grid->at(0).size();
    std::vector<std::vector<int>> next_grid(nrows, std::vector<int>(ncols, DEAD));

    // all from program->grid to next_grid
    copy_grid(*program->grid, next_grid);

    // Determine rows per process
    int rows_per_process = nrows / program->processes;
    int start_row = program->process_id * rows_per_process;
    int end_row = (program->process_id == program->processes - 1) ? nrows - 1 : start_row + rows_per_process - 1;

    // Constants for clarity in communication
    const int UP = 0;
    const int DOWN = 1;

    MPI_Request send_requests[2];
    MPI_Request recv_requests[2];

    bool is_first_process = program->process_id == 0;
    bool is_middle_process = program->process_id > 0 && program->process_id < program->processes - 1;
    bool is_last_process = program->process_id == program->processes - 1;

    // Main simulation loop
    if (program->process_id == 0) {
        print_grid(program);
        fprintf(stdout, "%c", '\n');
    }
    for (int step = 0; step < program->steps; ++step) {
//        DEBUG_PRINT_LITE("process_id: %d, step: %d, start_row: %d, end_row: %d\n", program->process_id, step, start_row, end_row);
        for (int row = start_row; row <= end_row; ++row) {
            for (int col = 0; col < ncols; ++col) {
                update_cell_state(program, row, col, next_grid);
            }
        }

        copy_grid(next_grid, *program->grid);

        if (program->process_id == 0) {  // First process
            // exchange data with process above
            int last_process = program->processes - 1;
            MPI_Isend(&((*program->grid)[start_row][0]), ncols, MPI_INT, last_process, DOWN, MPI_COMM_WORLD, &send_requests[UP]);
            MPI_Irecv(&((*program->grid)[nrows-1][0]), ncols, MPI_INT, last_process, UP, MPI_COMM_WORLD, &recv_requests[UP]);

            // exchange data with process below
            int next_process = program->process_id + 1;
            MPI_Isend(&((*program->grid)[end_row][0]), ncols, MPI_INT, next_process, UP, MPI_COMM_WORLD, &send_requests[DOWN]);
            MPI_Irecv(&((*program->grid)[end_row + 1][0]), ncols, MPI_INT, next_process, DOWN, MPI_COMM_WORLD, &recv_requests[DOWN]);
        } else if (program->process_id > 0 && program->process_id < program->processes - 1) {  // Middle processes
            // exchange data with process above
            int previous_process = program->process_id - 1;
            MPI_Isend(&((*program->grid)[start_row][0]), ncols, MPI_INT, previous_process, DOWN, MPI_COMM_WORLD, &send_requests[UP]);
            MPI_Irecv(&((*program->grid)[start_row - 1][0]), ncols, MPI_INT, previous_process, UP, MPI_COMM_WORLD, &recv_requests[UP]);

            // exchange data with process below
            int next_process = program->process_id + 1;
            MPI_Isend(&((*program->grid)[end_row][0]), ncols, MPI_INT, next_process, UP, MPI_COMM_WORLD, &send_requests[DOWN]);
            MPI_Irecv(&((*program->grid)[end_row + 1][0]), ncols, MPI_INT, next_process, DOWN, MPI_COMM_WORLD, &recv_requests[DOWN]);
        } else if (program->process_id == program->processes - 1) {  // Last process
            // exchange data with process above
            int previous_process = program->process_id - 1;
            MPI_Isend(&((*program->grid)[start_row][0]), ncols, MPI_INT, previous_process, DOWN, MPI_COMM_WORLD, &send_requests[UP]);
            MPI_Irecv(&((*program->grid)[start_row - 1][0]), ncols, MPI_INT, previous_process, UP, MPI_COMM_WORLD, &recv_requests[UP]);

            // exchange data with process below
            int first_process = 0;
            MPI_Isend(&((*program->grid)[end_row][0]), ncols, MPI_INT, first_process, UP, MPI_COMM_WORLD, &send_requests[DOWN]);
            MPI_Irecv(&((*program->grid)[0][0]), ncols, MPI_INT, first_process, DOWN, MPI_COMM_WORLD, &recv_requests[DOWN]);
        }

        MPI_Barrier(MPI_COMM_WORLD);
    }

    for (int row = start_row; row <= end_row; ++row) {
        fprintf(stdout, "%d: ", program->process_id);
        for (int col = 0; col < ncols; ++col) {
            fprintf(stdout, "%d", program->grid->at(row).at(col));
        }
        fprintf(stdout, "%c", '\n');
    }
}
//        if (is_first_process) {
//            auto end_row_address_current = &((*program->grid)[end_row][0]);
//            auto end_row_address_previous = &((*program->grid)[nrows - 1][0]);
//            MPI_Isend(end_row_address_current, ncols, MPI_INT, (program->process_id + 1) % program->processes, DOWN, MPI_COMM_WORLD, &send_requests[DOWN]);
//            MPI_Irecv(end_row_address_previous, ncols, MPI_INT, program->processes - 1, DOWN, MPI_COMM_WORLD, &recv_requests[DOWN]);
//
//            auto start_row_address_current = &((*program->grid)[start_row][0]);
//            auto index = (program->process_id + 1) * rows_per_process;
//            auto start_row_address_next = &((*program->grid)[index][0]);
//            MPI_Isend(start_row_address_current, ncols, MPI_INT, program->processes - 1, UP, MPI_COMM_WORLD, &send_requests[UP]);
//            MPI_Irecv(start_row_address_next, ncols, MPI_INT, (program->process_id + 1) % program->processes, UP, MPI_COMM_WORLD, &recv_requests[UP]);
//        } else if (is_middle_process) {
//            auto start_row_address_current = &((*program->grid)[start_row][0]);
//            auto start_row_address_previous = &((*program->grid)[start_row - 1][0]);
//            MPI_Isend(start_row_address_current, ncols, MPI_INT, 0, UP, MPI_COMM_WORLD, &send_requests[UP]);
//            MPI_Irecv(start_row_address_previous, ncols, MPI_INT, program->process_id - 1, UP, MPI_COMM_WORLD, &recv_requests[UP]);
//
//            auto end_row_address_current = &((*program->grid)[end_row][0]);
//            auto end_row_address_next = &((*program->grid)[0][0]);
//            MPI_Isend(end_row_address_current, ncols, MPI_INT, 0, DOWN, MPI_COMM_WORLD, &send_requests[DOWN]);
//            MPI_Irecv(end_row_address_next, ncols, MPI_INT, program->process_id + 1, DOWN, MPI_COMM_WORLD, &recv_requests[DOWN]);
//        } else if (is_last_process) {
//            auto start_row_address_current = &((*program->grid)[start_row][0]);
//            auto index = (program->process_id - 1) * rows_per_process;
//            auto start_row_address_previous = &((*program->grid)[index][0]);
//            // send first to prev process
//            MPI_Isend(start_row_address_current, ncols, MPI_INT, program->process_id - 1, UP, MPI_COMM_WORLD, &send_requests[UP]);
//            MPI_Irecv(start_row_address_previous, ncols, MPI_INT, program->process_id - 1, UP, MPI_COMM_WORLD, &recv_requests[UP]);
//
//            auto end_row_address_current = &((*program->grid)[end_row][0]);
//            auto end_row_address_next = &((*program->grid)[0][0]);
//            // send last to the 1st process
//            MPI_Isend(end_row_address_current, ncols, MPI_INT, 0, DOWN, MPI_COMM_WORLD, &send_requests[DOWN]);
//            MPI_Irecv(end_row_address_next, ncols, MPI_INT, 0, DOWN, MPI_COMM_WORLD, &recv_requests[DOWN]);
//        }


int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    auto program = new Program<void>(argc, argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &program->process_id);
    MPI_Comm_size(MPI_COMM_WORLD, &program->processes);

//    DEBUG_PRINT_LITE_PROCESS(0, "process_id: %d\n", program->process_id);
//    DEBUG_PRINT_LITE_PROCESS(0, "filename: %s\n", program->filename);
//    DEBUG_PRINT_LITE_PROCESS(0, "steps: %d\n", program->steps);
//    DEBUG_PRINT_LITE_PROCESS(0, "grid (rows x columns): %ld x %ld\n", program->grid->size(), program->grid->at(0).size());
    // print whole grid
//    if (program->process_id == 0) {
//        print_grid(program);
//    }
//    DEBUG_PRINT_LITE_PROCESS(0, "processes: %d\n", program->processes);
//    DEBUG_PRINT_LITE_PROCESS(0, "====================================%c", '\n');
    MPI_Barrier(MPI_COMM_WORLD);
    run_simulation(program);
    MPI_Barrier(MPI_COMM_WORLD);
//    if (program->process_id == 0) {
//        print_grid(program);
//    }
    MPI_Finalize();
    delete program;
    return 0;
}

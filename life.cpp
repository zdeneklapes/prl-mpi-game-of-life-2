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
        filename = argv[1];

        // Read file
        std::ifstream file_stream(filename);
        if (!file_stream.good()) {
            std::cerr << "File not found: " << filename << std::endl;
            die();
        } else {
            // Parse input
            std::string line;
            grid = new std::vector<std::vector<int>>();
            while (getline(file_stream, line)) {
                std::vector<int> row;
                for (char cell : line) {
                    if (cell == '0') {
                        row.push_back(0); // Dead cell
                    } else if (cell == '1') {
                        row.push_back(1); // Alive cell
                    }
                }
                grid->push_back(row);
            }
            file_stream.close();
        }

        // Parse steps
        try {
            steps = std::stoi(argv[2]);
            if (steps <= 0) {
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
    std::string input = std::string();
    int process_id = -1;
    int processes = -1;
    std::vector<std::vector<int>> *grid;
};

void run_simulation(Program<void> *program) {
    int nrows = program->grid->size();
    int ncols = program->grid->at(0).size();
    if (program->process_id == 0) {
        DEBUG_PRINT_LITE("nrows: %d\n", nrows);
        DEBUG_PRINT_LITE("ncols: %d\n", ncols);
    }
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    auto program = new Program<void>(argc, argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &program->process_id);
    MPI_Comm_size(MPI_COMM_WORLD, &program->processes);

    if (program->process_id == 0) {
        DEBUG_PRINT_LITE("filename: %s\n", program->filename);
        DEBUG_PRINT_LITE("steps: %d\n", program->steps);
        DEBUG_PRINT_LITE("input:\n%s\n", program->input.c_str());
        DEBUG_PRINT_LITE("process_id: %d\n", program->process_id);
        DEBUG_PRINT_LITE("processes: %d\n", program->processes);
    }

    run_simulation(program);

    MPI_Finalize();
    delete program;
    return 0;
}

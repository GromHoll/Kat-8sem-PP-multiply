#include<iostream>
#include<fstream>
#include<string>
#include<mpi.h>

using namespace std;

struct matrix {
    int x, y;
    double ** m;
};

matrix loadFromFile(string filename) {
    matrix m;

    ifstream in(filename.c_str());

    in >> m.x;
    in >> m.y;

    m.m = new double*[m.x];
    for(int i = 0; i < m.x; i++) {
        m.m[i] = new double[m.y];
    }

    for(int j = 0; j < m.y; j++) {
        for(int i = 0; i < m.x; i++) {
            in >> m.m[i][j];
        }
    }

    in.close();

    return m;
}

void writeToFile(string filename, matrix m) {
    ofstream out(filename.c_str());

    out << m.x << ' ' << m.y << endl;
    for(int j = 0; j < m.y; j++) {
        for(int i = 0; i < m.x; i++) {
            out << m.m[i][j] << "\t";
        }
        out << endl;
    }

    out.close();
}



int main(int argc, char *argv[]) {
    int myRank;
    int nodeNum;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nodeNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

    if(myRank == 0) {
        matrix A = loadFromFile("A.txt");
        matrix B = loadFromFile("B.txt");

        if(A.y == B.x) {
            matrix C;

            C.x = B.x;
            C.y = A.y;
            C.m = new double*[C.x];
            for(int i = 0; i < C.x; i++) {
                C.m[i] = new double[C.y];
            }

            MPI_Bcast(&C.x, 1, MPI_INT, 0, MPI_COMM_WORLD);

            for(int i = 0; i < C.x; i++) {
                MPI_Send(B.m[i], C.y, MPI::DOUBLE, i+1, 0, MPI_COMM_WORLD);
            }

            for(int i = 0; i < C.x; i++) {
                double row[C.x];
                for(int j = 0; j < C.x; j++) {
                    row[j] = A.m[j][i];
                }
                for(int j = 0; j < C.x; j++) {
                    MPI_Send(row, C.x, MPI::DOUBLE, j+1, 0, MPI_COMM_WORLD);
                }
                for(int j = 0; j < C.x; j++) {
                    MPI_Recv(&C.m[j][i], 1, MPI::DOUBLE, j+1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }

            writeToFile("C.txt", C);
        } else {
            cerr << "Invalid matrix dimention" << endl;
        }
    } else {
        double col[32];
        int size;
        MPI_Bcast(&size, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Recv(col, size, MPI::DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        
        for(int i = 0; i < size; i++) {
            double row[32];
            double res = 0;
            MPI_Recv(row, size, MPI::DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for(int j = 0; j < size; j++) {
                res += col[j]*row[j];
            }
            MPI_Send(&res, 1, MPI::DOUBLE, 0, 0, MPI_COMM_WORLD);
        }
	
    }

    MPI_Finalize();
    return 0;
}

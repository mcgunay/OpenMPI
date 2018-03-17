#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include <stdlib.h>

#define mpi_comm_world MPI_COMM_WORLD
const int MAX_STRING = 100;
const int MAX_NUMBERS = 100;
#define LMAX 255


// Define a vector type
typedef struct {
    int *array;
    size_t size;
    size_t capacity;
} Array;

int * intdup(int const*, size_t, int);
void printArray (double *, int, int );
//void printMatrix (double *, int , int , int);
void printMatrix ( double*, int,int, int );
void printResultArray (double *, int );


int main(void) {

    int com_sz;
    int my_rank;
    int my_sum = 0;
    int total_sum = 0;
    double ** receive_buffer;
    double **mat_buffer;
    double * matrix_buffer;
    double * vec_buffer;
    int* row_size, *column_size;
    int col;
    double** distribution_matrix;
    double *distribution_vec;
    Array myData;
    column_size = (int *) malloc(sizeof(int)); //For Calculating number size per core
    row_size = (int *) malloc(sizeof(int)); //For Calculating number size per core

    MPI_Init(NULL, NULL);
    MPI_Comm_size(mpi_comm_world, &com_sz);
    MPI_Comm_rank(mpi_comm_world, &my_rank);
    double *matrix;


    if(my_rank == 0) {

        printf("Total Number of Processor %d\n", com_sz);

        char **array = NULL;
        char *ln = NULL;
        size_t n = 0;
        ssize_t nchr = 0;
        size_t idx = 0;
        size_t it = 0;
        size_t lmax = LMAX;
        FILE *fp = NULL;
        if (!(fp = fopen ("../matrix.out", "r"))) { /
            fprintf (stderr, "error: file open failed '%s'.", "../matrix.out");
            return 1;
        }

        if (!(array = calloc (LMAX, sizeof *array))) {
            fprintf (stderr, "error: memory allocation failed.");
            return 1;
        }

        if((nchr = getline (&ln, &n, fp)) != -1){
            char num;
            int count = 0;

            while(sscanf(ln, "%c", &num) != -1) {
                if(num == ' ')
                    count++;
                ++ln;
            }

            *row_size = count ;
            *column_size = count;
            printf("column size = %d\n", *column_size);
        }

        if (fp) fclose (fp);


        double (*arr)[*column_size] = malloc(sizeof *arr * *row_size);

        FILE *file;
        file=fopen("../matrix.out", "r");

        for(int i = 0; i < *row_size; i++)
        {
            for(int j = 0; j < *column_size; j++)
            {
                if (!fscanf(file, "%lf", &arr[i][j]))
                    break;
                // mat[i][j] -= '0';
                //printf("%lf\n",arr[i][j]); //Use lf format specifier, \n is for new line
            }

        }
        fclose(file);

        double* vec = malloc(*column_size* sizeof(double));

        FILE *file_vec;
        file_vec = fopen("../vector.out", "r");

        for (int k = 0; k < *column_size ; ++k) {
            if(!fscanf(file_vec, "%lf", &vec[k]))
                break;

            //printf("%lf\n", vec[k]);
        }


        fclose(file_vec);


        distribution_matrix=malloc(*row_size*sizeof(double*));
        for(col=0;col<*column_size;++col)
            distribution_matrix[col]=malloc(*column_size*sizeof(double));

        distribution_vec = malloc(sizeof(double) * *row_size);


        for (int i1 = 0; i1 < com_sz; ++i1) {
            for (int i = 0,k = 0; k < *column_size/com_sz; ++k,i+=com_sz) {
                distribution_matrix[(i1**column_size/com_sz) +k] = arr[i1 + i];
            }
        }

        if (( matrix = malloc(*column_size**row_size*sizeof(double))) == NULL) {
            printf("Malloc error");
            exit(1);
        }

        for (int i = 0; i < *row_size; ++i) {
            for (int j = 0; j < *column_size; ++j) {
                matrix[i * (*row_size) + j] = distribution_matrix[i][j];
            }
        }


        for (int i1 = 0; i1 < com_sz; ++i1) {
            for (int i = 0,k = 0; k < *column_size/com_sz; ++k,i+=com_sz) {
                distribution_vec[(i1**column_size/com_sz) +k] = vec[i1 + i];
            }
        }

        //printMatrix(matrix, *column_size, *column_size, my_rank);
        printf("\n");
        printf("\n");
        printf("\n");
        //printArray(distribution_vec,*column_size,my_rank);


    }
    double starttime, endtime;

    if(0 == my_rank)
        starttime = MPI_Wtime();

    MPI_Bcast(column_size, 1, MPI_INT, 0, MPI_COMM_WORLD); // Core 0 broadcast number size so that other processes can use

    size_t rows, cols;
    rows = *column_size/com_sz;
    cols = *column_size;
    //double (*temp)[cols] = malloc(sizeof *temp * (rows));
    double* temp = malloc(sizeof(double) * *column_size * (*column_size/com_sz));
    MPI_Scatter(matrix, *column_size*(*column_size/com_sz), MPI_DOUBLE, // send one row, which contains p integers
                temp, *column_size*(*column_size/com_sz), MPI_DOUBLE, // receive one row, which contains p integers
                0, MPI_COMM_WORLD);

    vec_buffer = malloc(sizeof(double) * (*column_size/com_sz));
    MPI_Scatter(distribution_vec, *column_size/com_sz, MPI_DOUBLE, vec_buffer, *column_size/com_sz, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    //Elimination
    double* vec_bcast_buffer = malloc(sizeof(double));
    int divided_row_count = *column_size/com_sz - 1 ;


    for (int i1 = *column_size - 1; i1 >= 0 ; --i1) {

        if (my_rank == (i1 % com_sz)) {     //That means division for this process

            memcpy(vec_bcast_buffer, &vec_buffer[divided_row_count], sizeof(double));
            divided_row_count = divided_row_count - 1;

        }
        int root = i1 % com_sz;

        MPI_Bcast(vec_bcast_buffer, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);

        //Elimination
        for (int i = divided_row_count; i >= 0; --i) {//Each undivided row will be eliminatedi
            //Eliminate Vector
            vec_buffer[i] =
                    vec_buffer[i] - (temp[(i * *column_size) + i1] * *vec_bcast_buffer);

        }

    }


    //printMatrix(temp, *column_size/com_sz, *column_size, my_rank);
    //printArray(vec_buffer, *column_size/com_sz, my_rank);
    double* gatter_vector;

    if (0 == my_rank) {
        gatter_vector = malloc(sizeof(double) **column_size);
    }


    //Gatter Vector
    MPI_Gather(vec_buffer, *column_size / com_sz, MPI_DOUBLE, gatter_vector, *column_size / com_sz, MPI_DOUBLE, 0,
               MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    if(0 == my_rank){
        endtime   = MPI_Wtime();
        printf("That took %f seconds\n",endtime-starttime);
    }

    //Write to out files
    if(0 == my_rank){

        //Reorder the matrix

        double* ordered_vector = malloc(sizeof(double) **column_size);
        int number_of_row = *column_size / com_sz;


        int col = *column_size;
        int x = 0;
        for (int j = 0; j < *column_size - 1; ++j) {
            if(j % number_of_row == 0)
                x = j / number_of_row;
            else
                x = x + com_sz;

            ordered_vector[x] = gatter_vector[j];

        }

        ordered_vector[*column_size -1] = gatter_vector[*column_size - 1];



        FILE *fp_vec;
        char output_vec[]="vector_result.out";
        int n;

        fp_vec=fopen(output_vec,"wb");

        for(n=0;n<*column_size;++n) {

            fprintf(fp_vec,"%lf\n",ordered_vector[n]);

        }

        fclose(fp_vec);


    }

    MPI_Finalize();
    return 0;
}

// Functions for dynamic array
int * intdup(int const * src, size_t len, int startIndex)
{
    int * p = malloc(len * sizeof(int));
    memcpy(p, &src[startIndex], len * sizeof(int));
    return p;
}


void printArray (double *row, int nElements, int rank) {

    int i;
    for (i=0; i<nElements; ++i) {
        printf("process %d gets data %lf ", rank, row[i]);
    }
    printf("\n");
}

void printMatrix (double* mat, int rows,int col, int rank) {

    int i;
    for (i=0; i<rows; i++) {
        printArray(&mat[i*col], col, rank);
        printf("\n");
    }
}

void printResultArray (double *row, int nElements) {

    int i;
    for (i=0; i<nElements; ++i) {
        for (int j = 0; j < nElements; ++j) {
            printf("%lf ",  row[i*nElements + j]);

        }
        printf("\n");
    }
    printf("\n");
}




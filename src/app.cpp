#include <iostream>
#include <mutex>
#include <thread>
#include <cstring>
#include <vector>
#include <atomic>
#include <algorithm>
#include <chrono>
#include <sys/shm.h>
#include <sys/mman.h>
#include <sys/fcntl.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>
#include <unistd.h>
#include <string.h>

#define print(t) std::cerr << t
#define THREAD_COUNT 8
// Thread count is per field


template <typename T>
class ScalarField {
private:
    T* array_field;

    bool threadable;
    T* thread_source;

    void thread_cpt(int start_row, int row_count, bool thread_op_add) {
        if (thread_op_add)
            for(int i = start_row; i < (start_row + row_count); i++) array_field[i] += thread_source[i];
        else
            for(int i = start_row; i < (start_row + row_count); i++) array_field[i] -= thread_source[i];
    }

    void thread_compute(T* source, bool thread_op_add) {
        std::vector<std::thread> thread_list;

        int rows_per_thread = this->size_x / THREAD_COUNT;
        for (int i = 0; i < THREAD_COUNT - 1; i++) {
            int start_row = i*rows_per_thread;
            thread_list.push_back(std::thread([this, start_row, rows_per_thread, thread_op_add] {
                this->thread_cpt(start_row, rows_per_thread, thread_op_add);
            }));
        }
        int start_row = (THREAD_COUNT - 1)*rows_per_thread;
        rows_per_thread += this->size_x % THREAD_COUNT;
        thread_list.push_back(std::thread([this, start_row, rows_per_thread, thread_op_add] {
                this->thread_cpt(start_row, rows_per_thread, thread_op_add);
        }));

        std::for_each(thread_list.begin(), thread_list.end(), [](std::thread &t) { t.join(); });
    }

public:
    int size_x;
    int size_y;
    int size_z;

    ScalarField(int size_x, int size_y, int size_z, bool threadable = false) {
        this->size_x = size_x;
        this->size_y = size_y;
        this->size_z = size_z;
        this->threadable = threadable;
        this->array_field = new T[this->size()];
        memset(this->array_field, 0, this->size() * sizeof(T));
    }

    ~ScalarField() { 
        delete[] array_field;
    }

    T* arr() { return array_field; }

    ScalarField<T>* set_arr(T* source) { 
        delete[] array_field; 
        array_field = source; 
        return this; 
    }

    void debug_print() {
        for (int z = 0; z < size_z; z++) {
            std::cout << "[";
            for (int y = 0; y < size_y; y++) {
                std::cout << "    [";
                for (int x = 0; x < size_y; x++) {
                    std::cout << array_field[x + (size_x * (y + (size_y * z)))] << ", ";
                }
                std::cout << "],\n";
            }
            std::cout << "],\n";
        }
    }

    int size() const { return size_x * size_y * size_z; }

    ScalarField<T>* copy_from(T* source) {
        memcpy(this->array_field, source, sizeof(T)*this->size());
        return this;
    }

    static ScalarField<T>* from(T* source, int size_x, int size_y, int size_z) {
        ScalarField<T>* sf = new ScalarField<T>(size_x, size_y, size_z);
         return sf->copy_from(source);
    }

    static ScalarField<T>* from_ptr(T* source, int size_x, int size_y, int size_z) {
        ScalarField<T>* sf = new ScalarField<T>(size_x, size_y, size_z);
        return sf->set_arr(source);
    }

    ScalarField<T>* copy() const {
        return ScalarField<T>::from(array_field, size_x, size_y, size_z);
    }

    ScalarField<T>* slice(int x_from, int x_to, int y_from, int y_to, int z_from, int z_to) const {
        
        // This is kinda ugly... but idc
        // It allows negative indices, like in Python

        if (x_from < 0) return slice(this->size_x + x_from, x_to, y_from, y_to, z_from, z_to);
        if (x_to < 0) return slice(x_from, this->size_x + x_to, y_from, y_to, z_from, z_to);
        if (y_from < 0) return slice(x_from, x_to, this->size_y + y_from, y_to, z_from, z_to);
        if (y_to < 0) return slice(x_from, x_to, y_from, this->size_y + y_to, z_from, z_to);
        if (z_from < 0) return slice(x_from, x_to, y_from, y_to, this->size_z + z_from, z_to);
        if (z_to < 0) return slice(x_from, x_to, y_from, y_to, z_from, this->size_z + z_to);

        int new_size_x = (x_to - x_from) + 1;
        int new_size_y = (y_to - y_from) + 1;
        int new_size_z = (z_to - z_from) + 1;

        T* new_array = new T[new_size_x * new_size_y * new_size_z];
        size_t sx = sizeof(T)*new_size_x;

        for (int z = z_from; z <= z_to; z++) {
            for (int y = y_from; y <= y_to; y++) {
                int idx_dest = new_size_x * ((y - y_from) + (new_size_y * (z - z_from)));
                int idx_src = x_from + (size_x * (y + (size_y * z)));
                memcpy(&new_array[idx_dest], &array_field[idx_src], sx);
            }
        }

        ScalarField<T>* res = ScalarField<T>::from_ptr(new_array, new_size_x, new_size_y, new_size_z);
        return res;
    }

    ScalarField<T>* x(int y_index, int z_index) const {
        return slice(0, -1, y_index, y_index, z_index, z_index);        
    }
    ScalarField<T>* y(int y_index, int z_index) const {
        return slice(0, -1, y_index, y_index, z_index, z_index);        
    }
    ScalarField<T>* z(int y_index, int z_index) const {
        return slice(0, -1, y_index, y_index, z_index, z_index);
    }

    ScalarField<T>* prepend_x(T default_value) {
        T* new_array = new T[(size_x + 1) * size_y * size_z];
        size_t sx = sizeof(T)*size_x;

        for (int z = 0; z < size_z; z++) {
            for (int y = 0; y < size_y; y++) {
                int idx_dest = (size_x + 1) * (y + (size_y * z));
                int idx_src = size_x * (y + (size_y * z));
                memcpy(&new_array[idx_dest+1], &array_field[idx_src], sx);
                new_array[idx_dest] = default_value;
            }
        }

        delete[] array_field;
        array_field = new_array;
        size_x += 1;
        return this;
    }

    ScalarField<T>* prepend_y(T default_value) {
        T* new_array = new T[size_x * (size_y + 1) * size_z];
        size_t sx = sizeof(T)*size_x;

        for (int z = 0; z < size_z; z++) {
            memset(&new_array[(size_y + 1) * z], default_value, sx);
            for (int y = 1; y < size_y + 1; y++) {
                int idx_dest = size_x * (y + ((size_y + 1) * z));
                int idx_src = size_x * (y + (size_y * z));
                memcpy(&new_array[idx_dest], &array_field[idx_src], sx);
            }
        }

        delete[] array_field;
        array_field = new_array;
        size_y += 1;
        return this;
    }

    ScalarField<T>* prepend_z(T default_value) {
        T* new_array = new T[size_x * size_y * (size_z + 1)];
        size_t sx = sizeof(T)*size_x;

        memset(&new_array[0], default_value, sx*size_y);
        for (int z = 1; z < size_z + 1; z++) {
            for (int y = 0; y < size_y; y++) {
                int idx_dest = size_x * (y + (size_y * z));
                int idx_src = idx_dest - (size_y * size_x);
                memcpy(&new_array[idx_dest], &array_field[idx_src], sx);
            }
        }

        delete[] array_field;
        array_field = new_array;
        size_z += 1;
        return this;
    }

    ScalarField<T>* dx() const {
        ScalarField<T>* sl = slice(1, -1, 0, -1, 0, -1);
        ScalarField<T>* sl_t = slice(0, -2, 0, -1, 0, -1);
        sl->operator-=(*sl_t);
        sl->prepend_x(0);
        delete sl_t;
        return sl;
        //return &((*slice(1, -1, 0, -1, 0, -1)) - (*slice(0, -2, 0, -1, 0, -1)));
    }
    ScalarField<T>* dy() const {
        ScalarField<T>* sl = slice(0, -1, 1, -1, 0, -1);
        ScalarField<T>* sl_t = slice(0, -1, 0, -2, 0, -1);
        sl->operator-=(*sl_t);
        sl->prepend_y(0);
        delete sl_t;
        return sl;
        //return &((*slice(0, -1, 1, -1, 0, -1)) - (*slice(0, -1, 0, -2, 0, -1)));
    }
    ScalarField<T>* dz() const {
        ScalarField<T>* sl = slice(0, -1, 0, -1, 1, -1);
        ScalarField<T>* sl_t = slice(0, -1, 0, -1, 0, -2);
        sl->operator-=(*sl_t);
        sl->prepend_z(0);
        delete sl_t;
        return sl;
        //return &((*slice(0, -1, 0, -1, 1, -1)) - (*slice(0, -1, 0, -1, 0, -2)));
    }

    // FIXME something doesn't work in + and - operators... f seems to get destroyed when leaving this context.
    ScalarField<T>& operator+ (const ScalarField<T>& s) const {
        ScalarField<T>* f = this->copy();
        for(int i = 0; i < this->size(); i++) (*f)[i] += s[i];
        return *f;
    }

    ScalarField<T>& operator- (const ScalarField<T>& s) const {
        ScalarField<T>* f = this->copy();
        for(int i = 0; i < this->size(); i++) (*f)[i] -= s[i];
        return *f;
    }

    ScalarField<T>& operator+= (const ScalarField<T>& s) {
        if(threadable) 
            thread_compute(s.array_field, true);
        else {
            int sz = this->size();
            for(int i = 0; i < sz; i++) {
                array_field[i] += s[i];
            }
        }
        return *this;
    }

    ScalarField<T>& operator-= (const ScalarField<T>& s) {
        if(threadable) 
            thread_compute(s.array_field, false);
        else {
            int sz = this->size();
            for(int i = 0; i < sz; i++) array_field[i] -= s[i];
        }
        return *this;
    }

    T operator [](int i) const { return array_field[i]; } // Getter
    T& operator [](int i) { return array_field[i]; } // Setter
};

template <typename T>
class VectorField {
private:
    bool threadable;

public:
    ScalarField<T>* i;
    ScalarField<T>* j;
    ScalarField<T>* k;
    int size_x;
    int size_y;
    int size_z;

    VectorField(int size_x, int size_y, int size_z, bool threadable = false) {
        this->size_x = size_x;
        this->size_y = size_y;
        this->size_z = size_z;
        this->threadable = threadable;
        this->i = new ScalarField<T>(size_x, size_y, size_z, threadable);
        this->j = new ScalarField<T>(size_x, size_y, size_z, threadable);
        this->k = new ScalarField<T>(size_x, size_y, size_z, threadable);
    }

    ~VectorField() {
        delete i;
        delete j;
        delete k;
    }

    int size() { return size_x * size_y * size_z; }

    VectorField<T>* curl() {
        VectorField<T>* vf = this->copy();

        ScalarField<T>* kdy = vf->k->dy();
        ScalarField<T>* jdz = vf->j->dz();
        ScalarField<T>* idz = vf->i->dz();
        ScalarField<T>* kdx = vf->k->dx();
        ScalarField<T>* jdx = vf->j->dx();
        ScalarField<T>* idy = vf->i->dy();

        vf->i->operator+=(*kdy);
        vf->i->operator-=(*jdz);
        vf->j->operator+=(*idz);
        vf->j->operator-=(*kdx);
        vf->k->operator+=(*jdx);
        vf->k->operator-=(*idy);

        delete kdy;
        delete jdz;
        delete idz;
        delete kdx;
        delete jdx;
        delete idy;

        return vf;
    }

    VectorField<T>* copy() {
        VectorField<T>* vf = new VectorField<T>(size_x, size_y, size_z, threadable);
        vf->i = this->i->copy();
        vf->j = this->j->copy();
        vf->k = this->k->copy();
        return vf;
    }

};

///// imported
// Taille de la matrice de travail (un côté)
static const int MATRIX_SIZE = 100;
static const int BUFFER_SIZE = MATRIX_SIZE * MATRIX_SIZE * MATRIX_SIZE * sizeof(double) * 3 * 2;
// Tampon générique à utiliser pour créer le fichier
char buffer_[BUFFER_SIZE];

void wait_signal()
{
    // Attend une entrée (ligne complète avec \n) sur stdin.
    std::string msg;
    std::cin >> msg;
    std::cerr << "CPP: Got signal." << std::endl;
}

void ack_signal()
{
    // Répond avec un message vide.
    std::cout << "" << std::endl;
}
///// imported

void churros()
{
    int n = MATRIX_SIZE;
    // VectorField<double>* v = new VectorField<double>(n, n, n, true);

    // auto start = std::chrono::high_resolution_clock::now();
    // VectorField<double>* c = v->curl();
    // auto stop = std::chrono::high_resolution_clock::now();
    // auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

    // std::cout << "Time taken with threading: " << duration.count() << " microseconds" << std::endl;

    // VectorField<double>* v2 = new VectorField<double>(n, n, n, false);

    // start = std::chrono::high_resolution_clock::now();
    // VectorField<double>* c2 = v->curl();
    // stop = std::chrono::high_resolution_clock::now();
    // duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

    // std::cout << "Time taken without threading: " << duration.count() << " microseconds" << std::endl;

}

int main(int argc, char** argv) {
    // print("Hello World\n");
    int n = 100;
    //  pogne matrice

   if (argc < 2)
    {
        std::cerr << "Error : no shared file provided as argv[1]" << std::endl;
        return -1;
    }

    wait_signal();

    // Création d'un fichier "vide" (le fichier doit exister et être d'une
    // taille suffisante avant d'utiliser mmap).
    memset(buffer_, 0, BUFFER_SIZE);
    FILE* shm_f = fopen(argv[1], "w");
    fwrite(buffer_, sizeof(char), BUFFER_SIZE, shm_f);
    fclose(shm_f);

    

    // On signale que le fichier est prêt.
    std::cerr << "CPP:  File ready." << std::endl;
    ack_signal();

    // On ré-ouvre le fichier et le passe à mmap(...). Le fichier peut ensuite
    // être fermé sans problèmes (mmap y a toujours accès, jusqu'à munmap.)
    int shm_fd = open(argv[1], O_RDWR);
    void* shm_mmap = mmap(NULL, BUFFER_SIZE, PROT_WRITE | PROT_READ, MAP_SHARED, shm_fd, 0);
    close(shm_fd);

    if (shm_mmap == MAP_FAILED) {
        std::cerr << "ERROR SHM\n";
        perror(NULL);
        return -1;
    }

    // Pointeur format double qui représente la matrice partagée:
    double* mtx = (double*)shm_mmap;

    wait_signal();
    double courant_number = mtx[0];
    double source_pos_x = mtx[1];
    double source_pos_y = mtx[2];
    double source_pos_z = mtx[3];
    double source_pos_m = mtx[4];
    double source_val = mtx[5];

    ack_signal();

    // TODO SET A ZERO
    VectorField<double> *E = new VectorField<double>(MATRIX_SIZE,MATRIX_SIZE,MATRIX_SIZE);
    VectorField<double> *F = new VectorField<double>(MATRIX_SIZE,MATRIX_SIZE,MATRIX_SIZE);


    int i = 0;
    while (true)
    {
        wait_signal();
        
        E->i->operator[](i) = 100.0;
        // E += courant_number * F->curl();
        // F -= courant_number * E->curl();

        memcpy(&mtx[0], E->i->arr(), E->i->size()*sizeof(double));
        memcpy(&mtx[E->i->size()], E->j->arr(), E->j->size()*sizeof(double));
        memcpy(&mtx[E->i->size()+E->j->size()], E->k->arr(), E->k->size()*sizeof(double));

        std::cerr << "CPP: Work done." << std::endl;
        ack_signal();
        i++;
    }
    

    //c->i->debug_print();
    munmap(shm_mmap, BUFFER_SIZE);
    return 0;
}

/*

To take into consideration (VERY IMPORTANT):

When splitting the Vector Field up in pieces as to multithread / multiprocess, regions in common will not work as intended.

For example:

[
    1, 2, 3, 4, 5, 6, 7, 8,
    1, 2, 3, 4, 5, 6, 7, 8,
    1, 2, 3, 4, 5, 6, 7, 8,
    1, 2, 3, 4, 5, 6, 7, 8,
    1, 2, 3, 4, 5, 6, 7, 8,
    1, 2, 3, 4, 5, 6, 7, 8
]

Were I to split this array in 2, vertically, and then calculate the curl, whenever dy would be called, some information would be missing...


dy of [
    1, 2, 3, 4, 5, 6, 7, 8,
    1, 2, 3, 4, 5, 6, 7, 8,
    1, 2, 3, 4, 5, 6, 7, 8,
]
= [
    1, 2, 3, 4, 5, 6, 7, 8,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
]

Therefore, combining them into the final image would give:

[
    1, 2, 3, 4, 5, 6, 7, 8,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    1, 2, 3, 4, 5, 6, 7, 8,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
]

where the 4th column should be filled with zeroes...
To paliate this problem, we must give one more row on each side to each thread.

T1 = T2 =
[
    1, 2, 3, 4, 5, 6, 7, 8,
    1, 2, 3, 4, 5, 6, 7, 8,
    1, 2, 3, 4, 5, 6, 7, 8,
    1, 2, 3, 4, 5, 6, 7, 8,
]
dy of T1 = dy of T2 = [
    1, 2, 3, 4, 5, 6, 7, 8,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
]

We must then combine them back together by keeping the 4th row of T1 instead of T2, as to give:

[
    1, 2, 3, 4, 5, 6, 7, 8,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
]

*/
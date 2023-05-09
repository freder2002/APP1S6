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

    void thread_cpt(int start_row, int row_count, int thread_op) {
        if (thread_op == 0)
            for(int i = start_row; i < (start_row + row_count); i++) array_field[i] += thread_source[i];
        else if (thread_op == 1)
            for(int i = start_row; i < (start_row + row_count); i++) array_field[i] -= thread_source[i];
        else if (thread_op == 2)
            for(int i = start_row; i < (start_row + row_count); i++) array_field[i] *= thread_source[i];
    
    }

    void thread_compute(T* source, int thread_op) {
        std::vector<std::thread> thread_list;

        delete[] thread_source;
        thread_source = source;

        int rows_per_thread = this->size_x / THREAD_COUNT;
        for (int i = 0; i < THREAD_COUNT - 1; i++) {
            int start_row = i*rows_per_thread;
            thread_list.push_back(std::thread([this, start_row, rows_per_thread, thread_op] {
                this->thread_cpt(start_row, rows_per_thread, thread_op);
            }));
        }
        int start_row = (THREAD_COUNT - 1)*rows_per_thread;
        rows_per_thread += this->size_x % THREAD_COUNT;
        thread_list.push_back(std::thread([this, start_row, rows_per_thread, thread_op] {
                this->thread_cpt(start_row, rows_per_thread, thread_op);
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
        empty();
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

    void empty() {
        memset(this->array_field, 0, this->size() * sizeof(T));
    }

    void debug_print() {
        for (int z = 0; z < size_z; z++) {
            std::cerr << "[";
            for (int y = 0; y < size_y; y++) {
                std::cerr << "    [";
                for (int x = 0; x < size_y; x++) {
                    std::cerr << array_field[x + (size_x * (y + (size_y * z)))] << ", ";
                }
                std::cerr << "],\n";
            }
            std::cerr << "],\n";
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

    ScalarField<T>* prepend_x(T default_value, bool forward = true) {
        T* new_array = new T[(size_x + 1) * size_y * size_z];
        size_t sx = sizeof(T)*size_x;

        if (forward) {
            for (int z = 0; z < size_z; z++) {
                for (int y = 0; y < size_y; y++) {
                    int idx_dest = (size_x + 1) * (y + (size_y * z));
                    int idx_src = size_x * (y + (size_y * z));
                    memcpy(&new_array[idx_dest+1], &array_field[idx_src], sx);
                    new_array[idx_dest] = default_value;
                }
            }
        } else {
            for (int z = 0; z < size_z; z++) {
                for (int y = 0; y < size_y; y++) {
                    int idx_dest = (size_x + 1) * (y + (size_y * z));
                    int idx_src = size_x * (y + (size_y * z));
                    memcpy(&new_array[idx_dest], &array_field[idx_src], sx);
                    new_array[idx_dest + size_x] = default_value;
                }
            }
        }

        delete[] array_field;
        array_field = new_array;
        size_x += 1;
        return this;
    }

    ScalarField<T>* prepend_y(T default_value, bool forward = true) {
        T* new_array = new T[size_x * (size_y + 1) * size_z];
        size_t sx = sizeof(T)*size_x;

        if (forward) {
            for (int z = 0; z < size_z; z++) {
                memset(&new_array[(size_y + 1) * z], default_value, sx);
                for (int y = 1; y < size_y + 1; y++) {
                    int idx_dest = size_x * (y + ((size_y + 1) * z));
                    int idx_src = size_x * (y + (size_y * z));
                    memcpy(&new_array[idx_dest], &array_field[idx_src], sx);
                }
            }
        } else {
            for (int z = 0; z < size_z; z++) {
                for (int y = 0; y < size_y; y++) {
                    int idx_dest = size_x * (y + ((size_y + 1) * z));
                    int idx_src = size_x * (y + (size_y * z));
                    memcpy(&new_array[idx_dest], &array_field[idx_src], sx);
                }
                memset(&new_array[size_x * (size_y + ((size_y + 1) * z))], default_value, sx);
            }
        }

        delete[] array_field;
        array_field = new_array;
        size_y += 1;
        return this;
    }

    ScalarField<T>* prepend_z(T default_value, bool forward = true) {
        T* new_array = new T[size_x * size_y * (size_z + 1)];
        size_t sx = sizeof(T)*size_x;

        if(forward) {
            memset(&new_array[0], default_value, sx*size_y);
            for (int z = 1; z < size_z + 1; z++) {
                for (int y = 0; y < size_y; y++) {
                    int idx_dest = size_x * (y + (size_y * z));
                    int idx_src = idx_dest - (size_y * size_x);
                    memcpy(&new_array[idx_dest], &array_field[idx_src], sx);
                }
            }
        } else {
            for (int z = 0; z < size_z; z++) {
                for (int y = 0; y < size_y; y++) {
                    int idx_dest = size_x * (y + (size_y * z));
                    int idx_src = idx_dest - (size_y * size_x);
                    memcpy(&new_array[idx_dest], &array_field[idx_src], sx);
                }
            }
            memset(&new_array[size_z], default_value, sx*size_y);
        }

        delete[] array_field;
        array_field = new_array;
        size_z += 1;
        return this;
    }

    ScalarField<T>* dx(bool forward = true) const {
        ScalarField<T>* sl = slice(1, -1, 0, -1, 0, -1);
        ScalarField<T>* sl_t = slice(0, -2, 0, -1, 0, -1);
        sl->operator-=(*sl_t);
        sl->prepend_x(0, forward);
        delete sl_t;
        return sl;
        //return &((*slice(1, -1, 0, -1, 0, -1)) - (*slice(0, -2, 0, -1, 0, -1)));
    }
    ScalarField<T>* dy(bool forward = true) const {
        ScalarField<T>* sl = slice(0, -1, 1, -1, 0, -1);
        ScalarField<T>* sl_t = slice(0, -1, 0, -2, 0, -1);
        sl->operator-=(*sl_t);
        sl->prepend_y(0, forward);
        delete sl_t;
        return sl;
        //return &((*slice(0, -1, 1, -1, 0, -1)) - (*slice(0, -1, 0, -2, 0, -1)));
    }
    ScalarField<T>* dz(bool forward = true) const {
        ScalarField<T>* sl = slice(0, -1, 0, -1, 1, -1);
        ScalarField<T>* sl_t = slice(0, -1, 0, -1, 0, -2);
        sl->operator-=(*sl_t);
        sl->prepend_z(0, forward);
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
            thread_compute(s.array_field, 0);
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
            thread_compute(s.array_field, 1);
        else {
            int sz = this->size();
            for(int i = 0; i < sz; i++) array_field[i] -= s[i];
        }
        return *this;
    }

    ScalarField<T>& operator*= (const ScalarField<T>& s) {
        if(threadable) 
            thread_compute(s.array_field, 2);
        else {
            int sz = this->size();
            for(int i = 0; i < sz; i++) array_field[i] *= s[i];
        }
        return *this;
    }

    ScalarField<T>& operator+= (const T& s) {
        if(threadable) {
            T* source = new T[this->size()];
            memset(source, s, sizeof(T)*this->size());
            thread_compute(source, 0);
        }
        else {
            int sz = this->size();
            for(int i = 0; i < sz; i++) {
                array_field[i] += s;
            }
        }
        return *this;
    }

    ScalarField<T>& operator-= (const T& s) {
        if(threadable) {
            T* source = new T[this->size()];
            memset(source, s, sizeof(T)*this->size());
            thread_compute(source, 1);
        }
        else {
            int sz = this->size();
            for(int i = 0; i < sz; i++) {
                array_field[i] -= s;
            }
        }
        return *this;
    }

    ScalarField<T>& operator*= (const T& s) {
        if(threadable) {
            T* source = new T[this->size()];
            memset(source, s, sizeof(T)*this->size());
            thread_compute(source, 2);
        }
        else {
            int sz = this->size();
            for(int i = 0; i < sz; i++) {
                array_field[i] *= s;
            }
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

    VectorField<T>* curl(bool forward = true) {
        VectorField<T>* vf = new VectorField<T>(size_x, size_y, size_z, threadable);

        ScalarField<T>* kdy = this->k->dy(forward);
        ScalarField<T>* jdz = this->j->dz(forward);
        ScalarField<T>* idz = this->i->dz(forward);
        ScalarField<T>* kdx = this->k->dx(forward);
        ScalarField<T>* jdx = this->j->dx(forward);
        ScalarField<T>* idy = this->i->dy(forward);

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

    VectorField<T>& operator+= (const VectorField<T>& s) {
        this->i->operator+=(*s.i);
        this->j->operator+=(*s.j);
        this->k->operator+=(*s.k);
        return *this;
    }

    VectorField<T>& operator-= (const VectorField<T>& s) {
        this->i->operator-=(*s.i);
        this->j->operator-=(*s.j);
        this->k->operator-=(*s.k);
        return *this;
    }

    VectorField<T>& operator*= (const VectorField<T>& s) {
        this->i->operator*=(*s.i);
        this->j->operator*=(*s.j);
        this->k->operator*=(*s.k);
        return *this;
    }

    VectorField<T>& operator+= (const T& s) {
        this->i->operator+=(s);
        this->j->operator+=(s);
        this->k->operator+=(s);
        return *this;
    }
    
    VectorField<T>& operator-= (const T& s) {
        this->i->operator-=(s);
        this->j->operator-=(s);
        this->k->operator-=(s);
        return *this;
    }

    VectorField<T>& operator*= (const T& s) {
        this->i->operator*=(s);
        this->j->operator*=(s);
        this->k->operator*=(s);
        return *this;
    }

    void empty() {
        this->i->empty();
        this->j->empty();
        this->k->empty();
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

        courant_number = mtx[0];
        source_pos_x = mtx[1];
        source_pos_y = mtx[2];
        source_pos_z = mtx[3];
        source_pos_m = mtx[4];
        source_val = mtx[5];

        VectorField<double>* tmp_F = F->curl();
        tmp_F->operator*=(courant_number);
        E->operator+=(*tmp_F);

        int source_pos = source_pos_x + (E->size_x * (source_pos_y + (E->size_y * source_pos_z)));
        E->i->operator[](source_pos) += source_val;

        VectorField<double>* tmp_E = E->curl(false);
        tmp_E->operator*=(courant_number);
        F->operator-=(*tmp_E);
        
        delete tmp_E;
        delete tmp_F;

        // E->i->operator[](i) = 100.0;
        // E += courant_number * F->curl();
        // F -= courant_number * E->curl();
        

        memcpy(&mtx[0], E->i->arr(), E->i->size()*sizeof(double));
        memcpy(&mtx[E->i->size()], E->j->arr(), E->j->size()*sizeof(double));
        memcpy(&mtx[E->i->size()+E->j->size()], E->k->arr(), E->k->size()*sizeof(double));
        memcpy(&mtx[E->i->size()+E->j->size()+E->k->size()], F->i->arr(), F->i->size()*sizeof(double));
        memcpy(&mtx[E->i->size()+E->j->size()+E->k->size()+F->i->size()], F->j->arr(), F->j->size()*sizeof(double));
        memcpy(&mtx[E->i->size()+E->j->size()+E->k->size()+F->i->size()+F->j->size()], F->k->arr(), F->k->size()*sizeof(double));

        double* arr = E->j->arr();
        int ind = 505050;
        for (int i = 0; i < 100; i++)
            std::cerr << "[" << mtx[E->i->size()+i] << "," << arr[i] << "], ";
        std::cerr << std::endl;

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
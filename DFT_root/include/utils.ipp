#include "utils.hpp"
#include <slate/slate.hh>
#include <blas.hh>
#include <mpi.h>

/* ---- Vector class ----*/
template <typename scalar_t>
Vector<scalar_t>::Vector(int m, int nb, int p, slate::Target target) {
    vector = slate::Matrix < scalar_t >(m, 1, nb, p, 1, MPI_COMM_WORLD);
    vector.insertLocalTiles(target);
};

template <typename scalar_t>
Vector<scalar_t>::Vector(int m, int nb, int p) {
    Vector(m, nb, p, slate::HostTask);
};

template <typename scalar_t>
Vector<scalar_t>::~Vector() {
    vector.clear();
};

template <typename scalar_t>
void Vector<scalar_t>::Set(int index, scalar_t value){
    int mt = vector.mt();
    for (int64_t i = 0; i < mt; i++){
        if ( i*(mt+1) <= index && index < i*(mt+1)+vector.tileMb(i) ){
            if (vector.tileIsLocal( i, 0)) {
                auto tile = vector( i, 0);
                auto tiledata = vector.data();
                tiledata[index-i*(mt+1)] = value;
                return;
            }
        }
    }
};

template <typename scalar_t>
scalar_t Vector<scalar_t>::Get(int index){
    int mt = vector.mt();
    for (int64_t i = 0; i < mt; i++){
        if ( i*(mt+1) <= index && index < i*(mt+1)+vector.tileMb(i) ){
            if (vector.tileIsLocal( i, 0)) {
                auto tile = vector( i, 0 );
                auto tiledata = vector.data();
                return tiledata[index-i*(mt+1)];
            }
        }
    }
};

template <typename scalar_t>
slate::Matrix<scalar_t>* Vector<scalar_t>::Get(){
    return *vector;
};

template <typename scalar_t>
void Vector<scalar_t>::Multiply(slate::Matrix<scalar_t> matrix, scalar_t a, scalar_t b, slate::Option options){
    slate::multiply(matrix, a, b,  vector, vector, options);
};

template <typename scalar_t>
void Vector<scalar_t>::Add(Vector<scalar_t> vector, scalar_t a, scalar_t b, slate::Option options){
    slate::add(vector.Get(), a, vector, b, options);
};


/* ---- Dual Vector class ----*/
template <typename scalar_t>
DualVector<scalar_t>::DualVector(int m, int nb, int p, slate::Target target) {
    vector_up = slate::Matrix < scalar_t >(m, 1, nb, p, 1, MPI_COMM_WORLD);
    vector_up.insertLocalTiles(target);
    vector_down = slate::Matrix < scalar_t >(m, 1, nb, p, 1, MPI_COMM_WORLD);
    vector_down.insertLocalTiles(target);
};

template <typename scalar_t>
DualVector<scalar_t>::DualVector(int m, int nb, int p) {
    DualVector(m, nb, p, slate::HostTask);
};

template <typename scalar_t>
DualVector<scalar_t>::~DualVector() {
    vector_up.clear();
    vector_down.clear();
};

template <typename scalar_t>
void DualVector<scalar_t>::SetUp(int index, scalar_t value){
    int mt = vector_up.mt();
    for (int64_t i = 0; i < mt; i++){
        if ( i*(mt+1) <= index && index < i*(mt+1)+vector_up.tileMb(i) ){
            if (vector_up.tileIsLocal( i, 0)) {
                auto tile = vector_up( i, 0);
                auto tiledata = vector_up.data();
                tiledata[index-i*(mt+1)] = value;
                return;
            }
        }
    }
};

template <typename scalar_t>
scalar_t DualVector<scalar_t>::GetUp(int index){
    int mt = vector_up.mt();
    for (int64_t i = 0; i < mt; i++){
        if ( i*(mt+1) <= index && index < i*(mt+1)+vector_up.tileMb(i) ){
            if (vector_up.tileIsLocal( i, 0)) {
                auto tile = vector_up( i, 0 );
                auto tiledata = vector_up.data();
                return tiledata[index-i*(mt+1)];
            }
        }
    }
};

template <typename scalar_t>
slate::Matrix<scalar_t>* DualVector<scalar_t>::GetUp(){
    return *vector_up;
};

template <typename scalar_t>
void DualVector<scalar_t>::SetDown(int index, scalar_t value){
    int mt = vector_down.mt();
    for (int64_t i = 0; i < mt; i++){
        if ( i*(mt+1) <= index && index < i*(mt+1)+vector_down.tileMb(i) ){
            if (vector_down.tileIsLocal( i, 0)) {
                auto tile = vector_down( i, 0);
                auto tiledata = vector_down.data();
                tiledata[index-i*(mt+1)] = value;
                return;
            }
        }
    }
};

template <typename scalar_t>
scalar_t DualVector<scalar_t>::GetDown(int index){
    int mt = vector_down.mt();
    for (int64_t i = 0; i < mt; i++){
        if ( i*(mt+1) <= index && index < i*(mt+1)+vector_down.tileMb(i) ){
            if (vector_down.tileIsLocal( i, 0)) {
                auto tile = vector_down( i, 0 );
                auto tiledata = vector_down.data();
                return tiledata[index-i*(mt+1)];
            }
        }
    }
};

template <typename scalar_t>
slate::Matrix<scalar_t>* DualVector<scalar_t>::GetDown(){
    return *vector_down;
};

template <typename scalar_t>
void DualVector<scalar_t>::Set(int index, scalar_t value){
    SetUp(index, value);
    SetDown(index, value);
}

template <typename scalar_t>
void DualVector<scalar_t>::MultiplyUp(slate::Matrix<scalar_t> matrix, scalar_t a, scalar_t b, slate::Option options){
    slate::multiply(matrix, a, b,  vector_up, vector_up, options);
};

template <typename scalar_t>
void DualVector<scalar_t>::MultiplyDown(slate::Matrix<scalar_t> matrix, scalar_t a, scalar_t b, slate::Option options){
    slate::multiply(matrix, a, b,  vector_down, vector_down, options);
};

template <typename scalar_t>
void DualVector<scalar_t>::Multiply(slate::Matrix<scalar_t> matrix, scalar_t a, scalar_t b, slate::Option options){
    MultiplyUp(matrix, a, b, options);
    MultiplyDown(matrix, a, b, options);
}

template <typename scalar_t>
void DualVector<scalar_t>::AddUp(Vector<scalar_t> vector, scalar_t a, scalar_t b, slate::Option options){
    slate::add(vector.Get(), a, vector_up, b, options);
};

template <typename scalar_t>
void DualVector<scalar_t>::AddDown(Vector<scalar_t> vector, scalar_t a, scalar_t b, slate::Option options){
    slate::add(vector.Get(), a, vector_down, b, options);
};

template <typename scalar_t>
void DualVector<scalar_t>::Add(DualVector<scalar_t> vector, scalar_t a, scalar_t b, slate::Option options){
    slate::add(vector.GetUp(), a, vector_up, b, options);
    slate::add(vector.GetDown(), a, vector_down, b, options);
};

template <typename scalar_t>
void DualVector<scalar_t>::Add(Vector<scalar_t> vector, scalar_t a, scalar_t b, slate::Option options){
    slate::add(vector.Get(), a, vector_up, b, options);
    slate::add(vector.Get(), a, vector_down, b, options);
};


/* ---- Dual Matrix ----*/
template <typename scalar_t>
DualMatrix<scalar_t>::DualMatrix(int m, int nb, int p, slate::Target target) {
    matrix_up = slate::Matrix < scalar_t >(m, m, nb, p, 1, MPI_COMM_WORLD);
    matrix_up.insertLocalTiles(target);
    matrix_down = slate::Matrix < scalar_t >(m, m, nb, p, 1, MPI_COMM_WORLD);
    matrix_down.insertLocalTiles(target);
};

template <typename scalar_t>
DualMatrix<scalar_t>::DualMatrix(int m, int nb, int p) {
    DualVector(m, nb, p, slate::HostTask);
};

template <typename scalar_t>
DualMatrix<scalar_t>::~DualMatrix() {
    matrix_up.clear();
    matrix_down.clear();
};

template <typename scalar_t>
void DualMatrix<scalar_t>::SetUp(int row, int col, scalar_t value){
    int mt = matrix_up.mt();
    int nt = matrix_up.nt();
    for (int64_t i = 0; i < mt; i++){
        for (int64_t j = 0; j < nt; ++j) {
            if (i * (mt + 1) <= row && row < i * (mt + 1) + matrix_up.tileMb(i) &&
                j * (nt + 1) <= col && col < i * (nt + 1) + matrix_up.tileNb(j)) {
                if (matrix_up.tileIsLocal(i, j)) {
                    auto tile = matrix_up(i, j);
                    auto tiledata = matrix_up.data();
                    tiledata[(row - i * (mt + 1)) +
                             (col - j * (nt + 1))*tile.stride()] = value;
                    return;
                }
            }
        }
    }
};

template <typename scalar_t>
scalar_t DualMatrix<scalar_t>::GetUp(int row, int col){
    int mt = matrix_up.mt();
    int nt = matrix_up.nt();
    for (int64_t i = 0; i < mt; i++){
        for (int64_t j = 0; j < nt; ++j) {
            if (i * (mt + 1) <= row && row < i * (mt + 1) + matrix_up.tileMb(i) &&
                j * (nt + 1) <= col && col < i * (nt + 1) + matrix_up.tileNb(j)) {
                if (matrix_up.tileIsLocal(i, j)) {
                    auto tile = matrix_up(i, j);
                    auto tiledata = matrix_up.data();
                    return tiledata[(row - i * (mt + 1)) +
                             (col - j * (nt + 1))*tile.stride()];
                }
            }
        }
    }
};

template <typename scalar_t>
slate::Matrix<scalar_t>* DualMatrix<scalar_t>::GetUp(){
    return *matrix_up;
};

template <typename scalar_t>
void DualMatrix<scalar_t>::SetDown(int row, int col, scalar_t value){
    int mt = matrix_down.mt();
    int nt = matrix_down.nt();
    for (int64_t i = 0; i < mt; i++){
        for (int64_t j = 0; j < nt; ++j) {
            if (i * (mt + 1) <= row && row < i * (mt + 1) + matrix_down.tileMb(i) &&
                j * (nt + 1) <= col && col < i * (nt + 1) + matrix_down.tileNb(j)) {
                if (matrix_down.tileIsLocal(i, j)) {
                    auto tile = matrix_down(i, j);
                    auto tiledata = matrix_down.data();
                    tiledata[(row - i * (mt + 1)) +
                             (col - j * (nt + 1))*tile.stride()] = value;
                    return;
                }
            }
        }
    }
};

template <typename scalar_t>
scalar_t DualMatrix<scalar_t>::GetDown(int row, int col){
    int mt = matrix_down.mt();
    int nt = matrix_down.nt();
    for (int64_t i = 0; i < mt; i++){
        for (int64_t j = 0; j < nt; ++j) {
            if (i * (mt + 1) <= row && row < i * (mt + 1) + matrix_down.tileMb(i) &&
                j * (nt + 1) <= col && col < i * (nt + 1) + matrix_down.tileNb(j)) {
                if (matrix_down.tileIsLocal(i, j)) {
                    auto tile = matrix_down(i, j);
                    auto tiledata = matrix_down.data();
                    return tiledata[(row - i * (mt + 1)) +
                                    (col - j * (nt + 1))*tile.stride()];
                }
            }
        }
    }
};

template <typename scalar_t>
slate::Matrix<scalar_t>* DualMatrix<scalar_t>::GetDown(){
    return *matrix_down;
};

template <typename scalar_t>
void DualMatrix<scalar_t>::Set(int row, int col, scalar_t value){
    SetUp(row, col, value);
    SetDown(row, col, value);
}

template <typename scalar_t>
void DualMatrix<scalar_t>::MultiplyUp(slate::Matrix<scalar_t> matrix, scalar_t a, scalar_t b, slate::Option options){
    slate::multiply(matrix, a, b,  matrix_up, matrix_up, options);
};

template <typename scalar_t>
void DualMatrix<scalar_t>::MultiplyDown(slate::Matrix<scalar_t> matrix, scalar_t a, scalar_t b, slate::Option options){
    slate::multiply(matrix, a, b,  matrix_down, matrix_down, options);
};

template <typename scalar_t>
void DualMatrix<scalar_t>::Multiply(slate::Matrix<scalar_t> matrix, scalar_t a, scalar_t b, slate::Option options){
    MultiplyUp(matrix, a, b, options);
    MultiplyDown(matrix, a, b, options);
}

template <typename scalar_t>
void DualMatrix<scalar_t>::AddUp(slate::Matrix<scalar_t> matrix, scalar_t a, scalar_t b, slate::Option options){
    slate::add(matrix, a, matrix_up, b, options);
};

template <typename scalar_t>
void DualMatrix<scalar_t>::AddDown(slate::Matrix<scalar_t> matrix, scalar_t a, scalar_t b, slate::Option options){
    slate::add(matrix, a, matrix_down, b, options);
};

template <typename scalar_t>
void DualMatrix<scalar_t>::Add(DualMatrix<scalar_t> matrix, scalar_t a, scalar_t b, slate::Option options){
    slate::add(matrix.GetUp(), a, matrix_up, b, options);
    slate::add(matrix.GetDown(), a, matrix_down, b, options);
};

template <typename scalar_t>
void DualMatrix<scalar_t>::Add(slate::Matrix<scalar_t> matrix, scalar_t a, scalar_t b, slate::Option options){
    slate::add(matrix, a, matrix_up, b, options);
    slate::add(matrix, a, matrix_down, b, options);
};

template <typename scalar_t>
DualVector<scalar_t>::EmptyLike(){
    DualVector<scalar_t> vec;
    vec.GetUp() = vector_up.emptyLike();
    vec.GetDown() = vector_down.emptyLike();
    return vec;
};
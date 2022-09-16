#include "slate/slate.hh"
#include <blas.hh>
#include <mpi.h>

#pragma once

template <typename scalar_t>
class Vector{
public:
    Vector(int m,int nb, int p);
    Vector(int m,int nb, int p, slate::Target target);
    ~Vector();

    void Set(int index, scalar_t value);

    scalar_t Get(int index);
    slate::Matrix<scalar_t>* Get();

    void Multiply(slate::Matrix<scalar_t> matrix, scalar_t a = 1, scalar_t b = 1, slate::Option options = {
            // Set execution target to cpu Devices
            { slate :: Option :: Target , slate :: Target :: HostTask }});

    void Add(Vector<scalar_t> vector, scalar_t a = 1, scalar_t b = 1, slate::Option options = {
            // Set execution target to cpu Devices
            { slate :: Option :: Target , slate :: Target :: HostTask }});

private:
    slate::Matrix<scalar_t> vector;
};

template <typename scalar_t>
class DualVector{
public:
    DualVector(int m,int nb, int p);
    DualVector(int m,int nb, int p, slate::Target target);
    ~DualVector();

    void Set(int index, scalar_t value);
    void SetUp(int index, scalar_t value);
    void SetDown(int index, scalar_t value);

    scalar_t GetUp(int index);
    scalar_t GetDown(int index);
    slate::Matrix<scalar_t>* GetUp();
    slate::Matrix<scalar_t>* GetDown();

    void Multiply(slate::Matrix<scalar_t> matrix, scalar_t a = 1, scalar_t b = 1, slate::Option options = {
            // Set execution target to cpu Devices
            { slate :: Option :: Target , slate :: Target :: HostTask }});

    void MultiplyUp(slate::Matrix<scalar_t> matrix, scalar_t a = 1, scalar_t b = 1, slate::Option options = {
            // Set execution target to cpu Devices
            { slate :: Option :: Target , slate :: Target :: HostTask }});

    void MultiplyDown(slate::Matrix<scalar_t> matrix, scalar_t a = 1, scalar_t b = 1, slate::Option options = {
            // Set execution target to cpu Devices
            { slate :: Option :: Target , slate :: Target :: HostTask }});

    void AddUp(Vector<scalar_t> vector, scalar_t a = 1, scalar_t b = 1, slate::Option options = {
            // Set execution target to cpu Devices
            { slate :: Option :: Target , slate :: Target :: HostTask }});
    void AddDown(Vector<scalar_t> vector, scalar_t a = 1, scalar_t b = 1, slate::Option options = {
            // Set execution target to cpu Devices
            { slate :: Option :: Target , slate :: Target :: HostTask }});
    void Add(DualVector<scalar_t> vector, scalar_t a = 1, scalar_t b = 1, slate::Option options = {
            // Set execution target to cpu Devices
            { slate :: Option :: Target , slate :: Target :: HostTask }});
    void Add(Vector<scalar_t> vector, scalar_t a = 1, scalar_t b = 1, slate::Option options = {
            // Set execution target to cpu Devices
            { slate :: Option :: Target , slate :: Target :: HostTask }});

private:
    slate::Matrix<scalar_t> vector_up;
    slate::Matrix<scalar_t> vector_down;
};

template <typename scalar_t>
class DualMatrix{
public:
    DualMatrix(int m,int nb, int p);
    DualMatrix(int m,int nb, int p, slate::Target target);
    ~DualMatrix();

    void Set(int row, int col, scalar_t value);
    void SetUp(int row, int col, scalar_t value);
    void SetDown(int row, int col, scalar_t value);

    scalar_t GetUp(int row, int col);
    scalar_t GetDown(int row, int col);
    slate::Matrix<scalar_t>* GetUp();
    slate::Matrix<scalar_t>* GetDown();

    void Multiply(slate::Matrix<scalar_t> matrix, scalar_t a = 1, scalar_t b = 1, slate::Option options = {
            // Set execution target to cpu Devices
            { slate :: Option :: Target , slate :: Target :: HostTask }});

    void MultiplyUp(slate::Matrix<scalar_t> matrix, scalar_t a = 1, scalar_t b = 1, slate::Option options = {
            // Set execution target to cpu Devices
            { slate :: Option :: Target , slate :: Target :: HostTask }});

    void MultiplyDown(slate::Matrix<scalar_t> matrix, scalar_t a = 1, scalar_t b = 1, slate::Option options = {
            // Set execution target to cpu Devices
            { slate :: Option :: Target , slate :: Target :: HostTask }});

    void AddUp(slate::Matrix<scalar_t> matrix, scalar_t a = 1, scalar_t b = 1, slate::Option options = {
            // Set execution target to cpu Devices
            { slate :: Option :: Target , slate :: Target :: HostTask }});
    void AddDown(slate::Matrix<scalar_t> matrix, scalar_t a = 1, scalar_t b = 1, slate::Option options = {
            // Set execution target to cpu Devices
            { slate :: Option :: Target , slate :: Target :: HostTask }});
    void Add(DualMatrix<scalar_t> matrix, scalar_t a = 1, scalar_t b = 1, slate::Option options = {
            // Set execution target to cpu Devices
            { slate :: Option :: Target , slate :: Target :: HostTask }});
    void Add(slate::Matrix<scalar_t> matrix, scalar_t a = 1, scalar_t b = 1, slate::Option options = {
            // Set execution target to cpu Devices
            { slate :: Option :: Target , slate :: Target :: HostTask }});

    DualVector<scalar_t> EmptyLike();

private:
    slate::Matrix<scalar_t> matrix_up;
    slate::Matrix<scalar_t> matrix_down;
};

#include "utils.ipp"


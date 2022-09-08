#include "slate/slate.hh"
#ifndef DFT_UTILS_H
#define DFT_UTILS_H

#ifndef Vector_H
#define Vector_H
template <class scalar_t>
class Vector {
public:
    Vector(int m,int nb,int p, slate::Target target = slate::Target::HostTask);

    void Set(int index, scalar_t value);

    scalar_t Get(int index);

    void Multiply(slate::Matrix<scalar_t> matrix, scalar_t a = 1, scalar_t b = 1, slate::Option options = {
            // Set execution target to cpu Devices
            { slate :: Option :: Target , slate :: Target :: HostTask }});

private:
    slate::Matrix<scalar_t> vector_ptr;
};
#endif // Vector_H

#endif //DFT_UTILS_H

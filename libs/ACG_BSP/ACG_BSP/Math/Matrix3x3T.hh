#ifndef ACG_MATH_MATRIX3X3T_HH_
#define ACG_MATH_MATRIX3X3T_HH_

#include <array>
#include <ostream>
#include <ACG_BSP/Math/VectorT.hh>
#include <algorithm>
#include <cmath>
#include <utility>

#if defined(_MSC_VER) && _MSC_VER < 1900
#define constexpr
typedef unsigned char uint_fast8_t;
#endif

namespace ACG {

/**
 * Small, kinda fast 3x3 matrix class.
 */
template<typename Scalar>
class Matrix3x3T {
    public:
        typedef typename OpenMesh::VectorT<Scalar, 3> Vec3;

        static Matrix3x3T<Scalar> fromColumns(Vec3 c1, Vec3 c2, Vec3 c3) {
            return Matrix3x3T<Scalar> {{
                c1[0], c2[0], c3[0],
                c1[1], c2[1], c3[1],
                c1[2], c2[2], c3[2],
            }};
        }

        static Matrix3x3T<Scalar> fromRows(Vec3 r1, Vec3 r2, Vec3 r3) {
            return Matrix3x3T<Scalar> {{
                r1[0], r1[1], r1[2],
                r2[0], r2[1], r2[2],
                r3[0], r3[1], r3[2]
            }};
        }

        static constexpr Matrix3x3T<Scalar> identity() {
            return {{
                1, 0, 0,
                0, 1, 0,
                0, 0, 1
            }};
        }

        static constexpr Matrix3x3T<Scalar> zero() {
            return {{
                0, 0, 0,
                0, 0, 0,
                0, 0, 0
            }};
        }

    public:
        Matrix3x3T() = default;

        /**
         * Initialize matrix from array in row major format.
         */
        constexpr Matrix3x3T(std::array<Scalar, 9> row_major) :
                values_(std::move(row_major)) {}

        constexpr bool operator==(const Matrix3x3T &rhs) const {
            return values_ == rhs.values_;
        }

        //Matrix3x3T(std::initializer_list<Scalar> row_major) {
        //    static_assert(row_major.size() == 9, "Need exactly 9 values.");
        //    std::copy(row_major.begin(), row_major.end(), this->begin());
        //}

        /// Map row/column index to linearized index.
        constexpr static uint_fast8_t indexof(uint_fast8_t r, uint_fast8_t c) {
            return r*3+c;
        }

        constexpr const Scalar &operator() (uint_fast8_t r, uint_fast8_t c) const {
            return values_[indexof(r, c)];
        }

        Scalar &operator() (uint_fast8_t r, uint_fast8_t c) {
            return values_[indexof(r, c)];
        }

        /// Linearized row major access.
        constexpr const Scalar &operator[] (uint_fast8_t i) const {
            return values_[i];
        }

        /// Linearized row major access.
        Scalar &operator[] (uint_fast8_t i) {
            return values_[i];
        }

        Vec3 getRow(uint_fast8_t r) const {
            return Vec3((*this)(r,0), (*this)(r,1), (*this)(r,2));
        }
        void setRow(uint_fast8_t r, const Vec3 &v) {
            for (uint_fast8_t c = 0; c < 3; ++c) {
                (*this)(r, c) = v[c];
            }
        }
        Vec3 getCol(uint_fast8_t c) const {
            return Vec3((*this)(0,c), (*this)(1,c), (*this)(2,c));
        }
        void setCol(uint_fast8_t c, const Vec3 &v) {
            for (uint_fast8_t r = 0; r < 3; ++r) {
                (*this)(r, c) = v[r];
            }
        }

        constexpr Matrix3x3T operator*(const Matrix3x3T &rhs) const {
            return Matrix3x3T {{
                (*this)(0, 0) * rhs(0, 0) + (*this)(0, 1) * rhs(1, 0) + (*this)(0, 2) * rhs(2, 0),
                (*this)(0, 0) * rhs(0, 1) + (*this)(0, 1) * rhs(1, 1) + (*this)(0, 2) * rhs(2, 1),
                (*this)(0, 0) * rhs(0, 2) + (*this)(0, 1) * rhs(1, 2) + (*this)(0, 2) * rhs(2, 2),

                (*this)(1, 0) * rhs(0, 0) + (*this)(1, 1) * rhs(1, 0) + (*this)(1, 2) * rhs(2, 0),
                (*this)(1, 0) * rhs(0, 1) + (*this)(1, 1) * rhs(1, 1) + (*this)(1, 2) * rhs(2, 1),
                (*this)(1, 0) * rhs(0, 2) + (*this)(1, 1) * rhs(1, 2) + (*this)(1, 2) * rhs(2, 2),

                (*this)(2, 0) * rhs(0, 0) + (*this)(2, 1) * rhs(1, 0) + (*this)(2, 2) * rhs(2, 0),
                (*this)(2, 0) * rhs(0, 1) + (*this)(2, 1) * rhs(1, 1) + (*this)(2, 2) * rhs(2, 1),
                (*this)(2, 0) * rhs(0, 2) + (*this)(2, 1) * rhs(1, 2) + (*this)(2, 2) * rhs(2, 2),
            }};
        }

        template<typename OtherScalar>
        constexpr auto operator*(const VectorT<OtherScalar,3> &rhs) const
        -> OpenMesh::VectorT<decltype(std::declval<Scalar>() * rhs[0]), 3>
        {
            return {
                (*this)(0, 0) * rhs[0] + (*this)(0, 1) * rhs[1] + (*this)(0, 2) * rhs[2],
                (*this)(1, 0) * rhs[0] + (*this)(1, 1) * rhs[1] + (*this)(1, 2) * rhs[2],
                (*this)(2, 0) * rhs[0] + (*this)(2, 1) * rhs[1] + (*this)(2, 2) * rhs[2]
            };
        }

        template<typename OtherScalar>
        constexpr friend auto operator*(VectorT<OtherScalar,3> v, const Matrix3x3T &rhs)
        -> OpenMesh::VectorT<decltype(std::declval<Scalar>() * v[0]), 3>
        {
            return {
                rhs(0, 0) * v[0] + rhs(0, 1) * v[1] + rhs(0, 2) * v[2],
                rhs(1, 0) * v[0] + rhs(1, 1) * v[1] + rhs(1, 2) * v[2],
                rhs(2, 0) * v[0] + rhs(2, 1) * v[1] + rhs(2, 2) * v[2]
            };
        }

        constexpr Matrix3x3T operator*(Scalar c) const {
            return Matrix3x3T {{
                (*this)[0] * c, (*this)[1] * c, (*this)[2] * c,
                (*this)[3] * c, (*this)[4] * c, (*this)[5] * c,
                (*this)[6] * c, (*this)[7] * c, (*this)[8] * c,
            }};
        }

        constexpr friend Matrix3x3T operator*(Scalar c, const Matrix3x3T &rhs) {
            return Matrix3x3T {{
                rhs[0] * c, rhs[1] * c, rhs[2] * c,
                rhs[3] * c, rhs[4] * c, rhs[5] * c,
                rhs[6] * c, rhs[7] * c, rhs[8] * c,
            }};
        }

        constexpr Matrix3x3T operator+ (const Matrix3x3T &rhs) const {
            return Matrix3x3T {{
                (*this)[0] + rhs[0], (*this)[1] + rhs[1], (*this)[2] + rhs[2],
                (*this)[3] + rhs[3], (*this)[4] + rhs[4], (*this)[5] + rhs[5],
                (*this)[6] + rhs[6], (*this)[7] + rhs[7], (*this)[8] + rhs[8],
            }};
        }

        constexpr Matrix3x3T operator- (const Matrix3x3T &rhs) const {
            return Matrix3x3T {{
                (*this)[0] - rhs[0], (*this)[1] - rhs[1], (*this)[2] - rhs[2],
                (*this)[3] - rhs[3], (*this)[4] - rhs[4], (*this)[5] - rhs[5],
                (*this)[6] - rhs[6], (*this)[7] - rhs[7], (*this)[8] - rhs[8],
            }};
        }

        constexpr Matrix3x3T operator- () const {
            return Matrix3x3T {{
                -values_[0], -values_[1], -values_[2],
                -values_[3], -values_[4], -values_[5],
                -values_[6], -values_[7], -values_[8]
            }};
        }

        const Matrix3x3T &operator*=(const Matrix3x3T &rhs) {
            (*this) = operator*(rhs);
            return *this;
        }

        constexpr Scalar det() const {
            return (*this)(0, 0) * ((*this)(1, 1) * (*this)(2, 2) - (*this)(2, 1) * (*this)(1, 2)) -
                   (*this)(0, 1) * ((*this)(1, 0) * (*this)(2, 2) - (*this)(1, 2) * (*this)(2, 0)) +
                   (*this)(0, 2) * ((*this)(1, 0) * (*this)(2, 1) - (*this)(1, 1) * (*this)(2, 0));

            /*
            return (*this)(0, 0) * (*this)(1, 1) * (*this)(2, 2) +
                   (*this)(1, 0) * (*this)(2, 1) * (*this)(0, 2) +
                   (*this)(2, 0) * (*this)(0, 1) * (*this)(1, 2) -
                   (*this)(0, 0) * (*this)(2, 1) * (*this)(1, 2) -
                   (*this)(2, 0) * (*this)(1, 1) * (*this)(0, 2) -
                   (*this)(1, 0) * (*this)(0, 1) * (*this)(2, 2)
             */
        }

        constexpr Scalar trace() const {
            return (*this)[0] + (*this)[4] + (*this)[8];
        }

        void transpose() {
            std::swap(values_[1], values_[3]);
            std::swap(values_[2], values_[6]);
            std::swap(values_[5], values_[7]);
        }

        constexpr Matrix3x3T transposed() const {
            return Matrix3x3T {{
                values_[0], values_[3], values_[6],
                values_[1], values_[4], values_[7],
                values_[2], values_[5], values_[8],
            }};
        }

        void invert() {
            *this = inverse();
        }

        Matrix3x3T inverse() const {
            const Scalar invdet = 1.0 / det();
            return Matrix3x3T {{
                ((*this)(1, 1) * (*this)(2, 2) - (*this)(2, 1) * (*this)(1, 2)) * invdet,
                ((*this)(0, 2) * (*this)(2, 1) - (*this)(0, 1) * (*this)(2, 2)) * invdet,
                ((*this)(0, 1) * (*this)(1, 2) - (*this)(0, 2) * (*this)(1, 1)) * invdet,
                ((*this)(1, 2) * (*this)(2, 0) - (*this)(1, 0) * (*this)(2, 2)) * invdet,
                ((*this)(0, 0) * (*this)(2, 2) - (*this)(0, 2) * (*this)(2, 0)) * invdet,
                ((*this)(1, 0) * (*this)(0, 2) - (*this)(0, 0) * (*this)(1, 2)) * invdet,
                ((*this)(1, 0) * (*this)(2, 1) - (*this)(2, 0) * (*this)(1, 1)) * invdet,
                ((*this)(2, 0) * (*this)(0, 1) - (*this)(0, 0) * (*this)(2, 1)) * invdet,
                ((*this)(0, 0) * (*this)(1, 1) - (*this)(1, 0) * (*this)(0, 1)) * invdet,
            }};
        }

        constexpr Scalar frobeniusSquared() const {
            return std::inner_product(
                    values_.begin(), values_.end(), values_.begin(), Scalar(0.0));
        }

        constexpr double frobenius() const {
            return std::sqrt(frobeniusSquared());
        }


        friend
        std::ostream &operator<< (std::ostream &os, const Matrix3x3T &m) {
            os << "[[" << m[0] << ", " << m[1] << ", " << m[2] << "], "
                  "[" << m[3] << ", " << m[4] << ", " << m[5] << "], "
                  "[" << m[6] << ", " << m[7] << ", " << m[8] << "]]";
            return os;
        }

    private:
        std::array<Scalar, 9> values_;
};

typedef Matrix3x3T<float> Matrix3x3f;
typedef Matrix3x3T<double> Matrix3x3d;

} /* namespace ACG */

#if defined(_MSC_VER) && _MSC_VER < 1900
#undef constexpr
#endif

#endif /* ACG_MATH_MATRIX3X3T_HH_ */

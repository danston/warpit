/*! FILE Vector2D.hpp

 Class Vector2D.
 
 The original code is due to:
 Author: Frank Firsching
 Date: 17.03.2001
 
 Reimplementation is due to:
 Author: Dmitry Anisimov
 Date: 18.06.2013
 Mail: danston@ymail.com
 Web: http://www.anisimovdmitry.com
 
*/

/*! TROUBLES

 With all the questions related to the class, please let me know about them by the email
 provided at the top of this file.
 
*/

#ifndef Vector2D_HPP
#define Vector2D_HPP

// STL includes.
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

// warpit namespace.
namespace warpit {

    template<typename T> class Vector2D;

    //! Some predefined types for Vector2D ----->

    //! Vector2D with int
    typedef Vector2D<int> Vector2Di;

    //! Vector2D with float
    typedef Vector2D<float> Vector2Df;

    //! Vector2D with double
    typedef Vector2D<double> Vector2Dd;
    typedef    Vector2Dd      Vector2d;

    //! Compute the Pi value up to T precision
    template<typename T> inline T Pi()
    {
        Vector2D<T>vStart = Vector2D<T>(T(1), T(0));
        Vector2D<T>vEnd = Vector2D<T>(T(-1), T(0));
        return vStart.unsignedAngleRad(vEnd);
    }

    //! Cotangent of the angle between two vectors.
    template<typename T> inline T cotangent(const Vector2D<T> &a, const Vector2D<T> &b)
    {
        return (a.scalarProduct(b)) / std::fabs(a.crossProduct(b));
    }

    //! Compute the determinant of three vectors in the plane.
    template<typename T> inline T determinant(const Vector2D<T> &a, const Vector2D<T> &b, const Vector2D<T> &c)
    {
       return  (a.x - c.x)*(b.y - c.y) - (b.x - c.x)*(a.y - c.y);
    }

    //! This class provides functions to handle 2-dimensional vectors
    template<typename T>
        class Vector2D
    {
    public:
        //! typedefs
        typedef T Type;

        //! x coordinate of the vector
        T x;
        //! y coordinate of the vector
        T y;

        //! Default constructor
        Vector2D() { x = T(0); y = T(0); }
        //! Initialize Vector2D with one explicit coordinate
        Vector2D(const T _x) { x = _x; y = _x; }
        //! Initialize Vector2D with two explicit coordinates
        Vector2D(const T _x, const T _y) { x = _x; y = _y; }
        //! Copy-constructor
        Vector2D(const Vector2D &v) { x = v.x; y = v.y; }
        //! Initialize Vector2D with a pointer to an array
        Vector2D(T* a) { x = a[0]; y = a[1]; }
        //! Initialize Vector2D with a standard C++ vector
        Vector2D(const std::vector<T> &v)
        {
            const int size = v.size();
            if(size == 0) { x = T(0); y = T(0); }
            if(size == 1) { x = v[0]; y = v[0]; }
            if(size >= 2) { x = v[0]; y = v[1]; }
        }

        //! Destructor
        ~Vector2D() { }

        //! Access to vector-coordinates with []-operator for a constant object
        inline T operator[](const int i) const;
        //! Access to vector-coordinates with []-operator for a manipulation
        inline T& operator[](const int i);

        //! Assign a value to a Vector2D
        inline Vector2D& operator=(const T scalar);
        //! Assign the values of a Vector2D to another one
        inline Vector2D& operator=(const Vector2D &v);
        //! Assign the values of a standard C++ array to a Vector2D
        inline Vector2D& operator=(T* a);
        //! Assign the values of a standard C++ vector to a Vector2D
        inline Vector2D& operator=(const std::vector<T> &v);

        //! Check for equality
        inline bool operator==(const Vector2D &v) const;
        //! Check for inequality
        inline bool operator!=(const Vector2D &v) const;

        //! Addition of 2 vectors
        inline Vector2D operator+(const Vector2D &v) const;
        //! Subtraction of 2 vectors
        inline Vector2D operator-(const Vector2D &v) const;
        //! Negate the vector
        inline Vector2D operator-() const;
        //! Add vector v to the operand
        inline void operator+=(const Vector2D &v);
        //! Subtract vector v from the operand
        inline void operator-=(const Vector2D &v);

        //! Multiply a vector by a scalar
        inline Vector2D operator*(const T scalar) const;
        //! Multiply the operand by a scalar
        inline void operator*=(const T scalar);
        //! Divide a vector by a scalar
        inline Vector2D operator/(const T scalar) const;
        //! Divide the operand by a scalar
        inline void operator/=(const T scalar);

        //! Calculate the scalar-product of two vectors
        inline T operator*(const Vector2D &v) const;
        //! Calculate the scalar-product of two vectors
        inline T scalarProduct(const Vector2D &v) const;

        //! Calculate the cross-product of two vectors
        inline T operator%(const Vector2D &v) const;
        //! Calculate the cross-product of two vectors
        inline T crossProduct(const Vector2D &v) const;

        //! Return the signed angle between the current vector and v in radians
        inline T signedAngleRad(const Vector2D &v) const;
        //! Return the angle between the current vector and v in radians
        inline T unsignedAngleRad(const Vector2D &v) const;
        //! Return the angle between the current vector and v in degrees
        inline T angleDeg(const Vector2D &v) const;

        //! Return the squared length of the vector
        inline T squaredLength() const;
        //! Return the length of the vector
        inline T length() const;
        //! Return the L1-Norm of the vector (Taxicab norm or Manhattan norm)
        inline T L1Norm() const;
        //! Return the L2-Norm of the vector (Euclidean norm)
        inline T L2Norm() const;
        //! Return the LInfinity-Norm of the vector (Maximum norm, Uniform norm, or Supremum norm)
        inline T LInfNorm() const;

        //! Return the index of the smallest coordinate of the vector
        inline int indexOfTheSmallestCoordinate() const;
        //! Return the value of the smallest coordinate of the vector
        inline T theSmallestCoordinate() const;
        //! Return the index of the biggest coordinate of the vector
        inline int indexOfTheBiggestCoordinate() const;
        //! Return the value of the biggest coordinate of the vector
        inline T theBiggestCoordinate() const;

        //! Set the length of the vector to one, but the orientation stays the same
        inline void normalize();
        //! Return the normalized vector as a copy
        inline Vector2D normalized() const;

        //! Reflect the current vector accross a given v through the origin
        inline void reflect(const Vector2D &v);
        //! Return a vector reflected accross a given v through the origin
        Vector2D reflected(const Vector2D &v) const;

        //! Rotate the current vector by the given angle
        inline void rotate(const T angle);
        //! Return a vector obtained by rotation of the current vector by the given angle
        inline Vector2D rotated(const T angle) const;

        //! Horizontal shear of the current vector by the given factor
        inline void horizontalShear(const T factor);
        //! Return a vector obtained by horizontal shear of the current vector by the given factor
        inline Vector2D horizontallySheared(const T factor) const;

        //! Vertical shear of the current vector by the given factor
        inline void verticalShear(const T factor);
        //! Return a vector obtained by vertical shear of the current vector by the given factor
        inline Vector2D verticallySheared(const T factor) const;

        //! Scale the current vector by the given vector
        inline void scale(const Vector2D &v);
        //! Return a vector obtained by scaling the current vector by the given vector
        inline Vector2D scaled(const Vector2D &v) const;

        //! Translate the current vector by the given vector
        inline void translate(const Vector2D &v);
        //! Return a vector obtained by translating the current vector by the given vector
        inline Vector2D translated(const Vector2D &v) const;

        //! Set all the coordinates of the vector to their absolute values
        inline void abs();
    };

    /// Operator [] ---------->

    template<typename T> inline T Vector2D<T>::operator[](const int i) const
    {
        assert( ( i >= 0 ) && ( i < 2 ) );
        return (&x)[i];
    }

    template<typename T> inline T& Vector2D<T>::operator[](const int i)
    {
        assert( ( i >= 0 ) && ( i < 2 ));
        return (&x)[i];
    }

    /// Operator = ---------->

    template<typename T> inline Vector2D<T>& Vector2D<T>::operator=(const T scalar)
    {
        x = y = scalar;
        return *this;
    }

    template<typename T> inline Vector2D<T>& Vector2D<T>::operator=(const Vector2D &v)
    {
        x = v.x; y = v.y;
        return *this;
    }

    template<typename T> inline Vector2D<T>& Vector2D<T>::operator=(T* a)
    {
        x = a[0]; y = a[1];
        return *this;
    }

    template<typename T> inline Vector2D<T>& Vector2D<T>::operator=(const std::vector<T> &v)
    {
        const int size = v.size();
        if(size == 0) { x = T(0); y = T(0); }
        if(size == 1) { x = v[0]; y = v[0]; }
        if(size >= 2) { x = v[0]; y = v[1]; }
        return *this;
    }

    /// Operators == and != ---------->

    template<typename T> inline bool Vector2D<T>::operator==(const Vector2D &v) const
    {
        return ( fabs(x - v.x) <= std::numeric_limits<T>::epsilon() && fabs(y - v.y) <= std::numeric_limits<T>::epsilon() );
    }

    template<typename T> inline bool Vector2D<T>::operator!=(const Vector2D &v) const
    {
        return ( fabs(x - v.x) > std::numeric_limits<T>::epsilon() || fabs(y - v.y) > std::numeric_limits<T>::epsilon() );
    }

    /// Operators - and + ---------->

    template<typename T> inline Vector2D<T> Vector2D<T>::operator+(const Vector2D &v) const
    {
        return Vector2D<T>(x + v.x, y + v.y);
    }

    template<typename T> inline Vector2D<T> Vector2D<T>::operator-(const Vector2D &v) const
    {
        return Vector2D<T>(x - v.x, y - v.y);
    }

    template<typename T> inline Vector2D<T> Vector2D<T>::operator-() const
    {
        return Vector2D<T>(-x, -y);
    }

    template<typename T> inline void Vector2D<T>::operator+=(const Vector2D &v)
    {
        x += v.x; y += v.y;
    }

    template<typename T> inline void Vector2D<T>::operator-=(const Vector2D &v)
    {
        x -= v.x; y -= v.y;
    }

    /// Operators / and * ---------->

    template<typename T> inline Vector2D<T> Vector2D<T>::operator*(const T scalar) const
    {
        return Vector2D<T>(scalar*x, scalar*y);
    }

    template<typename T, typename S> inline Vector2D<T> operator*(const S scalar, const Vector2D<T> &v)
    {
        return v*scalar;
    }

    template<typename T> inline void Vector2D<T>::operator*=(const T scalar)
    {
        x *= scalar; y *= scalar;
    }

    template<typename T> inline Vector2D<T> Vector2D<T>::operator/(const T scalar) const
    {
        assert( fabs(scalar) > std::numeric_limits<T>::epsilon() );
        return Vector2D<T>(x / scalar, y / scalar);
    }

    template<typename T> inline void Vector2D<T>::operator/=(const T scalar)
    {
        assert( fabs(scalar) > std::numeric_limits<T>::epsilon() );
        x /= scalar; y /= scalar;
    }

    /// Scalar product ---------->

    template<typename T> inline T Vector2D<T>::operator*(const Vector2D &v) const
    {
        return this->scalarProduct(v);
    }

    template<typename T> inline T Vector2D<T>::scalarProduct(const Vector2D &v) const
    {
        return (x*v.x + y*v.y);
    }

    /// Cross product ---------->

    template<typename T> inline T Vector2D<T>::operator%(const Vector2D &v) const
    {
        return this->crossProduct(v);
    }

    template<typename T> inline T Vector2D<T>::crossProduct(const Vector2D &v) const
    {
        return (x*v.y - y*v.x);
    }

    /// Angles ---------->

    template<typename T> inline T Vector2D<T>::signedAngleRad(const Vector2D &v) const
    {
        return atan2(this->crossProduct(v), this->scalarProduct(v));
    }

    template<typename T> inline T Vector2D<T>::unsignedAngleRad(const Vector2D &v) const
    {
        return acos(this->scalarProduct(v) / (this->length()*v.length()));
    }

    template<typename T> inline T Vector2D<T>::angleDeg(const Vector2D &v) const
    {
        return this->unsignedAngleRad(v)*(T(180) / Pi<T>());
    }

    /// Norms ---------->

    template<typename T> inline T Vector2D<T>::squaredLength() const
    {
        return (x*x + y*y);
    }

    template<typename T> inline T Vector2D<T>::length() const
    {
        return sqrt(this->squaredLength());
    }

    template<typename T> inline T Vector2D<T>::L1Norm() const
    {
        return (fabs(x) + fabs(y));
    }

    template<typename T> inline T Vector2D<T>::L2Norm() const
    {
        return this->length();
    }

    template<typename T> inline T Vector2D<T>::LInfNorm() const
    {
        return max(fabs(x), fabs(y));
    }

    /// The smallest and the biggest coordinate ---------->

    template<typename T> inline int Vector2D<T>::indexOfTheSmallestCoordinate() const
    {
        return ( (x < y) ? 0 : 1 );
    }

    template<typename T> inline T Vector2D<T>::theSmallestCoordinate() const
    {
        return (&x)[this->indexOfTheSmallestCoordinate()];
    }

    template<typename T> inline int Vector2D<T>::indexOfTheBiggestCoordinate() const
    {
        return ( (x < y) ? 1 : 0 );
    }

    template<typename T> inline T Vector2D<T>::theBiggestCoordinate() const
    {
        return (&x)[this->indexOfTheBiggestCoordinate()];
    }

    /// Normalization ---------->

    template<typename T> inline void Vector2D<T>::normalize()
    {
        (*this) = this->normalized();
    }

    template<typename T> inline Vector2D<T> Vector2D<T>::normalized() const
    {
        T thisLength = this->length();
        if(fabs(thisLength) >  std::numeric_limits<T>::epsilon()) return Vector2D<T>(x / thisLength, y / thisLength);
        else return Vector2D<T>(T(0), T(1));
    }

    /// Reflection ---------->

    template<typename T> inline void Vector2D<T>::reflect(const Vector2D &v)
    {
        (*this) = this->reflected(v);
    }

    template<typename T> inline Vector2D<T> Vector2D<T>::reflected(const Vector2D &v) const
    {
        return (T(2)*(((*this)*v) / (v*v))*v - (*this));
    }

    /// Rotation ---------->

    template<typename T> inline void Vector2D<T>::rotate(const T angle)
    {
        (*this) = this->rotated(angle);
    }

    template<typename T> inline Vector2D<T> Vector2D<T>::rotated(const T angle) const
    {
        return Vector2D<T>(std::cos(angle) * x - std::sin(angle) * y, std::sin(angle) * x + std::cos(angle) * y);
    }

    /// Horizontal shear ---------->

    template<typename T> inline void Vector2D<T>::horizontalShear(const T factor)
    {
        (*this) = this->horizontallySheared(factor);
    }

    template<typename T> inline Vector2D<T> Vector2D<T>::horizontallySheared(const T factor) const
    {
        return Vector2D<T>(x + factor * y, y);
    }

    /// Vertical shear ---------->

    template<typename T> inline void Vector2D<T>::verticalShear(const T factor)
    {
        (*this) = this->verticallySheared(factor);
    }

    template<typename T> inline Vector2D<T> Vector2D<T>::verticallySheared(const T factor) const
    {
        return Vector2D<T>(x, factor * x + y);
    }

    /// Scaling ---------->

    template<typename T> inline void Vector2D<T>::scale(const Vector2D &v)
    {
        (*this) = this->scaled(v);
    }

    template<typename T> inline Vector2D<T> Vector2D<T>::scaled(const Vector2D &v) const
    {
        return Vector2D<T>(v.x * x, v.y * y);
    }

    /// Translation ---------->

    template<typename T> inline void Vector2D<T>::translate(const Vector2D &v)
    {
        (*this) = this->translated(v);
    }

    template<typename T> inline Vector2D<T> Vector2D<T>::translated(const Vector2D &v) const
    {
        return Vector2D<T>(x + v.x, y + v.y);
    }

    /// Absolute value ---------->

    template<typename T> inline void Vector2D<T>::abs()
    {
        x = fabs(x); y = fabs(y);
    }

    /// Read from or write to stream ---------->

    //! Write a 2D-vector to a stream
    template<typename T> inline std::ostream& operator<<(std::ostream &ostr, const Vector2D<T> &v)
    {
        return (ostr << "(" << v.x << ", " << v.y << ")");
    }

    //! Read a 2D-vector from a stream
    template<typename T> inline std::istream& operator>>(std::istream &istr, Vector2D<T> &v)
    {
        return (istr >> v.x >> v.y);
    }

} // namespace warpit

#endif // Vector2D_HPP

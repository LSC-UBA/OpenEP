#ifndef VECTORFIELD_H_
#define VECTORFIELD_H_

#include "scalar_field.h"

class VectorField
{
private:
    ScalarField * scalar_component;
    /* The mesh can be shared with other scalar field objects */
    const Mesh * mesh; 
    int ii;
    int jj;
    int kk;

    void save_csv(const std::string &filename);

    void save_vtk(const std::string &filename);

public:

    VectorField(const Mesh &mesh0);

    VectorField(const Mesh &mesh0, long double val);

    VectorField(const Mesh &mesh0, const ScalarField & vx, const ScalarField & vy, const ScalarField & vz);

    void set(const int coord, const int i, const int j, const int k, long double val);

    void set_coordinate(const int coord, const ScalarField & v);

    ScalarField & get_coordinate(const int coord);

    long double get(const int coord, const int i, const int j, const int k) const;

    const Mesh& get_mesh() const;

    void print() const;

    void save_components_vtk(const std::string &filename, const std::string &save_format);
    
    void save(const std::string &filename, const std::string &save_format);
    
    ScalarField get_norm();
    
    ScalarField & get_component ( int i );

    VectorField& operator =(const VectorField & other);

    long double & operator()(const int coord, const int i, const int j, const int k) const;

    VectorField operator +=(const VectorField & other) const;

    VectorField operator -=(const VectorField & other) const;

    ~VectorField();

};


#endif


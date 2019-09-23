#ifndef SCALARFIELD_H_
#define SCALARFIELD_H_

#include <fstream>
#include <sstream>

#include "mesh.h"

class ScalarField
{
private:
    int ii;
    int jj;
    int kk;
    int n;
    long double * p_mat;
    /* The mesh can be shared with other scalar field objects */
    const Mesh * mesh; 

    void save_vtk(const std::string &filename);
    
    void save_csv(const std::string &filename);
    

public:

    ScalarField();

    ScalarField(const Mesh &mesh0);

    ScalarField(const Mesh &mesh0, long double val);
    
    void init(const Mesh &mesh0, long double val);
    
    void init(long double val);

    void set(const int i, const int j, const int k, long double val);

    long double get(const int i, const int j, const int k) const;
    
    void save(const std::string &filename, const std::string &save_format);
    
    const Mesh& get_mesh() const;
    
    void print() const;

    long double get_norm() const;

    // It assums that ii, jj and kk doesn't change
    ScalarField& operator =(const ScalarField & other);

    long double& operator()(const int i, const int j, const int k) const;

    ScalarField operator +(const ScalarField & other) const;
    
    ScalarField operator *(const ScalarField & other) const;

    ScalarField operator -(const ScalarField & other) const;

    ~ScalarField();

};

#endif

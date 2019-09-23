#include "vector_field.h"

void VectorField::save_csv(const std::string &filename)
{
    std::ofstream file;
    // Open file
    file.open(filename.c_str());

    // Write content
    file << "X,Y,Z,I,J,K,Vx,Vy,Vz" << std::endl;
    for (int k = 0; k < kk; k++)
    {
        for (int j = 0; j < jj; j++)
        {
            for (int i = 0; i < ii; i++)
            {
                file    << mesh->get_x(i) << "," << mesh->get_y(j) << "," << mesh->get_z(k)
                        << "," << i << "," << j << "," << k << "," 
                        << get(0,i,j,k) << "," << get(1,i,j,k) << "," << get(2,i,j,k) << std::endl;
            }
        }
    }
    // Close file name
    file.close();
}


void VectorField::save_vtk(const std::string &filename)
{

    std::ofstream file;

    // Open file
    file.open(filename.c_str());

    // Write content
    file << "# vtk DataFile Version 2.0\n";
    file << "Comment goes here\n";
    file << "ASCII\n";
    file << "\n";
    file << "DATASET STRUCTURED_POINTS\n";
    file << "DIMENSIONS  " << ii << " " << jj  << " " << kk  << std::endl;
    file << "\n";
    file << "ORIGIN 0.000   0.000   0.000\n";
    file << "SPACING    "   << mesh->get_dx() <<" "<< mesh->get_dy()
                            << " " << mesh->get_dz() << "\n";
    file << "\n";
    file << "POINT_DATA " << ii * jj * kk << std::endl;
    file << "VECTORS vectors float\n";
    // file << "LOOKUP_TABLE default\n";
    file << "\n";

    
    for (int k = 0; k < kk; k++)
    {
        for (int j = 0; j < jj; j++)
        {
            for (int i = 0; i < ii; i++)
            {
                file << get(0,i,j,k) << " ";
                file << get(1,i,j,k) << " ";
                file << get(2,i,j,k) << " ";
                file << std::endl;
            }
        }
    }

    // Close file name
    file.close();
}

VectorField::VectorField(const Mesh &mesh0)
{
    mesh = &mesh0;
    ii = mesh->get_ii();
    jj = mesh->get_jj();
    kk = mesh->get_kk();
            
    scalar_component = new ScalarField[3];
    scalar_component[0].init(mesh0, 0.);
    scalar_component[1].init(mesh0, 0.);
    scalar_component[2].init(mesh0, 0.);
}

VectorField::VectorField(const Mesh &mesh0, long double val)
{
    mesh = &mesh0;
    ii = mesh->get_ii();
    jj = mesh->get_jj();
    kk = mesh->get_kk();
            
    scalar_component = new ScalarField[3];
    scalar_component[0].init(mesh0, val);
    scalar_component[1].init(mesh0, val);
    scalar_component[2].init(mesh0, val);
}

VectorField::VectorField(const Mesh &mesh0, const ScalarField & vx, const ScalarField & vy, const ScalarField & vz)
{
    mesh = &mesh0;
    ii = mesh->get_ii();
    jj = mesh->get_jj();
    kk = mesh->get_kk();
            
    scalar_component = new ScalarField[3];
    scalar_component[0] = vx;
    scalar_component[1] = vy;
    scalar_component[2] = vz;
}

void VectorField::set(const int coord, const int i, const int j, const int k, long double val)
{
    scalar_component[coord].set(i, j, k, val);
}

void VectorField::set_coordinate(const int coord, const ScalarField & v)
{
    scalar_component[coord] = v;
}

ScalarField & VectorField::get_coordinate(const int coord)
{
    return scalar_component[coord];
}

long double VectorField::get(const int coord, const int i, const int j, const int k) const
{
    return scalar_component[coord].get(i, j, k);        
}

const Mesh& VectorField::get_mesh() const
{
    return *mesh;
}

void VectorField::print() const
{
    scalar_component[0].print();
    scalar_component[1].print();
    scalar_component[2].print();  
}

void VectorField::save_components_vtk(const std::string &filename, const std::string &save_format)
{
    scalar_component[0].save( "x_" + filename, save_format);
    scalar_component[1].save( "y_" + filename, save_format);
    scalar_component[2].save( "z_" + filename, save_format);
}

void VectorField::save(const std::string &filename, const std::string &save_format)
{

    if ( save_format == "vtk" )
    {
        save_vtk( filename + ".vtk" );
    }
    else if ( save_format == "csv" )
    {
        save_csv( filename + ".csv" );
    }
    else
    {
        std::cout << "Invalid save format.";
        std::cout.flush();
        exit(EXIT_FAILURE);
    }

}

ScalarField VectorField::get_norm()
{
    ScalarField norm(*mesh, 0.);
            
    long double val;
    for (int k = 0; k < kk; k++)
    {
        for (int j = 0; j < jj; j++)
        {
            for (int i = 0; i < ii; i++)
            {
                val =   sqrt(   pow( get(0,i,j,k), 2) + 
                        pow( get(1,i,j,k), 2) +
                        pow( get(2,i,j,k), 2) );
                norm.set(i, j, k, val);
            }
        }
    }
    return norm;
}

ScalarField & VectorField::get_component ( int i )
{
    return scalar_component[i];
}

VectorField& VectorField::operator =(const VectorField & other)
{
    scalar_component[0] = (ScalarField) other.scalar_component[0];
    scalar_component[1] = (ScalarField) other.scalar_component[1];
    scalar_component[2] = (ScalarField) other.scalar_component[2];
    return *this;
}

long double VectorField::operator()(const int coord, const int i, const int j, const int k) const
{
    return get(coord, i, j, k);
}

VectorField VectorField::operator +=(const VectorField & other) const
{
    scalar_component[0] = scalar_component[0] + other.scalar_component[0];
    scalar_component[1] = scalar_component[1] + other.scalar_component[1];
    scalar_component[2] = scalar_component[2] + other.scalar_component[2];
    return *this;
}

VectorField VectorField::operator -=(const VectorField & other) const
{
    scalar_component[0] =  scalar_component[0] - other.scalar_component[0];
    scalar_component[1] =  scalar_component[1] - other.scalar_component[1];
    scalar_component[2] =  scalar_component[2] - other.scalar_component[2];
    return *this;
}

VectorField::~VectorField()
{
    if (scalar_component != NULL)
        delete[] scalar_component;
    // The mesh can be shared with other scalar field objects
    // ergo the mesh will not be deleted.
}


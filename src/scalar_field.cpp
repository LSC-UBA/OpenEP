#include "scalar_field.h"

void ScalarField::save_vtk(const std::string &filename)
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
    file << "DIMENSIONS  " << ii << " " << jj << " " << kk << std::endl;
    file << "\n";
    file << "ORIGIN    0.000   0.000   0.000\n";
    file << "SPACING    "   << mesh->get_dx() <<" "<< mesh->get_dy()
                            << " " << mesh->get_dz() << "\n";
    file << "\n";
    file << "POINT_DATA " << ii * jj * kk << std::endl;
    file << "SCALARS scalars float\n";
    file << "LOOKUP_TABLE default\n";
    file << "\n";

    for (int k = 0; k < kk; k++)
    {
        for (int j = 0; j < jj; j++)
        {
            for (int i = 0; i < ii; i++)
            {
                file << p_mat[i + j * ii + ii * jj * k] << " ";
            }
            file << std::endl;
        }
    }

    // Close file name
    file.close();
}

void ScalarField::save_csv(const std::string &filename)
{
    std::ofstream file;
    long double aux;
    // Open file
    file.open(filename.c_str());

    // Write content
    file << "X,Y,Z,I,J,K,S" << std::endl;
    for (int k = 0; k < kk; k++)
    {
        for (int j = 0; j < jj; j++)
        {
            for (int i = 0; i < ii; i++)
            {
                aux = p_mat[i + j * ii + ii * jj * k];
                file    << mesh->get_x(i) << "," 
                        << mesh->get_y(j) << ","
                        << mesh->get_z(k)
                        << "," << i << "," << j << "," << k 
                        << "," << aux  << std::endl;
            }
        }
    }
    // Close file name
    file.close();

}

ScalarField::ScalarField()
{
    mesh = NULL;
    p_mat = NULL;
    n = 0;
}

ScalarField::ScalarField(const Mesh &mesh0)
{
    mesh = &mesh0;
    ii = mesh->get_ii();
    jj = mesh->get_jj();
    kk = mesh->get_kk();
    n = ii * jj * kk;
    p_mat = new long double[n];
    init(0.);
}

ScalarField::ScalarField(const Mesh &mesh0, long double val)
{
    init(mesh0, val);
}

void ScalarField::init(const Mesh &mesh0, long double val)
{
    mesh = &mesh0;
    ii = mesh->get_ii();
    jj = mesh->get_jj();
    kk = mesh->get_kk();
    n = ii * jj * kk;
    p_mat = new long double[n];
    init(0.);
    for (int l = 0; l < n; l++)
        p_mat[l] = val;
}

void ScalarField::init(long double val)
{
    for (int l = 0; l < n; l++)
        p_mat[l] = val;
}

void ScalarField::set(const int i, const int j, const int k, long double val)
{
    p_mat[i + j * ii + ii * jj * k] = val;
}

long double ScalarField::get(const int i, const int j, const int k) const
{
    return p_mat[i + j * ii + ii * jj * k];

}

void ScalarField::save(const std::string &filename, const std::string &save_format)
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
        exit(EXIT_FAILURE);
    }

}

const Mesh& ScalarField::get_mesh() const
{
    return *mesh;
}

void ScalarField::print() const
{
    for (int i = 0; i < ii; i++)
    {
        for (int j = 0; j < jj; j++)
        {
            for (int k = 0; k < kk; k++)
            {
                std::cout   << mesh->get_x(i) << "," << mesh->get_y(j) << ","
                            << mesh->get_z(k) << ","
                            << p_mat[i + j * ii + ii * jj * k] << std::endl;
            }
        }
    }
}   

long double ScalarField::get_norm() const
{
    long double max = fabs(p_mat[0]);
    for (int l = 1; l < n; l++)
    {
        if (max < fabs(p_mat[l]))
            max = fabs(p_mat[l]);
    }
    return max;
}

// It assums that ii, jj and kk doesn't change
ScalarField& ScalarField::operator =(const ScalarField & other)
{

    mesh = &other.get_mesh();
    ii = mesh->get_ii();
    jj = mesh->get_jj();
    kk = mesh->get_kk();
    if (p_mat == NULL)
    {
        int n = ii * jj * kk;
        p_mat = new long double[n];
    }
    #pragma omp parallel for collapse(3) schedule(static)
    for (int i = 0; i < ii; i++)
    {
        for (int j = 0; j < jj; j++)
        {
            for (int k = 0; k < kk; k++)
            {
                p_mat[i + j * ii + ii * jj * k] = other(i, j, k);
            }
        }
    }
    return *this;
}

long double & ScalarField::operator()(const int i, const int j, const int k) const
{
    return p_mat[i+j*ii+ii*jj*k];
}

ScalarField ScalarField::operator +(const ScalarField & other) const
{
    ScalarField add(*mesh, 0);
    #pragma omp parallel for collapse(3) schedule(static)
    for (int i = 0; i < ii; i++)
    {
        for (int j = 0; j < jj; j++)
        {
            for (int k = 0; k < kk; k++)
            {
                add.set(i, j, k,
                        p_mat[i + j * ii + ii * jj * k]
                                - other(i, j, k));
            }
        }
    }
    return add;
}

ScalarField ScalarField::operator *(const ScalarField & other) const
{
    ScalarField mul(*mesh, 1);
    #pragma omp parallel for collapse(3) schedule(static)
    for (int i = 0; i < ii; i++)
    {
        for (int j = 0; j < jj; j++)
        {
            for (int k = 0; k < kk; k++)
            {
                mul.set(i, j, k,
                        p_mat[i + j * ii + ii * jj * k]
                                * other(i, j, k));
            }
        }
    }
    return mul;
}

ScalarField ScalarField::operator -(const ScalarField & other) const
{
    ScalarField diff(*mesh, 0);
    #pragma omp parallel for collapse(3) schedule(static)
    for (int i = 0; i < ii; i++)
    {
        for (int j = 0; j < jj; j++)
        {
            for (int k = 0; k < kk; k++)
            {
                diff.set(i, j, k,
                        p_mat[i + j * ii + ii * jj * k]
                                - other(i, j, k));
            }
        }
    }
    return diff;
}

ScalarField::~ScalarField()
{
    if (p_mat != NULL)
        delete[] p_mat;
    // The mesh can be shared with other scalar field objects
    // ergo the mesh will not be deleted.
}


#include "mesh.h"

/* Mesh class keeps information about geometry and domain discretization */

Mesh::Mesh( long double x_max0, long double y_max0, long double z_max0,
            long double electrode_length0, long double electrode_width0, 
            long double electrode_thickness0, long double gap_anode_cathode0,
            long double gap_elect_elect0, int no_electrodes0, int resolution0,
            std::string electrode_type0)

{

    /* Calc. anodes and cathodes positions */

    x_max = x_max0;
    y_max = y_max0;
    z_max = z_max0;

    electrode_length = electrode_length0;
    electrode_width = electrode_width0;
    electrode_thickness = electrode_thickness0;
    gap_anode_cathode = gap_anode_cathode0;
    gap_elect_elect = gap_elect_elect0; 
    resolution = resolution0;
    if(electrode_type0 == "plates"){
        no_electrodes = 1;
        e_type = PLATES;
    } 
    else if(electrode_type0 == "needles")
    {
        no_electrodes = no_electrodes0;
        e_type = NEEDLES;
    }
    else
    {
        std::cout << "Invalid type of electrode. Only 'plates' or 'needles' are valid."
        << std::endl;
        exit(EXIT_FAILURE);
    }

    dx = gap_anode_cathode / resolution;
    dy = dx;
    dz = dx;

    ii = int( x_max / dx );
    jj = int( y_max / dy );
    kk = int( z_max / dz );

    /*  No. of elect. nodes in z-axis */
    int nodes_long = electrode_length / dz + 1; 
    int nodes_width = electrode_width / dy;
    int nodes_thickness = electrode_thickness / dx;

    anode_number = nodes_width * nodes_long * nodes_thickness * no_electrodes; 
    cathode_number = nodes_width * nodes_long * nodes_thickness * no_electrodes;

    /* Anode positions */
    anodes = new int * [anode_number]; 
    /* Cathode positions */
    cathodes = new int * [cathode_number];

    int electrode_gap_nodes = gap_anode_cathode / dx + 1;

    for( int l = 0; l < anode_number; l++ )
        anodes[l] = new int[3];

    for( int l = 0; l < cathode_number; l++ )
        cathodes[l] = new int[3];

    int l = 0;
    int offset_anode_i = ii / 2. - electrode_gap_nodes / 2.; 
    int offset_cathode_i = ii / 2. + electrode_gap_nodes / 2. + 1;

    if (e_type == NEEDLES)
    {                              
        int nodes_elect_zone = nodes_width * no_electrodes 
                             + gap_elect_elect / dy 
                             * (no_electrodes - 1);
        int start_j = (jj - nodes_elect_zone) / 2;

        for(int ne = 0; ne < no_electrodes; ne++)
        {
            for(int i = 0; i < nodes_thickness; i++)
            {
                for(int j = 0; j < nodes_width; j++)
                {
                      for(int k = 0; k < nodes_long; k++)
                      {
                            anodes[l][0] = offset_anode_i - i;
                            anodes[l][1] = start_j + j;
                            anodes[l][2] = k;

                            cathodes[l][0] = offset_cathode_i + i;
                            cathodes[l][1] = start_j + j;
                            cathodes[l][2] = k;

                            l++;
                      }
                }
            }
            start_j += gap_elect_elect / dy + nodes_width;
        }
    }
    if (e_type == PLATES)
    {
        int offset_j = ( jj - nodes_width ) / 2.;

        for(int i = 0; i < nodes_thickness; i++)
        {
              for(int j = 0; j < nodes_width; j++)
              {
                    for(int k = 0; k < nodes_long; k++)
                    {
                          anodes[l][0] = offset_anode_i - i;
                          anodes[l][1] = offset_j + j;
                          anodes[l][2] = k;

                          cathodes[l][0] = offset_cathode_i + i;
                          cathodes[l][1] = offset_j + j;
                          cathodes[l][2] = k;

                          l++;
                    }
              }
        }
    }

    // Init. p_is_electrode and p_is_elec_iso_material

    int nn = ii * jj * kk;

    p_is_electrode = new bool [nn];
    p_is_elec_iso_material = new bool [nn];

    for ( int l = 0; l < nn; l++)
    {
        p_is_electrode[l] = false;
        p_is_elec_iso_material[l] = false;
    }

    int i, j, k;

    for( int l = 0; l < anode_number; l++ )
    {
        i = anodes[l][0];
        j = anodes[l][1];
        k = anodes[l][2];
        p_is_electrode[i + j * ii + ii * jj * k] = true;
    }

    for( int l = 0; l < cathode_number; l++ )
    {
        i = cathodes[l][0];
        j = cathodes[l][1];
        k = cathodes[l][2];
        p_is_electrode[i + j * ii + ii * jj * k] = true;
    }

    // Init. p_is_elec_iso_material if it is needed

}

long double Mesh::get_x_max() const
{
    return x_max;
}

long double Mesh::get_y_max() const
{
    return y_max;
}

long double Mesh::get_z_max() const
{
    return z_max;
}

long double Mesh::get_dx() const
{
    return dx;
}

long double Mesh::get_dy() const
{
    return dy;
}

long double Mesh::get_dz() const
{
    return dz;
}

int Mesh::get_ii() const
{
    return ii;
}

int Mesh::get_jj() const
{
    return jj;
}

int Mesh::get_kk() const
{
    return kk;
}

long double Mesh::get_x(const int i) const
{
    return dx * i;
}

long double Mesh::get_y(const int j) const
{
    return dy * j;
}

long double Mesh::get_z(const int k) const
{
    return dz * k;
}

int ** Mesh::get_anodes() const
{
    return anodes;
}

int ** Mesh::get_cathodes() const
{
    return cathodes;
}

int Mesh::get_anode_number() const
{
    return anode_number;
}

int Mesh::get_cathode_number() const
{
    return cathode_number;
}

bool Mesh::is_electrode(int i, int j, int k) const
{
    return p_is_electrode[i + j * ii + ii * jj * k];
}

bool Mesh::is_elec_isolating_material(int i, int j, int k) const
{
    return p_is_elec_iso_material[i + j * ii + ii * jj * k];
}

Mesh::~Mesh()
{
}

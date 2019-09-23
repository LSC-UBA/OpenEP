#ifndef MESH_H_
#define MESH_H_

#include <math.h> 
#include <iostream>
#include <string>
#include <cstdlib>

/* Mesh class keeps information about geometry and domain discretization */

class Mesh
{

private:

    /* Domain parameters */

    /* Maximum length of the x axis [m] */
    long double x_max; 
    /* Maximum length of the y axis [m] */
    long double y_max;
    /* Maximum length of the z axis [m] */
    long double z_max; 
    enum electrode_type {PLATES, NEEDLES};
    /* Type of electrode */
    electrode_type e_type; 
    /* Electrode length [m] */
    long double electrode_length;
    /* Electrode width [m] */ 
    long double electrode_width; 
    /* Electrode thickness [m] */
    long double electrode_thickness; 
    /* Gap between anode and cathode [m] */
    long double gap_anode_cathode; 
    /* Gap between anode-anode or cathode-cathode [m] */
    long double gap_elect_elect;
    /* No. of electrodes */
    int no_electrodes; 
    
    /* Numerical parameters */
    
    /* No. of nodes between the anode and cathode */
    int resolution;
    /* Space between nodes in x axis [m] */
    long double dx;
    /* Space between nodes in y axis [m] */
    long double dy;
    /* Space between nodes in z axis [m] */
    long double dz; 
    /* Number of divisions in x axis */
    int ii;
    /* Number of divisions in y axis */
    int jj; 
    /* Number of divisions in z axis */
    int kk; 
    /* Anodes index positions */
    int ** anodes; 
    /* No. of anodes */
    int anode_number;
    /* Cathodes index positions */ 
    int ** cathodes;
    /* No. of cathodes */
    int cathode_number; 
    /* Is electrode Mask */
    bool * p_is_electrode; 
    /* Is electrical isolating material mask */
    bool * p_is_elec_iso_material; 


public:

    Mesh(   long double x_max0, long double y_max0, long double z_max0,
            long double electrode_length0, long double electrode_width0, 
            long double electrode_thickness0, long double gap_anode_cathode0,
            long double gap_elect_elect0, int no_electrodes0, int resolution0,
            std::string electrode_type0);

    long double get_x_max() const;

    long double get_y_max() const;

    long double get_z_max() const;

    long double get_dx() const;

    long double get_dy() const;
    
    long double get_dz() const;

    int get_ii() const;

    int get_jj() const;
    
    int get_kk() const;

    long double get_x(const int i) const;

    long double get_y(const int j) const;

    long double get_z(const int k) const;
    
    int ** get_anodes() const;
    
    int ** get_cathodes() const;
    
    int get_anode_number() const;
    
    int get_cathode_number() const;

    bool is_electrode(int i, int j, int k) const;
    
    bool is_elec_isolating_material(int i, int j, int k) const;
    
    ~Mesh();
    
};

#endif


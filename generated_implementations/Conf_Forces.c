#ifdef _WIN32
#    define Conf_Forces_API __declspec(dllexport)
#else
#    define Conf_Forces_API
#endif
#include <math.h>
#include <stdint.h>
#include <stddef.h>

Conf_Forces_API void calc_Nodal_C(size_t num_elem_nodal,double CF_out[][3],double CF[][3],int64_t inverse[]){
    for (size_t i=0;i<num_elem_nodal;i++){
        for (size_t j=0;j<3;j++){
            CF_out[inverse[i]][j]+=CF[i][j];
        }
    }
}


void CPE4_dynamic(double *rst,double *coord,double rho,double *Element_U,double *Element_V,double *Element_A,double *S,double PENER,double SENER,double *res_0){
    double r = rst[0];
    double s = rst[1];
    double t = rst[2];
    double *U = Element_U;
    double *V = Element_V;
    double *A = Element_A;
    double S11 = S[0];
    double S22 = S[1];
    double S33 = S[2];
    double S12 = S[3];
    double S13 = S[4];
    double S23 = S[5];
    
    double x0 = (1.0/4.0)*r;
    double x1 = x0 + 1.0/4.0;
    double x2 = -x1;
    double x3 = x0 - 1.0/4.0;
    double x4 = -x3;
    double x5 = x1*coord[6] + x2*coord[3] + x3*coord[0] + x4*coord[9];
    double x6 = (1.0/4.0)*s;
    double x7 = x6 + 1.0/4.0;
    double x8 = -x7;
    double x9 = x6 - 1.0/4.0;
    double x10 = -x9;
    double x11 = x10*coord[4] + x7*coord[7] + x8*coord[10] + x9*coord[1];
    double x12 = x1*coord[7] + x2*coord[4] + x3*coord[1] + x4*coord[10];
    double x13 = x10*coord[3] + x7*coord[6] + x8*coord[9] + x9*coord[0];
    double x14 = -x11*x5 + x12*x13;
    double x15 = 1.0/x14;
    double x16 = x15*x3;
    double x17 = -x5;
    double x18 = x15*x9;
    double x19 = x13*x16 + x17*x18;
    double x20 = x1*U[7] + x2*U[4] + x3*U[1] + x4*U[10];
    double x21 = -x11;
    double x22 = x15*x21;
    double x23 = x20*x22;
    double x24 = x10*U[4] + x7*U[7] + x8*U[10] + x9*U[1];
    double x25 = x12*x15;
    double x26 = x24*x25;
    double x27 = -1.0*x23 - 1.0*x26;
    double x28 = x1*U[6] + x2*U[3] + x3*U[0] + x4*U[9];
    double x29 = x22*x28;
    double x30 = x10*U[3] + x7*U[6] + x8*U[9] + x9*U[0];
    double x31 = x25*x30;
    double x32 = 1.0*x29 + 1.0*x31 + 1.0;
    double x33 = S11*x27 + S12*x32;
    double x34 = x29 + x31;
    double x35 = S12*x27 + S22*x32;
    double x36 = x23 + x26;
    double x37 = S13*x27 + S23*x32;
    double x38 = x1*U[8] + x2*U[5] + x3*U[2] + x4*U[11];
    double x39 = x10*U[5] + x7*U[8] + x8*U[11] + x9*U[2];
    double x40 = x22*x38 + x25*x39;
    double x41 = -x33*x34 - x35*x36 - x37*x40;
    double x42 = x12*x18 + x16*x21;
    double x43 = x13*x15;
    double x44 = x20*x43;
    double x45 = x15*x17;
    double x46 = x24*x45;
    double x47 = 1.0*x44 + 1.0*x46 + 1.0;
    double x48 = x28*x43;
    double x49 = x30*x45;
    double x50 = -1.0*x48 - 1.0*x49;
    double x51 = S11*x47 + S12*x50;
    double x52 = S12*x47 + S22*x50;
    double x53 = S13*x47 + S23*x50;
    double x54 = -x0;
    double x55 = s*x0;
    double x56 = x55 - x6;
    double x57 = x54 + x56 + 1.0/4.0;
    double x58 = -x3 - x56;
    double x59 = -x54 - x55 - x9;
    double x60 = x1 + x55 + x6;
    double x61 = x57*V[0] + x58*V[9] + x59*V[3] + x60*V[6];
    double x62 = x57*V[1] + x58*V[10] + x59*V[4] + x60*V[7];
    double x63 = x57*V[2] + x58*V[11] + x59*V[5] + x60*V[8];
    double x64 = -1.0*PENER - 1.0*SENER + 0.5*rho*(x61*x61 + x62*x62 + x63*x63);
    double x65 = -x34*x51 - x36*x52 - x40*x53 - x64;
    double x66 = x57*A[0] + x58*A[9] + x59*A[3] + x60*A[6];
    double x67 = x57*A[1] + x58*A[10] + x59*A[4] + x60*A[7];
    double x68 = x57*A[2] + x58*A[11] + x59*A[5] + x60*A[8];
    double x69 = x1*V[6] + x2*V[3] + x3*V[0] + x4*V[9];
    double x70 = x10*V[3] + x7*V[6] + x8*V[9] + x9*V[0];
    double x71 = x1*V[7] + x2*V[4] + x3*V[1] + x4*V[10];
    double x72 = x10*V[4] + x7*V[7] + x8*V[10] + x9*V[1];
    double x73 = x1*V[8] + x2*V[5] + x3*V[2] + x4*V[11];
    double x74 = x10*V[5] + x7*V[8] + x8*V[11] + x9*V[2];
    double x75 = x34*x66 + x36*x67 + x40*x68 + x61*(x22*x69 + x25*x70) + x62*(x22*x71 + x25*x72) + x63*(x22*x73 + x25*x74);
    double x76 = rho*x57;
    double x77 = x48 + x49;
    double x78 = x44 + x46;
    double x79 = x38*x43 + x39*x45;
    double x80 = -x33*x77 - x35*x78 - x37*x79 - x64;
    double x81 = -x51*x77 - x52*x78 - x53*x79;
    double x82 = x61*(x43*x69 + x45*x70) + x62*(x43*x71 + x45*x72) + x63*(x43*x73 + x45*x74) + x66*x77 + x67*x78 + x68*x79;
    double x83 = x10*x45 + x2*x43;
    double x84 = x10*x25 + x2*x22;
    double x85 = rho*x59;
    double x86 = x1*x43 + x45*x7;
    double x87 = x1*x22 + x25*x7;
    double x88 = rho*x60;
    double x89 = x22*x4 + x25*x8;
    double x90 = x4*x43 + x45*x8;
    double x91 = rho*x58;
    
    res_0[0] = x14*(x19*x41 + x42*x65 - x75*x76);
    res_0[1] = x14*(x19*x80 + x42*x81 - x76*x82);
    res_0[2] = 0;
    res_0[3] = x14*(x41*x83 + x65*x84 - x75*x85);
    res_0[4] = x14*(x80*x83 + x81*x84 - x82*x85);
    res_0[5] = 0;
    res_0[6] = x14*(x41*x86 + x65*x87 - x75*x88);
    res_0[7] = x14*(x80*x86 + x81*x87 - x82*x88);
    res_0[8] = 0;
    res_0[9] = x14*(x41*x90 + x65*x89 - x75*x91);
    res_0[10] = x14*(x80*x90 + x81*x89 - x82*x91);
    res_0[11] = 0;
}

Conf_Forces_API void Integration_CPE4_dynamic(size_t num_elem,double Coords[][4][3],double *rho,double Element_U[][4][3],double Element_V[][4][3],double Element_A[][4][3],double S[][4][6],
    double PENER[][4],double SENER[][4],double Conf_Force[][4][3]){
    
    //Inputs:
    //    Coords.........Element_Nodal_Coordinates  [num_elem][numNodes,3]
    //    Element_U......Element_Nodal_Displacements[num_elem][numNodes,3]
    //    Element_V......Element_Nodal_Velocity     [num_elem][numNodes,3]
    //    Element_A......Element_Nodal_Acceleration [num_elem][numNodes,3]
    //    S..............Stress at integration points [num_elem][num_int_points,6]
    //    PENER..........Plastic strain energy [num_elem][num_int_points]
    //    SENER..........Elastic strain energy [num_elem][num_int_points]
    //    
    double weights[4]={1.0,1.0,1.0,1.0,};
    double int_points[4][3]={
        {-0.5773502691896258,-0.5773502691896258,0.0,},
        {0.5773502691896258,-0.5773502691896258,0.0,},
        {-0.5773502691896258,0.5773502691896258,0.0,},
        {0.5773502691896258,0.5773502691896258,0.0,}};

    #pragma omp parallel for
    for (size_t i=0;i<num_elem;i++){
        double TMP[4][3];
        for (size_t j=0;j<4;j++){
            CPE4_dynamic((double*) int_points[j],(double*) Coords[i],rho[i],(double*) Element_U[i],(double*) Element_V[i],(double*) Element_A[i],(double*) S[i][j],PENER[i][j],SENER[i][j],(double*) TMP);
            for (size_t num_node=0;num_node<4;num_node++){
                for (size_t ii=0;ii<3;ii++){
                    Conf_Force[i][num_node][ii]+=TMP[num_node][ii]*weights[j];
                }
            }
        }
    }
}

void CPE4_static_mbf(double *rst,double *coord,double *Element_U,double *S,double PENER,double SENER,double *res_0){
    double r = rst[0];
    double s = rst[1];
    double t = rst[2];
    double *U = Element_U;
    double S11 = S[0];
    double S22 = S[1];
    double S33 = S[2];
    double S12 = S[3];
    double S13 = S[4];
    double S23 = S[5];
    
    double x0 = (1.0/4.0)*s;
    double x1 = x0 + 1.0/4.0;
    double x2 = -x1;
    double x3 = x0 - 1.0/4.0;
    double x4 = -x3;
    double x5 = x1*coord[6] + x2*coord[9] + x3*coord[0] + x4*coord[3];
    double x6 = (1.0/4.0)*r;
    double x7 = x6 - 1.0/4.0;
    double x8 = x6 + 1.0/4.0;
    double x9 = -x8;
    double x10 = -x7;
    double x11 = x10*coord[9] + x7*coord[0] + x8*coord[6] + x9*coord[3];
    double x12 = x1*coord[7] + x2*coord[10] + x3*coord[1] + x4*coord[4];
    double x13 = x10*coord[10] + x7*coord[1] + x8*coord[7] + x9*coord[4];
    double x14 = -x11*x12 + x13*x5;
    double x15 = 1.0/x14;
    double x16 = x15*x7;
    double x17 = -x11;
    double x18 = x15*x3;
    double x19 = x16*x5 + x17*x18;
    double x20 = x10*U[10] + x7*U[1] + x8*U[7] + x9*U[4];
    double x21 = -x12;
    double x22 = x15*x21;
    double x23 = x20*x22;
    double x24 = x1*U[7] + x2*U[10] + x3*U[1] + x4*U[4];
    double x25 = x13*x15;
    double x26 = x24*x25;
    double x27 = -1.0*x23 - 1.0*x26;
    double x28 = x10*U[9] + x7*U[0] + x8*U[6] + x9*U[3];
    double x29 = x22*x28;
    double x30 = x1*U[6] + x2*U[9] + x3*U[0] + x4*U[3];
    double x31 = x25*x30;
    double x32 = 1.0*x29 + 1.0*x31 + 1.0;
    double x33 = S11*x27 + S12*x32;
    double x34 = x29 + x31 + 1.0;
    double x35 = S12*x27 + S22*x32;
    double x36 = x23 + x26;
    double x37 = S13*x27;
    double x38 = S23*x32;
    double x39 = x37 + x38;
    double x40 = x10*U[11] + x7*U[2] + x8*U[8] + x9*U[5];
    double x41 = x1*U[8] + x2*U[11] + x3*U[2] + x4*U[5];
    double x42 = x22*x40 + x25*x41;
    double x43 = -x33*x34 - x35*x36 - x39*x42;
    double x44 = x13*x18 + x16*x21;
    double x45 = x15*x5;
    double x46 = x20*x45;
    double x47 = x15*x17;
    double x48 = x24*x47;
    double x49 = 1.0*x46 + 1.0*x48 + 1.0;
    double x50 = x28*x45;
    double x51 = x30*x47;
    double x52 = -1.0*x50 - 1.0*x51;
    double x53 = S11*x49 + S12*x52;
    double x54 = S12*x49 + S22*x52;
    double x55 = S13*x49;
    double x56 = S23*x52;
    double x57 = x55 + x56;
    double x58 = -1.0*PENER - 1.0*SENER;
    double x59 = -x34*x53 - x36*x54 - x42*x57 - x58;
    double x60 = x50 + x51;
    double x61 = x46 + x48 + 1.0;
    double x62 = x40*x45 + x41*x47;
    double x63 = -x33*x60 - x35*x61 - x39*x62 - x58;
    double x64 = -x53*x60 - x54*x61 - x57*x62;
    double x65 = -1.0*x37 - 1.0*x38;
    double x66 = -1.0*x55 - 1.0*x56;
    double x67 = x4*x47 + x45*x9;
    double x68 = x22*x9 + x25*x4;
    double x69 = x1*x47 + x45*x8;
    double x70 = x1*x25 + x22*x8;
    double x71 = x10*x22 + x2*x25;
    double x72 = x10*x45 + x2*x47;
    
    res_0[0] = x14*(x19*x43 + x44*x59);
    res_0[1] = x14*(x19*x63 + x44*x64);
    res_0[2] = x14*(x19*x65 + x44*x66);
    res_0[3] = x14*(x43*x67 + x59*x68);
    res_0[4] = x14*(x63*x67 + x64*x68);
    res_0[5] = x14*(x65*x67 + x66*x68);
    res_0[6] = x14*(x43*x69 + x59*x70);
    res_0[7] = x14*(x63*x69 + x64*x70);
    res_0[8] = x14*(x65*x69 + x66*x70);
    res_0[9] = x14*(x43*x72 + x59*x71);
    res_0[10] = x14*(x63*x72 + x64*x71);
    res_0[11] = x14*(x65*x72 + x66*x71);
}

Conf_Forces_API void Integration_CPE4_static_mbf(size_t num_elem,double Coords[][4][3],double Element_U[][4][3],double S[][4][6],
    double PENER[][4],double SENER[][4],double Conf_Force[][4][3]){
    
    //Inputs:
    //    Coords.........Element_Nodal_Coordinates [num_elem][numNodes,3]
    //    Element_U......Element_Nodal_Displacements [num_elem][numNodes,3]
    //    S..............Stress at integration points [num_elem][num_int_points,6]
    //    PENER..........Plastic strain energy [num_elem][num_int_points]
    //    SENER..........Elastic strain energy [num_elem][num_int_points]
    //    
    double weights[4]={1.0,1.0,1.0,1.0,};
    double int_points[4][3]={
        {-0.5773502691896258,-0.5773502691896258,0.0,},
        {0.5773502691896258,-0.5773502691896258,0.0,},
        {-0.5773502691896258,0.5773502691896258,0.0,},
        {0.5773502691896258,0.5773502691896258,0.0,}};

    #pragma omp parallel for
    for (size_t i=0;i<num_elem;i++){
        double TMP[4][3];
        for (size_t j=0;j<4;j++){
            CPE4_static_mbf((double*) int_points[j],(double*) Coords[i],(double*) Element_U[i],(double*) S[i][j],PENER[i][j],SENER[i][j],(double*) TMP);
            for (size_t num_node=0;num_node<4;num_node++){
                for (size_t ii=0;ii<3;ii++){
                    Conf_Force[i][num_node][ii]+=TMP[num_node][ii]*weights[j];
                }
            }
        }
    }
}

void CPE4_static_dbf(double *rst,double *coord,double *Element_U,double *S,double PENER,double SENER,double *res_0){
    double r = rst[0];
    double s = rst[1];
    double t = rst[2];
    double *U = Element_U;
    double S11 = S[0];
    double S22 = S[1];
    double S33 = S[2];
    double S12 = S[3];
    double S13 = S[4];
    double S23 = S[5];
    
    double x0 = (1.0/4.0)*s;
    double x1 = x0 + 1.0/4.0;
    double x2 = -x1;
    double x3 = x0 - 1.0/4.0;
    double x4 = -x3;
    double x5 = x1*coord[6] + x2*coord[9] + x3*coord[0] + x4*coord[3];
    double x6 = (1.0/4.0)*r;
    double x7 = x6 - 1.0/4.0;
    double x8 = x6 + 1.0/4.0;
    double x9 = -x8;
    double x10 = -x7;
    double x11 = x10*coord[9] + x7*coord[0] + x8*coord[6] + x9*coord[3];
    double x12 = x1*coord[7] + x2*coord[10] + x3*coord[1] + x4*coord[4];
    double x13 = x10*coord[10] + x7*coord[1] + x8*coord[7] + x9*coord[4];
    double x14 = -x11*x12 + x13*x5;
    double x15 = 1.0/x14;
    double x16 = x15*x7;
    double x17 = -x11;
    double x18 = x15*x3;
    double x19 = x16*x5 + x17*x18;
    double x20 = x10*U[10] + x7*U[1] + x8*U[7] + x9*U[4];
    double x21 = -x12;
    double x22 = x15*x21;
    double x23 = x20*x22;
    double x24 = x1*U[7] + x2*U[10] + x3*U[1] + x4*U[4];
    double x25 = x13*x15;
    double x26 = x24*x25;
    double x27 = -1.0*x23 - 1.0*x26;
    double x28 = x10*U[9] + x7*U[0] + x8*U[6] + x9*U[3];
    double x29 = x22*x28;
    double x30 = x1*U[6] + x2*U[9] + x3*U[0] + x4*U[3];
    double x31 = x25*x30;
    double x32 = 1.0*x29 + 1.0*x31 + 1.0;
    double x33 = S11*x27 + S12*x32;
    double x34 = x29 + x31;
    double x35 = S12*x27 + S22*x32;
    double x36 = x23 + x26;
    double x37 = S13*x27 + S23*x32;
    double x38 = x10*U[11] + x7*U[2] + x8*U[8] + x9*U[5];
    double x39 = x1*U[8] + x2*U[11] + x3*U[2] + x4*U[5];
    double x40 = x22*x38 + x25*x39;
    double x41 = -x33*x34 - x35*x36 - x37*x40;
    double x42 = x13*x18 + x16*x21;
    double x43 = x15*x5;
    double x44 = x20*x43;
    double x45 = x15*x17;
    double x46 = x24*x45;
    double x47 = 1.0*x44 + 1.0*x46 + 1.0;
    double x48 = x28*x43;
    double x49 = x30*x45;
    double x50 = -1.0*x48 - 1.0*x49;
    double x51 = S11*x47 + S12*x50;
    double x52 = S12*x47 + S22*x50;
    double x53 = S13*x47 + S23*x50;
    double x54 = -1.0*PENER - 1.0*SENER;
    double x55 = -x34*x51 - x36*x52 - x40*x53 - x54;
    double x56 = x48 + x49;
    double x57 = x44 + x46;
    double x58 = x38*x43 + x39*x45;
    double x59 = -x33*x56 - x35*x57 - x37*x58 - x54;
    double x60 = -x51*x56 - x52*x57 - x53*x58;
    double x61 = x4*x45 + x43*x9;
    double x62 = x22*x9 + x25*x4;
    double x63 = x1*x45 + x43*x8;
    double x64 = x1*x25 + x22*x8;
    double x65 = x10*x22 + x2*x25;
    double x66 = x10*x43 + x2*x45;
    
    res_0[0] = x14*(x19*x41 + x42*x55);
    res_0[1] = x14*(x19*x59 + x42*x60);
    res_0[2] = 0;
    res_0[3] = x14*(x41*x61 + x55*x62);
    res_0[4] = x14*(x59*x61 + x60*x62);
    res_0[5] = 0;
    res_0[6] = x14*(x41*x63 + x55*x64);
    res_0[7] = x14*(x59*x63 + x60*x64);
    res_0[8] = 0;
    res_0[9] = x14*(x41*x66 + x55*x65);
    res_0[10] = x14*(x59*x66 + x60*x65);
    res_0[11] = 0;
}

Conf_Forces_API void Integration_CPE4_static_dbf(size_t num_elem,double Coords[][4][3],double Element_U[][4][3],double S[][4][6],
    double PENER[][4],double SENER[][4],double Conf_Force[][4][3]){
    
    //Inputs:
    //    Coords.........Element_Nodal_Coordinates [num_elem][numNodes,3]
    //    Element_U......Element_Nodal_Displacements [num_elem][numNodes,3]
    //    S..............Stress at integration points [num_elem][num_int_points,6]
    //    PENER..........Plastic strain energy [num_elem][num_int_points]
    //    SENER..........Elastic strain energy [num_elem][num_int_points]
    //    
    double weights[4]={1.0,1.0,1.0,1.0,};
    double int_points[4][3]={
        {-0.5773502691896258,-0.5773502691896258,0.0,},
        {0.5773502691896258,-0.5773502691896258,0.0,},
        {-0.5773502691896258,0.5773502691896258,0.0,},
        {0.5773502691896258,0.5773502691896258,0.0,}};

    #pragma omp parallel for
    for (size_t i=0;i<num_elem;i++){
        double TMP[4][3];
        for (size_t j=0;j<4;j++){
            CPE4_static_dbf((double*) int_points[j],(double*) Coords[i],(double*) Element_U[i],(double*) S[i][j],PENER[i][j],SENER[i][j],(double*) TMP);
            for (size_t num_node=0;num_node<4;num_node++){
                for (size_t ii=0;ii<3;ii++){
                    Conf_Force[i][num_node][ii]+=TMP[num_node][ii]*weights[j];
                }
            }
        }
    }
}

void CPE4R_dynamic(double *rst,double *coord,double rho,double *Element_U,double *Element_V,double *Element_A,double *S,double PENER,double SENER,double *res_0){
    double r = rst[0];
    double s = rst[1];
    double t = rst[2];
    double *U = Element_U;
    double *V = Element_V;
    double *A = Element_A;
    double S11 = S[0];
    double S22 = S[1];
    double S33 = S[2];
    double S12 = S[3];
    double S13 = S[4];
    double S23 = S[5];
    
    double x0 = (1.0/4.0)*r;
    double x1 = x0 + 1.0/4.0;
    double x2 = -x1;
    double x3 = x0 - 1.0/4.0;
    double x4 = -x3;
    double x5 = x1*coord[6] + x2*coord[3] + x3*coord[0] + x4*coord[9];
    double x6 = (1.0/4.0)*s;
    double x7 = x6 + 1.0/4.0;
    double x8 = -x7;
    double x9 = x6 - 1.0/4.0;
    double x10 = -x9;
    double x11 = x10*coord[4] + x7*coord[7] + x8*coord[10] + x9*coord[1];
    double x12 = x1*coord[7] + x2*coord[4] + x3*coord[1] + x4*coord[10];
    double x13 = x10*coord[3] + x7*coord[6] + x8*coord[9] + x9*coord[0];
    double x14 = -x11*x5 + x12*x13;
    double x15 = 1.0/x14;
    double x16 = x15*x3;
    double x17 = -x5;
    double x18 = x15*x9;
    double x19 = x13*x16 + x17*x18;
    double x20 = x1*U[7] + x2*U[4] + x3*U[1] + x4*U[10];
    double x21 = -x11;
    double x22 = x15*x21;
    double x23 = x20*x22;
    double x24 = x10*U[4] + x7*U[7] + x8*U[10] + x9*U[1];
    double x25 = x12*x15;
    double x26 = x24*x25;
    double x27 = -1.0*x23 - 1.0*x26;
    double x28 = x1*U[6] + x2*U[3] + x3*U[0] + x4*U[9];
    double x29 = x22*x28;
    double x30 = x10*U[3] + x7*U[6] + x8*U[9] + x9*U[0];
    double x31 = x25*x30;
    double x32 = 1.0*x29 + 1.0*x31 + 1.0;
    double x33 = S11*x27 + S12*x32;
    double x34 = x29 + x31;
    double x35 = S12*x27 + S22*x32;
    double x36 = x23 + x26;
    double x37 = S13*x27 + S23*x32;
    double x38 = x1*U[8] + x2*U[5] + x3*U[2] + x4*U[11];
    double x39 = x10*U[5] + x7*U[8] + x8*U[11] + x9*U[2];
    double x40 = x22*x38 + x25*x39;
    double x41 = -x33*x34 - x35*x36 - x37*x40;
    double x42 = x12*x18 + x16*x21;
    double x43 = x13*x15;
    double x44 = x20*x43;
    double x45 = x15*x17;
    double x46 = x24*x45;
    double x47 = 1.0*x44 + 1.0*x46 + 1.0;
    double x48 = x28*x43;
    double x49 = x30*x45;
    double x50 = -1.0*x48 - 1.0*x49;
    double x51 = S11*x47 + S12*x50;
    double x52 = S12*x47 + S22*x50;
    double x53 = S13*x47 + S23*x50;
    double x54 = -x0;
    double x55 = s*x0;
    double x56 = x55 - x6;
    double x57 = x54 + x56 + 1.0/4.0;
    double x58 = -x3 - x56;
    double x59 = -x54 - x55 - x9;
    double x60 = x1 + x55 + x6;
    double x61 = x57*V[0] + x58*V[9] + x59*V[3] + x60*V[6];
    double x62 = x57*V[1] + x58*V[10] + x59*V[4] + x60*V[7];
    double x63 = x57*V[2] + x58*V[11] + x59*V[5] + x60*V[8];
    double x64 = -1.0*PENER - 1.0*SENER + 0.5*rho*(x61*x61 + x62*x62 + x63*x63);
    double x65 = -x34*x51 - x36*x52 - x40*x53 - x64;
    double x66 = x57*A[0] + x58*A[9] + x59*A[3] + x60*A[6];
    double x67 = x57*A[1] + x58*A[10] + x59*A[4] + x60*A[7];
    double x68 = x57*A[2] + x58*A[11] + x59*A[5] + x60*A[8];
    double x69 = x1*V[6] + x2*V[3] + x3*V[0] + x4*V[9];
    double x70 = x10*V[3] + x7*V[6] + x8*V[9] + x9*V[0];
    double x71 = x1*V[7] + x2*V[4] + x3*V[1] + x4*V[10];
    double x72 = x10*V[4] + x7*V[7] + x8*V[10] + x9*V[1];
    double x73 = x1*V[8] + x2*V[5] + x3*V[2] + x4*V[11];
    double x74 = x10*V[5] + x7*V[8] + x8*V[11] + x9*V[2];
    double x75 = x34*x66 + x36*x67 + x40*x68 + x61*(x22*x69 + x25*x70) + x62*(x22*x71 + x25*x72) + x63*(x22*x73 + x25*x74);
    double x76 = rho*x57;
    double x77 = x48 + x49;
    double x78 = x44 + x46;
    double x79 = x38*x43 + x39*x45;
    double x80 = -x33*x77 - x35*x78 - x37*x79 - x64;
    double x81 = -x51*x77 - x52*x78 - x53*x79;
    double x82 = x61*(x43*x69 + x45*x70) + x62*(x43*x71 + x45*x72) + x63*(x43*x73 + x45*x74) + x66*x77 + x67*x78 + x68*x79;
    double x83 = x10*x45 + x2*x43;
    double x84 = x10*x25 + x2*x22;
    double x85 = rho*x59;
    double x86 = x1*x43 + x45*x7;
    double x87 = x1*x22 + x25*x7;
    double x88 = rho*x60;
    double x89 = x22*x4 + x25*x8;
    double x90 = x4*x43 + x45*x8;
    double x91 = rho*x58;
    
    res_0[0] = x14*(x19*x41 + x42*x65 - x75*x76);
    res_0[1] = x14*(x19*x80 + x42*x81 - x76*x82);
    res_0[2] = 0;
    res_0[3] = x14*(x41*x83 + x65*x84 - x75*x85);
    res_0[4] = x14*(x80*x83 + x81*x84 - x82*x85);
    res_0[5] = 0;
    res_0[6] = x14*(x41*x86 + x65*x87 - x75*x88);
    res_0[7] = x14*(x80*x86 + x81*x87 - x82*x88);
    res_0[8] = 0;
    res_0[9] = x14*(x41*x90 + x65*x89 - x75*x91);
    res_0[10] = x14*(x80*x90 + x81*x89 - x82*x91);
    res_0[11] = 0;
}

Conf_Forces_API void Integration_CPE4R_dynamic(size_t num_elem,double Coords[][4][3],double *rho,double Element_U[][4][3],double Element_V[][4][3],double Element_A[][4][3],double S[][1][6],
    double PENER[][1],double SENER[][1],double Conf_Force[][4][3]){
    
    //Inputs:
    //    Coords.........Element_Nodal_Coordinates  [num_elem][numNodes,3]
    //    Element_U......Element_Nodal_Displacements[num_elem][numNodes,3]
    //    Element_V......Element_Nodal_Velocity     [num_elem][numNodes,3]
    //    Element_A......Element_Nodal_Acceleration [num_elem][numNodes,3]
    //    S..............Stress at integration points [num_elem][num_int_points,6]
    //    PENER..........Plastic strain energy [num_elem][num_int_points]
    //    SENER..........Elastic strain energy [num_elem][num_int_points]
    //    
    double weights[1]={4.0,};
    double int_points[1][3]={
        {0.0,0.0,0.0,}};

    #pragma omp parallel for
    for (size_t i=0;i<num_elem;i++){
        double TMP[4][3];
        for (size_t j=0;j<1;j++){
            CPE4R_dynamic((double*) int_points[j],(double*) Coords[i],rho[i],(double*) Element_U[i],(double*) Element_V[i],(double*) Element_A[i],(double*) S[i][j],PENER[i][j],SENER[i][j],(double*) TMP);
            for (size_t num_node=0;num_node<4;num_node++){
                for (size_t ii=0;ii<3;ii++){
                    Conf_Force[i][num_node][ii]+=TMP[num_node][ii]*weights[j];
                }
            }
        }
    }
}

void CPE4R_static_mbf(double *rst,double *coord,double *Element_U,double *S,double PENER,double SENER,double *res_0){
    double r = rst[0];
    double s = rst[1];
    double t = rst[2];
    double *U = Element_U;
    double S11 = S[0];
    double S22 = S[1];
    double S33 = S[2];
    double S12 = S[3];
    double S13 = S[4];
    double S23 = S[5];
    
    double x0 = (1.0/4.0)*s;
    double x1 = x0 + 1.0/4.0;
    double x2 = -x1;
    double x3 = x0 - 1.0/4.0;
    double x4 = -x3;
    double x5 = x1*coord[6] + x2*coord[9] + x3*coord[0] + x4*coord[3];
    double x6 = (1.0/4.0)*r;
    double x7 = x6 - 1.0/4.0;
    double x8 = x6 + 1.0/4.0;
    double x9 = -x8;
    double x10 = -x7;
    double x11 = x10*coord[9] + x7*coord[0] + x8*coord[6] + x9*coord[3];
    double x12 = x1*coord[7] + x2*coord[10] + x3*coord[1] + x4*coord[4];
    double x13 = x10*coord[10] + x7*coord[1] + x8*coord[7] + x9*coord[4];
    double x14 = -x11*x12 + x13*x5;
    double x15 = 1.0/x14;
    double x16 = x15*x7;
    double x17 = -x11;
    double x18 = x15*x3;
    double x19 = x16*x5 + x17*x18;
    double x20 = x10*U[10] + x7*U[1] + x8*U[7] + x9*U[4];
    double x21 = -x12;
    double x22 = x15*x21;
    double x23 = x20*x22;
    double x24 = x1*U[7] + x2*U[10] + x3*U[1] + x4*U[4];
    double x25 = x13*x15;
    double x26 = x24*x25;
    double x27 = -1.0*x23 - 1.0*x26;
    double x28 = x10*U[9] + x7*U[0] + x8*U[6] + x9*U[3];
    double x29 = x22*x28;
    double x30 = x1*U[6] + x2*U[9] + x3*U[0] + x4*U[3];
    double x31 = x25*x30;
    double x32 = 1.0*x29 + 1.0*x31 + 1.0;
    double x33 = S11*x27 + S12*x32;
    double x34 = x29 + x31 + 1.0;
    double x35 = S12*x27 + S22*x32;
    double x36 = x23 + x26;
    double x37 = S13*x27;
    double x38 = S23*x32;
    double x39 = x37 + x38;
    double x40 = x10*U[11] + x7*U[2] + x8*U[8] + x9*U[5];
    double x41 = x1*U[8] + x2*U[11] + x3*U[2] + x4*U[5];
    double x42 = x22*x40 + x25*x41;
    double x43 = -x33*x34 - x35*x36 - x39*x42;
    double x44 = x13*x18 + x16*x21;
    double x45 = x15*x5;
    double x46 = x20*x45;
    double x47 = x15*x17;
    double x48 = x24*x47;
    double x49 = 1.0*x46 + 1.0*x48 + 1.0;
    double x50 = x28*x45;
    double x51 = x30*x47;
    double x52 = -1.0*x50 - 1.0*x51;
    double x53 = S11*x49 + S12*x52;
    double x54 = S12*x49 + S22*x52;
    double x55 = S13*x49;
    double x56 = S23*x52;
    double x57 = x55 + x56;
    double x58 = -1.0*PENER - 1.0*SENER;
    double x59 = -x34*x53 - x36*x54 - x42*x57 - x58;
    double x60 = x50 + x51;
    double x61 = x46 + x48 + 1.0;
    double x62 = x40*x45 + x41*x47;
    double x63 = -x33*x60 - x35*x61 - x39*x62 - x58;
    double x64 = -x53*x60 - x54*x61 - x57*x62;
    double x65 = -1.0*x37 - 1.0*x38;
    double x66 = -1.0*x55 - 1.0*x56;
    double x67 = x4*x47 + x45*x9;
    double x68 = x22*x9 + x25*x4;
    double x69 = x1*x47 + x45*x8;
    double x70 = x1*x25 + x22*x8;
    double x71 = x10*x22 + x2*x25;
    double x72 = x10*x45 + x2*x47;
    
    res_0[0] = x14*(x19*x43 + x44*x59);
    res_0[1] = x14*(x19*x63 + x44*x64);
    res_0[2] = x14*(x19*x65 + x44*x66);
    res_0[3] = x14*(x43*x67 + x59*x68);
    res_0[4] = x14*(x63*x67 + x64*x68);
    res_0[5] = x14*(x65*x67 + x66*x68);
    res_0[6] = x14*(x43*x69 + x59*x70);
    res_0[7] = x14*(x63*x69 + x64*x70);
    res_0[8] = x14*(x65*x69 + x66*x70);
    res_0[9] = x14*(x43*x72 + x59*x71);
    res_0[10] = x14*(x63*x72 + x64*x71);
    res_0[11] = x14*(x65*x72 + x66*x71);
}

Conf_Forces_API void Integration_CPE4R_static_mbf(size_t num_elem,double Coords[][4][3],double Element_U[][4][3],double S[][1][6],
    double PENER[][1],double SENER[][1],double Conf_Force[][4][3]){
    
    //Inputs:
    //    Coords.........Element_Nodal_Coordinates [num_elem][numNodes,3]
    //    Element_U......Element_Nodal_Displacements [num_elem][numNodes,3]
    //    S..............Stress at integration points [num_elem][num_int_points,6]
    //    PENER..........Plastic strain energy [num_elem][num_int_points]
    //    SENER..........Elastic strain energy [num_elem][num_int_points]
    //    
    double weights[1]={4.0,};
    double int_points[1][3]={
        {0.0,0.0,0.0,}};

    #pragma omp parallel for
    for (size_t i=0;i<num_elem;i++){
        double TMP[4][3];
        for (size_t j=0;j<1;j++){
            CPE4R_static_mbf((double*) int_points[j],(double*) Coords[i],(double*) Element_U[i],(double*) S[i][j],PENER[i][j],SENER[i][j],(double*) TMP);
            for (size_t num_node=0;num_node<4;num_node++){
                for (size_t ii=0;ii<3;ii++){
                    Conf_Force[i][num_node][ii]+=TMP[num_node][ii]*weights[j];
                }
            }
        }
    }
}

void CPE4R_static_dbf(double *rst,double *coord,double *Element_U,double *S,double PENER,double SENER,double *res_0){
    double r = rst[0];
    double s = rst[1];
    double t = rst[2];
    double *U = Element_U;
    double S11 = S[0];
    double S22 = S[1];
    double S33 = S[2];
    double S12 = S[3];
    double S13 = S[4];
    double S23 = S[5];
    
    double x0 = (1.0/4.0)*s;
    double x1 = x0 + 1.0/4.0;
    double x2 = -x1;
    double x3 = x0 - 1.0/4.0;
    double x4 = -x3;
    double x5 = x1*coord[6] + x2*coord[9] + x3*coord[0] + x4*coord[3];
    double x6 = (1.0/4.0)*r;
    double x7 = x6 - 1.0/4.0;
    double x8 = x6 + 1.0/4.0;
    double x9 = -x8;
    double x10 = -x7;
    double x11 = x10*coord[9] + x7*coord[0] + x8*coord[6] + x9*coord[3];
    double x12 = x1*coord[7] + x2*coord[10] + x3*coord[1] + x4*coord[4];
    double x13 = x10*coord[10] + x7*coord[1] + x8*coord[7] + x9*coord[4];
    double x14 = -x11*x12 + x13*x5;
    double x15 = 1.0/x14;
    double x16 = x15*x7;
    double x17 = -x11;
    double x18 = x15*x3;
    double x19 = x16*x5 + x17*x18;
    double x20 = x10*U[10] + x7*U[1] + x8*U[7] + x9*U[4];
    double x21 = -x12;
    double x22 = x15*x21;
    double x23 = x20*x22;
    double x24 = x1*U[7] + x2*U[10] + x3*U[1] + x4*U[4];
    double x25 = x13*x15;
    double x26 = x24*x25;
    double x27 = -1.0*x23 - 1.0*x26;
    double x28 = x10*U[9] + x7*U[0] + x8*U[6] + x9*U[3];
    double x29 = x22*x28;
    double x30 = x1*U[6] + x2*U[9] + x3*U[0] + x4*U[3];
    double x31 = x25*x30;
    double x32 = 1.0*x29 + 1.0*x31 + 1.0;
    double x33 = S11*x27 + S12*x32;
    double x34 = x29 + x31;
    double x35 = S12*x27 + S22*x32;
    double x36 = x23 + x26;
    double x37 = S13*x27 + S23*x32;
    double x38 = x10*U[11] + x7*U[2] + x8*U[8] + x9*U[5];
    double x39 = x1*U[8] + x2*U[11] + x3*U[2] + x4*U[5];
    double x40 = x22*x38 + x25*x39;
    double x41 = -x33*x34 - x35*x36 - x37*x40;
    double x42 = x13*x18 + x16*x21;
    double x43 = x15*x5;
    double x44 = x20*x43;
    double x45 = x15*x17;
    double x46 = x24*x45;
    double x47 = 1.0*x44 + 1.0*x46 + 1.0;
    double x48 = x28*x43;
    double x49 = x30*x45;
    double x50 = -1.0*x48 - 1.0*x49;
    double x51 = S11*x47 + S12*x50;
    double x52 = S12*x47 + S22*x50;
    double x53 = S13*x47 + S23*x50;
    double x54 = -1.0*PENER - 1.0*SENER;
    double x55 = -x34*x51 - x36*x52 - x40*x53 - x54;
    double x56 = x48 + x49;
    double x57 = x44 + x46;
    double x58 = x38*x43 + x39*x45;
    double x59 = -x33*x56 - x35*x57 - x37*x58 - x54;
    double x60 = -x51*x56 - x52*x57 - x53*x58;
    double x61 = x4*x45 + x43*x9;
    double x62 = x22*x9 + x25*x4;
    double x63 = x1*x45 + x43*x8;
    double x64 = x1*x25 + x22*x8;
    double x65 = x10*x22 + x2*x25;
    double x66 = x10*x43 + x2*x45;
    
    res_0[0] = x14*(x19*x41 + x42*x55);
    res_0[1] = x14*(x19*x59 + x42*x60);
    res_0[2] = 0;
    res_0[3] = x14*(x41*x61 + x55*x62);
    res_0[4] = x14*(x59*x61 + x60*x62);
    res_0[5] = 0;
    res_0[6] = x14*(x41*x63 + x55*x64);
    res_0[7] = x14*(x59*x63 + x60*x64);
    res_0[8] = 0;
    res_0[9] = x14*(x41*x66 + x55*x65);
    res_0[10] = x14*(x59*x66 + x60*x65);
    res_0[11] = 0;
}

Conf_Forces_API void Integration_CPE4R_static_dbf(size_t num_elem,double Coords[][4][3],double Element_U[][4][3],double S[][1][6],
    double PENER[][1],double SENER[][1],double Conf_Force[][4][3]){
    
    //Inputs:
    //    Coords.........Element_Nodal_Coordinates [num_elem][numNodes,3]
    //    Element_U......Element_Nodal_Displacements [num_elem][numNodes,3]
    //    S..............Stress at integration points [num_elem][num_int_points,6]
    //    PENER..........Plastic strain energy [num_elem][num_int_points]
    //    SENER..........Elastic strain energy [num_elem][num_int_points]
    //    
    double weights[1]={4.0,};
    double int_points[1][3]={
        {0.0,0.0,0.0,}};

    #pragma omp parallel for
    for (size_t i=0;i<num_elem;i++){
        double TMP[4][3];
        for (size_t j=0;j<1;j++){
            CPE4R_static_dbf((double*) int_points[j],(double*) Coords[i],(double*) Element_U[i],(double*) S[i][j],PENER[i][j],SENER[i][j],(double*) TMP);
            for (size_t num_node=0;num_node<4;num_node++){
                for (size_t ii=0;ii<3;ii++){
                    Conf_Force[i][num_node][ii]+=TMP[num_node][ii]*weights[j];
                }
            }
        }
    }
}

void CPE8_dynamic(double *rst,double *coord,double rho,double *Element_U,double *Element_V,double *Element_A,double *S,double PENER,double SENER,double *res_0){
    double r = rst[0];
    double s = rst[1];
    double t = rst[2];
    double *U = Element_U;
    double *V = Element_V;
    double *A = Element_A;
    double S11 = S[0];
    double S22 = S[1];
    double S33 = S[2];
    double S12 = S[3];
    double S13 = S[4];
    double S23 = S[5];
    
    double x0 = r*r;
    double x1 = (1.0/2.0)*x0;
    double x2 = x1 - 1.0/2.0;
    double x3 = -x2;
    double x4 = r*s;
    double x5 = -s + x4;
    double x6 = -s - x4;
    double x7 = (1.0/2.0)*s;
    double x8 = -x7;
    double x9 = (1.0/4.0)*r;
    double x10 = (1.0/4.0)*x0;
    double x11 = (1.0/2.0)*x4;
    double x12 = x10 - x11;
    double x13 = -x12 - x8 - x9;
    double x14 = -x9;
    double x15 = x12 + x14 + x7;
    double x16 = x10 + x11;
    double x17 = -x14 - x16 - x8;
    double x18 = x16 + x7 + x9;
    double x19 = x13*coord[3] + x15*coord[9] + x17*coord[0] + x18*coord[6] + x2*coord[12] + x3*coord[18] + x5*coord[21] + x6*coord[15];
    double x20 = s*s;
    double x21 = (1.0/2.0)*x20;
    double x22 = x21 - 1.0/2.0;
    double x23 = -x22;
    double x24 = -r + x4;
    double x25 = -r - x4;
    double x26 = (1.0/2.0)*r;
    double x27 = (1.0/4.0)*s;
    double x28 = -x27;
    double x29 = (1.0/4.0)*x20;
    double x30 = -x29;
    double x31 = x11 + x30;
    double x32 = x26 + x28 + x31;
    double x33 = -x26;
    double x34 = -x27 - x31 - x33;
    double x35 = x11 + x29;
    double x36 = -x28 - x33 - x35;
    double x37 = x26 + x27 + x35;
    double x38 = x22*coord[22] + x23*coord[16] + x24*coord[13] + x25*coord[19] + x32*coord[10] + x34*coord[4] + x36*coord[1] + x37*coord[7];
    double x39 = x13*coord[4] + x15*coord[10] + x17*coord[1] + x18*coord[7] + x2*coord[13] + x3*coord[19] + x5*coord[22] + x6*coord[16];
    double x40 = x22*coord[21] + x23*coord[15] + x24*coord[12] + x25*coord[18] + x32*coord[9] + x34*coord[3] + x36*coord[0] + x37*coord[6];
    double x41 = -x19*x38 + x39*x40;
    double x42 = 1.0/x41;
    double x43 = x39*x42;
    double x44 = -x38*x42;
    double x45 = x17*x44 + x36*x43;
    double x46 = x13*U[4] + x15*U[10] + x17*U[1] + x18*U[7] + x2*U[13] + x3*U[19] + x5*U[22] + x6*U[16];
    double x47 = x40*x42;
    double x48 = x46*x47;
    double x49 = x22*U[22] + x23*U[16] + x24*U[13] + x25*U[19] + x32*U[10] + x34*U[4] + x36*U[1] + x37*U[7];
    double x50 = -x19*x42;
    double x51 = x49*x50;
    double x52 = 1.0*x48 + 1.0*x51 + 1.0;
    double x53 = x13*U[3] + x15*U[9] + x17*U[0] + x18*U[6] + x2*U[12] + x3*U[18] + x5*U[21] + x6*U[15];
    double x54 = x47*x53;
    double x55 = x22*U[21] + x23*U[15] + x24*U[12] + x25*U[18] + x32*U[9] + x34*U[3] + x36*U[0] + x37*U[6];
    double x56 = x50*x55;
    double x57 = -1.0*x54 - 1.0*x56;
    double x58 = S11*x52 + S12*x57;
    double x59 = x44*x53;
    double x60 = x43*x55;
    double x61 = x59 + x60;
    double x62 = S12*x52 + S22*x57;
    double x63 = x44*x46;
    double x64 = x43*x49;
    double x65 = x63 + x64;
    double x66 = S13*x52 + S23*x57;
    double x67 = x13*U[5] + x15*U[11] + x17*U[2] + x18*U[8] + x2*U[14] + x3*U[20] + x5*U[23] + x6*U[17];
    double x68 = x22*U[23] + x23*U[17] + x24*U[14] + x25*U[20] + x32*U[11] + x34*U[5] + x36*U[2] + x37*U[8];
    double x69 = x43*x68 + x44*x67;
    double x70 = r*x21 + x33;
    double x71 = -x22 - x70;
    double x72 = s*x1 + x8;
    double x73 = -x2 - x72;
    double x74 = -x21 + x70 + 1.0/2.0;
    double x75 = -x1 + x72 + 1.0/2.0;
    double x76 = (1.0/4.0)*x4;
    double x77 = -x76;
    double x78 = x20*x9;
    double x79 = s*x10;
    double x80 = -x78 + x79;
    double x81 = x10 + x29 - 1.0/4.0;
    double x82 = x77 + x80 + x81;
    double x83 = -x10 + x30 + 1.0/4.0;
    double x84 = -x76 - x80 - x83;
    double x85 = x78 + x79;
    double x86 = -x77 - x83 - x85;
    double x87 = x76 + x81 + x85;
    double x88 = x71*V[15] + x73*V[18] + x74*V[21] + x75*V[12] + x82*V[9] + x84*V[3] + x86*V[0] + x87*V[6];
    double x89 = x71*V[16] + x73*V[19] + x74*V[22] + x75*V[13] + x82*V[10] + x84*V[4] + x86*V[1] + x87*V[7];
    double x90 = x71*V[17] + x73*V[20] + x74*V[23] + x75*V[14] + x82*V[11] + x84*V[5] + x86*V[2] + x87*V[8];
    double x91 = -1.0*PENER - 1.0*SENER + 0.5*rho*(x88*x88 + x89*x89 + x90*x90);
    double x92 = -x58*x61 - x62*x65 - x66*x69 - x91;
    double x93 = x17*x47 + x36*x50;
    double x94 = -1.0*x63 - 1.0*x64;
    double x95 = 1.0*x59 + 1.0*x60 + 1.0;
    double x96 = S11*x94 + S12*x95;
    double x97 = S12*x94 + S22*x95;
    double x98 = S13*x94 + S23*x95;
    double x99 = -x61*x96 - x65*x97 - x69*x98;
    double x100 = x71*A[15] + x73*A[18] + x74*A[21] + x75*A[12] + x82*A[9] + x84*A[3] + x86*A[0] + x87*A[6];
    double x101 = x71*A[16] + x73*A[19] + x74*A[22] + x75*A[13] + x82*A[10] + x84*A[4] + x86*A[1] + x87*A[7];
    double x102 = x71*A[17] + x73*A[20] + x74*A[23] + x75*A[14] + x82*A[11] + x84*A[5] + x86*A[2] + x87*A[8];
    double x103 = x13*V[3] + x15*V[9] + x17*V[0] + x18*V[6] + x2*V[12] + x3*V[18] + x5*V[21] + x6*V[15];
    double x104 = x22*V[21] + x23*V[15] + x24*V[12] + x25*V[18] + x32*V[9] + x34*V[3] + x36*V[0] + x37*V[6];
    double x105 = x13*V[4] + x15*V[10] + x17*V[1] + x18*V[7] + x2*V[13] + x3*V[19] + x5*V[22] + x6*V[16];
    double x106 = x22*V[22] + x23*V[16] + x24*V[13] + x25*V[19] + x32*V[10] + x34*V[4] + x36*V[1] + x37*V[7];
    double x107 = x13*V[5] + x15*V[11] + x17*V[2] + x18*V[8] + x2*V[14] + x3*V[20] + x5*V[23] + x6*V[17];
    double x108 = x22*V[23] + x23*V[17] + x24*V[14] + x25*V[20] + x32*V[11] + x34*V[5] + x36*V[2] + x37*V[8];
    double x109 = x100*x61 + x101*x65 + x102*x69 + x88*(x103*x44 + x104*x43) + x89*(x105*x44 + x106*x43) + x90*(x107*x44 + x108*x43);
    double x110 = rho*x86;
    double x111 = x54 + x56;
    double x112 = x48 + x51;
    double x113 = x47*x67 + x50*x68;
    double x114 = -x111*x58 - x112*x62 - x113*x66;
    double x115 = -x111*x96 - x112*x97 - x113*x98 - x91;
    double x116 = x100*x111 + x101*x112 + x102*x113 + x88*(x103*x47 + x104*x50) + x89*(x105*x47 + x106*x50) + x90*(x107*x47 + x108*x50);
    double x117 = x13*x44 + x34*x43;
    double x118 = x13*x47 + x34*x50;
    double x119 = rho*x84;
    double x120 = x18*x44 + x37*x43;
    double x121 = x18*x47 + x37*x50;
    double x122 = rho*x87;
    double x123 = x15*x44 + x32*x43;
    double x124 = x15*x47 + x32*x50;
    double x125 = rho*x82;
    double x126 = x2*x47 + x24*x50;
    double x127 = x2*x44 + x24*x43;
    double x128 = rho*x75;
    double x129 = x23*x43 + x44*x6;
    double x130 = x23*x50 + x47*x6;
    double x131 = rho*x71;
    double x132 = x25*x50 + x3*x47;
    double x133 = x25*x43 + x3*x44;
    double x134 = rho*x73;
    double x135 = x22*x43 + x44*x5;
    double x136 = x22*x50 + x47*x5;
    double x137 = rho*x74;
    
    res_0[0] = x41*(-x109*x110 + x45*x92 + x93*x99);
    res_0[1] = x41*(-x110*x116 + x114*x45 + x115*x93);
    res_0[2] = 0;
    res_0[3] = x41*(-x109*x119 + x117*x92 + x118*x99);
    res_0[4] = x41*(x114*x117 + x115*x118 - x116*x119);
    res_0[5] = 0;
    res_0[6] = x41*(-x109*x122 + x120*x92 + x121*x99);
    res_0[7] = x41*(x114*x120 + x115*x121 - x116*x122);
    res_0[8] = 0;
    res_0[9] = x41*(-x109*x125 + x123*x92 + x124*x99);
    res_0[10] = x41*(x114*x123 + x115*x124 - x116*x125);
    res_0[11] = 0;
    res_0[12] = x41*(-x109*x128 + x126*x99 + x127*x92);
    res_0[13] = x41*(x114*x127 + x115*x126 - x116*x128);
    res_0[14] = 0;
    res_0[15] = x41*(-x109*x131 + x129*x92 + x130*x99);
    res_0[16] = x41*(x114*x129 + x115*x130 - x116*x131);
    res_0[17] = 0;
    res_0[18] = x41*(-x109*x134 + x132*x99 + x133*x92);
    res_0[19] = x41*(x114*x133 + x115*x132 - x116*x134);
    res_0[20] = 0;
    res_0[21] = x41*(-x109*x137 + x135*x92 + x136*x99);
    res_0[22] = x41*(x114*x135 + x115*x136 - x116*x137);
    res_0[23] = 0;
}

Conf_Forces_API void Integration_CPE8_dynamic(size_t num_elem,double Coords[][8][3],double *rho,double Element_U[][8][3],double Element_V[][8][3],double Element_A[][8][3],double S[][9][6],
    double PENER[][9],double SENER[][9],double Conf_Force[][8][3]){
    
    //Inputs:
    //    Coords.........Element_Nodal_Coordinates  [num_elem][numNodes,3]
    //    Element_U......Element_Nodal_Displacements[num_elem][numNodes,3]
    //    Element_V......Element_Nodal_Velocity     [num_elem][numNodes,3]
    //    Element_A......Element_Nodal_Acceleration [num_elem][numNodes,3]
    //    S..............Stress at integration points [num_elem][num_int_points,6]
    //    PENER..........Plastic strain energy [num_elem][num_int_points]
    //    SENER..........Elastic strain energy [num_elem][num_int_points]
    //    
    double weights[9]={0.308641975308642,0.49382716049382713,0.308641975308642,0.49382716049382713,0.7901234567901234,0.49382716049382713,0.308641975308642,0.49382716049382713,0.308641975308642,};
    double int_points[9][3]={
        {-0.7745966692414834,-0.7745966692414834,0.0,},
        {0.0,-0.7745966692414834,0.0,},
        {0.7745966692414834,-0.7745966692414834,0.0,},
        {-0.7745966692414834,0.0,0.0,},
        {0.0,0.0,0.0,},
        {0.7745966692414834,0.0,0.0,},
        {-0.7745966692414834,0.7745966692414834,0.0,},
        {0.0,0.7745966692414834,0.0,},
        {0.7745966692414834,0.7745966692414834,0.0,}};

    #pragma omp parallel for
    for (size_t i=0;i<num_elem;i++){
        double TMP[8][3];
        for (size_t j=0;j<9;j++){
            CPE8_dynamic((double*) int_points[j],(double*) Coords[i],rho[i],(double*) Element_U[i],(double*) Element_V[i],(double*) Element_A[i],(double*) S[i][j],PENER[i][j],SENER[i][j],(double*) TMP);
            for (size_t num_node=0;num_node<8;num_node++){
                for (size_t ii=0;ii<3;ii++){
                    Conf_Force[i][num_node][ii]+=TMP[num_node][ii]*weights[j];
                }
            }
        }
    }
}

void CPE8_static_mbf(double *rst,double *coord,double *Element_U,double *S,double PENER,double SENER,double *res_0){
    double r = rst[0];
    double s = rst[1];
    double t = rst[2];
    double *U = Element_U;
    double S11 = S[0];
    double S22 = S[1];
    double S33 = S[2];
    double S12 = S[3];
    double S13 = S[4];
    double S23 = S[5];
    
    double x0 = (1.0/4.0)*s;
    double x1 = -x0;
    double x2 = -1.0/2.0*r;
    double x3 = r*s;
    double x4 = (1.0/2.0)*x3;
    double x5 = s*s;
    double x6 = (1.0/4.0)*x5;
    double x7 = x4 + x6;
    double x8 = -x1 - x2 - x7;
    double x9 = r*r;
    double x10 = (1.0/2.0)*x9 - 1.0/2.0;
    double x11 = -x10;
    double x12 = -s + x3;
    double x13 = -s - x3;
    double x14 = -1.0/2.0*s;
    double x15 = (1.0/4.0)*r;
    double x16 = (1.0/4.0)*x9;
    double x17 = x16 - x4;
    double x18 = -x14 - x15 - x17;
    double x19 = (1.0/2.0)*s;
    double x20 = -x15;
    double x21 = x17 + x19 + x20;
    double x22 = x16 + x4;
    double x23 = -x14 - x20 - x22;
    double x24 = x15 + x19 + x22;
    double x25 = x10*coord[13] + x11*coord[19] + x12*coord[22] + x13*coord[16] + x18*coord[4] + x21*coord[10] + x23*coord[1] + x24*coord[7];
    double x26 = x10*coord[12] + x11*coord[18] + x12*coord[21] + x13*coord[15] + x18*coord[3] + x21*coord[9] + x23*coord[0] + x24*coord[6];
    double x27 = (1.0/2.0)*x5 - 1.0/2.0;
    double x28 = -x27;
    double x29 = -r + x3;
    double x30 = -r - x3;
    double x31 = (1.0/2.0)*r;
    double x32 = x4 - x6;
    double x33 = x1 + x31 + x32;
    double x34 = -x0 - x2 - x32;
    double x35 = x0 + x31 + x7;
    double x36 = x27*coord[22] + x28*coord[16] + x29*coord[13] + x30*coord[19] + x33*coord[10] + x34*coord[4] + x35*coord[7] + x8*coord[1];
    double x37 = x27*coord[21] + x28*coord[15] + x29*coord[12] + x30*coord[18] + x33*coord[9] + x34*coord[3] + x35*coord[6] + x8*coord[0];
    double x38 = x25*x37 - x26*x36;
    double x39 = 1.0/x38;
    double x40 = x25*x39;
    double x41 = -x36*x39;
    double x42 = x23*x41 + x40*x8;
    double x43 = x10*U[13] + x11*U[19] + x12*U[22] + x13*U[16] + x18*U[4] + x21*U[10] + x23*U[1] + x24*U[7];
    double x44 = x37*x39;
    double x45 = x43*x44;
    double x46 = x27*U[22] + x28*U[16] + x29*U[13] + x30*U[19] + x33*U[10] + x34*U[4] + x35*U[7] + x8*U[1];
    double x47 = -x26*x39;
    double x48 = x46*x47;
    double x49 = 1.0*x45 + 1.0*x48 + 1.0;
    double x50 = x10*U[12] + x11*U[18] + x12*U[21] + x13*U[15] + x18*U[3] + x21*U[9] + x23*U[0] + x24*U[6];
    double x51 = x44*x50;
    double x52 = x27*U[21] + x28*U[15] + x29*U[12] + x30*U[18] + x33*U[9] + x34*U[3] + x35*U[6] + x8*U[0];
    double x53 = x47*x52;
    double x54 = -1.0*x51 - 1.0*x53;
    double x55 = S11*x49 + S12*x54;
    double x56 = x41*x50;
    double x57 = x40*x52;
    double x58 = x56 + x57 + 1.0;
    double x59 = S12*x49 + S22*x54;
    double x60 = x41*x43;
    double x61 = x40*x46;
    double x62 = x60 + x61;
    double x63 = S13*x49;
    double x64 = S23*x54;
    double x65 = x63 + x64;
    double x66 = x10*U[14] + x11*U[20] + x12*U[23] + x13*U[17] + x18*U[5] + x21*U[11] + x23*U[2] + x24*U[8];
    double x67 = x27*U[23] + x28*U[17] + x29*U[14] + x30*U[20] + x33*U[11] + x34*U[5] + x35*U[8] + x8*U[2];
    double x68 = x40*x67 + x41*x66;
    double x69 = -1.0*PENER - 1.0*SENER;
    double x70 = -x55*x58 - x59*x62 - x65*x68 - x69;
    double x71 = x23*x44 + x47*x8;
    double x72 = -1.0*x60 - 1.0*x61;
    double x73 = 1.0*x56 + 1.0*x57 + 1.0;
    double x74 = S11*x72 + S12*x73;
    double x75 = S12*x72 + S22*x73;
    double x76 = S13*x72;
    double x77 = S23*x73;
    double x78 = x76 + x77;
    double x79 = -x58*x74 - x62*x75 - x68*x78;
    double x80 = x51 + x53;
    double x81 = x45 + x48 + 1.0;
    double x82 = x44*x66 + x47*x67;
    double x83 = -x55*x80 - x59*x81 - x65*x82;
    double x84 = -x69 - x74*x80 - x75*x81 - x78*x82;
    double x85 = -1.0*x76 - 1.0*x77;
    double x86 = -1.0*x63 - 1.0*x64;
    double x87 = x18*x41 + x34*x40;
    double x88 = x18*x44 + x34*x47;
    double x89 = x24*x41 + x35*x40;
    double x90 = x24*x44 + x35*x47;
    double x91 = x21*x41 + x33*x40;
    double x92 = x21*x44 + x33*x47;
    double x93 = x10*x44 + x29*x47;
    double x94 = x10*x41 + x29*x40;
    double x95 = x13*x41 + x28*x40;
    double x96 = x13*x44 + x28*x47;
    double x97 = x11*x44 + x30*x47;
    double x98 = x11*x41 + x30*x40;
    double x99 = x12*x41 + x27*x40;
    double x100 = x12*x44 + x27*x47;
    
    res_0[0] = x38*(x42*x70 + x71*x79);
    res_0[1] = x38*(x42*x83 + x71*x84);
    res_0[2] = x38*(x42*x86 + x71*x85);
    res_0[3] = x38*(x70*x87 + x79*x88);
    res_0[4] = x38*(x83*x87 + x84*x88);
    res_0[5] = x38*(x85*x88 + x86*x87);
    res_0[6] = x38*(x70*x89 + x79*x90);
    res_0[7] = x38*(x83*x89 + x84*x90);
    res_0[8] = x38*(x85*x90 + x86*x89);
    res_0[9] = x38*(x70*x91 + x79*x92);
    res_0[10] = x38*(x83*x91 + x84*x92);
    res_0[11] = x38*(x85*x92 + x86*x91);
    res_0[12] = x38*(x70*x94 + x79*x93);
    res_0[13] = x38*(x83*x94 + x84*x93);
    res_0[14] = x38*(x85*x93 + x86*x94);
    res_0[15] = x38*(x70*x95 + x79*x96);
    res_0[16] = x38*(x83*x95 + x84*x96);
    res_0[17] = x38*(x85*x96 + x86*x95);
    res_0[18] = x38*(x70*x98 + x79*x97);
    res_0[19] = x38*(x83*x98 + x84*x97);
    res_0[20] = x38*(x85*x97 + x86*x98);
    res_0[21] = x38*(x100*x79 + x70*x99);
    res_0[22] = x38*(x100*x84 + x83*x99);
    res_0[23] = x38*(x100*x85 + x86*x99);
}

Conf_Forces_API void Integration_CPE8_static_mbf(size_t num_elem,double Coords[][8][3],double Element_U[][8][3],double S[][9][6],
    double PENER[][9],double SENER[][9],double Conf_Force[][8][3]){
    
    //Inputs:
    //    Coords.........Element_Nodal_Coordinates [num_elem][numNodes,3]
    //    Element_U......Element_Nodal_Displacements [num_elem][numNodes,3]
    //    S..............Stress at integration points [num_elem][num_int_points,6]
    //    PENER..........Plastic strain energy [num_elem][num_int_points]
    //    SENER..........Elastic strain energy [num_elem][num_int_points]
    //    
    double weights[9]={0.308641975308642,0.49382716049382713,0.308641975308642,0.49382716049382713,0.7901234567901234,0.49382716049382713,0.308641975308642,0.49382716049382713,0.308641975308642,};
    double int_points[9][3]={
        {-0.7745966692414834,-0.7745966692414834,0.0,},
        {0.0,-0.7745966692414834,0.0,},
        {0.7745966692414834,-0.7745966692414834,0.0,},
        {-0.7745966692414834,0.0,0.0,},
        {0.0,0.0,0.0,},
        {0.7745966692414834,0.0,0.0,},
        {-0.7745966692414834,0.7745966692414834,0.0,},
        {0.0,0.7745966692414834,0.0,},
        {0.7745966692414834,0.7745966692414834,0.0,}};

    #pragma omp parallel for
    for (size_t i=0;i<num_elem;i++){
        double TMP[8][3];
        for (size_t j=0;j<9;j++){
            CPE8_static_mbf((double*) int_points[j],(double*) Coords[i],(double*) Element_U[i],(double*) S[i][j],PENER[i][j],SENER[i][j],(double*) TMP);
            for (size_t num_node=0;num_node<8;num_node++){
                for (size_t ii=0;ii<3;ii++){
                    Conf_Force[i][num_node][ii]+=TMP[num_node][ii]*weights[j];
                }
            }
        }
    }
}

void CPE8_static_dbf(double *rst,double *coord,double *Element_U,double *S,double PENER,double SENER,double *res_0){
    double r = rst[0];
    double s = rst[1];
    double t = rst[2];
    double *U = Element_U;
    double S11 = S[0];
    double S22 = S[1];
    double S33 = S[2];
    double S12 = S[3];
    double S13 = S[4];
    double S23 = S[5];
    
    double x0 = (1.0/4.0)*s;
    double x1 = -x0;
    double x2 = -1.0/2.0*r;
    double x3 = r*s;
    double x4 = (1.0/2.0)*x3;
    double x5 = s*s;
    double x6 = (1.0/4.0)*x5;
    double x7 = x4 + x6;
    double x8 = -x1 - x2 - x7;
    double x9 = r*r;
    double x10 = (1.0/2.0)*x9 - 1.0/2.0;
    double x11 = -x10;
    double x12 = -s + x3;
    double x13 = -s - x3;
    double x14 = -1.0/2.0*s;
    double x15 = (1.0/4.0)*r;
    double x16 = (1.0/4.0)*x9;
    double x17 = x16 - x4;
    double x18 = -x14 - x15 - x17;
    double x19 = (1.0/2.0)*s;
    double x20 = -x15;
    double x21 = x17 + x19 + x20;
    double x22 = x16 + x4;
    double x23 = -x14 - x20 - x22;
    double x24 = x15 + x19 + x22;
    double x25 = x10*coord[13] + x11*coord[19] + x12*coord[22] + x13*coord[16] + x18*coord[4] + x21*coord[10] + x23*coord[1] + x24*coord[7];
    double x26 = x10*coord[12] + x11*coord[18] + x12*coord[21] + x13*coord[15] + x18*coord[3] + x21*coord[9] + x23*coord[0] + x24*coord[6];
    double x27 = (1.0/2.0)*x5 - 1.0/2.0;
    double x28 = -x27;
    double x29 = -r + x3;
    double x30 = -r - x3;
    double x31 = (1.0/2.0)*r;
    double x32 = x4 - x6;
    double x33 = x1 + x31 + x32;
    double x34 = -x0 - x2 - x32;
    double x35 = x0 + x31 + x7;
    double x36 = x27*coord[22] + x28*coord[16] + x29*coord[13] + x30*coord[19] + x33*coord[10] + x34*coord[4] + x35*coord[7] + x8*coord[1];
    double x37 = x27*coord[21] + x28*coord[15] + x29*coord[12] + x30*coord[18] + x33*coord[9] + x34*coord[3] + x35*coord[6] + x8*coord[0];
    double x38 = x25*x37 - x26*x36;
    double x39 = 1.0/x38;
    double x40 = x25*x39;
    double x41 = -x36*x39;
    double x42 = x23*x41 + x40*x8;
    double x43 = x10*U[13] + x11*U[19] + x12*U[22] + x13*U[16] + x18*U[4] + x21*U[10] + x23*U[1] + x24*U[7];
    double x44 = x37*x39;
    double x45 = x43*x44;
    double x46 = x27*U[22] + x28*U[16] + x29*U[13] + x30*U[19] + x33*U[10] + x34*U[4] + x35*U[7] + x8*U[1];
    double x47 = -x26*x39;
    double x48 = x46*x47;
    double x49 = 1.0*x45 + 1.0*x48 + 1.0;
    double x50 = x10*U[12] + x11*U[18] + x12*U[21] + x13*U[15] + x18*U[3] + x21*U[9] + x23*U[0] + x24*U[6];
    double x51 = x44*x50;
    double x52 = x27*U[21] + x28*U[15] + x29*U[12] + x30*U[18] + x33*U[9] + x34*U[3] + x35*U[6] + x8*U[0];
    double x53 = x47*x52;
    double x54 = -1.0*x51 - 1.0*x53;
    double x55 = S11*x49 + S12*x54;
    double x56 = x41*x50;
    double x57 = x40*x52;
    double x58 = x56 + x57;
    double x59 = S12*x49 + S22*x54;
    double x60 = x41*x43;
    double x61 = x40*x46;
    double x62 = x60 + x61;
    double x63 = S13*x49 + S23*x54;
    double x64 = x10*U[14] + x11*U[20] + x12*U[23] + x13*U[17] + x18*U[5] + x21*U[11] + x23*U[2] + x24*U[8];
    double x65 = x27*U[23] + x28*U[17] + x29*U[14] + x30*U[20] + x33*U[11] + x34*U[5] + x35*U[8] + x8*U[2];
    double x66 = x40*x65 + x41*x64;
    double x67 = -1.0*PENER - 1.0*SENER;
    double x68 = -x55*x58 - x59*x62 - x63*x66 - x67;
    double x69 = x23*x44 + x47*x8;
    double x70 = -1.0*x60 - 1.0*x61;
    double x71 = 1.0*x56 + 1.0*x57 + 1.0;
    double x72 = S11*x70 + S12*x71;
    double x73 = S12*x70 + S22*x71;
    double x74 = S13*x70 + S23*x71;
    double x75 = -x58*x72 - x62*x73 - x66*x74;
    double x76 = x51 + x53;
    double x77 = x45 + x48;
    double x78 = x44*x64 + x47*x65;
    double x79 = -x55*x76 - x59*x77 - x63*x78;
    double x80 = -x67 - x72*x76 - x73*x77 - x74*x78;
    double x81 = x18*x41 + x34*x40;
    double x82 = x18*x44 + x34*x47;
    double x83 = x24*x41 + x35*x40;
    double x84 = x24*x44 + x35*x47;
    double x85 = x21*x41 + x33*x40;
    double x86 = x21*x44 + x33*x47;
    double x87 = x10*x44 + x29*x47;
    double x88 = x10*x41 + x29*x40;
    double x89 = x13*x41 + x28*x40;
    double x90 = x13*x44 + x28*x47;
    double x91 = x11*x44 + x30*x47;
    double x92 = x11*x41 + x30*x40;
    double x93 = x12*x41 + x27*x40;
    double x94 = x12*x44 + x27*x47;
    
    res_0[0] = x38*(x42*x68 + x69*x75);
    res_0[1] = x38*(x42*x79 + x69*x80);
    res_0[2] = 0;
    res_0[3] = x38*(x68*x81 + x75*x82);
    res_0[4] = x38*(x79*x81 + x80*x82);
    res_0[5] = 0;
    res_0[6] = x38*(x68*x83 + x75*x84);
    res_0[7] = x38*(x79*x83 + x80*x84);
    res_0[8] = 0;
    res_0[9] = x38*(x68*x85 + x75*x86);
    res_0[10] = x38*(x79*x85 + x80*x86);
    res_0[11] = 0;
    res_0[12] = x38*(x68*x88 + x75*x87);
    res_0[13] = x38*(x79*x88 + x80*x87);
    res_0[14] = 0;
    res_0[15] = x38*(x68*x89 + x75*x90);
    res_0[16] = x38*(x79*x89 + x80*x90);
    res_0[17] = 0;
    res_0[18] = x38*(x68*x92 + x75*x91);
    res_0[19] = x38*(x79*x92 + x80*x91);
    res_0[20] = 0;
    res_0[21] = x38*(x68*x93 + x75*x94);
    res_0[22] = x38*(x79*x93 + x80*x94);
    res_0[23] = 0;
}

Conf_Forces_API void Integration_CPE8_static_dbf(size_t num_elem,double Coords[][8][3],double Element_U[][8][3],double S[][9][6],
    double PENER[][9],double SENER[][9],double Conf_Force[][8][3]){
    
    //Inputs:
    //    Coords.........Element_Nodal_Coordinates [num_elem][numNodes,3]
    //    Element_U......Element_Nodal_Displacements [num_elem][numNodes,3]
    //    S..............Stress at integration points [num_elem][num_int_points,6]
    //    PENER..........Plastic strain energy [num_elem][num_int_points]
    //    SENER..........Elastic strain energy [num_elem][num_int_points]
    //    
    double weights[9]={0.308641975308642,0.49382716049382713,0.308641975308642,0.49382716049382713,0.7901234567901234,0.49382716049382713,0.308641975308642,0.49382716049382713,0.308641975308642,};
    double int_points[9][3]={
        {-0.7745966692414834,-0.7745966692414834,0.0,},
        {0.0,-0.7745966692414834,0.0,},
        {0.7745966692414834,-0.7745966692414834,0.0,},
        {-0.7745966692414834,0.0,0.0,},
        {0.0,0.0,0.0,},
        {0.7745966692414834,0.0,0.0,},
        {-0.7745966692414834,0.7745966692414834,0.0,},
        {0.0,0.7745966692414834,0.0,},
        {0.7745966692414834,0.7745966692414834,0.0,}};

    #pragma omp parallel for
    for (size_t i=0;i<num_elem;i++){
        double TMP[8][3];
        for (size_t j=0;j<9;j++){
            CPE8_static_dbf((double*) int_points[j],(double*) Coords[i],(double*) Element_U[i],(double*) S[i][j],PENER[i][j],SENER[i][j],(double*) TMP);
            for (size_t num_node=0;num_node<8;num_node++){
                for (size_t ii=0;ii<3;ii++){
                    Conf_Force[i][num_node][ii]+=TMP[num_node][ii]*weights[j];
                }
            }
        }
    }
}

void CPE8R_dynamic(double *rst,double *coord,double rho,double *Element_U,double *Element_V,double *Element_A,double *S,double PENER,double SENER,double *res_0){
    double r = rst[0];
    double s = rst[1];
    double t = rst[2];
    double *U = Element_U;
    double *V = Element_V;
    double *A = Element_A;
    double S11 = S[0];
    double S22 = S[1];
    double S33 = S[2];
    double S12 = S[3];
    double S13 = S[4];
    double S23 = S[5];
    
    double x0 = r*r;
    double x1 = (1.0/2.0)*x0;
    double x2 = x1 - 1.0/2.0;
    double x3 = -x2;
    double x4 = r*s;
    double x5 = -s + x4;
    double x6 = -s - x4;
    double x7 = (1.0/2.0)*s;
    double x8 = -x7;
    double x9 = (1.0/4.0)*r;
    double x10 = (1.0/4.0)*x0;
    double x11 = (1.0/2.0)*x4;
    double x12 = x10 - x11;
    double x13 = -x12 - x8 - x9;
    double x14 = -x9;
    double x15 = x12 + x14 + x7;
    double x16 = x10 + x11;
    double x17 = -x14 - x16 - x8;
    double x18 = x16 + x7 + x9;
    double x19 = x13*coord[3] + x15*coord[9] + x17*coord[0] + x18*coord[6] + x2*coord[12] + x3*coord[18] + x5*coord[21] + x6*coord[15];
    double x20 = s*s;
    double x21 = (1.0/2.0)*x20;
    double x22 = x21 - 1.0/2.0;
    double x23 = -x22;
    double x24 = -r + x4;
    double x25 = -r - x4;
    double x26 = (1.0/2.0)*r;
    double x27 = (1.0/4.0)*s;
    double x28 = -x27;
    double x29 = (1.0/4.0)*x20;
    double x30 = -x29;
    double x31 = x11 + x30;
    double x32 = x26 + x28 + x31;
    double x33 = -x26;
    double x34 = -x27 - x31 - x33;
    double x35 = x11 + x29;
    double x36 = -x28 - x33 - x35;
    double x37 = x26 + x27 + x35;
    double x38 = x22*coord[22] + x23*coord[16] + x24*coord[13] + x25*coord[19] + x32*coord[10] + x34*coord[4] + x36*coord[1] + x37*coord[7];
    double x39 = x13*coord[4] + x15*coord[10] + x17*coord[1] + x18*coord[7] + x2*coord[13] + x3*coord[19] + x5*coord[22] + x6*coord[16];
    double x40 = x22*coord[21] + x23*coord[15] + x24*coord[12] + x25*coord[18] + x32*coord[9] + x34*coord[3] + x36*coord[0] + x37*coord[6];
    double x41 = -x19*x38 + x39*x40;
    double x42 = 1.0/x41;
    double x43 = x39*x42;
    double x44 = -x38*x42;
    double x45 = x17*x44 + x36*x43;
    double x46 = x13*U[4] + x15*U[10] + x17*U[1] + x18*U[7] + x2*U[13] + x3*U[19] + x5*U[22] + x6*U[16];
    double x47 = x40*x42;
    double x48 = x46*x47;
    double x49 = x22*U[22] + x23*U[16] + x24*U[13] + x25*U[19] + x32*U[10] + x34*U[4] + x36*U[1] + x37*U[7];
    double x50 = -x19*x42;
    double x51 = x49*x50;
    double x52 = 1.0*x48 + 1.0*x51 + 1.0;
    double x53 = x13*U[3] + x15*U[9] + x17*U[0] + x18*U[6] + x2*U[12] + x3*U[18] + x5*U[21] + x6*U[15];
    double x54 = x47*x53;
    double x55 = x22*U[21] + x23*U[15] + x24*U[12] + x25*U[18] + x32*U[9] + x34*U[3] + x36*U[0] + x37*U[6];
    double x56 = x50*x55;
    double x57 = -1.0*x54 - 1.0*x56;
    double x58 = S11*x52 + S12*x57;
    double x59 = x44*x53;
    double x60 = x43*x55;
    double x61 = x59 + x60;
    double x62 = S12*x52 + S22*x57;
    double x63 = x44*x46;
    double x64 = x43*x49;
    double x65 = x63 + x64;
    double x66 = S13*x52 + S23*x57;
    double x67 = x13*U[5] + x15*U[11] + x17*U[2] + x18*U[8] + x2*U[14] + x3*U[20] + x5*U[23] + x6*U[17];
    double x68 = x22*U[23] + x23*U[17] + x24*U[14] + x25*U[20] + x32*U[11] + x34*U[5] + x36*U[2] + x37*U[8];
    double x69 = x43*x68 + x44*x67;
    double x70 = r*x21 + x33;
    double x71 = -x22 - x70;
    double x72 = s*x1 + x8;
    double x73 = -x2 - x72;
    double x74 = -x21 + x70 + 1.0/2.0;
    double x75 = -x1 + x72 + 1.0/2.0;
    double x76 = (1.0/4.0)*x4;
    double x77 = -x76;
    double x78 = x20*x9;
    double x79 = s*x10;
    double x80 = -x78 + x79;
    double x81 = x10 + x29 - 1.0/4.0;
    double x82 = x77 + x80 + x81;
    double x83 = -x10 + x30 + 1.0/4.0;
    double x84 = -x76 - x80 - x83;
    double x85 = x78 + x79;
    double x86 = -x77 - x83 - x85;
    double x87 = x76 + x81 + x85;
    double x88 = x71*V[15] + x73*V[18] + x74*V[21] + x75*V[12] + x82*V[9] + x84*V[3] + x86*V[0] + x87*V[6];
    double x89 = x71*V[16] + x73*V[19] + x74*V[22] + x75*V[13] + x82*V[10] + x84*V[4] + x86*V[1] + x87*V[7];
    double x90 = x71*V[17] + x73*V[20] + x74*V[23] + x75*V[14] + x82*V[11] + x84*V[5] + x86*V[2] + x87*V[8];
    double x91 = -1.0*PENER - 1.0*SENER + 0.5*rho*(x88*x88 + x89*x89 + x90*x90);
    double x92 = -x58*x61 - x62*x65 - x66*x69 - x91;
    double x93 = x17*x47 + x36*x50;
    double x94 = -1.0*x63 - 1.0*x64;
    double x95 = 1.0*x59 + 1.0*x60 + 1.0;
    double x96 = S11*x94 + S12*x95;
    double x97 = S12*x94 + S22*x95;
    double x98 = S13*x94 + S23*x95;
    double x99 = -x61*x96 - x65*x97 - x69*x98;
    double x100 = x71*A[15] + x73*A[18] + x74*A[21] + x75*A[12] + x82*A[9] + x84*A[3] + x86*A[0] + x87*A[6];
    double x101 = x71*A[16] + x73*A[19] + x74*A[22] + x75*A[13] + x82*A[10] + x84*A[4] + x86*A[1] + x87*A[7];
    double x102 = x71*A[17] + x73*A[20] + x74*A[23] + x75*A[14] + x82*A[11] + x84*A[5] + x86*A[2] + x87*A[8];
    double x103 = x13*V[3] + x15*V[9] + x17*V[0] + x18*V[6] + x2*V[12] + x3*V[18] + x5*V[21] + x6*V[15];
    double x104 = x22*V[21] + x23*V[15] + x24*V[12] + x25*V[18] + x32*V[9] + x34*V[3] + x36*V[0] + x37*V[6];
    double x105 = x13*V[4] + x15*V[10] + x17*V[1] + x18*V[7] + x2*V[13] + x3*V[19] + x5*V[22] + x6*V[16];
    double x106 = x22*V[22] + x23*V[16] + x24*V[13] + x25*V[19] + x32*V[10] + x34*V[4] + x36*V[1] + x37*V[7];
    double x107 = x13*V[5] + x15*V[11] + x17*V[2] + x18*V[8] + x2*V[14] + x3*V[20] + x5*V[23] + x6*V[17];
    double x108 = x22*V[23] + x23*V[17] + x24*V[14] + x25*V[20] + x32*V[11] + x34*V[5] + x36*V[2] + x37*V[8];
    double x109 = x100*x61 + x101*x65 + x102*x69 + x88*(x103*x44 + x104*x43) + x89*(x105*x44 + x106*x43) + x90*(x107*x44 + x108*x43);
    double x110 = rho*x86;
    double x111 = x54 + x56;
    double x112 = x48 + x51;
    double x113 = x47*x67 + x50*x68;
    double x114 = -x111*x58 - x112*x62 - x113*x66;
    double x115 = -x111*x96 - x112*x97 - x113*x98 - x91;
    double x116 = x100*x111 + x101*x112 + x102*x113 + x88*(x103*x47 + x104*x50) + x89*(x105*x47 + x106*x50) + x90*(x107*x47 + x108*x50);
    double x117 = x13*x44 + x34*x43;
    double x118 = x13*x47 + x34*x50;
    double x119 = rho*x84;
    double x120 = x18*x44 + x37*x43;
    double x121 = x18*x47 + x37*x50;
    double x122 = rho*x87;
    double x123 = x15*x44 + x32*x43;
    double x124 = x15*x47 + x32*x50;
    double x125 = rho*x82;
    double x126 = x2*x47 + x24*x50;
    double x127 = x2*x44 + x24*x43;
    double x128 = rho*x75;
    double x129 = x23*x43 + x44*x6;
    double x130 = x23*x50 + x47*x6;
    double x131 = rho*x71;
    double x132 = x25*x50 + x3*x47;
    double x133 = x25*x43 + x3*x44;
    double x134 = rho*x73;
    double x135 = x22*x43 + x44*x5;
    double x136 = x22*x50 + x47*x5;
    double x137 = rho*x74;
    
    res_0[0] = x41*(-x109*x110 + x45*x92 + x93*x99);
    res_0[1] = x41*(-x110*x116 + x114*x45 + x115*x93);
    res_0[2] = 0;
    res_0[3] = x41*(-x109*x119 + x117*x92 + x118*x99);
    res_0[4] = x41*(x114*x117 + x115*x118 - x116*x119);
    res_0[5] = 0;
    res_0[6] = x41*(-x109*x122 + x120*x92 + x121*x99);
    res_0[7] = x41*(x114*x120 + x115*x121 - x116*x122);
    res_0[8] = 0;
    res_0[9] = x41*(-x109*x125 + x123*x92 + x124*x99);
    res_0[10] = x41*(x114*x123 + x115*x124 - x116*x125);
    res_0[11] = 0;
    res_0[12] = x41*(-x109*x128 + x126*x99 + x127*x92);
    res_0[13] = x41*(x114*x127 + x115*x126 - x116*x128);
    res_0[14] = 0;
    res_0[15] = x41*(-x109*x131 + x129*x92 + x130*x99);
    res_0[16] = x41*(x114*x129 + x115*x130 - x116*x131);
    res_0[17] = 0;
    res_0[18] = x41*(-x109*x134 + x132*x99 + x133*x92);
    res_0[19] = x41*(x114*x133 + x115*x132 - x116*x134);
    res_0[20] = 0;
    res_0[21] = x41*(-x109*x137 + x135*x92 + x136*x99);
    res_0[22] = x41*(x114*x135 + x115*x136 - x116*x137);
    res_0[23] = 0;
}

Conf_Forces_API void Integration_CPE8R_dynamic(size_t num_elem,double Coords[][8][3],double *rho,double Element_U[][8][3],double Element_V[][8][3],double Element_A[][8][3],double S[][4][6],
    double PENER[][4],double SENER[][4],double Conf_Force[][8][3]){
    
    //Inputs:
    //    Coords.........Element_Nodal_Coordinates  [num_elem][numNodes,3]
    //    Element_U......Element_Nodal_Displacements[num_elem][numNodes,3]
    //    Element_V......Element_Nodal_Velocity     [num_elem][numNodes,3]
    //    Element_A......Element_Nodal_Acceleration [num_elem][numNodes,3]
    //    S..............Stress at integration points [num_elem][num_int_points,6]
    //    PENER..........Plastic strain energy [num_elem][num_int_points]
    //    SENER..........Elastic strain energy [num_elem][num_int_points]
    //    
    double weights[4]={1.0,1.0,1.0,1.0,};
    double int_points[4][3]={
        {-0.5773502691896258,-0.5773502691896258,0.0,},
        {0.5773502691896258,-0.5773502691896258,0.0,},
        {-0.5773502691896258,0.5773502691896258,0.0,},
        {0.5773502691896258,0.5773502691896258,0.0,}};

    #pragma omp parallel for
    for (size_t i=0;i<num_elem;i++){
        double TMP[8][3];
        for (size_t j=0;j<4;j++){
            CPE8R_dynamic((double*) int_points[j],(double*) Coords[i],rho[i],(double*) Element_U[i],(double*) Element_V[i],(double*) Element_A[i],(double*) S[i][j],PENER[i][j],SENER[i][j],(double*) TMP);
            for (size_t num_node=0;num_node<8;num_node++){
                for (size_t ii=0;ii<3;ii++){
                    Conf_Force[i][num_node][ii]+=TMP[num_node][ii]*weights[j];
                }
            }
        }
    }
}

void CPE8R_static_mbf(double *rst,double *coord,double *Element_U,double *S,double PENER,double SENER,double *res_0){
    double r = rst[0];
    double s = rst[1];
    double t = rst[2];
    double *U = Element_U;
    double S11 = S[0];
    double S22 = S[1];
    double S33 = S[2];
    double S12 = S[3];
    double S13 = S[4];
    double S23 = S[5];
    
    double x0 = (1.0/4.0)*s;
    double x1 = -x0;
    double x2 = -1.0/2.0*r;
    double x3 = r*s;
    double x4 = (1.0/2.0)*x3;
    double x5 = s*s;
    double x6 = (1.0/4.0)*x5;
    double x7 = x4 + x6;
    double x8 = -x1 - x2 - x7;
    double x9 = r*r;
    double x10 = (1.0/2.0)*x9 - 1.0/2.0;
    double x11 = -x10;
    double x12 = -s + x3;
    double x13 = -s - x3;
    double x14 = -1.0/2.0*s;
    double x15 = (1.0/4.0)*r;
    double x16 = (1.0/4.0)*x9;
    double x17 = x16 - x4;
    double x18 = -x14 - x15 - x17;
    double x19 = (1.0/2.0)*s;
    double x20 = -x15;
    double x21 = x17 + x19 + x20;
    double x22 = x16 + x4;
    double x23 = -x14 - x20 - x22;
    double x24 = x15 + x19 + x22;
    double x25 = x10*coord[13] + x11*coord[19] + x12*coord[22] + x13*coord[16] + x18*coord[4] + x21*coord[10] + x23*coord[1] + x24*coord[7];
    double x26 = x10*coord[12] + x11*coord[18] + x12*coord[21] + x13*coord[15] + x18*coord[3] + x21*coord[9] + x23*coord[0] + x24*coord[6];
    double x27 = (1.0/2.0)*x5 - 1.0/2.0;
    double x28 = -x27;
    double x29 = -r + x3;
    double x30 = -r - x3;
    double x31 = (1.0/2.0)*r;
    double x32 = x4 - x6;
    double x33 = x1 + x31 + x32;
    double x34 = -x0 - x2 - x32;
    double x35 = x0 + x31 + x7;
    double x36 = x27*coord[22] + x28*coord[16] + x29*coord[13] + x30*coord[19] + x33*coord[10] + x34*coord[4] + x35*coord[7] + x8*coord[1];
    double x37 = x27*coord[21] + x28*coord[15] + x29*coord[12] + x30*coord[18] + x33*coord[9] + x34*coord[3] + x35*coord[6] + x8*coord[0];
    double x38 = x25*x37 - x26*x36;
    double x39 = 1.0/x38;
    double x40 = x25*x39;
    double x41 = -x36*x39;
    double x42 = x23*x41 + x40*x8;
    double x43 = x10*U[13] + x11*U[19] + x12*U[22] + x13*U[16] + x18*U[4] + x21*U[10] + x23*U[1] + x24*U[7];
    double x44 = x37*x39;
    double x45 = x43*x44;
    double x46 = x27*U[22] + x28*U[16] + x29*U[13] + x30*U[19] + x33*U[10] + x34*U[4] + x35*U[7] + x8*U[1];
    double x47 = -x26*x39;
    double x48 = x46*x47;
    double x49 = 1.0*x45 + 1.0*x48 + 1.0;
    double x50 = x10*U[12] + x11*U[18] + x12*U[21] + x13*U[15] + x18*U[3] + x21*U[9] + x23*U[0] + x24*U[6];
    double x51 = x44*x50;
    double x52 = x27*U[21] + x28*U[15] + x29*U[12] + x30*U[18] + x33*U[9] + x34*U[3] + x35*U[6] + x8*U[0];
    double x53 = x47*x52;
    double x54 = -1.0*x51 - 1.0*x53;
    double x55 = S11*x49 + S12*x54;
    double x56 = x41*x50;
    double x57 = x40*x52;
    double x58 = x56 + x57 + 1.0;
    double x59 = S12*x49 + S22*x54;
    double x60 = x41*x43;
    double x61 = x40*x46;
    double x62 = x60 + x61;
    double x63 = S13*x49;
    double x64 = S23*x54;
    double x65 = x63 + x64;
    double x66 = x10*U[14] + x11*U[20] + x12*U[23] + x13*U[17] + x18*U[5] + x21*U[11] + x23*U[2] + x24*U[8];
    double x67 = x27*U[23] + x28*U[17] + x29*U[14] + x30*U[20] + x33*U[11] + x34*U[5] + x35*U[8] + x8*U[2];
    double x68 = x40*x67 + x41*x66;
    double x69 = -1.0*PENER - 1.0*SENER;
    double x70 = -x55*x58 - x59*x62 - x65*x68 - x69;
    double x71 = x23*x44 + x47*x8;
    double x72 = -1.0*x60 - 1.0*x61;
    double x73 = 1.0*x56 + 1.0*x57 + 1.0;
    double x74 = S11*x72 + S12*x73;
    double x75 = S12*x72 + S22*x73;
    double x76 = S13*x72;
    double x77 = S23*x73;
    double x78 = x76 + x77;
    double x79 = -x58*x74 - x62*x75 - x68*x78;
    double x80 = x51 + x53;
    double x81 = x45 + x48 + 1.0;
    double x82 = x44*x66 + x47*x67;
    double x83 = -x55*x80 - x59*x81 - x65*x82;
    double x84 = -x69 - x74*x80 - x75*x81 - x78*x82;
    double x85 = -1.0*x76 - 1.0*x77;
    double x86 = -1.0*x63 - 1.0*x64;
    double x87 = x18*x41 + x34*x40;
    double x88 = x18*x44 + x34*x47;
    double x89 = x24*x41 + x35*x40;
    double x90 = x24*x44 + x35*x47;
    double x91 = x21*x41 + x33*x40;
    double x92 = x21*x44 + x33*x47;
    double x93 = x10*x44 + x29*x47;
    double x94 = x10*x41 + x29*x40;
    double x95 = x13*x41 + x28*x40;
    double x96 = x13*x44 + x28*x47;
    double x97 = x11*x44 + x30*x47;
    double x98 = x11*x41 + x30*x40;
    double x99 = x12*x41 + x27*x40;
    double x100 = x12*x44 + x27*x47;
    
    res_0[0] = x38*(x42*x70 + x71*x79);
    res_0[1] = x38*(x42*x83 + x71*x84);
    res_0[2] = x38*(x42*x86 + x71*x85);
    res_0[3] = x38*(x70*x87 + x79*x88);
    res_0[4] = x38*(x83*x87 + x84*x88);
    res_0[5] = x38*(x85*x88 + x86*x87);
    res_0[6] = x38*(x70*x89 + x79*x90);
    res_0[7] = x38*(x83*x89 + x84*x90);
    res_0[8] = x38*(x85*x90 + x86*x89);
    res_0[9] = x38*(x70*x91 + x79*x92);
    res_0[10] = x38*(x83*x91 + x84*x92);
    res_0[11] = x38*(x85*x92 + x86*x91);
    res_0[12] = x38*(x70*x94 + x79*x93);
    res_0[13] = x38*(x83*x94 + x84*x93);
    res_0[14] = x38*(x85*x93 + x86*x94);
    res_0[15] = x38*(x70*x95 + x79*x96);
    res_0[16] = x38*(x83*x95 + x84*x96);
    res_0[17] = x38*(x85*x96 + x86*x95);
    res_0[18] = x38*(x70*x98 + x79*x97);
    res_0[19] = x38*(x83*x98 + x84*x97);
    res_0[20] = x38*(x85*x97 + x86*x98);
    res_0[21] = x38*(x100*x79 + x70*x99);
    res_0[22] = x38*(x100*x84 + x83*x99);
    res_0[23] = x38*(x100*x85 + x86*x99);
}

Conf_Forces_API void Integration_CPE8R_static_mbf(size_t num_elem,double Coords[][8][3],double Element_U[][8][3],double S[][4][6],
    double PENER[][4],double SENER[][4],double Conf_Force[][8][3]){
    
    //Inputs:
    //    Coords.........Element_Nodal_Coordinates [num_elem][numNodes,3]
    //    Element_U......Element_Nodal_Displacements [num_elem][numNodes,3]
    //    S..............Stress at integration points [num_elem][num_int_points,6]
    //    PENER..........Plastic strain energy [num_elem][num_int_points]
    //    SENER..........Elastic strain energy [num_elem][num_int_points]
    //    
    double weights[4]={1.0,1.0,1.0,1.0,};
    double int_points[4][3]={
        {-0.5773502691896258,-0.5773502691896258,0.0,},
        {0.5773502691896258,-0.5773502691896258,0.0,},
        {-0.5773502691896258,0.5773502691896258,0.0,},
        {0.5773502691896258,0.5773502691896258,0.0,}};

    #pragma omp parallel for
    for (size_t i=0;i<num_elem;i++){
        double TMP[8][3];
        for (size_t j=0;j<4;j++){
            CPE8R_static_mbf((double*) int_points[j],(double*) Coords[i],(double*) Element_U[i],(double*) S[i][j],PENER[i][j],SENER[i][j],(double*) TMP);
            for (size_t num_node=0;num_node<8;num_node++){
                for (size_t ii=0;ii<3;ii++){
                    Conf_Force[i][num_node][ii]+=TMP[num_node][ii]*weights[j];
                }
            }
        }
    }
}

void CPE8R_static_dbf(double *rst,double *coord,double *Element_U,double *S,double PENER,double SENER,double *res_0){
    double r = rst[0];
    double s = rst[1];
    double t = rst[2];
    double *U = Element_U;
    double S11 = S[0];
    double S22 = S[1];
    double S33 = S[2];
    double S12 = S[3];
    double S13 = S[4];
    double S23 = S[5];
    
    double x0 = (1.0/4.0)*s;
    double x1 = -x0;
    double x2 = -1.0/2.0*r;
    double x3 = r*s;
    double x4 = (1.0/2.0)*x3;
    double x5 = s*s;
    double x6 = (1.0/4.0)*x5;
    double x7 = x4 + x6;
    double x8 = -x1 - x2 - x7;
    double x9 = r*r;
    double x10 = (1.0/2.0)*x9 - 1.0/2.0;
    double x11 = -x10;
    double x12 = -s + x3;
    double x13 = -s - x3;
    double x14 = -1.0/2.0*s;
    double x15 = (1.0/4.0)*r;
    double x16 = (1.0/4.0)*x9;
    double x17 = x16 - x4;
    double x18 = -x14 - x15 - x17;
    double x19 = (1.0/2.0)*s;
    double x20 = -x15;
    double x21 = x17 + x19 + x20;
    double x22 = x16 + x4;
    double x23 = -x14 - x20 - x22;
    double x24 = x15 + x19 + x22;
    double x25 = x10*coord[13] + x11*coord[19] + x12*coord[22] + x13*coord[16] + x18*coord[4] + x21*coord[10] + x23*coord[1] + x24*coord[7];
    double x26 = x10*coord[12] + x11*coord[18] + x12*coord[21] + x13*coord[15] + x18*coord[3] + x21*coord[9] + x23*coord[0] + x24*coord[6];
    double x27 = (1.0/2.0)*x5 - 1.0/2.0;
    double x28 = -x27;
    double x29 = -r + x3;
    double x30 = -r - x3;
    double x31 = (1.0/2.0)*r;
    double x32 = x4 - x6;
    double x33 = x1 + x31 + x32;
    double x34 = -x0 - x2 - x32;
    double x35 = x0 + x31 + x7;
    double x36 = x27*coord[22] + x28*coord[16] + x29*coord[13] + x30*coord[19] + x33*coord[10] + x34*coord[4] + x35*coord[7] + x8*coord[1];
    double x37 = x27*coord[21] + x28*coord[15] + x29*coord[12] + x30*coord[18] + x33*coord[9] + x34*coord[3] + x35*coord[6] + x8*coord[0];
    double x38 = x25*x37 - x26*x36;
    double x39 = 1.0/x38;
    double x40 = x25*x39;
    double x41 = -x36*x39;
    double x42 = x23*x41 + x40*x8;
    double x43 = x10*U[13] + x11*U[19] + x12*U[22] + x13*U[16] + x18*U[4] + x21*U[10] + x23*U[1] + x24*U[7];
    double x44 = x37*x39;
    double x45 = x43*x44;
    double x46 = x27*U[22] + x28*U[16] + x29*U[13] + x30*U[19] + x33*U[10] + x34*U[4] + x35*U[7] + x8*U[1];
    double x47 = -x26*x39;
    double x48 = x46*x47;
    double x49 = 1.0*x45 + 1.0*x48 + 1.0;
    double x50 = x10*U[12] + x11*U[18] + x12*U[21] + x13*U[15] + x18*U[3] + x21*U[9] + x23*U[0] + x24*U[6];
    double x51 = x44*x50;
    double x52 = x27*U[21] + x28*U[15] + x29*U[12] + x30*U[18] + x33*U[9] + x34*U[3] + x35*U[6] + x8*U[0];
    double x53 = x47*x52;
    double x54 = -1.0*x51 - 1.0*x53;
    double x55 = S11*x49 + S12*x54;
    double x56 = x41*x50;
    double x57 = x40*x52;
    double x58 = x56 + x57;
    double x59 = S12*x49 + S22*x54;
    double x60 = x41*x43;
    double x61 = x40*x46;
    double x62 = x60 + x61;
    double x63 = S13*x49 + S23*x54;
    double x64 = x10*U[14] + x11*U[20] + x12*U[23] + x13*U[17] + x18*U[5] + x21*U[11] + x23*U[2] + x24*U[8];
    double x65 = x27*U[23] + x28*U[17] + x29*U[14] + x30*U[20] + x33*U[11] + x34*U[5] + x35*U[8] + x8*U[2];
    double x66 = x40*x65 + x41*x64;
    double x67 = -1.0*PENER - 1.0*SENER;
    double x68 = -x55*x58 - x59*x62 - x63*x66 - x67;
    double x69 = x23*x44 + x47*x8;
    double x70 = -1.0*x60 - 1.0*x61;
    double x71 = 1.0*x56 + 1.0*x57 + 1.0;
    double x72 = S11*x70 + S12*x71;
    double x73 = S12*x70 + S22*x71;
    double x74 = S13*x70 + S23*x71;
    double x75 = -x58*x72 - x62*x73 - x66*x74;
    double x76 = x51 + x53;
    double x77 = x45 + x48;
    double x78 = x44*x64 + x47*x65;
    double x79 = -x55*x76 - x59*x77 - x63*x78;
    double x80 = -x67 - x72*x76 - x73*x77 - x74*x78;
    double x81 = x18*x41 + x34*x40;
    double x82 = x18*x44 + x34*x47;
    double x83 = x24*x41 + x35*x40;
    double x84 = x24*x44 + x35*x47;
    double x85 = x21*x41 + x33*x40;
    double x86 = x21*x44 + x33*x47;
    double x87 = x10*x44 + x29*x47;
    double x88 = x10*x41 + x29*x40;
    double x89 = x13*x41 + x28*x40;
    double x90 = x13*x44 + x28*x47;
    double x91 = x11*x44 + x30*x47;
    double x92 = x11*x41 + x30*x40;
    double x93 = x12*x41 + x27*x40;
    double x94 = x12*x44 + x27*x47;
    
    res_0[0] = x38*(x42*x68 + x69*x75);
    res_0[1] = x38*(x42*x79 + x69*x80);
    res_0[2] = 0;
    res_0[3] = x38*(x68*x81 + x75*x82);
    res_0[4] = x38*(x79*x81 + x80*x82);
    res_0[5] = 0;
    res_0[6] = x38*(x68*x83 + x75*x84);
    res_0[7] = x38*(x79*x83 + x80*x84);
    res_0[8] = 0;
    res_0[9] = x38*(x68*x85 + x75*x86);
    res_0[10] = x38*(x79*x85 + x80*x86);
    res_0[11] = 0;
    res_0[12] = x38*(x68*x88 + x75*x87);
    res_0[13] = x38*(x79*x88 + x80*x87);
    res_0[14] = 0;
    res_0[15] = x38*(x68*x89 + x75*x90);
    res_0[16] = x38*(x79*x89 + x80*x90);
    res_0[17] = 0;
    res_0[18] = x38*(x68*x92 + x75*x91);
    res_0[19] = x38*(x79*x92 + x80*x91);
    res_0[20] = 0;
    res_0[21] = x38*(x68*x93 + x75*x94);
    res_0[22] = x38*(x79*x93 + x80*x94);
    res_0[23] = 0;
}

Conf_Forces_API void Integration_CPE8R_static_dbf(size_t num_elem,double Coords[][8][3],double Element_U[][8][3],double S[][4][6],
    double PENER[][4],double SENER[][4],double Conf_Force[][8][3]){
    
    //Inputs:
    //    Coords.........Element_Nodal_Coordinates [num_elem][numNodes,3]
    //    Element_U......Element_Nodal_Displacements [num_elem][numNodes,3]
    //    S..............Stress at integration points [num_elem][num_int_points,6]
    //    PENER..........Plastic strain energy [num_elem][num_int_points]
    //    SENER..........Elastic strain energy [num_elem][num_int_points]
    //    
    double weights[4]={1.0,1.0,1.0,1.0,};
    double int_points[4][3]={
        {-0.5773502691896258,-0.5773502691896258,0.0,},
        {0.5773502691896258,-0.5773502691896258,0.0,},
        {-0.5773502691896258,0.5773502691896258,0.0,},
        {0.5773502691896258,0.5773502691896258,0.0,}};

    #pragma omp parallel for
    for (size_t i=0;i<num_elem;i++){
        double TMP[8][3];
        for (size_t j=0;j<4;j++){
            CPE8R_static_dbf((double*) int_points[j],(double*) Coords[i],(double*) Element_U[i],(double*) S[i][j],PENER[i][j],SENER[i][j],(double*) TMP);
            for (size_t num_node=0;num_node<8;num_node++){
                for (size_t ii=0;ii<3;ii++){
                    Conf_Force[i][num_node][ii]+=TMP[num_node][ii]*weights[j];
                }
            }
        }
    }
}

void CPE3_dynamic(double *rst,double *coord,double rho,double *Element_U,double *Element_V,double *Element_A,double *S,double PENER,double SENER,double *res_0){
    double r = rst[0];
    double s = rst[1];
    double t = rst[2];
    double *U = Element_U;
    double *V = Element_V;
    double *A = Element_A;
    double S11 = S[0];
    double S22 = S[1];
    double S33 = S[2];
    double S12 = S[3];
    double S13 = S[4];
    double S23 = S[5];
    
    double x0 = -coord[0] + coord[3];
    double x1 = -coord[1] + coord[7];
    double x2 = coord[0] - coord[6];
    double x3 = coord[1] - coord[4];
    double x4 = x0*x1 - x2*x3;
    double x5 = 1.0/x4;
    double x6 = x0*x5;
    double x7 = x2*x5;
    double x8 = -x6 - x7;
    double x9 = -U[1] + U[4];
    double x10 = x1*x5;
    double x11 = x10*x9;
    double x12 = -U[1] + U[7];
    double x13 = x3*x5;
    double x14 = x12*x13;
    double x15 = -1.0*x11 - 1.0*x14;
    double x16 = -U[0] + U[3];
    double x17 = x10*x16;
    double x18 = -U[0] + U[6];
    double x19 = x13*x18;
    double x20 = 1.0*x17 + 1.0*x19 + 1.0;
    double x21 = S11*x15 + S12*x20;
    double x22 = x17 + x19;
    double x23 = S12*x15 + S22*x20;
    double x24 = x11 + x14;
    double x25 = S13*x15 + S23*x20;
    double x26 = -U[2] + U[5];
    double x27 = -U[2] + U[8];
    double x28 = x10*x26 + x13*x27;
    double x29 = -x21*x22 - x23*x24 - x25*x28;
    double x30 = -x10 - x13;
    double x31 = x7*x9;
    double x32 = x12*x6;
    double x33 = 1.0*x31 + 1.0*x32 + 1.0;
    double x34 = x16*x7;
    double x35 = x18*x6;
    double x36 = -1.0*x34 - 1.0*x35;
    double x37 = S11*x33 + S12*x36;
    double x38 = S12*x33 + S22*x36;
    double x39 = S13*x33 + S23*x36;
    double x40 = -r - s + 1;
    double x41 = r*V[3] + s*V[6] + x40*V[0];
    double x42 = r*V[4] + s*V[7] + x40*V[1];
    double x43 = r*V[5] + s*V[8] + x40*V[2];
    double x44 = -1.0*PENER - 1.0*SENER + 0.5*rho*(x41*x41 + x42*x42 + x43*x43);
    double x45 = -x22*x37 - x24*x38 - x28*x39 - x44;
    double x46 = r*A[3] + s*A[6] + x40*A[0];
    double x47 = r*A[4] + s*A[7] + x40*A[1];
    double x48 = r*A[5] + s*A[8] + x40*A[2];
    double x49 = -V[0] + V[3];
    double x50 = -V[0] + V[6];
    double x51 = -V[1] + V[4];
    double x52 = -V[1] + V[7];
    double x53 = -V[2] + V[5];
    double x54 = -V[2] + V[8];
    double x55 = x22*x46 + x24*x47 + x28*x48 + x41*(x10*x49 + x13*x50) + x42*(x10*x51 + x13*x52) + x43*(x10*x53 + x13*x54);
    double x56 = rho*x40;
    double x57 = x34 + x35;
    double x58 = x31 + x32;
    double x59 = x26*x7 + x27*x6;
    double x60 = -x21*x57 - x23*x58 - x25*x59 - x44;
    double x61 = -x37*x57 - x38*x58 - x39*x59;
    double x62 = x41*(x49*x7 + x50*x6) + x42*(x51*x7 + x52*x6) + x43*(x53*x7 + x54*x6) + x46*x57 + x47*x58 + x48*x59;
    double x63 = r*rho;
    double x64 = rho*s;
    
    res_0[0] = x4*(x29*x8 + x30*x45 - x55*x56);
    res_0[1] = x4*(x30*x61 - x56*x62 + x60*x8);
    res_0[2] = 0;
    res_0[3] = x4*(x10*x45 + x29*x7 - x55*x63);
    res_0[4] = x4*(x10*x61 + x60*x7 - x62*x63);
    res_0[5] = 0;
    res_0[6] = x4*(x13*x45 + x29*x6 - x55*x64);
    res_0[7] = x4*(x13*x61 + x6*x60 - x62*x64);
    res_0[8] = 0;
}

Conf_Forces_API void Integration_CPE3_dynamic(size_t num_elem,double Coords[][3][3],double *rho,double Element_U[][3][3],double Element_V[][3][3],double Element_A[][3][3],double S[][1][6],
    double PENER[][1],double SENER[][1],double Conf_Force[][3][3]){
    
    //Inputs:
    //    Coords.........Element_Nodal_Coordinates  [num_elem][numNodes,3]
    //    Element_U......Element_Nodal_Displacements[num_elem][numNodes,3]
    //    Element_V......Element_Nodal_Velocity     [num_elem][numNodes,3]
    //    Element_A......Element_Nodal_Acceleration [num_elem][numNodes,3]
    //    S..............Stress at integration points [num_elem][num_int_points,6]
    //    PENER..........Plastic strain energy [num_elem][num_int_points]
    //    SENER..........Elastic strain energy [num_elem][num_int_points]
    //    
    double weights[1]={0.5,};
    double int_points[1][3]={
        {0.3333333333333333,0.3333333333333333,0.0,}};

    #pragma omp parallel for
    for (size_t i=0;i<num_elem;i++){
        double TMP[3][3];
        for (size_t j=0;j<1;j++){
            CPE3_dynamic((double*) int_points[j],(double*) Coords[i],rho[i],(double*) Element_U[i],(double*) Element_V[i],(double*) Element_A[i],(double*) S[i][j],PENER[i][j],SENER[i][j],(double*) TMP);
            for (size_t num_node=0;num_node<3;num_node++){
                for (size_t ii=0;ii<3;ii++){
                    Conf_Force[i][num_node][ii]+=TMP[num_node][ii]*weights[j];
                }
            }
        }
    }
}

void CPE3_static_mbf(double *rst,double *coord,double *Element_U,double *S,double PENER,double SENER,double *res_0){
    double r = rst[0];
    double s = rst[1];
    double t = rst[2];
    double *U = Element_U;
    double S11 = S[0];
    double S22 = S[1];
    double S33 = S[2];
    double S12 = S[3];
    double S13 = S[4];
    double S23 = S[5];
    
    double x0 = -coord[0] + coord[3];
    double x1 = -coord[1] + coord[7];
    double x2 = coord[0] - coord[6];
    double x3 = coord[1] - coord[4];
    double x4 = x0*x1 - x2*x3;
    double x5 = 1.0/x4;
    double x6 = x0*x5;
    double x7 = x2*x5;
    double x8 = -x6 - x7;
    double x9 = -U[1] + U[4];
    double x10 = x1*x5;
    double x11 = x10*x9;
    double x12 = -U[1] + U[7];
    double x13 = x3*x5;
    double x14 = x12*x13;
    double x15 = -1.0*x11 - 1.0*x14;
    double x16 = -U[0] + U[3];
    double x17 = x10*x16;
    double x18 = -U[0] + U[6];
    double x19 = x13*x18;
    double x20 = 1.0*x17 + 1.0*x19 + 1.0;
    double x21 = S11*x15 + S12*x20;
    double x22 = x17 + x19 + 1.0;
    double x23 = S12*x15 + S22*x20;
    double x24 = x11 + x14;
    double x25 = S13*x15;
    double x26 = S23*x20;
    double x27 = x25 + x26;
    double x28 = -U[2] + U[5];
    double x29 = -U[2] + U[8];
    double x30 = x10*x28 + x13*x29;
    double x31 = -x21*x22 - x23*x24 - x27*x30;
    double x32 = -x10 - x13;
    double x33 = x7*x9;
    double x34 = x12*x6;
    double x35 = 1.0*x33 + 1.0*x34 + 1.0;
    double x36 = x16*x7;
    double x37 = x18*x6;
    double x38 = -1.0*x36 - 1.0*x37;
    double x39 = S11*x35 + S12*x38;
    double x40 = S12*x35 + S22*x38;
    double x41 = S13*x35;
    double x42 = S23*x38;
    double x43 = x41 + x42;
    double x44 = -1.0*PENER - 1.0*SENER;
    double x45 = -x22*x39 - x24*x40 - x30*x43 - x44;
    double x46 = x36 + x37;
    double x47 = x33 + x34 + 1.0;
    double x48 = x28*x7 + x29*x6;
    double x49 = -x21*x46 - x23*x47 - x27*x48 - x44;
    double x50 = -x39*x46 - x40*x47 - x43*x48;
    double x51 = -1.0*x25 - 1.0*x26;
    double x52 = -1.0*x41 - 1.0*x42;
    
    res_0[0] = x4*(x31*x8 + x32*x45);
    res_0[1] = x4*(x32*x50 + x49*x8);
    res_0[2] = x4*(x32*x52 + x51*x8);
    res_0[3] = x4*(x10*x45 + x31*x7);
    res_0[4] = x4*(x10*x50 + x49*x7);
    res_0[5] = x4*(x10*x52 + x51*x7);
    res_0[6] = x4*(x13*x45 + x31*x6);
    res_0[7] = x4*(x13*x50 + x49*x6);
    res_0[8] = x4*(x13*x52 + x51*x6);
}

Conf_Forces_API void Integration_CPE3_static_mbf(size_t num_elem,double Coords[][3][3],double Element_U[][3][3],double S[][1][6],
    double PENER[][1],double SENER[][1],double Conf_Force[][3][3]){
    
    //Inputs:
    //    Coords.........Element_Nodal_Coordinates [num_elem][numNodes,3]
    //    Element_U......Element_Nodal_Displacements [num_elem][numNodes,3]
    //    S..............Stress at integration points [num_elem][num_int_points,6]
    //    PENER..........Plastic strain energy [num_elem][num_int_points]
    //    SENER..........Elastic strain energy [num_elem][num_int_points]
    //    
    double weights[1]={0.5,};
    double int_points[1][3]={
        {0.3333333333333333,0.3333333333333333,0.0,}};

    #pragma omp parallel for
    for (size_t i=0;i<num_elem;i++){
        double TMP[3][3];
        for (size_t j=0;j<1;j++){
            CPE3_static_mbf((double*) int_points[j],(double*) Coords[i],(double*) Element_U[i],(double*) S[i][j],PENER[i][j],SENER[i][j],(double*) TMP);
            for (size_t num_node=0;num_node<3;num_node++){
                for (size_t ii=0;ii<3;ii++){
                    Conf_Force[i][num_node][ii]+=TMP[num_node][ii]*weights[j];
                }
            }
        }
    }
}

void CPE3_static_dbf(double *rst,double *coord,double *Element_U,double *S,double PENER,double SENER,double *res_0){
    double r = rst[0];
    double s = rst[1];
    double t = rst[2];
    double *U = Element_U;
    double S11 = S[0];
    double S22 = S[1];
    double S33 = S[2];
    double S12 = S[3];
    double S13 = S[4];
    double S23 = S[5];
    
    double x0 = -coord[0] + coord[3];
    double x1 = -coord[1] + coord[7];
    double x2 = coord[0] - coord[6];
    double x3 = coord[1] - coord[4];
    double x4 = x0*x1 - x2*x3;
    double x5 = 1.0/x4;
    double x6 = x0*x5;
    double x7 = x2*x5;
    double x8 = -x6 - x7;
    double x9 = -U[1] + U[4];
    double x10 = x1*x5;
    double x11 = x10*x9;
    double x12 = -U[1] + U[7];
    double x13 = x3*x5;
    double x14 = x12*x13;
    double x15 = -1.0*x11 - 1.0*x14;
    double x16 = -U[0] + U[3];
    double x17 = x10*x16;
    double x18 = -U[0] + U[6];
    double x19 = x13*x18;
    double x20 = 1.0*x17 + 1.0*x19 + 1.0;
    double x21 = S11*x15 + S12*x20;
    double x22 = x17 + x19;
    double x23 = S12*x15 + S22*x20;
    double x24 = x11 + x14;
    double x25 = S13*x15 + S23*x20;
    double x26 = -U[2] + U[5];
    double x27 = -U[2] + U[8];
    double x28 = x10*x26 + x13*x27;
    double x29 = -x21*x22 - x23*x24 - x25*x28;
    double x30 = -x10 - x13;
    double x31 = x7*x9;
    double x32 = x12*x6;
    double x33 = 1.0*x31 + 1.0*x32 + 1.0;
    double x34 = x16*x7;
    double x35 = x18*x6;
    double x36 = -1.0*x34 - 1.0*x35;
    double x37 = S11*x33 + S12*x36;
    double x38 = S12*x33 + S22*x36;
    double x39 = S13*x33 + S23*x36;
    double x40 = -1.0*PENER - 1.0*SENER;
    double x41 = -x22*x37 - x24*x38 - x28*x39 - x40;
    double x42 = x34 + x35;
    double x43 = x31 + x32;
    double x44 = x26*x7 + x27*x6;
    double x45 = -x21*x42 - x23*x43 - x25*x44 - x40;
    double x46 = -x37*x42 - x38*x43 - x39*x44;
    
    res_0[0] = x4*(x29*x8 + x30*x41);
    res_0[1] = x4*(x30*x46 + x45*x8);
    res_0[2] = 0;
    res_0[3] = x4*(x10*x41 + x29*x7);
    res_0[4] = x4*(x10*x46 + x45*x7);
    res_0[5] = 0;
    res_0[6] = x4*(x13*x41 + x29*x6);
    res_0[7] = x4*(x13*x46 + x45*x6);
    res_0[8] = 0;
}

Conf_Forces_API void Integration_CPE3_static_dbf(size_t num_elem,double Coords[][3][3],double Element_U[][3][3],double S[][1][6],
    double PENER[][1],double SENER[][1],double Conf_Force[][3][3]){
    
    //Inputs:
    //    Coords.........Element_Nodal_Coordinates [num_elem][numNodes,3]
    //    Element_U......Element_Nodal_Displacements [num_elem][numNodes,3]
    //    S..............Stress at integration points [num_elem][num_int_points,6]
    //    PENER..........Plastic strain energy [num_elem][num_int_points]
    //    SENER..........Elastic strain energy [num_elem][num_int_points]
    //    
    double weights[1]={0.5,};
    double int_points[1][3]={
        {0.3333333333333333,0.3333333333333333,0.0,}};

    #pragma omp parallel for
    for (size_t i=0;i<num_elem;i++){
        double TMP[3][3];
        for (size_t j=0;j<1;j++){
            CPE3_static_dbf((double*) int_points[j],(double*) Coords[i],(double*) Element_U[i],(double*) S[i][j],PENER[i][j],SENER[i][j],(double*) TMP);
            for (size_t num_node=0;num_node<3;num_node++){
                for (size_t ii=0;ii<3;ii++){
                    Conf_Force[i][num_node][ii]+=TMP[num_node][ii]*weights[j];
                }
            }
        }
    }
}

void CPE6_dynamic(double *rst,double *coord,double rho,double *Element_U,double *Element_V,double *Element_A,double *S,double PENER,double SENER,double *res_0){
    double r = rst[0];
    double s = rst[1];
    double t = rst[2];
    double *U = Element_U;
    double *V = Element_V;
    double *A = Element_A;
    double S11 = S[0];
    double S22 = S[1];
    double S33 = S[2];
    double S12 = S[3];
    double S13 = S[4];
    double S23 = S[5];
    
    double x0 = 4*r;
    double x1 = x0 - 1;
    double x2 = 4*s;
    double x3 = x0 + x2 - 3;
    double x4 = x3*coord[1];
    double x5 = -8*r - x2 + 4;
    double x6 = x1*coord[4] + x2*coord[13] - x2*coord[16] + x4 + x5*coord[10];
    double x7 = x2 - 1;
    double x8 = x3*coord[0];
    double x9 = -8*s - x0 + 4;
    double x10 = -x0*coord[9] + x0*coord[12] + x7*coord[6] + x8 + x9*coord[15];
    double x11 = x1*coord[3] + x2*coord[12] - x2*coord[15] + x5*coord[9] + x8;
    double x12 = -x0*coord[10] + x0*coord[13] + x4 + x7*coord[7] + x9*coord[16];
    double x13 = -x10*x6 + x11*x12;
    double x14 = 1.0/x13;
    double x15 = x14*x3;
    double x16 = -x10;
    double x17 = x11*x15 + x15*x16;
    double x18 = x3*U[1];
    double x19 = x1*U[4] + x18 + x2*U[13] - x2*U[16] + x5*U[10];
    double x20 = x12*x14;
    double x21 = x19*x20;
    double x22 = -x0*U[10] + x0*U[13] + x18 + x7*U[7] + x9*U[16];
    double x23 = -x6;
    double x24 = x14*x23;
    double x25 = x22*x24;
    double x26 = -1.0*x21 - 1.0*x25;
    double x27 = x3*U[0];
    double x28 = x1*U[3] + x2*U[12] - x2*U[15] + x27 + x5*U[9];
    double x29 = x20*x28;
    double x30 = -x0*U[9] + x0*U[12] + x27 + x7*U[6] + x9*U[15];
    double x31 = x24*x30;
    double x32 = 1.0*x29 + 1.0*x31 + 1.0;
    double x33 = S11*x26 + S12*x32;
    double x34 = x29 + x31;
    double x35 = S12*x26 + S22*x32;
    double x36 = x21 + x25;
    double x37 = S13*x26 + S23*x32;
    double x38 = x3*U[2];
    double x39 = x1*U[5] + x2*U[14] - x2*U[17] + x38 + x5*U[11];
    double x40 = -x0*U[11] + x0*U[14] + x38 + x7*U[8] + x9*U[17];
    double x41 = x20*x39 + x24*x40;
    double x42 = -x33*x34 - x35*x36 - x37*x41;
    double x43 = x12*x15 + x15*x23;
    double x44 = x14*x16;
    double x45 = x19*x44;
    double x46 = x11*x14;
    double x47 = x22*x46;
    double x48 = 1.0*x45 + 1.0*x47 + 1.0;
    double x49 = x28*x44;
    double x50 = x30*x46;
    double x51 = -1.0*x49 - 1.0*x50;
    double x52 = S11*x48 + S12*x51;
    double x53 = S12*x48 + S22*x51;
    double x54 = S13*x48 + S23*x51;
    double x55 = r*r;
    double x56 = 2*x55;
    double x57 = -r + x56;
    double x58 = s*s;
    double x59 = 2*x58;
    double x60 = -s + x59;
    double x61 = s*x0;
    double x62 = x0 - 4*x55 - x61;
    double x63 = x2 - 4*x58 - x61;
    double x64 = -3*r - 3*s + x56 + x59 + x61 + 1;
    double x65 = x57*V[3] + x60*V[6] + x61*V[12] + x62*V[9] + x63*V[15] + x64*V[0];
    double x66 = x57*V[4] + x60*V[7] + x61*V[13] + x62*V[10] + x63*V[16] + x64*V[1];
    double x67 = x57*V[5] + x60*V[8] + x61*V[14] + x62*V[11] + x63*V[17] + x64*V[2];
    double x68 = -1.0*PENER - 1.0*SENER + 0.5*rho*(x65*x65 + x66*x66 + x67*x67);
    double x69 = -x34*x52 - x36*x53 - x41*x54 - x68;
    double x70 = x57*A[3] + x60*A[6] + x61*A[12] + x62*A[9] + x63*A[15] + x64*A[0];
    double x71 = x57*A[4] + x60*A[7] + x61*A[13] + x62*A[10] + x63*A[16] + x64*A[1];
    double x72 = x57*A[5] + x60*A[8] + x61*A[14] + x62*A[11] + x63*A[17] + x64*A[2];
    double x73 = x3*V[0];
    double x74 = x1*V[3] + x2*V[12] - x2*V[15] + x5*V[9] + x73;
    double x75 = -x0*V[9] + x0*V[12] + x7*V[6] + x73 + x9*V[15];
    double x76 = x3*V[1];
    double x77 = x1*V[4] + x2*V[13] - x2*V[16] + x5*V[10] + x76;
    double x78 = -x0*V[10] + x0*V[13] + x7*V[7] + x76 + x9*V[16];
    double x79 = x3*V[2];
    double x80 = x1*V[5] + x2*V[14] - x2*V[17] + x5*V[11] + x79;
    double x81 = -x0*V[11] + x0*V[14] + x7*V[8] + x79 + x9*V[17];
    double x82 = x34*x70 + x36*x71 + x41*x72 + x65*(x20*x74 + x24*x75) + x66*(x20*x77 + x24*x78) + x67*(x20*x80 + x24*x81);
    double x83 = rho*x64;
    double x84 = x49 + x50;
    double x85 = x45 + x47;
    double x86 = x39*x44 + x40*x46;
    double x87 = -x33*x84 - x35*x85 - x37*x86 - x68;
    double x88 = -x52*x84 - x53*x85 - x54*x86;
    double x89 = x65*(x44*x74 + x46*x75) + x66*(x44*x77 + x46*x78) + x67*(x44*x80 + x46*x81) + x70*x84 + x71*x85 + x72*x86;
    double x90 = rho*x57;
    double x91 = x1*x44;
    double x92 = x1*x20;
    double x93 = rho*x60;
    double x94 = x46*x7;
    double x95 = x24*x7;
    double x96 = x0*x24;
    double x97 = x12*x14*x5 - x96;
    double x98 = x0*x46;
    double x99 = x14*x16*x5 - x98;
    double x100 = rho*x62;
    double x101 = x2*x44;
    double x102 = x101 + x98;
    double x103 = x2*x20;
    double x104 = x103 + x96;
    double x105 = rho*x61;
    double x106 = -x101 + x11*x14*x9;
    double x107 = -x103 + x14*x23*x9;
    double x108 = rho*x63;
    
    res_0[0] = x13*(x17*x42 + x43*x69 - x82*x83);
    res_0[1] = x13*(x17*x87 + x43*x88 - x83*x89);
    res_0[2] = 0;
    res_0[3] = x13*(x42*x91 + x69*x92 - x82*x90);
    res_0[4] = x13*(x87*x91 + x88*x92 - x89*x90);
    res_0[5] = 0;
    res_0[6] = x13*(x42*x94 + x69*x95 - x82*x93);
    res_0[7] = x13*(x87*x94 + x88*x95 - x89*x93);
    res_0[8] = 0;
    res_0[9] = x13*(-x100*x82 + x42*x99 + x69*x97);
    res_0[10] = x13*(-x100*x89 + x87*x99 + x88*x97);
    res_0[11] = 0;
    res_0[12] = x13*(x102*x42 + x104*x69 - x105*x82);
    res_0[13] = x13*(x102*x87 + x104*x88 - x105*x89);
    res_0[14] = 0;
    res_0[15] = x13*(x106*x42 + x107*x69 - x108*x82);
    res_0[16] = x13*(x106*x87 + x107*x88 - x108*x89);
    res_0[17] = 0;
}

Conf_Forces_API void Integration_CPE6_dynamic(size_t num_elem,double Coords[][6][3],double *rho,double Element_U[][6][3],double Element_V[][6][3],double Element_A[][6][3],double S[][3][6],
    double PENER[][3],double SENER[][3],double Conf_Force[][6][3]){
    
    //Inputs:
    //    Coords.........Element_Nodal_Coordinates  [num_elem][numNodes,3]
    //    Element_U......Element_Nodal_Displacements[num_elem][numNodes,3]
    //    Element_V......Element_Nodal_Velocity     [num_elem][numNodes,3]
    //    Element_A......Element_Nodal_Acceleration [num_elem][numNodes,3]
    //    S..............Stress at integration points [num_elem][num_int_points,6]
    //    PENER..........Plastic strain energy [num_elem][num_int_points]
    //    SENER..........Elastic strain energy [num_elem][num_int_points]
    //    
    double weights[3]={0.16666666666666666,0.16666666666666666,0.16666666666666666,};
    double int_points[3][3]={
        {0.16666666666666666,0.16666666666666666,0.0,},
        {0.6666666666666666,0.16666666666666666,0.0,},
        {0.16666666666666666,0.6666666666666666,0.0,}};

    #pragma omp parallel for
    for (size_t i=0;i<num_elem;i++){
        double TMP[6][3];
        for (size_t j=0;j<3;j++){
            CPE6_dynamic((double*) int_points[j],(double*) Coords[i],rho[i],(double*) Element_U[i],(double*) Element_V[i],(double*) Element_A[i],(double*) S[i][j],PENER[i][j],SENER[i][j],(double*) TMP);
            for (size_t num_node=0;num_node<6;num_node++){
                for (size_t ii=0;ii<3;ii++){
                    Conf_Force[i][num_node][ii]+=TMP[num_node][ii]*weights[j];
                }
            }
        }
    }
}

void CPE6_static_mbf(double *rst,double *coord,double *Element_U,double *S,double PENER,double SENER,double *res_0){
    double r = rst[0];
    double s = rst[1];
    double t = rst[2];
    double *U = Element_U;
    double S11 = S[0];
    double S22 = S[1];
    double S33 = S[2];
    double S12 = S[3];
    double S13 = S[4];
    double S23 = S[5];
    
    double x0 = 4*r;
    double x1 = x0 - 1;
    double x2 = 4*s;
    double x3 = x0 + x2 - 3;
    double x4 = x3*coord[0];
    double x5 = -8*r - x2 + 4;
    double x6 = x1*coord[3] + x2*coord[12] - x2*coord[15] + x4 + x5*coord[9];
    double x7 = x3*coord[1];
    double x8 = x1*coord[4] + x2*coord[13] - x2*coord[16] + x5*coord[10] + x7;
    double x9 = x2 - 1;
    double x10 = -8*s - x0 + 4;
    double x11 = -x0*coord[9] + x0*coord[12] + x10*coord[15] + x4 + x9*coord[6];
    double x12 = -x0*coord[10] + x0*coord[13] + x10*coord[16] + x7 + x9*coord[7];
    double x13 = -x11*x8 + x12*x6;
    double x14 = 1.0/x13;
    double x15 = x14*x3;
    double x16 = -x11;
    double x17 = x15*x16 + x15*x6;
    double x18 = x3*U[1];
    double x19 = x1*U[4] + x18 + x2*U[13] - x2*U[16] + x5*U[10];
    double x20 = x12*x14;
    double x21 = x19*x20;
    double x22 = -x0*U[10] + x0*U[13] + x10*U[16] + x18 + x9*U[7];
    double x23 = -x8;
    double x24 = x14*x23;
    double x25 = x22*x24;
    double x26 = -1.0*x21 - 1.0*x25;
    double x27 = x3*U[0];
    double x28 = x1*U[3] + x2*U[12] - x2*U[15] + x27 + x5*U[9];
    double x29 = x20*x28;
    double x30 = -x0*U[9] + x0*U[12] + x10*U[15] + x27 + x9*U[6];
    double x31 = x24*x30;
    double x32 = 1.0*x29 + 1.0*x31 + 1.0;
    double x33 = S11*x26 + S12*x32;
    double x34 = x29 + x31 + 1.0;
    double x35 = S12*x26 + S22*x32;
    double x36 = x21 + x25;
    double x37 = S13*x26;
    double x38 = S23*x32;
    double x39 = x37 + x38;
    double x40 = x3*U[2];
    double x41 = x1*U[5] + x2*U[14] - x2*U[17] + x40 + x5*U[11];
    double x42 = -x0*U[11] + x0*U[14] + x10*U[17] + x40 + x9*U[8];
    double x43 = x20*x41 + x24*x42;
    double x44 = -x33*x34 - x35*x36 - x39*x43;
    double x45 = x12*x15 + x15*x23;
    double x46 = x14*x16;
    double x47 = x19*x46;
    double x48 = x14*x6;
    double x49 = x22*x48;
    double x50 = 1.0*x47 + 1.0*x49 + 1.0;
    double x51 = x28*x46;
    double x52 = x30*x48;
    double x53 = -1.0*x51 - 1.0*x52;
    double x54 = S11*x50 + S12*x53;
    double x55 = S12*x50 + S22*x53;
    double x56 = S13*x50;
    double x57 = S23*x53;
    double x58 = x56 + x57;
    double x59 = -1.0*PENER - 1.0*SENER;
    double x60 = -x34*x54 - x36*x55 - x43*x58 - x59;
    double x61 = x51 + x52;
    double x62 = x47 + x49 + 1.0;
    double x63 = x41*x46 + x42*x48;
    double x64 = -x33*x61 - x35*x62 - x39*x63 - x59;
    double x65 = -x54*x61 - x55*x62 - x58*x63;
    double x66 = -1.0*x37 - 1.0*x38;
    double x67 = -1.0*x56 - 1.0*x57;
    double x68 = x1*x46;
    double x69 = x1*x20;
    double x70 = x48*x9;
    double x71 = x24*x9;
    double x72 = x0*x24;
    double x73 = x12*x14*x5 - x72;
    double x74 = x0*x48;
    double x75 = x14*x16*x5 - x74;
    double x76 = x2*x46;
    double x77 = x74 + x76;
    double x78 = x2*x20;
    double x79 = x72 + x78;
    double x80 = x10*x14*x6 - x76;
    double x81 = x10*x14*x23 - x78;
    
    res_0[0] = x13*(x17*x44 + x45*x60);
    res_0[1] = x13*(x17*x64 + x45*x65);
    res_0[2] = x13*(x17*x66 + x45*x67);
    res_0[3] = x13*(x44*x68 + x60*x69);
    res_0[4] = x13*(x64*x68 + x65*x69);
    res_0[5] = x13*(x66*x68 + x67*x69);
    res_0[6] = x13*(x44*x70 + x60*x71);
    res_0[7] = x13*(x64*x70 + x65*x71);
    res_0[8] = x13*(x66*x70 + x67*x71);
    res_0[9] = x13*(x44*x75 + x60*x73);
    res_0[10] = x13*(x64*x75 + x65*x73);
    res_0[11] = x13*(x66*x75 + x67*x73);
    res_0[12] = x13*(x44*x77 + x60*x79);
    res_0[13] = x13*(x64*x77 + x65*x79);
    res_0[14] = x13*(x66*x77 + x67*x79);
    res_0[15] = x13*(x44*x80 + x60*x81);
    res_0[16] = x13*(x64*x80 + x65*x81);
    res_0[17] = x13*(x66*x80 + x67*x81);
}

Conf_Forces_API void Integration_CPE6_static_mbf(size_t num_elem,double Coords[][6][3],double Element_U[][6][3],double S[][3][6],
    double PENER[][3],double SENER[][3],double Conf_Force[][6][3]){
    
    //Inputs:
    //    Coords.........Element_Nodal_Coordinates [num_elem][numNodes,3]
    //    Element_U......Element_Nodal_Displacements [num_elem][numNodes,3]
    //    S..............Stress at integration points [num_elem][num_int_points,6]
    //    PENER..........Plastic strain energy [num_elem][num_int_points]
    //    SENER..........Elastic strain energy [num_elem][num_int_points]
    //    
    double weights[3]={0.16666666666666666,0.16666666666666666,0.16666666666666666,};
    double int_points[3][3]={
        {0.16666666666666666,0.16666666666666666,0.0,},
        {0.6666666666666666,0.16666666666666666,0.0,},
        {0.16666666666666666,0.6666666666666666,0.0,}};

    #pragma omp parallel for
    for (size_t i=0;i<num_elem;i++){
        double TMP[6][3];
        for (size_t j=0;j<3;j++){
            CPE6_static_mbf((double*) int_points[j],(double*) Coords[i],(double*) Element_U[i],(double*) S[i][j],PENER[i][j],SENER[i][j],(double*) TMP);
            for (size_t num_node=0;num_node<6;num_node++){
                for (size_t ii=0;ii<3;ii++){
                    Conf_Force[i][num_node][ii]+=TMP[num_node][ii]*weights[j];
                }
            }
        }
    }
}

void CPE6_static_dbf(double *rst,double *coord,double *Element_U,double *S,double PENER,double SENER,double *res_0){
    double r = rst[0];
    double s = rst[1];
    double t = rst[2];
    double *U = Element_U;
    double S11 = S[0];
    double S22 = S[1];
    double S33 = S[2];
    double S12 = S[3];
    double S13 = S[4];
    double S23 = S[5];
    
    double x0 = 4*r;
    double x1 = x0 - 1;
    double x2 = 4*s;
    double x3 = x0 + x2 - 3;
    double x4 = x3*coord[0];
    double x5 = -8*r - x2 + 4;
    double x6 = x1*coord[3] + x2*coord[12] - x2*coord[15] + x4 + x5*coord[9];
    double x7 = x3*coord[1];
    double x8 = x1*coord[4] + x2*coord[13] - x2*coord[16] + x5*coord[10] + x7;
    double x9 = x2 - 1;
    double x10 = -8*s - x0 + 4;
    double x11 = -x0*coord[9] + x0*coord[12] + x10*coord[15] + x4 + x9*coord[6];
    double x12 = -x0*coord[10] + x0*coord[13] + x10*coord[16] + x7 + x9*coord[7];
    double x13 = -x11*x8 + x12*x6;
    double x14 = 1.0/x13;
    double x15 = x14*x3;
    double x16 = -x11;
    double x17 = x15*x16 + x15*x6;
    double x18 = x3*U[1];
    double x19 = x1*U[4] + x18 + x2*U[13] - x2*U[16] + x5*U[10];
    double x20 = x12*x14;
    double x21 = x19*x20;
    double x22 = -x0*U[10] + x0*U[13] + x10*U[16] + x18 + x9*U[7];
    double x23 = -x8;
    double x24 = x14*x23;
    double x25 = x22*x24;
    double x26 = -1.0*x21 - 1.0*x25;
    double x27 = x3*U[0];
    double x28 = x1*U[3] + x2*U[12] - x2*U[15] + x27 + x5*U[9];
    double x29 = x20*x28;
    double x30 = -x0*U[9] + x0*U[12] + x10*U[15] + x27 + x9*U[6];
    double x31 = x24*x30;
    double x32 = 1.0*x29 + 1.0*x31 + 1.0;
    double x33 = S11*x26 + S12*x32;
    double x34 = x29 + x31;
    double x35 = S12*x26 + S22*x32;
    double x36 = x21 + x25;
    double x37 = S13*x26 + S23*x32;
    double x38 = x3*U[2];
    double x39 = x1*U[5] + x2*U[14] - x2*U[17] + x38 + x5*U[11];
    double x40 = -x0*U[11] + x0*U[14] + x10*U[17] + x38 + x9*U[8];
    double x41 = x20*x39 + x24*x40;
    double x42 = -x33*x34 - x35*x36 - x37*x41;
    double x43 = x12*x15 + x15*x23;
    double x44 = x14*x16;
    double x45 = x19*x44;
    double x46 = x14*x6;
    double x47 = x22*x46;
    double x48 = 1.0*x45 + 1.0*x47 + 1.0;
    double x49 = x28*x44;
    double x50 = x30*x46;
    double x51 = -1.0*x49 - 1.0*x50;
    double x52 = S11*x48 + S12*x51;
    double x53 = S12*x48 + S22*x51;
    double x54 = S13*x48 + S23*x51;
    double x55 = -1.0*PENER - 1.0*SENER;
    double x56 = -x34*x52 - x36*x53 - x41*x54 - x55;
    double x57 = x49 + x50;
    double x58 = x45 + x47;
    double x59 = x39*x44 + x40*x46;
    double x60 = -x33*x57 - x35*x58 - x37*x59 - x55;
    double x61 = -x52*x57 - x53*x58 - x54*x59;
    double x62 = x1*x44;
    double x63 = x1*x20;
    double x64 = x46*x9;
    double x65 = x24*x9;
    double x66 = x0*x24;
    double x67 = x12*x14*x5 - x66;
    double x68 = x0*x46;
    double x69 = x14*x16*x5 - x68;
    double x70 = x2*x44;
    double x71 = x68 + x70;
    double x72 = x2*x20;
    double x73 = x66 + x72;
    double x74 = x10*x14*x6 - x70;
    double x75 = x10*x14*x23 - x72;
    
    res_0[0] = x13*(x17*x42 + x43*x56);
    res_0[1] = x13*(x17*x60 + x43*x61);
    res_0[2] = 0;
    res_0[3] = x13*(x42*x62 + x56*x63);
    res_0[4] = x13*(x60*x62 + x61*x63);
    res_0[5] = 0;
    res_0[6] = x13*(x42*x64 + x56*x65);
    res_0[7] = x13*(x60*x64 + x61*x65);
    res_0[8] = 0;
    res_0[9] = x13*(x42*x69 + x56*x67);
    res_0[10] = x13*(x60*x69 + x61*x67);
    res_0[11] = 0;
    res_0[12] = x13*(x42*x71 + x56*x73);
    res_0[13] = x13*(x60*x71 + x61*x73);
    res_0[14] = 0;
    res_0[15] = x13*(x42*x74 + x56*x75);
    res_0[16] = x13*(x60*x74 + x61*x75);
    res_0[17] = 0;
}

Conf_Forces_API void Integration_CPE6_static_dbf(size_t num_elem,double Coords[][6][3],double Element_U[][6][3],double S[][3][6],
    double PENER[][3],double SENER[][3],double Conf_Force[][6][3]){
    
    //Inputs:
    //    Coords.........Element_Nodal_Coordinates [num_elem][numNodes,3]
    //    Element_U......Element_Nodal_Displacements [num_elem][numNodes,3]
    //    S..............Stress at integration points [num_elem][num_int_points,6]
    //    PENER..........Plastic strain energy [num_elem][num_int_points]
    //    SENER..........Elastic strain energy [num_elem][num_int_points]
    //    
    double weights[3]={0.16666666666666666,0.16666666666666666,0.16666666666666666,};
    double int_points[3][3]={
        {0.16666666666666666,0.16666666666666666,0.0,},
        {0.6666666666666666,0.16666666666666666,0.0,},
        {0.16666666666666666,0.6666666666666666,0.0,}};

    #pragma omp parallel for
    for (size_t i=0;i<num_elem;i++){
        double TMP[6][3];
        for (size_t j=0;j<3;j++){
            CPE6_static_dbf((double*) int_points[j],(double*) Coords[i],(double*) Element_U[i],(double*) S[i][j],PENER[i][j],SENER[i][j],(double*) TMP);
            for (size_t num_node=0;num_node<6;num_node++){
                for (size_t ii=0;ii<3;ii++){
                    Conf_Force[i][num_node][ii]+=TMP[num_node][ii]*weights[j];
                }
            }
        }
    }
}

void C3D6_dynamic(double *rst,double *coord,double rho,double *Element_U,double *Element_V,double *Element_A,double *S,double PENER,double SENER,double *res_0){
    double r = rst[0];
    double s = rst[1];
    double t = rst[2];
    double *U = Element_U;
    double *V = Element_V;
    double *A = Element_A;
    double S11 = S[0];
    double S22 = S[1];
    double S33 = S[2];
    double S12 = S[3];
    double S13 = S[4];
    double S23 = S[5];
    
    double x0 = (1.0/2.0)*t;
    double x1 = x0 + 1.0/2.0;
    double x2 = x0 - 1.0/2.0;
    double x3 = -x2;
    double x4 = -x1;
    double x5 = x2*coord[0] + x4*coord[9];
    double x6 = x1*coord[12] + x3*coord[3] + x5;
    double x7 = x2*coord[1] + x4*coord[10];
    double x8 = x1*coord[16] + x3*coord[7] + x7;
    double x9 = x6*x8;
    double x10 = x1*coord[13] + x3*coord[4] + x7;
    double x11 = x1*coord[15] + x3*coord[6] + x5;
    double x12 = x10*x11;
    double x13 = -x12 + x9;
    double x14 = (1.0/2.0)*r;
    double x15 = (1.0/2.0)*s;
    double x16 = x14 + x15 - 1.0/2.0;
    double x17 = -x16;
    double x18 = -x14*coord[5] + x14*coord[14] - x15*coord[8] + x15*coord[17] + x16*coord[2] + x17*coord[11];
    double x19 = x2*coord[2] + x4*coord[11];
    double x20 = x1*coord[17] + x19 + x3*coord[8];
    double x21 = -x14*coord[3] + x14*coord[12] - x15*coord[6] + x15*coord[15] + x16*coord[0] + x17*coord[9];
    double x22 = x1*coord[14] + x19 + x3*coord[5];
    double x23 = -x14*coord[4] + x14*coord[13] - x15*coord[7] + x15*coord[16] + x16*coord[1] + x17*coord[10];
    double x24 = x11*x23;
    double x25 = x23*x6;
    double x26 = x21*x8;
    double x27 = x10*x20*x21 - x12*x18 + x18*x9 - x20*x25 + x22*x24 - x22*x26;
    double x28 = 1.0/x27;
    double x29 = x28*(-x14*U[4] + x14*U[13] - x15*U[7] + x15*U[16] + x16*U[1] + x17*U[10]);
    double x30 = x10*x21 - x25;
    double x31 = x2*U[1] + x4*U[10];
    double x32 = x28*(x1*U[16] + x3*U[7] + x31);
    double x33 = x24 - x26;
    double x34 = x28*(x1*U[13] + x3*U[4] + x31);
    double x35 = x13*x29 + x30*x32 + x33*x34;
    double x36 = x10*x20 - x22*x8;
    double x37 = x28*(-x14*U[5] + x14*U[14] - x15*U[8] + x15*U[17] + x16*U[2] + x17*U[11]);
    double x38 = -x10*x18 + x22*x23;
    double x39 = x2*U[2] + x4*U[11];
    double x40 = x28*(x1*U[17] + x3*U[8] + x39);
    double x41 = x18*x8 - x20*x23;
    double x42 = x28*(x1*U[14] + x3*U[5] + x39);
    double x43 = x36*x37 + x38*x40 + x41*x42;
    double x44 = x29*x36 + x32*x38 + x34*x41;
    double x45 = x13*x37 + x30*x40 + x33*x42;
    double x46 = x45 + 1.0;
    double x47 = x35*x43 - x44*x46;
    double x48 = x28*(-x14*U[3] + x14*U[12] - x15*U[6] + x15*U[15] + x16*U[0] + x17*U[9]);
    double x49 = x2*U[0] + x4*U[9];
    double x50 = x28*(x1*U[15] + x3*U[6] + x49);
    double x51 = x28*(x1*U[12] + x3*U[3] + x49);
    double x52 = x13*x48 + x30*x50 + x33*x51;
    double x53 = x36*x48 + x38*x50 + x41*x51;
    double x54 = x53 + 1.0;
    double x55 = -x43*x52 + x46*x54;
    double x56 = -x35*x54 + x44*x52;
    double x57 = S11*x47 + S12*x55 + S13*x56;
    double x58 = S12*x47 + S22*x55 + S23*x56;
    double x59 = S13*x47 + S23*x55 + S33*x56;
    double x60 = -x43*x59 - x44*x58 - x53*x57;
    double x61 = x18*x6 - x21*x22;
    double x62 = x2*x28;
    double x63 = -x11*x18 + x20*x21;
    double x64 = x11*x22 - x20*x6;
    double x65 = x16*x28;
    double x66 = x61*x62 + x62*x63 + x64*x65;
    double x67 = x37*x64 + x40*x61 + x42*x63;
    double x68 = x29*x64 + x32*x61 + x34*x63;
    double x69 = x68 + 1.0;
    double x70 = -x43*x69 + x44*x67;
    double x71 = x48*x64 + x50*x61 + x51*x63;
    double x72 = x43*x71 - x54*x67;
    double x73 = -x44*x71 + x54*x69;
    double x74 = S11*x70 + S12*x72 + S13*x73;
    double x75 = S12*x70 + S22*x72 + S23*x73;
    double x76 = S13*x70 + S23*x72 + S33*x73;
    double x77 = -x43*x76 - x44*x75 - x53*x74;
    double x78 = x13*x65 + x30*x62 + x33*x62;
    double x79 = x36*x65 + x38*x62 + x41*x62;
    double x80 = -x35*x67 + x46*x69;
    double x81 = -x46*x71 + x52*x67;
    double x82 = x35*x71 - x52*x69;
    double x83 = S11*x80 + S12*x81 + S13*x82;
    double x84 = S12*x80 + S22*x81 + S23*x82;
    double x85 = S13*x80 + S23*x81 + S33*x82;
    double x86 = r*x0;
    double x87 = x14 + x86;
    double x88 = -x14 + x86;
    double x89 = -x88;
    double x90 = s*x0;
    double x91 = x15 + x90;
    double x92 = -x15 + x90;
    double x93 = -x92;
    double x94 = -x0;
    double x95 = -x87 - x91 - x94 + 1.0/2.0;
    double x96 = x88 + x92 + x94 + 1.0/2.0;
    double x97 = x87*V[12] + x89*V[3] + x91*V[15] + x93*V[6] + x95*V[9] + x96*V[0];
    double x98 = x87*V[13] + x89*V[4] + x91*V[16] + x93*V[7] + x95*V[10] + x96*V[1];
    double x99 = x87*V[14] + x89*V[5] + x91*V[17] + x93*V[8] + x95*V[11] + x96*V[2];
    double x100 = -1.0*PENER - 1.0*SENER + 0.5*rho*(x97*x97 + x98*x98 + x99*x99);
    double x101 = -x100 - x43*x85 - x44*x84 - x53*x83;
    double x102 = x87*A[12] + x89*A[3] + x91*A[15] + x93*A[6] + x95*A[9] + x96*A[0];
    double x103 = x87*A[13] + x89*A[4] + x91*A[16] + x93*A[7] + x95*A[10] + x96*A[1];
    double x104 = x87*A[14] + x89*A[5] + x91*A[17] + x93*A[8] + x95*A[11] + x96*A[2];
    double x105 = -x14*V[3] + x14*V[12] - x15*V[6] + x15*V[15] + x16*V[0] + x17*V[9];
    double x106 = x28*x36;
    double x107 = x2*V[0] + x4*V[9];
    double x108 = x1*V[15] + x107 + x3*V[6];
    double x109 = x28*x38;
    double x110 = x1*V[12] + x107 + x3*V[3];
    double x111 = x28*x41;
    double x112 = -x14*V[4] + x14*V[13] - x15*V[7] + x15*V[16] + x16*V[1] + x17*V[10];
    double x113 = x2*V[1] + x4*V[10];
    double x114 = x1*V[16] + x113 + x3*V[7];
    double x115 = x1*V[13] + x113 + x3*V[4];
    double x116 = -x14*V[5] + x14*V[14] - x15*V[8] + x15*V[17] + x16*V[2] + x17*V[11];
    double x117 = x2*V[2] + x4*V[11];
    double x118 = x1*V[17] + x117 + x3*V[8];
    double x119 = x1*V[14] + x117 + x3*V[5];
    double x120 = x102*x53 + x103*x44 + x104*x43 + x97*(x105*x106 + x108*x109 + x110*x111) + x98*(x106*x112 + x109*x114 + x111*x115) + x99*(x106*x116 + x109*x118 + x111*x119);
    double x121 = rho*x96;
    double x122 = -x67*x76 - x68*x75 - x71*x74;
    double x123 = -x67*x85 - x68*x84 - x71*x83;
    double x124 = -x100 - x57*x71 - x58*x68 - x59*x67;
    double x125 = x28*x61;
    double x126 = x28*x64;
    double x127 = x28*x63;
    double x128 = x102*x71 + x103*x68 + x104*x67 + x97*(x105*x126 + x108*x125 + x110*x127) + x98*(x112*x126 + x114*x125 + x115*x127) + x99*(x116*x126 + x118*x125 + x119*x127);
    double x129 = -x35*x58 - x45*x59 - x52*x57;
    double x130 = -x35*x84 - x45*x85 - x52*x83;
    double x131 = -x100 - x35*x75 - x45*x76 - x52*x74;
    double x132 = x13*x28;
    double x133 = x28*x30;
    double x134 = x28*x33;
    double x135 = x102*x52 + x103*x35 + x104*x45 + x97*(x105*x132 + x108*x133 + x110*x134) + x98*(x112*x132 + x114*x133 + x115*x134) + x99*(x116*x132 + x118*x133 + x119*x134);
    double x136 = x132*x14;
    double x137 = -x136 + x28*x3*x33;
    double x138 = x106*x14;
    double x139 = -x138 + x28*x3*x41;
    double x140 = x126*x14;
    double x141 = -x140 + x28*x3*x63;
    double x142 = rho*x89;
    double x143 = x126*x15;
    double x144 = -x143 + x28*x3*x61;
    double x145 = x132*x15;
    double x146 = -x145 + x28*x3*x30;
    double x147 = x106*x15;
    double x148 = -x147 + x28*x3*x38;
    double x149 = rho*x93;
    double x150 = x125*x4 + x126*x17 + x127*x4;
    double x151 = x132*x17 + x133*x4 + x134*x4;
    double x152 = x106*x17 + x109*x4 + x111*x4;
    double x153 = rho*x95;
    double x154 = x1*x134 + x136;
    double x155 = x1*x111 + x138;
    double x156 = x1*x127 + x140;
    double x157 = rho*x87;
    double x158 = x1*x125 + x143;
    double x159 = x1*x133 + x145;
    double x160 = x1*x109 + x147;
    double x161 = rho*x91;
    
    res_0[0] = x27*(x101*x79 - x120*x121 + x60*x66 + x77*x78);
    res_0[1] = x27*(-x121*x128 + x122*x78 + x123*x79 + x124*x66);
    res_0[2] = x27*(-x121*x135 + x129*x66 + x130*x79 + x131*x78);
    res_0[3] = x27*(x101*x139 - x120*x142 + x137*x77 + x141*x60);
    res_0[4] = x27*(x122*x137 + x123*x139 + x124*x141 - x128*x142);
    res_0[5] = x27*(x129*x141 + x130*x139 + x131*x137 - x135*x142);
    res_0[6] = x27*(x101*x148 - x120*x149 + x144*x60 + x146*x77);
    res_0[7] = x27*(x122*x146 + x123*x148 + x124*x144 - x128*x149);
    res_0[8] = x27*(x129*x144 + x130*x148 + x131*x146 - x135*x149);
    res_0[9] = x27*(x101*x152 - x120*x153 + x150*x60 + x151*x77);
    res_0[10] = x27*(x122*x151 + x123*x152 + x124*x150 - x128*x153);
    res_0[11] = x27*(x129*x150 + x130*x152 + x131*x151 - x135*x153);
    res_0[12] = x27*(x101*x155 - x120*x157 + x154*x77 + x156*x60);
    res_0[13] = x27*(x122*x154 + x123*x155 + x124*x156 - x128*x157);
    res_0[14] = x27*(x129*x156 + x130*x155 + x131*x154 - x135*x157);
    res_0[15] = x27*(x101*x160 - x120*x161 + x158*x60 + x159*x77);
    res_0[16] = x27*(x122*x159 + x123*x160 + x124*x158 - x128*x161);
    res_0[17] = x27*(x129*x158 + x130*x160 + x131*x159 - x135*x161);
}

Conf_Forces_API void Integration_C3D6_dynamic(size_t num_elem,double Coords[][6][3],double *rho,double Element_U[][6][3],double Element_V[][6][3],double Element_A[][6][3],double S[][2][6],
    double PENER[][2],double SENER[][2],double Conf_Force[][6][3]){
    
    //Inputs:
    //    Coords.........Element_Nodal_Coordinates  [num_elem][numNodes,3]
    //    Element_U......Element_Nodal_Displacements[num_elem][numNodes,3]
    //    Element_V......Element_Nodal_Velocity     [num_elem][numNodes,3]
    //    Element_A......Element_Nodal_Acceleration [num_elem][numNodes,3]
    //    S..............Stress at integration points [num_elem][num_int_points,6]
    //    PENER..........Plastic strain energy [num_elem][num_int_points]
    //    SENER..........Elastic strain energy [num_elem][num_int_points]
    //    
    double weights[2]={0.5,0.5,};
    double int_points[2][3]={
        {0.3333333333333333,0.3333333333333333,-0.5773502691896258,},
        {0.3333333333333333,0.3333333333333333,0.5773502691896258,}};

    #pragma omp parallel for
    for (size_t i=0;i<num_elem;i++){
        double TMP[6][3];
        for (size_t j=0;j<2;j++){
            C3D6_dynamic((double*) int_points[j],(double*) Coords[i],rho[i],(double*) Element_U[i],(double*) Element_V[i],(double*) Element_A[i],(double*) S[i][j],PENER[i][j],SENER[i][j],(double*) TMP);
            for (size_t num_node=0;num_node<6;num_node++){
                for (size_t ii=0;ii<3;ii++){
                    Conf_Force[i][num_node][ii]+=TMP[num_node][ii]*weights[j];
                }
            }
        }
    }
}

void C3D6_static_mbf(double *rst,double *coord,double *Element_U,double *S,double PENER,double SENER,double *res_0){
    double r = rst[0];
    double s = rst[1];
    double t = rst[2];
    double *U = Element_U;
    double S11 = S[0];
    double S22 = S[1];
    double S33 = S[2];
    double S12 = S[3];
    double S13 = S[4];
    double S23 = S[5];
    
    double x0 = (1.0/2.0)*t;
    double x1 = x0 + 1.0/2.0;
    double x2 = x0 - 1.0/2.0;
    double x3 = -x2;
    double x4 = -x1;
    double x5 = x2*coord[0] + x4*coord[9];
    double x6 = x1*coord[12] + x3*coord[3] + x5;
    double x7 = x2*coord[1] + x4*coord[10];
    double x8 = x1*coord[16] + x3*coord[7] + x7;
    double x9 = x6*x8;
    double x10 = x1*coord[13] + x3*coord[4] + x7;
    double x11 = x1*coord[15] + x3*coord[6] + x5;
    double x12 = x10*x11;
    double x13 = -x12 + x9;
    double x14 = (1.0/2.0)*r;
    double x15 = (1.0/2.0)*s;
    double x16 = x14 + x15 - 1.0/2.0;
    double x17 = -x16;
    double x18 = -x14*coord[5] + x14*coord[14] - x15*coord[8] + x15*coord[17] + x16*coord[2] + x17*coord[11];
    double x19 = x2*coord[2] + x4*coord[11];
    double x20 = x1*coord[17] + x19 + x3*coord[8];
    double x21 = -x14*coord[3] + x14*coord[12] - x15*coord[6] + x15*coord[15] + x16*coord[0] + x17*coord[9];
    double x22 = x1*coord[14] + x19 + x3*coord[5];
    double x23 = -x14*coord[4] + x14*coord[13] - x15*coord[7] + x15*coord[16] + x16*coord[1] + x17*coord[10];
    double x24 = x11*x23;
    double x25 = x23*x6;
    double x26 = x21*x8;
    double x27 = x10*x20*x21 - x12*x18 + x18*x9 - x20*x25 + x22*x24 - x22*x26;
    double x28 = 1.0/x27;
    double x29 = x28*(-x14*U[4] + x14*U[13] - x15*U[7] + x15*U[16] + x16*U[1] + x17*U[10]);
    double x30 = x10*x21 - x25;
    double x31 = x2*U[1] + x4*U[10];
    double x32 = x28*(x1*U[16] + x3*U[7] + x31);
    double x33 = x24 - x26;
    double x34 = x28*(x1*U[13] + x3*U[4] + x31);
    double x35 = x13*x29 + x30*x32 + x33*x34;
    double x36 = x10*x20 - x22*x8;
    double x37 = x28*(-x14*U[5] + x14*U[14] - x15*U[8] + x15*U[17] + x16*U[2] + x17*U[11]);
    double x38 = -x10*x18 + x22*x23;
    double x39 = x2*U[2] + x4*U[11];
    double x40 = x28*(x1*U[17] + x3*U[8] + x39);
    double x41 = x18*x8 - x20*x23;
    double x42 = x28*(x1*U[14] + x3*U[5] + x39);
    double x43 = x36*x37 + x38*x40 + x41*x42;
    double x44 = x29*x36 + x32*x38 + x34*x41;
    double x45 = x13*x37 + x30*x40 + x33*x42 + 1.0;
    double x46 = x35*x43 - x44*x45;
    double x47 = x28*(-x14*U[3] + x14*U[12] - x15*U[6] + x15*U[15] + x16*U[0] + x17*U[9]);
    double x48 = x2*U[0] + x4*U[9];
    double x49 = x28*(x1*U[15] + x3*U[6] + x48);
    double x50 = x28*(x1*U[12] + x3*U[3] + x48);
    double x51 = x13*x47 + x30*x49 + x33*x50;
    double x52 = x36*x47 + x38*x49 + x41*x50 + 1.0;
    double x53 = -x43*x51 + x45*x52;
    double x54 = -x35*x52 + x44*x51;
    double x55 = S11*x46 + S12*x53 + S13*x54;
    double x56 = S12*x46 + S22*x53 + S23*x54;
    double x57 = S13*x46 + S23*x53 + S33*x54;
    double x58 = -x43*x57 - x44*x56 - x52*x55;
    double x59 = x18*x6 - x21*x22;
    double x60 = x2*x28;
    double x61 = -x11*x18 + x20*x21;
    double x62 = x11*x22 - x20*x6;
    double x63 = x16*x28;
    double x64 = x59*x60 + x60*x61 + x62*x63;
    double x65 = x37*x62 + x40*x59 + x42*x61;
    double x66 = x29*x62 + x32*x59 + x34*x61 + 1.0;
    double x67 = -x43*x66 + x44*x65;
    double x68 = x47*x62 + x49*x59 + x50*x61;
    double x69 = x43*x68 - x52*x65;
    double x70 = -x44*x68 + x52*x66;
    double x71 = S11*x67 + S12*x69 + S13*x70;
    double x72 = S12*x67 + S22*x69 + S23*x70;
    double x73 = S13*x67 + S23*x69 + S33*x70;
    double x74 = -x43*x73 - x44*x72 - x52*x71;
    double x75 = x13*x63 + x30*x60 + x33*x60;
    double x76 = x36*x63 + x38*x60 + x41*x60;
    double x77 = -x35*x65 + x45*x66;
    double x78 = -x45*x68 + x51*x65;
    double x79 = x35*x68 - x51*x66;
    double x80 = S11*x77 + S12*x78 + S13*x79;
    double x81 = S12*x77 + S22*x78 + S23*x79;
    double x82 = S13*x77 + S23*x78 + S33*x79;
    double x83 = -1.0*PENER - 1.0*SENER;
    double x84 = -x43*x82 - x44*x81 - x52*x80 - x83;
    double x85 = -x65*x73 - x66*x72 - x68*x71;
    double x86 = -x65*x82 - x66*x81 - x68*x80;
    double x87 = -x55*x68 - x56*x66 - x57*x65 - x83;
    double x88 = -x35*x56 - x45*x57 - x51*x55;
    double x89 = -x35*x81 - x45*x82 - x51*x80;
    double x90 = -x35*x72 - x45*x73 - x51*x71 - x83;
    double x91 = x14*x28;
    double x92 = x13*x91;
    double x93 = x28*x3*x33 - x92;
    double x94 = x36*x91;
    double x95 = x28*x3*x41 - x94;
    double x96 = x62*x91;
    double x97 = x28*x3*x61 - x96;
    double x98 = x15*x28;
    double x99 = x62*x98;
    double x100 = x28*x3*x59 - x99;
    double x101 = x13*x98;
    double x102 = -x101 + x28*x3*x30;
    double x103 = x36*x98;
    double x104 = -x103 + x28*x3*x38;
    double x105 = x28*x4;
    double x106 = x17*x28;
    double x107 = x105*x59 + x105*x61 + x106*x62;
    double x108 = x105*x30 + x105*x33 + x106*x13;
    double x109 = x105*x38 + x105*x41 + x106*x36;
    double x110 = x1*x28;
    double x111 = x110*x33 + x92;
    double x112 = x110*x41 + x94;
    double x113 = x110*x61 + x96;
    double x114 = x110*x59 + x99;
    double x115 = x101 + x110*x30;
    double x116 = x103 + x110*x38;
    
    res_0[0] = x27*(x58*x64 + x74*x75 + x76*x84);
    res_0[1] = x27*(x64*x87 + x75*x85 + x76*x86);
    res_0[2] = x27*(x64*x88 + x75*x90 + x76*x89);
    res_0[3] = x27*(x58*x97 + x74*x93 + x84*x95);
    res_0[4] = x27*(x85*x93 + x86*x95 + x87*x97);
    res_0[5] = x27*(x88*x97 + x89*x95 + x90*x93);
    res_0[6] = x27*(x100*x58 + x102*x74 + x104*x84);
    res_0[7] = x27*(x100*x87 + x102*x85 + x104*x86);
    res_0[8] = x27*(x100*x88 + x102*x90 + x104*x89);
    res_0[9] = x27*(x107*x58 + x108*x74 + x109*x84);
    res_0[10] = x27*(x107*x87 + x108*x85 + x109*x86);
    res_0[11] = x27*(x107*x88 + x108*x90 + x109*x89);
    res_0[12] = x27*(x111*x74 + x112*x84 + x113*x58);
    res_0[13] = x27*(x111*x85 + x112*x86 + x113*x87);
    res_0[14] = x27*(x111*x90 + x112*x89 + x113*x88);
    res_0[15] = x27*(x114*x58 + x115*x74 + x116*x84);
    res_0[16] = x27*(x114*x87 + x115*x85 + x116*x86);
    res_0[17] = x27*(x114*x88 + x115*x90 + x116*x89);
}

Conf_Forces_API void Integration_C3D6_static_mbf(size_t num_elem,double Coords[][6][3],double Element_U[][6][3],double S[][2][6],
    double PENER[][2],double SENER[][2],double Conf_Force[][6][3]){
    
    //Inputs:
    //    Coords.........Element_Nodal_Coordinates [num_elem][numNodes,3]
    //    Element_U......Element_Nodal_Displacements [num_elem][numNodes,3]
    //    S..............Stress at integration points [num_elem][num_int_points,6]
    //    PENER..........Plastic strain energy [num_elem][num_int_points]
    //    SENER..........Elastic strain energy [num_elem][num_int_points]
    //    
    double weights[2]={0.5,0.5,};
    double int_points[2][3]={
        {0.3333333333333333,0.3333333333333333,-0.5773502691896258,},
        {0.3333333333333333,0.3333333333333333,0.5773502691896258,}};

    #pragma omp parallel for
    for (size_t i=0;i<num_elem;i++){
        double TMP[6][3];
        for (size_t j=0;j<2;j++){
            C3D6_static_mbf((double*) int_points[j],(double*) Coords[i],(double*) Element_U[i],(double*) S[i][j],PENER[i][j],SENER[i][j],(double*) TMP);
            for (size_t num_node=0;num_node<6;num_node++){
                for (size_t ii=0;ii<3;ii++){
                    Conf_Force[i][num_node][ii]+=TMP[num_node][ii]*weights[j];
                }
            }
        }
    }
}

void C3D6_static_dbf(double *rst,double *coord,double *Element_U,double *S,double PENER,double SENER,double *res_0){
    double r = rst[0];
    double s = rst[1];
    double t = rst[2];
    double *U = Element_U;
    double S11 = S[0];
    double S22 = S[1];
    double S33 = S[2];
    double S12 = S[3];
    double S13 = S[4];
    double S23 = S[5];
    
    double x0 = (1.0/2.0)*t;
    double x1 = x0 + 1.0/2.0;
    double x2 = x0 - 1.0/2.0;
    double x3 = -x2;
    double x4 = -x1;
    double x5 = x2*coord[0] + x4*coord[9];
    double x6 = x1*coord[12] + x3*coord[3] + x5;
    double x7 = x2*coord[1] + x4*coord[10];
    double x8 = x1*coord[16] + x3*coord[7] + x7;
    double x9 = x6*x8;
    double x10 = x1*coord[13] + x3*coord[4] + x7;
    double x11 = x1*coord[15] + x3*coord[6] + x5;
    double x12 = x10*x11;
    double x13 = -x12 + x9;
    double x14 = (1.0/2.0)*r;
    double x15 = (1.0/2.0)*s;
    double x16 = x14 + x15 - 1.0/2.0;
    double x17 = -x16;
    double x18 = -x14*coord[5] + x14*coord[14] - x15*coord[8] + x15*coord[17] + x16*coord[2] + x17*coord[11];
    double x19 = x2*coord[2] + x4*coord[11];
    double x20 = x1*coord[17] + x19 + x3*coord[8];
    double x21 = -x14*coord[3] + x14*coord[12] - x15*coord[6] + x15*coord[15] + x16*coord[0] + x17*coord[9];
    double x22 = x1*coord[14] + x19 + x3*coord[5];
    double x23 = -x14*coord[4] + x14*coord[13] - x15*coord[7] + x15*coord[16] + x16*coord[1] + x17*coord[10];
    double x24 = x11*x23;
    double x25 = x23*x6;
    double x26 = x21*x8;
    double x27 = x10*x20*x21 - x12*x18 + x18*x9 - x20*x25 + x22*x24 - x22*x26;
    double x28 = 1.0/x27;
    double x29 = x28*(-x14*U[4] + x14*U[13] - x15*U[7] + x15*U[16] + x16*U[1] + x17*U[10]);
    double x30 = x10*x21 - x25;
    double x31 = x2*U[1] + x4*U[10];
    double x32 = x28*(x1*U[16] + x3*U[7] + x31);
    double x33 = x24 - x26;
    double x34 = x28*(x1*U[13] + x3*U[4] + x31);
    double x35 = x13*x29 + x30*x32 + x33*x34;
    double x36 = x10*x20 - x22*x8;
    double x37 = x28*(-x14*U[5] + x14*U[14] - x15*U[8] + x15*U[17] + x16*U[2] + x17*U[11]);
    double x38 = -x10*x18 + x22*x23;
    double x39 = x2*U[2] + x4*U[11];
    double x40 = x28*(x1*U[17] + x3*U[8] + x39);
    double x41 = x18*x8 - x20*x23;
    double x42 = x28*(x1*U[14] + x3*U[5] + x39);
    double x43 = x36*x37 + x38*x40 + x41*x42;
    double x44 = x29*x36 + x32*x38 + x34*x41;
    double x45 = x13*x37 + x30*x40 + x33*x42;
    double x46 = x45 + 1.0;
    double x47 = x35*x43 - x44*x46;
    double x48 = x28*(-x14*U[3] + x14*U[12] - x15*U[6] + x15*U[15] + x16*U[0] + x17*U[9]);
    double x49 = x2*U[0] + x4*U[9];
    double x50 = x28*(x1*U[15] + x3*U[6] + x49);
    double x51 = x28*(x1*U[12] + x3*U[3] + x49);
    double x52 = x13*x48 + x30*x50 + x33*x51;
    double x53 = x36*x48 + x38*x50 + x41*x51;
    double x54 = x53 + 1.0;
    double x55 = -x43*x52 + x46*x54;
    double x56 = -x35*x54 + x44*x52;
    double x57 = S11*x47 + S12*x55 + S13*x56;
    double x58 = S12*x47 + S22*x55 + S23*x56;
    double x59 = S13*x47 + S23*x55 + S33*x56;
    double x60 = -x43*x59 - x44*x58 - x53*x57;
    double x61 = x18*x6 - x21*x22;
    double x62 = x2*x28;
    double x63 = -x11*x18 + x20*x21;
    double x64 = x11*x22 - x20*x6;
    double x65 = x16*x28;
    double x66 = x61*x62 + x62*x63 + x64*x65;
    double x67 = x37*x64 + x40*x61 + x42*x63;
    double x68 = x29*x64 + x32*x61 + x34*x63;
    double x69 = x68 + 1.0;
    double x70 = -x43*x69 + x44*x67;
    double x71 = x48*x64 + x50*x61 + x51*x63;
    double x72 = x43*x71 - x54*x67;
    double x73 = -x44*x71 + x54*x69;
    double x74 = S11*x70 + S12*x72 + S13*x73;
    double x75 = S12*x70 + S22*x72 + S23*x73;
    double x76 = S13*x70 + S23*x72 + S33*x73;
    double x77 = -x43*x76 - x44*x75 - x53*x74;
    double x78 = x13*x65 + x30*x62 + x33*x62;
    double x79 = x36*x65 + x38*x62 + x41*x62;
    double x80 = -x35*x67 + x46*x69;
    double x81 = -x46*x71 + x52*x67;
    double x82 = x35*x71 - x52*x69;
    double x83 = S11*x80 + S12*x81 + S13*x82;
    double x84 = S12*x80 + S22*x81 + S23*x82;
    double x85 = S13*x80 + S23*x81 + S33*x82;
    double x86 = -1.0*PENER - 1.0*SENER;
    double x87 = -x43*x85 - x44*x84 - x53*x83 - x86;
    double x88 = -x67*x76 - x68*x75 - x71*x74;
    double x89 = -x67*x85 - x68*x84 - x71*x83;
    double x90 = -x57*x71 - x58*x68 - x59*x67 - x86;
    double x91 = -x35*x58 - x45*x59 - x52*x57;
    double x92 = -x35*x84 - x45*x85 - x52*x83;
    double x93 = -x35*x75 - x45*x76 - x52*x74 - x86;
    double x94 = x14*x28;
    double x95 = x13*x94;
    double x96 = x28*x3*x33 - x95;
    double x97 = x36*x94;
    double x98 = x28*x3*x41 - x97;
    double x99 = x64*x94;
    double x100 = x28*x3*x63 - x99;
    double x101 = x15*x28;
    double x102 = x101*x64;
    double x103 = -x102 + x28*x3*x61;
    double x104 = x101*x13;
    double x105 = -x104 + x28*x3*x30;
    double x106 = x101*x36;
    double x107 = -x106 + x28*x3*x38;
    double x108 = x28*x4;
    double x109 = x17*x28;
    double x110 = x108*x61 + x108*x63 + x109*x64;
    double x111 = x108*x30 + x108*x33 + x109*x13;
    double x112 = x108*x38 + x108*x41 + x109*x36;
    double x113 = x1*x28;
    double x114 = x113*x33 + x95;
    double x115 = x113*x41 + x97;
    double x116 = x113*x63 + x99;
    double x117 = x102 + x113*x61;
    double x118 = x104 + x113*x30;
    double x119 = x106 + x113*x38;
    
    res_0[0] = x27*(x60*x66 + x77*x78 + x79*x87);
    res_0[1] = x27*(x66*x90 + x78*x88 + x79*x89);
    res_0[2] = x27*(x66*x91 + x78*x93 + x79*x92);
    res_0[3] = x27*(x100*x60 + x77*x96 + x87*x98);
    res_0[4] = x27*(x100*x90 + x88*x96 + x89*x98);
    res_0[5] = x27*(x100*x91 + x92*x98 + x93*x96);
    res_0[6] = x27*(x103*x60 + x105*x77 + x107*x87);
    res_0[7] = x27*(x103*x90 + x105*x88 + x107*x89);
    res_0[8] = x27*(x103*x91 + x105*x93 + x107*x92);
    res_0[9] = x27*(x110*x60 + x111*x77 + x112*x87);
    res_0[10] = x27*(x110*x90 + x111*x88 + x112*x89);
    res_0[11] = x27*(x110*x91 + x111*x93 + x112*x92);
    res_0[12] = x27*(x114*x77 + x115*x87 + x116*x60);
    res_0[13] = x27*(x114*x88 + x115*x89 + x116*x90);
    res_0[14] = x27*(x114*x93 + x115*x92 + x116*x91);
    res_0[15] = x27*(x117*x60 + x118*x77 + x119*x87);
    res_0[16] = x27*(x117*x90 + x118*x88 + x119*x89);
    res_0[17] = x27*(x117*x91 + x118*x93 + x119*x92);
}

Conf_Forces_API void Integration_C3D6_static_dbf(size_t num_elem,double Coords[][6][3],double Element_U[][6][3],double S[][2][6],
    double PENER[][2],double SENER[][2],double Conf_Force[][6][3]){
    
    //Inputs:
    //    Coords.........Element_Nodal_Coordinates [num_elem][numNodes,3]
    //    Element_U......Element_Nodal_Displacements [num_elem][numNodes,3]
    //    S..............Stress at integration points [num_elem][num_int_points,6]
    //    PENER..........Plastic strain energy [num_elem][num_int_points]
    //    SENER..........Elastic strain energy [num_elem][num_int_points]
    //    
    double weights[2]={0.5,0.5,};
    double int_points[2][3]={
        {0.3333333333333333,0.3333333333333333,-0.5773502691896258,},
        {0.3333333333333333,0.3333333333333333,0.5773502691896258,}};

    #pragma omp parallel for
    for (size_t i=0;i<num_elem;i++){
        double TMP[6][3];
        for (size_t j=0;j<2;j++){
            C3D6_static_dbf((double*) int_points[j],(double*) Coords[i],(double*) Element_U[i],(double*) S[i][j],PENER[i][j],SENER[i][j],(double*) TMP);
            for (size_t num_node=0;num_node<6;num_node++){
                for (size_t ii=0;ii<3;ii++){
                    Conf_Force[i][num_node][ii]+=TMP[num_node][ii]*weights[j];
                }
            }
        }
    }
}

void C3D8_dynamic(double *rst,double *coord,double rho,double *Element_U,double *Element_V,double *Element_A,double *S,double PENER,double SENER,double *res_0){
    double r = rst[0];
    double s = rst[1];
    double t = rst[2];
    double *U = Element_U;
    double *V = Element_V;
    double *A = Element_A;
    double S11 = S[0];
    double S22 = S[1];
    double S33 = S[2];
    double S12 = S[3];
    double S13 = S[4];
    double S23 = S[5];
    
    double x0 = (1.0/8.0)*s;
    double x1 = (1.0/8.0)*r;
    double x2 = s*x1;
    double x3 = x0 + x2;
    double x4 = x1 + 1.0/8.0;
    double x5 = x3 + x4;
    double x6 = -x5;
    double x7 = -x1;
    double x8 = x7 - 1.0/8.0;
    double x9 = x3 + x8;
    double x10 = -x0;
    double x11 = x10 + x2;
    double x12 = x1 - 1.0/8.0;
    double x13 = x11 + x12;
    double x14 = x7 + 1.0/8.0;
    double x15 = x11 + x14;
    double x16 = -x15;
    double x17 = -x13;
    double x18 = -x9;
    double x19 = x13*coord[9] + x15*coord[12] + x16*coord[0] + x17*coord[21] + x18*coord[15] + x5*coord[18] + x6*coord[6] + x9*coord[3];
    double x20 = (1.0/8.0)*t;
    double x21 = t*x1;
    double x22 = x20 + x21;
    double x23 = x22 + x4;
    double x24 = -x23;
    double x25 = x22 + x8;
    double x26 = -x20;
    double x27 = x21 + x26;
    double x28 = x12 + x27;
    double x29 = x14 + x27;
    double x30 = -x29;
    double x31 = -x28;
    double x32 = -x25;
    double x33 = x23*coord[20] + x24*coord[17] + x25*coord[5] + x28*coord[14] + x29*coord[11] + x30*coord[2] + x31*coord[23] + x32*coord[8];
    double x34 = x13*coord[11] + x15*coord[14] + x16*coord[2] + x17*coord[23] + x18*coord[17] + x5*coord[20] + x6*coord[8] + x9*coord[5];
    double x35 = x23*coord[18] + x24*coord[15] + x25*coord[3] + x28*coord[12] + x29*coord[9] + x30*coord[0] + x31*coord[21] + x32*coord[6];
    double x36 = x34*x35;
    double x37 = x19*x33 - x36;
    double x38 = x23*coord[19] + x24*coord[16] + x25*coord[4] + x28*coord[13] + x29*coord[10] + x30*coord[1] + x31*coord[22] + x32*coord[7];
    double x39 = t*x0;
    double x40 = x20 + x39;
    double x41 = x0 + x40 + 1.0/8.0;
    double x42 = -x41;
    double x43 = x10 + x40 - 1.0/8.0;
    double x44 = x26 + x39;
    double x45 = x0 + x44 - 1.0/8.0;
    double x46 = x10 + x44;
    double x47 = x46 + 1.0/8.0;
    double x48 = -x47;
    double x49 = -x45;
    double x50 = -x43;
    double x51 = x41*coord[20] + x42*coord[23] + x43*coord[11] + x45*coord[14] + x47*coord[5] + x48*coord[2] + x49*coord[17] + x50*coord[8];
    double x52 = x19*x51;
    double x53 = x13*coord[10] + x15*coord[13] + x16*coord[1] + x17*coord[22] + x18*coord[16] + x5*coord[19] + x6*coord[7] + x9*coord[4];
    double x54 = x41*coord[18] + x42*coord[21] + x43*coord[9] + x45*coord[12] + x47*coord[3] + x48*coord[0] + x49*coord[15] + x50*coord[6];
    double x55 = x33*x54;
    double x56 = x41*coord[19] + x42*coord[22] + x43*coord[10] + x45*coord[13] + x47*coord[4] + x48*coord[1] + x49*coord[16] + x50*coord[7];
    double x57 = x19*x33*x56 + x34*x38*x54 + x35*x51*x53 - x36*x56 - x38*x52 - x53*x55;
    double x58 = 1.0/x57;
    double x59 = x58*(x41*U[20] + x42*U[23] + x43*U[11] + x45*U[14] + x47*U[5] + x48*U[2] + x49*U[17] + x50*U[8]);
    double x60 = x34*x54 - x52;
    double x61 = x58*(x23*U[20] + x24*U[17] + x25*U[5] + x28*U[14] + x29*U[11] + x30*U[2] + x31*U[23] + x32*U[8]);
    double x62 = x35*x51 - x55;
    double x63 = x58*(x13*U[11] + x15*U[14] + x16*U[2] + x17*U[23] + x18*U[17] + x5*U[20] + x6*U[8] + x9*U[5]);
    double x64 = x37*x59 + x60*x61 + x62*x63;
    double x65 = -x34*x56 + x51*x53;
    double x66 = x58*(x23*U[19] + x24*U[16] + x25*U[4] + x28*U[13] + x29*U[10] + x30*U[1] + x31*U[22] + x32*U[7]);
    double x67 = -x33*x53 + x34*x38;
    double x68 = x58*(x41*U[19] + x42*U[22] + x43*U[10] + x45*U[13] + x47*U[4] + x48*U[1] + x49*U[16] + x50*U[7]);
    double x69 = x33*x56 - x38*x51;
    double x70 = x58*(x13*U[10] + x15*U[13] + x16*U[1] + x17*U[22] + x18*U[16] + x5*U[19] + x6*U[7] + x9*U[4]);
    double x71 = x65*x66 + x67*x68 + x69*x70;
    double x72 = x59*x67 + x61*x65 + x63*x69;
    double x73 = x37*x68 + x60*x66 + x62*x70;
    double x74 = x73 + 1.0;
    double x75 = x64*x71 - x72*x74;
    double x76 = x58*(x41*U[18] + x42*U[21] + x43*U[9] + x45*U[12] + x47*U[3] + x48*U[0] + x49*U[15] + x50*U[6]);
    double x77 = x58*(x23*U[18] + x24*U[15] + x25*U[3] + x28*U[12] + x29*U[9] + x30*U[0] + x31*U[21] + x32*U[6]);
    double x78 = x58*(x13*U[9] + x15*U[12] + x16*U[0] + x17*U[21] + x18*U[15] + x5*U[18] + x6*U[6] + x9*U[3]);
    double x79 = x37*x76 + x60*x77 + x62*x78;
    double x80 = x65*x77 + x67*x76 + x69*x78;
    double x81 = x80 + 1.0;
    double x82 = -x64*x81 + x72*x79;
    double x83 = -x71*x79 + x74*x81;
    double x84 = S11*x75 + S12*x82 + S13*x83;
    double x85 = S12*x75 + S22*x82 + S23*x83;
    double x86 = S13*x75 + S23*x82 + S33*x83;
    double x87 = -x71*x85 - x72*x86 - x80*x84;
    double x88 = x19*x56 - x53*x54;
    double x89 = x30*x58;
    double x90 = -x19*x38 + x35*x53;
    double x91 = x48*x58;
    double x92 = -x35*x56 + x38*x54;
    double x93 = x16*x58;
    double x94 = x88*x89 + x90*x91 + x92*x93;
    double x95 = x66*x88 + x68*x90 + x70*x92;
    double x96 = x59*x90 + x61*x88 + x63*x92;
    double x97 = x96 + 1.0;
    double x98 = -x71*x97 + x72*x95;
    double x99 = x76*x90 + x77*x88 + x78*x92;
    double x100 = -x72*x99 + x81*x97;
    double x101 = x71*x99 - x81*x95;
    double x102 = S11*x98 + S12*x100 + S13*x101;
    double x103 = S12*x98 + S22*x100 + S23*x101;
    double x104 = S13*x98 + S23*x100 + S33*x101;
    double x105 = -x102*x80 - x103*x71 - x104*x72;
    double x106 = x37*x91 + x60*x89 + x62*x93;
    double x107 = x65*x89 + x67*x91 + x69*x93;
    double x108 = -x64*x95 + x74*x97;
    double x109 = x64*x99 - x79*x97;
    double x110 = -x74*x99 + x79*x95;
    double x111 = S11*x108 + S12*x109 + S13*x110;
    double x112 = S12*x108 + S22*x109 + S23*x110;
    double x113 = S13*x108 + S23*x109 + S33*x110;
    double x114 = -x2;
    double x115 = -x21;
    double x116 = x114 + x115;
    double x117 = t*x2;
    double x118 = x117 - x39;
    double x119 = x118 + x20;
    double x120 = -x0 - x116 - x119 - x12;
    double x121 = x115 + x119 + x15;
    double x122 = x0 + x114 + x118 + x29;
    double x123 = -x118 - x13 - x27;
    double x124 = x116 + x117 + x4 + x46;
    double x125 = -x115 - x117 - x44 - x9;
    double x126 = x117 + x39;
    double x127 = -x10 - x114 - x126 - x25;
    double x128 = x126 + x22 + x5;
    double x129 = x120*V[0] + x121*V[12] + x122*V[9] + x123*V[21] + x124*V[3] + x125*V[15] + x127*V[6] + x128*V[18];
    double x130 = x120*V[1] + x121*V[13] + x122*V[10] + x123*V[22] + x124*V[4] + x125*V[16] + x127*V[7] + x128*V[19];
    double x131 = x120*V[2] + x121*V[14] + x122*V[11] + x123*V[23] + x124*V[5] + x125*V[17] + x127*V[8] + x128*V[20];
    double x132 = -1.0*PENER - 1.0*SENER + 0.5*rho*(x129*x129 + x130*x130 + x131*x131);
    double x133 = -x111*x80 - x112*x71 - x113*x72 - x132;
    double x134 = x120*A[0] + x121*A[12] + x122*A[9] + x123*A[21] + x124*A[3] + x125*A[15] + x127*A[6] + x128*A[18];
    double x135 = x120*A[1] + x121*A[13] + x122*A[10] + x123*A[22] + x124*A[4] + x125*A[16] + x127*A[7] + x128*A[19];
    double x136 = x120*A[2] + x121*A[14] + x122*A[11] + x123*A[23] + x124*A[5] + x125*A[17] + x127*A[8] + x128*A[20];
    double x137 = x23*V[18] + x24*V[15] + x25*V[3] + x28*V[12] + x29*V[9] + x30*V[0] + x31*V[21] + x32*V[6];
    double x138 = x58*x65;
    double x139 = x41*V[18] + x42*V[21] + x43*V[9] + x45*V[12] + x47*V[3] + x48*V[0] + x49*V[15] + x50*V[6];
    double x140 = x58*x67;
    double x141 = x13*V[9] + x15*V[12] + x16*V[0] + x17*V[21] + x18*V[15] + x5*V[18] + x6*V[6] + x9*V[3];
    double x142 = x58*x69;
    double x143 = x23*V[19] + x24*V[16] + x25*V[4] + x28*V[13] + x29*V[10] + x30*V[1] + x31*V[22] + x32*V[7];
    double x144 = x41*V[19] + x42*V[22] + x43*V[10] + x45*V[13] + x47*V[4] + x48*V[1] + x49*V[16] + x50*V[7];
    double x145 = x13*V[10] + x15*V[13] + x16*V[1] + x17*V[22] + x18*V[16] + x5*V[19] + x6*V[7] + x9*V[4];
    double x146 = x23*V[20] + x24*V[17] + x25*V[5] + x28*V[14] + x29*V[11] + x30*V[2] + x31*V[23] + x32*V[8];
    double x147 = x41*V[20] + x42*V[23] + x43*V[11] + x45*V[14] + x47*V[5] + x48*V[2] + x49*V[17] + x50*V[8];
    double x148 = x13*V[11] + x15*V[14] + x16*V[2] + x17*V[23] + x18*V[17] + x5*V[20] + x6*V[8] + x9*V[5];
    double x149 = x129*(x137*x138 + x139*x140 + x141*x142) + x130*(x138*x143 + x140*x144 + x142*x145) + x131*(x138*x146 + x140*x147 + x142*x148) + x134*x80 + x135*x71 + x136*x72;
    double x150 = rho*x120;
    double x151 = -x64*x86 - x73*x85 - x79*x84;
    double x152 = -x111*x79 - x112*x73 - x113*x64;
    double x153 = -x102*x79 - x103*x73 - x104*x64 - x132;
    double x154 = x37*x58;
    double x155 = x58*x60;
    double x156 = x58*x62;
    double x157 = x129*(x137*x155 + x139*x154 + x141*x156) + x130*(x143*x155 + x144*x154 + x145*x156) + x131*(x146*x155 + x147*x154 + x148*x156) + x134*x79 + x135*x73 + x136*x64;
    double x158 = -x102*x99 - x103*x95 - x104*x96;
    double x159 = -x111*x99 - x112*x95 - x113*x96;
    double x160 = -x132 - x84*x99 - x85*x95 - x86*x96;
    double x161 = x58*x88;
    double x162 = x58*x90;
    double x163 = x58*x92;
    double x164 = x129*(x137*x161 + x139*x162 + x141*x163) + x130*(x143*x161 + x144*x162 + x145*x163) + x131*(x146*x161 + x147*x162 + x148*x163) + x134*x99 + x135*x95 + x136*x96;
    double x165 = x161*x25 + x162*x47 + x163*x9;
    double x166 = x154*x47 + x155*x25 + x156*x9;
    double x167 = x138*x25 + x140*x47 + x142*x9;
    double x168 = rho*x124;
    double x169 = x161*x32 + x162*x50 + x163*x6;
    double x170 = x154*x50 + x155*x32 + x156*x6;
    double x171 = x138*x32 + x140*x50 + x142*x6;
    double x172 = rho*x127;
    double x173 = x13*x163 + x161*x29 + x162*x43;
    double x174 = x13*x156 + x154*x43 + x155*x29;
    double x175 = x13*x142 + x138*x29 + x140*x43;
    double x176 = rho*x122;
    double x177 = x15*x163 + x161*x28 + x162*x45;
    double x178 = x15*x156 + x154*x45 + x155*x28;
    double x179 = x138*x28 + x140*x45 + x142*x15;
    double x180 = rho*x121;
    double x181 = x161*x24 + x162*x49 + x163*x18;
    double x182 = x154*x49 + x155*x24 + x156*x18;
    double x183 = x138*x24 + x140*x49 + x142*x18;
    double x184 = rho*x125;
    double x185 = x161*x23 + x162*x41 + x163*x5;
    double x186 = x154*x41 + x155*x23 + x156*x5;
    double x187 = x138*x23 + x140*x41 + x142*x5;
    double x188 = rho*x128;
    double x189 = x161*x31 + x162*x42 + x163*x17;
    double x190 = x154*x42 + x155*x31 + x156*x17;
    double x191 = x138*x31 + x140*x42 + x142*x17;
    double x192 = rho*x123;
    
    res_0[0] = x57*(x105*x106 + x107*x133 - x149*x150 + x87*x94);
    res_0[1] = x57*(x106*x153 + x107*x152 - x150*x157 + x151*x94);
    res_0[2] = x57*(x106*x158 + x107*x159 - x150*x164 + x160*x94);
    res_0[3] = x57*(x105*x166 + x133*x167 - x149*x168 + x165*x87);
    res_0[4] = x57*(x151*x165 + x152*x167 + x153*x166 - x157*x168);
    res_0[5] = x57*(x158*x166 + x159*x167 + x160*x165 - x164*x168);
    res_0[6] = x57*(x105*x170 + x133*x171 - x149*x172 + x169*x87);
    res_0[7] = x57*(x151*x169 + x152*x171 + x153*x170 - x157*x172);
    res_0[8] = x57*(x158*x170 + x159*x171 + x160*x169 - x164*x172);
    res_0[9] = x57*(x105*x174 + x133*x175 - x149*x176 + x173*x87);
    res_0[10] = x57*(x151*x173 + x152*x175 + x153*x174 - x157*x176);
    res_0[11] = x57*(x158*x174 + x159*x175 + x160*x173 - x164*x176);
    res_0[12] = x57*(x105*x178 + x133*x179 - x149*x180 + x177*x87);
    res_0[13] = x57*(x151*x177 + x152*x179 + x153*x178 - x157*x180);
    res_0[14] = x57*(x158*x178 + x159*x179 + x160*x177 - x164*x180);
    res_0[15] = x57*(x105*x182 + x133*x183 - x149*x184 + x181*x87);
    res_0[16] = x57*(x151*x181 + x152*x183 + x153*x182 - x157*x184);
    res_0[17] = x57*(x158*x182 + x159*x183 + x160*x181 - x164*x184);
    res_0[18] = x57*(x105*x186 + x133*x187 - x149*x188 + x185*x87);
    res_0[19] = x57*(x151*x185 + x152*x187 + x153*x186 - x157*x188);
    res_0[20] = x57*(x158*x186 + x159*x187 + x160*x185 - x164*x188);
    res_0[21] = x57*(x105*x190 + x133*x191 - x149*x192 + x189*x87);
    res_0[22] = x57*(x151*x189 + x152*x191 + x153*x190 - x157*x192);
    res_0[23] = x57*(x158*x190 + x159*x191 + x160*x189 - x164*x192);
}

Conf_Forces_API void Integration_C3D8_dynamic(size_t num_elem,double Coords[][8][3],double *rho,double Element_U[][8][3],double Element_V[][8][3],double Element_A[][8][3],double S[][8][6],
    double PENER[][8],double SENER[][8],double Conf_Force[][8][3]){
    
    //Inputs:
    //    Coords.........Element_Nodal_Coordinates  [num_elem][numNodes,3]
    //    Element_U......Element_Nodal_Displacements[num_elem][numNodes,3]
    //    Element_V......Element_Nodal_Velocity     [num_elem][numNodes,3]
    //    Element_A......Element_Nodal_Acceleration [num_elem][numNodes,3]
    //    S..............Stress at integration points [num_elem][num_int_points,6]
    //    PENER..........Plastic strain energy [num_elem][num_int_points]
    //    SENER..........Elastic strain energy [num_elem][num_int_points]
    //    
    double weights[8]={1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,};
    double int_points[8][3]={
        {-0.5773502691896258,-0.5773502691896258,-0.5773502691896258,},
        {0.5773502691896258,-0.5773502691896258,-0.5773502691896258,},
        {-0.5773502691896258,0.5773502691896258,-0.5773502691896258,},
        {0.5773502691896258,0.5773502691896258,-0.5773502691896258,},
        {-0.5773502691896258,-0.5773502691896258,0.5773502691896258,},
        {0.5773502691896258,-0.5773502691896258,0.5773502691896258,},
        {-0.5773502691896258,0.5773502691896258,0.5773502691896258,},
        {0.5773502691896258,0.5773502691896258,0.5773502691896258,}};

    #pragma omp parallel for
    for (size_t i=0;i<num_elem;i++){
        double TMP[8][3];
        for (size_t j=0;j<8;j++){
            C3D8_dynamic((double*) int_points[j],(double*) Coords[i],rho[i],(double*) Element_U[i],(double*) Element_V[i],(double*) Element_A[i],(double*) S[i][j],PENER[i][j],SENER[i][j],(double*) TMP);
            for (size_t num_node=0;num_node<8;num_node++){
                for (size_t ii=0;ii<3;ii++){
                    Conf_Force[i][num_node][ii]+=TMP[num_node][ii]*weights[j];
                }
            }
        }
    }
}

void C3D8_static_mbf(double *rst,double *coord,double *Element_U,double *S,double PENER,double SENER,double *res_0){
    double r = rst[0];
    double s = rst[1];
    double t = rst[2];
    double *U = Element_U;
    double S11 = S[0];
    double S22 = S[1];
    double S33 = S[2];
    double S12 = S[3];
    double S13 = S[4];
    double S23 = S[5];
    
    double x0 = (1.0/8.0)*s;
    double x1 = (1.0/8.0)*r;
    double x2 = s*x1;
    double x3 = x0 + x2;
    double x4 = x1 + 1.0/8.0;
    double x5 = x3 + x4;
    double x6 = -x5;
    double x7 = -x1;
    double x8 = x7 - 1.0/8.0;
    double x9 = x3 + x8;
    double x10 = -x0;
    double x11 = x10 + x2;
    double x12 = x1 - 1.0/8.0;
    double x13 = x11 + x12;
    double x14 = x7 + 1.0/8.0;
    double x15 = x11 + x14;
    double x16 = -x15;
    double x17 = -x13;
    double x18 = -x9;
    double x19 = x13*coord[9] + x15*coord[12] + x16*coord[0] + x17*coord[21] + x18*coord[15] + x5*coord[18] + x6*coord[6] + x9*coord[3];
    double x20 = (1.0/8.0)*t;
    double x21 = t*x1;
    double x22 = x20 + x21;
    double x23 = x22 + x4;
    double x24 = -x23;
    double x25 = x22 + x8;
    double x26 = -x20;
    double x27 = x21 + x26;
    double x28 = x12 + x27;
    double x29 = x14 + x27;
    double x30 = -x29;
    double x31 = -x28;
    double x32 = -x25;
    double x33 = x23*coord[20] + x24*coord[17] + x25*coord[5] + x28*coord[14] + x29*coord[11] + x30*coord[2] + x31*coord[23] + x32*coord[8];
    double x34 = x13*coord[11] + x15*coord[14] + x16*coord[2] + x17*coord[23] + x18*coord[17] + x5*coord[20] + x6*coord[8] + x9*coord[5];
    double x35 = x23*coord[18] + x24*coord[15] + x25*coord[3] + x28*coord[12] + x29*coord[9] + x30*coord[0] + x31*coord[21] + x32*coord[6];
    double x36 = x34*x35;
    double x37 = x19*x33 - x36;
    double x38 = x23*coord[19] + x24*coord[16] + x25*coord[4] + x28*coord[13] + x29*coord[10] + x30*coord[1] + x31*coord[22] + x32*coord[7];
    double x39 = t*x0;
    double x40 = x20 + x39;
    double x41 = x0 + x40 + 1.0/8.0;
    double x42 = -x41;
    double x43 = x10 + x40 - 1.0/8.0;
    double x44 = x26 + x39;
    double x45 = x0 + x44 - 1.0/8.0;
    double x46 = x10 + x44 + 1.0/8.0;
    double x47 = -x46;
    double x48 = -x45;
    double x49 = -x43;
    double x50 = x41*coord[20] + x42*coord[23] + x43*coord[11] + x45*coord[14] + x46*coord[5] + x47*coord[2] + x48*coord[17] + x49*coord[8];
    double x51 = x19*x50;
    double x52 = x13*coord[10] + x15*coord[13] + x16*coord[1] + x17*coord[22] + x18*coord[16] + x5*coord[19] + x6*coord[7] + x9*coord[4];
    double x53 = x41*coord[18] + x42*coord[21] + x43*coord[9] + x45*coord[12] + x46*coord[3] + x47*coord[0] + x48*coord[15] + x49*coord[6];
    double x54 = x33*x53;
    double x55 = x41*coord[19] + x42*coord[22] + x43*coord[10] + x45*coord[13] + x46*coord[4] + x47*coord[1] + x48*coord[16] + x49*coord[7];
    double x56 = x19*x33*x55 + x34*x38*x53 + x35*x50*x52 - x36*x55 - x38*x51 - x52*x54;
    double x57 = 1.0/x56;
    double x58 = x57*(x41*U[20] + x42*U[23] + x43*U[11] + x45*U[14] + x46*U[5] + x47*U[2] + x48*U[17] + x49*U[8]);
    double x59 = x34*x53 - x51;
    double x60 = x57*(x23*U[20] + x24*U[17] + x25*U[5] + x28*U[14] + x29*U[11] + x30*U[2] + x31*U[23] + x32*U[8]);
    double x61 = x35*x50 - x54;
    double x62 = x57*(x13*U[11] + x15*U[14] + x16*U[2] + x17*U[23] + x18*U[17] + x5*U[20] + x6*U[8] + x9*U[5]);
    double x63 = x37*x58 + x59*x60 + x61*x62;
    double x64 = -x34*x55 + x50*x52;
    double x65 = x57*(x23*U[19] + x24*U[16] + x25*U[4] + x28*U[13] + x29*U[10] + x30*U[1] + x31*U[22] + x32*U[7]);
    double x66 = -x33*x52 + x34*x38;
    double x67 = x57*(x41*U[19] + x42*U[22] + x43*U[10] + x45*U[13] + x46*U[4] + x47*U[1] + x48*U[16] + x49*U[7]);
    double x68 = x33*x55 - x38*x50;
    double x69 = x57*(x13*U[10] + x15*U[13] + x16*U[1] + x17*U[22] + x18*U[16] + x5*U[19] + x6*U[7] + x9*U[4]);
    double x70 = x64*x65 + x66*x67 + x68*x69;
    double x71 = x58*x66 + x60*x64 + x62*x68;
    double x72 = x37*x67 + x59*x65 + x61*x69 + 1.0;
    double x73 = x63*x70 - x71*x72;
    double x74 = x57*(x41*U[18] + x42*U[21] + x43*U[9] + x45*U[12] + x46*U[3] + x47*U[0] + x48*U[15] + x49*U[6]);
    double x75 = x57*(x23*U[18] + x24*U[15] + x25*U[3] + x28*U[12] + x29*U[9] + x30*U[0] + x31*U[21] + x32*U[6]);
    double x76 = x57*(x13*U[9] + x15*U[12] + x16*U[0] + x17*U[21] + x18*U[15] + x5*U[18] + x6*U[6] + x9*U[3]);
    double x77 = x37*x74 + x59*x75 + x61*x76;
    double x78 = x64*x75 + x66*x74 + x68*x76 + 1.0;
    double x79 = -x63*x78 + x71*x77;
    double x80 = -x70*x77 + x72*x78;
    double x81 = S11*x73 + S12*x79 + S13*x80;
    double x82 = S12*x73 + S22*x79 + S23*x80;
    double x83 = S13*x73 + S23*x79 + S33*x80;
    double x84 = -x70*x82 - x71*x83 - x78*x81;
    double x85 = x19*x55 - x52*x53;
    double x86 = x30*x57;
    double x87 = -x19*x38 + x35*x52;
    double x88 = x47*x57;
    double x89 = -x35*x55 + x38*x53;
    double x90 = x16*x57;
    double x91 = x85*x86 + x87*x88 + x89*x90;
    double x92 = x65*x85 + x67*x87 + x69*x89;
    double x93 = x58*x87 + x60*x85 + x62*x89 + 1.0;
    double x94 = -x70*x93 + x71*x92;
    double x95 = x74*x87 + x75*x85 + x76*x89;
    double x96 = -x71*x95 + x78*x93;
    double x97 = x70*x95 - x78*x92;
    double x98 = S11*x94 + S12*x96 + S13*x97;
    double x99 = S12*x94 + S22*x96 + S23*x97;
    double x100 = S13*x94 + S23*x96 + S33*x97;
    double x101 = -x100*x71 - x70*x99 - x78*x98;
    double x102 = x37*x88 + x59*x86 + x61*x90;
    double x103 = x64*x86 + x66*x88 + x68*x90;
    double x104 = -x63*x92 + x72*x93;
    double x105 = x63*x95 - x77*x93;
    double x106 = -x72*x95 + x77*x92;
    double x107 = S11*x104 + S12*x105 + S13*x106;
    double x108 = S12*x104 + S22*x105 + S23*x106;
    double x109 = S13*x104 + S23*x105 + S33*x106;
    double x110 = -1.0*PENER - 1.0*SENER;
    double x111 = -x107*x78 - x108*x70 - x109*x71 - x110;
    double x112 = -x63*x83 - x72*x82 - x77*x81;
    double x113 = -x107*x77 - x108*x72 - x109*x63;
    double x114 = -x100*x63 - x110 - x72*x99 - x77*x98;
    double x115 = -x100*x93 - x92*x99 - x95*x98;
    double x116 = -x107*x95 - x108*x92 - x109*x93;
    double x117 = -x110 - x81*x95 - x82*x92 - x83*x93;
    double x118 = x25*x57;
    double x119 = x46*x57;
    double x120 = x57*x9;
    double x121 = x118*x85 + x119*x87 + x120*x89;
    double x122 = x118*x59 + x119*x37 + x120*x61;
    double x123 = x118*x64 + x119*x66 + x120*x68;
    double x124 = x32*x57;
    double x125 = x49*x57;
    double x126 = x57*x6;
    double x127 = x124*x85 + x125*x87 + x126*x89;
    double x128 = x124*x59 + x125*x37 + x126*x61;
    double x129 = x124*x64 + x125*x66 + x126*x68;
    double x130 = x29*x57;
    double x131 = x43*x57;
    double x132 = x13*x57;
    double x133 = x130*x85 + x131*x87 + x132*x89;
    double x134 = x130*x59 + x131*x37 + x132*x61;
    double x135 = x130*x64 + x131*x66 + x132*x68;
    double x136 = x28*x57;
    double x137 = x45*x57;
    double x138 = x15*x57;
    double x139 = x136*x85 + x137*x87 + x138*x89;
    double x140 = x136*x59 + x137*x37 + x138*x61;
    double x141 = x136*x64 + x137*x66 + x138*x68;
    double x142 = x24*x57;
    double x143 = x48*x57;
    double x144 = x18*x57;
    double x145 = x142*x85 + x143*x87 + x144*x89;
    double x146 = x142*x59 + x143*x37 + x144*x61;
    double x147 = x142*x64 + x143*x66 + x144*x68;
    double x148 = x23*x57;
    double x149 = x41*x57;
    double x150 = x5*x57;
    double x151 = x148*x85 + x149*x87 + x150*x89;
    double x152 = x148*x59 + x149*x37 + x150*x61;
    double x153 = x148*x64 + x149*x66 + x150*x68;
    double x154 = x31*x57;
    double x155 = x42*x57;
    double x156 = x17*x57;
    double x157 = x154*x85 + x155*x87 + x156*x89;
    double x158 = x154*x59 + x155*x37 + x156*x61;
    double x159 = x154*x64 + x155*x66 + x156*x68;
    
    res_0[0] = x56*(x101*x102 + x103*x111 + x84*x91);
    res_0[1] = x56*(x102*x114 + x103*x113 + x112*x91);
    res_0[2] = x56*(x102*x115 + x103*x116 + x117*x91);
    res_0[3] = x56*(x101*x122 + x111*x123 + x121*x84);
    res_0[4] = x56*(x112*x121 + x113*x123 + x114*x122);
    res_0[5] = x56*(x115*x122 + x116*x123 + x117*x121);
    res_0[6] = x56*(x101*x128 + x111*x129 + x127*x84);
    res_0[7] = x56*(x112*x127 + x113*x129 + x114*x128);
    res_0[8] = x56*(x115*x128 + x116*x129 + x117*x127);
    res_0[9] = x56*(x101*x134 + x111*x135 + x133*x84);
    res_0[10] = x56*(x112*x133 + x113*x135 + x114*x134);
    res_0[11] = x56*(x115*x134 + x116*x135 + x117*x133);
    res_0[12] = x56*(x101*x140 + x111*x141 + x139*x84);
    res_0[13] = x56*(x112*x139 + x113*x141 + x114*x140);
    res_0[14] = x56*(x115*x140 + x116*x141 + x117*x139);
    res_0[15] = x56*(x101*x146 + x111*x147 + x145*x84);
    res_0[16] = x56*(x112*x145 + x113*x147 + x114*x146);
    res_0[17] = x56*(x115*x146 + x116*x147 + x117*x145);
    res_0[18] = x56*(x101*x152 + x111*x153 + x151*x84);
    res_0[19] = x56*(x112*x151 + x113*x153 + x114*x152);
    res_0[20] = x56*(x115*x152 + x116*x153 + x117*x151);
    res_0[21] = x56*(x101*x158 + x111*x159 + x157*x84);
    res_0[22] = x56*(x112*x157 + x113*x159 + x114*x158);
    res_0[23] = x56*(x115*x158 + x116*x159 + x117*x157);
}

Conf_Forces_API void Integration_C3D8_static_mbf(size_t num_elem,double Coords[][8][3],double Element_U[][8][3],double S[][8][6],
    double PENER[][8],double SENER[][8],double Conf_Force[][8][3]){
    
    //Inputs:
    //    Coords.........Element_Nodal_Coordinates [num_elem][numNodes,3]
    //    Element_U......Element_Nodal_Displacements [num_elem][numNodes,3]
    //    S..............Stress at integration points [num_elem][num_int_points,6]
    //    PENER..........Plastic strain energy [num_elem][num_int_points]
    //    SENER..........Elastic strain energy [num_elem][num_int_points]
    //    
    double weights[8]={1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,};
    double int_points[8][3]={
        {-0.5773502691896258,-0.5773502691896258,-0.5773502691896258,},
        {0.5773502691896258,-0.5773502691896258,-0.5773502691896258,},
        {-0.5773502691896258,0.5773502691896258,-0.5773502691896258,},
        {0.5773502691896258,0.5773502691896258,-0.5773502691896258,},
        {-0.5773502691896258,-0.5773502691896258,0.5773502691896258,},
        {0.5773502691896258,-0.5773502691896258,0.5773502691896258,},
        {-0.5773502691896258,0.5773502691896258,0.5773502691896258,},
        {0.5773502691896258,0.5773502691896258,0.5773502691896258,}};

    #pragma omp parallel for
    for (size_t i=0;i<num_elem;i++){
        double TMP[8][3];
        for (size_t j=0;j<8;j++){
            C3D8_static_mbf((double*) int_points[j],(double*) Coords[i],(double*) Element_U[i],(double*) S[i][j],PENER[i][j],SENER[i][j],(double*) TMP);
            for (size_t num_node=0;num_node<8;num_node++){
                for (size_t ii=0;ii<3;ii++){
                    Conf_Force[i][num_node][ii]+=TMP[num_node][ii]*weights[j];
                }
            }
        }
    }
}

void C3D8_static_dbf(double *rst,double *coord,double *Element_U,double *S,double PENER,double SENER,double *res_0){
    double r = rst[0];
    double s = rst[1];
    double t = rst[2];
    double *U = Element_U;
    double S11 = S[0];
    double S22 = S[1];
    double S33 = S[2];
    double S12 = S[3];
    double S13 = S[4];
    double S23 = S[5];
    
    double x0 = (1.0/8.0)*s;
    double x1 = (1.0/8.0)*r;
    double x2 = s*x1;
    double x3 = x0 + x2;
    double x4 = x1 + 1.0/8.0;
    double x5 = x3 + x4;
    double x6 = -x5;
    double x7 = -x1;
    double x8 = x7 - 1.0/8.0;
    double x9 = x3 + x8;
    double x10 = -x0;
    double x11 = x10 + x2;
    double x12 = x1 - 1.0/8.0;
    double x13 = x11 + x12;
    double x14 = x7 + 1.0/8.0;
    double x15 = x11 + x14;
    double x16 = -x15;
    double x17 = -x13;
    double x18 = -x9;
    double x19 = x13*coord[9] + x15*coord[12] + x16*coord[0] + x17*coord[21] + x18*coord[15] + x5*coord[18] + x6*coord[6] + x9*coord[3];
    double x20 = (1.0/8.0)*t;
    double x21 = t*x1;
    double x22 = x20 + x21;
    double x23 = x22 + x4;
    double x24 = -x23;
    double x25 = x22 + x8;
    double x26 = -x20;
    double x27 = x21 + x26;
    double x28 = x12 + x27;
    double x29 = x14 + x27;
    double x30 = -x29;
    double x31 = -x28;
    double x32 = -x25;
    double x33 = x23*coord[20] + x24*coord[17] + x25*coord[5] + x28*coord[14] + x29*coord[11] + x30*coord[2] + x31*coord[23] + x32*coord[8];
    double x34 = x13*coord[11] + x15*coord[14] + x16*coord[2] + x17*coord[23] + x18*coord[17] + x5*coord[20] + x6*coord[8] + x9*coord[5];
    double x35 = x23*coord[18] + x24*coord[15] + x25*coord[3] + x28*coord[12] + x29*coord[9] + x30*coord[0] + x31*coord[21] + x32*coord[6];
    double x36 = x34*x35;
    double x37 = x19*x33 - x36;
    double x38 = x23*coord[19] + x24*coord[16] + x25*coord[4] + x28*coord[13] + x29*coord[10] + x30*coord[1] + x31*coord[22] + x32*coord[7];
    double x39 = t*x0;
    double x40 = x20 + x39;
    double x41 = x0 + x40 + 1.0/8.0;
    double x42 = -x41;
    double x43 = x10 + x40 - 1.0/8.0;
    double x44 = x26 + x39;
    double x45 = x0 + x44 - 1.0/8.0;
    double x46 = x10 + x44 + 1.0/8.0;
    double x47 = -x46;
    double x48 = -x45;
    double x49 = -x43;
    double x50 = x41*coord[20] + x42*coord[23] + x43*coord[11] + x45*coord[14] + x46*coord[5] + x47*coord[2] + x48*coord[17] + x49*coord[8];
    double x51 = x19*x50;
    double x52 = x13*coord[10] + x15*coord[13] + x16*coord[1] + x17*coord[22] + x18*coord[16] + x5*coord[19] + x6*coord[7] + x9*coord[4];
    double x53 = x41*coord[18] + x42*coord[21] + x43*coord[9] + x45*coord[12] + x46*coord[3] + x47*coord[0] + x48*coord[15] + x49*coord[6];
    double x54 = x33*x53;
    double x55 = x41*coord[19] + x42*coord[22] + x43*coord[10] + x45*coord[13] + x46*coord[4] + x47*coord[1] + x48*coord[16] + x49*coord[7];
    double x56 = x19*x33*x55 + x34*x38*x53 + x35*x50*x52 - x36*x55 - x38*x51 - x52*x54;
    double x57 = 1.0/x56;
    double x58 = x57*(x41*U[20] + x42*U[23] + x43*U[11] + x45*U[14] + x46*U[5] + x47*U[2] + x48*U[17] + x49*U[8]);
    double x59 = x34*x53 - x51;
    double x60 = x57*(x23*U[20] + x24*U[17] + x25*U[5] + x28*U[14] + x29*U[11] + x30*U[2] + x31*U[23] + x32*U[8]);
    double x61 = x35*x50 - x54;
    double x62 = x57*(x13*U[11] + x15*U[14] + x16*U[2] + x17*U[23] + x18*U[17] + x5*U[20] + x6*U[8] + x9*U[5]);
    double x63 = x37*x58 + x59*x60 + x61*x62;
    double x64 = -x34*x55 + x50*x52;
    double x65 = x57*(x23*U[19] + x24*U[16] + x25*U[4] + x28*U[13] + x29*U[10] + x30*U[1] + x31*U[22] + x32*U[7]);
    double x66 = -x33*x52 + x34*x38;
    double x67 = x57*(x41*U[19] + x42*U[22] + x43*U[10] + x45*U[13] + x46*U[4] + x47*U[1] + x48*U[16] + x49*U[7]);
    double x68 = x33*x55 - x38*x50;
    double x69 = x57*(x13*U[10] + x15*U[13] + x16*U[1] + x17*U[22] + x18*U[16] + x5*U[19] + x6*U[7] + x9*U[4]);
    double x70 = x64*x65 + x66*x67 + x68*x69;
    double x71 = x58*x66 + x60*x64 + x62*x68;
    double x72 = x37*x67 + x59*x65 + x61*x69;
    double x73 = x72 + 1.0;
    double x74 = x63*x70 - x71*x73;
    double x75 = x57*(x41*U[18] + x42*U[21] + x43*U[9] + x45*U[12] + x46*U[3] + x47*U[0] + x48*U[15] + x49*U[6]);
    double x76 = x57*(x23*U[18] + x24*U[15] + x25*U[3] + x28*U[12] + x29*U[9] + x30*U[0] + x31*U[21] + x32*U[6]);
    double x77 = x57*(x13*U[9] + x15*U[12] + x16*U[0] + x17*U[21] + x18*U[15] + x5*U[18] + x6*U[6] + x9*U[3]);
    double x78 = x37*x75 + x59*x76 + x61*x77;
    double x79 = x64*x76 + x66*x75 + x68*x77;
    double x80 = x79 + 1.0;
    double x81 = -x63*x80 + x71*x78;
    double x82 = -x70*x78 + x73*x80;
    double x83 = S11*x74 + S12*x81 + S13*x82;
    double x84 = S12*x74 + S22*x81 + S23*x82;
    double x85 = S13*x74 + S23*x81 + S33*x82;
    double x86 = -x70*x84 - x71*x85 - x79*x83;
    double x87 = x19*x55 - x52*x53;
    double x88 = x30*x57;
    double x89 = -x19*x38 + x35*x52;
    double x90 = x47*x57;
    double x91 = -x35*x55 + x38*x53;
    double x92 = x16*x57;
    double x93 = x87*x88 + x89*x90 + x91*x92;
    double x94 = x65*x87 + x67*x89 + x69*x91;
    double x95 = x58*x89 + x60*x87 + x62*x91;
    double x96 = x95 + 1.0;
    double x97 = -x70*x96 + x71*x94;
    double x98 = x75*x89 + x76*x87 + x77*x91;
    double x99 = -x71*x98 + x80*x96;
    double x100 = x70*x98 - x80*x94;
    double x101 = S11*x97 + S12*x99 + S13*x100;
    double x102 = S12*x97 + S22*x99 + S23*x100;
    double x103 = S13*x97 + S23*x99 + S33*x100;
    double x104 = -x101*x79 - x102*x70 - x103*x71;
    double x105 = x37*x90 + x59*x88 + x61*x92;
    double x106 = x64*x88 + x66*x90 + x68*x92;
    double x107 = -x63*x94 + x73*x96;
    double x108 = x63*x98 - x78*x96;
    double x109 = -x73*x98 + x78*x94;
    double x110 = S11*x107 + S12*x108 + S13*x109;
    double x111 = S12*x107 + S22*x108 + S23*x109;
    double x112 = S13*x107 + S23*x108 + S33*x109;
    double x113 = -1.0*PENER - 1.0*SENER;
    double x114 = -x110*x79 - x111*x70 - x112*x71 - x113;
    double x115 = -x63*x85 - x72*x84 - x78*x83;
    double x116 = -x110*x78 - x111*x72 - x112*x63;
    double x117 = -x101*x78 - x102*x72 - x103*x63 - x113;
    double x118 = -x101*x98 - x102*x94 - x103*x95;
    double x119 = -x110*x98 - x111*x94 - x112*x95;
    double x120 = -x113 - x83*x98 - x84*x94 - x85*x95;
    double x121 = x25*x57;
    double x122 = x46*x57;
    double x123 = x57*x9;
    double x124 = x121*x87 + x122*x89 + x123*x91;
    double x125 = x121*x59 + x122*x37 + x123*x61;
    double x126 = x121*x64 + x122*x66 + x123*x68;
    double x127 = x32*x57;
    double x128 = x49*x57;
    double x129 = x57*x6;
    double x130 = x127*x87 + x128*x89 + x129*x91;
    double x131 = x127*x59 + x128*x37 + x129*x61;
    double x132 = x127*x64 + x128*x66 + x129*x68;
    double x133 = x29*x57;
    double x134 = x43*x57;
    double x135 = x13*x57;
    double x136 = x133*x87 + x134*x89 + x135*x91;
    double x137 = x133*x59 + x134*x37 + x135*x61;
    double x138 = x133*x64 + x134*x66 + x135*x68;
    double x139 = x28*x57;
    double x140 = x45*x57;
    double x141 = x15*x57;
    double x142 = x139*x87 + x140*x89 + x141*x91;
    double x143 = x139*x59 + x140*x37 + x141*x61;
    double x144 = x139*x64 + x140*x66 + x141*x68;
    double x145 = x24*x57;
    double x146 = x48*x57;
    double x147 = x18*x57;
    double x148 = x145*x87 + x146*x89 + x147*x91;
    double x149 = x145*x59 + x146*x37 + x147*x61;
    double x150 = x145*x64 + x146*x66 + x147*x68;
    double x151 = x23*x57;
    double x152 = x41*x57;
    double x153 = x5*x57;
    double x154 = x151*x87 + x152*x89 + x153*x91;
    double x155 = x151*x59 + x152*x37 + x153*x61;
    double x156 = x151*x64 + x152*x66 + x153*x68;
    double x157 = x31*x57;
    double x158 = x42*x57;
    double x159 = x17*x57;
    double x160 = x157*x87 + x158*x89 + x159*x91;
    double x161 = x157*x59 + x158*x37 + x159*x61;
    double x162 = x157*x64 + x158*x66 + x159*x68;
    
    res_0[0] = x56*(x104*x105 + x106*x114 + x86*x93);
    res_0[1] = x56*(x105*x117 + x106*x116 + x115*x93);
    res_0[2] = x56*(x105*x118 + x106*x119 + x120*x93);
    res_0[3] = x56*(x104*x125 + x114*x126 + x124*x86);
    res_0[4] = x56*(x115*x124 + x116*x126 + x117*x125);
    res_0[5] = x56*(x118*x125 + x119*x126 + x120*x124);
    res_0[6] = x56*(x104*x131 + x114*x132 + x130*x86);
    res_0[7] = x56*(x115*x130 + x116*x132 + x117*x131);
    res_0[8] = x56*(x118*x131 + x119*x132 + x120*x130);
    res_0[9] = x56*(x104*x137 + x114*x138 + x136*x86);
    res_0[10] = x56*(x115*x136 + x116*x138 + x117*x137);
    res_0[11] = x56*(x118*x137 + x119*x138 + x120*x136);
    res_0[12] = x56*(x104*x143 + x114*x144 + x142*x86);
    res_0[13] = x56*(x115*x142 + x116*x144 + x117*x143);
    res_0[14] = x56*(x118*x143 + x119*x144 + x120*x142);
    res_0[15] = x56*(x104*x149 + x114*x150 + x148*x86);
    res_0[16] = x56*(x115*x148 + x116*x150 + x117*x149);
    res_0[17] = x56*(x118*x149 + x119*x150 + x120*x148);
    res_0[18] = x56*(x104*x155 + x114*x156 + x154*x86);
    res_0[19] = x56*(x115*x154 + x116*x156 + x117*x155);
    res_0[20] = x56*(x118*x155 + x119*x156 + x120*x154);
    res_0[21] = x56*(x104*x161 + x114*x162 + x160*x86);
    res_0[22] = x56*(x115*x160 + x116*x162 + x117*x161);
    res_0[23] = x56*(x118*x161 + x119*x162 + x120*x160);
}

Conf_Forces_API void Integration_C3D8_static_dbf(size_t num_elem,double Coords[][8][3],double Element_U[][8][3],double S[][8][6],
    double PENER[][8],double SENER[][8],double Conf_Force[][8][3]){
    
    //Inputs:
    //    Coords.........Element_Nodal_Coordinates [num_elem][numNodes,3]
    //    Element_U......Element_Nodal_Displacements [num_elem][numNodes,3]
    //    S..............Stress at integration points [num_elem][num_int_points,6]
    //    PENER..........Plastic strain energy [num_elem][num_int_points]
    //    SENER..........Elastic strain energy [num_elem][num_int_points]
    //    
    double weights[8]={1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,};
    double int_points[8][3]={
        {-0.5773502691896258,-0.5773502691896258,-0.5773502691896258,},
        {0.5773502691896258,-0.5773502691896258,-0.5773502691896258,},
        {-0.5773502691896258,0.5773502691896258,-0.5773502691896258,},
        {0.5773502691896258,0.5773502691896258,-0.5773502691896258,},
        {-0.5773502691896258,-0.5773502691896258,0.5773502691896258,},
        {0.5773502691896258,-0.5773502691896258,0.5773502691896258,},
        {-0.5773502691896258,0.5773502691896258,0.5773502691896258,},
        {0.5773502691896258,0.5773502691896258,0.5773502691896258,}};

    #pragma omp parallel for
    for (size_t i=0;i<num_elem;i++){
        double TMP[8][3];
        for (size_t j=0;j<8;j++){
            C3D8_static_dbf((double*) int_points[j],(double*) Coords[i],(double*) Element_U[i],(double*) S[i][j],PENER[i][j],SENER[i][j],(double*) TMP);
            for (size_t num_node=0;num_node<8;num_node++){
                for (size_t ii=0;ii<3;ii++){
                    Conf_Force[i][num_node][ii]+=TMP[num_node][ii]*weights[j];
                }
            }
        }
    }
}

void C3D8R_dynamic(double *rst,double *coord,double rho,double *Element_U,double *Element_V,double *Element_A,double *S,double PENER,double SENER,double *res_0){
    double r = rst[0];
    double s = rst[1];
    double t = rst[2];
    double *U = Element_U;
    double *V = Element_V;
    double *A = Element_A;
    double S11 = S[0];
    double S22 = S[1];
    double S33 = S[2];
    double S12 = S[3];
    double S13 = S[4];
    double S23 = S[5];
    
    double x0 = (1.0/8.0)*s;
    double x1 = (1.0/8.0)*r;
    double x2 = s*x1;
    double x3 = x0 + x2;
    double x4 = x1 + 1.0/8.0;
    double x5 = x3 + x4;
    double x6 = -x5;
    double x7 = -x1;
    double x8 = x7 - 1.0/8.0;
    double x9 = x3 + x8;
    double x10 = -x0;
    double x11 = x10 + x2;
    double x12 = x1 - 1.0/8.0;
    double x13 = x11 + x12;
    double x14 = x7 + 1.0/8.0;
    double x15 = x11 + x14;
    double x16 = -x15;
    double x17 = -x13;
    double x18 = -x9;
    double x19 = x13*coord[9] + x15*coord[12] + x16*coord[0] + x17*coord[21] + x18*coord[15] + x5*coord[18] + x6*coord[6] + x9*coord[3];
    double x20 = (1.0/8.0)*t;
    double x21 = t*x1;
    double x22 = x20 + x21;
    double x23 = x22 + x4;
    double x24 = -x23;
    double x25 = x22 + x8;
    double x26 = -x20;
    double x27 = x21 + x26;
    double x28 = x12 + x27;
    double x29 = x14 + x27;
    double x30 = -x29;
    double x31 = -x28;
    double x32 = -x25;
    double x33 = x23*coord[20] + x24*coord[17] + x25*coord[5] + x28*coord[14] + x29*coord[11] + x30*coord[2] + x31*coord[23] + x32*coord[8];
    double x34 = x13*coord[11] + x15*coord[14] + x16*coord[2] + x17*coord[23] + x18*coord[17] + x5*coord[20] + x6*coord[8] + x9*coord[5];
    double x35 = x23*coord[18] + x24*coord[15] + x25*coord[3] + x28*coord[12] + x29*coord[9] + x30*coord[0] + x31*coord[21] + x32*coord[6];
    double x36 = x34*x35;
    double x37 = x19*x33 - x36;
    double x38 = x23*coord[19] + x24*coord[16] + x25*coord[4] + x28*coord[13] + x29*coord[10] + x30*coord[1] + x31*coord[22] + x32*coord[7];
    double x39 = t*x0;
    double x40 = x20 + x39;
    double x41 = x0 + x40 + 1.0/8.0;
    double x42 = -x41;
    double x43 = x10 + x40 - 1.0/8.0;
    double x44 = x26 + x39;
    double x45 = x0 + x44 - 1.0/8.0;
    double x46 = x10 + x44;
    double x47 = x46 + 1.0/8.0;
    double x48 = -x47;
    double x49 = -x45;
    double x50 = -x43;
    double x51 = x41*coord[20] + x42*coord[23] + x43*coord[11] + x45*coord[14] + x47*coord[5] + x48*coord[2] + x49*coord[17] + x50*coord[8];
    double x52 = x19*x51;
    double x53 = x13*coord[10] + x15*coord[13] + x16*coord[1] + x17*coord[22] + x18*coord[16] + x5*coord[19] + x6*coord[7] + x9*coord[4];
    double x54 = x41*coord[18] + x42*coord[21] + x43*coord[9] + x45*coord[12] + x47*coord[3] + x48*coord[0] + x49*coord[15] + x50*coord[6];
    double x55 = x33*x54;
    double x56 = x41*coord[19] + x42*coord[22] + x43*coord[10] + x45*coord[13] + x47*coord[4] + x48*coord[1] + x49*coord[16] + x50*coord[7];
    double x57 = x19*x33*x56 + x34*x38*x54 + x35*x51*x53 - x36*x56 - x38*x52 - x53*x55;
    double x58 = 1.0/x57;
    double x59 = x58*(x41*U[20] + x42*U[23] + x43*U[11] + x45*U[14] + x47*U[5] + x48*U[2] + x49*U[17] + x50*U[8]);
    double x60 = x34*x54 - x52;
    double x61 = x58*(x23*U[20] + x24*U[17] + x25*U[5] + x28*U[14] + x29*U[11] + x30*U[2] + x31*U[23] + x32*U[8]);
    double x62 = x35*x51 - x55;
    double x63 = x58*(x13*U[11] + x15*U[14] + x16*U[2] + x17*U[23] + x18*U[17] + x5*U[20] + x6*U[8] + x9*U[5]);
    double x64 = x37*x59 + x60*x61 + x62*x63;
    double x65 = -x34*x56 + x51*x53;
    double x66 = x58*(x23*U[19] + x24*U[16] + x25*U[4] + x28*U[13] + x29*U[10] + x30*U[1] + x31*U[22] + x32*U[7]);
    double x67 = -x33*x53 + x34*x38;
    double x68 = x58*(x41*U[19] + x42*U[22] + x43*U[10] + x45*U[13] + x47*U[4] + x48*U[1] + x49*U[16] + x50*U[7]);
    double x69 = x33*x56 - x38*x51;
    double x70 = x58*(x13*U[10] + x15*U[13] + x16*U[1] + x17*U[22] + x18*U[16] + x5*U[19] + x6*U[7] + x9*U[4]);
    double x71 = x65*x66 + x67*x68 + x69*x70;
    double x72 = x59*x67 + x61*x65 + x63*x69;
    double x73 = x37*x68 + x60*x66 + x62*x70;
    double x74 = x73 + 1.0;
    double x75 = x64*x71 - x72*x74;
    double x76 = x58*(x41*U[18] + x42*U[21] + x43*U[9] + x45*U[12] + x47*U[3] + x48*U[0] + x49*U[15] + x50*U[6]);
    double x77 = x58*(x23*U[18] + x24*U[15] + x25*U[3] + x28*U[12] + x29*U[9] + x30*U[0] + x31*U[21] + x32*U[6]);
    double x78 = x58*(x13*U[9] + x15*U[12] + x16*U[0] + x17*U[21] + x18*U[15] + x5*U[18] + x6*U[6] + x9*U[3]);
    double x79 = x37*x76 + x60*x77 + x62*x78;
    double x80 = x65*x77 + x67*x76 + x69*x78;
    double x81 = x80 + 1.0;
    double x82 = -x64*x81 + x72*x79;
    double x83 = -x71*x79 + x74*x81;
    double x84 = S11*x75 + S12*x82 + S13*x83;
    double x85 = S12*x75 + S22*x82 + S23*x83;
    double x86 = S13*x75 + S23*x82 + S33*x83;
    double x87 = -x71*x85 - x72*x86 - x80*x84;
    double x88 = x19*x56 - x53*x54;
    double x89 = x30*x58;
    double x90 = -x19*x38 + x35*x53;
    double x91 = x48*x58;
    double x92 = -x35*x56 + x38*x54;
    double x93 = x16*x58;
    double x94 = x88*x89 + x90*x91 + x92*x93;
    double x95 = x66*x88 + x68*x90 + x70*x92;
    double x96 = x59*x90 + x61*x88 + x63*x92;
    double x97 = x96 + 1.0;
    double x98 = -x71*x97 + x72*x95;
    double x99 = x76*x90 + x77*x88 + x78*x92;
    double x100 = -x72*x99 + x81*x97;
    double x101 = x71*x99 - x81*x95;
    double x102 = S11*x98 + S12*x100 + S13*x101;
    double x103 = S12*x98 + S22*x100 + S23*x101;
    double x104 = S13*x98 + S23*x100 + S33*x101;
    double x105 = -x102*x80 - x103*x71 - x104*x72;
    double x106 = x37*x91 + x60*x89 + x62*x93;
    double x107 = x65*x89 + x67*x91 + x69*x93;
    double x108 = -x64*x95 + x74*x97;
    double x109 = x64*x99 - x79*x97;
    double x110 = -x74*x99 + x79*x95;
    double x111 = S11*x108 + S12*x109 + S13*x110;
    double x112 = S12*x108 + S22*x109 + S23*x110;
    double x113 = S13*x108 + S23*x109 + S33*x110;
    double x114 = -x2;
    double x115 = -x21;
    double x116 = x114 + x115;
    double x117 = t*x2;
    double x118 = x117 - x39;
    double x119 = x118 + x20;
    double x120 = -x0 - x116 - x119 - x12;
    double x121 = x115 + x119 + x15;
    double x122 = x0 + x114 + x118 + x29;
    double x123 = -x118 - x13 - x27;
    double x124 = x116 + x117 + x4 + x46;
    double x125 = -x115 - x117 - x44 - x9;
    double x126 = x117 + x39;
    double x127 = -x10 - x114 - x126 - x25;
    double x128 = x126 + x22 + x5;
    double x129 = x120*V[0] + x121*V[12] + x122*V[9] + x123*V[21] + x124*V[3] + x125*V[15] + x127*V[6] + x128*V[18];
    double x130 = x120*V[1] + x121*V[13] + x122*V[10] + x123*V[22] + x124*V[4] + x125*V[16] + x127*V[7] + x128*V[19];
    double x131 = x120*V[2] + x121*V[14] + x122*V[11] + x123*V[23] + x124*V[5] + x125*V[17] + x127*V[8] + x128*V[20];
    double x132 = -1.0*PENER - 1.0*SENER + 0.5*rho*(x129*x129 + x130*x130 + x131*x131);
    double x133 = -x111*x80 - x112*x71 - x113*x72 - x132;
    double x134 = x120*A[0] + x121*A[12] + x122*A[9] + x123*A[21] + x124*A[3] + x125*A[15] + x127*A[6] + x128*A[18];
    double x135 = x120*A[1] + x121*A[13] + x122*A[10] + x123*A[22] + x124*A[4] + x125*A[16] + x127*A[7] + x128*A[19];
    double x136 = x120*A[2] + x121*A[14] + x122*A[11] + x123*A[23] + x124*A[5] + x125*A[17] + x127*A[8] + x128*A[20];
    double x137 = x23*V[18] + x24*V[15] + x25*V[3] + x28*V[12] + x29*V[9] + x30*V[0] + x31*V[21] + x32*V[6];
    double x138 = x58*x65;
    double x139 = x41*V[18] + x42*V[21] + x43*V[9] + x45*V[12] + x47*V[3] + x48*V[0] + x49*V[15] + x50*V[6];
    double x140 = x58*x67;
    double x141 = x13*V[9] + x15*V[12] + x16*V[0] + x17*V[21] + x18*V[15] + x5*V[18] + x6*V[6] + x9*V[3];
    double x142 = x58*x69;
    double x143 = x23*V[19] + x24*V[16] + x25*V[4] + x28*V[13] + x29*V[10] + x30*V[1] + x31*V[22] + x32*V[7];
    double x144 = x41*V[19] + x42*V[22] + x43*V[10] + x45*V[13] + x47*V[4] + x48*V[1] + x49*V[16] + x50*V[7];
    double x145 = x13*V[10] + x15*V[13] + x16*V[1] + x17*V[22] + x18*V[16] + x5*V[19] + x6*V[7] + x9*V[4];
    double x146 = x23*V[20] + x24*V[17] + x25*V[5] + x28*V[14] + x29*V[11] + x30*V[2] + x31*V[23] + x32*V[8];
    double x147 = x41*V[20] + x42*V[23] + x43*V[11] + x45*V[14] + x47*V[5] + x48*V[2] + x49*V[17] + x50*V[8];
    double x148 = x13*V[11] + x15*V[14] + x16*V[2] + x17*V[23] + x18*V[17] + x5*V[20] + x6*V[8] + x9*V[5];
    double x149 = x129*(x137*x138 + x139*x140 + x141*x142) + x130*(x138*x143 + x140*x144 + x142*x145) + x131*(x138*x146 + x140*x147 + x142*x148) + x134*x80 + x135*x71 + x136*x72;
    double x150 = rho*x120;
    double x151 = -x64*x86 - x73*x85 - x79*x84;
    double x152 = -x111*x79 - x112*x73 - x113*x64;
    double x153 = -x102*x79 - x103*x73 - x104*x64 - x132;
    double x154 = x37*x58;
    double x155 = x58*x60;
    double x156 = x58*x62;
    double x157 = x129*(x137*x155 + x139*x154 + x141*x156) + x130*(x143*x155 + x144*x154 + x145*x156) + x131*(x146*x155 + x147*x154 + x148*x156) + x134*x79 + x135*x73 + x136*x64;
    double x158 = -x102*x99 - x103*x95 - x104*x96;
    double x159 = -x111*x99 - x112*x95 - x113*x96;
    double x160 = -x132 - x84*x99 - x85*x95 - x86*x96;
    double x161 = x58*x88;
    double x162 = x58*x90;
    double x163 = x58*x92;
    double x164 = x129*(x137*x161 + x139*x162 + x141*x163) + x130*(x143*x161 + x144*x162 + x145*x163) + x131*(x146*x161 + x147*x162 + x148*x163) + x134*x99 + x135*x95 + x136*x96;
    double x165 = x161*x25 + x162*x47 + x163*x9;
    double x166 = x154*x47 + x155*x25 + x156*x9;
    double x167 = x138*x25 + x140*x47 + x142*x9;
    double x168 = rho*x124;
    double x169 = x161*x32 + x162*x50 + x163*x6;
    double x170 = x154*x50 + x155*x32 + x156*x6;
    double x171 = x138*x32 + x140*x50 + x142*x6;
    double x172 = rho*x127;
    double x173 = x13*x163 + x161*x29 + x162*x43;
    double x174 = x13*x156 + x154*x43 + x155*x29;
    double x175 = x13*x142 + x138*x29 + x140*x43;
    double x176 = rho*x122;
    double x177 = x15*x163 + x161*x28 + x162*x45;
    double x178 = x15*x156 + x154*x45 + x155*x28;
    double x179 = x138*x28 + x140*x45 + x142*x15;
    double x180 = rho*x121;
    double x181 = x161*x24 + x162*x49 + x163*x18;
    double x182 = x154*x49 + x155*x24 + x156*x18;
    double x183 = x138*x24 + x140*x49 + x142*x18;
    double x184 = rho*x125;
    double x185 = x161*x23 + x162*x41 + x163*x5;
    double x186 = x154*x41 + x155*x23 + x156*x5;
    double x187 = x138*x23 + x140*x41 + x142*x5;
    double x188 = rho*x128;
    double x189 = x161*x31 + x162*x42 + x163*x17;
    double x190 = x154*x42 + x155*x31 + x156*x17;
    double x191 = x138*x31 + x140*x42 + x142*x17;
    double x192 = rho*x123;
    
    res_0[0] = x57*(x105*x106 + x107*x133 - x149*x150 + x87*x94);
    res_0[1] = x57*(x106*x153 + x107*x152 - x150*x157 + x151*x94);
    res_0[2] = x57*(x106*x158 + x107*x159 - x150*x164 + x160*x94);
    res_0[3] = x57*(x105*x166 + x133*x167 - x149*x168 + x165*x87);
    res_0[4] = x57*(x151*x165 + x152*x167 + x153*x166 - x157*x168);
    res_0[5] = x57*(x158*x166 + x159*x167 + x160*x165 - x164*x168);
    res_0[6] = x57*(x105*x170 + x133*x171 - x149*x172 + x169*x87);
    res_0[7] = x57*(x151*x169 + x152*x171 + x153*x170 - x157*x172);
    res_0[8] = x57*(x158*x170 + x159*x171 + x160*x169 - x164*x172);
    res_0[9] = x57*(x105*x174 + x133*x175 - x149*x176 + x173*x87);
    res_0[10] = x57*(x151*x173 + x152*x175 + x153*x174 - x157*x176);
    res_0[11] = x57*(x158*x174 + x159*x175 + x160*x173 - x164*x176);
    res_0[12] = x57*(x105*x178 + x133*x179 - x149*x180 + x177*x87);
    res_0[13] = x57*(x151*x177 + x152*x179 + x153*x178 - x157*x180);
    res_0[14] = x57*(x158*x178 + x159*x179 + x160*x177 - x164*x180);
    res_0[15] = x57*(x105*x182 + x133*x183 - x149*x184 + x181*x87);
    res_0[16] = x57*(x151*x181 + x152*x183 + x153*x182 - x157*x184);
    res_0[17] = x57*(x158*x182 + x159*x183 + x160*x181 - x164*x184);
    res_0[18] = x57*(x105*x186 + x133*x187 - x149*x188 + x185*x87);
    res_0[19] = x57*(x151*x185 + x152*x187 + x153*x186 - x157*x188);
    res_0[20] = x57*(x158*x186 + x159*x187 + x160*x185 - x164*x188);
    res_0[21] = x57*(x105*x190 + x133*x191 - x149*x192 + x189*x87);
    res_0[22] = x57*(x151*x189 + x152*x191 + x153*x190 - x157*x192);
    res_0[23] = x57*(x158*x190 + x159*x191 + x160*x189 - x164*x192);
}

Conf_Forces_API void Integration_C3D8R_dynamic(size_t num_elem,double Coords[][8][3],double *rho,double Element_U[][8][3],double Element_V[][8][3],double Element_A[][8][3],double S[][1][6],
    double PENER[][1],double SENER[][1],double Conf_Force[][8][3]){
    
    //Inputs:
    //    Coords.........Element_Nodal_Coordinates  [num_elem][numNodes,3]
    //    Element_U......Element_Nodal_Displacements[num_elem][numNodes,3]
    //    Element_V......Element_Nodal_Velocity     [num_elem][numNodes,3]
    //    Element_A......Element_Nodal_Acceleration [num_elem][numNodes,3]
    //    S..............Stress at integration points [num_elem][num_int_points,6]
    //    PENER..........Plastic strain energy [num_elem][num_int_points]
    //    SENER..........Elastic strain energy [num_elem][num_int_points]
    //    
    double weights[1]={8.0,};
    double int_points[1][3]={
        {0.0,0.0,0.0,}};

    #pragma omp parallel for
    for (size_t i=0;i<num_elem;i++){
        double TMP[8][3];
        for (size_t j=0;j<1;j++){
            C3D8R_dynamic((double*) int_points[j],(double*) Coords[i],rho[i],(double*) Element_U[i],(double*) Element_V[i],(double*) Element_A[i],(double*) S[i][j],PENER[i][j],SENER[i][j],(double*) TMP);
            for (size_t num_node=0;num_node<8;num_node++){
                for (size_t ii=0;ii<3;ii++){
                    Conf_Force[i][num_node][ii]+=TMP[num_node][ii]*weights[j];
                }
            }
        }
    }
}

void C3D8R_static_mbf(double *rst,double *coord,double *Element_U,double *S,double PENER,double SENER,double *res_0){
    double r = rst[0];
    double s = rst[1];
    double t = rst[2];
    double *U = Element_U;
    double S11 = S[0];
    double S22 = S[1];
    double S33 = S[2];
    double S12 = S[3];
    double S13 = S[4];
    double S23 = S[5];
    
    double x0 = (1.0/8.0)*s;
    double x1 = (1.0/8.0)*r;
    double x2 = s*x1;
    double x3 = x0 + x2;
    double x4 = x1 + 1.0/8.0;
    double x5 = x3 + x4;
    double x6 = -x5;
    double x7 = -x1;
    double x8 = x7 - 1.0/8.0;
    double x9 = x3 + x8;
    double x10 = -x0;
    double x11 = x10 + x2;
    double x12 = x1 - 1.0/8.0;
    double x13 = x11 + x12;
    double x14 = x7 + 1.0/8.0;
    double x15 = x11 + x14;
    double x16 = -x15;
    double x17 = -x13;
    double x18 = -x9;
    double x19 = x13*coord[9] + x15*coord[12] + x16*coord[0] + x17*coord[21] + x18*coord[15] + x5*coord[18] + x6*coord[6] + x9*coord[3];
    double x20 = (1.0/8.0)*t;
    double x21 = t*x1;
    double x22 = x20 + x21;
    double x23 = x22 + x4;
    double x24 = -x23;
    double x25 = x22 + x8;
    double x26 = -x20;
    double x27 = x21 + x26;
    double x28 = x12 + x27;
    double x29 = x14 + x27;
    double x30 = -x29;
    double x31 = -x28;
    double x32 = -x25;
    double x33 = x23*coord[20] + x24*coord[17] + x25*coord[5] + x28*coord[14] + x29*coord[11] + x30*coord[2] + x31*coord[23] + x32*coord[8];
    double x34 = x13*coord[11] + x15*coord[14] + x16*coord[2] + x17*coord[23] + x18*coord[17] + x5*coord[20] + x6*coord[8] + x9*coord[5];
    double x35 = x23*coord[18] + x24*coord[15] + x25*coord[3] + x28*coord[12] + x29*coord[9] + x30*coord[0] + x31*coord[21] + x32*coord[6];
    double x36 = x34*x35;
    double x37 = x19*x33 - x36;
    double x38 = x23*coord[19] + x24*coord[16] + x25*coord[4] + x28*coord[13] + x29*coord[10] + x30*coord[1] + x31*coord[22] + x32*coord[7];
    double x39 = t*x0;
    double x40 = x20 + x39;
    double x41 = x0 + x40 + 1.0/8.0;
    double x42 = -x41;
    double x43 = x10 + x40 - 1.0/8.0;
    double x44 = x26 + x39;
    double x45 = x0 + x44 - 1.0/8.0;
    double x46 = x10 + x44 + 1.0/8.0;
    double x47 = -x46;
    double x48 = -x45;
    double x49 = -x43;
    double x50 = x41*coord[20] + x42*coord[23] + x43*coord[11] + x45*coord[14] + x46*coord[5] + x47*coord[2] + x48*coord[17] + x49*coord[8];
    double x51 = x19*x50;
    double x52 = x13*coord[10] + x15*coord[13] + x16*coord[1] + x17*coord[22] + x18*coord[16] + x5*coord[19] + x6*coord[7] + x9*coord[4];
    double x53 = x41*coord[18] + x42*coord[21] + x43*coord[9] + x45*coord[12] + x46*coord[3] + x47*coord[0] + x48*coord[15] + x49*coord[6];
    double x54 = x33*x53;
    double x55 = x41*coord[19] + x42*coord[22] + x43*coord[10] + x45*coord[13] + x46*coord[4] + x47*coord[1] + x48*coord[16] + x49*coord[7];
    double x56 = x19*x33*x55 + x34*x38*x53 + x35*x50*x52 - x36*x55 - x38*x51 - x52*x54;
    double x57 = 1.0/x56;
    double x58 = x57*(x41*U[20] + x42*U[23] + x43*U[11] + x45*U[14] + x46*U[5] + x47*U[2] + x48*U[17] + x49*U[8]);
    double x59 = x34*x53 - x51;
    double x60 = x57*(x23*U[20] + x24*U[17] + x25*U[5] + x28*U[14] + x29*U[11] + x30*U[2] + x31*U[23] + x32*U[8]);
    double x61 = x35*x50 - x54;
    double x62 = x57*(x13*U[11] + x15*U[14] + x16*U[2] + x17*U[23] + x18*U[17] + x5*U[20] + x6*U[8] + x9*U[5]);
    double x63 = x37*x58 + x59*x60 + x61*x62;
    double x64 = -x34*x55 + x50*x52;
    double x65 = x57*(x23*U[19] + x24*U[16] + x25*U[4] + x28*U[13] + x29*U[10] + x30*U[1] + x31*U[22] + x32*U[7]);
    double x66 = -x33*x52 + x34*x38;
    double x67 = x57*(x41*U[19] + x42*U[22] + x43*U[10] + x45*U[13] + x46*U[4] + x47*U[1] + x48*U[16] + x49*U[7]);
    double x68 = x33*x55 - x38*x50;
    double x69 = x57*(x13*U[10] + x15*U[13] + x16*U[1] + x17*U[22] + x18*U[16] + x5*U[19] + x6*U[7] + x9*U[4]);
    double x70 = x64*x65 + x66*x67 + x68*x69;
    double x71 = x58*x66 + x60*x64 + x62*x68;
    double x72 = x37*x67 + x59*x65 + x61*x69 + 1.0;
    double x73 = x63*x70 - x71*x72;
    double x74 = x57*(x41*U[18] + x42*U[21] + x43*U[9] + x45*U[12] + x46*U[3] + x47*U[0] + x48*U[15] + x49*U[6]);
    double x75 = x57*(x23*U[18] + x24*U[15] + x25*U[3] + x28*U[12] + x29*U[9] + x30*U[0] + x31*U[21] + x32*U[6]);
    double x76 = x57*(x13*U[9] + x15*U[12] + x16*U[0] + x17*U[21] + x18*U[15] + x5*U[18] + x6*U[6] + x9*U[3]);
    double x77 = x37*x74 + x59*x75 + x61*x76;
    double x78 = x64*x75 + x66*x74 + x68*x76 + 1.0;
    double x79 = -x63*x78 + x71*x77;
    double x80 = -x70*x77 + x72*x78;
    double x81 = S11*x73 + S12*x79 + S13*x80;
    double x82 = S12*x73 + S22*x79 + S23*x80;
    double x83 = S13*x73 + S23*x79 + S33*x80;
    double x84 = -x70*x82 - x71*x83 - x78*x81;
    double x85 = x19*x55 - x52*x53;
    double x86 = x30*x57;
    double x87 = -x19*x38 + x35*x52;
    double x88 = x47*x57;
    double x89 = -x35*x55 + x38*x53;
    double x90 = x16*x57;
    double x91 = x85*x86 + x87*x88 + x89*x90;
    double x92 = x65*x85 + x67*x87 + x69*x89;
    double x93 = x58*x87 + x60*x85 + x62*x89 + 1.0;
    double x94 = -x70*x93 + x71*x92;
    double x95 = x74*x87 + x75*x85 + x76*x89;
    double x96 = -x71*x95 + x78*x93;
    double x97 = x70*x95 - x78*x92;
    double x98 = S11*x94 + S12*x96 + S13*x97;
    double x99 = S12*x94 + S22*x96 + S23*x97;
    double x100 = S13*x94 + S23*x96 + S33*x97;
    double x101 = -x100*x71 - x70*x99 - x78*x98;
    double x102 = x37*x88 + x59*x86 + x61*x90;
    double x103 = x64*x86 + x66*x88 + x68*x90;
    double x104 = -x63*x92 + x72*x93;
    double x105 = x63*x95 - x77*x93;
    double x106 = -x72*x95 + x77*x92;
    double x107 = S11*x104 + S12*x105 + S13*x106;
    double x108 = S12*x104 + S22*x105 + S23*x106;
    double x109 = S13*x104 + S23*x105 + S33*x106;
    double x110 = -1.0*PENER - 1.0*SENER;
    double x111 = -x107*x78 - x108*x70 - x109*x71 - x110;
    double x112 = -x63*x83 - x72*x82 - x77*x81;
    double x113 = -x107*x77 - x108*x72 - x109*x63;
    double x114 = -x100*x63 - x110 - x72*x99 - x77*x98;
    double x115 = -x100*x93 - x92*x99 - x95*x98;
    double x116 = -x107*x95 - x108*x92 - x109*x93;
    double x117 = -x110 - x81*x95 - x82*x92 - x83*x93;
    double x118 = x25*x57;
    double x119 = x46*x57;
    double x120 = x57*x9;
    double x121 = x118*x85 + x119*x87 + x120*x89;
    double x122 = x118*x59 + x119*x37 + x120*x61;
    double x123 = x118*x64 + x119*x66 + x120*x68;
    double x124 = x32*x57;
    double x125 = x49*x57;
    double x126 = x57*x6;
    double x127 = x124*x85 + x125*x87 + x126*x89;
    double x128 = x124*x59 + x125*x37 + x126*x61;
    double x129 = x124*x64 + x125*x66 + x126*x68;
    double x130 = x29*x57;
    double x131 = x43*x57;
    double x132 = x13*x57;
    double x133 = x130*x85 + x131*x87 + x132*x89;
    double x134 = x130*x59 + x131*x37 + x132*x61;
    double x135 = x130*x64 + x131*x66 + x132*x68;
    double x136 = x28*x57;
    double x137 = x45*x57;
    double x138 = x15*x57;
    double x139 = x136*x85 + x137*x87 + x138*x89;
    double x140 = x136*x59 + x137*x37 + x138*x61;
    double x141 = x136*x64 + x137*x66 + x138*x68;
    double x142 = x24*x57;
    double x143 = x48*x57;
    double x144 = x18*x57;
    double x145 = x142*x85 + x143*x87 + x144*x89;
    double x146 = x142*x59 + x143*x37 + x144*x61;
    double x147 = x142*x64 + x143*x66 + x144*x68;
    double x148 = x23*x57;
    double x149 = x41*x57;
    double x150 = x5*x57;
    double x151 = x148*x85 + x149*x87 + x150*x89;
    double x152 = x148*x59 + x149*x37 + x150*x61;
    double x153 = x148*x64 + x149*x66 + x150*x68;
    double x154 = x31*x57;
    double x155 = x42*x57;
    double x156 = x17*x57;
    double x157 = x154*x85 + x155*x87 + x156*x89;
    double x158 = x154*x59 + x155*x37 + x156*x61;
    double x159 = x154*x64 + x155*x66 + x156*x68;
    
    res_0[0] = x56*(x101*x102 + x103*x111 + x84*x91);
    res_0[1] = x56*(x102*x114 + x103*x113 + x112*x91);
    res_0[2] = x56*(x102*x115 + x103*x116 + x117*x91);
    res_0[3] = x56*(x101*x122 + x111*x123 + x121*x84);
    res_0[4] = x56*(x112*x121 + x113*x123 + x114*x122);
    res_0[5] = x56*(x115*x122 + x116*x123 + x117*x121);
    res_0[6] = x56*(x101*x128 + x111*x129 + x127*x84);
    res_0[7] = x56*(x112*x127 + x113*x129 + x114*x128);
    res_0[8] = x56*(x115*x128 + x116*x129 + x117*x127);
    res_0[9] = x56*(x101*x134 + x111*x135 + x133*x84);
    res_0[10] = x56*(x112*x133 + x113*x135 + x114*x134);
    res_0[11] = x56*(x115*x134 + x116*x135 + x117*x133);
    res_0[12] = x56*(x101*x140 + x111*x141 + x139*x84);
    res_0[13] = x56*(x112*x139 + x113*x141 + x114*x140);
    res_0[14] = x56*(x115*x140 + x116*x141 + x117*x139);
    res_0[15] = x56*(x101*x146 + x111*x147 + x145*x84);
    res_0[16] = x56*(x112*x145 + x113*x147 + x114*x146);
    res_0[17] = x56*(x115*x146 + x116*x147 + x117*x145);
    res_0[18] = x56*(x101*x152 + x111*x153 + x151*x84);
    res_0[19] = x56*(x112*x151 + x113*x153 + x114*x152);
    res_0[20] = x56*(x115*x152 + x116*x153 + x117*x151);
    res_0[21] = x56*(x101*x158 + x111*x159 + x157*x84);
    res_0[22] = x56*(x112*x157 + x113*x159 + x114*x158);
    res_0[23] = x56*(x115*x158 + x116*x159 + x117*x157);
}

Conf_Forces_API void Integration_C3D8R_static_mbf(size_t num_elem,double Coords[][8][3],double Element_U[][8][3],double S[][1][6],
    double PENER[][1],double SENER[][1],double Conf_Force[][8][3]){
    
    //Inputs:
    //    Coords.........Element_Nodal_Coordinates [num_elem][numNodes,3]
    //    Element_U......Element_Nodal_Displacements [num_elem][numNodes,3]
    //    S..............Stress at integration points [num_elem][num_int_points,6]
    //    PENER..........Plastic strain energy [num_elem][num_int_points]
    //    SENER..........Elastic strain energy [num_elem][num_int_points]
    //    
    double weights[1]={8.0,};
    double int_points[1][3]={
        {0.0,0.0,0.0,}};

    #pragma omp parallel for
    for (size_t i=0;i<num_elem;i++){
        double TMP[8][3];
        for (size_t j=0;j<1;j++){
            C3D8R_static_mbf((double*) int_points[j],(double*) Coords[i],(double*) Element_U[i],(double*) S[i][j],PENER[i][j],SENER[i][j],(double*) TMP);
            for (size_t num_node=0;num_node<8;num_node++){
                for (size_t ii=0;ii<3;ii++){
                    Conf_Force[i][num_node][ii]+=TMP[num_node][ii]*weights[j];
                }
            }
        }
    }
}

void C3D8R_static_dbf(double *rst,double *coord,double *Element_U,double *S,double PENER,double SENER,double *res_0){
    double r = rst[0];
    double s = rst[1];
    double t = rst[2];
    double *U = Element_U;
    double S11 = S[0];
    double S22 = S[1];
    double S33 = S[2];
    double S12 = S[3];
    double S13 = S[4];
    double S23 = S[5];
    
    double x0 = (1.0/8.0)*s;
    double x1 = (1.0/8.0)*r;
    double x2 = s*x1;
    double x3 = x0 + x2;
    double x4 = x1 + 1.0/8.0;
    double x5 = x3 + x4;
    double x6 = -x5;
    double x7 = -x1;
    double x8 = x7 - 1.0/8.0;
    double x9 = x3 + x8;
    double x10 = -x0;
    double x11 = x10 + x2;
    double x12 = x1 - 1.0/8.0;
    double x13 = x11 + x12;
    double x14 = x7 + 1.0/8.0;
    double x15 = x11 + x14;
    double x16 = -x15;
    double x17 = -x13;
    double x18 = -x9;
    double x19 = x13*coord[9] + x15*coord[12] + x16*coord[0] + x17*coord[21] + x18*coord[15] + x5*coord[18] + x6*coord[6] + x9*coord[3];
    double x20 = (1.0/8.0)*t;
    double x21 = t*x1;
    double x22 = x20 + x21;
    double x23 = x22 + x4;
    double x24 = -x23;
    double x25 = x22 + x8;
    double x26 = -x20;
    double x27 = x21 + x26;
    double x28 = x12 + x27;
    double x29 = x14 + x27;
    double x30 = -x29;
    double x31 = -x28;
    double x32 = -x25;
    double x33 = x23*coord[20] + x24*coord[17] + x25*coord[5] + x28*coord[14] + x29*coord[11] + x30*coord[2] + x31*coord[23] + x32*coord[8];
    double x34 = x13*coord[11] + x15*coord[14] + x16*coord[2] + x17*coord[23] + x18*coord[17] + x5*coord[20] + x6*coord[8] + x9*coord[5];
    double x35 = x23*coord[18] + x24*coord[15] + x25*coord[3] + x28*coord[12] + x29*coord[9] + x30*coord[0] + x31*coord[21] + x32*coord[6];
    double x36 = x34*x35;
    double x37 = x19*x33 - x36;
    double x38 = x23*coord[19] + x24*coord[16] + x25*coord[4] + x28*coord[13] + x29*coord[10] + x30*coord[1] + x31*coord[22] + x32*coord[7];
    double x39 = t*x0;
    double x40 = x20 + x39;
    double x41 = x0 + x40 + 1.0/8.0;
    double x42 = -x41;
    double x43 = x10 + x40 - 1.0/8.0;
    double x44 = x26 + x39;
    double x45 = x0 + x44 - 1.0/8.0;
    double x46 = x10 + x44 + 1.0/8.0;
    double x47 = -x46;
    double x48 = -x45;
    double x49 = -x43;
    double x50 = x41*coord[20] + x42*coord[23] + x43*coord[11] + x45*coord[14] + x46*coord[5] + x47*coord[2] + x48*coord[17] + x49*coord[8];
    double x51 = x19*x50;
    double x52 = x13*coord[10] + x15*coord[13] + x16*coord[1] + x17*coord[22] + x18*coord[16] + x5*coord[19] + x6*coord[7] + x9*coord[4];
    double x53 = x41*coord[18] + x42*coord[21] + x43*coord[9] + x45*coord[12] + x46*coord[3] + x47*coord[0] + x48*coord[15] + x49*coord[6];
    double x54 = x33*x53;
    double x55 = x41*coord[19] + x42*coord[22] + x43*coord[10] + x45*coord[13] + x46*coord[4] + x47*coord[1] + x48*coord[16] + x49*coord[7];
    double x56 = x19*x33*x55 + x34*x38*x53 + x35*x50*x52 - x36*x55 - x38*x51 - x52*x54;
    double x57 = 1.0/x56;
    double x58 = x57*(x41*U[20] + x42*U[23] + x43*U[11] + x45*U[14] + x46*U[5] + x47*U[2] + x48*U[17] + x49*U[8]);
    double x59 = x34*x53 - x51;
    double x60 = x57*(x23*U[20] + x24*U[17] + x25*U[5] + x28*U[14] + x29*U[11] + x30*U[2] + x31*U[23] + x32*U[8]);
    double x61 = x35*x50 - x54;
    double x62 = x57*(x13*U[11] + x15*U[14] + x16*U[2] + x17*U[23] + x18*U[17] + x5*U[20] + x6*U[8] + x9*U[5]);
    double x63 = x37*x58 + x59*x60 + x61*x62;
    double x64 = -x34*x55 + x50*x52;
    double x65 = x57*(x23*U[19] + x24*U[16] + x25*U[4] + x28*U[13] + x29*U[10] + x30*U[1] + x31*U[22] + x32*U[7]);
    double x66 = -x33*x52 + x34*x38;
    double x67 = x57*(x41*U[19] + x42*U[22] + x43*U[10] + x45*U[13] + x46*U[4] + x47*U[1] + x48*U[16] + x49*U[7]);
    double x68 = x33*x55 - x38*x50;
    double x69 = x57*(x13*U[10] + x15*U[13] + x16*U[1] + x17*U[22] + x18*U[16] + x5*U[19] + x6*U[7] + x9*U[4]);
    double x70 = x64*x65 + x66*x67 + x68*x69;
    double x71 = x58*x66 + x60*x64 + x62*x68;
    double x72 = x37*x67 + x59*x65 + x61*x69;
    double x73 = x72 + 1.0;
    double x74 = x63*x70 - x71*x73;
    double x75 = x57*(x41*U[18] + x42*U[21] + x43*U[9] + x45*U[12] + x46*U[3] + x47*U[0] + x48*U[15] + x49*U[6]);
    double x76 = x57*(x23*U[18] + x24*U[15] + x25*U[3] + x28*U[12] + x29*U[9] + x30*U[0] + x31*U[21] + x32*U[6]);
    double x77 = x57*(x13*U[9] + x15*U[12] + x16*U[0] + x17*U[21] + x18*U[15] + x5*U[18] + x6*U[6] + x9*U[3]);
    double x78 = x37*x75 + x59*x76 + x61*x77;
    double x79 = x64*x76 + x66*x75 + x68*x77;
    double x80 = x79 + 1.0;
    double x81 = -x63*x80 + x71*x78;
    double x82 = -x70*x78 + x73*x80;
    double x83 = S11*x74 + S12*x81 + S13*x82;
    double x84 = S12*x74 + S22*x81 + S23*x82;
    double x85 = S13*x74 + S23*x81 + S33*x82;
    double x86 = -x70*x84 - x71*x85 - x79*x83;
    double x87 = x19*x55 - x52*x53;
    double x88 = x30*x57;
    double x89 = -x19*x38 + x35*x52;
    double x90 = x47*x57;
    double x91 = -x35*x55 + x38*x53;
    double x92 = x16*x57;
    double x93 = x87*x88 + x89*x90 + x91*x92;
    double x94 = x65*x87 + x67*x89 + x69*x91;
    double x95 = x58*x89 + x60*x87 + x62*x91;
    double x96 = x95 + 1.0;
    double x97 = -x70*x96 + x71*x94;
    double x98 = x75*x89 + x76*x87 + x77*x91;
    double x99 = -x71*x98 + x80*x96;
    double x100 = x70*x98 - x80*x94;
    double x101 = S11*x97 + S12*x99 + S13*x100;
    double x102 = S12*x97 + S22*x99 + S23*x100;
    double x103 = S13*x97 + S23*x99 + S33*x100;
    double x104 = -x101*x79 - x102*x70 - x103*x71;
    double x105 = x37*x90 + x59*x88 + x61*x92;
    double x106 = x64*x88 + x66*x90 + x68*x92;
    double x107 = -x63*x94 + x73*x96;
    double x108 = x63*x98 - x78*x96;
    double x109 = -x73*x98 + x78*x94;
    double x110 = S11*x107 + S12*x108 + S13*x109;
    double x111 = S12*x107 + S22*x108 + S23*x109;
    double x112 = S13*x107 + S23*x108 + S33*x109;
    double x113 = -1.0*PENER - 1.0*SENER;
    double x114 = -x110*x79 - x111*x70 - x112*x71 - x113;
    double x115 = -x63*x85 - x72*x84 - x78*x83;
    double x116 = -x110*x78 - x111*x72 - x112*x63;
    double x117 = -x101*x78 - x102*x72 - x103*x63 - x113;
    double x118 = -x101*x98 - x102*x94 - x103*x95;
    double x119 = -x110*x98 - x111*x94 - x112*x95;
    double x120 = -x113 - x83*x98 - x84*x94 - x85*x95;
    double x121 = x25*x57;
    double x122 = x46*x57;
    double x123 = x57*x9;
    double x124 = x121*x87 + x122*x89 + x123*x91;
    double x125 = x121*x59 + x122*x37 + x123*x61;
    double x126 = x121*x64 + x122*x66 + x123*x68;
    double x127 = x32*x57;
    double x128 = x49*x57;
    double x129 = x57*x6;
    double x130 = x127*x87 + x128*x89 + x129*x91;
    double x131 = x127*x59 + x128*x37 + x129*x61;
    double x132 = x127*x64 + x128*x66 + x129*x68;
    double x133 = x29*x57;
    double x134 = x43*x57;
    double x135 = x13*x57;
    double x136 = x133*x87 + x134*x89 + x135*x91;
    double x137 = x133*x59 + x134*x37 + x135*x61;
    double x138 = x133*x64 + x134*x66 + x135*x68;
    double x139 = x28*x57;
    double x140 = x45*x57;
    double x141 = x15*x57;
    double x142 = x139*x87 + x140*x89 + x141*x91;
    double x143 = x139*x59 + x140*x37 + x141*x61;
    double x144 = x139*x64 + x140*x66 + x141*x68;
    double x145 = x24*x57;
    double x146 = x48*x57;
    double x147 = x18*x57;
    double x148 = x145*x87 + x146*x89 + x147*x91;
    double x149 = x145*x59 + x146*x37 + x147*x61;
    double x150 = x145*x64 + x146*x66 + x147*x68;
    double x151 = x23*x57;
    double x152 = x41*x57;
    double x153 = x5*x57;
    double x154 = x151*x87 + x152*x89 + x153*x91;
    double x155 = x151*x59 + x152*x37 + x153*x61;
    double x156 = x151*x64 + x152*x66 + x153*x68;
    double x157 = x31*x57;
    double x158 = x42*x57;
    double x159 = x17*x57;
    double x160 = x157*x87 + x158*x89 + x159*x91;
    double x161 = x157*x59 + x158*x37 + x159*x61;
    double x162 = x157*x64 + x158*x66 + x159*x68;
    
    res_0[0] = x56*(x104*x105 + x106*x114 + x86*x93);
    res_0[1] = x56*(x105*x117 + x106*x116 + x115*x93);
    res_0[2] = x56*(x105*x118 + x106*x119 + x120*x93);
    res_0[3] = x56*(x104*x125 + x114*x126 + x124*x86);
    res_0[4] = x56*(x115*x124 + x116*x126 + x117*x125);
    res_0[5] = x56*(x118*x125 + x119*x126 + x120*x124);
    res_0[6] = x56*(x104*x131 + x114*x132 + x130*x86);
    res_0[7] = x56*(x115*x130 + x116*x132 + x117*x131);
    res_0[8] = x56*(x118*x131 + x119*x132 + x120*x130);
    res_0[9] = x56*(x104*x137 + x114*x138 + x136*x86);
    res_0[10] = x56*(x115*x136 + x116*x138 + x117*x137);
    res_0[11] = x56*(x118*x137 + x119*x138 + x120*x136);
    res_0[12] = x56*(x104*x143 + x114*x144 + x142*x86);
    res_0[13] = x56*(x115*x142 + x116*x144 + x117*x143);
    res_0[14] = x56*(x118*x143 + x119*x144 + x120*x142);
    res_0[15] = x56*(x104*x149 + x114*x150 + x148*x86);
    res_0[16] = x56*(x115*x148 + x116*x150 + x117*x149);
    res_0[17] = x56*(x118*x149 + x119*x150 + x120*x148);
    res_0[18] = x56*(x104*x155 + x114*x156 + x154*x86);
    res_0[19] = x56*(x115*x154 + x116*x156 + x117*x155);
    res_0[20] = x56*(x118*x155 + x119*x156 + x120*x154);
    res_0[21] = x56*(x104*x161 + x114*x162 + x160*x86);
    res_0[22] = x56*(x115*x160 + x116*x162 + x117*x161);
    res_0[23] = x56*(x118*x161 + x119*x162 + x120*x160);
}

Conf_Forces_API void Integration_C3D8R_static_dbf(size_t num_elem,double Coords[][8][3],double Element_U[][8][3],double S[][1][6],
    double PENER[][1],double SENER[][1],double Conf_Force[][8][3]){
    
    //Inputs:
    //    Coords.........Element_Nodal_Coordinates [num_elem][numNodes,3]
    //    Element_U......Element_Nodal_Displacements [num_elem][numNodes,3]
    //    S..............Stress at integration points [num_elem][num_int_points,6]
    //    PENER..........Plastic strain energy [num_elem][num_int_points]
    //    SENER..........Elastic strain energy [num_elem][num_int_points]
    //    
    double weights[1]={8.0,};
    double int_points[1][3]={
        {0.0,0.0,0.0,}};

    #pragma omp parallel for
    for (size_t i=0;i<num_elem;i++){
        double TMP[8][3];
        for (size_t j=0;j<1;j++){
            C3D8R_static_dbf((double*) int_points[j],(double*) Coords[i],(double*) Element_U[i],(double*) S[i][j],PENER[i][j],SENER[i][j],(double*) TMP);
            for (size_t num_node=0;num_node<8;num_node++){
                for (size_t ii=0;ii<3;ii++){
                    Conf_Force[i][num_node][ii]+=TMP[num_node][ii]*weights[j];
                }
            }
        }
    }
}

void C3D20_dynamic(double *rst,double *coord,double rho,double *Element_U,double *Element_V,double *Element_A,double *S,double PENER,double SENER,double *res_0){
    double r = rst[0];
    double s = rst[1];
    double t = rst[2];
    double *U = Element_U;
    double *V = Element_V;
    double *A = Element_A;
    double S11 = S[0];
    double S22 = S[1];
    double S33 = S[2];
    double S12 = S[3];
    double S13 = S[4];
    double S23 = S[5];
    
    double x0 = (1.0/4.0)*r;
    double x1 = -x0;
    double x2 = s*s;
    double x3 = x0*x2;
    double x4 = x1 + x3;
    double x5 = (1.0/4.0)*x2;
    double x6 = x5 - 1.0/4.0;
    double x7 = x4 + x6;
    double x8 = (1.0/4.0)*s;
    double x9 = -x8;
    double x10 = r*r;
    double x11 = x10*x8;
    double x12 = x11 + x9;
    double x13 = (1.0/4.0)*x10;
    double x14 = x13 - 1.0/4.0;
    double x15 = x12 + x14;
    double x16 = 1.0/4.0 - x5;
    double x17 = x16 + x4;
    double x18 = -x17;
    double x19 = 1.0/4.0 - x13;
    double x20 = x12 + x19;
    double x21 = -x20;
    double x22 = -x15;
    double x23 = -x7;
    double x24 = (1.0/2.0)*t;
    double x25 = r*x24;
    double x26 = -x25;
    double x27 = s*x24;
    double x28 = s*x25;
    double x29 = -x27 + x28;
    double x30 = -x24 - x26 - x29;
    double x31 = -x24;
    double x32 = x25 + x29 + x31;
    double x33 = x27 + x28;
    double x34 = x26 + x31 + x33;
    double x35 = -x24 - x25 - x33;
    double x36 = (1.0/8.0)*r;
    double x37 = s*x36;
    double x38 = -x37;
    double x39 = t*x8;
    double x40 = t*x0;
    double x41 = s*x40;
    double x42 = -x41;
    double x43 = x39 + x42;
    double x44 = x38 + x43;
    double x45 = (1.0/4.0)*t;
    double x46 = -x45;
    double x47 = (1.0/8.0)*x2;
    double x48 = r*x47;
    double x49 = (1.0/8.0)*x10;
    double x50 = s*x49;
    double x51 = x48 + x50;
    double x52 = x46 + x51;
    double x53 = 1.0/8.0 - x49;
    double x54 = -x47;
    double x55 = x40 + x54;
    double x56 = x53 + x55;
    double x57 = -x44 - x52 - x56;
    double x58 = -x40;
    double x59 = x45 + x58;
    double x60 = -x48;
    double x61 = x50 + x60;
    double x62 = x49 - 1.0/8.0;
    double x63 = x47 + x62;
    double x64 = x61 + x63;
    double x65 = x44 + x59 + x64;
    double x66 = x41 + x46;
    double x67 = x37 + x61;
    double x68 = x53 + x54;
    double x69 = -x39 - x58 - x66 - x67 - x68;
    double x70 = x41 + x51;
    double x71 = x37 + x63;
    double x72 = x39 + x40 + x45 + x70 + x71;
    double x73 = -x39;
    double x74 = x38 + x73;
    double x75 = x59 + x68 + x70 + x74;
    double x76 = -x40 - x64 - x66 - x74;
    double x77 = x42 + x73;
    double x78 = x45 + x56 + x67 + x77;
    double x79 = -x52 - x58 - x71 - x77;
    double x80 = x15*coord[30] + x17*coord[45] + x18*coord[33] + x20*coord[36] + x21*coord[24] + x22*coord[42] + x23*coord[39] + x30*coord[48] + x32*coord[57] + x34*coord[51] + x35*coord[54] + x57*coord[12] + x65*coord[21] + x69*coord[15] + x7*coord[27] + x72*coord[18] + x75*coord[0] + x76*coord[9] + x78*coord[3] + x79*coord[6];
    double x81 = t*t;
    double x82 = x0*x81;
    double x83 = x1 + x82;
    double x84 = (1.0/4.0)*x81;
    double x85 = x84 - 1.0/4.0;
    double x86 = x83 + x85;
    double x87 = t*x13;
    double x88 = x46 + x87;
    double x89 = x14 + x88;
    double x90 = 1.0/4.0 - x84;
    double x91 = x83 + x90;
    double x92 = -x91;
    double x93 = x19 + x88;
    double x94 = -x93;
    double x95 = -x89;
    double x96 = -x86;
    double x97 = (1.0/2.0)*s;
    double x98 = r*x97;
    double x99 = -x98;
    double x100 = -x29 - x97 - x99;
    double x101 = -x97;
    double x102 = x101 + x29 + x98;
    double x103 = x101 + x33 + x99;
    double x104 = -x33 - x97 - x98;
    double x105 = t*x36;
    double x106 = -x105;
    double x107 = x106 + x43;
    double x108 = (1.0/8.0)*x81;
    double x109 = r*x108;
    double x110 = t*x49;
    double x111 = x109 + x110;
    double x112 = x111 + x9;
    double x113 = -x108;
    double x114 = s*x0;
    double x115 = x113 + x114;
    double x116 = x115 + x53;
    double x117 = -x107 - x112 - x116;
    double x118 = -x114;
    double x119 = x118 + x8;
    double x120 = -x109;
    double x121 = x108 + x120;
    double x122 = x110 + x121;
    double x123 = x122 + x62;
    double x124 = x107 + x119 + x123;
    double x125 = x41 + x9;
    double x126 = x105 + x118;
    double x127 = x113 + x53;
    double x128 = x110 + x120;
    double x129 = -x125 - x126 - x127 - x128 - x39;
    double x130 = x111 + x41;
    double x131 = x105 + x8;
    double x132 = x108 + x62;
    double x133 = x114 + x130 + x131 + x132 + x39;
    double x134 = x106 + x73;
    double x135 = x119 + x127 + x130 + x134;
    double x136 = -x114 - x123 - x125 - x134;
    double x137 = x116 + x128 + x131 + x77;
    double x138 = -x112 - x126 - x132 - x77;
    double x139 = x100*coord[35] + x102*coord[47] + x103*coord[29] + x104*coord[41] + x117*coord[11] + x124*coord[23] + x129*coord[8] + x133*coord[20] + x135*coord[2] + x136*coord[14] + x137*coord[5] + x138*coord[17] + x86*coord[53] + x89*coord[38] + x91*coord[59] + x92*coord[50] + x93*coord[32] + x94*coord[26] + x95*coord[44] + x96*coord[56];
    double x140 = x15*coord[32] + x17*coord[47] + x18*coord[35] + x20*coord[38] + x21*coord[26] + x22*coord[44] + x23*coord[41] + x30*coord[50] + x32*coord[59] + x34*coord[53] + x35*coord[56] + x57*coord[14] + x65*coord[23] + x69*coord[17] + x7*coord[29] + x72*coord[20] + x75*coord[2] + x76*coord[11] + x78*coord[5] + x79*coord[8];
    double x141 = x100*coord[33] + x102*coord[45] + x103*coord[27] + x104*coord[39] + x117*coord[9] + x124*coord[21] + x129*coord[6] + x133*coord[18] + x135*coord[0] + x136*coord[12] + x137*coord[3] + x138*coord[15] + x86*coord[51] + x89*coord[36] + x91*coord[57] + x92*coord[48] + x93*coord[30] + x94*coord[24] + x95*coord[42] + x96*coord[54];
    double x142 = x140*x141;
    double x143 = x139*x80 - x142;
    double x144 = x100*coord[34] + x102*coord[46] + x103*coord[28] + x104*coord[40] + x117*coord[10] + x124*coord[22] + x129*coord[7] + x133*coord[19] + x135*coord[1] + x136*coord[13] + x137*coord[4] + x138*coord[16] + x86*coord[52] + x89*coord[37] + x91*coord[58] + x92*coord[49] + x93*coord[31] + x94*coord[25] + x95*coord[43] + x96*coord[55];
    double x145 = x8*x81;
    double x146 = x145 + x9;
    double x147 = x146 + x85;
    double x148 = t*x5;
    double x149 = x148 + x46;
    double x150 = x149 + x6;
    double x151 = x146 + x90;
    double x152 = -x151;
    double x153 = x149 + x16;
    double x154 = -x153;
    double x155 = -x150;
    double x156 = -x147;
    double x157 = (1.0/2.0)*r;
    double x158 = x26 + x28;
    double x159 = -x157 - x158 - x99;
    double x160 = -x157;
    double x161 = x158 + x160 + x98;
    double x162 = x25 + x28;
    double x163 = x160 + x162 + x99;
    double x164 = -x157 - x162 - x98;
    double x165 = t*x47;
    double x166 = -x165;
    double x167 = x108 + x47;
    double x168 = x166 + x167;
    double x169 = (1.0/8.0)*s;
    double x170 = t*x169;
    double x171 = x170 + x41 + x58;
    double x172 = s*x108;
    double x173 = -x172;
    double x174 = x118 + x173;
    double x175 = x0 - 1.0/8.0;
    double x176 = x168 + x171 + x174 + x175;
    double x177 = x1 + 1.0/8.0;
    double x178 = x166 + x172;
    double x179 = -x115 - x171 - x177 - x178 - x54;
    double x180 = x165 + x41;
    double x181 = x170 + x180;
    double x182 = -x113 - x174 - x177 - x181 - x55;
    double x183 = x114 + x167;
    double x184 = x172 + x40;
    double x185 = x175 + x181 + x183 + x184;
    double x186 = -x170;
    double x187 = x118 + x186;
    double x188 = x180 + x58;
    double x189 = x0 + 1.0/8.0;
    double x190 = x113 + x172 + x187 + x188 + x189 + x54;
    double x191 = x1 - 1.0/8.0;
    double x192 = x173 + x186;
    double x193 = -x183 - x188 - x191 - x192;
    double x194 = -x168 - x184 - x187 - x191 - x41;
    double x195 = x115 + x166 + x189 + x192 + x41 + x55;
    double x196 = x147*coord[59] + x150*coord[47] + x151*coord[53] + x152*coord[50] + x153*coord[29] + x154*coord[35] + x155*coord[41] + x156*coord[56] + x159*coord[26] + x161*coord[38] + x163*coord[32] + x164*coord[44] + x176*coord[5] + x179*coord[17] + x182*coord[8] + x185*coord[20] + x190*coord[2] + x193*coord[14] + x194*coord[11] + x195*coord[23];
    double x197 = x196*x80;
    double x198 = x15*coord[31] + x17*coord[46] + x18*coord[34] + x20*coord[37] + x21*coord[25] + x22*coord[43] + x23*coord[40] + x30*coord[49] + x32*coord[58] + x34*coord[52] + x35*coord[55] + x57*coord[13] + x65*coord[22] + x69*coord[16] + x7*coord[28] + x72*coord[19] + x75*coord[1] + x76*coord[10] + x78*coord[4] + x79*coord[7];
    double x199 = x147*coord[57] + x150*coord[45] + x151*coord[51] + x152*coord[48] + x153*coord[27] + x154*coord[33] + x155*coord[39] + x156*coord[54] + x159*coord[24] + x161*coord[36] + x163*coord[30] + x164*coord[42] + x176*coord[3] + x179*coord[15] + x182*coord[6] + x185*coord[18] + x190*coord[0] + x193*coord[12] + x194*coord[9] + x195*coord[21];
    double x200 = x139*x199;
    double x201 = x147*coord[58] + x150*coord[46] + x151*coord[52] + x152*coord[49] + x153*coord[28] + x154*coord[34] + x155*coord[40] + x156*coord[55] + x159*coord[25] + x161*coord[37] + x163*coord[31] + x164*coord[43] + x176*coord[4] + x179*coord[16] + x182*coord[7] + x185*coord[19] + x190*coord[1] + x193*coord[13] + x194*coord[10] + x195*coord[22];
    double x202 = x139*x201*x80 + x140*x144*x199 + x141*x196*x198 - x142*x201 - x144*x197 - x198*x200;
    double x203 = 1.0/x202;
    double x204 = x203*(x147*U[59] + x150*U[47] + x151*U[53] + x152*U[50] + x153*U[29] + x154*U[35] + x155*U[41] + x156*U[56] + x159*U[26] + x161*U[38] + x163*U[32] + x164*U[44] + x176*U[5] + x179*U[17] + x182*U[8] + x185*U[20] + x190*U[2] + x193*U[14] + x194*U[11] + x195*U[23]);
    double x205 = x140*x199 - x197;
    double x206 = x203*(x100*U[35] + x102*U[47] + x103*U[29] + x104*U[41] + x117*U[11] + x124*U[23] + x129*U[8] + x133*U[20] + x135*U[2] + x136*U[14] + x137*U[5] + x138*U[17] + x86*U[53] + x89*U[38] + x91*U[59] + x92*U[50] + x93*U[32] + x94*U[26] + x95*U[44] + x96*U[56]);
    double x207 = x141*x196 - x200;
    double x208 = x203*(x15*U[32] + x17*U[47] + x18*U[35] + x20*U[38] + x21*U[26] + x22*U[44] + x23*U[41] + x30*U[50] + x32*U[59] + x34*U[53] + x35*U[56] + x57*U[14] + x65*U[23] + x69*U[17] + x7*U[29] + x72*U[20] + x75*U[2] + x76*U[11] + x78*U[5] + x79*U[8]);
    double x209 = x143*x204 + x205*x206 + x207*x208;
    double x210 = -x140*x201 + x196*x198;
    double x211 = x203*(x100*U[34] + x102*U[46] + x103*U[28] + x104*U[40] + x117*U[10] + x124*U[22] + x129*U[7] + x133*U[19] + x135*U[1] + x136*U[13] + x137*U[4] + x138*U[16] + x86*U[52] + x89*U[37] + x91*U[58] + x92*U[49] + x93*U[31] + x94*U[25] + x95*U[43] + x96*U[55]);
    double x212 = -x139*x198 + x140*x144;
    double x213 = x203*(x147*U[58] + x150*U[46] + x151*U[52] + x152*U[49] + x153*U[28] + x154*U[34] + x155*U[40] + x156*U[55] + x159*U[25] + x161*U[37] + x163*U[31] + x164*U[43] + x176*U[4] + x179*U[16] + x182*U[7] + x185*U[19] + x190*U[1] + x193*U[13] + x194*U[10] + x195*U[22]);
    double x214 = x139*x201 - x144*x196;
    double x215 = x203*(x15*U[31] + x17*U[46] + x18*U[34] + x20*U[37] + x21*U[25] + x22*U[43] + x23*U[40] + x30*U[49] + x32*U[58] + x34*U[52] + x35*U[55] + x57*U[13] + x65*U[22] + x69*U[16] + x7*U[28] + x72*U[19] + x75*U[1] + x76*U[10] + x78*U[4] + x79*U[7]);
    double x216 = x210*x211 + x212*x213 + x214*x215;
    double x217 = x204*x212 + x206*x210 + x208*x214;
    double x218 = x143*x213 + x205*x211 + x207*x215;
    double x219 = x218 + 1.0;
    double x220 = x209*x216 - x217*x219;
    double x221 = x203*(x147*U[57] + x150*U[45] + x151*U[51] + x152*U[48] + x153*U[27] + x154*U[33] + x155*U[39] + x156*U[54] + x159*U[24] + x161*U[36] + x163*U[30] + x164*U[42] + x176*U[3] + x179*U[15] + x182*U[6] + x185*U[18] + x190*U[0] + x193*U[12] + x194*U[9] + x195*U[21]);
    double x222 = x203*(x100*U[33] + x102*U[45] + x103*U[27] + x104*U[39] + x117*U[9] + x124*U[21] + x129*U[6] + x133*U[18] + x135*U[0] + x136*U[12] + x137*U[3] + x138*U[15] + x86*U[51] + x89*U[36] + x91*U[57] + x92*U[48] + x93*U[30] + x94*U[24] + x95*U[42] + x96*U[54]);
    double x223 = x203*(x15*U[30] + x17*U[45] + x18*U[33] + x20*U[36] + x21*U[24] + x22*U[42] + x23*U[39] + x30*U[48] + x32*U[57] + x34*U[51] + x35*U[54] + x57*U[12] + x65*U[21] + x69*U[15] + x7*U[27] + x72*U[18] + x75*U[0] + x76*U[9] + x78*U[3] + x79*U[6]);
    double x224 = x143*x221 + x205*x222 + x207*x223;
    double x225 = x210*x222 + x212*x221 + x214*x223;
    double x226 = x225 + 1.0;
    double x227 = -x209*x226 + x217*x224;
    double x228 = -x216*x224 + x219*x226;
    double x229 = S11*x220 + S12*x227 + S13*x228;
    double x230 = S12*x220 + S22*x227 + S23*x228;
    double x231 = S13*x220 + S23*x227 + S33*x228;
    double x232 = -x216*x230 - x217*x231 - x225*x229;
    double x233 = -x198*x199 + x201*x80;
    double x234 = x135*x203;
    double x235 = x141*x198 - x144*x80;
    double x236 = x190*x203;
    double x237 = -x141*x201 + x144*x199;
    double x238 = x203*x75;
    double x239 = x233*x234 + x235*x236 + x237*x238;
    double x240 = x211*x233 + x213*x235 + x215*x237;
    double x241 = x204*x235 + x206*x233 + x208*x237;
    double x242 = x241 + 1.0;
    double x243 = -x216*x242 + x217*x240;
    double x244 = x221*x235 + x222*x233 + x223*x237;
    double x245 = -x217*x244 + x226*x242;
    double x246 = x216*x244 - x226*x240;
    double x247 = S11*x243 + S12*x245 + S13*x246;
    double x248 = S12*x243 + S22*x245 + S23*x246;
    double x249 = S13*x243 + S23*x245 + S33*x246;
    double x250 = -x216*x248 - x217*x249 - x225*x247;
    double x251 = x143*x236 + x205*x234 + x207*x238;
    double x252 = x210*x234 + x212*x236 + x214*x238;
    double x253 = -x209*x240 + x219*x242;
    double x254 = x209*x244 - x224*x242;
    double x255 = -x219*x244 + x224*x240;
    double x256 = S11*x253 + S12*x254 + S13*x255;
    double x257 = S12*x253 + S22*x254 + S23*x255;
    double x258 = S13*x253 + S23*x254 + S33*x255;
    double x259 = s*x82 + x118;
    double x260 = -x145 + x259 + x8;
    double x261 = x0 - x82;
    double x262 = -x260 - x261 - x85;
    double x263 = t*x3 + x58;
    double x264 = -x148 + x263 + x45;
    double x265 = x0 - x3;
    double x266 = -x264 - x265 - x6;
    double x267 = x17 + x264;
    double x268 = x260 + x91;
    double x269 = t*x11 + x73;
    double x270 = x269 + x45 - x87;
    double x271 = -x11 + x8;
    double x272 = -x14 - x270 - x271;
    double x273 = x20 + x270;
    double x274 = x151 + x259 + x261;
    double x275 = x269 + x271 + x93;
    double x276 = x153 + x263 + x265;
    double x277 = -x15 - x269 - x88;
    double x278 = -x149 - x263 - x7;
    double x279 = -x146 - x259 - x86;
    double x280 = -x36;
    double x281 = t*x37;
    double x282 = t*x50;
    double x283 = x280 + x281 + x282 + x49 - 1.0/4.0;
    double x284 = -x169;
    double x285 = t*x48;
    double x286 = x284 + x285;
    double x287 = (1.0/8.0)*t;
    double x288 = -x287;
    double x289 = s*x109;
    double x290 = x165 + x288 + x289;
    double x291 = x172 + x51;
    double x292 = x111 + x167;
    double x293 = x283 + x286 + x290 + x291 + x292;
    double x294 = -x285;
    double x295 = -x282;
    double x296 = x294 + x295 + x49 - 1.0/4.0;
    double x297 = -x281;
    double x298 = x280 + x297;
    double x299 = -x110;
    double x300 = x109 + x168 + x287 + x299;
    double x301 = x284 + x289 + x291 + x296 + x298 + x300;
    double x302 = -x289;
    double x303 = x295 + x302 + x49 - 1.0/4.0;
    double x304 = x169 + x173 - x50;
    double x305 = x304 + x48;
    double x306 = x165 + x288;
    double x307 = x285 + x292 + x298 + x303 + x305 + x306;
    double x308 = x294 + x302;
    double x309 = x283 + x300 + x305 + x308;
    double x310 = x36 + x47 + x61;
    double x311 = x282 + x297 + x49 - 1.0/4.0;
    double x312 = x122 + x172 + x284 + x306 + x308 + x310 + x311;
    double x313 = x121 + x287 + x299;
    double x314 = x178 + x281 + x286 + x303 + x310 + x313;
    double x315 = x304 + x36 + x47 + x60;
    double x316 = x122 + x281 + x290 + x296 + x315;
    double x317 = x166 + x285 + x289 + x311 + x313 + x315;
    double x318 = x262*V[48] + x266*V[33] + x267*V[45] + x268*V[57] + x272*V[24] + x273*V[36] + x274*V[51] + x275*V[30] + x276*V[27] + x277*V[42] + x278*V[39] + x279*V[54] + x293*V[18] + x301*V[6] + x307*V[15] + x309*V[3] + x312*V[21] + x314*V[9] + x316*V[12] + x317*V[0];
    double x319 = x262*V[49] + x266*V[34] + x267*V[46] + x268*V[58] + x272*V[25] + x273*V[37] + x274*V[52] + x275*V[31] + x276*V[28] + x277*V[43] + x278*V[40] + x279*V[55] + x293*V[19] + x301*V[7] + x307*V[16] + x309*V[4] + x312*V[22] + x314*V[10] + x316*V[13] + x317*V[1];
    double x320 = x262*V[50] + x266*V[35] + x267*V[47] + x268*V[59] + x272*V[26] + x273*V[38] + x274*V[53] + x275*V[32] + x276*V[29] + x277*V[44] + x278*V[41] + x279*V[56] + x293*V[20] + x301*V[8] + x307*V[17] + x309*V[5] + x312*V[23] + x314*V[11] + x316*V[14] + x317*V[2];
    double x321 = -1.0*PENER - 1.0*SENER + 0.5*rho*(x318*x318 + x319*x319 + x320*x320);
    double x322 = -x216*x257 - x217*x258 - x225*x256 - x321;
    double x323 = x262*A[48] + x266*A[33] + x267*A[45] + x268*A[57] + x272*A[24] + x273*A[36] + x274*A[51] + x275*A[30] + x276*A[27] + x277*A[42] + x278*A[39] + x279*A[54] + x293*A[18] + x301*A[6] + x307*A[15] + x309*A[3] + x312*A[21] + x314*A[9] + x316*A[12] + x317*A[0];
    double x324 = x262*A[49] + x266*A[34] + x267*A[46] + x268*A[58] + x272*A[25] + x273*A[37] + x274*A[52] + x275*A[31] + x276*A[28] + x277*A[43] + x278*A[40] + x279*A[55] + x293*A[19] + x301*A[7] + x307*A[16] + x309*A[4] + x312*A[22] + x314*A[10] + x316*A[13] + x317*A[1];
    double x325 = x262*A[50] + x266*A[35] + x267*A[47] + x268*A[59] + x272*A[26] + x273*A[38] + x274*A[53] + x275*A[32] + x276*A[29] + x277*A[44] + x278*A[41] + x279*A[56] + x293*A[20] + x301*A[8] + x307*A[17] + x309*A[5] + x312*A[23] + x314*A[11] + x316*A[14] + x317*A[2];
    double x326 = x100*V[33] + x102*V[45] + x103*V[27] + x104*V[39] + x117*V[9] + x124*V[21] + x129*V[6] + x133*V[18] + x135*V[0] + x136*V[12] + x137*V[3] + x138*V[15] + x86*V[51] + x89*V[36] + x91*V[57] + x92*V[48] + x93*V[30] + x94*V[24] + x95*V[42] + x96*V[54];
    double x327 = x203*x210;
    double x328 = x147*V[57] + x150*V[45] + x151*V[51] + x152*V[48] + x153*V[27] + x154*V[33] + x155*V[39] + x156*V[54] + x159*V[24] + x161*V[36] + x163*V[30] + x164*V[42] + x176*V[3] + x179*V[15] + x182*V[6] + x185*V[18] + x190*V[0] + x193*V[12] + x194*V[9] + x195*V[21];
    double x329 = x203*x212;
    double x330 = x15*V[30] + x17*V[45] + x18*V[33] + x20*V[36] + x21*V[24] + x22*V[42] + x23*V[39] + x30*V[48] + x32*V[57] + x34*V[51] + x35*V[54] + x57*V[12] + x65*V[21] + x69*V[15] + x7*V[27] + x72*V[18] + x75*V[0] + x76*V[9] + x78*V[3] + x79*V[6];
    double x331 = x203*x214;
    double x332 = x100*V[34] + x102*V[46] + x103*V[28] + x104*V[40] + x117*V[10] + x124*V[22] + x129*V[7] + x133*V[19] + x135*V[1] + x136*V[13] + x137*V[4] + x138*V[16] + x86*V[52] + x89*V[37] + x91*V[58] + x92*V[49] + x93*V[31] + x94*V[25] + x95*V[43] + x96*V[55];
    double x333 = x147*V[58] + x150*V[46] + x151*V[52] + x152*V[49] + x153*V[28] + x154*V[34] + x155*V[40] + x156*V[55] + x159*V[25] + x161*V[37] + x163*V[31] + x164*V[43] + x176*V[4] + x179*V[16] + x182*V[7] + x185*V[19] + x190*V[1] + x193*V[13] + x194*V[10] + x195*V[22];
    double x334 = x15*V[31] + x17*V[46] + x18*V[34] + x20*V[37] + x21*V[25] + x22*V[43] + x23*V[40] + x30*V[49] + x32*V[58] + x34*V[52] + x35*V[55] + x57*V[13] + x65*V[22] + x69*V[16] + x7*V[28] + x72*V[19] + x75*V[1] + x76*V[10] + x78*V[4] + x79*V[7];
    double x335 = x100*V[35] + x102*V[47] + x103*V[29] + x104*V[41] + x117*V[11] + x124*V[23] + x129*V[8] + x133*V[20] + x135*V[2] + x136*V[14] + x137*V[5] + x138*V[17] + x86*V[53] + x89*V[38] + x91*V[59] + x92*V[50] + x93*V[32] + x94*V[26] + x95*V[44] + x96*V[56];
    double x336 = x147*V[59] + x150*V[47] + x151*V[53] + x152*V[50] + x153*V[29] + x154*V[35] + x155*V[41] + x156*V[56] + x159*V[26] + x161*V[38] + x163*V[32] + x164*V[44] + x176*V[5] + x179*V[17] + x182*V[8] + x185*V[20] + x190*V[2] + x193*V[14] + x194*V[11] + x195*V[23];
    double x337 = x15*V[32] + x17*V[47] + x18*V[35] + x20*V[38] + x21*V[26] + x22*V[44] + x23*V[41] + x30*V[50] + x32*V[59] + x34*V[53] + x35*V[56] + x57*V[14] + x65*V[23] + x69*V[17] + x7*V[29] + x72*V[20] + x75*V[2] + x76*V[11] + x78*V[5] + x79*V[8];
    double x338 = x216*x324 + x217*x325 + x225*x323 + x318*(x326*x327 + x328*x329 + x330*x331) + x319*(x327*x332 + x329*x333 + x331*x334) + x320*(x327*x335 + x329*x336 + x331*x337);
    double x339 = rho*x317;
    double x340 = -x209*x231 - x218*x230 - x224*x229;
    double x341 = -x209*x258 - x218*x257 - x224*x256;
    double x342 = -x209*x249 - x218*x248 - x224*x247 - x321;
    double x343 = x143*x203;
    double x344 = x203*x205;
    double x345 = x203*x207;
    double x346 = x209*x325 + x218*x324 + x224*x323 + x318*(x326*x344 + x328*x343 + x330*x345) + x319*(x332*x344 + x333*x343 + x334*x345) + x320*(x335*x344 + x336*x343 + x337*x345);
    double x347 = -x240*x248 - x241*x249 - x244*x247;
    double x348 = -x240*x257 - x241*x258 - x244*x256;
    double x349 = -x229*x244 - x230*x240 - x231*x241 - x321;
    double x350 = x203*x233;
    double x351 = x203*x235;
    double x352 = x203*x237;
    double x353 = x240*x324 + x241*x325 + x244*x323 + x318*(x326*x350 + x328*x351 + x330*x352) + x319*(x332*x350 + x333*x351 + x334*x352) + x320*(x335*x350 + x336*x351 + x337*x352);
    double x354 = x137*x350 + x176*x351 + x352*x78;
    double x355 = x137*x344 + x176*x343 + x345*x78;
    double x356 = x137*x327 + x176*x329 + x331*x78;
    double x357 = rho*x309;
    double x358 = x129*x350 + x182*x351 + x352*x79;
    double x359 = x129*x344 + x182*x343 + x345*x79;
    double x360 = x129*x327 + x182*x329 + x331*x79;
    double x361 = rho*x301;
    double x362 = x117*x350 + x194*x351 + x352*x76;
    double x363 = x117*x344 + x194*x343 + x345*x76;
    double x364 = x117*x327 + x194*x329 + x331*x76;
    double x365 = rho*x314;
    double x366 = x136*x350 + x193*x351 + x352*x57;
    double x367 = x136*x344 + x193*x343 + x345*x57;
    double x368 = x136*x327 + x193*x329 + x331*x57;
    double x369 = rho*x316;
    double x370 = x138*x350 + x179*x351 + x352*x69;
    double x371 = x138*x344 + x179*x343 + x345*x69;
    double x372 = x138*x327 + x179*x329 + x331*x69;
    double x373 = rho*x307;
    double x374 = x133*x350 + x185*x351 + x352*x72;
    double x375 = x133*x344 + x185*x343 + x345*x72;
    double x376 = x133*x327 + x185*x329 + x331*x72;
    double x377 = rho*x293;
    double x378 = x124*x350 + x195*x351 + x352*x65;
    double x379 = x124*x344 + x195*x343 + x345*x65;
    double x380 = x124*x327 + x195*x329 + x331*x65;
    double x381 = rho*x312;
    double x382 = x159*x351 + x21*x352 + x350*x94;
    double x383 = x159*x343 + x21*x345 + x344*x94;
    double x384 = x159*x329 + x21*x331 + x327*x94;
    double x385 = rho*x272;
    double x386 = x103*x350 + x153*x351 + x352*x7;
    double x387 = x103*x344 + x153*x343 + x345*x7;
    double x388 = x103*x327 + x153*x329 + x331*x7;
    double x389 = rho*x276;
    double x390 = x15*x352 + x163*x351 + x350*x93;
    double x391 = x15*x345 + x163*x343 + x344*x93;
    double x392 = x15*x331 + x163*x329 + x327*x93;
    double x393 = rho*x275;
    double x394 = x100*x350 + x154*x351 + x18*x352;
    double x395 = x100*x344 + x154*x343 + x18*x345;
    double x396 = x100*x327 + x154*x329 + x18*x331;
    double x397 = rho*x266;
    double x398 = x161*x351 + x20*x352 + x350*x89;
    double x399 = x161*x343 + x20*x345 + x344*x89;
    double x400 = x161*x329 + x20*x331 + x327*x89;
    double x401 = rho*x273;
    double x402 = x104*x350 + x155*x351 + x23*x352;
    double x403 = x104*x344 + x155*x343 + x23*x345;
    double x404 = x104*x327 + x155*x329 + x23*x331;
    double x405 = rho*x278;
    double x406 = x164*x351 + x22*x352 + x350*x95;
    double x407 = x164*x343 + x22*x345 + x344*x95;
    double x408 = x164*x329 + x22*x331 + x327*x95;
    double x409 = rho*x277;
    double x410 = x102*x350 + x150*x351 + x17*x352;
    double x411 = x102*x344 + x150*x343 + x17*x345;
    double x412 = x102*x327 + x150*x329 + x17*x331;
    double x413 = rho*x267;
    double x414 = x152*x351 + x30*x352 + x350*x92;
    double x415 = x152*x343 + x30*x345 + x344*x92;
    double x416 = x152*x329 + x30*x331 + x327*x92;
    double x417 = rho*x262;
    double x418 = x151*x351 + x34*x352 + x350*x86;
    double x419 = x151*x343 + x34*x345 + x344*x86;
    double x420 = x151*x329 + x327*x86 + x331*x34;
    double x421 = rho*x274;
    double x422 = x156*x351 + x35*x352 + x350*x96;
    double x423 = x156*x343 + x344*x96 + x345*x35;
    double x424 = x156*x329 + x327*x96 + x331*x35;
    double x425 = rho*x279;
    double x426 = x147*x351 + x32*x352 + x350*x91;
    double x427 = x147*x343 + x32*x345 + x344*x91;
    double x428 = x147*x329 + x32*x331 + x327*x91;
    double x429 = rho*x268;
    
    res_0[0] = x202*(x232*x239 + x250*x251 + x252*x322 - x338*x339);
    res_0[1] = x202*(x239*x340 + x251*x342 + x252*x341 - x339*x346);
    res_0[2] = x202*(x239*x349 + x251*x347 + x252*x348 - x339*x353);
    res_0[3] = x202*(x232*x354 + x250*x355 + x322*x356 - x338*x357);
    res_0[4] = x202*(x340*x354 + x341*x356 + x342*x355 - x346*x357);
    res_0[5] = x202*(x347*x355 + x348*x356 + x349*x354 - x353*x357);
    res_0[6] = x202*(x232*x358 + x250*x359 + x322*x360 - x338*x361);
    res_0[7] = x202*(x340*x358 + x341*x360 + x342*x359 - x346*x361);
    res_0[8] = x202*(x347*x359 + x348*x360 + x349*x358 - x353*x361);
    res_0[9] = x202*(x232*x362 + x250*x363 + x322*x364 - x338*x365);
    res_0[10] = x202*(x340*x362 + x341*x364 + x342*x363 - x346*x365);
    res_0[11] = x202*(x347*x363 + x348*x364 + x349*x362 - x353*x365);
    res_0[12] = x202*(x232*x366 + x250*x367 + x322*x368 - x338*x369);
    res_0[13] = x202*(x340*x366 + x341*x368 + x342*x367 - x346*x369);
    res_0[14] = x202*(x347*x367 + x348*x368 + x349*x366 - x353*x369);
    res_0[15] = x202*(x232*x370 + x250*x371 + x322*x372 - x338*x373);
    res_0[16] = x202*(x340*x370 + x341*x372 + x342*x371 - x346*x373);
    res_0[17] = x202*(x347*x371 + x348*x372 + x349*x370 - x353*x373);
    res_0[18] = x202*(x232*x374 + x250*x375 + x322*x376 - x338*x377);
    res_0[19] = x202*(x340*x374 + x341*x376 + x342*x375 - x346*x377);
    res_0[20] = x202*(x347*x375 + x348*x376 + x349*x374 - x353*x377);
    res_0[21] = x202*(x232*x378 + x250*x379 + x322*x380 - x338*x381);
    res_0[22] = x202*(x340*x378 + x341*x380 + x342*x379 - x346*x381);
    res_0[23] = x202*(x347*x379 + x348*x380 + x349*x378 - x353*x381);
    res_0[24] = x202*(x232*x382 + x250*x383 + x322*x384 - x338*x385);
    res_0[25] = x202*(x340*x382 + x341*x384 + x342*x383 - x346*x385);
    res_0[26] = x202*(x347*x383 + x348*x384 + x349*x382 - x353*x385);
    res_0[27] = x202*(x232*x386 + x250*x387 + x322*x388 - x338*x389);
    res_0[28] = x202*(x340*x386 + x341*x388 + x342*x387 - x346*x389);
    res_0[29] = x202*(x347*x387 + x348*x388 + x349*x386 - x353*x389);
    res_0[30] = x202*(x232*x390 + x250*x391 + x322*x392 - x338*x393);
    res_0[31] = x202*(x340*x390 + x341*x392 + x342*x391 - x346*x393);
    res_0[32] = x202*(x347*x391 + x348*x392 + x349*x390 - x353*x393);
    res_0[33] = x202*(x232*x394 + x250*x395 + x322*x396 - x338*x397);
    res_0[34] = x202*(x340*x394 + x341*x396 + x342*x395 - x346*x397);
    res_0[35] = x202*(x347*x395 + x348*x396 + x349*x394 - x353*x397);
    res_0[36] = x202*(x232*x398 + x250*x399 + x322*x400 - x338*x401);
    res_0[37] = x202*(x340*x398 + x341*x400 + x342*x399 - x346*x401);
    res_0[38] = x202*(x347*x399 + x348*x400 + x349*x398 - x353*x401);
    res_0[39] = x202*(x232*x402 + x250*x403 + x322*x404 - x338*x405);
    res_0[40] = x202*(x340*x402 + x341*x404 + x342*x403 - x346*x405);
    res_0[41] = x202*(x347*x403 + x348*x404 + x349*x402 - x353*x405);
    res_0[42] = x202*(x232*x406 + x250*x407 + x322*x408 - x338*x409);
    res_0[43] = x202*(x340*x406 + x341*x408 + x342*x407 - x346*x409);
    res_0[44] = x202*(x347*x407 + x348*x408 + x349*x406 - x353*x409);
    res_0[45] = x202*(x232*x410 + x250*x411 + x322*x412 - x338*x413);
    res_0[46] = x202*(x340*x410 + x341*x412 + x342*x411 - x346*x413);
    res_0[47] = x202*(x347*x411 + x348*x412 + x349*x410 - x353*x413);
    res_0[48] = x202*(x232*x414 + x250*x415 + x322*x416 - x338*x417);
    res_0[49] = x202*(x340*x414 + x341*x416 + x342*x415 - x346*x417);
    res_0[50] = x202*(x347*x415 + x348*x416 + x349*x414 - x353*x417);
    res_0[51] = x202*(x232*x418 + x250*x419 + x322*x420 - x338*x421);
    res_0[52] = x202*(x340*x418 + x341*x420 + x342*x419 - x346*x421);
    res_0[53] = x202*(x347*x419 + x348*x420 + x349*x418 - x353*x421);
    res_0[54] = x202*(x232*x422 + x250*x423 + x322*x424 - x338*x425);
    res_0[55] = x202*(x340*x422 + x341*x424 + x342*x423 - x346*x425);
    res_0[56] = x202*(x347*x423 + x348*x424 + x349*x422 - x353*x425);
    res_0[57] = x202*(x232*x426 + x250*x427 + x322*x428 - x338*x429);
    res_0[58] = x202*(x340*x426 + x341*x428 + x342*x427 - x346*x429);
    res_0[59] = x202*(x347*x427 + x348*x428 + x349*x426 - x353*x429);
}

Conf_Forces_API void Integration_C3D20_dynamic(size_t num_elem,double Coords[][20][3],double *rho,double Element_U[][20][3],double Element_V[][20][3],double Element_A[][20][3],double S[][27][6],
    double PENER[][27],double SENER[][27],double Conf_Force[][20][3]){
    
    //Inputs:
    //    Coords.........Element_Nodal_Coordinates  [num_elem][numNodes,3]
    //    Element_U......Element_Nodal_Displacements[num_elem][numNodes,3]
    //    Element_V......Element_Nodal_Velocity     [num_elem][numNodes,3]
    //    Element_A......Element_Nodal_Acceleration [num_elem][numNodes,3]
    //    S..............Stress at integration points [num_elem][num_int_points,6]
    //    PENER..........Plastic strain energy [num_elem][num_int_points]
    //    SENER..........Elastic strain energy [num_elem][num_int_points]
    //    
    double weights[27]={0.1714677640603567,0.27434842249657065,0.1714677640603567,0.27434842249657065,0.438957475994513,0.27434842249657065,0.1714677640603567,0.27434842249657065,0.1714677640603567,0.27434842249657065,0.43895747599451296,0.27434842249657065,0.43895747599451296,0.7023319615912208,0.43895747599451296,0.27434842249657065,0.43895747599451296,0.27434842249657065,0.1714677640603567,0.27434842249657065,0.1714677640603567,0.27434842249657065,0.438957475994513,0.27434842249657065,0.1714677640603567,0.27434842249657065,0.1714677640603567,};
    double int_points[27][3]={
        {-0.7745966692414834,-0.7745966692414834,-0.7745966692414834,},
        {0.0,-0.7745966692414834,-0.7745966692414834,},
        {0.7745966692414834,-0.7745966692414834,-0.7745966692414834,},
        {-0.7745966692414834,0.0,-0.7745966692414834,},
        {0.0,0.0,-0.7745966692414834,},
        {0.7745966692414834,0.0,-0.7745966692414834,},
        {-0.7745966692414834,0.7745966692414834,-0.7745966692414834,},
        {0.0,0.7745966692414834,-0.7745966692414834,},
        {0.7745966692414834,0.7745966692414834,-0.7745966692414834,},
        {-0.7745966692414834,-0.7745966692414834,0.0,},
        {0.0,-0.7745966692414834,0.0,},
        {0.7745966692414834,-0.7745966692414834,0.0,},
        {-0.7745966692414834,0.0,0.0,},
        {0.0,0.0,0.0,},
        {0.7745966692414834,0.0,0.0,},
        {-0.7745966692414834,0.7745966692414834,0.0,},
        {0.0,0.7745966692414834,0.0,},
        {0.7745966692414834,0.7745966692414834,0.0,},
        {-0.7745966692414834,-0.7745966692414834,0.7745966692414834,},
        {0.0,-0.7745966692414834,0.7745966692414834,},
        {0.7745966692414834,-0.7745966692414834,0.7745966692414834,},
        {-0.7745966692414834,0.0,0.7745966692414834,},
        {0.0,0.0,0.7745966692414834,},
        {0.7745966692414834,0.0,0.7745966692414834,},
        {-0.7745966692414834,0.7745966692414834,0.7745966692414834,},
        {0.0,0.7745966692414834,0.7745966692414834,},
        {0.7745966692414834,0.7745966692414834,0.7745966692414834,}};

    #pragma omp parallel for
    for (size_t i=0;i<num_elem;i++){
        double TMP[20][3];
        for (size_t j=0;j<27;j++){
            C3D20_dynamic((double*) int_points[j],(double*) Coords[i],rho[i],(double*) Element_U[i],(double*) Element_V[i],(double*) Element_A[i],(double*) S[i][j],PENER[i][j],SENER[i][j],(double*) TMP);
            for (size_t num_node=0;num_node<20;num_node++){
                for (size_t ii=0;ii<3;ii++){
                    Conf_Force[i][num_node][ii]+=TMP[num_node][ii]*weights[j];
                }
            }
        }
    }
}

void C3D20_static_mbf(double *rst,double *coord,double *Element_U,double *S,double PENER,double SENER,double *res_0){
    double r = rst[0];
    double s = rst[1];
    double t = rst[2];
    double *U = Element_U;
    double S11 = S[0];
    double S22 = S[1];
    double S33 = S[2];
    double S12 = S[3];
    double S13 = S[4];
    double S23 = S[5];
    
    double x0 = (1.0/4.0)*r;
    double x1 = -x0;
    double x2 = s*s;
    double x3 = x0*x2 + x1;
    double x4 = (1.0/4.0)*x2;
    double x5 = x4 - 1.0/4.0;
    double x6 = x3 + x5;
    double x7 = (1.0/4.0)*s;
    double x8 = -x7;
    double x9 = r*r;
    double x10 = x7*x9 + x8;
    double x11 = (1.0/4.0)*x9;
    double x12 = x11 - 1.0/4.0;
    double x13 = x10 + x12;
    double x14 = 1.0/4.0 - x4;
    double x15 = x14 + x3;
    double x16 = -x15;
    double x17 = 1.0/4.0 - x11;
    double x18 = x10 + x17;
    double x19 = -x18;
    double x20 = -x13;
    double x21 = -x6;
    double x22 = (1.0/2.0)*t;
    double x23 = r*x22;
    double x24 = -x23;
    double x25 = s*x22;
    double x26 = s*x23;
    double x27 = -x25 + x26;
    double x28 = -x22 - x24 - x27;
    double x29 = -x22;
    double x30 = x23 + x27 + x29;
    double x31 = x25 + x26;
    double x32 = x24 + x29 + x31;
    double x33 = -x22 - x23 - x31;
    double x34 = (1.0/8.0)*r;
    double x35 = s*x34;
    double x36 = -x35;
    double x37 = (1.0/8.0)*x9;
    double x38 = s*x37;
    double x39 = t*x7;
    double x40 = t*x0;
    double x41 = s*x40;
    double x42 = -x41;
    double x43 = x39 + x42;
    double x44 = x36 + x38 + x43;
    double x45 = (1.0/4.0)*t;
    double x46 = -x45;
    double x47 = (1.0/8.0)*x2;
    double x48 = r*x47;
    double x49 = x46 + x48;
    double x50 = 1.0/8.0 - x37;
    double x51 = -x47;
    double x52 = x40 + x51;
    double x53 = x50 + x52;
    double x54 = -x44 - x49 - x53;
    double x55 = -x40;
    double x56 = -x48;
    double x57 = x55 + x56;
    double x58 = x37 - 1.0/8.0;
    double x59 = x45 + x47 + x58;
    double x60 = x44 + x57 + x59;
    double x61 = x35 + x38;
    double x62 = x39 + x41;
    double x63 = x50 + x62;
    double x64 = -x46 - x51 - x57 - x61 - x63;
    double x65 = x40 + x48 + x59 + x61 + x62;
    double x66 = -x39;
    double x67 = x45 + x66;
    double x68 = x41 + x55;
    double x69 = x51 + x68;
    double x70 = x36 + x38;
    double x71 = x48 + x50 + x67 + x69 + x70;
    double x72 = x41 + x66;
    double x73 = x47 + x58;
    double x74 = -x40 - x46 - x56 - x70 - x72 - x73;
    double x75 = x42 + x61;
    double x76 = x53 + x56 + x67 + x75;
    double x77 = -x49 - x55 - x66 - x73 - x75;
    double x78 = x13*coord[30] + x15*coord[45] + x16*coord[33] + x18*coord[36] + x19*coord[24] + x20*coord[42] + x21*coord[39] + x28*coord[48] + x30*coord[57] + x32*coord[51] + x33*coord[54] + x54*coord[12] + x6*coord[27] + x60*coord[21] + x64*coord[15] + x65*coord[18] + x71*coord[0] + x74*coord[9] + x76*coord[3] + x77*coord[6];
    double x79 = t*t;
    double x80 = x0*x79 + x1;
    double x81 = (1.0/4.0)*x79;
    double x82 = x81 - 1.0/4.0;
    double x83 = x80 + x82;
    double x84 = t*x11 + x46;
    double x85 = x12 + x84;
    double x86 = 1.0/4.0 - x81;
    double x87 = x80 + x86;
    double x88 = -x87;
    double x89 = x17 + x84;
    double x90 = -x89;
    double x91 = -x85;
    double x92 = -x83;
    double x93 = (1.0/2.0)*s;
    double x94 = r*x93;
    double x95 = -x94;
    double x96 = -x27 - x93 - x95;
    double x97 = -x93;
    double x98 = x27 + x94 + x97;
    double x99 = x31 + x95 + x97;
    double x100 = -x31 - x93 - x94;
    double x101 = t*x34;
    double x102 = -x101;
    double x103 = t*x37;
    double x104 = x102 + x103 + x43;
    double x105 = (1.0/8.0)*x79;
    double x106 = r*x105;
    double x107 = x106 + x8;
    double x108 = -x105;
    double x109 = s*x0;
    double x110 = x108 + x109;
    double x111 = x110 + x50;
    double x112 = -x104 - x107 - x111;
    double x113 = -x109;
    double x114 = -x106;
    double x115 = x113 + x114;
    double x116 = x105 + x58 + x7;
    double x117 = x104 + x115 + x116;
    double x118 = x101 + x103;
    double x119 = -x108 - x115 - x118 - x63 - x8;
    double x120 = x106 + x109 + x116 + x118 + x62;
    double x121 = x108 + x113;
    double x122 = x102 + x103 + x72;
    double x123 = x106 + x121 + x122 + x50 + x7;
    double x124 = x105 + x58;
    double x125 = -x109 - x114 - x122 - x124 - x8;
    double x126 = x118 + x42 + x66;
    double x127 = x111 + x114 + x126 + x7;
    double x128 = -x107 - x113 - x124 - x126;
    double x129 = x100*coord[41] + x112*coord[11] + x117*coord[23] + x119*coord[8] + x120*coord[20] + x123*coord[2] + x125*coord[14] + x127*coord[5] + x128*coord[17] + x83*coord[53] + x85*coord[38] + x87*coord[59] + x88*coord[50] + x89*coord[32] + x90*coord[26] + x91*coord[44] + x92*coord[56] + x96*coord[35] + x98*coord[47] + x99*coord[29];
    double x130 = x13*coord[32] + x15*coord[47] + x16*coord[35] + x18*coord[38] + x19*coord[26] + x20*coord[44] + x21*coord[41] + x28*coord[50] + x30*coord[59] + x32*coord[53] + x33*coord[56] + x54*coord[14] + x6*coord[29] + x60*coord[23] + x64*coord[17] + x65*coord[20] + x71*coord[2] + x74*coord[11] + x76*coord[5] + x77*coord[8];
    double x131 = x100*coord[39] + x112*coord[9] + x117*coord[21] + x119*coord[6] + x120*coord[18] + x123*coord[0] + x125*coord[12] + x127*coord[3] + x128*coord[15] + x83*coord[51] + x85*coord[36] + x87*coord[57] + x88*coord[48] + x89*coord[30] + x90*coord[24] + x91*coord[42] + x92*coord[54] + x96*coord[33] + x98*coord[45] + x99*coord[27];
    double x132 = x130*x131;
    double x133 = x129*x78 - x132;
    double x134 = x100*coord[40] + x112*coord[10] + x117*coord[22] + x119*coord[7] + x120*coord[19] + x123*coord[1] + x125*coord[13] + x127*coord[4] + x128*coord[16] + x83*coord[52] + x85*coord[37] + x87*coord[58] + x88*coord[49] + x89*coord[31] + x90*coord[25] + x91*coord[43] + x92*coord[55] + x96*coord[34] + x98*coord[46] + x99*coord[28];
    double x135 = x7*x79 + x8;
    double x136 = x135 + x82;
    double x137 = t*x4 + x46;
    double x138 = x137 + x5;
    double x139 = x135 + x86;
    double x140 = -x139;
    double x141 = x137 + x14;
    double x142 = -x141;
    double x143 = -x138;
    double x144 = -x136;
    double x145 = (1.0/2.0)*r;
    double x146 = x24 + x26;
    double x147 = -x145 - x146 - x95;
    double x148 = -x145;
    double x149 = x146 + x148 + x94;
    double x150 = x23 + x26;
    double x151 = x148 + x150 + x95;
    double x152 = -x145 - x150 - x94;
    double x153 = t*x47;
    double x154 = -x153;
    double x155 = (1.0/8.0)*s*t;
    double x156 = x154 + x155;
    double x157 = s*x105;
    double x158 = -x157;
    double x159 = x0 + x158;
    double x160 = x105 + x47 - 1.0/8.0;
    double x161 = x160 + x68;
    double x162 = x113 + x156 + x159 + x161;
    double x163 = x1 + 1.0/8.0;
    double x164 = x157 + x69;
    double x165 = -x110 - x156 - x163 - x164;
    double x166 = x121 + x153;
    double x167 = x41 + x52;
    double x168 = -x155 - x158 - x163 - x166 - x167;
    double x169 = x109 + x153;
    double x170 = x157 + x160 + x40 + x41;
    double x171 = x0 + x155 + x169 + x170;
    double x172 = -x155;
    double x173 = x172 + 1.0/8.0;
    double x174 = x0 + x164 + x166 + x173;
    double x175 = x1 + x172;
    double x176 = -x158 - x161 - x169 - x175;
    double x177 = -x113 - x154 - x170 - x175;
    double x178 = x110 + x154 + x159 + x167 + x173;
    double x179 = x136*coord[59] + x138*coord[47] + x139*coord[53] + x140*coord[50] + x141*coord[29] + x142*coord[35] + x143*coord[41] + x144*coord[56] + x147*coord[26] + x149*coord[38] + x151*coord[32] + x152*coord[44] + x162*coord[5] + x165*coord[17] + x168*coord[8] + x171*coord[20] + x174*coord[2] + x176*coord[14] + x177*coord[11] + x178*coord[23];
    double x180 = x179*x78;
    double x181 = x13*coord[31] + x15*coord[46] + x16*coord[34] + x18*coord[37] + x19*coord[25] + x20*coord[43] + x21*coord[40] + x28*coord[49] + x30*coord[58] + x32*coord[52] + x33*coord[55] + x54*coord[13] + x6*coord[28] + x60*coord[22] + x64*coord[16] + x65*coord[19] + x71*coord[1] + x74*coord[10] + x76*coord[4] + x77*coord[7];
    double x182 = x136*coord[57] + x138*coord[45] + x139*coord[51] + x140*coord[48] + x141*coord[27] + x142*coord[33] + x143*coord[39] + x144*coord[54] + x147*coord[24] + x149*coord[36] + x151*coord[30] + x152*coord[42] + x162*coord[3] + x165*coord[15] + x168*coord[6] + x171*coord[18] + x174*coord[0] + x176*coord[12] + x177*coord[9] + x178*coord[21];
    double x183 = x129*x182;
    double x184 = x136*coord[58] + x138*coord[46] + x139*coord[52] + x140*coord[49] + x141*coord[28] + x142*coord[34] + x143*coord[40] + x144*coord[55] + x147*coord[25] + x149*coord[37] + x151*coord[31] + x152*coord[43] + x162*coord[4] + x165*coord[16] + x168*coord[7] + x171*coord[19] + x174*coord[1] + x176*coord[13] + x177*coord[10] + x178*coord[22];
    double x185 = x129*x184*x78 + x130*x134*x182 + x131*x179*x181 - x132*x184 - x134*x180 - x181*x183;
    double x186 = 1.0/x185;
    double x187 = x186*(x136*U[59] + x138*U[47] + x139*U[53] + x140*U[50] + x141*U[29] + x142*U[35] + x143*U[41] + x144*U[56] + x147*U[26] + x149*U[38] + x151*U[32] + x152*U[44] + x162*U[5] + x165*U[17] + x168*U[8] + x171*U[20] + x174*U[2] + x176*U[14] + x177*U[11] + x178*U[23]);
    double x188 = x130*x182 - x180;
    double x189 = x186*(x100*U[41] + x112*U[11] + x117*U[23] + x119*U[8] + x120*U[20] + x123*U[2] + x125*U[14] + x127*U[5] + x128*U[17] + x83*U[53] + x85*U[38] + x87*U[59] + x88*U[50] + x89*U[32] + x90*U[26] + x91*U[44] + x92*U[56] + x96*U[35] + x98*U[47] + x99*U[29]);
    double x190 = x131*x179 - x183;
    double x191 = x186*(x13*U[32] + x15*U[47] + x16*U[35] + x18*U[38] + x19*U[26] + x20*U[44] + x21*U[41] + x28*U[50] + x30*U[59] + x32*U[53] + x33*U[56] + x54*U[14] + x6*U[29] + x60*U[23] + x64*U[17] + x65*U[20] + x71*U[2] + x74*U[11] + x76*U[5] + x77*U[8]);
    double x192 = x133*x187 + x188*x189 + x190*x191;
    double x193 = -x130*x184 + x179*x181;
    double x194 = x186*(x100*U[40] + x112*U[10] + x117*U[22] + x119*U[7] + x120*U[19] + x123*U[1] + x125*U[13] + x127*U[4] + x128*U[16] + x83*U[52] + x85*U[37] + x87*U[58] + x88*U[49] + x89*U[31] + x90*U[25] + x91*U[43] + x92*U[55] + x96*U[34] + x98*U[46] + x99*U[28]);
    double x195 = -x129*x181 + x130*x134;
    double x196 = x186*(x136*U[58] + x138*U[46] + x139*U[52] + x140*U[49] + x141*U[28] + x142*U[34] + x143*U[40] + x144*U[55] + x147*U[25] + x149*U[37] + x151*U[31] + x152*U[43] + x162*U[4] + x165*U[16] + x168*U[7] + x171*U[19] + x174*U[1] + x176*U[13] + x177*U[10] + x178*U[22]);
    double x197 = x129*x184 - x134*x179;
    double x198 = x186*(x13*U[31] + x15*U[46] + x16*U[34] + x18*U[37] + x19*U[25] + x20*U[43] + x21*U[40] + x28*U[49] + x30*U[58] + x32*U[52] + x33*U[55] + x54*U[13] + x6*U[28] + x60*U[22] + x64*U[16] + x65*U[19] + x71*U[1] + x74*U[10] + x76*U[4] + x77*U[7]);
    double x199 = x193*x194 + x195*x196 + x197*x198;
    double x200 = x187*x195 + x189*x193 + x191*x197;
    double x201 = x133*x196 + x188*x194 + x190*x198 + 1.0;
    double x202 = x192*x199 - x200*x201;
    double x203 = x186*(x136*U[57] + x138*U[45] + x139*U[51] + x140*U[48] + x141*U[27] + x142*U[33] + x143*U[39] + x144*U[54] + x147*U[24] + x149*U[36] + x151*U[30] + x152*U[42] + x162*U[3] + x165*U[15] + x168*U[6] + x171*U[18] + x174*U[0] + x176*U[12] + x177*U[9] + x178*U[21]);
    double x204 = x186*(x100*U[39] + x112*U[9] + x117*U[21] + x119*U[6] + x120*U[18] + x123*U[0] + x125*U[12] + x127*U[3] + x128*U[15] + x83*U[51] + x85*U[36] + x87*U[57] + x88*U[48] + x89*U[30] + x90*U[24] + x91*U[42] + x92*U[54] + x96*U[33] + x98*U[45] + x99*U[27]);
    double x205 = x186*(x13*U[30] + x15*U[45] + x16*U[33] + x18*U[36] + x19*U[24] + x20*U[42] + x21*U[39] + x28*U[48] + x30*U[57] + x32*U[51] + x33*U[54] + x54*U[12] + x6*U[27] + x60*U[21] + x64*U[15] + x65*U[18] + x71*U[0] + x74*U[9] + x76*U[3] + x77*U[6]);
    double x206 = x133*x203 + x188*x204 + x190*x205;
    double x207 = x193*x204 + x195*x203 + x197*x205 + 1.0;
    double x208 = -x192*x207 + x200*x206;
    double x209 = -x199*x206 + x201*x207;
    double x210 = S11*x202 + S12*x208 + S13*x209;
    double x211 = S12*x202 + S22*x208 + S23*x209;
    double x212 = S13*x202 + S23*x208 + S33*x209;
    double x213 = -x199*x211 - x200*x212 - x207*x210;
    double x214 = -x181*x182 + x184*x78;
    double x215 = x123*x186;
    double x216 = x131*x181 - x134*x78;
    double x217 = x174*x186;
    double x218 = -x131*x184 + x134*x182;
    double x219 = x186*x71;
    double x220 = x214*x215 + x216*x217 + x218*x219;
    double x221 = x194*x214 + x196*x216 + x198*x218;
    double x222 = x187*x216 + x189*x214 + x191*x218 + 1.0;
    double x223 = -x199*x222 + x200*x221;
    double x224 = x203*x216 + x204*x214 + x205*x218;
    double x225 = -x200*x224 + x207*x222;
    double x226 = x199*x224 - x207*x221;
    double x227 = S11*x223 + S12*x225 + S13*x226;
    double x228 = S12*x223 + S22*x225 + S23*x226;
    double x229 = S13*x223 + S23*x225 + S33*x226;
    double x230 = -x199*x228 - x200*x229 - x207*x227;
    double x231 = x133*x217 + x188*x215 + x190*x219;
    double x232 = x193*x215 + x195*x217 + x197*x219;
    double x233 = -x192*x221 + x201*x222;
    double x234 = x192*x224 - x206*x222;
    double x235 = -x201*x224 + x206*x221;
    double x236 = S11*x233 + S12*x234 + S13*x235;
    double x237 = S12*x233 + S22*x234 + S23*x235;
    double x238 = S13*x233 + S23*x234 + S33*x235;
    double x239 = -1.0*PENER - 1.0*SENER;
    double x240 = -x199*x237 - x200*x238 - x207*x236 - x239;
    double x241 = -x192*x212 - x201*x211 - x206*x210;
    double x242 = -x192*x238 - x201*x237 - x206*x236;
    double x243 = -x192*x229 - x201*x228 - x206*x227 - x239;
    double x244 = -x221*x228 - x222*x229 - x224*x227;
    double x245 = -x221*x237 - x222*x238 - x224*x236;
    double x246 = -x210*x224 - x211*x221 - x212*x222 - x239;
    double x247 = x127*x186;
    double x248 = x162*x186;
    double x249 = x186*x76;
    double x250 = x214*x247 + x216*x248 + x218*x249;
    double x251 = x133*x248 + x188*x247 + x190*x249;
    double x252 = x193*x247 + x195*x248 + x197*x249;
    double x253 = x119*x186;
    double x254 = x168*x186;
    double x255 = x186*x77;
    double x256 = x214*x253 + x216*x254 + x218*x255;
    double x257 = x133*x254 + x188*x253 + x190*x255;
    double x258 = x193*x253 + x195*x254 + x197*x255;
    double x259 = x112*x186;
    double x260 = x177*x186;
    double x261 = x186*x74;
    double x262 = x214*x259 + x216*x260 + x218*x261;
    double x263 = x133*x260 + x188*x259 + x190*x261;
    double x264 = x193*x259 + x195*x260 + x197*x261;
    double x265 = x125*x186;
    double x266 = x176*x186;
    double x267 = x186*x54;
    double x268 = x214*x265 + x216*x266 + x218*x267;
    double x269 = x133*x266 + x188*x265 + x190*x267;
    double x270 = x193*x265 + x195*x266 + x197*x267;
    double x271 = x128*x186;
    double x272 = x165*x186;
    double x273 = x186*x64;
    double x274 = x214*x271 + x216*x272 + x218*x273;
    double x275 = x133*x272 + x188*x271 + x190*x273;
    double x276 = x193*x271 + x195*x272 + x197*x273;
    double x277 = x120*x186;
    double x278 = x171*x186;
    double x279 = x186*x65;
    double x280 = x214*x277 + x216*x278 + x218*x279;
    double x281 = x133*x278 + x188*x277 + x190*x279;
    double x282 = x193*x277 + x195*x278 + x197*x279;
    double x283 = x117*x186;
    double x284 = x178*x186;
    double x285 = x186*x60;
    double x286 = x214*x283 + x216*x284 + x218*x285;
    double x287 = x133*x284 + x188*x283 + x190*x285;
    double x288 = x193*x283 + x195*x284 + x197*x285;
    double x289 = x186*x90;
    double x290 = x147*x186;
    double x291 = x186*x19;
    double x292 = x214*x289 + x216*x290 + x218*x291;
    double x293 = x133*x290 + x188*x289 + x190*x291;
    double x294 = x193*x289 + x195*x290 + x197*x291;
    double x295 = x186*x99;
    double x296 = x141*x186;
    double x297 = x186*x6;
    double x298 = x214*x295 + x216*x296 + x218*x297;
    double x299 = x133*x296 + x188*x295 + x190*x297;
    double x300 = x193*x295 + x195*x296 + x197*x297;
    double x301 = x186*x89;
    double x302 = x151*x186;
    double x303 = x13*x186;
    double x304 = x214*x301 + x216*x302 + x218*x303;
    double x305 = x133*x302 + x188*x301 + x190*x303;
    double x306 = x193*x301 + x195*x302 + x197*x303;
    double x307 = x186*x96;
    double x308 = x142*x186;
    double x309 = x16*x186;
    double x310 = x214*x307 + x216*x308 + x218*x309;
    double x311 = x133*x308 + x188*x307 + x190*x309;
    double x312 = x193*x307 + x195*x308 + x197*x309;
    double x313 = x186*x85;
    double x314 = x149*x186;
    double x315 = x18*x186;
    double x316 = x214*x313 + x216*x314 + x218*x315;
    double x317 = x133*x314 + x188*x313 + x190*x315;
    double x318 = x193*x313 + x195*x314 + x197*x315;
    double x319 = x100*x186;
    double x320 = x143*x186;
    double x321 = x186*x21;
    double x322 = x214*x319 + x216*x320 + x218*x321;
    double x323 = x133*x320 + x188*x319 + x190*x321;
    double x324 = x193*x319 + x195*x320 + x197*x321;
    double x325 = x186*x91;
    double x326 = x152*x186;
    double x327 = x186*x20;
    double x328 = x214*x325 + x216*x326 + x218*x327;
    double x329 = x133*x326 + x188*x325 + x190*x327;
    double x330 = x193*x325 + x195*x326 + x197*x327;
    double x331 = x186*x98;
    double x332 = x138*x186;
    double x333 = x15*x186;
    double x334 = x214*x331 + x216*x332 + x218*x333;
    double x335 = x133*x332 + x188*x331 + x190*x333;
    double x336 = x193*x331 + x195*x332 + x197*x333;
    double x337 = x186*x88;
    double x338 = x140*x186;
    double x339 = x186*x28;
    double x340 = x214*x337 + x216*x338 + x218*x339;
    double x341 = x133*x338 + x188*x337 + x190*x339;
    double x342 = x193*x337 + x195*x338 + x197*x339;
    double x343 = x186*x83;
    double x344 = x139*x186;
    double x345 = x186*x32;
    double x346 = x214*x343 + x216*x344 + x218*x345;
    double x347 = x133*x344 + x188*x343 + x190*x345;
    double x348 = x193*x343 + x195*x344 + x197*x345;
    double x349 = x186*x92;
    double x350 = x144*x186;
    double x351 = x186*x33;
    double x352 = x214*x349 + x216*x350 + x218*x351;
    double x353 = x133*x350 + x188*x349 + x190*x351;
    double x354 = x193*x349 + x195*x350 + x197*x351;
    double x355 = x186*x87;
    double x356 = x136*x186;
    double x357 = x186*x30;
    double x358 = x214*x355 + x216*x356 + x218*x357;
    double x359 = x133*x356 + x188*x355 + x190*x357;
    double x360 = x193*x355 + x195*x356 + x197*x357;
    
    res_0[0] = x185*(x213*x220 + x230*x231 + x232*x240);
    res_0[1] = x185*(x220*x241 + x231*x243 + x232*x242);
    res_0[2] = x185*(x220*x246 + x231*x244 + x232*x245);
    res_0[3] = x185*(x213*x250 + x230*x251 + x240*x252);
    res_0[4] = x185*(x241*x250 + x242*x252 + x243*x251);
    res_0[5] = x185*(x244*x251 + x245*x252 + x246*x250);
    res_0[6] = x185*(x213*x256 + x230*x257 + x240*x258);
    res_0[7] = x185*(x241*x256 + x242*x258 + x243*x257);
    res_0[8] = x185*(x244*x257 + x245*x258 + x246*x256);
    res_0[9] = x185*(x213*x262 + x230*x263 + x240*x264);
    res_0[10] = x185*(x241*x262 + x242*x264 + x243*x263);
    res_0[11] = x185*(x244*x263 + x245*x264 + x246*x262);
    res_0[12] = x185*(x213*x268 + x230*x269 + x240*x270);
    res_0[13] = x185*(x241*x268 + x242*x270 + x243*x269);
    res_0[14] = x185*(x244*x269 + x245*x270 + x246*x268);
    res_0[15] = x185*(x213*x274 + x230*x275 + x240*x276);
    res_0[16] = x185*(x241*x274 + x242*x276 + x243*x275);
    res_0[17] = x185*(x244*x275 + x245*x276 + x246*x274);
    res_0[18] = x185*(x213*x280 + x230*x281 + x240*x282);
    res_0[19] = x185*(x241*x280 + x242*x282 + x243*x281);
    res_0[20] = x185*(x244*x281 + x245*x282 + x246*x280);
    res_0[21] = x185*(x213*x286 + x230*x287 + x240*x288);
    res_0[22] = x185*(x241*x286 + x242*x288 + x243*x287);
    res_0[23] = x185*(x244*x287 + x245*x288 + x246*x286);
    res_0[24] = x185*(x213*x292 + x230*x293 + x240*x294);
    res_0[25] = x185*(x241*x292 + x242*x294 + x243*x293);
    res_0[26] = x185*(x244*x293 + x245*x294 + x246*x292);
    res_0[27] = x185*(x213*x298 + x230*x299 + x240*x300);
    res_0[28] = x185*(x241*x298 + x242*x300 + x243*x299);
    res_0[29] = x185*(x244*x299 + x245*x300 + x246*x298);
    res_0[30] = x185*(x213*x304 + x230*x305 + x240*x306);
    res_0[31] = x185*(x241*x304 + x242*x306 + x243*x305);
    res_0[32] = x185*(x244*x305 + x245*x306 + x246*x304);
    res_0[33] = x185*(x213*x310 + x230*x311 + x240*x312);
    res_0[34] = x185*(x241*x310 + x242*x312 + x243*x311);
    res_0[35] = x185*(x244*x311 + x245*x312 + x246*x310);
    res_0[36] = x185*(x213*x316 + x230*x317 + x240*x318);
    res_0[37] = x185*(x241*x316 + x242*x318 + x243*x317);
    res_0[38] = x185*(x244*x317 + x245*x318 + x246*x316);
    res_0[39] = x185*(x213*x322 + x230*x323 + x240*x324);
    res_0[40] = x185*(x241*x322 + x242*x324 + x243*x323);
    res_0[41] = x185*(x244*x323 + x245*x324 + x246*x322);
    res_0[42] = x185*(x213*x328 + x230*x329 + x240*x330);
    res_0[43] = x185*(x241*x328 + x242*x330 + x243*x329);
    res_0[44] = x185*(x244*x329 + x245*x330 + x246*x328);
    res_0[45] = x185*(x213*x334 + x230*x335 + x240*x336);
    res_0[46] = x185*(x241*x334 + x242*x336 + x243*x335);
    res_0[47] = x185*(x244*x335 + x245*x336 + x246*x334);
    res_0[48] = x185*(x213*x340 + x230*x341 + x240*x342);
    res_0[49] = x185*(x241*x340 + x242*x342 + x243*x341);
    res_0[50] = x185*(x244*x341 + x245*x342 + x246*x340);
    res_0[51] = x185*(x213*x346 + x230*x347 + x240*x348);
    res_0[52] = x185*(x241*x346 + x242*x348 + x243*x347);
    res_0[53] = x185*(x244*x347 + x245*x348 + x246*x346);
    res_0[54] = x185*(x213*x352 + x230*x353 + x240*x354);
    res_0[55] = x185*(x241*x352 + x242*x354 + x243*x353);
    res_0[56] = x185*(x244*x353 + x245*x354 + x246*x352);
    res_0[57] = x185*(x213*x358 + x230*x359 + x240*x360);
    res_0[58] = x185*(x241*x358 + x242*x360 + x243*x359);
    res_0[59] = x185*(x244*x359 + x245*x360 + x246*x358);
}

Conf_Forces_API void Integration_C3D20_static_mbf(size_t num_elem,double Coords[][20][3],double Element_U[][20][3],double S[][27][6],
    double PENER[][27],double SENER[][27],double Conf_Force[][20][3]){
    
    //Inputs:
    //    Coords.........Element_Nodal_Coordinates [num_elem][numNodes,3]
    //    Element_U......Element_Nodal_Displacements [num_elem][numNodes,3]
    //    S..............Stress at integration points [num_elem][num_int_points,6]
    //    PENER..........Plastic strain energy [num_elem][num_int_points]
    //    SENER..........Elastic strain energy [num_elem][num_int_points]
    //    
    double weights[27]={0.1714677640603567,0.27434842249657065,0.1714677640603567,0.27434842249657065,0.438957475994513,0.27434842249657065,0.1714677640603567,0.27434842249657065,0.1714677640603567,0.27434842249657065,0.43895747599451296,0.27434842249657065,0.43895747599451296,0.7023319615912208,0.43895747599451296,0.27434842249657065,0.43895747599451296,0.27434842249657065,0.1714677640603567,0.27434842249657065,0.1714677640603567,0.27434842249657065,0.438957475994513,0.27434842249657065,0.1714677640603567,0.27434842249657065,0.1714677640603567,};
    double int_points[27][3]={
        {-0.7745966692414834,-0.7745966692414834,-0.7745966692414834,},
        {0.0,-0.7745966692414834,-0.7745966692414834,},
        {0.7745966692414834,-0.7745966692414834,-0.7745966692414834,},
        {-0.7745966692414834,0.0,-0.7745966692414834,},
        {0.0,0.0,-0.7745966692414834,},
        {0.7745966692414834,0.0,-0.7745966692414834,},
        {-0.7745966692414834,0.7745966692414834,-0.7745966692414834,},
        {0.0,0.7745966692414834,-0.7745966692414834,},
        {0.7745966692414834,0.7745966692414834,-0.7745966692414834,},
        {-0.7745966692414834,-0.7745966692414834,0.0,},
        {0.0,-0.7745966692414834,0.0,},
        {0.7745966692414834,-0.7745966692414834,0.0,},
        {-0.7745966692414834,0.0,0.0,},
        {0.0,0.0,0.0,},
        {0.7745966692414834,0.0,0.0,},
        {-0.7745966692414834,0.7745966692414834,0.0,},
        {0.0,0.7745966692414834,0.0,},
        {0.7745966692414834,0.7745966692414834,0.0,},
        {-0.7745966692414834,-0.7745966692414834,0.7745966692414834,},
        {0.0,-0.7745966692414834,0.7745966692414834,},
        {0.7745966692414834,-0.7745966692414834,0.7745966692414834,},
        {-0.7745966692414834,0.0,0.7745966692414834,},
        {0.0,0.0,0.7745966692414834,},
        {0.7745966692414834,0.0,0.7745966692414834,},
        {-0.7745966692414834,0.7745966692414834,0.7745966692414834,},
        {0.0,0.7745966692414834,0.7745966692414834,},
        {0.7745966692414834,0.7745966692414834,0.7745966692414834,}};

    #pragma omp parallel for
    for (size_t i=0;i<num_elem;i++){
        double TMP[20][3];
        for (size_t j=0;j<27;j++){
            C3D20_static_mbf((double*) int_points[j],(double*) Coords[i],(double*) Element_U[i],(double*) S[i][j],PENER[i][j],SENER[i][j],(double*) TMP);
            for (size_t num_node=0;num_node<20;num_node++){
                for (size_t ii=0;ii<3;ii++){
                    Conf_Force[i][num_node][ii]+=TMP[num_node][ii]*weights[j];
                }
            }
        }
    }
}

void C3D20_static_dbf(double *rst,double *coord,double *Element_U,double *S,double PENER,double SENER,double *res_0){
    double r = rst[0];
    double s = rst[1];
    double t = rst[2];
    double *U = Element_U;
    double S11 = S[0];
    double S22 = S[1];
    double S33 = S[2];
    double S12 = S[3];
    double S13 = S[4];
    double S23 = S[5];
    
    double x0 = (1.0/4.0)*r;
    double x1 = -x0;
    double x2 = s*s;
    double x3 = x0*x2 + x1;
    double x4 = (1.0/4.0)*x2;
    double x5 = x4 - 1.0/4.0;
    double x6 = x3 + x5;
    double x7 = (1.0/4.0)*s;
    double x8 = -x7;
    double x9 = r*r;
    double x10 = x7*x9 + x8;
    double x11 = (1.0/4.0)*x9;
    double x12 = x11 - 1.0/4.0;
    double x13 = x10 + x12;
    double x14 = 1.0/4.0 - x4;
    double x15 = x14 + x3;
    double x16 = -x15;
    double x17 = 1.0/4.0 - x11;
    double x18 = x10 + x17;
    double x19 = -x18;
    double x20 = -x13;
    double x21 = -x6;
    double x22 = (1.0/2.0)*t;
    double x23 = r*x22;
    double x24 = -x23;
    double x25 = s*x22;
    double x26 = s*x23;
    double x27 = -x25 + x26;
    double x28 = -x22 - x24 - x27;
    double x29 = -x22;
    double x30 = x23 + x27 + x29;
    double x31 = x25 + x26;
    double x32 = x24 + x29 + x31;
    double x33 = -x22 - x23 - x31;
    double x34 = (1.0/8.0)*r;
    double x35 = s*x34;
    double x36 = -x35;
    double x37 = (1.0/8.0)*x9;
    double x38 = s*x37;
    double x39 = t*x7;
    double x40 = t*x0;
    double x41 = s*x40;
    double x42 = -x41;
    double x43 = x39 + x42;
    double x44 = x36 + x38 + x43;
    double x45 = (1.0/4.0)*t;
    double x46 = -x45;
    double x47 = (1.0/8.0)*x2;
    double x48 = r*x47;
    double x49 = x46 + x48;
    double x50 = 1.0/8.0 - x37;
    double x51 = -x47;
    double x52 = x40 + x51;
    double x53 = x50 + x52;
    double x54 = -x44 - x49 - x53;
    double x55 = -x40;
    double x56 = -x48;
    double x57 = x55 + x56;
    double x58 = x37 - 1.0/8.0;
    double x59 = x45 + x47 + x58;
    double x60 = x44 + x57 + x59;
    double x61 = x35 + x38;
    double x62 = x39 + x41;
    double x63 = x50 + x62;
    double x64 = -x46 - x51 - x57 - x61 - x63;
    double x65 = x40 + x48 + x59 + x61 + x62;
    double x66 = -x39;
    double x67 = x45 + x66;
    double x68 = x41 + x55;
    double x69 = x51 + x68;
    double x70 = x36 + x38;
    double x71 = x48 + x50 + x67 + x69 + x70;
    double x72 = x41 + x66;
    double x73 = x47 + x58;
    double x74 = -x40 - x46 - x56 - x70 - x72 - x73;
    double x75 = x42 + x61;
    double x76 = x53 + x56 + x67 + x75;
    double x77 = -x49 - x55 - x66 - x73 - x75;
    double x78 = x13*coord[30] + x15*coord[45] + x16*coord[33] + x18*coord[36] + x19*coord[24] + x20*coord[42] + x21*coord[39] + x28*coord[48] + x30*coord[57] + x32*coord[51] + x33*coord[54] + x54*coord[12] + x6*coord[27] + x60*coord[21] + x64*coord[15] + x65*coord[18] + x71*coord[0] + x74*coord[9] + x76*coord[3] + x77*coord[6];
    double x79 = t*t;
    double x80 = x0*x79 + x1;
    double x81 = (1.0/4.0)*x79;
    double x82 = x81 - 1.0/4.0;
    double x83 = x80 + x82;
    double x84 = t*x11 + x46;
    double x85 = x12 + x84;
    double x86 = 1.0/4.0 - x81;
    double x87 = x80 + x86;
    double x88 = -x87;
    double x89 = x17 + x84;
    double x90 = -x89;
    double x91 = -x85;
    double x92 = -x83;
    double x93 = (1.0/2.0)*s;
    double x94 = r*x93;
    double x95 = -x94;
    double x96 = -x27 - x93 - x95;
    double x97 = -x93;
    double x98 = x27 + x94 + x97;
    double x99 = x31 + x95 + x97;
    double x100 = -x31 - x93 - x94;
    double x101 = t*x34;
    double x102 = -x101;
    double x103 = t*x37;
    double x104 = x102 + x103 + x43;
    double x105 = (1.0/8.0)*x79;
    double x106 = r*x105;
    double x107 = x106 + x8;
    double x108 = -x105;
    double x109 = s*x0;
    double x110 = x108 + x109;
    double x111 = x110 + x50;
    double x112 = -x104 - x107 - x111;
    double x113 = -x109;
    double x114 = -x106;
    double x115 = x113 + x114;
    double x116 = x105 + x58 + x7;
    double x117 = x104 + x115 + x116;
    double x118 = x101 + x103;
    double x119 = -x108 - x115 - x118 - x63 - x8;
    double x120 = x106 + x109 + x116 + x118 + x62;
    double x121 = x108 + x113;
    double x122 = x102 + x103 + x72;
    double x123 = x106 + x121 + x122 + x50 + x7;
    double x124 = x105 + x58;
    double x125 = -x109 - x114 - x122 - x124 - x8;
    double x126 = x118 + x42 + x66;
    double x127 = x111 + x114 + x126 + x7;
    double x128 = -x107 - x113 - x124 - x126;
    double x129 = x100*coord[41] + x112*coord[11] + x117*coord[23] + x119*coord[8] + x120*coord[20] + x123*coord[2] + x125*coord[14] + x127*coord[5] + x128*coord[17] + x83*coord[53] + x85*coord[38] + x87*coord[59] + x88*coord[50] + x89*coord[32] + x90*coord[26] + x91*coord[44] + x92*coord[56] + x96*coord[35] + x98*coord[47] + x99*coord[29];
    double x130 = x13*coord[32] + x15*coord[47] + x16*coord[35] + x18*coord[38] + x19*coord[26] + x20*coord[44] + x21*coord[41] + x28*coord[50] + x30*coord[59] + x32*coord[53] + x33*coord[56] + x54*coord[14] + x6*coord[29] + x60*coord[23] + x64*coord[17] + x65*coord[20] + x71*coord[2] + x74*coord[11] + x76*coord[5] + x77*coord[8];
    double x131 = x100*coord[39] + x112*coord[9] + x117*coord[21] + x119*coord[6] + x120*coord[18] + x123*coord[0] + x125*coord[12] + x127*coord[3] + x128*coord[15] + x83*coord[51] + x85*coord[36] + x87*coord[57] + x88*coord[48] + x89*coord[30] + x90*coord[24] + x91*coord[42] + x92*coord[54] + x96*coord[33] + x98*coord[45] + x99*coord[27];
    double x132 = x130*x131;
    double x133 = x129*x78 - x132;
    double x134 = x100*coord[40] + x112*coord[10] + x117*coord[22] + x119*coord[7] + x120*coord[19] + x123*coord[1] + x125*coord[13] + x127*coord[4] + x128*coord[16] + x83*coord[52] + x85*coord[37] + x87*coord[58] + x88*coord[49] + x89*coord[31] + x90*coord[25] + x91*coord[43] + x92*coord[55] + x96*coord[34] + x98*coord[46] + x99*coord[28];
    double x135 = x7*x79 + x8;
    double x136 = x135 + x82;
    double x137 = t*x4 + x46;
    double x138 = x137 + x5;
    double x139 = x135 + x86;
    double x140 = -x139;
    double x141 = x137 + x14;
    double x142 = -x141;
    double x143 = -x138;
    double x144 = -x136;
    double x145 = (1.0/2.0)*r;
    double x146 = x24 + x26;
    double x147 = -x145 - x146 - x95;
    double x148 = -x145;
    double x149 = x146 + x148 + x94;
    double x150 = x23 + x26;
    double x151 = x148 + x150 + x95;
    double x152 = -x145 - x150 - x94;
    double x153 = t*x47;
    double x154 = -x153;
    double x155 = (1.0/8.0)*s*t;
    double x156 = x154 + x155;
    double x157 = s*x105;
    double x158 = -x157;
    double x159 = x0 + x158;
    double x160 = x105 + x47 - 1.0/8.0;
    double x161 = x160 + x68;
    double x162 = x113 + x156 + x159 + x161;
    double x163 = x1 + 1.0/8.0;
    double x164 = x157 + x69;
    double x165 = -x110 - x156 - x163 - x164;
    double x166 = x121 + x153;
    double x167 = x41 + x52;
    double x168 = -x155 - x158 - x163 - x166 - x167;
    double x169 = x109 + x153;
    double x170 = x157 + x160 + x40 + x41;
    double x171 = x0 + x155 + x169 + x170;
    double x172 = -x155;
    double x173 = x172 + 1.0/8.0;
    double x174 = x0 + x164 + x166 + x173;
    double x175 = x1 + x172;
    double x176 = -x158 - x161 - x169 - x175;
    double x177 = -x113 - x154 - x170 - x175;
    double x178 = x110 + x154 + x159 + x167 + x173;
    double x179 = x136*coord[59] + x138*coord[47] + x139*coord[53] + x140*coord[50] + x141*coord[29] + x142*coord[35] + x143*coord[41] + x144*coord[56] + x147*coord[26] + x149*coord[38] + x151*coord[32] + x152*coord[44] + x162*coord[5] + x165*coord[17] + x168*coord[8] + x171*coord[20] + x174*coord[2] + x176*coord[14] + x177*coord[11] + x178*coord[23];
    double x180 = x179*x78;
    double x181 = x13*coord[31] + x15*coord[46] + x16*coord[34] + x18*coord[37] + x19*coord[25] + x20*coord[43] + x21*coord[40] + x28*coord[49] + x30*coord[58] + x32*coord[52] + x33*coord[55] + x54*coord[13] + x6*coord[28] + x60*coord[22] + x64*coord[16] + x65*coord[19] + x71*coord[1] + x74*coord[10] + x76*coord[4] + x77*coord[7];
    double x182 = x136*coord[57] + x138*coord[45] + x139*coord[51] + x140*coord[48] + x141*coord[27] + x142*coord[33] + x143*coord[39] + x144*coord[54] + x147*coord[24] + x149*coord[36] + x151*coord[30] + x152*coord[42] + x162*coord[3] + x165*coord[15] + x168*coord[6] + x171*coord[18] + x174*coord[0] + x176*coord[12] + x177*coord[9] + x178*coord[21];
    double x183 = x129*x182;
    double x184 = x136*coord[58] + x138*coord[46] + x139*coord[52] + x140*coord[49] + x141*coord[28] + x142*coord[34] + x143*coord[40] + x144*coord[55] + x147*coord[25] + x149*coord[37] + x151*coord[31] + x152*coord[43] + x162*coord[4] + x165*coord[16] + x168*coord[7] + x171*coord[19] + x174*coord[1] + x176*coord[13] + x177*coord[10] + x178*coord[22];
    double x185 = x129*x184*x78 + x130*x134*x182 + x131*x179*x181 - x132*x184 - x134*x180 - x181*x183;
    double x186 = 1.0/x185;
    double x187 = x186*(x136*U[59] + x138*U[47] + x139*U[53] + x140*U[50] + x141*U[29] + x142*U[35] + x143*U[41] + x144*U[56] + x147*U[26] + x149*U[38] + x151*U[32] + x152*U[44] + x162*U[5] + x165*U[17] + x168*U[8] + x171*U[20] + x174*U[2] + x176*U[14] + x177*U[11] + x178*U[23]);
    double x188 = x130*x182 - x180;
    double x189 = x186*(x100*U[41] + x112*U[11] + x117*U[23] + x119*U[8] + x120*U[20] + x123*U[2] + x125*U[14] + x127*U[5] + x128*U[17] + x83*U[53] + x85*U[38] + x87*U[59] + x88*U[50] + x89*U[32] + x90*U[26] + x91*U[44] + x92*U[56] + x96*U[35] + x98*U[47] + x99*U[29]);
    double x190 = x131*x179 - x183;
    double x191 = x186*(x13*U[32] + x15*U[47] + x16*U[35] + x18*U[38] + x19*U[26] + x20*U[44] + x21*U[41] + x28*U[50] + x30*U[59] + x32*U[53] + x33*U[56] + x54*U[14] + x6*U[29] + x60*U[23] + x64*U[17] + x65*U[20] + x71*U[2] + x74*U[11] + x76*U[5] + x77*U[8]);
    double x192 = x133*x187 + x188*x189 + x190*x191;
    double x193 = -x130*x184 + x179*x181;
    double x194 = x186*(x100*U[40] + x112*U[10] + x117*U[22] + x119*U[7] + x120*U[19] + x123*U[1] + x125*U[13] + x127*U[4] + x128*U[16] + x83*U[52] + x85*U[37] + x87*U[58] + x88*U[49] + x89*U[31] + x90*U[25] + x91*U[43] + x92*U[55] + x96*U[34] + x98*U[46] + x99*U[28]);
    double x195 = -x129*x181 + x130*x134;
    double x196 = x186*(x136*U[58] + x138*U[46] + x139*U[52] + x140*U[49] + x141*U[28] + x142*U[34] + x143*U[40] + x144*U[55] + x147*U[25] + x149*U[37] + x151*U[31] + x152*U[43] + x162*U[4] + x165*U[16] + x168*U[7] + x171*U[19] + x174*U[1] + x176*U[13] + x177*U[10] + x178*U[22]);
    double x197 = x129*x184 - x134*x179;
    double x198 = x186*(x13*U[31] + x15*U[46] + x16*U[34] + x18*U[37] + x19*U[25] + x20*U[43] + x21*U[40] + x28*U[49] + x30*U[58] + x32*U[52] + x33*U[55] + x54*U[13] + x6*U[28] + x60*U[22] + x64*U[16] + x65*U[19] + x71*U[1] + x74*U[10] + x76*U[4] + x77*U[7]);
    double x199 = x193*x194 + x195*x196 + x197*x198;
    double x200 = x187*x195 + x189*x193 + x191*x197;
    double x201 = x133*x196 + x188*x194 + x190*x198;
    double x202 = x201 + 1.0;
    double x203 = x192*x199 - x200*x202;
    double x204 = x186*(x136*U[57] + x138*U[45] + x139*U[51] + x140*U[48] + x141*U[27] + x142*U[33] + x143*U[39] + x144*U[54] + x147*U[24] + x149*U[36] + x151*U[30] + x152*U[42] + x162*U[3] + x165*U[15] + x168*U[6] + x171*U[18] + x174*U[0] + x176*U[12] + x177*U[9] + x178*U[21]);
    double x205 = x186*(x100*U[39] + x112*U[9] + x117*U[21] + x119*U[6] + x120*U[18] + x123*U[0] + x125*U[12] + x127*U[3] + x128*U[15] + x83*U[51] + x85*U[36] + x87*U[57] + x88*U[48] + x89*U[30] + x90*U[24] + x91*U[42] + x92*U[54] + x96*U[33] + x98*U[45] + x99*U[27]);
    double x206 = x186*(x13*U[30] + x15*U[45] + x16*U[33] + x18*U[36] + x19*U[24] + x20*U[42] + x21*U[39] + x28*U[48] + x30*U[57] + x32*U[51] + x33*U[54] + x54*U[12] + x6*U[27] + x60*U[21] + x64*U[15] + x65*U[18] + x71*U[0] + x74*U[9] + x76*U[3] + x77*U[6]);
    double x207 = x133*x204 + x188*x205 + x190*x206;
    double x208 = x193*x205 + x195*x204 + x197*x206;
    double x209 = x208 + 1.0;
    double x210 = -x192*x209 + x200*x207;
    double x211 = -x199*x207 + x202*x209;
    double x212 = S11*x203 + S12*x210 + S13*x211;
    double x213 = S12*x203 + S22*x210 + S23*x211;
    double x214 = S13*x203 + S23*x210 + S33*x211;
    double x215 = -x199*x213 - x200*x214 - x208*x212;
    double x216 = -x181*x182 + x184*x78;
    double x217 = x123*x186;
    double x218 = x131*x181 - x134*x78;
    double x219 = x174*x186;
    double x220 = -x131*x184 + x134*x182;
    double x221 = x186*x71;
    double x222 = x216*x217 + x218*x219 + x220*x221;
    double x223 = x194*x216 + x196*x218 + x198*x220;
    double x224 = x187*x218 + x189*x216 + x191*x220;
    double x225 = x224 + 1.0;
    double x226 = -x199*x225 + x200*x223;
    double x227 = x204*x218 + x205*x216 + x206*x220;
    double x228 = -x200*x227 + x209*x225;
    double x229 = x199*x227 - x209*x223;
    double x230 = S11*x226 + S12*x228 + S13*x229;
    double x231 = S12*x226 + S22*x228 + S23*x229;
    double x232 = S13*x226 + S23*x228 + S33*x229;
    double x233 = -x199*x231 - x200*x232 - x208*x230;
    double x234 = x133*x219 + x188*x217 + x190*x221;
    double x235 = x193*x217 + x195*x219 + x197*x221;
    double x236 = -x192*x223 + x202*x225;
    double x237 = x192*x227 - x207*x225;
    double x238 = -x202*x227 + x207*x223;
    double x239 = S11*x236 + S12*x237 + S13*x238;
    double x240 = S12*x236 + S22*x237 + S23*x238;
    double x241 = S13*x236 + S23*x237 + S33*x238;
    double x242 = -1.0*PENER - 1.0*SENER;
    double x243 = -x199*x240 - x200*x241 - x208*x239 - x242;
    double x244 = -x192*x214 - x201*x213 - x207*x212;
    double x245 = -x192*x241 - x201*x240 - x207*x239;
    double x246 = -x192*x232 - x201*x231 - x207*x230 - x242;
    double x247 = -x223*x231 - x224*x232 - x227*x230;
    double x248 = -x223*x240 - x224*x241 - x227*x239;
    double x249 = -x212*x227 - x213*x223 - x214*x224 - x242;
    double x250 = x127*x186;
    double x251 = x162*x186;
    double x252 = x186*x76;
    double x253 = x216*x250 + x218*x251 + x220*x252;
    double x254 = x133*x251 + x188*x250 + x190*x252;
    double x255 = x193*x250 + x195*x251 + x197*x252;
    double x256 = x119*x186;
    double x257 = x168*x186;
    double x258 = x186*x77;
    double x259 = x216*x256 + x218*x257 + x220*x258;
    double x260 = x133*x257 + x188*x256 + x190*x258;
    double x261 = x193*x256 + x195*x257 + x197*x258;
    double x262 = x112*x186;
    double x263 = x177*x186;
    double x264 = x186*x74;
    double x265 = x216*x262 + x218*x263 + x220*x264;
    double x266 = x133*x263 + x188*x262 + x190*x264;
    double x267 = x193*x262 + x195*x263 + x197*x264;
    double x268 = x125*x186;
    double x269 = x176*x186;
    double x270 = x186*x54;
    double x271 = x216*x268 + x218*x269 + x220*x270;
    double x272 = x133*x269 + x188*x268 + x190*x270;
    double x273 = x193*x268 + x195*x269 + x197*x270;
    double x274 = x128*x186;
    double x275 = x165*x186;
    double x276 = x186*x64;
    double x277 = x216*x274 + x218*x275 + x220*x276;
    double x278 = x133*x275 + x188*x274 + x190*x276;
    double x279 = x193*x274 + x195*x275 + x197*x276;
    double x280 = x120*x186;
    double x281 = x171*x186;
    double x282 = x186*x65;
    double x283 = x216*x280 + x218*x281 + x220*x282;
    double x284 = x133*x281 + x188*x280 + x190*x282;
    double x285 = x193*x280 + x195*x281 + x197*x282;
    double x286 = x117*x186;
    double x287 = x178*x186;
    double x288 = x186*x60;
    double x289 = x216*x286 + x218*x287 + x220*x288;
    double x290 = x133*x287 + x188*x286 + x190*x288;
    double x291 = x193*x286 + x195*x287 + x197*x288;
    double x292 = x186*x90;
    double x293 = x147*x186;
    double x294 = x186*x19;
    double x295 = x216*x292 + x218*x293 + x220*x294;
    double x296 = x133*x293 + x188*x292 + x190*x294;
    double x297 = x193*x292 + x195*x293 + x197*x294;
    double x298 = x186*x99;
    double x299 = x141*x186;
    double x300 = x186*x6;
    double x301 = x216*x298 + x218*x299 + x220*x300;
    double x302 = x133*x299 + x188*x298 + x190*x300;
    double x303 = x193*x298 + x195*x299 + x197*x300;
    double x304 = x186*x89;
    double x305 = x151*x186;
    double x306 = x13*x186;
    double x307 = x216*x304 + x218*x305 + x220*x306;
    double x308 = x133*x305 + x188*x304 + x190*x306;
    double x309 = x193*x304 + x195*x305 + x197*x306;
    double x310 = x186*x96;
    double x311 = x142*x186;
    double x312 = x16*x186;
    double x313 = x216*x310 + x218*x311 + x220*x312;
    double x314 = x133*x311 + x188*x310 + x190*x312;
    double x315 = x193*x310 + x195*x311 + x197*x312;
    double x316 = x186*x85;
    double x317 = x149*x186;
    double x318 = x18*x186;
    double x319 = x216*x316 + x218*x317 + x220*x318;
    double x320 = x133*x317 + x188*x316 + x190*x318;
    double x321 = x193*x316 + x195*x317 + x197*x318;
    double x322 = x100*x186;
    double x323 = x143*x186;
    double x324 = x186*x21;
    double x325 = x216*x322 + x218*x323 + x220*x324;
    double x326 = x133*x323 + x188*x322 + x190*x324;
    double x327 = x193*x322 + x195*x323 + x197*x324;
    double x328 = x186*x91;
    double x329 = x152*x186;
    double x330 = x186*x20;
    double x331 = x216*x328 + x218*x329 + x220*x330;
    double x332 = x133*x329 + x188*x328 + x190*x330;
    double x333 = x193*x328 + x195*x329 + x197*x330;
    double x334 = x186*x98;
    double x335 = x138*x186;
    double x336 = x15*x186;
    double x337 = x216*x334 + x218*x335 + x220*x336;
    double x338 = x133*x335 + x188*x334 + x190*x336;
    double x339 = x193*x334 + x195*x335 + x197*x336;
    double x340 = x186*x88;
    double x341 = x140*x186;
    double x342 = x186*x28;
    double x343 = x216*x340 + x218*x341 + x220*x342;
    double x344 = x133*x341 + x188*x340 + x190*x342;
    double x345 = x193*x340 + x195*x341 + x197*x342;
    double x346 = x186*x83;
    double x347 = x139*x186;
    double x348 = x186*x32;
    double x349 = x216*x346 + x218*x347 + x220*x348;
    double x350 = x133*x347 + x188*x346 + x190*x348;
    double x351 = x193*x346 + x195*x347 + x197*x348;
    double x352 = x186*x92;
    double x353 = x144*x186;
    double x354 = x186*x33;
    double x355 = x216*x352 + x218*x353 + x220*x354;
    double x356 = x133*x353 + x188*x352 + x190*x354;
    double x357 = x193*x352 + x195*x353 + x197*x354;
    double x358 = x186*x87;
    double x359 = x136*x186;
    double x360 = x186*x30;
    double x361 = x216*x358 + x218*x359 + x220*x360;
    double x362 = x133*x359 + x188*x358 + x190*x360;
    double x363 = x193*x358 + x195*x359 + x197*x360;
    
    res_0[0] = x185*(x215*x222 + x233*x234 + x235*x243);
    res_0[1] = x185*(x222*x244 + x234*x246 + x235*x245);
    res_0[2] = x185*(x222*x249 + x234*x247 + x235*x248);
    res_0[3] = x185*(x215*x253 + x233*x254 + x243*x255);
    res_0[4] = x185*(x244*x253 + x245*x255 + x246*x254);
    res_0[5] = x185*(x247*x254 + x248*x255 + x249*x253);
    res_0[6] = x185*(x215*x259 + x233*x260 + x243*x261);
    res_0[7] = x185*(x244*x259 + x245*x261 + x246*x260);
    res_0[8] = x185*(x247*x260 + x248*x261 + x249*x259);
    res_0[9] = x185*(x215*x265 + x233*x266 + x243*x267);
    res_0[10] = x185*(x244*x265 + x245*x267 + x246*x266);
    res_0[11] = x185*(x247*x266 + x248*x267 + x249*x265);
    res_0[12] = x185*(x215*x271 + x233*x272 + x243*x273);
    res_0[13] = x185*(x244*x271 + x245*x273 + x246*x272);
    res_0[14] = x185*(x247*x272 + x248*x273 + x249*x271);
    res_0[15] = x185*(x215*x277 + x233*x278 + x243*x279);
    res_0[16] = x185*(x244*x277 + x245*x279 + x246*x278);
    res_0[17] = x185*(x247*x278 + x248*x279 + x249*x277);
    res_0[18] = x185*(x215*x283 + x233*x284 + x243*x285);
    res_0[19] = x185*(x244*x283 + x245*x285 + x246*x284);
    res_0[20] = x185*(x247*x284 + x248*x285 + x249*x283);
    res_0[21] = x185*(x215*x289 + x233*x290 + x243*x291);
    res_0[22] = x185*(x244*x289 + x245*x291 + x246*x290);
    res_0[23] = x185*(x247*x290 + x248*x291 + x249*x289);
    res_0[24] = x185*(x215*x295 + x233*x296 + x243*x297);
    res_0[25] = x185*(x244*x295 + x245*x297 + x246*x296);
    res_0[26] = x185*(x247*x296 + x248*x297 + x249*x295);
    res_0[27] = x185*(x215*x301 + x233*x302 + x243*x303);
    res_0[28] = x185*(x244*x301 + x245*x303 + x246*x302);
    res_0[29] = x185*(x247*x302 + x248*x303 + x249*x301);
    res_0[30] = x185*(x215*x307 + x233*x308 + x243*x309);
    res_0[31] = x185*(x244*x307 + x245*x309 + x246*x308);
    res_0[32] = x185*(x247*x308 + x248*x309 + x249*x307);
    res_0[33] = x185*(x215*x313 + x233*x314 + x243*x315);
    res_0[34] = x185*(x244*x313 + x245*x315 + x246*x314);
    res_0[35] = x185*(x247*x314 + x248*x315 + x249*x313);
    res_0[36] = x185*(x215*x319 + x233*x320 + x243*x321);
    res_0[37] = x185*(x244*x319 + x245*x321 + x246*x320);
    res_0[38] = x185*(x247*x320 + x248*x321 + x249*x319);
    res_0[39] = x185*(x215*x325 + x233*x326 + x243*x327);
    res_0[40] = x185*(x244*x325 + x245*x327 + x246*x326);
    res_0[41] = x185*(x247*x326 + x248*x327 + x249*x325);
    res_0[42] = x185*(x215*x331 + x233*x332 + x243*x333);
    res_0[43] = x185*(x244*x331 + x245*x333 + x246*x332);
    res_0[44] = x185*(x247*x332 + x248*x333 + x249*x331);
    res_0[45] = x185*(x215*x337 + x233*x338 + x243*x339);
    res_0[46] = x185*(x244*x337 + x245*x339 + x246*x338);
    res_0[47] = x185*(x247*x338 + x248*x339 + x249*x337);
    res_0[48] = x185*(x215*x343 + x233*x344 + x243*x345);
    res_0[49] = x185*(x244*x343 + x245*x345 + x246*x344);
    res_0[50] = x185*(x247*x344 + x248*x345 + x249*x343);
    res_0[51] = x185*(x215*x349 + x233*x350 + x243*x351);
    res_0[52] = x185*(x244*x349 + x245*x351 + x246*x350);
    res_0[53] = x185*(x247*x350 + x248*x351 + x249*x349);
    res_0[54] = x185*(x215*x355 + x233*x356 + x243*x357);
    res_0[55] = x185*(x244*x355 + x245*x357 + x246*x356);
    res_0[56] = x185*(x247*x356 + x248*x357 + x249*x355);
    res_0[57] = x185*(x215*x361 + x233*x362 + x243*x363);
    res_0[58] = x185*(x244*x361 + x245*x363 + x246*x362);
    res_0[59] = x185*(x247*x362 + x248*x363 + x249*x361);
}

Conf_Forces_API void Integration_C3D20_static_dbf(size_t num_elem,double Coords[][20][3],double Element_U[][20][3],double S[][27][6],
    double PENER[][27],double SENER[][27],double Conf_Force[][20][3]){
    
    //Inputs:
    //    Coords.........Element_Nodal_Coordinates [num_elem][numNodes,3]
    //    Element_U......Element_Nodal_Displacements [num_elem][numNodes,3]
    //    S..............Stress at integration points [num_elem][num_int_points,6]
    //    PENER..........Plastic strain energy [num_elem][num_int_points]
    //    SENER..........Elastic strain energy [num_elem][num_int_points]
    //    
    double weights[27]={0.1714677640603567,0.27434842249657065,0.1714677640603567,0.27434842249657065,0.438957475994513,0.27434842249657065,0.1714677640603567,0.27434842249657065,0.1714677640603567,0.27434842249657065,0.43895747599451296,0.27434842249657065,0.43895747599451296,0.7023319615912208,0.43895747599451296,0.27434842249657065,0.43895747599451296,0.27434842249657065,0.1714677640603567,0.27434842249657065,0.1714677640603567,0.27434842249657065,0.438957475994513,0.27434842249657065,0.1714677640603567,0.27434842249657065,0.1714677640603567,};
    double int_points[27][3]={
        {-0.7745966692414834,-0.7745966692414834,-0.7745966692414834,},
        {0.0,-0.7745966692414834,-0.7745966692414834,},
        {0.7745966692414834,-0.7745966692414834,-0.7745966692414834,},
        {-0.7745966692414834,0.0,-0.7745966692414834,},
        {0.0,0.0,-0.7745966692414834,},
        {0.7745966692414834,0.0,-0.7745966692414834,},
        {-0.7745966692414834,0.7745966692414834,-0.7745966692414834,},
        {0.0,0.7745966692414834,-0.7745966692414834,},
        {0.7745966692414834,0.7745966692414834,-0.7745966692414834,},
        {-0.7745966692414834,-0.7745966692414834,0.0,},
        {0.0,-0.7745966692414834,0.0,},
        {0.7745966692414834,-0.7745966692414834,0.0,},
        {-0.7745966692414834,0.0,0.0,},
        {0.0,0.0,0.0,},
        {0.7745966692414834,0.0,0.0,},
        {-0.7745966692414834,0.7745966692414834,0.0,},
        {0.0,0.7745966692414834,0.0,},
        {0.7745966692414834,0.7745966692414834,0.0,},
        {-0.7745966692414834,-0.7745966692414834,0.7745966692414834,},
        {0.0,-0.7745966692414834,0.7745966692414834,},
        {0.7745966692414834,-0.7745966692414834,0.7745966692414834,},
        {-0.7745966692414834,0.0,0.7745966692414834,},
        {0.0,0.0,0.7745966692414834,},
        {0.7745966692414834,0.0,0.7745966692414834,},
        {-0.7745966692414834,0.7745966692414834,0.7745966692414834,},
        {0.0,0.7745966692414834,0.7745966692414834,},
        {0.7745966692414834,0.7745966692414834,0.7745966692414834,}};

    #pragma omp parallel for
    for (size_t i=0;i<num_elem;i++){
        double TMP[20][3];
        for (size_t j=0;j<27;j++){
            C3D20_static_dbf((double*) int_points[j],(double*) Coords[i],(double*) Element_U[i],(double*) S[i][j],PENER[i][j],SENER[i][j],(double*) TMP);
            for (size_t num_node=0;num_node<20;num_node++){
                for (size_t ii=0;ii<3;ii++){
                    Conf_Force[i][num_node][ii]+=TMP[num_node][ii]*weights[j];
                }
            }
        }
    }
}

void C3D20R_dynamic(double *rst,double *coord,double rho,double *Element_U,double *Element_V,double *Element_A,double *S,double PENER,double SENER,double *res_0){
    double r = rst[0];
    double s = rst[1];
    double t = rst[2];
    double *U = Element_U;
    double *V = Element_V;
    double *A = Element_A;
    double S11 = S[0];
    double S22 = S[1];
    double S33 = S[2];
    double S12 = S[3];
    double S13 = S[4];
    double S23 = S[5];
    
    double x0 = (1.0/4.0)*r;
    double x1 = -x0;
    double x2 = s*s;
    double x3 = x0*x2;
    double x4 = x1 + x3;
    double x5 = (1.0/4.0)*x2;
    double x6 = x5 - 1.0/4.0;
    double x7 = x4 + x6;
    double x8 = (1.0/4.0)*s;
    double x9 = -x8;
    double x10 = r*r;
    double x11 = x10*x8;
    double x12 = x11 + x9;
    double x13 = (1.0/4.0)*x10;
    double x14 = x13 - 1.0/4.0;
    double x15 = x12 + x14;
    double x16 = 1.0/4.0 - x5;
    double x17 = x16 + x4;
    double x18 = -x17;
    double x19 = 1.0/4.0 - x13;
    double x20 = x12 + x19;
    double x21 = -x20;
    double x22 = -x15;
    double x23 = -x7;
    double x24 = (1.0/2.0)*t;
    double x25 = r*x24;
    double x26 = -x25;
    double x27 = s*x24;
    double x28 = s*x25;
    double x29 = -x27 + x28;
    double x30 = -x24 - x26 - x29;
    double x31 = -x24;
    double x32 = x25 + x29 + x31;
    double x33 = x27 + x28;
    double x34 = x26 + x31 + x33;
    double x35 = -x24 - x25 - x33;
    double x36 = (1.0/8.0)*r;
    double x37 = s*x36;
    double x38 = -x37;
    double x39 = t*x8;
    double x40 = t*x0;
    double x41 = s*x40;
    double x42 = -x41;
    double x43 = x39 + x42;
    double x44 = x38 + x43;
    double x45 = (1.0/4.0)*t;
    double x46 = -x45;
    double x47 = (1.0/8.0)*x2;
    double x48 = r*x47;
    double x49 = (1.0/8.0)*x10;
    double x50 = s*x49;
    double x51 = x48 + x50;
    double x52 = x46 + x51;
    double x53 = 1.0/8.0 - x49;
    double x54 = -x47;
    double x55 = x40 + x54;
    double x56 = x53 + x55;
    double x57 = -x44 - x52 - x56;
    double x58 = -x40;
    double x59 = x45 + x58;
    double x60 = -x48;
    double x61 = x50 + x60;
    double x62 = x49 - 1.0/8.0;
    double x63 = x47 + x62;
    double x64 = x61 + x63;
    double x65 = x44 + x59 + x64;
    double x66 = x41 + x46;
    double x67 = x37 + x61;
    double x68 = x53 + x54;
    double x69 = -x39 - x58 - x66 - x67 - x68;
    double x70 = x41 + x51;
    double x71 = x37 + x63;
    double x72 = x39 + x40 + x45 + x70 + x71;
    double x73 = -x39;
    double x74 = x38 + x73;
    double x75 = x59 + x68 + x70 + x74;
    double x76 = -x40 - x64 - x66 - x74;
    double x77 = x42 + x73;
    double x78 = x45 + x56 + x67 + x77;
    double x79 = -x52 - x58 - x71 - x77;
    double x80 = x15*coord[30] + x17*coord[45] + x18*coord[33] + x20*coord[36] + x21*coord[24] + x22*coord[42] + x23*coord[39] + x30*coord[48] + x32*coord[57] + x34*coord[51] + x35*coord[54] + x57*coord[12] + x65*coord[21] + x69*coord[15] + x7*coord[27] + x72*coord[18] + x75*coord[0] + x76*coord[9] + x78*coord[3] + x79*coord[6];
    double x81 = t*t;
    double x82 = x0*x81;
    double x83 = x1 + x82;
    double x84 = (1.0/4.0)*x81;
    double x85 = x84 - 1.0/4.0;
    double x86 = x83 + x85;
    double x87 = t*x13;
    double x88 = x46 + x87;
    double x89 = x14 + x88;
    double x90 = 1.0/4.0 - x84;
    double x91 = x83 + x90;
    double x92 = -x91;
    double x93 = x19 + x88;
    double x94 = -x93;
    double x95 = -x89;
    double x96 = -x86;
    double x97 = (1.0/2.0)*s;
    double x98 = r*x97;
    double x99 = -x98;
    double x100 = -x29 - x97 - x99;
    double x101 = -x97;
    double x102 = x101 + x29 + x98;
    double x103 = x101 + x33 + x99;
    double x104 = -x33 - x97 - x98;
    double x105 = t*x36;
    double x106 = -x105;
    double x107 = x106 + x43;
    double x108 = (1.0/8.0)*x81;
    double x109 = r*x108;
    double x110 = t*x49;
    double x111 = x109 + x110;
    double x112 = x111 + x9;
    double x113 = -x108;
    double x114 = s*x0;
    double x115 = x113 + x114;
    double x116 = x115 + x53;
    double x117 = -x107 - x112 - x116;
    double x118 = -x114;
    double x119 = x118 + x8;
    double x120 = -x109;
    double x121 = x108 + x120;
    double x122 = x110 + x121;
    double x123 = x122 + x62;
    double x124 = x107 + x119 + x123;
    double x125 = x41 + x9;
    double x126 = x105 + x118;
    double x127 = x113 + x53;
    double x128 = x110 + x120;
    double x129 = -x125 - x126 - x127 - x128 - x39;
    double x130 = x111 + x41;
    double x131 = x105 + x8;
    double x132 = x108 + x62;
    double x133 = x114 + x130 + x131 + x132 + x39;
    double x134 = x106 + x73;
    double x135 = x119 + x127 + x130 + x134;
    double x136 = -x114 - x123 - x125 - x134;
    double x137 = x116 + x128 + x131 + x77;
    double x138 = -x112 - x126 - x132 - x77;
    double x139 = x100*coord[35] + x102*coord[47] + x103*coord[29] + x104*coord[41] + x117*coord[11] + x124*coord[23] + x129*coord[8] + x133*coord[20] + x135*coord[2] + x136*coord[14] + x137*coord[5] + x138*coord[17] + x86*coord[53] + x89*coord[38] + x91*coord[59] + x92*coord[50] + x93*coord[32] + x94*coord[26] + x95*coord[44] + x96*coord[56];
    double x140 = x15*coord[32] + x17*coord[47] + x18*coord[35] + x20*coord[38] + x21*coord[26] + x22*coord[44] + x23*coord[41] + x30*coord[50] + x32*coord[59] + x34*coord[53] + x35*coord[56] + x57*coord[14] + x65*coord[23] + x69*coord[17] + x7*coord[29] + x72*coord[20] + x75*coord[2] + x76*coord[11] + x78*coord[5] + x79*coord[8];
    double x141 = x100*coord[33] + x102*coord[45] + x103*coord[27] + x104*coord[39] + x117*coord[9] + x124*coord[21] + x129*coord[6] + x133*coord[18] + x135*coord[0] + x136*coord[12] + x137*coord[3] + x138*coord[15] + x86*coord[51] + x89*coord[36] + x91*coord[57] + x92*coord[48] + x93*coord[30] + x94*coord[24] + x95*coord[42] + x96*coord[54];
    double x142 = x140*x141;
    double x143 = x139*x80 - x142;
    double x144 = x100*coord[34] + x102*coord[46] + x103*coord[28] + x104*coord[40] + x117*coord[10] + x124*coord[22] + x129*coord[7] + x133*coord[19] + x135*coord[1] + x136*coord[13] + x137*coord[4] + x138*coord[16] + x86*coord[52] + x89*coord[37] + x91*coord[58] + x92*coord[49] + x93*coord[31] + x94*coord[25] + x95*coord[43] + x96*coord[55];
    double x145 = x8*x81;
    double x146 = x145 + x9;
    double x147 = x146 + x85;
    double x148 = t*x5;
    double x149 = x148 + x46;
    double x150 = x149 + x6;
    double x151 = x146 + x90;
    double x152 = -x151;
    double x153 = x149 + x16;
    double x154 = -x153;
    double x155 = -x150;
    double x156 = -x147;
    double x157 = (1.0/2.0)*r;
    double x158 = x26 + x28;
    double x159 = -x157 - x158 - x99;
    double x160 = -x157;
    double x161 = x158 + x160 + x98;
    double x162 = x25 + x28;
    double x163 = x160 + x162 + x99;
    double x164 = -x157 - x162 - x98;
    double x165 = t*x47;
    double x166 = -x165;
    double x167 = x108 + x47;
    double x168 = x166 + x167;
    double x169 = (1.0/8.0)*s;
    double x170 = t*x169;
    double x171 = x170 + x41 + x58;
    double x172 = s*x108;
    double x173 = -x172;
    double x174 = x118 + x173;
    double x175 = x0 - 1.0/8.0;
    double x176 = x168 + x171 + x174 + x175;
    double x177 = x1 + 1.0/8.0;
    double x178 = x166 + x172;
    double x179 = -x115 - x171 - x177 - x178 - x54;
    double x180 = x165 + x41;
    double x181 = x170 + x180;
    double x182 = -x113 - x174 - x177 - x181 - x55;
    double x183 = x114 + x167;
    double x184 = x172 + x40;
    double x185 = x175 + x181 + x183 + x184;
    double x186 = -x170;
    double x187 = x118 + x186;
    double x188 = x180 + x58;
    double x189 = x0 + 1.0/8.0;
    double x190 = x113 + x172 + x187 + x188 + x189 + x54;
    double x191 = x1 - 1.0/8.0;
    double x192 = x173 + x186;
    double x193 = -x183 - x188 - x191 - x192;
    double x194 = -x168 - x184 - x187 - x191 - x41;
    double x195 = x115 + x166 + x189 + x192 + x41 + x55;
    double x196 = x147*coord[59] + x150*coord[47] + x151*coord[53] + x152*coord[50] + x153*coord[29] + x154*coord[35] + x155*coord[41] + x156*coord[56] + x159*coord[26] + x161*coord[38] + x163*coord[32] + x164*coord[44] + x176*coord[5] + x179*coord[17] + x182*coord[8] + x185*coord[20] + x190*coord[2] + x193*coord[14] + x194*coord[11] + x195*coord[23];
    double x197 = x196*x80;
    double x198 = x15*coord[31] + x17*coord[46] + x18*coord[34] + x20*coord[37] + x21*coord[25] + x22*coord[43] + x23*coord[40] + x30*coord[49] + x32*coord[58] + x34*coord[52] + x35*coord[55] + x57*coord[13] + x65*coord[22] + x69*coord[16] + x7*coord[28] + x72*coord[19] + x75*coord[1] + x76*coord[10] + x78*coord[4] + x79*coord[7];
    double x199 = x147*coord[57] + x150*coord[45] + x151*coord[51] + x152*coord[48] + x153*coord[27] + x154*coord[33] + x155*coord[39] + x156*coord[54] + x159*coord[24] + x161*coord[36] + x163*coord[30] + x164*coord[42] + x176*coord[3] + x179*coord[15] + x182*coord[6] + x185*coord[18] + x190*coord[0] + x193*coord[12] + x194*coord[9] + x195*coord[21];
    double x200 = x139*x199;
    double x201 = x147*coord[58] + x150*coord[46] + x151*coord[52] + x152*coord[49] + x153*coord[28] + x154*coord[34] + x155*coord[40] + x156*coord[55] + x159*coord[25] + x161*coord[37] + x163*coord[31] + x164*coord[43] + x176*coord[4] + x179*coord[16] + x182*coord[7] + x185*coord[19] + x190*coord[1] + x193*coord[13] + x194*coord[10] + x195*coord[22];
    double x202 = x139*x201*x80 + x140*x144*x199 + x141*x196*x198 - x142*x201 - x144*x197 - x198*x200;
    double x203 = 1.0/x202;
    double x204 = x203*(x147*U[59] + x150*U[47] + x151*U[53] + x152*U[50] + x153*U[29] + x154*U[35] + x155*U[41] + x156*U[56] + x159*U[26] + x161*U[38] + x163*U[32] + x164*U[44] + x176*U[5] + x179*U[17] + x182*U[8] + x185*U[20] + x190*U[2] + x193*U[14] + x194*U[11] + x195*U[23]);
    double x205 = x140*x199 - x197;
    double x206 = x203*(x100*U[35] + x102*U[47] + x103*U[29] + x104*U[41] + x117*U[11] + x124*U[23] + x129*U[8] + x133*U[20] + x135*U[2] + x136*U[14] + x137*U[5] + x138*U[17] + x86*U[53] + x89*U[38] + x91*U[59] + x92*U[50] + x93*U[32] + x94*U[26] + x95*U[44] + x96*U[56]);
    double x207 = x141*x196 - x200;
    double x208 = x203*(x15*U[32] + x17*U[47] + x18*U[35] + x20*U[38] + x21*U[26] + x22*U[44] + x23*U[41] + x30*U[50] + x32*U[59] + x34*U[53] + x35*U[56] + x57*U[14] + x65*U[23] + x69*U[17] + x7*U[29] + x72*U[20] + x75*U[2] + x76*U[11] + x78*U[5] + x79*U[8]);
    double x209 = x143*x204 + x205*x206 + x207*x208;
    double x210 = -x140*x201 + x196*x198;
    double x211 = x203*(x100*U[34] + x102*U[46] + x103*U[28] + x104*U[40] + x117*U[10] + x124*U[22] + x129*U[7] + x133*U[19] + x135*U[1] + x136*U[13] + x137*U[4] + x138*U[16] + x86*U[52] + x89*U[37] + x91*U[58] + x92*U[49] + x93*U[31] + x94*U[25] + x95*U[43] + x96*U[55]);
    double x212 = -x139*x198 + x140*x144;
    double x213 = x203*(x147*U[58] + x150*U[46] + x151*U[52] + x152*U[49] + x153*U[28] + x154*U[34] + x155*U[40] + x156*U[55] + x159*U[25] + x161*U[37] + x163*U[31] + x164*U[43] + x176*U[4] + x179*U[16] + x182*U[7] + x185*U[19] + x190*U[1] + x193*U[13] + x194*U[10] + x195*U[22]);
    double x214 = x139*x201 - x144*x196;
    double x215 = x203*(x15*U[31] + x17*U[46] + x18*U[34] + x20*U[37] + x21*U[25] + x22*U[43] + x23*U[40] + x30*U[49] + x32*U[58] + x34*U[52] + x35*U[55] + x57*U[13] + x65*U[22] + x69*U[16] + x7*U[28] + x72*U[19] + x75*U[1] + x76*U[10] + x78*U[4] + x79*U[7]);
    double x216 = x210*x211 + x212*x213 + x214*x215;
    double x217 = x204*x212 + x206*x210 + x208*x214;
    double x218 = x143*x213 + x205*x211 + x207*x215;
    double x219 = x218 + 1.0;
    double x220 = x209*x216 - x217*x219;
    double x221 = x203*(x147*U[57] + x150*U[45] + x151*U[51] + x152*U[48] + x153*U[27] + x154*U[33] + x155*U[39] + x156*U[54] + x159*U[24] + x161*U[36] + x163*U[30] + x164*U[42] + x176*U[3] + x179*U[15] + x182*U[6] + x185*U[18] + x190*U[0] + x193*U[12] + x194*U[9] + x195*U[21]);
    double x222 = x203*(x100*U[33] + x102*U[45] + x103*U[27] + x104*U[39] + x117*U[9] + x124*U[21] + x129*U[6] + x133*U[18] + x135*U[0] + x136*U[12] + x137*U[3] + x138*U[15] + x86*U[51] + x89*U[36] + x91*U[57] + x92*U[48] + x93*U[30] + x94*U[24] + x95*U[42] + x96*U[54]);
    double x223 = x203*(x15*U[30] + x17*U[45] + x18*U[33] + x20*U[36] + x21*U[24] + x22*U[42] + x23*U[39] + x30*U[48] + x32*U[57] + x34*U[51] + x35*U[54] + x57*U[12] + x65*U[21] + x69*U[15] + x7*U[27] + x72*U[18] + x75*U[0] + x76*U[9] + x78*U[3] + x79*U[6]);
    double x224 = x143*x221 + x205*x222 + x207*x223;
    double x225 = x210*x222 + x212*x221 + x214*x223;
    double x226 = x225 + 1.0;
    double x227 = -x209*x226 + x217*x224;
    double x228 = -x216*x224 + x219*x226;
    double x229 = S11*x220 + S12*x227 + S13*x228;
    double x230 = S12*x220 + S22*x227 + S23*x228;
    double x231 = S13*x220 + S23*x227 + S33*x228;
    double x232 = -x216*x230 - x217*x231 - x225*x229;
    double x233 = -x198*x199 + x201*x80;
    double x234 = x135*x203;
    double x235 = x141*x198 - x144*x80;
    double x236 = x190*x203;
    double x237 = -x141*x201 + x144*x199;
    double x238 = x203*x75;
    double x239 = x233*x234 + x235*x236 + x237*x238;
    double x240 = x211*x233 + x213*x235 + x215*x237;
    double x241 = x204*x235 + x206*x233 + x208*x237;
    double x242 = x241 + 1.0;
    double x243 = -x216*x242 + x217*x240;
    double x244 = x221*x235 + x222*x233 + x223*x237;
    double x245 = -x217*x244 + x226*x242;
    double x246 = x216*x244 - x226*x240;
    double x247 = S11*x243 + S12*x245 + S13*x246;
    double x248 = S12*x243 + S22*x245 + S23*x246;
    double x249 = S13*x243 + S23*x245 + S33*x246;
    double x250 = -x216*x248 - x217*x249 - x225*x247;
    double x251 = x143*x236 + x205*x234 + x207*x238;
    double x252 = x210*x234 + x212*x236 + x214*x238;
    double x253 = -x209*x240 + x219*x242;
    double x254 = x209*x244 - x224*x242;
    double x255 = -x219*x244 + x224*x240;
    double x256 = S11*x253 + S12*x254 + S13*x255;
    double x257 = S12*x253 + S22*x254 + S23*x255;
    double x258 = S13*x253 + S23*x254 + S33*x255;
    double x259 = s*x82 + x118;
    double x260 = -x145 + x259 + x8;
    double x261 = x0 - x82;
    double x262 = -x260 - x261 - x85;
    double x263 = t*x3 + x58;
    double x264 = -x148 + x263 + x45;
    double x265 = x0 - x3;
    double x266 = -x264 - x265 - x6;
    double x267 = x17 + x264;
    double x268 = x260 + x91;
    double x269 = t*x11 + x73;
    double x270 = x269 + x45 - x87;
    double x271 = -x11 + x8;
    double x272 = -x14 - x270 - x271;
    double x273 = x20 + x270;
    double x274 = x151 + x259 + x261;
    double x275 = x269 + x271 + x93;
    double x276 = x153 + x263 + x265;
    double x277 = -x15 - x269 - x88;
    double x278 = -x149 - x263 - x7;
    double x279 = -x146 - x259 - x86;
    double x280 = -x36;
    double x281 = t*x37;
    double x282 = t*x50;
    double x283 = x280 + x281 + x282 + x49 - 1.0/4.0;
    double x284 = -x169;
    double x285 = t*x48;
    double x286 = x284 + x285;
    double x287 = (1.0/8.0)*t;
    double x288 = -x287;
    double x289 = s*x109;
    double x290 = x165 + x288 + x289;
    double x291 = x172 + x51;
    double x292 = x111 + x167;
    double x293 = x283 + x286 + x290 + x291 + x292;
    double x294 = -x285;
    double x295 = -x282;
    double x296 = x294 + x295 + x49 - 1.0/4.0;
    double x297 = -x281;
    double x298 = x280 + x297;
    double x299 = -x110;
    double x300 = x109 + x168 + x287 + x299;
    double x301 = x284 + x289 + x291 + x296 + x298 + x300;
    double x302 = -x289;
    double x303 = x295 + x302 + x49 - 1.0/4.0;
    double x304 = x169 + x173 - x50;
    double x305 = x304 + x48;
    double x306 = x165 + x288;
    double x307 = x285 + x292 + x298 + x303 + x305 + x306;
    double x308 = x294 + x302;
    double x309 = x283 + x300 + x305 + x308;
    double x310 = x36 + x47 + x61;
    double x311 = x282 + x297 + x49 - 1.0/4.0;
    double x312 = x122 + x172 + x284 + x306 + x308 + x310 + x311;
    double x313 = x121 + x287 + x299;
    double x314 = x178 + x281 + x286 + x303 + x310 + x313;
    double x315 = x304 + x36 + x47 + x60;
    double x316 = x122 + x281 + x290 + x296 + x315;
    double x317 = x166 + x285 + x289 + x311 + x313 + x315;
    double x318 = x262*V[48] + x266*V[33] + x267*V[45] + x268*V[57] + x272*V[24] + x273*V[36] + x274*V[51] + x275*V[30] + x276*V[27] + x277*V[42] + x278*V[39] + x279*V[54] + x293*V[18] + x301*V[6] + x307*V[15] + x309*V[3] + x312*V[21] + x314*V[9] + x316*V[12] + x317*V[0];
    double x319 = x262*V[49] + x266*V[34] + x267*V[46] + x268*V[58] + x272*V[25] + x273*V[37] + x274*V[52] + x275*V[31] + x276*V[28] + x277*V[43] + x278*V[40] + x279*V[55] + x293*V[19] + x301*V[7] + x307*V[16] + x309*V[4] + x312*V[22] + x314*V[10] + x316*V[13] + x317*V[1];
    double x320 = x262*V[50] + x266*V[35] + x267*V[47] + x268*V[59] + x272*V[26] + x273*V[38] + x274*V[53] + x275*V[32] + x276*V[29] + x277*V[44] + x278*V[41] + x279*V[56] + x293*V[20] + x301*V[8] + x307*V[17] + x309*V[5] + x312*V[23] + x314*V[11] + x316*V[14] + x317*V[2];
    double x321 = -1.0*PENER - 1.0*SENER + 0.5*rho*(x318*x318 + x319*x319 + x320*x320);
    double x322 = -x216*x257 - x217*x258 - x225*x256 - x321;
    double x323 = x262*A[48] + x266*A[33] + x267*A[45] + x268*A[57] + x272*A[24] + x273*A[36] + x274*A[51] + x275*A[30] + x276*A[27] + x277*A[42] + x278*A[39] + x279*A[54] + x293*A[18] + x301*A[6] + x307*A[15] + x309*A[3] + x312*A[21] + x314*A[9] + x316*A[12] + x317*A[0];
    double x324 = x262*A[49] + x266*A[34] + x267*A[46] + x268*A[58] + x272*A[25] + x273*A[37] + x274*A[52] + x275*A[31] + x276*A[28] + x277*A[43] + x278*A[40] + x279*A[55] + x293*A[19] + x301*A[7] + x307*A[16] + x309*A[4] + x312*A[22] + x314*A[10] + x316*A[13] + x317*A[1];
    double x325 = x262*A[50] + x266*A[35] + x267*A[47] + x268*A[59] + x272*A[26] + x273*A[38] + x274*A[53] + x275*A[32] + x276*A[29] + x277*A[44] + x278*A[41] + x279*A[56] + x293*A[20] + x301*A[8] + x307*A[17] + x309*A[5] + x312*A[23] + x314*A[11] + x316*A[14] + x317*A[2];
    double x326 = x100*V[33] + x102*V[45] + x103*V[27] + x104*V[39] + x117*V[9] + x124*V[21] + x129*V[6] + x133*V[18] + x135*V[0] + x136*V[12] + x137*V[3] + x138*V[15] + x86*V[51] + x89*V[36] + x91*V[57] + x92*V[48] + x93*V[30] + x94*V[24] + x95*V[42] + x96*V[54];
    double x327 = x203*x210;
    double x328 = x147*V[57] + x150*V[45] + x151*V[51] + x152*V[48] + x153*V[27] + x154*V[33] + x155*V[39] + x156*V[54] + x159*V[24] + x161*V[36] + x163*V[30] + x164*V[42] + x176*V[3] + x179*V[15] + x182*V[6] + x185*V[18] + x190*V[0] + x193*V[12] + x194*V[9] + x195*V[21];
    double x329 = x203*x212;
    double x330 = x15*V[30] + x17*V[45] + x18*V[33] + x20*V[36] + x21*V[24] + x22*V[42] + x23*V[39] + x30*V[48] + x32*V[57] + x34*V[51] + x35*V[54] + x57*V[12] + x65*V[21] + x69*V[15] + x7*V[27] + x72*V[18] + x75*V[0] + x76*V[9] + x78*V[3] + x79*V[6];
    double x331 = x203*x214;
    double x332 = x100*V[34] + x102*V[46] + x103*V[28] + x104*V[40] + x117*V[10] + x124*V[22] + x129*V[7] + x133*V[19] + x135*V[1] + x136*V[13] + x137*V[4] + x138*V[16] + x86*V[52] + x89*V[37] + x91*V[58] + x92*V[49] + x93*V[31] + x94*V[25] + x95*V[43] + x96*V[55];
    double x333 = x147*V[58] + x150*V[46] + x151*V[52] + x152*V[49] + x153*V[28] + x154*V[34] + x155*V[40] + x156*V[55] + x159*V[25] + x161*V[37] + x163*V[31] + x164*V[43] + x176*V[4] + x179*V[16] + x182*V[7] + x185*V[19] + x190*V[1] + x193*V[13] + x194*V[10] + x195*V[22];
    double x334 = x15*V[31] + x17*V[46] + x18*V[34] + x20*V[37] + x21*V[25] + x22*V[43] + x23*V[40] + x30*V[49] + x32*V[58] + x34*V[52] + x35*V[55] + x57*V[13] + x65*V[22] + x69*V[16] + x7*V[28] + x72*V[19] + x75*V[1] + x76*V[10] + x78*V[4] + x79*V[7];
    double x335 = x100*V[35] + x102*V[47] + x103*V[29] + x104*V[41] + x117*V[11] + x124*V[23] + x129*V[8] + x133*V[20] + x135*V[2] + x136*V[14] + x137*V[5] + x138*V[17] + x86*V[53] + x89*V[38] + x91*V[59] + x92*V[50] + x93*V[32] + x94*V[26] + x95*V[44] + x96*V[56];
    double x336 = x147*V[59] + x150*V[47] + x151*V[53] + x152*V[50] + x153*V[29] + x154*V[35] + x155*V[41] + x156*V[56] + x159*V[26] + x161*V[38] + x163*V[32] + x164*V[44] + x176*V[5] + x179*V[17] + x182*V[8] + x185*V[20] + x190*V[2] + x193*V[14] + x194*V[11] + x195*V[23];
    double x337 = x15*V[32] + x17*V[47] + x18*V[35] + x20*V[38] + x21*V[26] + x22*V[44] + x23*V[41] + x30*V[50] + x32*V[59] + x34*V[53] + x35*V[56] + x57*V[14] + x65*V[23] + x69*V[17] + x7*V[29] + x72*V[20] + x75*V[2] + x76*V[11] + x78*V[5] + x79*V[8];
    double x338 = x216*x324 + x217*x325 + x225*x323 + x318*(x326*x327 + x328*x329 + x330*x331) + x319*(x327*x332 + x329*x333 + x331*x334) + x320*(x327*x335 + x329*x336 + x331*x337);
    double x339 = rho*x317;
    double x340 = -x209*x231 - x218*x230 - x224*x229;
    double x341 = -x209*x258 - x218*x257 - x224*x256;
    double x342 = -x209*x249 - x218*x248 - x224*x247 - x321;
    double x343 = x143*x203;
    double x344 = x203*x205;
    double x345 = x203*x207;
    double x346 = x209*x325 + x218*x324 + x224*x323 + x318*(x326*x344 + x328*x343 + x330*x345) + x319*(x332*x344 + x333*x343 + x334*x345) + x320*(x335*x344 + x336*x343 + x337*x345);
    double x347 = -x240*x248 - x241*x249 - x244*x247;
    double x348 = -x240*x257 - x241*x258 - x244*x256;
    double x349 = -x229*x244 - x230*x240 - x231*x241 - x321;
    double x350 = x203*x233;
    double x351 = x203*x235;
    double x352 = x203*x237;
    double x353 = x240*x324 + x241*x325 + x244*x323 + x318*(x326*x350 + x328*x351 + x330*x352) + x319*(x332*x350 + x333*x351 + x334*x352) + x320*(x335*x350 + x336*x351 + x337*x352);
    double x354 = x137*x350 + x176*x351 + x352*x78;
    double x355 = x137*x344 + x176*x343 + x345*x78;
    double x356 = x137*x327 + x176*x329 + x331*x78;
    double x357 = rho*x309;
    double x358 = x129*x350 + x182*x351 + x352*x79;
    double x359 = x129*x344 + x182*x343 + x345*x79;
    double x360 = x129*x327 + x182*x329 + x331*x79;
    double x361 = rho*x301;
    double x362 = x117*x350 + x194*x351 + x352*x76;
    double x363 = x117*x344 + x194*x343 + x345*x76;
    double x364 = x117*x327 + x194*x329 + x331*x76;
    double x365 = rho*x314;
    double x366 = x136*x350 + x193*x351 + x352*x57;
    double x367 = x136*x344 + x193*x343 + x345*x57;
    double x368 = x136*x327 + x193*x329 + x331*x57;
    double x369 = rho*x316;
    double x370 = x138*x350 + x179*x351 + x352*x69;
    double x371 = x138*x344 + x179*x343 + x345*x69;
    double x372 = x138*x327 + x179*x329 + x331*x69;
    double x373 = rho*x307;
    double x374 = x133*x350 + x185*x351 + x352*x72;
    double x375 = x133*x344 + x185*x343 + x345*x72;
    double x376 = x133*x327 + x185*x329 + x331*x72;
    double x377 = rho*x293;
    double x378 = x124*x350 + x195*x351 + x352*x65;
    double x379 = x124*x344 + x195*x343 + x345*x65;
    double x380 = x124*x327 + x195*x329 + x331*x65;
    double x381 = rho*x312;
    double x382 = x159*x351 + x21*x352 + x350*x94;
    double x383 = x159*x343 + x21*x345 + x344*x94;
    double x384 = x159*x329 + x21*x331 + x327*x94;
    double x385 = rho*x272;
    double x386 = x103*x350 + x153*x351 + x352*x7;
    double x387 = x103*x344 + x153*x343 + x345*x7;
    double x388 = x103*x327 + x153*x329 + x331*x7;
    double x389 = rho*x276;
    double x390 = x15*x352 + x163*x351 + x350*x93;
    double x391 = x15*x345 + x163*x343 + x344*x93;
    double x392 = x15*x331 + x163*x329 + x327*x93;
    double x393 = rho*x275;
    double x394 = x100*x350 + x154*x351 + x18*x352;
    double x395 = x100*x344 + x154*x343 + x18*x345;
    double x396 = x100*x327 + x154*x329 + x18*x331;
    double x397 = rho*x266;
    double x398 = x161*x351 + x20*x352 + x350*x89;
    double x399 = x161*x343 + x20*x345 + x344*x89;
    double x400 = x161*x329 + x20*x331 + x327*x89;
    double x401 = rho*x273;
    double x402 = x104*x350 + x155*x351 + x23*x352;
    double x403 = x104*x344 + x155*x343 + x23*x345;
    double x404 = x104*x327 + x155*x329 + x23*x331;
    double x405 = rho*x278;
    double x406 = x164*x351 + x22*x352 + x350*x95;
    double x407 = x164*x343 + x22*x345 + x344*x95;
    double x408 = x164*x329 + x22*x331 + x327*x95;
    double x409 = rho*x277;
    double x410 = x102*x350 + x150*x351 + x17*x352;
    double x411 = x102*x344 + x150*x343 + x17*x345;
    double x412 = x102*x327 + x150*x329 + x17*x331;
    double x413 = rho*x267;
    double x414 = x152*x351 + x30*x352 + x350*x92;
    double x415 = x152*x343 + x30*x345 + x344*x92;
    double x416 = x152*x329 + x30*x331 + x327*x92;
    double x417 = rho*x262;
    double x418 = x151*x351 + x34*x352 + x350*x86;
    double x419 = x151*x343 + x34*x345 + x344*x86;
    double x420 = x151*x329 + x327*x86 + x331*x34;
    double x421 = rho*x274;
    double x422 = x156*x351 + x35*x352 + x350*x96;
    double x423 = x156*x343 + x344*x96 + x345*x35;
    double x424 = x156*x329 + x327*x96 + x331*x35;
    double x425 = rho*x279;
    double x426 = x147*x351 + x32*x352 + x350*x91;
    double x427 = x147*x343 + x32*x345 + x344*x91;
    double x428 = x147*x329 + x32*x331 + x327*x91;
    double x429 = rho*x268;
    
    res_0[0] = x202*(x232*x239 + x250*x251 + x252*x322 - x338*x339);
    res_0[1] = x202*(x239*x340 + x251*x342 + x252*x341 - x339*x346);
    res_0[2] = x202*(x239*x349 + x251*x347 + x252*x348 - x339*x353);
    res_0[3] = x202*(x232*x354 + x250*x355 + x322*x356 - x338*x357);
    res_0[4] = x202*(x340*x354 + x341*x356 + x342*x355 - x346*x357);
    res_0[5] = x202*(x347*x355 + x348*x356 + x349*x354 - x353*x357);
    res_0[6] = x202*(x232*x358 + x250*x359 + x322*x360 - x338*x361);
    res_0[7] = x202*(x340*x358 + x341*x360 + x342*x359 - x346*x361);
    res_0[8] = x202*(x347*x359 + x348*x360 + x349*x358 - x353*x361);
    res_0[9] = x202*(x232*x362 + x250*x363 + x322*x364 - x338*x365);
    res_0[10] = x202*(x340*x362 + x341*x364 + x342*x363 - x346*x365);
    res_0[11] = x202*(x347*x363 + x348*x364 + x349*x362 - x353*x365);
    res_0[12] = x202*(x232*x366 + x250*x367 + x322*x368 - x338*x369);
    res_0[13] = x202*(x340*x366 + x341*x368 + x342*x367 - x346*x369);
    res_0[14] = x202*(x347*x367 + x348*x368 + x349*x366 - x353*x369);
    res_0[15] = x202*(x232*x370 + x250*x371 + x322*x372 - x338*x373);
    res_0[16] = x202*(x340*x370 + x341*x372 + x342*x371 - x346*x373);
    res_0[17] = x202*(x347*x371 + x348*x372 + x349*x370 - x353*x373);
    res_0[18] = x202*(x232*x374 + x250*x375 + x322*x376 - x338*x377);
    res_0[19] = x202*(x340*x374 + x341*x376 + x342*x375 - x346*x377);
    res_0[20] = x202*(x347*x375 + x348*x376 + x349*x374 - x353*x377);
    res_0[21] = x202*(x232*x378 + x250*x379 + x322*x380 - x338*x381);
    res_0[22] = x202*(x340*x378 + x341*x380 + x342*x379 - x346*x381);
    res_0[23] = x202*(x347*x379 + x348*x380 + x349*x378 - x353*x381);
    res_0[24] = x202*(x232*x382 + x250*x383 + x322*x384 - x338*x385);
    res_0[25] = x202*(x340*x382 + x341*x384 + x342*x383 - x346*x385);
    res_0[26] = x202*(x347*x383 + x348*x384 + x349*x382 - x353*x385);
    res_0[27] = x202*(x232*x386 + x250*x387 + x322*x388 - x338*x389);
    res_0[28] = x202*(x340*x386 + x341*x388 + x342*x387 - x346*x389);
    res_0[29] = x202*(x347*x387 + x348*x388 + x349*x386 - x353*x389);
    res_0[30] = x202*(x232*x390 + x250*x391 + x322*x392 - x338*x393);
    res_0[31] = x202*(x340*x390 + x341*x392 + x342*x391 - x346*x393);
    res_0[32] = x202*(x347*x391 + x348*x392 + x349*x390 - x353*x393);
    res_0[33] = x202*(x232*x394 + x250*x395 + x322*x396 - x338*x397);
    res_0[34] = x202*(x340*x394 + x341*x396 + x342*x395 - x346*x397);
    res_0[35] = x202*(x347*x395 + x348*x396 + x349*x394 - x353*x397);
    res_0[36] = x202*(x232*x398 + x250*x399 + x322*x400 - x338*x401);
    res_0[37] = x202*(x340*x398 + x341*x400 + x342*x399 - x346*x401);
    res_0[38] = x202*(x347*x399 + x348*x400 + x349*x398 - x353*x401);
    res_0[39] = x202*(x232*x402 + x250*x403 + x322*x404 - x338*x405);
    res_0[40] = x202*(x340*x402 + x341*x404 + x342*x403 - x346*x405);
    res_0[41] = x202*(x347*x403 + x348*x404 + x349*x402 - x353*x405);
    res_0[42] = x202*(x232*x406 + x250*x407 + x322*x408 - x338*x409);
    res_0[43] = x202*(x340*x406 + x341*x408 + x342*x407 - x346*x409);
    res_0[44] = x202*(x347*x407 + x348*x408 + x349*x406 - x353*x409);
    res_0[45] = x202*(x232*x410 + x250*x411 + x322*x412 - x338*x413);
    res_0[46] = x202*(x340*x410 + x341*x412 + x342*x411 - x346*x413);
    res_0[47] = x202*(x347*x411 + x348*x412 + x349*x410 - x353*x413);
    res_0[48] = x202*(x232*x414 + x250*x415 + x322*x416 - x338*x417);
    res_0[49] = x202*(x340*x414 + x341*x416 + x342*x415 - x346*x417);
    res_0[50] = x202*(x347*x415 + x348*x416 + x349*x414 - x353*x417);
    res_0[51] = x202*(x232*x418 + x250*x419 + x322*x420 - x338*x421);
    res_0[52] = x202*(x340*x418 + x341*x420 + x342*x419 - x346*x421);
    res_0[53] = x202*(x347*x419 + x348*x420 + x349*x418 - x353*x421);
    res_0[54] = x202*(x232*x422 + x250*x423 + x322*x424 - x338*x425);
    res_0[55] = x202*(x340*x422 + x341*x424 + x342*x423 - x346*x425);
    res_0[56] = x202*(x347*x423 + x348*x424 + x349*x422 - x353*x425);
    res_0[57] = x202*(x232*x426 + x250*x427 + x322*x428 - x338*x429);
    res_0[58] = x202*(x340*x426 + x341*x428 + x342*x427 - x346*x429);
    res_0[59] = x202*(x347*x427 + x348*x428 + x349*x426 - x353*x429);
}

Conf_Forces_API void Integration_C3D20R_dynamic(size_t num_elem,double Coords[][20][3],double *rho,double Element_U[][20][3],double Element_V[][20][3],double Element_A[][20][3],double S[][8][6],
    double PENER[][8],double SENER[][8],double Conf_Force[][20][3]){
    
    //Inputs:
    //    Coords.........Element_Nodal_Coordinates  [num_elem][numNodes,3]
    //    Element_U......Element_Nodal_Displacements[num_elem][numNodes,3]
    //    Element_V......Element_Nodal_Velocity     [num_elem][numNodes,3]
    //    Element_A......Element_Nodal_Acceleration [num_elem][numNodes,3]
    //    S..............Stress at integration points [num_elem][num_int_points,6]
    //    PENER..........Plastic strain energy [num_elem][num_int_points]
    //    SENER..........Elastic strain energy [num_elem][num_int_points]
    //    
    double weights[8]={1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,};
    double int_points[8][3]={
        {-0.5773502691896258,-0.5773502691896258,-0.5773502691896258,},
        {0.5773502691896258,-0.5773502691896258,-0.5773502691896258,},
        {-0.5773502691896258,0.5773502691896258,-0.5773502691896258,},
        {0.5773502691896258,0.5773502691896258,-0.5773502691896258,},
        {-0.5773502691896258,-0.5773502691896258,0.5773502691896258,},
        {0.5773502691896258,-0.5773502691896258,0.5773502691896258,},
        {-0.5773502691896258,0.5773502691896258,0.5773502691896258,},
        {0.5773502691896258,0.5773502691896258,0.5773502691896258,}};

    #pragma omp parallel for
    for (size_t i=0;i<num_elem;i++){
        double TMP[20][3];
        for (size_t j=0;j<8;j++){
            C3D20R_dynamic((double*) int_points[j],(double*) Coords[i],rho[i],(double*) Element_U[i],(double*) Element_V[i],(double*) Element_A[i],(double*) S[i][j],PENER[i][j],SENER[i][j],(double*) TMP);
            for (size_t num_node=0;num_node<20;num_node++){
                for (size_t ii=0;ii<3;ii++){
                    Conf_Force[i][num_node][ii]+=TMP[num_node][ii]*weights[j];
                }
            }
        }
    }
}

void C3D20R_static_mbf(double *rst,double *coord,double *Element_U,double *S,double PENER,double SENER,double *res_0){
    double r = rst[0];
    double s = rst[1];
    double t = rst[2];
    double *U = Element_U;
    double S11 = S[0];
    double S22 = S[1];
    double S33 = S[2];
    double S12 = S[3];
    double S13 = S[4];
    double S23 = S[5];
    
    double x0 = (1.0/4.0)*r;
    double x1 = -x0;
    double x2 = s*s;
    double x3 = x0*x2 + x1;
    double x4 = (1.0/4.0)*x2;
    double x5 = x4 - 1.0/4.0;
    double x6 = x3 + x5;
    double x7 = (1.0/4.0)*s;
    double x8 = -x7;
    double x9 = r*r;
    double x10 = x7*x9 + x8;
    double x11 = (1.0/4.0)*x9;
    double x12 = x11 - 1.0/4.0;
    double x13 = x10 + x12;
    double x14 = 1.0/4.0 - x4;
    double x15 = x14 + x3;
    double x16 = -x15;
    double x17 = 1.0/4.0 - x11;
    double x18 = x10 + x17;
    double x19 = -x18;
    double x20 = -x13;
    double x21 = -x6;
    double x22 = (1.0/2.0)*t;
    double x23 = r*x22;
    double x24 = -x23;
    double x25 = s*x22;
    double x26 = s*x23;
    double x27 = -x25 + x26;
    double x28 = -x22 - x24 - x27;
    double x29 = -x22;
    double x30 = x23 + x27 + x29;
    double x31 = x25 + x26;
    double x32 = x24 + x29 + x31;
    double x33 = -x22 - x23 - x31;
    double x34 = (1.0/8.0)*r;
    double x35 = s*x34;
    double x36 = -x35;
    double x37 = (1.0/8.0)*x9;
    double x38 = s*x37;
    double x39 = t*x7;
    double x40 = t*x0;
    double x41 = s*x40;
    double x42 = -x41;
    double x43 = x39 + x42;
    double x44 = x36 + x38 + x43;
    double x45 = (1.0/4.0)*t;
    double x46 = -x45;
    double x47 = (1.0/8.0)*x2;
    double x48 = r*x47;
    double x49 = x46 + x48;
    double x50 = 1.0/8.0 - x37;
    double x51 = -x47;
    double x52 = x40 + x51;
    double x53 = x50 + x52;
    double x54 = -x44 - x49 - x53;
    double x55 = -x40;
    double x56 = -x48;
    double x57 = x55 + x56;
    double x58 = x37 - 1.0/8.0;
    double x59 = x45 + x47 + x58;
    double x60 = x44 + x57 + x59;
    double x61 = x35 + x38;
    double x62 = x39 + x41;
    double x63 = x50 + x62;
    double x64 = -x46 - x51 - x57 - x61 - x63;
    double x65 = x40 + x48 + x59 + x61 + x62;
    double x66 = -x39;
    double x67 = x45 + x66;
    double x68 = x41 + x55;
    double x69 = x51 + x68;
    double x70 = x36 + x38;
    double x71 = x48 + x50 + x67 + x69 + x70;
    double x72 = x41 + x66;
    double x73 = x47 + x58;
    double x74 = -x40 - x46 - x56 - x70 - x72 - x73;
    double x75 = x42 + x61;
    double x76 = x53 + x56 + x67 + x75;
    double x77 = -x49 - x55 - x66 - x73 - x75;
    double x78 = x13*coord[30] + x15*coord[45] + x16*coord[33] + x18*coord[36] + x19*coord[24] + x20*coord[42] + x21*coord[39] + x28*coord[48] + x30*coord[57] + x32*coord[51] + x33*coord[54] + x54*coord[12] + x6*coord[27] + x60*coord[21] + x64*coord[15] + x65*coord[18] + x71*coord[0] + x74*coord[9] + x76*coord[3] + x77*coord[6];
    double x79 = t*t;
    double x80 = x0*x79 + x1;
    double x81 = (1.0/4.0)*x79;
    double x82 = x81 - 1.0/4.0;
    double x83 = x80 + x82;
    double x84 = t*x11 + x46;
    double x85 = x12 + x84;
    double x86 = 1.0/4.0 - x81;
    double x87 = x80 + x86;
    double x88 = -x87;
    double x89 = x17 + x84;
    double x90 = -x89;
    double x91 = -x85;
    double x92 = -x83;
    double x93 = (1.0/2.0)*s;
    double x94 = r*x93;
    double x95 = -x94;
    double x96 = -x27 - x93 - x95;
    double x97 = -x93;
    double x98 = x27 + x94 + x97;
    double x99 = x31 + x95 + x97;
    double x100 = -x31 - x93 - x94;
    double x101 = t*x34;
    double x102 = -x101;
    double x103 = t*x37;
    double x104 = x102 + x103 + x43;
    double x105 = (1.0/8.0)*x79;
    double x106 = r*x105;
    double x107 = x106 + x8;
    double x108 = -x105;
    double x109 = s*x0;
    double x110 = x108 + x109;
    double x111 = x110 + x50;
    double x112 = -x104 - x107 - x111;
    double x113 = -x109;
    double x114 = -x106;
    double x115 = x113 + x114;
    double x116 = x105 + x58 + x7;
    double x117 = x104 + x115 + x116;
    double x118 = x101 + x103;
    double x119 = -x108 - x115 - x118 - x63 - x8;
    double x120 = x106 + x109 + x116 + x118 + x62;
    double x121 = x108 + x113;
    double x122 = x102 + x103 + x72;
    double x123 = x106 + x121 + x122 + x50 + x7;
    double x124 = x105 + x58;
    double x125 = -x109 - x114 - x122 - x124 - x8;
    double x126 = x118 + x42 + x66;
    double x127 = x111 + x114 + x126 + x7;
    double x128 = -x107 - x113 - x124 - x126;
    double x129 = x100*coord[41] + x112*coord[11] + x117*coord[23] + x119*coord[8] + x120*coord[20] + x123*coord[2] + x125*coord[14] + x127*coord[5] + x128*coord[17] + x83*coord[53] + x85*coord[38] + x87*coord[59] + x88*coord[50] + x89*coord[32] + x90*coord[26] + x91*coord[44] + x92*coord[56] + x96*coord[35] + x98*coord[47] + x99*coord[29];
    double x130 = x13*coord[32] + x15*coord[47] + x16*coord[35] + x18*coord[38] + x19*coord[26] + x20*coord[44] + x21*coord[41] + x28*coord[50] + x30*coord[59] + x32*coord[53] + x33*coord[56] + x54*coord[14] + x6*coord[29] + x60*coord[23] + x64*coord[17] + x65*coord[20] + x71*coord[2] + x74*coord[11] + x76*coord[5] + x77*coord[8];
    double x131 = x100*coord[39] + x112*coord[9] + x117*coord[21] + x119*coord[6] + x120*coord[18] + x123*coord[0] + x125*coord[12] + x127*coord[3] + x128*coord[15] + x83*coord[51] + x85*coord[36] + x87*coord[57] + x88*coord[48] + x89*coord[30] + x90*coord[24] + x91*coord[42] + x92*coord[54] + x96*coord[33] + x98*coord[45] + x99*coord[27];
    double x132 = x130*x131;
    double x133 = x129*x78 - x132;
    double x134 = x100*coord[40] + x112*coord[10] + x117*coord[22] + x119*coord[7] + x120*coord[19] + x123*coord[1] + x125*coord[13] + x127*coord[4] + x128*coord[16] + x83*coord[52] + x85*coord[37] + x87*coord[58] + x88*coord[49] + x89*coord[31] + x90*coord[25] + x91*coord[43] + x92*coord[55] + x96*coord[34] + x98*coord[46] + x99*coord[28];
    double x135 = x7*x79 + x8;
    double x136 = x135 + x82;
    double x137 = t*x4 + x46;
    double x138 = x137 + x5;
    double x139 = x135 + x86;
    double x140 = -x139;
    double x141 = x137 + x14;
    double x142 = -x141;
    double x143 = -x138;
    double x144 = -x136;
    double x145 = (1.0/2.0)*r;
    double x146 = x24 + x26;
    double x147 = -x145 - x146 - x95;
    double x148 = -x145;
    double x149 = x146 + x148 + x94;
    double x150 = x23 + x26;
    double x151 = x148 + x150 + x95;
    double x152 = -x145 - x150 - x94;
    double x153 = t*x47;
    double x154 = -x153;
    double x155 = (1.0/8.0)*s*t;
    double x156 = x154 + x155;
    double x157 = s*x105;
    double x158 = -x157;
    double x159 = x0 + x158;
    double x160 = x105 + x47 - 1.0/8.0;
    double x161 = x160 + x68;
    double x162 = x113 + x156 + x159 + x161;
    double x163 = x1 + 1.0/8.0;
    double x164 = x157 + x69;
    double x165 = -x110 - x156 - x163 - x164;
    double x166 = x121 + x153;
    double x167 = x41 + x52;
    double x168 = -x155 - x158 - x163 - x166 - x167;
    double x169 = x109 + x153;
    double x170 = x157 + x160 + x40 + x41;
    double x171 = x0 + x155 + x169 + x170;
    double x172 = -x155;
    double x173 = x172 + 1.0/8.0;
    double x174 = x0 + x164 + x166 + x173;
    double x175 = x1 + x172;
    double x176 = -x158 - x161 - x169 - x175;
    double x177 = -x113 - x154 - x170 - x175;
    double x178 = x110 + x154 + x159 + x167 + x173;
    double x179 = x136*coord[59] + x138*coord[47] + x139*coord[53] + x140*coord[50] + x141*coord[29] + x142*coord[35] + x143*coord[41] + x144*coord[56] + x147*coord[26] + x149*coord[38] + x151*coord[32] + x152*coord[44] + x162*coord[5] + x165*coord[17] + x168*coord[8] + x171*coord[20] + x174*coord[2] + x176*coord[14] + x177*coord[11] + x178*coord[23];
    double x180 = x179*x78;
    double x181 = x13*coord[31] + x15*coord[46] + x16*coord[34] + x18*coord[37] + x19*coord[25] + x20*coord[43] + x21*coord[40] + x28*coord[49] + x30*coord[58] + x32*coord[52] + x33*coord[55] + x54*coord[13] + x6*coord[28] + x60*coord[22] + x64*coord[16] + x65*coord[19] + x71*coord[1] + x74*coord[10] + x76*coord[4] + x77*coord[7];
    double x182 = x136*coord[57] + x138*coord[45] + x139*coord[51] + x140*coord[48] + x141*coord[27] + x142*coord[33] + x143*coord[39] + x144*coord[54] + x147*coord[24] + x149*coord[36] + x151*coord[30] + x152*coord[42] + x162*coord[3] + x165*coord[15] + x168*coord[6] + x171*coord[18] + x174*coord[0] + x176*coord[12] + x177*coord[9] + x178*coord[21];
    double x183 = x129*x182;
    double x184 = x136*coord[58] + x138*coord[46] + x139*coord[52] + x140*coord[49] + x141*coord[28] + x142*coord[34] + x143*coord[40] + x144*coord[55] + x147*coord[25] + x149*coord[37] + x151*coord[31] + x152*coord[43] + x162*coord[4] + x165*coord[16] + x168*coord[7] + x171*coord[19] + x174*coord[1] + x176*coord[13] + x177*coord[10] + x178*coord[22];
    double x185 = x129*x184*x78 + x130*x134*x182 + x131*x179*x181 - x132*x184 - x134*x180 - x181*x183;
    double x186 = 1.0/x185;
    double x187 = x186*(x136*U[59] + x138*U[47] + x139*U[53] + x140*U[50] + x141*U[29] + x142*U[35] + x143*U[41] + x144*U[56] + x147*U[26] + x149*U[38] + x151*U[32] + x152*U[44] + x162*U[5] + x165*U[17] + x168*U[8] + x171*U[20] + x174*U[2] + x176*U[14] + x177*U[11] + x178*U[23]);
    double x188 = x130*x182 - x180;
    double x189 = x186*(x100*U[41] + x112*U[11] + x117*U[23] + x119*U[8] + x120*U[20] + x123*U[2] + x125*U[14] + x127*U[5] + x128*U[17] + x83*U[53] + x85*U[38] + x87*U[59] + x88*U[50] + x89*U[32] + x90*U[26] + x91*U[44] + x92*U[56] + x96*U[35] + x98*U[47] + x99*U[29]);
    double x190 = x131*x179 - x183;
    double x191 = x186*(x13*U[32] + x15*U[47] + x16*U[35] + x18*U[38] + x19*U[26] + x20*U[44] + x21*U[41] + x28*U[50] + x30*U[59] + x32*U[53] + x33*U[56] + x54*U[14] + x6*U[29] + x60*U[23] + x64*U[17] + x65*U[20] + x71*U[2] + x74*U[11] + x76*U[5] + x77*U[8]);
    double x192 = x133*x187 + x188*x189 + x190*x191;
    double x193 = -x130*x184 + x179*x181;
    double x194 = x186*(x100*U[40] + x112*U[10] + x117*U[22] + x119*U[7] + x120*U[19] + x123*U[1] + x125*U[13] + x127*U[4] + x128*U[16] + x83*U[52] + x85*U[37] + x87*U[58] + x88*U[49] + x89*U[31] + x90*U[25] + x91*U[43] + x92*U[55] + x96*U[34] + x98*U[46] + x99*U[28]);
    double x195 = -x129*x181 + x130*x134;
    double x196 = x186*(x136*U[58] + x138*U[46] + x139*U[52] + x140*U[49] + x141*U[28] + x142*U[34] + x143*U[40] + x144*U[55] + x147*U[25] + x149*U[37] + x151*U[31] + x152*U[43] + x162*U[4] + x165*U[16] + x168*U[7] + x171*U[19] + x174*U[1] + x176*U[13] + x177*U[10] + x178*U[22]);
    double x197 = x129*x184 - x134*x179;
    double x198 = x186*(x13*U[31] + x15*U[46] + x16*U[34] + x18*U[37] + x19*U[25] + x20*U[43] + x21*U[40] + x28*U[49] + x30*U[58] + x32*U[52] + x33*U[55] + x54*U[13] + x6*U[28] + x60*U[22] + x64*U[16] + x65*U[19] + x71*U[1] + x74*U[10] + x76*U[4] + x77*U[7]);
    double x199 = x193*x194 + x195*x196 + x197*x198;
    double x200 = x187*x195 + x189*x193 + x191*x197;
    double x201 = x133*x196 + x188*x194 + x190*x198 + 1.0;
    double x202 = x192*x199 - x200*x201;
    double x203 = x186*(x136*U[57] + x138*U[45] + x139*U[51] + x140*U[48] + x141*U[27] + x142*U[33] + x143*U[39] + x144*U[54] + x147*U[24] + x149*U[36] + x151*U[30] + x152*U[42] + x162*U[3] + x165*U[15] + x168*U[6] + x171*U[18] + x174*U[0] + x176*U[12] + x177*U[9] + x178*U[21]);
    double x204 = x186*(x100*U[39] + x112*U[9] + x117*U[21] + x119*U[6] + x120*U[18] + x123*U[0] + x125*U[12] + x127*U[3] + x128*U[15] + x83*U[51] + x85*U[36] + x87*U[57] + x88*U[48] + x89*U[30] + x90*U[24] + x91*U[42] + x92*U[54] + x96*U[33] + x98*U[45] + x99*U[27]);
    double x205 = x186*(x13*U[30] + x15*U[45] + x16*U[33] + x18*U[36] + x19*U[24] + x20*U[42] + x21*U[39] + x28*U[48] + x30*U[57] + x32*U[51] + x33*U[54] + x54*U[12] + x6*U[27] + x60*U[21] + x64*U[15] + x65*U[18] + x71*U[0] + x74*U[9] + x76*U[3] + x77*U[6]);
    double x206 = x133*x203 + x188*x204 + x190*x205;
    double x207 = x193*x204 + x195*x203 + x197*x205 + 1.0;
    double x208 = -x192*x207 + x200*x206;
    double x209 = -x199*x206 + x201*x207;
    double x210 = S11*x202 + S12*x208 + S13*x209;
    double x211 = S12*x202 + S22*x208 + S23*x209;
    double x212 = S13*x202 + S23*x208 + S33*x209;
    double x213 = -x199*x211 - x200*x212 - x207*x210;
    double x214 = -x181*x182 + x184*x78;
    double x215 = x123*x186;
    double x216 = x131*x181 - x134*x78;
    double x217 = x174*x186;
    double x218 = -x131*x184 + x134*x182;
    double x219 = x186*x71;
    double x220 = x214*x215 + x216*x217 + x218*x219;
    double x221 = x194*x214 + x196*x216 + x198*x218;
    double x222 = x187*x216 + x189*x214 + x191*x218 + 1.0;
    double x223 = -x199*x222 + x200*x221;
    double x224 = x203*x216 + x204*x214 + x205*x218;
    double x225 = -x200*x224 + x207*x222;
    double x226 = x199*x224 - x207*x221;
    double x227 = S11*x223 + S12*x225 + S13*x226;
    double x228 = S12*x223 + S22*x225 + S23*x226;
    double x229 = S13*x223 + S23*x225 + S33*x226;
    double x230 = -x199*x228 - x200*x229 - x207*x227;
    double x231 = x133*x217 + x188*x215 + x190*x219;
    double x232 = x193*x215 + x195*x217 + x197*x219;
    double x233 = -x192*x221 + x201*x222;
    double x234 = x192*x224 - x206*x222;
    double x235 = -x201*x224 + x206*x221;
    double x236 = S11*x233 + S12*x234 + S13*x235;
    double x237 = S12*x233 + S22*x234 + S23*x235;
    double x238 = S13*x233 + S23*x234 + S33*x235;
    double x239 = -1.0*PENER - 1.0*SENER;
    double x240 = -x199*x237 - x200*x238 - x207*x236 - x239;
    double x241 = -x192*x212 - x201*x211 - x206*x210;
    double x242 = -x192*x238 - x201*x237 - x206*x236;
    double x243 = -x192*x229 - x201*x228 - x206*x227 - x239;
    double x244 = -x221*x228 - x222*x229 - x224*x227;
    double x245 = -x221*x237 - x222*x238 - x224*x236;
    double x246 = -x210*x224 - x211*x221 - x212*x222 - x239;
    double x247 = x127*x186;
    double x248 = x162*x186;
    double x249 = x186*x76;
    double x250 = x214*x247 + x216*x248 + x218*x249;
    double x251 = x133*x248 + x188*x247 + x190*x249;
    double x252 = x193*x247 + x195*x248 + x197*x249;
    double x253 = x119*x186;
    double x254 = x168*x186;
    double x255 = x186*x77;
    double x256 = x214*x253 + x216*x254 + x218*x255;
    double x257 = x133*x254 + x188*x253 + x190*x255;
    double x258 = x193*x253 + x195*x254 + x197*x255;
    double x259 = x112*x186;
    double x260 = x177*x186;
    double x261 = x186*x74;
    double x262 = x214*x259 + x216*x260 + x218*x261;
    double x263 = x133*x260 + x188*x259 + x190*x261;
    double x264 = x193*x259 + x195*x260 + x197*x261;
    double x265 = x125*x186;
    double x266 = x176*x186;
    double x267 = x186*x54;
    double x268 = x214*x265 + x216*x266 + x218*x267;
    double x269 = x133*x266 + x188*x265 + x190*x267;
    double x270 = x193*x265 + x195*x266 + x197*x267;
    double x271 = x128*x186;
    double x272 = x165*x186;
    double x273 = x186*x64;
    double x274 = x214*x271 + x216*x272 + x218*x273;
    double x275 = x133*x272 + x188*x271 + x190*x273;
    double x276 = x193*x271 + x195*x272 + x197*x273;
    double x277 = x120*x186;
    double x278 = x171*x186;
    double x279 = x186*x65;
    double x280 = x214*x277 + x216*x278 + x218*x279;
    double x281 = x133*x278 + x188*x277 + x190*x279;
    double x282 = x193*x277 + x195*x278 + x197*x279;
    double x283 = x117*x186;
    double x284 = x178*x186;
    double x285 = x186*x60;
    double x286 = x214*x283 + x216*x284 + x218*x285;
    double x287 = x133*x284 + x188*x283 + x190*x285;
    double x288 = x193*x283 + x195*x284 + x197*x285;
    double x289 = x186*x90;
    double x290 = x147*x186;
    double x291 = x186*x19;
    double x292 = x214*x289 + x216*x290 + x218*x291;
    double x293 = x133*x290 + x188*x289 + x190*x291;
    double x294 = x193*x289 + x195*x290 + x197*x291;
    double x295 = x186*x99;
    double x296 = x141*x186;
    double x297 = x186*x6;
    double x298 = x214*x295 + x216*x296 + x218*x297;
    double x299 = x133*x296 + x188*x295 + x190*x297;
    double x300 = x193*x295 + x195*x296 + x197*x297;
    double x301 = x186*x89;
    double x302 = x151*x186;
    double x303 = x13*x186;
    double x304 = x214*x301 + x216*x302 + x218*x303;
    double x305 = x133*x302 + x188*x301 + x190*x303;
    double x306 = x193*x301 + x195*x302 + x197*x303;
    double x307 = x186*x96;
    double x308 = x142*x186;
    double x309 = x16*x186;
    double x310 = x214*x307 + x216*x308 + x218*x309;
    double x311 = x133*x308 + x188*x307 + x190*x309;
    double x312 = x193*x307 + x195*x308 + x197*x309;
    double x313 = x186*x85;
    double x314 = x149*x186;
    double x315 = x18*x186;
    double x316 = x214*x313 + x216*x314 + x218*x315;
    double x317 = x133*x314 + x188*x313 + x190*x315;
    double x318 = x193*x313 + x195*x314 + x197*x315;
    double x319 = x100*x186;
    double x320 = x143*x186;
    double x321 = x186*x21;
    double x322 = x214*x319 + x216*x320 + x218*x321;
    double x323 = x133*x320 + x188*x319 + x190*x321;
    double x324 = x193*x319 + x195*x320 + x197*x321;
    double x325 = x186*x91;
    double x326 = x152*x186;
    double x327 = x186*x20;
    double x328 = x214*x325 + x216*x326 + x218*x327;
    double x329 = x133*x326 + x188*x325 + x190*x327;
    double x330 = x193*x325 + x195*x326 + x197*x327;
    double x331 = x186*x98;
    double x332 = x138*x186;
    double x333 = x15*x186;
    double x334 = x214*x331 + x216*x332 + x218*x333;
    double x335 = x133*x332 + x188*x331 + x190*x333;
    double x336 = x193*x331 + x195*x332 + x197*x333;
    double x337 = x186*x88;
    double x338 = x140*x186;
    double x339 = x186*x28;
    double x340 = x214*x337 + x216*x338 + x218*x339;
    double x341 = x133*x338 + x188*x337 + x190*x339;
    double x342 = x193*x337 + x195*x338 + x197*x339;
    double x343 = x186*x83;
    double x344 = x139*x186;
    double x345 = x186*x32;
    double x346 = x214*x343 + x216*x344 + x218*x345;
    double x347 = x133*x344 + x188*x343 + x190*x345;
    double x348 = x193*x343 + x195*x344 + x197*x345;
    double x349 = x186*x92;
    double x350 = x144*x186;
    double x351 = x186*x33;
    double x352 = x214*x349 + x216*x350 + x218*x351;
    double x353 = x133*x350 + x188*x349 + x190*x351;
    double x354 = x193*x349 + x195*x350 + x197*x351;
    double x355 = x186*x87;
    double x356 = x136*x186;
    double x357 = x186*x30;
    double x358 = x214*x355 + x216*x356 + x218*x357;
    double x359 = x133*x356 + x188*x355 + x190*x357;
    double x360 = x193*x355 + x195*x356 + x197*x357;
    
    res_0[0] = x185*(x213*x220 + x230*x231 + x232*x240);
    res_0[1] = x185*(x220*x241 + x231*x243 + x232*x242);
    res_0[2] = x185*(x220*x246 + x231*x244 + x232*x245);
    res_0[3] = x185*(x213*x250 + x230*x251 + x240*x252);
    res_0[4] = x185*(x241*x250 + x242*x252 + x243*x251);
    res_0[5] = x185*(x244*x251 + x245*x252 + x246*x250);
    res_0[6] = x185*(x213*x256 + x230*x257 + x240*x258);
    res_0[7] = x185*(x241*x256 + x242*x258 + x243*x257);
    res_0[8] = x185*(x244*x257 + x245*x258 + x246*x256);
    res_0[9] = x185*(x213*x262 + x230*x263 + x240*x264);
    res_0[10] = x185*(x241*x262 + x242*x264 + x243*x263);
    res_0[11] = x185*(x244*x263 + x245*x264 + x246*x262);
    res_0[12] = x185*(x213*x268 + x230*x269 + x240*x270);
    res_0[13] = x185*(x241*x268 + x242*x270 + x243*x269);
    res_0[14] = x185*(x244*x269 + x245*x270 + x246*x268);
    res_0[15] = x185*(x213*x274 + x230*x275 + x240*x276);
    res_0[16] = x185*(x241*x274 + x242*x276 + x243*x275);
    res_0[17] = x185*(x244*x275 + x245*x276 + x246*x274);
    res_0[18] = x185*(x213*x280 + x230*x281 + x240*x282);
    res_0[19] = x185*(x241*x280 + x242*x282 + x243*x281);
    res_0[20] = x185*(x244*x281 + x245*x282 + x246*x280);
    res_0[21] = x185*(x213*x286 + x230*x287 + x240*x288);
    res_0[22] = x185*(x241*x286 + x242*x288 + x243*x287);
    res_0[23] = x185*(x244*x287 + x245*x288 + x246*x286);
    res_0[24] = x185*(x213*x292 + x230*x293 + x240*x294);
    res_0[25] = x185*(x241*x292 + x242*x294 + x243*x293);
    res_0[26] = x185*(x244*x293 + x245*x294 + x246*x292);
    res_0[27] = x185*(x213*x298 + x230*x299 + x240*x300);
    res_0[28] = x185*(x241*x298 + x242*x300 + x243*x299);
    res_0[29] = x185*(x244*x299 + x245*x300 + x246*x298);
    res_0[30] = x185*(x213*x304 + x230*x305 + x240*x306);
    res_0[31] = x185*(x241*x304 + x242*x306 + x243*x305);
    res_0[32] = x185*(x244*x305 + x245*x306 + x246*x304);
    res_0[33] = x185*(x213*x310 + x230*x311 + x240*x312);
    res_0[34] = x185*(x241*x310 + x242*x312 + x243*x311);
    res_0[35] = x185*(x244*x311 + x245*x312 + x246*x310);
    res_0[36] = x185*(x213*x316 + x230*x317 + x240*x318);
    res_0[37] = x185*(x241*x316 + x242*x318 + x243*x317);
    res_0[38] = x185*(x244*x317 + x245*x318 + x246*x316);
    res_0[39] = x185*(x213*x322 + x230*x323 + x240*x324);
    res_0[40] = x185*(x241*x322 + x242*x324 + x243*x323);
    res_0[41] = x185*(x244*x323 + x245*x324 + x246*x322);
    res_0[42] = x185*(x213*x328 + x230*x329 + x240*x330);
    res_0[43] = x185*(x241*x328 + x242*x330 + x243*x329);
    res_0[44] = x185*(x244*x329 + x245*x330 + x246*x328);
    res_0[45] = x185*(x213*x334 + x230*x335 + x240*x336);
    res_0[46] = x185*(x241*x334 + x242*x336 + x243*x335);
    res_0[47] = x185*(x244*x335 + x245*x336 + x246*x334);
    res_0[48] = x185*(x213*x340 + x230*x341 + x240*x342);
    res_0[49] = x185*(x241*x340 + x242*x342 + x243*x341);
    res_0[50] = x185*(x244*x341 + x245*x342 + x246*x340);
    res_0[51] = x185*(x213*x346 + x230*x347 + x240*x348);
    res_0[52] = x185*(x241*x346 + x242*x348 + x243*x347);
    res_0[53] = x185*(x244*x347 + x245*x348 + x246*x346);
    res_0[54] = x185*(x213*x352 + x230*x353 + x240*x354);
    res_0[55] = x185*(x241*x352 + x242*x354 + x243*x353);
    res_0[56] = x185*(x244*x353 + x245*x354 + x246*x352);
    res_0[57] = x185*(x213*x358 + x230*x359 + x240*x360);
    res_0[58] = x185*(x241*x358 + x242*x360 + x243*x359);
    res_0[59] = x185*(x244*x359 + x245*x360 + x246*x358);
}

Conf_Forces_API void Integration_C3D20R_static_mbf(size_t num_elem,double Coords[][20][3],double Element_U[][20][3],double S[][8][6],
    double PENER[][8],double SENER[][8],double Conf_Force[][20][3]){
    
    //Inputs:
    //    Coords.........Element_Nodal_Coordinates [num_elem][numNodes,3]
    //    Element_U......Element_Nodal_Displacements [num_elem][numNodes,3]
    //    S..............Stress at integration points [num_elem][num_int_points,6]
    //    PENER..........Plastic strain energy [num_elem][num_int_points]
    //    SENER..........Elastic strain energy [num_elem][num_int_points]
    //    
    double weights[8]={1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,};
    double int_points[8][3]={
        {-0.5773502691896258,-0.5773502691896258,-0.5773502691896258,},
        {0.5773502691896258,-0.5773502691896258,-0.5773502691896258,},
        {-0.5773502691896258,0.5773502691896258,-0.5773502691896258,},
        {0.5773502691896258,0.5773502691896258,-0.5773502691896258,},
        {-0.5773502691896258,-0.5773502691896258,0.5773502691896258,},
        {0.5773502691896258,-0.5773502691896258,0.5773502691896258,},
        {-0.5773502691896258,0.5773502691896258,0.5773502691896258,},
        {0.5773502691896258,0.5773502691896258,0.5773502691896258,}};

    #pragma omp parallel for
    for (size_t i=0;i<num_elem;i++){
        double TMP[20][3];
        for (size_t j=0;j<8;j++){
            C3D20R_static_mbf((double*) int_points[j],(double*) Coords[i],(double*) Element_U[i],(double*) S[i][j],PENER[i][j],SENER[i][j],(double*) TMP);
            for (size_t num_node=0;num_node<20;num_node++){
                for (size_t ii=0;ii<3;ii++){
                    Conf_Force[i][num_node][ii]+=TMP[num_node][ii]*weights[j];
                }
            }
        }
    }
}

void C3D20R_static_dbf(double *rst,double *coord,double *Element_U,double *S,double PENER,double SENER,double *res_0){
    double r = rst[0];
    double s = rst[1];
    double t = rst[2];
    double *U = Element_U;
    double S11 = S[0];
    double S22 = S[1];
    double S33 = S[2];
    double S12 = S[3];
    double S13 = S[4];
    double S23 = S[5];
    
    double x0 = (1.0/4.0)*r;
    double x1 = -x0;
    double x2 = s*s;
    double x3 = x0*x2 + x1;
    double x4 = (1.0/4.0)*x2;
    double x5 = x4 - 1.0/4.0;
    double x6 = x3 + x5;
    double x7 = (1.0/4.0)*s;
    double x8 = -x7;
    double x9 = r*r;
    double x10 = x7*x9 + x8;
    double x11 = (1.0/4.0)*x9;
    double x12 = x11 - 1.0/4.0;
    double x13 = x10 + x12;
    double x14 = 1.0/4.0 - x4;
    double x15 = x14 + x3;
    double x16 = -x15;
    double x17 = 1.0/4.0 - x11;
    double x18 = x10 + x17;
    double x19 = -x18;
    double x20 = -x13;
    double x21 = -x6;
    double x22 = (1.0/2.0)*t;
    double x23 = r*x22;
    double x24 = -x23;
    double x25 = s*x22;
    double x26 = s*x23;
    double x27 = -x25 + x26;
    double x28 = -x22 - x24 - x27;
    double x29 = -x22;
    double x30 = x23 + x27 + x29;
    double x31 = x25 + x26;
    double x32 = x24 + x29 + x31;
    double x33 = -x22 - x23 - x31;
    double x34 = (1.0/8.0)*r;
    double x35 = s*x34;
    double x36 = -x35;
    double x37 = (1.0/8.0)*x9;
    double x38 = s*x37;
    double x39 = t*x7;
    double x40 = t*x0;
    double x41 = s*x40;
    double x42 = -x41;
    double x43 = x39 + x42;
    double x44 = x36 + x38 + x43;
    double x45 = (1.0/4.0)*t;
    double x46 = -x45;
    double x47 = (1.0/8.0)*x2;
    double x48 = r*x47;
    double x49 = x46 + x48;
    double x50 = 1.0/8.0 - x37;
    double x51 = -x47;
    double x52 = x40 + x51;
    double x53 = x50 + x52;
    double x54 = -x44 - x49 - x53;
    double x55 = -x40;
    double x56 = -x48;
    double x57 = x55 + x56;
    double x58 = x37 - 1.0/8.0;
    double x59 = x45 + x47 + x58;
    double x60 = x44 + x57 + x59;
    double x61 = x35 + x38;
    double x62 = x39 + x41;
    double x63 = x50 + x62;
    double x64 = -x46 - x51 - x57 - x61 - x63;
    double x65 = x40 + x48 + x59 + x61 + x62;
    double x66 = -x39;
    double x67 = x45 + x66;
    double x68 = x41 + x55;
    double x69 = x51 + x68;
    double x70 = x36 + x38;
    double x71 = x48 + x50 + x67 + x69 + x70;
    double x72 = x41 + x66;
    double x73 = x47 + x58;
    double x74 = -x40 - x46 - x56 - x70 - x72 - x73;
    double x75 = x42 + x61;
    double x76 = x53 + x56 + x67 + x75;
    double x77 = -x49 - x55 - x66 - x73 - x75;
    double x78 = x13*coord[30] + x15*coord[45] + x16*coord[33] + x18*coord[36] + x19*coord[24] + x20*coord[42] + x21*coord[39] + x28*coord[48] + x30*coord[57] + x32*coord[51] + x33*coord[54] + x54*coord[12] + x6*coord[27] + x60*coord[21] + x64*coord[15] + x65*coord[18] + x71*coord[0] + x74*coord[9] + x76*coord[3] + x77*coord[6];
    double x79 = t*t;
    double x80 = x0*x79 + x1;
    double x81 = (1.0/4.0)*x79;
    double x82 = x81 - 1.0/4.0;
    double x83 = x80 + x82;
    double x84 = t*x11 + x46;
    double x85 = x12 + x84;
    double x86 = 1.0/4.0 - x81;
    double x87 = x80 + x86;
    double x88 = -x87;
    double x89 = x17 + x84;
    double x90 = -x89;
    double x91 = -x85;
    double x92 = -x83;
    double x93 = (1.0/2.0)*s;
    double x94 = r*x93;
    double x95 = -x94;
    double x96 = -x27 - x93 - x95;
    double x97 = -x93;
    double x98 = x27 + x94 + x97;
    double x99 = x31 + x95 + x97;
    double x100 = -x31 - x93 - x94;
    double x101 = t*x34;
    double x102 = -x101;
    double x103 = t*x37;
    double x104 = x102 + x103 + x43;
    double x105 = (1.0/8.0)*x79;
    double x106 = r*x105;
    double x107 = x106 + x8;
    double x108 = -x105;
    double x109 = s*x0;
    double x110 = x108 + x109;
    double x111 = x110 + x50;
    double x112 = -x104 - x107 - x111;
    double x113 = -x109;
    double x114 = -x106;
    double x115 = x113 + x114;
    double x116 = x105 + x58 + x7;
    double x117 = x104 + x115 + x116;
    double x118 = x101 + x103;
    double x119 = -x108 - x115 - x118 - x63 - x8;
    double x120 = x106 + x109 + x116 + x118 + x62;
    double x121 = x108 + x113;
    double x122 = x102 + x103 + x72;
    double x123 = x106 + x121 + x122 + x50 + x7;
    double x124 = x105 + x58;
    double x125 = -x109 - x114 - x122 - x124 - x8;
    double x126 = x118 + x42 + x66;
    double x127 = x111 + x114 + x126 + x7;
    double x128 = -x107 - x113 - x124 - x126;
    double x129 = x100*coord[41] + x112*coord[11] + x117*coord[23] + x119*coord[8] + x120*coord[20] + x123*coord[2] + x125*coord[14] + x127*coord[5] + x128*coord[17] + x83*coord[53] + x85*coord[38] + x87*coord[59] + x88*coord[50] + x89*coord[32] + x90*coord[26] + x91*coord[44] + x92*coord[56] + x96*coord[35] + x98*coord[47] + x99*coord[29];
    double x130 = x13*coord[32] + x15*coord[47] + x16*coord[35] + x18*coord[38] + x19*coord[26] + x20*coord[44] + x21*coord[41] + x28*coord[50] + x30*coord[59] + x32*coord[53] + x33*coord[56] + x54*coord[14] + x6*coord[29] + x60*coord[23] + x64*coord[17] + x65*coord[20] + x71*coord[2] + x74*coord[11] + x76*coord[5] + x77*coord[8];
    double x131 = x100*coord[39] + x112*coord[9] + x117*coord[21] + x119*coord[6] + x120*coord[18] + x123*coord[0] + x125*coord[12] + x127*coord[3] + x128*coord[15] + x83*coord[51] + x85*coord[36] + x87*coord[57] + x88*coord[48] + x89*coord[30] + x90*coord[24] + x91*coord[42] + x92*coord[54] + x96*coord[33] + x98*coord[45] + x99*coord[27];
    double x132 = x130*x131;
    double x133 = x129*x78 - x132;
    double x134 = x100*coord[40] + x112*coord[10] + x117*coord[22] + x119*coord[7] + x120*coord[19] + x123*coord[1] + x125*coord[13] + x127*coord[4] + x128*coord[16] + x83*coord[52] + x85*coord[37] + x87*coord[58] + x88*coord[49] + x89*coord[31] + x90*coord[25] + x91*coord[43] + x92*coord[55] + x96*coord[34] + x98*coord[46] + x99*coord[28];
    double x135 = x7*x79 + x8;
    double x136 = x135 + x82;
    double x137 = t*x4 + x46;
    double x138 = x137 + x5;
    double x139 = x135 + x86;
    double x140 = -x139;
    double x141 = x137 + x14;
    double x142 = -x141;
    double x143 = -x138;
    double x144 = -x136;
    double x145 = (1.0/2.0)*r;
    double x146 = x24 + x26;
    double x147 = -x145 - x146 - x95;
    double x148 = -x145;
    double x149 = x146 + x148 + x94;
    double x150 = x23 + x26;
    double x151 = x148 + x150 + x95;
    double x152 = -x145 - x150 - x94;
    double x153 = t*x47;
    double x154 = -x153;
    double x155 = (1.0/8.0)*s*t;
    double x156 = x154 + x155;
    double x157 = s*x105;
    double x158 = -x157;
    double x159 = x0 + x158;
    double x160 = x105 + x47 - 1.0/8.0;
    double x161 = x160 + x68;
    double x162 = x113 + x156 + x159 + x161;
    double x163 = x1 + 1.0/8.0;
    double x164 = x157 + x69;
    double x165 = -x110 - x156 - x163 - x164;
    double x166 = x121 + x153;
    double x167 = x41 + x52;
    double x168 = -x155 - x158 - x163 - x166 - x167;
    double x169 = x109 + x153;
    double x170 = x157 + x160 + x40 + x41;
    double x171 = x0 + x155 + x169 + x170;
    double x172 = -x155;
    double x173 = x172 + 1.0/8.0;
    double x174 = x0 + x164 + x166 + x173;
    double x175 = x1 + x172;
    double x176 = -x158 - x161 - x169 - x175;
    double x177 = -x113 - x154 - x170 - x175;
    double x178 = x110 + x154 + x159 + x167 + x173;
    double x179 = x136*coord[59] + x138*coord[47] + x139*coord[53] + x140*coord[50] + x141*coord[29] + x142*coord[35] + x143*coord[41] + x144*coord[56] + x147*coord[26] + x149*coord[38] + x151*coord[32] + x152*coord[44] + x162*coord[5] + x165*coord[17] + x168*coord[8] + x171*coord[20] + x174*coord[2] + x176*coord[14] + x177*coord[11] + x178*coord[23];
    double x180 = x179*x78;
    double x181 = x13*coord[31] + x15*coord[46] + x16*coord[34] + x18*coord[37] + x19*coord[25] + x20*coord[43] + x21*coord[40] + x28*coord[49] + x30*coord[58] + x32*coord[52] + x33*coord[55] + x54*coord[13] + x6*coord[28] + x60*coord[22] + x64*coord[16] + x65*coord[19] + x71*coord[1] + x74*coord[10] + x76*coord[4] + x77*coord[7];
    double x182 = x136*coord[57] + x138*coord[45] + x139*coord[51] + x140*coord[48] + x141*coord[27] + x142*coord[33] + x143*coord[39] + x144*coord[54] + x147*coord[24] + x149*coord[36] + x151*coord[30] + x152*coord[42] + x162*coord[3] + x165*coord[15] + x168*coord[6] + x171*coord[18] + x174*coord[0] + x176*coord[12] + x177*coord[9] + x178*coord[21];
    double x183 = x129*x182;
    double x184 = x136*coord[58] + x138*coord[46] + x139*coord[52] + x140*coord[49] + x141*coord[28] + x142*coord[34] + x143*coord[40] + x144*coord[55] + x147*coord[25] + x149*coord[37] + x151*coord[31] + x152*coord[43] + x162*coord[4] + x165*coord[16] + x168*coord[7] + x171*coord[19] + x174*coord[1] + x176*coord[13] + x177*coord[10] + x178*coord[22];
    double x185 = x129*x184*x78 + x130*x134*x182 + x131*x179*x181 - x132*x184 - x134*x180 - x181*x183;
    double x186 = 1.0/x185;
    double x187 = x186*(x136*U[59] + x138*U[47] + x139*U[53] + x140*U[50] + x141*U[29] + x142*U[35] + x143*U[41] + x144*U[56] + x147*U[26] + x149*U[38] + x151*U[32] + x152*U[44] + x162*U[5] + x165*U[17] + x168*U[8] + x171*U[20] + x174*U[2] + x176*U[14] + x177*U[11] + x178*U[23]);
    double x188 = x130*x182 - x180;
    double x189 = x186*(x100*U[41] + x112*U[11] + x117*U[23] + x119*U[8] + x120*U[20] + x123*U[2] + x125*U[14] + x127*U[5] + x128*U[17] + x83*U[53] + x85*U[38] + x87*U[59] + x88*U[50] + x89*U[32] + x90*U[26] + x91*U[44] + x92*U[56] + x96*U[35] + x98*U[47] + x99*U[29]);
    double x190 = x131*x179 - x183;
    double x191 = x186*(x13*U[32] + x15*U[47] + x16*U[35] + x18*U[38] + x19*U[26] + x20*U[44] + x21*U[41] + x28*U[50] + x30*U[59] + x32*U[53] + x33*U[56] + x54*U[14] + x6*U[29] + x60*U[23] + x64*U[17] + x65*U[20] + x71*U[2] + x74*U[11] + x76*U[5] + x77*U[8]);
    double x192 = x133*x187 + x188*x189 + x190*x191;
    double x193 = -x130*x184 + x179*x181;
    double x194 = x186*(x100*U[40] + x112*U[10] + x117*U[22] + x119*U[7] + x120*U[19] + x123*U[1] + x125*U[13] + x127*U[4] + x128*U[16] + x83*U[52] + x85*U[37] + x87*U[58] + x88*U[49] + x89*U[31] + x90*U[25] + x91*U[43] + x92*U[55] + x96*U[34] + x98*U[46] + x99*U[28]);
    double x195 = -x129*x181 + x130*x134;
    double x196 = x186*(x136*U[58] + x138*U[46] + x139*U[52] + x140*U[49] + x141*U[28] + x142*U[34] + x143*U[40] + x144*U[55] + x147*U[25] + x149*U[37] + x151*U[31] + x152*U[43] + x162*U[4] + x165*U[16] + x168*U[7] + x171*U[19] + x174*U[1] + x176*U[13] + x177*U[10] + x178*U[22]);
    double x197 = x129*x184 - x134*x179;
    double x198 = x186*(x13*U[31] + x15*U[46] + x16*U[34] + x18*U[37] + x19*U[25] + x20*U[43] + x21*U[40] + x28*U[49] + x30*U[58] + x32*U[52] + x33*U[55] + x54*U[13] + x6*U[28] + x60*U[22] + x64*U[16] + x65*U[19] + x71*U[1] + x74*U[10] + x76*U[4] + x77*U[7]);
    double x199 = x193*x194 + x195*x196 + x197*x198;
    double x200 = x187*x195 + x189*x193 + x191*x197;
    double x201 = x133*x196 + x188*x194 + x190*x198;
    double x202 = x201 + 1.0;
    double x203 = x192*x199 - x200*x202;
    double x204 = x186*(x136*U[57] + x138*U[45] + x139*U[51] + x140*U[48] + x141*U[27] + x142*U[33] + x143*U[39] + x144*U[54] + x147*U[24] + x149*U[36] + x151*U[30] + x152*U[42] + x162*U[3] + x165*U[15] + x168*U[6] + x171*U[18] + x174*U[0] + x176*U[12] + x177*U[9] + x178*U[21]);
    double x205 = x186*(x100*U[39] + x112*U[9] + x117*U[21] + x119*U[6] + x120*U[18] + x123*U[0] + x125*U[12] + x127*U[3] + x128*U[15] + x83*U[51] + x85*U[36] + x87*U[57] + x88*U[48] + x89*U[30] + x90*U[24] + x91*U[42] + x92*U[54] + x96*U[33] + x98*U[45] + x99*U[27]);
    double x206 = x186*(x13*U[30] + x15*U[45] + x16*U[33] + x18*U[36] + x19*U[24] + x20*U[42] + x21*U[39] + x28*U[48] + x30*U[57] + x32*U[51] + x33*U[54] + x54*U[12] + x6*U[27] + x60*U[21] + x64*U[15] + x65*U[18] + x71*U[0] + x74*U[9] + x76*U[3] + x77*U[6]);
    double x207 = x133*x204 + x188*x205 + x190*x206;
    double x208 = x193*x205 + x195*x204 + x197*x206;
    double x209 = x208 + 1.0;
    double x210 = -x192*x209 + x200*x207;
    double x211 = -x199*x207 + x202*x209;
    double x212 = S11*x203 + S12*x210 + S13*x211;
    double x213 = S12*x203 + S22*x210 + S23*x211;
    double x214 = S13*x203 + S23*x210 + S33*x211;
    double x215 = -x199*x213 - x200*x214 - x208*x212;
    double x216 = -x181*x182 + x184*x78;
    double x217 = x123*x186;
    double x218 = x131*x181 - x134*x78;
    double x219 = x174*x186;
    double x220 = -x131*x184 + x134*x182;
    double x221 = x186*x71;
    double x222 = x216*x217 + x218*x219 + x220*x221;
    double x223 = x194*x216 + x196*x218 + x198*x220;
    double x224 = x187*x218 + x189*x216 + x191*x220;
    double x225 = x224 + 1.0;
    double x226 = -x199*x225 + x200*x223;
    double x227 = x204*x218 + x205*x216 + x206*x220;
    double x228 = -x200*x227 + x209*x225;
    double x229 = x199*x227 - x209*x223;
    double x230 = S11*x226 + S12*x228 + S13*x229;
    double x231 = S12*x226 + S22*x228 + S23*x229;
    double x232 = S13*x226 + S23*x228 + S33*x229;
    double x233 = -x199*x231 - x200*x232 - x208*x230;
    double x234 = x133*x219 + x188*x217 + x190*x221;
    double x235 = x193*x217 + x195*x219 + x197*x221;
    double x236 = -x192*x223 + x202*x225;
    double x237 = x192*x227 - x207*x225;
    double x238 = -x202*x227 + x207*x223;
    double x239 = S11*x236 + S12*x237 + S13*x238;
    double x240 = S12*x236 + S22*x237 + S23*x238;
    double x241 = S13*x236 + S23*x237 + S33*x238;
    double x242 = -1.0*PENER - 1.0*SENER;
    double x243 = -x199*x240 - x200*x241 - x208*x239 - x242;
    double x244 = -x192*x214 - x201*x213 - x207*x212;
    double x245 = -x192*x241 - x201*x240 - x207*x239;
    double x246 = -x192*x232 - x201*x231 - x207*x230 - x242;
    double x247 = -x223*x231 - x224*x232 - x227*x230;
    double x248 = -x223*x240 - x224*x241 - x227*x239;
    double x249 = -x212*x227 - x213*x223 - x214*x224 - x242;
    double x250 = x127*x186;
    double x251 = x162*x186;
    double x252 = x186*x76;
    double x253 = x216*x250 + x218*x251 + x220*x252;
    double x254 = x133*x251 + x188*x250 + x190*x252;
    double x255 = x193*x250 + x195*x251 + x197*x252;
    double x256 = x119*x186;
    double x257 = x168*x186;
    double x258 = x186*x77;
    double x259 = x216*x256 + x218*x257 + x220*x258;
    double x260 = x133*x257 + x188*x256 + x190*x258;
    double x261 = x193*x256 + x195*x257 + x197*x258;
    double x262 = x112*x186;
    double x263 = x177*x186;
    double x264 = x186*x74;
    double x265 = x216*x262 + x218*x263 + x220*x264;
    double x266 = x133*x263 + x188*x262 + x190*x264;
    double x267 = x193*x262 + x195*x263 + x197*x264;
    double x268 = x125*x186;
    double x269 = x176*x186;
    double x270 = x186*x54;
    double x271 = x216*x268 + x218*x269 + x220*x270;
    double x272 = x133*x269 + x188*x268 + x190*x270;
    double x273 = x193*x268 + x195*x269 + x197*x270;
    double x274 = x128*x186;
    double x275 = x165*x186;
    double x276 = x186*x64;
    double x277 = x216*x274 + x218*x275 + x220*x276;
    double x278 = x133*x275 + x188*x274 + x190*x276;
    double x279 = x193*x274 + x195*x275 + x197*x276;
    double x280 = x120*x186;
    double x281 = x171*x186;
    double x282 = x186*x65;
    double x283 = x216*x280 + x218*x281 + x220*x282;
    double x284 = x133*x281 + x188*x280 + x190*x282;
    double x285 = x193*x280 + x195*x281 + x197*x282;
    double x286 = x117*x186;
    double x287 = x178*x186;
    double x288 = x186*x60;
    double x289 = x216*x286 + x218*x287 + x220*x288;
    double x290 = x133*x287 + x188*x286 + x190*x288;
    double x291 = x193*x286 + x195*x287 + x197*x288;
    double x292 = x186*x90;
    double x293 = x147*x186;
    double x294 = x186*x19;
    double x295 = x216*x292 + x218*x293 + x220*x294;
    double x296 = x133*x293 + x188*x292 + x190*x294;
    double x297 = x193*x292 + x195*x293 + x197*x294;
    double x298 = x186*x99;
    double x299 = x141*x186;
    double x300 = x186*x6;
    double x301 = x216*x298 + x218*x299 + x220*x300;
    double x302 = x133*x299 + x188*x298 + x190*x300;
    double x303 = x193*x298 + x195*x299 + x197*x300;
    double x304 = x186*x89;
    double x305 = x151*x186;
    double x306 = x13*x186;
    double x307 = x216*x304 + x218*x305 + x220*x306;
    double x308 = x133*x305 + x188*x304 + x190*x306;
    double x309 = x193*x304 + x195*x305 + x197*x306;
    double x310 = x186*x96;
    double x311 = x142*x186;
    double x312 = x16*x186;
    double x313 = x216*x310 + x218*x311 + x220*x312;
    double x314 = x133*x311 + x188*x310 + x190*x312;
    double x315 = x193*x310 + x195*x311 + x197*x312;
    double x316 = x186*x85;
    double x317 = x149*x186;
    double x318 = x18*x186;
    double x319 = x216*x316 + x218*x317 + x220*x318;
    double x320 = x133*x317 + x188*x316 + x190*x318;
    double x321 = x193*x316 + x195*x317 + x197*x318;
    double x322 = x100*x186;
    double x323 = x143*x186;
    double x324 = x186*x21;
    double x325 = x216*x322 + x218*x323 + x220*x324;
    double x326 = x133*x323 + x188*x322 + x190*x324;
    double x327 = x193*x322 + x195*x323 + x197*x324;
    double x328 = x186*x91;
    double x329 = x152*x186;
    double x330 = x186*x20;
    double x331 = x216*x328 + x218*x329 + x220*x330;
    double x332 = x133*x329 + x188*x328 + x190*x330;
    double x333 = x193*x328 + x195*x329 + x197*x330;
    double x334 = x186*x98;
    double x335 = x138*x186;
    double x336 = x15*x186;
    double x337 = x216*x334 + x218*x335 + x220*x336;
    double x338 = x133*x335 + x188*x334 + x190*x336;
    double x339 = x193*x334 + x195*x335 + x197*x336;
    double x340 = x186*x88;
    double x341 = x140*x186;
    double x342 = x186*x28;
    double x343 = x216*x340 + x218*x341 + x220*x342;
    double x344 = x133*x341 + x188*x340 + x190*x342;
    double x345 = x193*x340 + x195*x341 + x197*x342;
    double x346 = x186*x83;
    double x347 = x139*x186;
    double x348 = x186*x32;
    double x349 = x216*x346 + x218*x347 + x220*x348;
    double x350 = x133*x347 + x188*x346 + x190*x348;
    double x351 = x193*x346 + x195*x347 + x197*x348;
    double x352 = x186*x92;
    double x353 = x144*x186;
    double x354 = x186*x33;
    double x355 = x216*x352 + x218*x353 + x220*x354;
    double x356 = x133*x353 + x188*x352 + x190*x354;
    double x357 = x193*x352 + x195*x353 + x197*x354;
    double x358 = x186*x87;
    double x359 = x136*x186;
    double x360 = x186*x30;
    double x361 = x216*x358 + x218*x359 + x220*x360;
    double x362 = x133*x359 + x188*x358 + x190*x360;
    double x363 = x193*x358 + x195*x359 + x197*x360;
    
    res_0[0] = x185*(x215*x222 + x233*x234 + x235*x243);
    res_0[1] = x185*(x222*x244 + x234*x246 + x235*x245);
    res_0[2] = x185*(x222*x249 + x234*x247 + x235*x248);
    res_0[3] = x185*(x215*x253 + x233*x254 + x243*x255);
    res_0[4] = x185*(x244*x253 + x245*x255 + x246*x254);
    res_0[5] = x185*(x247*x254 + x248*x255 + x249*x253);
    res_0[6] = x185*(x215*x259 + x233*x260 + x243*x261);
    res_0[7] = x185*(x244*x259 + x245*x261 + x246*x260);
    res_0[8] = x185*(x247*x260 + x248*x261 + x249*x259);
    res_0[9] = x185*(x215*x265 + x233*x266 + x243*x267);
    res_0[10] = x185*(x244*x265 + x245*x267 + x246*x266);
    res_0[11] = x185*(x247*x266 + x248*x267 + x249*x265);
    res_0[12] = x185*(x215*x271 + x233*x272 + x243*x273);
    res_0[13] = x185*(x244*x271 + x245*x273 + x246*x272);
    res_0[14] = x185*(x247*x272 + x248*x273 + x249*x271);
    res_0[15] = x185*(x215*x277 + x233*x278 + x243*x279);
    res_0[16] = x185*(x244*x277 + x245*x279 + x246*x278);
    res_0[17] = x185*(x247*x278 + x248*x279 + x249*x277);
    res_0[18] = x185*(x215*x283 + x233*x284 + x243*x285);
    res_0[19] = x185*(x244*x283 + x245*x285 + x246*x284);
    res_0[20] = x185*(x247*x284 + x248*x285 + x249*x283);
    res_0[21] = x185*(x215*x289 + x233*x290 + x243*x291);
    res_0[22] = x185*(x244*x289 + x245*x291 + x246*x290);
    res_0[23] = x185*(x247*x290 + x248*x291 + x249*x289);
    res_0[24] = x185*(x215*x295 + x233*x296 + x243*x297);
    res_0[25] = x185*(x244*x295 + x245*x297 + x246*x296);
    res_0[26] = x185*(x247*x296 + x248*x297 + x249*x295);
    res_0[27] = x185*(x215*x301 + x233*x302 + x243*x303);
    res_0[28] = x185*(x244*x301 + x245*x303 + x246*x302);
    res_0[29] = x185*(x247*x302 + x248*x303 + x249*x301);
    res_0[30] = x185*(x215*x307 + x233*x308 + x243*x309);
    res_0[31] = x185*(x244*x307 + x245*x309 + x246*x308);
    res_0[32] = x185*(x247*x308 + x248*x309 + x249*x307);
    res_0[33] = x185*(x215*x313 + x233*x314 + x243*x315);
    res_0[34] = x185*(x244*x313 + x245*x315 + x246*x314);
    res_0[35] = x185*(x247*x314 + x248*x315 + x249*x313);
    res_0[36] = x185*(x215*x319 + x233*x320 + x243*x321);
    res_0[37] = x185*(x244*x319 + x245*x321 + x246*x320);
    res_0[38] = x185*(x247*x320 + x248*x321 + x249*x319);
    res_0[39] = x185*(x215*x325 + x233*x326 + x243*x327);
    res_0[40] = x185*(x244*x325 + x245*x327 + x246*x326);
    res_0[41] = x185*(x247*x326 + x248*x327 + x249*x325);
    res_0[42] = x185*(x215*x331 + x233*x332 + x243*x333);
    res_0[43] = x185*(x244*x331 + x245*x333 + x246*x332);
    res_0[44] = x185*(x247*x332 + x248*x333 + x249*x331);
    res_0[45] = x185*(x215*x337 + x233*x338 + x243*x339);
    res_0[46] = x185*(x244*x337 + x245*x339 + x246*x338);
    res_0[47] = x185*(x247*x338 + x248*x339 + x249*x337);
    res_0[48] = x185*(x215*x343 + x233*x344 + x243*x345);
    res_0[49] = x185*(x244*x343 + x245*x345 + x246*x344);
    res_0[50] = x185*(x247*x344 + x248*x345 + x249*x343);
    res_0[51] = x185*(x215*x349 + x233*x350 + x243*x351);
    res_0[52] = x185*(x244*x349 + x245*x351 + x246*x350);
    res_0[53] = x185*(x247*x350 + x248*x351 + x249*x349);
    res_0[54] = x185*(x215*x355 + x233*x356 + x243*x357);
    res_0[55] = x185*(x244*x355 + x245*x357 + x246*x356);
    res_0[56] = x185*(x247*x356 + x248*x357 + x249*x355);
    res_0[57] = x185*(x215*x361 + x233*x362 + x243*x363);
    res_0[58] = x185*(x244*x361 + x245*x363 + x246*x362);
    res_0[59] = x185*(x247*x362 + x248*x363 + x249*x361);
}

Conf_Forces_API void Integration_C3D20R_static_dbf(size_t num_elem,double Coords[][20][3],double Element_U[][20][3],double S[][8][6],
    double PENER[][8],double SENER[][8],double Conf_Force[][20][3]){
    
    //Inputs:
    //    Coords.........Element_Nodal_Coordinates [num_elem][numNodes,3]
    //    Element_U......Element_Nodal_Displacements [num_elem][numNodes,3]
    //    S..............Stress at integration points [num_elem][num_int_points,6]
    //    PENER..........Plastic strain energy [num_elem][num_int_points]
    //    SENER..........Elastic strain energy [num_elem][num_int_points]
    //    
    double weights[8]={1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,};
    double int_points[8][3]={
        {-0.5773502691896258,-0.5773502691896258,-0.5773502691896258,},
        {0.5773502691896258,-0.5773502691896258,-0.5773502691896258,},
        {-0.5773502691896258,0.5773502691896258,-0.5773502691896258,},
        {0.5773502691896258,0.5773502691896258,-0.5773502691896258,},
        {-0.5773502691896258,-0.5773502691896258,0.5773502691896258,},
        {0.5773502691896258,-0.5773502691896258,0.5773502691896258,},
        {-0.5773502691896258,0.5773502691896258,0.5773502691896258,},
        {0.5773502691896258,0.5773502691896258,0.5773502691896258,}};

    #pragma omp parallel for
    for (size_t i=0;i<num_elem;i++){
        double TMP[20][3];
        for (size_t j=0;j<8;j++){
            C3D20R_static_dbf((double*) int_points[j],(double*) Coords[i],(double*) Element_U[i],(double*) S[i][j],PENER[i][j],SENER[i][j],(double*) TMP);
            for (size_t num_node=0;num_node<20;num_node++){
                for (size_t ii=0;ii<3;ii++){
                    Conf_Force[i][num_node][ii]+=TMP[num_node][ii]*weights[j];
                }
            }
        }
    }
}

void C3D4_dynamic(double *rst,double *coord,double rho,double *Element_U,double *Element_V,double *Element_A,double *S,double PENER,double SENER,double *res_0){
    double r = rst[0];
    double s = rst[1];
    double t = rst[2];
    double *U = Element_U;
    double *V = Element_V;
    double *A = Element_A;
    double S11 = S[0];
    double S22 = S[1];
    double S33 = S[2];
    double S12 = S[3];
    double S13 = S[4];
    double S23 = S[5];
    
    double x0 = -coord[2] + coord[11];
    double x1 = -coord[0] + coord[3];
    double x2 = -coord[1] + coord[7];
    double x3 = x1*x2;
    double x4 = -coord[2] + coord[5];
    double x5 = -coord[0] + coord[6];
    double x6 = -coord[1] + coord[10];
    double x7 = x5*x6;
    double x8 = -coord[2] + coord[8];
    double x9 = -coord[0] + coord[9];
    double x10 = -coord[1] + coord[4];
    double x11 = x1*x6;
    double x12 = x10*x5;
    double x13 = x2*x9;
    double x14 = -x0*x12 + x0*x3 + x10*x8*x9 - x11*x8 - x13*x4 + x4*x7;
    double x15 = 1.0/x14;
    double x16 = x15*(-x12 + x3);
    double x17 = x15*(-x13 + x7);
    double x18 = x15*(x10*x9 - x11);
    double x19 = -x16 - x17 - x18;
    double x20 = -U[1] + U[4];
    double x21 = x15*(x0*x2 - x6*x8);
    double x22 = -U[1] + U[7];
    double x23 = x15*(-x0*x10 + x4*x6);
    double x24 = -U[1] + U[10];
    double x25 = x15*(x10*x8 - x2*x4);
    double x26 = x20*x21 + x22*x23 + x24*x25;
    double x27 = -U[2] + U[5];
    double x28 = x15*(-x0*x5 + x8*x9);
    double x29 = -U[2] + U[8];
    double x30 = x15*(x0*x1 - x4*x9);
    double x31 = -U[2] + U[11];
    double x32 = x15*(-x1*x8 + x4*x5);
    double x33 = x27*x28 + x29*x30 + x31*x32;
    double x34 = x21*x27 + x23*x29 + x25*x31;
    double x35 = x20*x28 + x22*x30 + x24*x32;
    double x36 = x35 + 1.0;
    double x37 = x26*x33 - x34*x36;
    double x38 = -U[0] + U[3];
    double x39 = -U[0] + U[6];
    double x40 = -U[0] + U[9];
    double x41 = x28*x38 + x30*x39 + x32*x40;
    double x42 = x21*x38 + x23*x39 + x25*x40;
    double x43 = x42 + 1.0;
    double x44 = -x33*x43 + x34*x41;
    double x45 = -x26*x41 + x36*x43;
    double x46 = S11*x37 + S12*x44 + S13*x45;
    double x47 = S12*x37 + S22*x44 + S23*x45;
    double x48 = S13*x37 + S23*x44 + S33*x45;
    double x49 = -x26*x47 - x34*x48 - x42*x46;
    double x50 = -x28 - x30 - x32;
    double x51 = x16*x24 + x17*x20 + x18*x22;
    double x52 = x16*x31 + x17*x27 + x18*x29;
    double x53 = x52 + 1.0;
    double x54 = -x26*x53 + x34*x51;
    double x55 = x16*x40 + x17*x38 + x18*x39;
    double x56 = -x34*x55 + x43*x53;
    double x57 = x26*x55 - x43*x51;
    double x58 = S11*x54 + S12*x56 + S13*x57;
    double x59 = S12*x54 + S22*x56 + S23*x57;
    double x60 = S13*x54 + S23*x56 + S33*x57;
    double x61 = -x26*x59 - x34*x60 - x42*x58;
    double x62 = -x21 - x23 - x25;
    double x63 = -x33*x51 + x36*x53;
    double x64 = x33*x55 - x41*x53;
    double x65 = -x36*x55 + x41*x51;
    double x66 = S11*x63 + S12*x64 + S13*x65;
    double x67 = S12*x63 + S22*x64 + S23*x65;
    double x68 = S13*x63 + S23*x64 + S33*x65;
    double x69 = -r - s - t + 1;
    double x70 = r*V[3] + s*V[6] + t*V[9] + x69*V[0];
    double x71 = r*V[4] + s*V[7] + t*V[10] + x69*V[1];
    double x72 = r*V[5] + s*V[8] + t*V[11] + x69*V[2];
    double x73 = -1.0*PENER - 1.0*SENER + 0.5*rho*(x70*x70 + x71*x71 + x72*x72);
    double x74 = -x26*x67 - x34*x68 - x42*x66 - x73;
    double x75 = r*A[3] + s*A[6] + t*A[9] + x69*A[0];
    double x76 = r*A[4] + s*A[7] + t*A[10] + x69*A[1];
    double x77 = r*A[5] + s*A[8] + t*A[11] + x69*A[2];
    double x78 = -V[0] + V[3];
    double x79 = -V[0] + V[6];
    double x80 = -V[0] + V[9];
    double x81 = -V[1] + V[4];
    double x82 = -V[1] + V[7];
    double x83 = -V[1] + V[10];
    double x84 = -V[2] + V[5];
    double x85 = -V[2] + V[8];
    double x86 = -V[2] + V[11];
    double x87 = x26*x76 + x34*x77 + x42*x75 + x70*(x21*x78 + x23*x79 + x25*x80) + x71*(x21*x81 + x23*x82 + x25*x83) + x72*(x21*x84 + x23*x85 + x25*x86);
    double x88 = rho*x69;
    double x89 = -x33*x48 - x35*x47 - x41*x46;
    double x90 = -x33*x60 - x35*x59 - x41*x58 - x73;
    double x91 = -x33*x68 - x35*x67 - x41*x66;
    double x92 = x33*x77 + x35*x76 + x41*x75 + x70*(x28*x78 + x30*x79 + x32*x80) + x71*(x28*x81 + x30*x82 + x32*x83) + x72*(x28*x84 + x30*x85 + x32*x86);
    double x93 = -x46*x55 - x47*x51 - x48*x52 - x73;
    double x94 = -x51*x59 - x52*x60 - x55*x58;
    double x95 = -x51*x67 - x52*x68 - x55*x66;
    double x96 = x51*x76 + x52*x77 + x55*x75 + x70*(x16*x80 + x17*x78 + x18*x79) + x71*(x16*x83 + x17*x81 + x18*x82) + x72*(x16*x86 + x17*x84 + x18*x85);
    double x97 = r*rho;
    double x98 = rho*s;
    double x99 = rho*t;
    
    res_0[0] = x14*(x19*x49 + x50*x61 + x62*x74 - x87*x88);
    res_0[1] = x14*(x19*x89 + x50*x90 + x62*x91 - x88*x92);
    res_0[2] = x14*(x19*x93 + x50*x94 + x62*x95 - x88*x96);
    res_0[3] = x14*(x17*x49 + x21*x74 + x28*x61 - x87*x97);
    res_0[4] = x14*(x17*x89 + x21*x91 + x28*x90 - x92*x97);
    res_0[5] = x14*(x17*x93 + x21*x95 + x28*x94 - x96*x97);
    res_0[6] = x14*(x18*x49 + x23*x74 + x30*x61 - x87*x98);
    res_0[7] = x14*(x18*x89 + x23*x91 + x30*x90 - x92*x98);
    res_0[8] = x14*(x18*x93 + x23*x95 + x30*x94 - x96*x98);
    res_0[9] = x14*(x16*x49 + x25*x74 + x32*x61 - x87*x99);
    res_0[10] = x14*(x16*x89 + x25*x91 + x32*x90 - x92*x99);
    res_0[11] = x14*(x16*x93 + x25*x95 + x32*x94 - x96*x99);
}

Conf_Forces_API void Integration_C3D4_dynamic(size_t num_elem,double Coords[][4][3],double *rho,double Element_U[][4][3],double Element_V[][4][3],double Element_A[][4][3],double S[][1][6],
    double PENER[][1],double SENER[][1],double Conf_Force[][4][3]){
    
    //Inputs:
    //    Coords.........Element_Nodal_Coordinates  [num_elem][numNodes,3]
    //    Element_U......Element_Nodal_Displacements[num_elem][numNodes,3]
    //    Element_V......Element_Nodal_Velocity     [num_elem][numNodes,3]
    //    Element_A......Element_Nodal_Acceleration [num_elem][numNodes,3]
    //    S..............Stress at integration points [num_elem][num_int_points,6]
    //    PENER..........Plastic strain energy [num_elem][num_int_points]
    //    SENER..........Elastic strain energy [num_elem][num_int_points]
    //    
    double weights[1]={0.16666666666666666,};
    double int_points[1][3]={
        {0.25,0.25,0.25,}};

    #pragma omp parallel for
    for (size_t i=0;i<num_elem;i++){
        double TMP[4][3];
        for (size_t j=0;j<1;j++){
            C3D4_dynamic((double*) int_points[j],(double*) Coords[i],rho[i],(double*) Element_U[i],(double*) Element_V[i],(double*) Element_A[i],(double*) S[i][j],PENER[i][j],SENER[i][j],(double*) TMP);
            for (size_t num_node=0;num_node<4;num_node++){
                for (size_t ii=0;ii<3;ii++){
                    Conf_Force[i][num_node][ii]+=TMP[num_node][ii]*weights[j];
                }
            }
        }
    }
}

void C3D4_static_mbf(double *rst,double *coord,double *Element_U,double *S,double PENER,double SENER,double *res_0){
    double r = rst[0];
    double s = rst[1];
    double t = rst[2];
    double *U = Element_U;
    double S11 = S[0];
    double S22 = S[1];
    double S33 = S[2];
    double S12 = S[3];
    double S13 = S[4];
    double S23 = S[5];
    
    double x0 = -coord[2] + coord[11];
    double x1 = -coord[0] + coord[3];
    double x2 = -coord[1] + coord[7];
    double x3 = x1*x2;
    double x4 = -coord[2] + coord[5];
    double x5 = -coord[0] + coord[6];
    double x6 = -coord[1] + coord[10];
    double x7 = x5*x6;
    double x8 = -coord[2] + coord[8];
    double x9 = -coord[0] + coord[9];
    double x10 = -coord[1] + coord[4];
    double x11 = x1*x6;
    double x12 = x10*x5;
    double x13 = x2*x9;
    double x14 = -x0*x12 + x0*x3 + x10*x8*x9 - x11*x8 - x13*x4 + x4*x7;
    double x15 = 1.0/x14;
    double x16 = x15*(-x12 + x3);
    double x17 = x15*(-x13 + x7);
    double x18 = x15*(x10*x9 - x11);
    double x19 = -x16 - x17 - x18;
    double x20 = -U[1] + U[4];
    double x21 = x15*(x0*x2 - x6*x8);
    double x22 = -U[1] + U[7];
    double x23 = x15*(-x0*x10 + x4*x6);
    double x24 = -U[1] + U[10];
    double x25 = x15*(x10*x8 - x2*x4);
    double x26 = x20*x21 + x22*x23 + x24*x25;
    double x27 = -U[2] + U[5];
    double x28 = x15*(-x0*x5 + x8*x9);
    double x29 = -U[2] + U[8];
    double x30 = x15*(x0*x1 - x4*x9);
    double x31 = -U[2] + U[11];
    double x32 = x15*(-x1*x8 + x4*x5);
    double x33 = x27*x28 + x29*x30 + x31*x32;
    double x34 = x21*x27 + x23*x29 + x25*x31;
    double x35 = x20*x28 + x22*x30 + x24*x32 + 1.0;
    double x36 = x26*x33 - x34*x35;
    double x37 = -U[0] + U[3];
    double x38 = -U[0] + U[6];
    double x39 = -U[0] + U[9];
    double x40 = x28*x37 + x30*x38 + x32*x39;
    double x41 = x21*x37 + x23*x38 + x25*x39 + 1.0;
    double x42 = -x33*x41 + x34*x40;
    double x43 = -x26*x40 + x35*x41;
    double x44 = S11*x36 + S12*x42 + S13*x43;
    double x45 = S12*x36 + S22*x42 + S23*x43;
    double x46 = S13*x36 + S23*x42 + S33*x43;
    double x47 = -x26*x45 - x34*x46 - x41*x44;
    double x48 = -x28 - x30 - x32;
    double x49 = x16*x24 + x17*x20 + x18*x22;
    double x50 = x16*x31 + x17*x27 + x18*x29 + 1.0;
    double x51 = -x26*x50 + x34*x49;
    double x52 = x16*x39 + x17*x37 + x18*x38;
    double x53 = -x34*x52 + x41*x50;
    double x54 = x26*x52 - x41*x49;
    double x55 = S11*x51 + S12*x53 + S13*x54;
    double x56 = S12*x51 + S22*x53 + S23*x54;
    double x57 = S13*x51 + S23*x53 + S33*x54;
    double x58 = -x26*x56 - x34*x57 - x41*x55;
    double x59 = -x21 - x23 - x25;
    double x60 = -x33*x49 + x35*x50;
    double x61 = x33*x52 - x40*x50;
    double x62 = -x35*x52 + x40*x49;
    double x63 = S11*x60 + S12*x61 + S13*x62;
    double x64 = S12*x60 + S22*x61 + S23*x62;
    double x65 = S13*x60 + S23*x61 + S33*x62;
    double x66 = -1.0*PENER - 1.0*SENER;
    double x67 = -x26*x64 - x34*x65 - x41*x63 - x66;
    double x68 = -x33*x46 - x35*x45 - x40*x44;
    double x69 = -x33*x57 - x35*x56 - x40*x55 - x66;
    double x70 = -x33*x65 - x35*x64 - x40*x63;
    double x71 = -x44*x52 - x45*x49 - x46*x50 - x66;
    double x72 = -x49*x56 - x50*x57 - x52*x55;
    double x73 = -x49*x64 - x50*x65 - x52*x63;
    
    res_0[0] = x14*(x19*x47 + x48*x58 + x59*x67);
    res_0[1] = x14*(x19*x68 + x48*x69 + x59*x70);
    res_0[2] = x14*(x19*x71 + x48*x72 + x59*x73);
    res_0[3] = x14*(x17*x47 + x21*x67 + x28*x58);
    res_0[4] = x14*(x17*x68 + x21*x70 + x28*x69);
    res_0[5] = x14*(x17*x71 + x21*x73 + x28*x72);
    res_0[6] = x14*(x18*x47 + x23*x67 + x30*x58);
    res_0[7] = x14*(x18*x68 + x23*x70 + x30*x69);
    res_0[8] = x14*(x18*x71 + x23*x73 + x30*x72);
    res_0[9] = x14*(x16*x47 + x25*x67 + x32*x58);
    res_0[10] = x14*(x16*x68 + x25*x70 + x32*x69);
    res_0[11] = x14*(x16*x71 + x25*x73 + x32*x72);
}

Conf_Forces_API void Integration_C3D4_static_mbf(size_t num_elem,double Coords[][4][3],double Element_U[][4][3],double S[][1][6],
    double PENER[][1],double SENER[][1],double Conf_Force[][4][3]){
    
    //Inputs:
    //    Coords.........Element_Nodal_Coordinates [num_elem][numNodes,3]
    //    Element_U......Element_Nodal_Displacements [num_elem][numNodes,3]
    //    S..............Stress at integration points [num_elem][num_int_points,6]
    //    PENER..........Plastic strain energy [num_elem][num_int_points]
    //    SENER..........Elastic strain energy [num_elem][num_int_points]
    //    
    double weights[1]={0.16666666666666666,};
    double int_points[1][3]={
        {0.25,0.25,0.25,}};

    #pragma omp parallel for
    for (size_t i=0;i<num_elem;i++){
        double TMP[4][3];
        for (size_t j=0;j<1;j++){
            C3D4_static_mbf((double*) int_points[j],(double*) Coords[i],(double*) Element_U[i],(double*) S[i][j],PENER[i][j],SENER[i][j],(double*) TMP);
            for (size_t num_node=0;num_node<4;num_node++){
                for (size_t ii=0;ii<3;ii++){
                    Conf_Force[i][num_node][ii]+=TMP[num_node][ii]*weights[j];
                }
            }
        }
    }
}

void C3D4_static_dbf(double *rst,double *coord,double *Element_U,double *S,double PENER,double SENER,double *res_0){
    double r = rst[0];
    double s = rst[1];
    double t = rst[2];
    double *U = Element_U;
    double S11 = S[0];
    double S22 = S[1];
    double S33 = S[2];
    double S12 = S[3];
    double S13 = S[4];
    double S23 = S[5];
    
    double x0 = -coord[2] + coord[11];
    double x1 = -coord[0] + coord[3];
    double x2 = -coord[1] + coord[7];
    double x3 = x1*x2;
    double x4 = -coord[2] + coord[5];
    double x5 = -coord[0] + coord[6];
    double x6 = -coord[1] + coord[10];
    double x7 = x5*x6;
    double x8 = -coord[2] + coord[8];
    double x9 = -coord[0] + coord[9];
    double x10 = -coord[1] + coord[4];
    double x11 = x1*x6;
    double x12 = x10*x5;
    double x13 = x2*x9;
    double x14 = -x0*x12 + x0*x3 + x10*x8*x9 - x11*x8 - x13*x4 + x4*x7;
    double x15 = 1.0/x14;
    double x16 = x15*(-x12 + x3);
    double x17 = x15*(-x13 + x7);
    double x18 = x15*(x10*x9 - x11);
    double x19 = -x16 - x17 - x18;
    double x20 = -U[1] + U[4];
    double x21 = x15*(x0*x2 - x6*x8);
    double x22 = -U[1] + U[7];
    double x23 = x15*(-x0*x10 + x4*x6);
    double x24 = -U[1] + U[10];
    double x25 = x15*(x10*x8 - x2*x4);
    double x26 = x20*x21 + x22*x23 + x24*x25;
    double x27 = -U[2] + U[5];
    double x28 = x15*(-x0*x5 + x8*x9);
    double x29 = -U[2] + U[8];
    double x30 = x15*(x0*x1 - x4*x9);
    double x31 = -U[2] + U[11];
    double x32 = x15*(-x1*x8 + x4*x5);
    double x33 = x27*x28 + x29*x30 + x31*x32;
    double x34 = x21*x27 + x23*x29 + x25*x31;
    double x35 = x20*x28 + x22*x30 + x24*x32;
    double x36 = x35 + 1.0;
    double x37 = x26*x33 - x34*x36;
    double x38 = -U[0] + U[3];
    double x39 = -U[0] + U[6];
    double x40 = -U[0] + U[9];
    double x41 = x28*x38 + x30*x39 + x32*x40;
    double x42 = x21*x38 + x23*x39 + x25*x40;
    double x43 = x42 + 1.0;
    double x44 = -x33*x43 + x34*x41;
    double x45 = -x26*x41 + x36*x43;
    double x46 = S11*x37 + S12*x44 + S13*x45;
    double x47 = S12*x37 + S22*x44 + S23*x45;
    double x48 = S13*x37 + S23*x44 + S33*x45;
    double x49 = -x26*x47 - x34*x48 - x42*x46;
    double x50 = -x28 - x30 - x32;
    double x51 = x16*x24 + x17*x20 + x18*x22;
    double x52 = x16*x31 + x17*x27 + x18*x29;
    double x53 = x52 + 1.0;
    double x54 = -x26*x53 + x34*x51;
    double x55 = x16*x40 + x17*x38 + x18*x39;
    double x56 = -x34*x55 + x43*x53;
    double x57 = x26*x55 - x43*x51;
    double x58 = S11*x54 + S12*x56 + S13*x57;
    double x59 = S12*x54 + S22*x56 + S23*x57;
    double x60 = S13*x54 + S23*x56 + S33*x57;
    double x61 = -x26*x59 - x34*x60 - x42*x58;
    double x62 = -x21 - x23 - x25;
    double x63 = -x33*x51 + x36*x53;
    double x64 = x33*x55 - x41*x53;
    double x65 = -x36*x55 + x41*x51;
    double x66 = S11*x63 + S12*x64 + S13*x65;
    double x67 = S12*x63 + S22*x64 + S23*x65;
    double x68 = S13*x63 + S23*x64 + S33*x65;
    double x69 = -1.0*PENER - 1.0*SENER;
    double x70 = -x26*x67 - x34*x68 - x42*x66 - x69;
    double x71 = -x33*x48 - x35*x47 - x41*x46;
    double x72 = -x33*x60 - x35*x59 - x41*x58 - x69;
    double x73 = -x33*x68 - x35*x67 - x41*x66;
    double x74 = -x46*x55 - x47*x51 - x48*x52 - x69;
    double x75 = -x51*x59 - x52*x60 - x55*x58;
    double x76 = -x51*x67 - x52*x68 - x55*x66;
    
    res_0[0] = x14*(x19*x49 + x50*x61 + x62*x70);
    res_0[1] = x14*(x19*x71 + x50*x72 + x62*x73);
    res_0[2] = x14*(x19*x74 + x50*x75 + x62*x76);
    res_0[3] = x14*(x17*x49 + x21*x70 + x28*x61);
    res_0[4] = x14*(x17*x71 + x21*x73 + x28*x72);
    res_0[5] = x14*(x17*x74 + x21*x76 + x28*x75);
    res_0[6] = x14*(x18*x49 + x23*x70 + x30*x61);
    res_0[7] = x14*(x18*x71 + x23*x73 + x30*x72);
    res_0[8] = x14*(x18*x74 + x23*x76 + x30*x75);
    res_0[9] = x14*(x16*x49 + x25*x70 + x32*x61);
    res_0[10] = x14*(x16*x71 + x25*x73 + x32*x72);
    res_0[11] = x14*(x16*x74 + x25*x76 + x32*x75);
}

Conf_Forces_API void Integration_C3D4_static_dbf(size_t num_elem,double Coords[][4][3],double Element_U[][4][3],double S[][1][6],
    double PENER[][1],double SENER[][1],double Conf_Force[][4][3]){
    
    //Inputs:
    //    Coords.........Element_Nodal_Coordinates [num_elem][numNodes,3]
    //    Element_U......Element_Nodal_Displacements [num_elem][numNodes,3]
    //    S..............Stress at integration points [num_elem][num_int_points,6]
    //    PENER..........Plastic strain energy [num_elem][num_int_points]
    //    SENER..........Elastic strain energy [num_elem][num_int_points]
    //    
    double weights[1]={0.16666666666666666,};
    double int_points[1][3]={
        {0.25,0.25,0.25,}};

    #pragma omp parallel for
    for (size_t i=0;i<num_elem;i++){
        double TMP[4][3];
        for (size_t j=0;j<1;j++){
            C3D4_static_dbf((double*) int_points[j],(double*) Coords[i],(double*) Element_U[i],(double*) S[i][j],PENER[i][j],SENER[i][j],(double*) TMP);
            for (size_t num_node=0;num_node<4;num_node++){
                for (size_t ii=0;ii<3;ii++){
                    Conf_Force[i][num_node][ii]+=TMP[num_node][ii]*weights[j];
                }
            }
        }
    }
}

void C3D10_dynamic(double *rst,double *coord,double rho,double *Element_U,double *Element_V,double *Element_A,double *S,double PENER,double SENER,double *res_0){
    double r = rst[0];
    double s = rst[1];
    double t = rst[2];
    double *U = Element_U;
    double *V = Element_V;
    double *A = Element_A;
    double S11 = S[0];
    double S22 = S[1];
    double S33 = S[2];
    double S12 = S[3];
    double S13 = S[4];
    double S23 = S[5];
    
    double x0 = 4*r;
    double x1 = x0 - 1;
    double x2 = 4*s;
    double x3 = 4*t;
    double x4 = x2 + x3;
    double x5 = -8*r - x4 + 4;
    double x6 = -x2*coord[19];
    double x7 = x0 + x4 - 3;
    double x8 = x7*coord[1];
    double x9 = -x3*coord[22] + x8;
    double x10 = x1*coord[4] + x2*coord[16] + x3*coord[25] + x5*coord[13] + x6 + x9;
    double x11 = x2 - 1;
    double x12 = x7*coord[0];
    double x13 = x0 - 4;
    double x14 = -8*s - x13 - x3;
    double x15 = -x0*coord[12];
    double x16 = -x3*coord[21];
    double x17 = x0*coord[15] + x11*coord[6] + x12 + x14*coord[18] + x15 + x16 + x3*coord[27];
    double x18 = x10*x17;
    double x19 = x12 - x2*coord[18];
    double x20 = x1*coord[3] + x16 + x19 + x2*coord[15] + x3*coord[24] + x5*coord[12];
    double x21 = -x0*coord[13];
    double x22 = x0*coord[16] + x11*coord[7] + x14*coord[19] + x21 + x3*coord[28] + x9;
    double x23 = -x18 + x20*x22;
    double x24 = x3 - 1;
    double x25 = -8*t - x13 - x2;
    double x26 = -x2*coord[20];
    double x27 = x7*coord[2];
    double x28 = -x0*coord[14] + x27;
    double x29 = x0*coord[26] + x2*coord[29] + x24*coord[11] + x25*coord[23] + x26 + x28;
    double x30 = -x3*coord[23];
    double x31 = x0*coord[17] + x11*coord[8] + x14*coord[20] + x28 + x3*coord[29] + x30;
    double x32 = x0*coord[24] + x15 + x19 + x2*coord[27] + x24*coord[9] + x25*coord[21];
    double x33 = x10*x32;
    double x34 = x1*coord[5] + x2*coord[17] + x26 + x27 + x3*coord[26] + x30 + x5*coord[14];
    double x35 = x0*coord[25] + x2*coord[28] + x21 + x24*coord[10] + x25*coord[22] + x6 + x8;
    double x36 = x17*x35;
    double x37 = x20*x35;
    double x38 = x22*x32;
    double x39 = -x18*x29 + x20*x22*x29 + x31*x33 - x31*x37 + x34*x36 - x34*x38;
    double x40 = 1.0/x39;
    double x41 = -x2*U[19];
    double x42 = x7*U[1];
    double x43 = -x0*U[13] + x42;
    double x44 = x40*(x0*U[25] + x2*U[28] + x24*U[10] + x25*U[22] + x41 + x43);
    double x45 = x33 - x37;
    double x46 = -x3*U[22];
    double x47 = x40*(x0*U[16] + x11*U[7] + x14*U[19] + x3*U[28] + x43 + x46);
    double x48 = x36 - x38;
    double x49 = x40*(x1*U[4] + x2*U[16] + x3*U[25] + x41 + x42 + x46 + x5*U[13]);
    double x50 = x23*x44 + x45*x47 + x48*x49;
    double x51 = x10*x31 - x22*x34;
    double x52 = -x2*U[20];
    double x53 = x7*U[2];
    double x54 = -x0*U[14] + x53;
    double x55 = x40*(x0*U[26] + x2*U[29] + x24*U[11] + x25*U[23] + x52 + x54);
    double x56 = -x10*x29 + x34*x35;
    double x57 = -x3*U[23];
    double x58 = x40*(x0*U[17] + x11*U[8] + x14*U[20] + x3*U[29] + x54 + x57);
    double x59 = x22*x29 - x31*x35;
    double x60 = x40*(x1*U[5] + x2*U[17] + x3*U[26] + x5*U[14] + x52 + x53 + x57);
    double x61 = x51*x55 + x56*x58 + x59*x60;
    double x62 = x44*x51 + x47*x56 + x49*x59;
    double x63 = x23*x55 + x45*x58 + x48*x60;
    double x64 = x63 + 1.0;
    double x65 = x50*x61 - x62*x64;
    double x66 = -x2*U[18];
    double x67 = x7*U[0];
    double x68 = -x0*U[12] + x67;
    double x69 = x40*(x0*U[24] + x2*U[27] + x24*U[9] + x25*U[21] + x66 + x68);
    double x70 = -x3*U[21];
    double x71 = x40*(x0*U[15] + x11*U[6] + x14*U[18] + x3*U[27] + x68 + x70);
    double x72 = x40*(x1*U[3] + x2*U[15] + x3*U[24] + x5*U[12] + x66 + x67 + x70);
    double x73 = x23*x69 + x45*x71 + x48*x72;
    double x74 = x51*x69 + x56*x71 + x59*x72;
    double x75 = x74 + 1.0;
    double x76 = -x61*x73 + x64*x75;
    double x77 = -x50*x75 + x62*x73;
    double x78 = S11*x65 + S12*x76 + S13*x77;
    double x79 = S12*x65 + S22*x76 + S23*x77;
    double x80 = S13*x65 + S23*x76 + S33*x77;
    double x81 = -x61*x80 - x62*x79 - x74*x78;
    double x82 = x20*x29 - x32*x34;
    double x83 = x40*x7;
    double x84 = x17*x34 - x20*x31;
    double x85 = -x17*x29 + x31*x32;
    double x86 = x82*x83 + x83*x84 + x83*x85;
    double x87 = x55*x84 + x58*x82 + x60*x85;
    double x88 = x44*x84 + x47*x82 + x49*x85;
    double x89 = x88 + 1.0;
    double x90 = -x61*x89 + x62*x87;
    double x91 = x69*x84 + x71*x82 + x72*x85;
    double x92 = x61*x91 - x75*x87;
    double x93 = -x62*x91 + x75*x89;
    double x94 = S11*x90 + S12*x92 + S13*x93;
    double x95 = S12*x90 + S22*x92 + S23*x93;
    double x96 = S13*x90 + S23*x92 + S33*x93;
    double x97 = -x61*x96 - x62*x95 - x74*x94;
    double x98 = x23*x83 + x45*x83 + x48*x83;
    double x99 = x51*x83 + x56*x83 + x59*x83;
    double x100 = -x50*x87 + x64*x89;
    double x101 = -x64*x91 + x73*x87;
    double x102 = x50*x91 - x73*x89;
    double x103 = S11*x100 + S12*x101 + S13*x102;
    double x104 = S12*x100 + S22*x101 + S23*x102;
    double x105 = S13*x100 + S23*x101 + S33*x102;
    double x106 = r*r;
    double x107 = 2*x106;
    double x108 = -r + x107;
    double x109 = s*s;
    double x110 = 2*x109;
    double x111 = -s + x110;
    double x112 = t*t;
    double x113 = 2*x112;
    double x114 = -t + x113;
    double x115 = s*x0;
    double x116 = t*x0;
    double x117 = x115 + x116;
    double x118 = x0 - 4*x106 - x117;
    double x119 = t*x2;
    double x120 = -4*x109 - x115 - x119 + x2;
    double x121 = -4*x112 - x116 - x119 + x3;
    double x122 = -3*r - 3*s - 3*t + x107 + x110 + x113 + x117 + x119 + 1;
    double x123 = x108*V[3] + x111*V[6] + x114*V[9] + x115*V[15] + x116*V[24] + x118*V[12] + x119*V[27] + x120*V[18] + x121*V[21] + x122*V[0];
    double x124 = x108*V[4] + x111*V[7] + x114*V[10] + x115*V[16] + x116*V[25] + x118*V[13] + x119*V[28] + x120*V[19] + x121*V[22] + x122*V[1];
    double x125 = x108*V[5] + x111*V[8] + x114*V[11] + x115*V[17] + x116*V[26] + x118*V[14] + x119*V[29] + x120*V[20] + x121*V[23] + x122*V[2];
    double x126 = -1.0*PENER - 1.0*SENER + 0.5*rho*(x123*x123 + x124*x124 + x125*x125);
    double x127 = -x103*x74 - x104*x62 - x105*x61 - x126;
    double x128 = x108*A[3] + x111*A[6] + x114*A[9] + x115*A[15] + x116*A[24] + x118*A[12] + x119*A[27] + x120*A[18] + x121*A[21] + x122*A[0];
    double x129 = x108*A[4] + x111*A[7] + x114*A[10] + x115*A[16] + x116*A[25] + x118*A[13] + x119*A[28] + x120*A[19] + x121*A[22] + x122*A[1];
    double x130 = x108*A[5] + x111*A[8] + x114*A[11] + x115*A[17] + x116*A[26] + x118*A[14] + x119*A[29] + x120*A[20] + x121*A[23] + x122*A[2];
    double x131 = -x2*V[18];
    double x132 = x7*V[0];
    double x133 = -x0*V[12] + x132;
    double x134 = x0*V[24] + x131 + x133 + x2*V[27] + x24*V[9] + x25*V[21];
    double x135 = x40*x51;
    double x136 = -x3*V[21];
    double x137 = x0*V[15] + x11*V[6] + x133 + x136 + x14*V[18] + x3*V[27];
    double x138 = x40*x56;
    double x139 = x1*V[3] + x131 + x132 + x136 + x2*V[15] + x3*V[24] + x5*V[12];
    double x140 = x40*x59;
    double x141 = -x2*V[19];
    double x142 = x7*V[1];
    double x143 = -x0*V[13] + x142;
    double x144 = x0*V[25] + x141 + x143 + x2*V[28] + x24*V[10] + x25*V[22];
    double x145 = -x3*V[22];
    double x146 = x0*V[16] + x11*V[7] + x14*V[19] + x143 + x145 + x3*V[28];
    double x147 = x1*V[4] + x141 + x142 + x145 + x2*V[16] + x3*V[25] + x5*V[13];
    double x148 = -x2*V[20];
    double x149 = x7*V[2];
    double x150 = -x0*V[14] + x149;
    double x151 = x0*V[26] + x148 + x150 + x2*V[29] + x24*V[11] + x25*V[23];
    double x152 = -x3*V[23];
    double x153 = x0*V[17] + x11*V[8] + x14*V[20] + x150 + x152 + x3*V[29];
    double x154 = x1*V[5] + x148 + x149 + x152 + x2*V[17] + x3*V[26] + x5*V[14];
    double x155 = x123*(x134*x135 + x137*x138 + x139*x140) + x124*(x135*x144 + x138*x146 + x140*x147) + x125*(x135*x151 + x138*x153 + x140*x154) + x128*x74 + x129*x62 + x130*x61;
    double x156 = rho*x122;
    double x157 = -x87*x96 - x88*x95 - x91*x94;
    double x158 = -x103*x91 - x104*x88 - x105*x87;
    double x159 = -x126 - x78*x91 - x79*x88 - x80*x87;
    double x160 = x40*x82;
    double x161 = x40*x84;
    double x162 = x40*x85;
    double x163 = x123*(x134*x161 + x137*x160 + x139*x162) + x124*(x144*x161 + x146*x160 + x147*x162) + x125*(x151*x161 + x153*x160 + x154*x162) + x128*x91 + x129*x88 + x130*x87;
    double x164 = -x50*x79 - x63*x80 - x73*x78;
    double x165 = -x103*x73 - x104*x50 - x105*x63;
    double x166 = -x126 - x50*x95 - x63*x96 - x73*x94;
    double x167 = x23*x40;
    double x168 = x40*x45;
    double x169 = x40*x48;
    double x170 = x123*(x134*x167 + x137*x168 + x139*x169) + x124*(x144*x167 + x146*x168 + x147*x169) + x125*(x151*x167 + x153*x168 + x154*x169) + x128*x73 + x129*x50 + x130*x63;
    double x171 = rho*x108;
    double x172 = x1*x169;
    double x173 = x1*x140;
    double x174 = x1*x162;
    double x175 = rho*x111;
    double x176 = x11*x160;
    double x177 = x11*x168;
    double x178 = x11*x138;
    double x179 = rho*x114;
    double x180 = x167*x24;
    double x181 = x135*x24;
    double x182 = x161*x24;
    double x183 = x0*x160;
    double x184 = x0*x161;
    double x185 = -x183 - x184 + x40*x5*x85;
    double x186 = x0*x167;
    double x187 = x0*x168;
    double x188 = -x186 - x187 + x40*x48*x5;
    double x189 = x0*x135;
    double x190 = x0*x138;
    double x191 = -x189 - x190 + x40*x5*x59;
    double x192 = rho*x118;
    double x193 = x162*x2;
    double x194 = x183 + x193;
    double x195 = x169*x2;
    double x196 = x187 + x195;
    double x197 = x140*x2;
    double x198 = x190 + x197;
    double x199 = rho*x115;
    double x200 = x161*x2;
    double x201 = x14*x40*x82 - x193 - x200;
    double x202 = x167*x2;
    double x203 = x14*x40*x45 - x195 - x202;
    double x204 = x135*x2;
    double x205 = x14*x40*x56 - x197 - x204;
    double x206 = rho*x120;
    double x207 = x160*x3;
    double x208 = x162*x3;
    double x209 = -x207 - x208 + x25*x40*x84;
    double x210 = x168*x3;
    double x211 = x169*x3;
    double x212 = -x210 - x211 + x23*x25*x40;
    double x213 = x138*x3;
    double x214 = x140*x3;
    double x215 = -x213 - x214 + x25*x40*x51;
    double x216 = rho*x121;
    double x217 = x186 + x211;
    double x218 = x189 + x214;
    double x219 = x184 + x208;
    double x220 = rho*x116;
    double x221 = x202 + x210;
    double x222 = x204 + x213;
    double x223 = x200 + x207;
    double x224 = rho*x119;
    
    res_0[0] = x39*(x127*x99 - x155*x156 + x81*x86 + x97*x98);
    res_0[1] = x39*(-x156*x163 + x157*x98 + x158*x99 + x159*x86);
    res_0[2] = x39*(-x156*x170 + x164*x86 + x165*x99 + x166*x98);
    res_0[3] = x39*(x127*x173 - x155*x171 + x172*x97 + x174*x81);
    res_0[4] = x39*(x157*x172 + x158*x173 + x159*x174 - x163*x171);
    res_0[5] = x39*(x164*x174 + x165*x173 + x166*x172 - x170*x171);
    res_0[6] = x39*(x127*x178 - x155*x175 + x176*x81 + x177*x97);
    res_0[7] = x39*(x157*x177 + x158*x178 + x159*x176 - x163*x175);
    res_0[8] = x39*(x164*x176 + x165*x178 + x166*x177 - x170*x175);
    res_0[9] = x39*(x127*x181 - x155*x179 + x180*x97 + x182*x81);
    res_0[10] = x39*(x157*x180 + x158*x181 + x159*x182 - x163*x179);
    res_0[11] = x39*(x164*x182 + x165*x181 + x166*x180 - x170*x179);
    res_0[12] = x39*(x127*x191 - x155*x192 + x185*x81 + x188*x97);
    res_0[13] = x39*(x157*x188 + x158*x191 + x159*x185 - x163*x192);
    res_0[14] = x39*(x164*x185 + x165*x191 + x166*x188 - x170*x192);
    res_0[15] = x39*(x127*x198 - x155*x199 + x194*x81 + x196*x97);
    res_0[16] = x39*(x157*x196 + x158*x198 + x159*x194 - x163*x199);
    res_0[17] = x39*(x164*x194 + x165*x198 + x166*x196 - x170*x199);
    res_0[18] = x39*(x127*x205 - x155*x206 + x201*x81 + x203*x97);
    res_0[19] = x39*(x157*x203 + x158*x205 + x159*x201 - x163*x206);
    res_0[20] = x39*(x164*x201 + x165*x205 + x166*x203 - x170*x206);
    res_0[21] = x39*(x127*x215 - x155*x216 + x209*x81 + x212*x97);
    res_0[22] = x39*(x157*x212 + x158*x215 + x159*x209 - x163*x216);
    res_0[23] = x39*(x164*x209 + x165*x215 + x166*x212 - x170*x216);
    res_0[24] = x39*(x127*x218 - x155*x220 + x217*x97 + x219*x81);
    res_0[25] = x39*(x157*x217 + x158*x218 + x159*x219 - x163*x220);
    res_0[26] = x39*(x164*x219 + x165*x218 + x166*x217 - x170*x220);
    res_0[27] = x39*(x127*x222 - x155*x224 + x221*x97 + x223*x81);
    res_0[28] = x39*(x157*x221 + x158*x222 + x159*x223 - x163*x224);
    res_0[29] = x39*(x164*x223 + x165*x222 + x166*x221 - x170*x224);
}

Conf_Forces_API void Integration_C3D10_dynamic(size_t num_elem,double Coords[][10][3],double *rho,double Element_U[][10][3],double Element_V[][10][3],double Element_A[][10][3],double S[][4][6],
    double PENER[][4],double SENER[][4],double Conf_Force[][10][3]){
    
    //Inputs:
    //    Coords.........Element_Nodal_Coordinates  [num_elem][numNodes,3]
    //    Element_U......Element_Nodal_Displacements[num_elem][numNodes,3]
    //    Element_V......Element_Nodal_Velocity     [num_elem][numNodes,3]
    //    Element_A......Element_Nodal_Acceleration [num_elem][numNodes,3]
    //    S..............Stress at integration points [num_elem][num_int_points,6]
    //    PENER..........Plastic strain energy [num_elem][num_int_points]
    //    SENER..........Elastic strain energy [num_elem][num_int_points]
    //    
    double weights[4]={0.041666666666666664,0.041666666666666664,0.041666666666666664,0.041666666666666664,};
    double int_points[4][3]={
        {0.1381966011250105,0.1381966011250105,0.1381966011250105,},
        {0.5854101966249685,0.1381966011250105,0.1381966011250105,},
        {0.1381966011250105,0.5854101966249685,0.1381966011250105,},
        {0.1381966011250105,0.1381966011250105,0.5854101966249685,}};

    #pragma omp parallel for
    for (size_t i=0;i<num_elem;i++){
        double TMP[10][3];
        for (size_t j=0;j<4;j++){
            C3D10_dynamic((double*) int_points[j],(double*) Coords[i],rho[i],(double*) Element_U[i],(double*) Element_V[i],(double*) Element_A[i],(double*) S[i][j],PENER[i][j],SENER[i][j],(double*) TMP);
            for (size_t num_node=0;num_node<10;num_node++){
                for (size_t ii=0;ii<3;ii++){
                    Conf_Force[i][num_node][ii]+=TMP[num_node][ii]*weights[j];
                }
            }
        }
    }
}

void C3D10_static_mbf(double *rst,double *coord,double *Element_U,double *S,double PENER,double SENER,double *res_0){
    double r = rst[0];
    double s = rst[1];
    double t = rst[2];
    double *U = Element_U;
    double S11 = S[0];
    double S22 = S[1];
    double S33 = S[2];
    double S12 = S[3];
    double S13 = S[4];
    double S23 = S[5];
    
    double x0 = 4*r;
    double x1 = x0 - 1;
    double x2 = 4*s;
    double x3 = 4*t;
    double x4 = x2 + x3;
    double x5 = -8*r - x4 + 4;
    double x6 = -x2*coord[19];
    double x7 = x0 + x4 - 3;
    double x8 = x7*coord[1];
    double x9 = -x3*coord[22] + x8;
    double x10 = x1*coord[4] + x2*coord[16] + x3*coord[25] + x5*coord[13] + x6 + x9;
    double x11 = x2 - 1;
    double x12 = x7*coord[0];
    double x13 = x0 - 4;
    double x14 = -8*s - x13 - x3;
    double x15 = -x0*coord[12];
    double x16 = -x3*coord[21];
    double x17 = x0*coord[15] + x11*coord[6] + x12 + x14*coord[18] + x15 + x16 + x3*coord[27];
    double x18 = x10*x17;
    double x19 = x12 - x2*coord[18];
    double x20 = x1*coord[3] + x16 + x19 + x2*coord[15] + x3*coord[24] + x5*coord[12];
    double x21 = -x0*coord[13];
    double x22 = x0*coord[16] + x11*coord[7] + x14*coord[19] + x21 + x3*coord[28] + x9;
    double x23 = -x18 + x20*x22;
    double x24 = x3 - 1;
    double x25 = -8*t - x13 - x2;
    double x26 = -x2*coord[20];
    double x27 = x7*coord[2];
    double x28 = -x0*coord[14] + x27;
    double x29 = x0*coord[26] + x2*coord[29] + x24*coord[11] + x25*coord[23] + x26 + x28;
    double x30 = -x3*coord[23];
    double x31 = x0*coord[17] + x11*coord[8] + x14*coord[20] + x28 + x3*coord[29] + x30;
    double x32 = x0*coord[24] + x15 + x19 + x2*coord[27] + x24*coord[9] + x25*coord[21];
    double x33 = x10*x32;
    double x34 = x1*coord[5] + x2*coord[17] + x26 + x27 + x3*coord[26] + x30 + x5*coord[14];
    double x35 = x0*coord[25] + x2*coord[28] + x21 + x24*coord[10] + x25*coord[22] + x6 + x8;
    double x36 = x17*x35;
    double x37 = x20*x35;
    double x38 = x22*x32;
    double x39 = -x18*x29 + x20*x22*x29 + x31*x33 - x31*x37 + x34*x36 - x34*x38;
    double x40 = 1.0/x39;
    double x41 = -x2*U[19];
    double x42 = x7*U[1];
    double x43 = -x0*U[13] + x42;
    double x44 = x40*(x0*U[25] + x2*U[28] + x24*U[10] + x25*U[22] + x41 + x43);
    double x45 = x33 - x37;
    double x46 = -x3*U[22];
    double x47 = x40*(x0*U[16] + x11*U[7] + x14*U[19] + x3*U[28] + x43 + x46);
    double x48 = x36 - x38;
    double x49 = x40*(x1*U[4] + x2*U[16] + x3*U[25] + x41 + x42 + x46 + x5*U[13]);
    double x50 = x23*x44 + x45*x47 + x48*x49;
    double x51 = x10*x31 - x22*x34;
    double x52 = -x2*U[20];
    double x53 = x7*U[2];
    double x54 = -x0*U[14] + x53;
    double x55 = x40*(x0*U[26] + x2*U[29] + x24*U[11] + x25*U[23] + x52 + x54);
    double x56 = -x10*x29 + x34*x35;
    double x57 = -x3*U[23];
    double x58 = x40*(x0*U[17] + x11*U[8] + x14*U[20] + x3*U[29] + x54 + x57);
    double x59 = x22*x29 - x31*x35;
    double x60 = x40*(x1*U[5] + x2*U[17] + x3*U[26] + x5*U[14] + x52 + x53 + x57);
    double x61 = x51*x55 + x56*x58 + x59*x60;
    double x62 = x44*x51 + x47*x56 + x49*x59;
    double x63 = x23*x55 + x45*x58 + x48*x60 + 1.0;
    double x64 = x50*x61 - x62*x63;
    double x65 = -x2*U[18];
    double x66 = x7*U[0];
    double x67 = -x0*U[12] + x66;
    double x68 = x40*(x0*U[24] + x2*U[27] + x24*U[9] + x25*U[21] + x65 + x67);
    double x69 = -x3*U[21];
    double x70 = x40*(x0*U[15] + x11*U[6] + x14*U[18] + x3*U[27] + x67 + x69);
    double x71 = x40*(x1*U[3] + x2*U[15] + x3*U[24] + x5*U[12] + x65 + x66 + x69);
    double x72 = x23*x68 + x45*x70 + x48*x71;
    double x73 = x51*x68 + x56*x70 + x59*x71 + 1.0;
    double x74 = -x61*x72 + x63*x73;
    double x75 = -x50*x73 + x62*x72;
    double x76 = S11*x64 + S12*x74 + S13*x75;
    double x77 = S12*x64 + S22*x74 + S23*x75;
    double x78 = S13*x64 + S23*x74 + S33*x75;
    double x79 = -x61*x78 - x62*x77 - x73*x76;
    double x80 = x20*x29 - x32*x34;
    double x81 = x40*x7;
    double x82 = x17*x34 - x20*x31;
    double x83 = -x17*x29 + x31*x32;
    double x84 = x80*x81 + x81*x82 + x81*x83;
    double x85 = x55*x82 + x58*x80 + x60*x83;
    double x86 = x44*x82 + x47*x80 + x49*x83 + 1.0;
    double x87 = -x61*x86 + x62*x85;
    double x88 = x68*x82 + x70*x80 + x71*x83;
    double x89 = x61*x88 - x73*x85;
    double x90 = -x62*x88 + x73*x86;
    double x91 = S11*x87 + S12*x89 + S13*x90;
    double x92 = S12*x87 + S22*x89 + S23*x90;
    double x93 = S13*x87 + S23*x89 + S33*x90;
    double x94 = -x61*x93 - x62*x92 - x73*x91;
    double x95 = x23*x81 + x45*x81 + x48*x81;
    double x96 = x51*x81 + x56*x81 + x59*x81;
    double x97 = -x50*x85 + x63*x86;
    double x98 = -x63*x88 + x72*x85;
    double x99 = x50*x88 - x72*x86;
    double x100 = S11*x97 + S12*x98 + S13*x99;
    double x101 = S12*x97 + S22*x98 + S23*x99;
    double x102 = S13*x97 + S23*x98 + S33*x99;
    double x103 = -1.0*PENER - 1.0*SENER;
    double x104 = -x100*x73 - x101*x62 - x102*x61 - x103;
    double x105 = -x85*x93 - x86*x92 - x88*x91;
    double x106 = -x100*x88 - x101*x86 - x102*x85;
    double x107 = -x103 - x76*x88 - x77*x86 - x78*x85;
    double x108 = -x50*x77 - x63*x78 - x72*x76;
    double x109 = -x100*x72 - x101*x50 - x102*x63;
    double x110 = -x103 - x50*x92 - x63*x93 - x72*x91;
    double x111 = x1*x40;
    double x112 = x111*x48;
    double x113 = x111*x59;
    double x114 = x111*x83;
    double x115 = x11*x40;
    double x116 = x115*x80;
    double x117 = x115*x45;
    double x118 = x115*x56;
    double x119 = x24*x40;
    double x120 = x119*x23;
    double x121 = x119*x51;
    double x122 = x119*x82;
    double x123 = x0*x40;
    double x124 = x123*x80;
    double x125 = x123*x82;
    double x126 = -x124 - x125 + x40*x5*x83;
    double x127 = x123*x23;
    double x128 = x123*x45;
    double x129 = -x127 - x128 + x40*x48*x5;
    double x130 = x123*x51;
    double x131 = x123*x56;
    double x132 = -x130 - x131 + x40*x5*x59;
    double x133 = x2*x40;
    double x134 = x133*x83;
    double x135 = x124 + x134;
    double x136 = x133*x48;
    double x137 = x128 + x136;
    double x138 = x133*x59;
    double x139 = x131 + x138;
    double x140 = x133*x82;
    double x141 = -x134 + x14*x40*x80 - x140;
    double x142 = x133*x23;
    double x143 = -x136 + x14*x40*x45 - x142;
    double x144 = x133*x51;
    double x145 = -x138 + x14*x40*x56 - x144;
    double x146 = x3*x40;
    double x147 = x146*x80;
    double x148 = x146*x83;
    double x149 = -x147 - x148 + x25*x40*x82;
    double x150 = x146*x45;
    double x151 = x146*x48;
    double x152 = -x150 - x151 + x23*x25*x40;
    double x153 = x146*x56;
    double x154 = x146*x59;
    double x155 = -x153 - x154 + x25*x40*x51;
    double x156 = x127 + x151;
    double x157 = x130 + x154;
    double x158 = x125 + x148;
    double x159 = x142 + x150;
    double x160 = x144 + x153;
    double x161 = x140 + x147;
    
    res_0[0] = x39*(x104*x96 + x79*x84 + x94*x95);
    res_0[1] = x39*(x105*x95 + x106*x96 + x107*x84);
    res_0[2] = x39*(x108*x84 + x109*x96 + x110*x95);
    res_0[3] = x39*(x104*x113 + x112*x94 + x114*x79);
    res_0[4] = x39*(x105*x112 + x106*x113 + x107*x114);
    res_0[5] = x39*(x108*x114 + x109*x113 + x110*x112);
    res_0[6] = x39*(x104*x118 + x116*x79 + x117*x94);
    res_0[7] = x39*(x105*x117 + x106*x118 + x107*x116);
    res_0[8] = x39*(x108*x116 + x109*x118 + x110*x117);
    res_0[9] = x39*(x104*x121 + x120*x94 + x122*x79);
    res_0[10] = x39*(x105*x120 + x106*x121 + x107*x122);
    res_0[11] = x39*(x108*x122 + x109*x121 + x110*x120);
    res_0[12] = x39*(x104*x132 + x126*x79 + x129*x94);
    res_0[13] = x39*(x105*x129 + x106*x132 + x107*x126);
    res_0[14] = x39*(x108*x126 + x109*x132 + x110*x129);
    res_0[15] = x39*(x104*x139 + x135*x79 + x137*x94);
    res_0[16] = x39*(x105*x137 + x106*x139 + x107*x135);
    res_0[17] = x39*(x108*x135 + x109*x139 + x110*x137);
    res_0[18] = x39*(x104*x145 + x141*x79 + x143*x94);
    res_0[19] = x39*(x105*x143 + x106*x145 + x107*x141);
    res_0[20] = x39*(x108*x141 + x109*x145 + x110*x143);
    res_0[21] = x39*(x104*x155 + x149*x79 + x152*x94);
    res_0[22] = x39*(x105*x152 + x106*x155 + x107*x149);
    res_0[23] = x39*(x108*x149 + x109*x155 + x110*x152);
    res_0[24] = x39*(x104*x157 + x156*x94 + x158*x79);
    res_0[25] = x39*(x105*x156 + x106*x157 + x107*x158);
    res_0[26] = x39*(x108*x158 + x109*x157 + x110*x156);
    res_0[27] = x39*(x104*x160 + x159*x94 + x161*x79);
    res_0[28] = x39*(x105*x159 + x106*x160 + x107*x161);
    res_0[29] = x39*(x108*x161 + x109*x160 + x110*x159);
}

Conf_Forces_API void Integration_C3D10_static_mbf(size_t num_elem,double Coords[][10][3],double Element_U[][10][3],double S[][4][6],
    double PENER[][4],double SENER[][4],double Conf_Force[][10][3]){
    
    //Inputs:
    //    Coords.........Element_Nodal_Coordinates [num_elem][numNodes,3]
    //    Element_U......Element_Nodal_Displacements [num_elem][numNodes,3]
    //    S..............Stress at integration points [num_elem][num_int_points,6]
    //    PENER..........Plastic strain energy [num_elem][num_int_points]
    //    SENER..........Elastic strain energy [num_elem][num_int_points]
    //    
    double weights[4]={0.041666666666666664,0.041666666666666664,0.041666666666666664,0.041666666666666664,};
    double int_points[4][3]={
        {0.1381966011250105,0.1381966011250105,0.1381966011250105,},
        {0.5854101966249685,0.1381966011250105,0.1381966011250105,},
        {0.1381966011250105,0.5854101966249685,0.1381966011250105,},
        {0.1381966011250105,0.1381966011250105,0.5854101966249685,}};

    #pragma omp parallel for
    for (size_t i=0;i<num_elem;i++){
        double TMP[10][3];
        for (size_t j=0;j<4;j++){
            C3D10_static_mbf((double*) int_points[j],(double*) Coords[i],(double*) Element_U[i],(double*) S[i][j],PENER[i][j],SENER[i][j],(double*) TMP);
            for (size_t num_node=0;num_node<10;num_node++){
                for (size_t ii=0;ii<3;ii++){
                    Conf_Force[i][num_node][ii]+=TMP[num_node][ii]*weights[j];
                }
            }
        }
    }
}

void C3D10_static_dbf(double *rst,double *coord,double *Element_U,double *S,double PENER,double SENER,double *res_0){
    double r = rst[0];
    double s = rst[1];
    double t = rst[2];
    double *U = Element_U;
    double S11 = S[0];
    double S22 = S[1];
    double S33 = S[2];
    double S12 = S[3];
    double S13 = S[4];
    double S23 = S[5];
    
    double x0 = 4*r;
    double x1 = x0 - 1;
    double x2 = 4*s;
    double x3 = 4*t;
    double x4 = x2 + x3;
    double x5 = -8*r - x4 + 4;
    double x6 = -x2*coord[19];
    double x7 = x0 + x4 - 3;
    double x8 = x7*coord[1];
    double x9 = -x3*coord[22] + x8;
    double x10 = x1*coord[4] + x2*coord[16] + x3*coord[25] + x5*coord[13] + x6 + x9;
    double x11 = x2 - 1;
    double x12 = x7*coord[0];
    double x13 = x0 - 4;
    double x14 = -8*s - x13 - x3;
    double x15 = -x0*coord[12];
    double x16 = -x3*coord[21];
    double x17 = x0*coord[15] + x11*coord[6] + x12 + x14*coord[18] + x15 + x16 + x3*coord[27];
    double x18 = x10*x17;
    double x19 = x12 - x2*coord[18];
    double x20 = x1*coord[3] + x16 + x19 + x2*coord[15] + x3*coord[24] + x5*coord[12];
    double x21 = -x0*coord[13];
    double x22 = x0*coord[16] + x11*coord[7] + x14*coord[19] + x21 + x3*coord[28] + x9;
    double x23 = -x18 + x20*x22;
    double x24 = x3 - 1;
    double x25 = -8*t - x13 - x2;
    double x26 = -x2*coord[20];
    double x27 = x7*coord[2];
    double x28 = -x0*coord[14] + x27;
    double x29 = x0*coord[26] + x2*coord[29] + x24*coord[11] + x25*coord[23] + x26 + x28;
    double x30 = -x3*coord[23];
    double x31 = x0*coord[17] + x11*coord[8] + x14*coord[20] + x28 + x3*coord[29] + x30;
    double x32 = x0*coord[24] + x15 + x19 + x2*coord[27] + x24*coord[9] + x25*coord[21];
    double x33 = x10*x32;
    double x34 = x1*coord[5] + x2*coord[17] + x26 + x27 + x3*coord[26] + x30 + x5*coord[14];
    double x35 = x0*coord[25] + x2*coord[28] + x21 + x24*coord[10] + x25*coord[22] + x6 + x8;
    double x36 = x17*x35;
    double x37 = x20*x35;
    double x38 = x22*x32;
    double x39 = -x18*x29 + x20*x22*x29 + x31*x33 - x31*x37 + x34*x36 - x34*x38;
    double x40 = 1.0/x39;
    double x41 = -x2*U[19];
    double x42 = x7*U[1];
    double x43 = -x0*U[13] + x42;
    double x44 = x40*(x0*U[25] + x2*U[28] + x24*U[10] + x25*U[22] + x41 + x43);
    double x45 = x33 - x37;
    double x46 = -x3*U[22];
    double x47 = x40*(x0*U[16] + x11*U[7] + x14*U[19] + x3*U[28] + x43 + x46);
    double x48 = x36 - x38;
    double x49 = x40*(x1*U[4] + x2*U[16] + x3*U[25] + x41 + x42 + x46 + x5*U[13]);
    double x50 = x23*x44 + x45*x47 + x48*x49;
    double x51 = x10*x31 - x22*x34;
    double x52 = -x2*U[20];
    double x53 = x7*U[2];
    double x54 = -x0*U[14] + x53;
    double x55 = x40*(x0*U[26] + x2*U[29] + x24*U[11] + x25*U[23] + x52 + x54);
    double x56 = -x10*x29 + x34*x35;
    double x57 = -x3*U[23];
    double x58 = x40*(x0*U[17] + x11*U[8] + x14*U[20] + x3*U[29] + x54 + x57);
    double x59 = x22*x29 - x31*x35;
    double x60 = x40*(x1*U[5] + x2*U[17] + x3*U[26] + x5*U[14] + x52 + x53 + x57);
    double x61 = x51*x55 + x56*x58 + x59*x60;
    double x62 = x44*x51 + x47*x56 + x49*x59;
    double x63 = x23*x55 + x45*x58 + x48*x60;
    double x64 = x63 + 1.0;
    double x65 = x50*x61 - x62*x64;
    double x66 = -x2*U[18];
    double x67 = x7*U[0];
    double x68 = -x0*U[12] + x67;
    double x69 = x40*(x0*U[24] + x2*U[27] + x24*U[9] + x25*U[21] + x66 + x68);
    double x70 = -x3*U[21];
    double x71 = x40*(x0*U[15] + x11*U[6] + x14*U[18] + x3*U[27] + x68 + x70);
    double x72 = x40*(x1*U[3] + x2*U[15] + x3*U[24] + x5*U[12] + x66 + x67 + x70);
    double x73 = x23*x69 + x45*x71 + x48*x72;
    double x74 = x51*x69 + x56*x71 + x59*x72;
    double x75 = x74 + 1.0;
    double x76 = -x61*x73 + x64*x75;
    double x77 = -x50*x75 + x62*x73;
    double x78 = S11*x65 + S12*x76 + S13*x77;
    double x79 = S12*x65 + S22*x76 + S23*x77;
    double x80 = S13*x65 + S23*x76 + S33*x77;
    double x81 = -x61*x80 - x62*x79 - x74*x78;
    double x82 = x20*x29 - x32*x34;
    double x83 = x40*x7;
    double x84 = x17*x34 - x20*x31;
    double x85 = -x17*x29 + x31*x32;
    double x86 = x82*x83 + x83*x84 + x83*x85;
    double x87 = x55*x84 + x58*x82 + x60*x85;
    double x88 = x44*x84 + x47*x82 + x49*x85;
    double x89 = x88 + 1.0;
    double x90 = -x61*x89 + x62*x87;
    double x91 = x69*x84 + x71*x82 + x72*x85;
    double x92 = x61*x91 - x75*x87;
    double x93 = -x62*x91 + x75*x89;
    double x94 = S11*x90 + S12*x92 + S13*x93;
    double x95 = S12*x90 + S22*x92 + S23*x93;
    double x96 = S13*x90 + S23*x92 + S33*x93;
    double x97 = -x61*x96 - x62*x95 - x74*x94;
    double x98 = x23*x83 + x45*x83 + x48*x83;
    double x99 = x51*x83 + x56*x83 + x59*x83;
    double x100 = -x50*x87 + x64*x89;
    double x101 = -x64*x91 + x73*x87;
    double x102 = x50*x91 - x73*x89;
    double x103 = S11*x100 + S12*x101 + S13*x102;
    double x104 = S12*x100 + S22*x101 + S23*x102;
    double x105 = S13*x100 + S23*x101 + S33*x102;
    double x106 = -1.0*PENER - 1.0*SENER;
    double x107 = -x103*x74 - x104*x62 - x105*x61 - x106;
    double x108 = -x87*x96 - x88*x95 - x91*x94;
    double x109 = -x103*x91 - x104*x88 - x105*x87;
    double x110 = -x106 - x78*x91 - x79*x88 - x80*x87;
    double x111 = -x50*x79 - x63*x80 - x73*x78;
    double x112 = -x103*x73 - x104*x50 - x105*x63;
    double x113 = -x106 - x50*x95 - x63*x96 - x73*x94;
    double x114 = x1*x40;
    double x115 = x114*x48;
    double x116 = x114*x59;
    double x117 = x114*x85;
    double x118 = x11*x40;
    double x119 = x118*x82;
    double x120 = x118*x45;
    double x121 = x118*x56;
    double x122 = x24*x40;
    double x123 = x122*x23;
    double x124 = x122*x51;
    double x125 = x122*x84;
    double x126 = x0*x40;
    double x127 = x126*x82;
    double x128 = x126*x84;
    double x129 = -x127 - x128 + x40*x5*x85;
    double x130 = x126*x23;
    double x131 = x126*x45;
    double x132 = -x130 - x131 + x40*x48*x5;
    double x133 = x126*x51;
    double x134 = x126*x56;
    double x135 = -x133 - x134 + x40*x5*x59;
    double x136 = x2*x40;
    double x137 = x136*x85;
    double x138 = x127 + x137;
    double x139 = x136*x48;
    double x140 = x131 + x139;
    double x141 = x136*x59;
    double x142 = x134 + x141;
    double x143 = x136*x84;
    double x144 = -x137 + x14*x40*x82 - x143;
    double x145 = x136*x23;
    double x146 = -x139 + x14*x40*x45 - x145;
    double x147 = x136*x51;
    double x148 = x14*x40*x56 - x141 - x147;
    double x149 = x3*x40;
    double x150 = x149*x82;
    double x151 = x149*x85;
    double x152 = -x150 - x151 + x25*x40*x84;
    double x153 = x149*x45;
    double x154 = x149*x48;
    double x155 = -x153 - x154 + x23*x25*x40;
    double x156 = x149*x56;
    double x157 = x149*x59;
    double x158 = -x156 - x157 + x25*x40*x51;
    double x159 = x130 + x154;
    double x160 = x133 + x157;
    double x161 = x128 + x151;
    double x162 = x145 + x153;
    double x163 = x147 + x156;
    double x164 = x143 + x150;
    
    res_0[0] = x39*(x107*x99 + x81*x86 + x97*x98);
    res_0[1] = x39*(x108*x98 + x109*x99 + x110*x86);
    res_0[2] = x39*(x111*x86 + x112*x99 + x113*x98);
    res_0[3] = x39*(x107*x116 + x115*x97 + x117*x81);
    res_0[4] = x39*(x108*x115 + x109*x116 + x110*x117);
    res_0[5] = x39*(x111*x117 + x112*x116 + x113*x115);
    res_0[6] = x39*(x107*x121 + x119*x81 + x120*x97);
    res_0[7] = x39*(x108*x120 + x109*x121 + x110*x119);
    res_0[8] = x39*(x111*x119 + x112*x121 + x113*x120);
    res_0[9] = x39*(x107*x124 + x123*x97 + x125*x81);
    res_0[10] = x39*(x108*x123 + x109*x124 + x110*x125);
    res_0[11] = x39*(x111*x125 + x112*x124 + x113*x123);
    res_0[12] = x39*(x107*x135 + x129*x81 + x132*x97);
    res_0[13] = x39*(x108*x132 + x109*x135 + x110*x129);
    res_0[14] = x39*(x111*x129 + x112*x135 + x113*x132);
    res_0[15] = x39*(x107*x142 + x138*x81 + x140*x97);
    res_0[16] = x39*(x108*x140 + x109*x142 + x110*x138);
    res_0[17] = x39*(x111*x138 + x112*x142 + x113*x140);
    res_0[18] = x39*(x107*x148 + x144*x81 + x146*x97);
    res_0[19] = x39*(x108*x146 + x109*x148 + x110*x144);
    res_0[20] = x39*(x111*x144 + x112*x148 + x113*x146);
    res_0[21] = x39*(x107*x158 + x152*x81 + x155*x97);
    res_0[22] = x39*(x108*x155 + x109*x158 + x110*x152);
    res_0[23] = x39*(x111*x152 + x112*x158 + x113*x155);
    res_0[24] = x39*(x107*x160 + x159*x97 + x161*x81);
    res_0[25] = x39*(x108*x159 + x109*x160 + x110*x161);
    res_0[26] = x39*(x111*x161 + x112*x160 + x113*x159);
    res_0[27] = x39*(x107*x163 + x162*x97 + x164*x81);
    res_0[28] = x39*(x108*x162 + x109*x163 + x110*x164);
    res_0[29] = x39*(x111*x164 + x112*x163 + x113*x162);
}

Conf_Forces_API void Integration_C3D10_static_dbf(size_t num_elem,double Coords[][10][3],double Element_U[][10][3],double S[][4][6],
    double PENER[][4],double SENER[][4],double Conf_Force[][10][3]){
    
    //Inputs:
    //    Coords.........Element_Nodal_Coordinates [num_elem][numNodes,3]
    //    Element_U......Element_Nodal_Displacements [num_elem][numNodes,3]
    //    S..............Stress at integration points [num_elem][num_int_points,6]
    //    PENER..........Plastic strain energy [num_elem][num_int_points]
    //    SENER..........Elastic strain energy [num_elem][num_int_points]
    //    
    double weights[4]={0.041666666666666664,0.041666666666666664,0.041666666666666664,0.041666666666666664,};
    double int_points[4][3]={
        {0.1381966011250105,0.1381966011250105,0.1381966011250105,},
        {0.5854101966249685,0.1381966011250105,0.1381966011250105,},
        {0.1381966011250105,0.5854101966249685,0.1381966011250105,},
        {0.1381966011250105,0.1381966011250105,0.5854101966249685,}};

    #pragma omp parallel for
    for (size_t i=0;i<num_elem;i++){
        double TMP[10][3];
        for (size_t j=0;j<4;j++){
            C3D10_static_dbf((double*) int_points[j],(double*) Coords[i],(double*) Element_U[i],(double*) S[i][j],PENER[i][j],SENER[i][j],(double*) TMP);
            for (size_t num_node=0;num_node<10;num_node++){
                for (size_t ii=0;ii<3;ii++){
                    Conf_Force[i][num_node][ii]+=TMP[num_node][ii]*weights[j];
                }
            }
        }
    }
}

void C3D10R_dynamic(double *rst,double *coord,double rho,double *Element_U,double *Element_V,double *Element_A,double *S,double PENER,double SENER,double *res_0){
    double r = rst[0];
    double s = rst[1];
    double t = rst[2];
    double *U = Element_U;
    double *V = Element_V;
    double *A = Element_A;
    double S11 = S[0];
    double S22 = S[1];
    double S33 = S[2];
    double S12 = S[3];
    double S13 = S[4];
    double S23 = S[5];
    
    double x0 = 4*r;
    double x1 = x0 - 1;
    double x2 = 4*s;
    double x3 = 4*t;
    double x4 = x2 + x3;
    double x5 = -8*r - x4 + 4;
    double x6 = -x2*coord[19];
    double x7 = x0 + x4 - 3;
    double x8 = x7*coord[1];
    double x9 = -x3*coord[22] + x8;
    double x10 = x1*coord[4] + x2*coord[16] + x3*coord[25] + x5*coord[13] + x6 + x9;
    double x11 = x2 - 1;
    double x12 = x7*coord[0];
    double x13 = x0 - 4;
    double x14 = -8*s - x13 - x3;
    double x15 = -x0*coord[12];
    double x16 = -x3*coord[21];
    double x17 = x0*coord[15] + x11*coord[6] + x12 + x14*coord[18] + x15 + x16 + x3*coord[27];
    double x18 = x10*x17;
    double x19 = x12 - x2*coord[18];
    double x20 = x1*coord[3] + x16 + x19 + x2*coord[15] + x3*coord[24] + x5*coord[12];
    double x21 = -x0*coord[13];
    double x22 = x0*coord[16] + x11*coord[7] + x14*coord[19] + x21 + x3*coord[28] + x9;
    double x23 = -x18 + x20*x22;
    double x24 = x3 - 1;
    double x25 = -8*t - x13 - x2;
    double x26 = -x2*coord[20];
    double x27 = x7*coord[2];
    double x28 = -x0*coord[14] + x27;
    double x29 = x0*coord[26] + x2*coord[29] + x24*coord[11] + x25*coord[23] + x26 + x28;
    double x30 = -x3*coord[23];
    double x31 = x0*coord[17] + x11*coord[8] + x14*coord[20] + x28 + x3*coord[29] + x30;
    double x32 = x0*coord[24] + x15 + x19 + x2*coord[27] + x24*coord[9] + x25*coord[21];
    double x33 = x10*x32;
    double x34 = x1*coord[5] + x2*coord[17] + x26 + x27 + x3*coord[26] + x30 + x5*coord[14];
    double x35 = x0*coord[25] + x2*coord[28] + x21 + x24*coord[10] + x25*coord[22] + x6 + x8;
    double x36 = x17*x35;
    double x37 = x20*x35;
    double x38 = x22*x32;
    double x39 = -x18*x29 + x20*x22*x29 + x31*x33 - x31*x37 + x34*x36 - x34*x38;
    double x40 = 1.0/x39;
    double x41 = -x2*U[19];
    double x42 = x7*U[1];
    double x43 = -x0*U[13] + x42;
    double x44 = x40*(x0*U[25] + x2*U[28] + x24*U[10] + x25*U[22] + x41 + x43);
    double x45 = x33 - x37;
    double x46 = -x3*U[22];
    double x47 = x40*(x0*U[16] + x11*U[7] + x14*U[19] + x3*U[28] + x43 + x46);
    double x48 = x36 - x38;
    double x49 = x40*(x1*U[4] + x2*U[16] + x3*U[25] + x41 + x42 + x46 + x5*U[13]);
    double x50 = x23*x44 + x45*x47 + x48*x49;
    double x51 = x10*x31 - x22*x34;
    double x52 = -x2*U[20];
    double x53 = x7*U[2];
    double x54 = -x0*U[14] + x53;
    double x55 = x40*(x0*U[26] + x2*U[29] + x24*U[11] + x25*U[23] + x52 + x54);
    double x56 = -x10*x29 + x34*x35;
    double x57 = -x3*U[23];
    double x58 = x40*(x0*U[17] + x11*U[8] + x14*U[20] + x3*U[29] + x54 + x57);
    double x59 = x22*x29 - x31*x35;
    double x60 = x40*(x1*U[5] + x2*U[17] + x3*U[26] + x5*U[14] + x52 + x53 + x57);
    double x61 = x51*x55 + x56*x58 + x59*x60;
    double x62 = x44*x51 + x47*x56 + x49*x59;
    double x63 = x23*x55 + x45*x58 + x48*x60;
    double x64 = x63 + 1.0;
    double x65 = x50*x61 - x62*x64;
    double x66 = -x2*U[18];
    double x67 = x7*U[0];
    double x68 = -x0*U[12] + x67;
    double x69 = x40*(x0*U[24] + x2*U[27] + x24*U[9] + x25*U[21] + x66 + x68);
    double x70 = -x3*U[21];
    double x71 = x40*(x0*U[15] + x11*U[6] + x14*U[18] + x3*U[27] + x68 + x70);
    double x72 = x40*(x1*U[3] + x2*U[15] + x3*U[24] + x5*U[12] + x66 + x67 + x70);
    double x73 = x23*x69 + x45*x71 + x48*x72;
    double x74 = x51*x69 + x56*x71 + x59*x72;
    double x75 = x74 + 1.0;
    double x76 = -x61*x73 + x64*x75;
    double x77 = -x50*x75 + x62*x73;
    double x78 = S11*x65 + S12*x76 + S13*x77;
    double x79 = S12*x65 + S22*x76 + S23*x77;
    double x80 = S13*x65 + S23*x76 + S33*x77;
    double x81 = -x61*x80 - x62*x79 - x74*x78;
    double x82 = x20*x29 - x32*x34;
    double x83 = x40*x7;
    double x84 = x17*x34 - x20*x31;
    double x85 = -x17*x29 + x31*x32;
    double x86 = x82*x83 + x83*x84 + x83*x85;
    double x87 = x55*x84 + x58*x82 + x60*x85;
    double x88 = x44*x84 + x47*x82 + x49*x85;
    double x89 = x88 + 1.0;
    double x90 = -x61*x89 + x62*x87;
    double x91 = x69*x84 + x71*x82 + x72*x85;
    double x92 = x61*x91 - x75*x87;
    double x93 = -x62*x91 + x75*x89;
    double x94 = S11*x90 + S12*x92 + S13*x93;
    double x95 = S12*x90 + S22*x92 + S23*x93;
    double x96 = S13*x90 + S23*x92 + S33*x93;
    double x97 = -x61*x96 - x62*x95 - x74*x94;
    double x98 = x23*x83 + x45*x83 + x48*x83;
    double x99 = x51*x83 + x56*x83 + x59*x83;
    double x100 = -x50*x87 + x64*x89;
    double x101 = -x64*x91 + x73*x87;
    double x102 = x50*x91 - x73*x89;
    double x103 = S11*x100 + S12*x101 + S13*x102;
    double x104 = S12*x100 + S22*x101 + S23*x102;
    double x105 = S13*x100 + S23*x101 + S33*x102;
    double x106 = r*r;
    double x107 = 2*x106;
    double x108 = -r + x107;
    double x109 = s*s;
    double x110 = 2*x109;
    double x111 = -s + x110;
    double x112 = t*t;
    double x113 = 2*x112;
    double x114 = -t + x113;
    double x115 = s*x0;
    double x116 = t*x0;
    double x117 = x115 + x116;
    double x118 = x0 - 4*x106 - x117;
    double x119 = t*x2;
    double x120 = -4*x109 - x115 - x119 + x2;
    double x121 = -4*x112 - x116 - x119 + x3;
    double x122 = -3*r - 3*s - 3*t + x107 + x110 + x113 + x117 + x119 + 1;
    double x123 = x108*V[3] + x111*V[6] + x114*V[9] + x115*V[15] + x116*V[24] + x118*V[12] + x119*V[27] + x120*V[18] + x121*V[21] + x122*V[0];
    double x124 = x108*V[4] + x111*V[7] + x114*V[10] + x115*V[16] + x116*V[25] + x118*V[13] + x119*V[28] + x120*V[19] + x121*V[22] + x122*V[1];
    double x125 = x108*V[5] + x111*V[8] + x114*V[11] + x115*V[17] + x116*V[26] + x118*V[14] + x119*V[29] + x120*V[20] + x121*V[23] + x122*V[2];
    double x126 = -1.0*PENER - 1.0*SENER + 0.5*rho*(x123*x123 + x124*x124 + x125*x125);
    double x127 = -x103*x74 - x104*x62 - x105*x61 - x126;
    double x128 = x108*A[3] + x111*A[6] + x114*A[9] + x115*A[15] + x116*A[24] + x118*A[12] + x119*A[27] + x120*A[18] + x121*A[21] + x122*A[0];
    double x129 = x108*A[4] + x111*A[7] + x114*A[10] + x115*A[16] + x116*A[25] + x118*A[13] + x119*A[28] + x120*A[19] + x121*A[22] + x122*A[1];
    double x130 = x108*A[5] + x111*A[8] + x114*A[11] + x115*A[17] + x116*A[26] + x118*A[14] + x119*A[29] + x120*A[20] + x121*A[23] + x122*A[2];
    double x131 = -x2*V[18];
    double x132 = x7*V[0];
    double x133 = -x0*V[12] + x132;
    double x134 = x0*V[24] + x131 + x133 + x2*V[27] + x24*V[9] + x25*V[21];
    double x135 = x40*x51;
    double x136 = -x3*V[21];
    double x137 = x0*V[15] + x11*V[6] + x133 + x136 + x14*V[18] + x3*V[27];
    double x138 = x40*x56;
    double x139 = x1*V[3] + x131 + x132 + x136 + x2*V[15] + x3*V[24] + x5*V[12];
    double x140 = x40*x59;
    double x141 = -x2*V[19];
    double x142 = x7*V[1];
    double x143 = -x0*V[13] + x142;
    double x144 = x0*V[25] + x141 + x143 + x2*V[28] + x24*V[10] + x25*V[22];
    double x145 = -x3*V[22];
    double x146 = x0*V[16] + x11*V[7] + x14*V[19] + x143 + x145 + x3*V[28];
    double x147 = x1*V[4] + x141 + x142 + x145 + x2*V[16] + x3*V[25] + x5*V[13];
    double x148 = -x2*V[20];
    double x149 = x7*V[2];
    double x150 = -x0*V[14] + x149;
    double x151 = x0*V[26] + x148 + x150 + x2*V[29] + x24*V[11] + x25*V[23];
    double x152 = -x3*V[23];
    double x153 = x0*V[17] + x11*V[8] + x14*V[20] + x150 + x152 + x3*V[29];
    double x154 = x1*V[5] + x148 + x149 + x152 + x2*V[17] + x3*V[26] + x5*V[14];
    double x155 = x123*(x134*x135 + x137*x138 + x139*x140) + x124*(x135*x144 + x138*x146 + x140*x147) + x125*(x135*x151 + x138*x153 + x140*x154) + x128*x74 + x129*x62 + x130*x61;
    double x156 = rho*x122;
    double x157 = -x87*x96 - x88*x95 - x91*x94;
    double x158 = -x103*x91 - x104*x88 - x105*x87;
    double x159 = -x126 - x78*x91 - x79*x88 - x80*x87;
    double x160 = x40*x82;
    double x161 = x40*x84;
    double x162 = x40*x85;
    double x163 = x123*(x134*x161 + x137*x160 + x139*x162) + x124*(x144*x161 + x146*x160 + x147*x162) + x125*(x151*x161 + x153*x160 + x154*x162) + x128*x91 + x129*x88 + x130*x87;
    double x164 = -x50*x79 - x63*x80 - x73*x78;
    double x165 = -x103*x73 - x104*x50 - x105*x63;
    double x166 = -x126 - x50*x95 - x63*x96 - x73*x94;
    double x167 = x23*x40;
    double x168 = x40*x45;
    double x169 = x40*x48;
    double x170 = x123*(x134*x167 + x137*x168 + x139*x169) + x124*(x144*x167 + x146*x168 + x147*x169) + x125*(x151*x167 + x153*x168 + x154*x169) + x128*x73 + x129*x50 + x130*x63;
    double x171 = rho*x108;
    double x172 = x1*x169;
    double x173 = x1*x140;
    double x174 = x1*x162;
    double x175 = rho*x111;
    double x176 = x11*x160;
    double x177 = x11*x168;
    double x178 = x11*x138;
    double x179 = rho*x114;
    double x180 = x167*x24;
    double x181 = x135*x24;
    double x182 = x161*x24;
    double x183 = x0*x160;
    double x184 = x0*x161;
    double x185 = -x183 - x184 + x40*x5*x85;
    double x186 = x0*x167;
    double x187 = x0*x168;
    double x188 = -x186 - x187 + x40*x48*x5;
    double x189 = x0*x135;
    double x190 = x0*x138;
    double x191 = -x189 - x190 + x40*x5*x59;
    double x192 = rho*x118;
    double x193 = x162*x2;
    double x194 = x183 + x193;
    double x195 = x169*x2;
    double x196 = x187 + x195;
    double x197 = x140*x2;
    double x198 = x190 + x197;
    double x199 = rho*x115;
    double x200 = x161*x2;
    double x201 = x14*x40*x82 - x193 - x200;
    double x202 = x167*x2;
    double x203 = x14*x40*x45 - x195 - x202;
    double x204 = x135*x2;
    double x205 = x14*x40*x56 - x197 - x204;
    double x206 = rho*x120;
    double x207 = x160*x3;
    double x208 = x162*x3;
    double x209 = -x207 - x208 + x25*x40*x84;
    double x210 = x168*x3;
    double x211 = x169*x3;
    double x212 = -x210 - x211 + x23*x25*x40;
    double x213 = x138*x3;
    double x214 = x140*x3;
    double x215 = -x213 - x214 + x25*x40*x51;
    double x216 = rho*x121;
    double x217 = x186 + x211;
    double x218 = x189 + x214;
    double x219 = x184 + x208;
    double x220 = rho*x116;
    double x221 = x202 + x210;
    double x222 = x204 + x213;
    double x223 = x200 + x207;
    double x224 = rho*x119;
    
    res_0[0] = x39*(x127*x99 - x155*x156 + x81*x86 + x97*x98);
    res_0[1] = x39*(-x156*x163 + x157*x98 + x158*x99 + x159*x86);
    res_0[2] = x39*(-x156*x170 + x164*x86 + x165*x99 + x166*x98);
    res_0[3] = x39*(x127*x173 - x155*x171 + x172*x97 + x174*x81);
    res_0[4] = x39*(x157*x172 + x158*x173 + x159*x174 - x163*x171);
    res_0[5] = x39*(x164*x174 + x165*x173 + x166*x172 - x170*x171);
    res_0[6] = x39*(x127*x178 - x155*x175 + x176*x81 + x177*x97);
    res_0[7] = x39*(x157*x177 + x158*x178 + x159*x176 - x163*x175);
    res_0[8] = x39*(x164*x176 + x165*x178 + x166*x177 - x170*x175);
    res_0[9] = x39*(x127*x181 - x155*x179 + x180*x97 + x182*x81);
    res_0[10] = x39*(x157*x180 + x158*x181 + x159*x182 - x163*x179);
    res_0[11] = x39*(x164*x182 + x165*x181 + x166*x180 - x170*x179);
    res_0[12] = x39*(x127*x191 - x155*x192 + x185*x81 + x188*x97);
    res_0[13] = x39*(x157*x188 + x158*x191 + x159*x185 - x163*x192);
    res_0[14] = x39*(x164*x185 + x165*x191 + x166*x188 - x170*x192);
    res_0[15] = x39*(x127*x198 - x155*x199 + x194*x81 + x196*x97);
    res_0[16] = x39*(x157*x196 + x158*x198 + x159*x194 - x163*x199);
    res_0[17] = x39*(x164*x194 + x165*x198 + x166*x196 - x170*x199);
    res_0[18] = x39*(x127*x205 - x155*x206 + x201*x81 + x203*x97);
    res_0[19] = x39*(x157*x203 + x158*x205 + x159*x201 - x163*x206);
    res_0[20] = x39*(x164*x201 + x165*x205 + x166*x203 - x170*x206);
    res_0[21] = x39*(x127*x215 - x155*x216 + x209*x81 + x212*x97);
    res_0[22] = x39*(x157*x212 + x158*x215 + x159*x209 - x163*x216);
    res_0[23] = x39*(x164*x209 + x165*x215 + x166*x212 - x170*x216);
    res_0[24] = x39*(x127*x218 - x155*x220 + x217*x97 + x219*x81);
    res_0[25] = x39*(x157*x217 + x158*x218 + x159*x219 - x163*x220);
    res_0[26] = x39*(x164*x219 + x165*x218 + x166*x217 - x170*x220);
    res_0[27] = x39*(x127*x222 - x155*x224 + x221*x97 + x223*x81);
    res_0[28] = x39*(x157*x221 + x158*x222 + x159*x223 - x163*x224);
    res_0[29] = x39*(x164*x223 + x165*x222 + x166*x221 - x170*x224);
}

Conf_Forces_API void Integration_C3D10R_dynamic(size_t num_elem,double Coords[][10][3],double *rho,double Element_U[][10][3],double Element_V[][10][3],double Element_A[][10][3],double S[][1][6],
    double PENER[][1],double SENER[][1],double Conf_Force[][10][3]){
    
    //Inputs:
    //    Coords.........Element_Nodal_Coordinates  [num_elem][numNodes,3]
    //    Element_U......Element_Nodal_Displacements[num_elem][numNodes,3]
    //    Element_V......Element_Nodal_Velocity     [num_elem][numNodes,3]
    //    Element_A......Element_Nodal_Acceleration [num_elem][numNodes,3]
    //    S..............Stress at integration points [num_elem][num_int_points,6]
    //    PENER..........Plastic strain energy [num_elem][num_int_points]
    //    SENER..........Elastic strain energy [num_elem][num_int_points]
    //    
    double weights[1]={0.16666666666666666,};
    double int_points[1][3]={
        {0.25,0.25,0.25,}};

    #pragma omp parallel for
    for (size_t i=0;i<num_elem;i++){
        double TMP[10][3];
        for (size_t j=0;j<1;j++){
            C3D10R_dynamic((double*) int_points[j],(double*) Coords[i],rho[i],(double*) Element_U[i],(double*) Element_V[i],(double*) Element_A[i],(double*) S[i][j],PENER[i][j],SENER[i][j],(double*) TMP);
            for (size_t num_node=0;num_node<10;num_node++){
                for (size_t ii=0;ii<3;ii++){
                    Conf_Force[i][num_node][ii]+=TMP[num_node][ii]*weights[j];
                }
            }
        }
    }
}

void C3D10R_static_mbf(double *rst,double *coord,double *Element_U,double *S,double PENER,double SENER,double *res_0){
    double r = rst[0];
    double s = rst[1];
    double t = rst[2];
    double *U = Element_U;
    double S11 = S[0];
    double S22 = S[1];
    double S33 = S[2];
    double S12 = S[3];
    double S13 = S[4];
    double S23 = S[5];
    
    double x0 = 4*r;
    double x1 = x0 - 1;
    double x2 = 4*s;
    double x3 = 4*t;
    double x4 = x2 + x3;
    double x5 = -8*r - x4 + 4;
    double x6 = -x2*coord[19];
    double x7 = x0 + x4 - 3;
    double x8 = x7*coord[1];
    double x9 = -x3*coord[22] + x8;
    double x10 = x1*coord[4] + x2*coord[16] + x3*coord[25] + x5*coord[13] + x6 + x9;
    double x11 = x2 - 1;
    double x12 = x7*coord[0];
    double x13 = x0 - 4;
    double x14 = -8*s - x13 - x3;
    double x15 = -x0*coord[12];
    double x16 = -x3*coord[21];
    double x17 = x0*coord[15] + x11*coord[6] + x12 + x14*coord[18] + x15 + x16 + x3*coord[27];
    double x18 = x10*x17;
    double x19 = x12 - x2*coord[18];
    double x20 = x1*coord[3] + x16 + x19 + x2*coord[15] + x3*coord[24] + x5*coord[12];
    double x21 = -x0*coord[13];
    double x22 = x0*coord[16] + x11*coord[7] + x14*coord[19] + x21 + x3*coord[28] + x9;
    double x23 = -x18 + x20*x22;
    double x24 = x3 - 1;
    double x25 = -8*t - x13 - x2;
    double x26 = -x2*coord[20];
    double x27 = x7*coord[2];
    double x28 = -x0*coord[14] + x27;
    double x29 = x0*coord[26] + x2*coord[29] + x24*coord[11] + x25*coord[23] + x26 + x28;
    double x30 = -x3*coord[23];
    double x31 = x0*coord[17] + x11*coord[8] + x14*coord[20] + x28 + x3*coord[29] + x30;
    double x32 = x0*coord[24] + x15 + x19 + x2*coord[27] + x24*coord[9] + x25*coord[21];
    double x33 = x10*x32;
    double x34 = x1*coord[5] + x2*coord[17] + x26 + x27 + x3*coord[26] + x30 + x5*coord[14];
    double x35 = x0*coord[25] + x2*coord[28] + x21 + x24*coord[10] + x25*coord[22] + x6 + x8;
    double x36 = x17*x35;
    double x37 = x20*x35;
    double x38 = x22*x32;
    double x39 = -x18*x29 + x20*x22*x29 + x31*x33 - x31*x37 + x34*x36 - x34*x38;
    double x40 = 1.0/x39;
    double x41 = -x2*U[19];
    double x42 = x7*U[1];
    double x43 = -x0*U[13] + x42;
    double x44 = x40*(x0*U[25] + x2*U[28] + x24*U[10] + x25*U[22] + x41 + x43);
    double x45 = x33 - x37;
    double x46 = -x3*U[22];
    double x47 = x40*(x0*U[16] + x11*U[7] + x14*U[19] + x3*U[28] + x43 + x46);
    double x48 = x36 - x38;
    double x49 = x40*(x1*U[4] + x2*U[16] + x3*U[25] + x41 + x42 + x46 + x5*U[13]);
    double x50 = x23*x44 + x45*x47 + x48*x49;
    double x51 = x10*x31 - x22*x34;
    double x52 = -x2*U[20];
    double x53 = x7*U[2];
    double x54 = -x0*U[14] + x53;
    double x55 = x40*(x0*U[26] + x2*U[29] + x24*U[11] + x25*U[23] + x52 + x54);
    double x56 = -x10*x29 + x34*x35;
    double x57 = -x3*U[23];
    double x58 = x40*(x0*U[17] + x11*U[8] + x14*U[20] + x3*U[29] + x54 + x57);
    double x59 = x22*x29 - x31*x35;
    double x60 = x40*(x1*U[5] + x2*U[17] + x3*U[26] + x5*U[14] + x52 + x53 + x57);
    double x61 = x51*x55 + x56*x58 + x59*x60;
    double x62 = x44*x51 + x47*x56 + x49*x59;
    double x63 = x23*x55 + x45*x58 + x48*x60 + 1.0;
    double x64 = x50*x61 - x62*x63;
    double x65 = -x2*U[18];
    double x66 = x7*U[0];
    double x67 = -x0*U[12] + x66;
    double x68 = x40*(x0*U[24] + x2*U[27] + x24*U[9] + x25*U[21] + x65 + x67);
    double x69 = -x3*U[21];
    double x70 = x40*(x0*U[15] + x11*U[6] + x14*U[18] + x3*U[27] + x67 + x69);
    double x71 = x40*(x1*U[3] + x2*U[15] + x3*U[24] + x5*U[12] + x65 + x66 + x69);
    double x72 = x23*x68 + x45*x70 + x48*x71;
    double x73 = x51*x68 + x56*x70 + x59*x71 + 1.0;
    double x74 = -x61*x72 + x63*x73;
    double x75 = -x50*x73 + x62*x72;
    double x76 = S11*x64 + S12*x74 + S13*x75;
    double x77 = S12*x64 + S22*x74 + S23*x75;
    double x78 = S13*x64 + S23*x74 + S33*x75;
    double x79 = -x61*x78 - x62*x77 - x73*x76;
    double x80 = x20*x29 - x32*x34;
    double x81 = x40*x7;
    double x82 = x17*x34 - x20*x31;
    double x83 = -x17*x29 + x31*x32;
    double x84 = x80*x81 + x81*x82 + x81*x83;
    double x85 = x55*x82 + x58*x80 + x60*x83;
    double x86 = x44*x82 + x47*x80 + x49*x83 + 1.0;
    double x87 = -x61*x86 + x62*x85;
    double x88 = x68*x82 + x70*x80 + x71*x83;
    double x89 = x61*x88 - x73*x85;
    double x90 = -x62*x88 + x73*x86;
    double x91 = S11*x87 + S12*x89 + S13*x90;
    double x92 = S12*x87 + S22*x89 + S23*x90;
    double x93 = S13*x87 + S23*x89 + S33*x90;
    double x94 = -x61*x93 - x62*x92 - x73*x91;
    double x95 = x23*x81 + x45*x81 + x48*x81;
    double x96 = x51*x81 + x56*x81 + x59*x81;
    double x97 = -x50*x85 + x63*x86;
    double x98 = -x63*x88 + x72*x85;
    double x99 = x50*x88 - x72*x86;
    double x100 = S11*x97 + S12*x98 + S13*x99;
    double x101 = S12*x97 + S22*x98 + S23*x99;
    double x102 = S13*x97 + S23*x98 + S33*x99;
    double x103 = -1.0*PENER - 1.0*SENER;
    double x104 = -x100*x73 - x101*x62 - x102*x61 - x103;
    double x105 = -x85*x93 - x86*x92 - x88*x91;
    double x106 = -x100*x88 - x101*x86 - x102*x85;
    double x107 = -x103 - x76*x88 - x77*x86 - x78*x85;
    double x108 = -x50*x77 - x63*x78 - x72*x76;
    double x109 = -x100*x72 - x101*x50 - x102*x63;
    double x110 = -x103 - x50*x92 - x63*x93 - x72*x91;
    double x111 = x1*x40;
    double x112 = x111*x48;
    double x113 = x111*x59;
    double x114 = x111*x83;
    double x115 = x11*x40;
    double x116 = x115*x80;
    double x117 = x115*x45;
    double x118 = x115*x56;
    double x119 = x24*x40;
    double x120 = x119*x23;
    double x121 = x119*x51;
    double x122 = x119*x82;
    double x123 = x0*x40;
    double x124 = x123*x80;
    double x125 = x123*x82;
    double x126 = -x124 - x125 + x40*x5*x83;
    double x127 = x123*x23;
    double x128 = x123*x45;
    double x129 = -x127 - x128 + x40*x48*x5;
    double x130 = x123*x51;
    double x131 = x123*x56;
    double x132 = -x130 - x131 + x40*x5*x59;
    double x133 = x2*x40;
    double x134 = x133*x83;
    double x135 = x124 + x134;
    double x136 = x133*x48;
    double x137 = x128 + x136;
    double x138 = x133*x59;
    double x139 = x131 + x138;
    double x140 = x133*x82;
    double x141 = -x134 + x14*x40*x80 - x140;
    double x142 = x133*x23;
    double x143 = -x136 + x14*x40*x45 - x142;
    double x144 = x133*x51;
    double x145 = -x138 + x14*x40*x56 - x144;
    double x146 = x3*x40;
    double x147 = x146*x80;
    double x148 = x146*x83;
    double x149 = -x147 - x148 + x25*x40*x82;
    double x150 = x146*x45;
    double x151 = x146*x48;
    double x152 = -x150 - x151 + x23*x25*x40;
    double x153 = x146*x56;
    double x154 = x146*x59;
    double x155 = -x153 - x154 + x25*x40*x51;
    double x156 = x127 + x151;
    double x157 = x130 + x154;
    double x158 = x125 + x148;
    double x159 = x142 + x150;
    double x160 = x144 + x153;
    double x161 = x140 + x147;
    
    res_0[0] = x39*(x104*x96 + x79*x84 + x94*x95);
    res_0[1] = x39*(x105*x95 + x106*x96 + x107*x84);
    res_0[2] = x39*(x108*x84 + x109*x96 + x110*x95);
    res_0[3] = x39*(x104*x113 + x112*x94 + x114*x79);
    res_0[4] = x39*(x105*x112 + x106*x113 + x107*x114);
    res_0[5] = x39*(x108*x114 + x109*x113 + x110*x112);
    res_0[6] = x39*(x104*x118 + x116*x79 + x117*x94);
    res_0[7] = x39*(x105*x117 + x106*x118 + x107*x116);
    res_0[8] = x39*(x108*x116 + x109*x118 + x110*x117);
    res_0[9] = x39*(x104*x121 + x120*x94 + x122*x79);
    res_0[10] = x39*(x105*x120 + x106*x121 + x107*x122);
    res_0[11] = x39*(x108*x122 + x109*x121 + x110*x120);
    res_0[12] = x39*(x104*x132 + x126*x79 + x129*x94);
    res_0[13] = x39*(x105*x129 + x106*x132 + x107*x126);
    res_0[14] = x39*(x108*x126 + x109*x132 + x110*x129);
    res_0[15] = x39*(x104*x139 + x135*x79 + x137*x94);
    res_0[16] = x39*(x105*x137 + x106*x139 + x107*x135);
    res_0[17] = x39*(x108*x135 + x109*x139 + x110*x137);
    res_0[18] = x39*(x104*x145 + x141*x79 + x143*x94);
    res_0[19] = x39*(x105*x143 + x106*x145 + x107*x141);
    res_0[20] = x39*(x108*x141 + x109*x145 + x110*x143);
    res_0[21] = x39*(x104*x155 + x149*x79 + x152*x94);
    res_0[22] = x39*(x105*x152 + x106*x155 + x107*x149);
    res_0[23] = x39*(x108*x149 + x109*x155 + x110*x152);
    res_0[24] = x39*(x104*x157 + x156*x94 + x158*x79);
    res_0[25] = x39*(x105*x156 + x106*x157 + x107*x158);
    res_0[26] = x39*(x108*x158 + x109*x157 + x110*x156);
    res_0[27] = x39*(x104*x160 + x159*x94 + x161*x79);
    res_0[28] = x39*(x105*x159 + x106*x160 + x107*x161);
    res_0[29] = x39*(x108*x161 + x109*x160 + x110*x159);
}

Conf_Forces_API void Integration_C3D10R_static_mbf(size_t num_elem,double Coords[][10][3],double Element_U[][10][3],double S[][1][6],
    double PENER[][1],double SENER[][1],double Conf_Force[][10][3]){
    
    //Inputs:
    //    Coords.........Element_Nodal_Coordinates [num_elem][numNodes,3]
    //    Element_U......Element_Nodal_Displacements [num_elem][numNodes,3]
    //    S..............Stress at integration points [num_elem][num_int_points,6]
    //    PENER..........Plastic strain energy [num_elem][num_int_points]
    //    SENER..........Elastic strain energy [num_elem][num_int_points]
    //    
    double weights[1]={0.16666666666666666,};
    double int_points[1][3]={
        {0.25,0.25,0.25,}};

    #pragma omp parallel for
    for (size_t i=0;i<num_elem;i++){
        double TMP[10][3];
        for (size_t j=0;j<1;j++){
            C3D10R_static_mbf((double*) int_points[j],(double*) Coords[i],(double*) Element_U[i],(double*) S[i][j],PENER[i][j],SENER[i][j],(double*) TMP);
            for (size_t num_node=0;num_node<10;num_node++){
                for (size_t ii=0;ii<3;ii++){
                    Conf_Force[i][num_node][ii]+=TMP[num_node][ii]*weights[j];
                }
            }
        }
    }
}

void C3D10R_static_dbf(double *rst,double *coord,double *Element_U,double *S,double PENER,double SENER,double *res_0){
    double r = rst[0];
    double s = rst[1];
    double t = rst[2];
    double *U = Element_U;
    double S11 = S[0];
    double S22 = S[1];
    double S33 = S[2];
    double S12 = S[3];
    double S13 = S[4];
    double S23 = S[5];
    
    double x0 = 4*r;
    double x1 = x0 - 1;
    double x2 = 4*s;
    double x3 = 4*t;
    double x4 = x2 + x3;
    double x5 = -8*r - x4 + 4;
    double x6 = -x2*coord[19];
    double x7 = x0 + x4 - 3;
    double x8 = x7*coord[1];
    double x9 = -x3*coord[22] + x8;
    double x10 = x1*coord[4] + x2*coord[16] + x3*coord[25] + x5*coord[13] + x6 + x9;
    double x11 = x2 - 1;
    double x12 = x7*coord[0];
    double x13 = x0 - 4;
    double x14 = -8*s - x13 - x3;
    double x15 = -x0*coord[12];
    double x16 = -x3*coord[21];
    double x17 = x0*coord[15] + x11*coord[6] + x12 + x14*coord[18] + x15 + x16 + x3*coord[27];
    double x18 = x10*x17;
    double x19 = x12 - x2*coord[18];
    double x20 = x1*coord[3] + x16 + x19 + x2*coord[15] + x3*coord[24] + x5*coord[12];
    double x21 = -x0*coord[13];
    double x22 = x0*coord[16] + x11*coord[7] + x14*coord[19] + x21 + x3*coord[28] + x9;
    double x23 = -x18 + x20*x22;
    double x24 = x3 - 1;
    double x25 = -8*t - x13 - x2;
    double x26 = -x2*coord[20];
    double x27 = x7*coord[2];
    double x28 = -x0*coord[14] + x27;
    double x29 = x0*coord[26] + x2*coord[29] + x24*coord[11] + x25*coord[23] + x26 + x28;
    double x30 = -x3*coord[23];
    double x31 = x0*coord[17] + x11*coord[8] + x14*coord[20] + x28 + x3*coord[29] + x30;
    double x32 = x0*coord[24] + x15 + x19 + x2*coord[27] + x24*coord[9] + x25*coord[21];
    double x33 = x10*x32;
    double x34 = x1*coord[5] + x2*coord[17] + x26 + x27 + x3*coord[26] + x30 + x5*coord[14];
    double x35 = x0*coord[25] + x2*coord[28] + x21 + x24*coord[10] + x25*coord[22] + x6 + x8;
    double x36 = x17*x35;
    double x37 = x20*x35;
    double x38 = x22*x32;
    double x39 = -x18*x29 + x20*x22*x29 + x31*x33 - x31*x37 + x34*x36 - x34*x38;
    double x40 = 1.0/x39;
    double x41 = -x2*U[19];
    double x42 = x7*U[1];
    double x43 = -x0*U[13] + x42;
    double x44 = x40*(x0*U[25] + x2*U[28] + x24*U[10] + x25*U[22] + x41 + x43);
    double x45 = x33 - x37;
    double x46 = -x3*U[22];
    double x47 = x40*(x0*U[16] + x11*U[7] + x14*U[19] + x3*U[28] + x43 + x46);
    double x48 = x36 - x38;
    double x49 = x40*(x1*U[4] + x2*U[16] + x3*U[25] + x41 + x42 + x46 + x5*U[13]);
    double x50 = x23*x44 + x45*x47 + x48*x49;
    double x51 = x10*x31 - x22*x34;
    double x52 = -x2*U[20];
    double x53 = x7*U[2];
    double x54 = -x0*U[14] + x53;
    double x55 = x40*(x0*U[26] + x2*U[29] + x24*U[11] + x25*U[23] + x52 + x54);
    double x56 = -x10*x29 + x34*x35;
    double x57 = -x3*U[23];
    double x58 = x40*(x0*U[17] + x11*U[8] + x14*U[20] + x3*U[29] + x54 + x57);
    double x59 = x22*x29 - x31*x35;
    double x60 = x40*(x1*U[5] + x2*U[17] + x3*U[26] + x5*U[14] + x52 + x53 + x57);
    double x61 = x51*x55 + x56*x58 + x59*x60;
    double x62 = x44*x51 + x47*x56 + x49*x59;
    double x63 = x23*x55 + x45*x58 + x48*x60;
    double x64 = x63 + 1.0;
    double x65 = x50*x61 - x62*x64;
    double x66 = -x2*U[18];
    double x67 = x7*U[0];
    double x68 = -x0*U[12] + x67;
    double x69 = x40*(x0*U[24] + x2*U[27] + x24*U[9] + x25*U[21] + x66 + x68);
    double x70 = -x3*U[21];
    double x71 = x40*(x0*U[15] + x11*U[6] + x14*U[18] + x3*U[27] + x68 + x70);
    double x72 = x40*(x1*U[3] + x2*U[15] + x3*U[24] + x5*U[12] + x66 + x67 + x70);
    double x73 = x23*x69 + x45*x71 + x48*x72;
    double x74 = x51*x69 + x56*x71 + x59*x72;
    double x75 = x74 + 1.0;
    double x76 = -x61*x73 + x64*x75;
    double x77 = -x50*x75 + x62*x73;
    double x78 = S11*x65 + S12*x76 + S13*x77;
    double x79 = S12*x65 + S22*x76 + S23*x77;
    double x80 = S13*x65 + S23*x76 + S33*x77;
    double x81 = -x61*x80 - x62*x79 - x74*x78;
    double x82 = x20*x29 - x32*x34;
    double x83 = x40*x7;
    double x84 = x17*x34 - x20*x31;
    double x85 = -x17*x29 + x31*x32;
    double x86 = x82*x83 + x83*x84 + x83*x85;
    double x87 = x55*x84 + x58*x82 + x60*x85;
    double x88 = x44*x84 + x47*x82 + x49*x85;
    double x89 = x88 + 1.0;
    double x90 = -x61*x89 + x62*x87;
    double x91 = x69*x84 + x71*x82 + x72*x85;
    double x92 = x61*x91 - x75*x87;
    double x93 = -x62*x91 + x75*x89;
    double x94 = S11*x90 + S12*x92 + S13*x93;
    double x95 = S12*x90 + S22*x92 + S23*x93;
    double x96 = S13*x90 + S23*x92 + S33*x93;
    double x97 = -x61*x96 - x62*x95 - x74*x94;
    double x98 = x23*x83 + x45*x83 + x48*x83;
    double x99 = x51*x83 + x56*x83 + x59*x83;
    double x100 = -x50*x87 + x64*x89;
    double x101 = -x64*x91 + x73*x87;
    double x102 = x50*x91 - x73*x89;
    double x103 = S11*x100 + S12*x101 + S13*x102;
    double x104 = S12*x100 + S22*x101 + S23*x102;
    double x105 = S13*x100 + S23*x101 + S33*x102;
    double x106 = -1.0*PENER - 1.0*SENER;
    double x107 = -x103*x74 - x104*x62 - x105*x61 - x106;
    double x108 = -x87*x96 - x88*x95 - x91*x94;
    double x109 = -x103*x91 - x104*x88 - x105*x87;
    double x110 = -x106 - x78*x91 - x79*x88 - x80*x87;
    double x111 = -x50*x79 - x63*x80 - x73*x78;
    double x112 = -x103*x73 - x104*x50 - x105*x63;
    double x113 = -x106 - x50*x95 - x63*x96 - x73*x94;
    double x114 = x1*x40;
    double x115 = x114*x48;
    double x116 = x114*x59;
    double x117 = x114*x85;
    double x118 = x11*x40;
    double x119 = x118*x82;
    double x120 = x118*x45;
    double x121 = x118*x56;
    double x122 = x24*x40;
    double x123 = x122*x23;
    double x124 = x122*x51;
    double x125 = x122*x84;
    double x126 = x0*x40;
    double x127 = x126*x82;
    double x128 = x126*x84;
    double x129 = -x127 - x128 + x40*x5*x85;
    double x130 = x126*x23;
    double x131 = x126*x45;
    double x132 = -x130 - x131 + x40*x48*x5;
    double x133 = x126*x51;
    double x134 = x126*x56;
    double x135 = -x133 - x134 + x40*x5*x59;
    double x136 = x2*x40;
    double x137 = x136*x85;
    double x138 = x127 + x137;
    double x139 = x136*x48;
    double x140 = x131 + x139;
    double x141 = x136*x59;
    double x142 = x134 + x141;
    double x143 = x136*x84;
    double x144 = -x137 + x14*x40*x82 - x143;
    double x145 = x136*x23;
    double x146 = -x139 + x14*x40*x45 - x145;
    double x147 = x136*x51;
    double x148 = x14*x40*x56 - x141 - x147;
    double x149 = x3*x40;
    double x150 = x149*x82;
    double x151 = x149*x85;
    double x152 = -x150 - x151 + x25*x40*x84;
    double x153 = x149*x45;
    double x154 = x149*x48;
    double x155 = -x153 - x154 + x23*x25*x40;
    double x156 = x149*x56;
    double x157 = x149*x59;
    double x158 = -x156 - x157 + x25*x40*x51;
    double x159 = x130 + x154;
    double x160 = x133 + x157;
    double x161 = x128 + x151;
    double x162 = x145 + x153;
    double x163 = x147 + x156;
    double x164 = x143 + x150;
    
    res_0[0] = x39*(x107*x99 + x81*x86 + x97*x98);
    res_0[1] = x39*(x108*x98 + x109*x99 + x110*x86);
    res_0[2] = x39*(x111*x86 + x112*x99 + x113*x98);
    res_0[3] = x39*(x107*x116 + x115*x97 + x117*x81);
    res_0[4] = x39*(x108*x115 + x109*x116 + x110*x117);
    res_0[5] = x39*(x111*x117 + x112*x116 + x113*x115);
    res_0[6] = x39*(x107*x121 + x119*x81 + x120*x97);
    res_0[7] = x39*(x108*x120 + x109*x121 + x110*x119);
    res_0[8] = x39*(x111*x119 + x112*x121 + x113*x120);
    res_0[9] = x39*(x107*x124 + x123*x97 + x125*x81);
    res_0[10] = x39*(x108*x123 + x109*x124 + x110*x125);
    res_0[11] = x39*(x111*x125 + x112*x124 + x113*x123);
    res_0[12] = x39*(x107*x135 + x129*x81 + x132*x97);
    res_0[13] = x39*(x108*x132 + x109*x135 + x110*x129);
    res_0[14] = x39*(x111*x129 + x112*x135 + x113*x132);
    res_0[15] = x39*(x107*x142 + x138*x81 + x140*x97);
    res_0[16] = x39*(x108*x140 + x109*x142 + x110*x138);
    res_0[17] = x39*(x111*x138 + x112*x142 + x113*x140);
    res_0[18] = x39*(x107*x148 + x144*x81 + x146*x97);
    res_0[19] = x39*(x108*x146 + x109*x148 + x110*x144);
    res_0[20] = x39*(x111*x144 + x112*x148 + x113*x146);
    res_0[21] = x39*(x107*x158 + x152*x81 + x155*x97);
    res_0[22] = x39*(x108*x155 + x109*x158 + x110*x152);
    res_0[23] = x39*(x111*x152 + x112*x158 + x113*x155);
    res_0[24] = x39*(x107*x160 + x159*x97 + x161*x81);
    res_0[25] = x39*(x108*x159 + x109*x160 + x110*x161);
    res_0[26] = x39*(x111*x161 + x112*x160 + x113*x159);
    res_0[27] = x39*(x107*x163 + x162*x97 + x164*x81);
    res_0[28] = x39*(x108*x162 + x109*x163 + x110*x164);
    res_0[29] = x39*(x111*x164 + x112*x163 + x113*x162);
}

Conf_Forces_API void Integration_C3D10R_static_dbf(size_t num_elem,double Coords[][10][3],double Element_U[][10][3],double S[][1][6],
    double PENER[][1],double SENER[][1],double Conf_Force[][10][3]){
    
    //Inputs:
    //    Coords.........Element_Nodal_Coordinates [num_elem][numNodes,3]
    //    Element_U......Element_Nodal_Displacements [num_elem][numNodes,3]
    //    S..............Stress at integration points [num_elem][num_int_points,6]
    //    PENER..........Plastic strain energy [num_elem][num_int_points]
    //    SENER..........Elastic strain energy [num_elem][num_int_points]
    //    
    double weights[1]={0.16666666666666666,};
    double int_points[1][3]={
        {0.25,0.25,0.25,}};

    #pragma omp parallel for
    for (size_t i=0;i<num_elem;i++){
        double TMP[10][3];
        for (size_t j=0;j<1;j++){
            C3D10R_static_dbf((double*) int_points[j],(double*) Coords[i],(double*) Element_U[i],(double*) S[i][j],PENER[i][j],SENER[i][j],(double*) TMP);
            for (size_t num_node=0;num_node<10;num_node++){
                for (size_t ii=0;ii<3;ii++){
                    Conf_Force[i][num_node][ii]+=TMP[num_node][ii]*weights[j];
                }
            }
        }
    }
}

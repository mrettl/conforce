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
    double x1 = -x0;
    double x2 = x1 - 1.0/4.0;
    double *x3 = coord;
    double x4 = x3[4];
    double x5 = x0 - 1.0/4.0;
    double x6 = x3[1];
    double x7 = x1 + 1.0/4.0;
    double x8 = x3[10];
    double x9 = x0 + 1.0/4.0;
    double x10 = x3[7];
    double x11 = x10*x9 + x2*x4 + x5*x6 + x7*x8;
    double x12 = (1.0/4.0)*s;
    double x13 = -x12;
    double x14 = x13 - 1.0/4.0;
    double x15 = x3[9];
    double x16 = x12 - 1.0/4.0;
    double x17 = x3[0];
    double x18 = x13 + 1.0/4.0;
    double x19 = x3[3];
    double x20 = x12 + 1.0/4.0;
    double x21 = x3[6];
    double x22 = x14*x15 + x16*x17 + x18*x19 + x20*x21;
    double x23 = x19*x2;
    double x24 = x17*x5;
    double x25 = x15*x7;
    double x26 = x21*x9;
    double x27 = x14*x8;
    double x28 = x16*x6;
    double x29 = x18*x4;
    double x30 = x10*x20;
    double x31 = x11*x22 - (x23 + x24 + x25 + x26)*(x27 + x28 + x29 + x30);
    double x32 = 1.0/x31;
    double x33 = x32*x5;
    double x34 = -x23 - x24 - x25 - x26;
    double x35 = x16*x32;
    double x36 = x22*x33 + x34*x35;
    double *x37 = U;
    double x38 = x37[4];
    double x39 = x37[1];
    double x40 = x37[10];
    double x41 = x37[7];
    double x42 = x2*x38 + x39*x5 + x40*x7 + x41*x9;
    double x43 = -x27 - x28 - x29 - x30;
    double x44 = x32*x43;
    double x45 = x42*x44;
    double x46 = x14*x40 + x16*x39 + x18*x38 + x20*x41;
    double x47 = x11*x32;
    double x48 = x46*x47;
    double x49 = -1.0*x45 - 1.0*x48;
    double x50 = x37[3];
    double x51 = x37[0];
    double x52 = x37[9];
    double x53 = x37[6];
    double x54 = x2*x50 + x5*x51 + x52*x7 + x53*x9;
    double x55 = x44*x54;
    double x56 = x14*x52 + x16*x51 + x18*x50 + x20*x53;
    double x57 = x47*x56;
    double x58 = 1.0*x55 + 1.0*x57 + 1.0;
    double x59 = S11*x49 + S12*x58;
    double x60 = x55 + x57;
    double x61 = S12*x49 + S22*x58;
    double x62 = x45 + x48;
    double x63 = S13*x49 + S23*x58;
    double x64 = x37[5];
    double x65 = x37[2];
    double x66 = x37[11];
    double x67 = x37[8];
    double x68 = x2*x64 + x5*x65 + x66*x7 + x67*x9;
    double x69 = x14*x66 + x16*x65 + x18*x64 + x20*x67;
    double x70 = x44*x68 + x47*x69;
    double x71 = -x59*x60 - x61*x62 - x63*x70;
    double x72 = x11*x35 + x33*x43;
    double x73 = x22*x32;
    double x74 = x42*x73;
    double x75 = x32*x34;
    double x76 = x46*x75;
    double x77 = 1.0*x74 + 1.0*x76 + 1.0;
    double x78 = x54*x73;
    double x79 = x56*x75;
    double x80 = -1.0*x78 - 1.0*x79;
    double x81 = S11*x77 + S12*x80;
    double x82 = S12*x77 + S22*x80;
    double x83 = S13*x77 + S23*x80;
    double x84 = s*x0;
    double x85 = x13 + x7 + x84;
    double *x86 = V;
    double x87 = x86[0];
    double x88 = -x84;
    double x89 = x12 + x7 + x88;
    double x90 = x86[9];
    double x91 = x13 + x88 + x9;
    double x92 = x86[3];
    double x93 = x12 + x84 + x9;
    double x94 = x86[6];
    double x95 = x85*x87 + x89*x90 + x91*x92 + x93*x94;
    double x96 = x86[1];
    double x97 = x86[10];
    double x98 = x86[4];
    double x99 = x86[7];
    double x100 = x85*x96 + x89*x97 + x91*x98 + x93*x99;
    double x101 = x86[2];
    double x102 = x86[11];
    double x103 = x86[5];
    double x104 = x86[8];
    double x105 = x101*x85 + x102*x89 + x103*x91 + x104*x93;
    double x106 = 1.0*PENER + 1.0*SENER - 0.5*rho*(x100*x100 + x105*x105 + x95*x95);
    double x107 = x106 - x60*x81 - x62*x82 - x70*x83;
    double *x108 = A;
    double x109 = x85*x108[0] + x89*x108[9] + x91*x108[3] + x93*x108[6];
    double x110 = x85*x108[1] + x89*x108[10] + x91*x108[4] + x93*x108[7];
    double x111 = x85*x108[2] + x89*x108[11] + x91*x108[5] + x93*x108[8];
    double x112 = x2*x92 + x5*x87 + x7*x90 + x9*x94;
    double x113 = x14*x90 + x16*x87 + x18*x92 + x20*x94;
    double x114 = x2*x98 + x5*x96 + x7*x97 + x9*x99;
    double x115 = x14*x97 + x16*x96 + x18*x98 + x20*x99;
    double x116 = x101*x5 + x102*x7 + x103*x2 + x104*x9;
    double x117 = x101*x16 + x102*x14 + x103*x18 + x104*x20;
    double x118 = x100*(x114*x44 + x115*x47) + x105*(x116*x44 + x117*x47) + x109*x60 + x110*x62 + x111*x70 + x95*(x112*x44 + x113*x47);
    double x119 = rho*x85;
    double x120 = x78 + x79;
    double x121 = x74 + x76;
    double x122 = x68*x73 + x69*x75;
    double x123 = x106 - x120*x59 - x121*x61 - x122*x63;
    double x124 = -x120*x81 - x121*x82 - x122*x83;
    double x125 = x100*(x114*x73 + x115*x75) + x105*(x116*x73 + x117*x75) + x109*x120 + x110*x121 + x111*x122 + x95*(x112*x73 + x113*x75);
    double x126 = x18*x75 + x2*x73;
    double x127 = x18*x47 + x2*x44;
    double x128 = rho*x91;
    double x129 = x20*x75 + x73*x9;
    double x130 = x20*x47 + x44*x9;
    double x131 = rho*x93;
    double x132 = x14*x47 + x44*x7;
    double x133 = x14*x75 + x7*x73;
    double x134 = rho*x89;
    
    res_0[0] = x31*(x107*x72 - x118*x119 + x36*x71);
    res_0[1] = x31*(-x119*x125 + x123*x36 + x124*x72);
    res_0[2] = 0;
    res_0[3] = x31*(x107*x127 - x118*x128 + x126*x71);
    res_0[4] = x31*(x123*x126 + x124*x127 - x125*x128);
    res_0[5] = 0;
    res_0[6] = x31*(x107*x130 - x118*x131 + x129*x71);
    res_0[7] = x31*(x123*x129 + x124*x130 - x125*x131);
    res_0[8] = 0;
    res_0[9] = x31*(x107*x132 - x118*x134 + x133*x71);
    res_0[10] = x31*(x123*x133 + x124*x132 - x125*x134);
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
    double x1 = -x0;
    double x2 = x1 - 1.0/4.0;
    double *x3 = coord;
    double x4 = x3[9];
    double x5 = x0 - 1.0/4.0;
    double x6 = x3[0];
    double x7 = x1 + 1.0/4.0;
    double x8 = x3[3];
    double x9 = x0 + 1.0/4.0;
    double x10 = x3[6];
    double x11 = x10*x9 + x2*x4 + x5*x6 + x7*x8;
    double x12 = (1.0/4.0)*r;
    double x13 = x12 - 1.0/4.0;
    double x14 = -x12;
    double x15 = x14 - 1.0/4.0;
    double x16 = x3[4];
    double x17 = x3[1];
    double x18 = x14 + 1.0/4.0;
    double x19 = x3[10];
    double x20 = x12 + 1.0/4.0;
    double x21 = x3[7];
    double x22 = x13*x17 + x15*x16 + x18*x19 + x20*x21;
    double x23 = x15*x8;
    double x24 = x13*x6;
    double x25 = x18*x4;
    double x26 = x10*x20;
    double x27 = x19*x2;
    double x28 = x17*x5;
    double x29 = x16*x7;
    double x30 = x21*x9;
    double x31 = x11*x22 - (x23 + x24 + x25 + x26)*(x27 + x28 + x29 + x30);
    double x32 = 1.0/x31;
    double x33 = x13*x32;
    double x34 = -x23 - x24 - x25 - x26;
    double x35 = x32*x5;
    double x36 = x11*x33 + x34*x35;
    double *x37 = U;
    double x38 = x37[4];
    double x39 = x37[1];
    double x40 = x37[10];
    double x41 = x37[7];
    double x42 = x13*x39 + x15*x38 + x18*x40 + x20*x41;
    double x43 = -x27 - x28 - x29 - x30;
    double x44 = x32*x43;
    double x45 = x42*x44;
    double x46 = x2*x40 + x38*x7 + x39*x5 + x41*x9;
    double x47 = x22*x32;
    double x48 = x46*x47;
    double x49 = -1.0*x45 - 1.0*x48;
    double x50 = x37[3];
    double x51 = x37[0];
    double x52 = x37[9];
    double x53 = x37[6];
    double x54 = x13*x51 + x15*x50 + x18*x52 + x20*x53;
    double x55 = x44*x54;
    double x56 = x2*x52 + x5*x51 + x50*x7 + x53*x9;
    double x57 = x47*x56;
    double x58 = 1.0*x55 + 1.0*x57 + 1.0;
    double x59 = S11*x49 + S12*x58;
    double x60 = x55 + x57 + 1.0;
    double x61 = S12*x49 + S22*x58;
    double x62 = x45 + x48;
    double x63 = S13*x49;
    double x64 = S23*x58;
    double x65 = x63 + x64;
    double x66 = x37[5];
    double x67 = x37[2];
    double x68 = x37[11];
    double x69 = x37[8];
    double x70 = x13*x67 + x15*x66 + x18*x68 + x20*x69;
    double x71 = x2*x68 + x5*x67 + x66*x7 + x69*x9;
    double x72 = x44*x70 + x47*x71;
    double x73 = -x59*x60 - x61*x62 - x65*x72;
    double x74 = x22*x35 + x33*x43;
    double x75 = x11*x32;
    double x76 = x42*x75;
    double x77 = x32*x34;
    double x78 = x46*x77;
    double x79 = 1.0*x76 + 1.0*x78 + 1.0;
    double x80 = x54*x75;
    double x81 = x56*x77;
    double x82 = -1.0*x80 - 1.0*x81;
    double x83 = S11*x79 + S12*x82;
    double x84 = S12*x79 + S22*x82;
    double x85 = S13*x79;
    double x86 = S23*x82;
    double x87 = x85 + x86;
    double x88 = 1.0*PENER + 1.0*SENER;
    double x89 = -x60*x83 - x62*x84 - x72*x87 + x88;
    double x90 = x80 + x81;
    double x91 = x76 + x78 + 1.0;
    double x92 = x70*x75 + x71*x77;
    double x93 = -x59*x90 - x61*x91 - x65*x92 + x88;
    double x94 = -x83*x90 - x84*x91 - x87*x92;
    double x95 = -1.0*x63 - 1.0*x64;
    double x96 = -1.0*x85 - 1.0*x86;
    double x97 = x15*x75 + x7*x77;
    double x98 = x15*x44 + x47*x7;
    double x99 = x20*x75 + x77*x9;
    double x100 = x20*x44 + x47*x9;
    double x101 = x18*x44 + x2*x47;
    double x102 = x18*x75 + x2*x77;
    
    res_0[0] = x31*(x36*x73 + x74*x89);
    res_0[1] = x31*(x36*x93 + x74*x94);
    res_0[2] = x31*(x36*x95 + x74*x96);
    res_0[3] = x31*(x73*x97 + x89*x98);
    res_0[4] = x31*(x93*x97 + x94*x98);
    res_0[5] = x31*(x95*x97 + x96*x98);
    res_0[6] = x31*(x100*x89 + x73*x99);
    res_0[7] = x31*(x100*x94 + x93*x99);
    res_0[8] = x31*(x100*x96 + x95*x99);
    res_0[9] = x31*(x101*x89 + x102*x73);
    res_0[10] = x31*(x101*x94 + x102*x93);
    res_0[11] = x31*(x101*x96 + x102*x95);
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
    double x1 = -x0;
    double x2 = x1 - 1.0/4.0;
    double *x3 = coord;
    double x4 = x3[9];
    double x5 = x0 - 1.0/4.0;
    double x6 = x3[0];
    double x7 = x1 + 1.0/4.0;
    double x8 = x3[3];
    double x9 = x0 + 1.0/4.0;
    double x10 = x3[6];
    double x11 = x10*x9 + x2*x4 + x5*x6 + x7*x8;
    double x12 = (1.0/4.0)*r;
    double x13 = x12 - 1.0/4.0;
    double x14 = -x12;
    double x15 = x14 - 1.0/4.0;
    double x16 = x3[4];
    double x17 = x3[1];
    double x18 = x14 + 1.0/4.0;
    double x19 = x3[10];
    double x20 = x12 + 1.0/4.0;
    double x21 = x3[7];
    double x22 = x13*x17 + x15*x16 + x18*x19 + x20*x21;
    double x23 = x15*x8;
    double x24 = x13*x6;
    double x25 = x18*x4;
    double x26 = x10*x20;
    double x27 = x19*x2;
    double x28 = x17*x5;
    double x29 = x16*x7;
    double x30 = x21*x9;
    double x31 = x11*x22 - (x23 + x24 + x25 + x26)*(x27 + x28 + x29 + x30);
    double x32 = 1.0/x31;
    double x33 = x13*x32;
    double x34 = -x23 - x24 - x25 - x26;
    double x35 = x32*x5;
    double x36 = x11*x33 + x34*x35;
    double *x37 = U;
    double x38 = x37[4];
    double x39 = x37[1];
    double x40 = x37[10];
    double x41 = x37[7];
    double x42 = x13*x39 + x15*x38 + x18*x40 + x20*x41;
    double x43 = -x27 - x28 - x29 - x30;
    double x44 = x32*x43;
    double x45 = x42*x44;
    double x46 = x2*x40 + x38*x7 + x39*x5 + x41*x9;
    double x47 = x22*x32;
    double x48 = x46*x47;
    double x49 = -1.0*x45 - 1.0*x48;
    double x50 = x37[3];
    double x51 = x37[0];
    double x52 = x37[9];
    double x53 = x37[6];
    double x54 = x13*x51 + x15*x50 + x18*x52 + x20*x53;
    double x55 = x44*x54;
    double x56 = x2*x52 + x5*x51 + x50*x7 + x53*x9;
    double x57 = x47*x56;
    double x58 = 1.0*x55 + 1.0*x57 + 1.0;
    double x59 = S11*x49 + S12*x58;
    double x60 = x55 + x57;
    double x61 = S12*x49 + S22*x58;
    double x62 = x45 + x48;
    double x63 = S13*x49 + S23*x58;
    double x64 = x37[5];
    double x65 = x37[2];
    double x66 = x37[11];
    double x67 = x37[8];
    double x68 = x13*x65 + x15*x64 + x18*x66 + x20*x67;
    double x69 = x2*x66 + x5*x65 + x64*x7 + x67*x9;
    double x70 = x44*x68 + x47*x69;
    double x71 = -x59*x60 - x61*x62 - x63*x70;
    double x72 = x22*x35 + x33*x43;
    double x73 = x11*x32;
    double x74 = x42*x73;
    double x75 = x32*x34;
    double x76 = x46*x75;
    double x77 = 1.0*x74 + 1.0*x76 + 1.0;
    double x78 = x54*x73;
    double x79 = x56*x75;
    double x80 = -1.0*x78 - 1.0*x79;
    double x81 = S11*x77 + S12*x80;
    double x82 = S12*x77 + S22*x80;
    double x83 = S13*x77 + S23*x80;
    double x84 = 1.0*PENER + 1.0*SENER;
    double x85 = -x60*x81 - x62*x82 - x70*x83 + x84;
    double x86 = x78 + x79;
    double x87 = x74 + x76;
    double x88 = x68*x73 + x69*x75;
    double x89 = -x59*x86 - x61*x87 - x63*x88 + x84;
    double x90 = -x81*x86 - x82*x87 - x83*x88;
    double x91 = x15*x73 + x7*x75;
    double x92 = x15*x44 + x47*x7;
    double x93 = x20*x73 + x75*x9;
    double x94 = x20*x44 + x47*x9;
    double x95 = x18*x44 + x2*x47;
    double x96 = x18*x73 + x2*x75;
    
    res_0[0] = x31*(x36*x71 + x72*x85);
    res_0[1] = x31*(x36*x89 + x72*x90);
    res_0[2] = 0;
    res_0[3] = x31*(x71*x91 + x85*x92);
    res_0[4] = x31*(x89*x91 + x90*x92);
    res_0[5] = 0;
    res_0[6] = x31*(x71*x93 + x85*x94);
    res_0[7] = x31*(x89*x93 + x90*x94);
    res_0[8] = 0;
    res_0[9] = x31*(x71*x96 + x85*x95);
    res_0[10] = x31*(x89*x96 + x90*x95);
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
    double x1 = -x0;
    double x2 = x1 - 1.0/4.0;
    double *x3 = coord;
    double x4 = x3[4];
    double x5 = x0 - 1.0/4.0;
    double x6 = x3[1];
    double x7 = x1 + 1.0/4.0;
    double x8 = x3[10];
    double x9 = x0 + 1.0/4.0;
    double x10 = x3[7];
    double x11 = x10*x9 + x2*x4 + x5*x6 + x7*x8;
    double x12 = (1.0/4.0)*s;
    double x13 = -x12;
    double x14 = x13 - 1.0/4.0;
    double x15 = x3[9];
    double x16 = x12 - 1.0/4.0;
    double x17 = x3[0];
    double x18 = x13 + 1.0/4.0;
    double x19 = x3[3];
    double x20 = x12 + 1.0/4.0;
    double x21 = x3[6];
    double x22 = x14*x15 + x16*x17 + x18*x19 + x20*x21;
    double x23 = x19*x2;
    double x24 = x17*x5;
    double x25 = x15*x7;
    double x26 = x21*x9;
    double x27 = x14*x8;
    double x28 = x16*x6;
    double x29 = x18*x4;
    double x30 = x10*x20;
    double x31 = x11*x22 - (x23 + x24 + x25 + x26)*(x27 + x28 + x29 + x30);
    double x32 = 1.0/x31;
    double x33 = x32*x5;
    double x34 = -x23 - x24 - x25 - x26;
    double x35 = x16*x32;
    double x36 = x22*x33 + x34*x35;
    double *x37 = U;
    double x38 = x37[4];
    double x39 = x37[1];
    double x40 = x37[10];
    double x41 = x37[7];
    double x42 = x2*x38 + x39*x5 + x40*x7 + x41*x9;
    double x43 = -x27 - x28 - x29 - x30;
    double x44 = x32*x43;
    double x45 = x42*x44;
    double x46 = x14*x40 + x16*x39 + x18*x38 + x20*x41;
    double x47 = x11*x32;
    double x48 = x46*x47;
    double x49 = -1.0*x45 - 1.0*x48;
    double x50 = x37[3];
    double x51 = x37[0];
    double x52 = x37[9];
    double x53 = x37[6];
    double x54 = x2*x50 + x5*x51 + x52*x7 + x53*x9;
    double x55 = x44*x54;
    double x56 = x14*x52 + x16*x51 + x18*x50 + x20*x53;
    double x57 = x47*x56;
    double x58 = 1.0*x55 + 1.0*x57 + 1.0;
    double x59 = S11*x49 + S12*x58;
    double x60 = x55 + x57;
    double x61 = S12*x49 + S22*x58;
    double x62 = x45 + x48;
    double x63 = S13*x49 + S23*x58;
    double x64 = x37[5];
    double x65 = x37[2];
    double x66 = x37[11];
    double x67 = x37[8];
    double x68 = x2*x64 + x5*x65 + x66*x7 + x67*x9;
    double x69 = x14*x66 + x16*x65 + x18*x64 + x20*x67;
    double x70 = x44*x68 + x47*x69;
    double x71 = -x59*x60 - x61*x62 - x63*x70;
    double x72 = x11*x35 + x33*x43;
    double x73 = x22*x32;
    double x74 = x42*x73;
    double x75 = x32*x34;
    double x76 = x46*x75;
    double x77 = 1.0*x74 + 1.0*x76 + 1.0;
    double x78 = x54*x73;
    double x79 = x56*x75;
    double x80 = -1.0*x78 - 1.0*x79;
    double x81 = S11*x77 + S12*x80;
    double x82 = S12*x77 + S22*x80;
    double x83 = S13*x77 + S23*x80;
    double x84 = s*x0;
    double x85 = x13 + x7 + x84;
    double *x86 = V;
    double x87 = x86[0];
    double x88 = -x84;
    double x89 = x12 + x7 + x88;
    double x90 = x86[9];
    double x91 = x13 + x88 + x9;
    double x92 = x86[3];
    double x93 = x12 + x84 + x9;
    double x94 = x86[6];
    double x95 = x85*x87 + x89*x90 + x91*x92 + x93*x94;
    double x96 = x86[1];
    double x97 = x86[10];
    double x98 = x86[4];
    double x99 = x86[7];
    double x100 = x85*x96 + x89*x97 + x91*x98 + x93*x99;
    double x101 = x86[2];
    double x102 = x86[11];
    double x103 = x86[5];
    double x104 = x86[8];
    double x105 = x101*x85 + x102*x89 + x103*x91 + x104*x93;
    double x106 = 1.0*PENER + 1.0*SENER - 0.5*rho*(x100*x100 + x105*x105 + x95*x95);
    double x107 = x106 - x60*x81 - x62*x82 - x70*x83;
    double *x108 = A;
    double x109 = x85*x108[0] + x89*x108[9] + x91*x108[3] + x93*x108[6];
    double x110 = x85*x108[1] + x89*x108[10] + x91*x108[4] + x93*x108[7];
    double x111 = x85*x108[2] + x89*x108[11] + x91*x108[5] + x93*x108[8];
    double x112 = x2*x92 + x5*x87 + x7*x90 + x9*x94;
    double x113 = x14*x90 + x16*x87 + x18*x92 + x20*x94;
    double x114 = x2*x98 + x5*x96 + x7*x97 + x9*x99;
    double x115 = x14*x97 + x16*x96 + x18*x98 + x20*x99;
    double x116 = x101*x5 + x102*x7 + x103*x2 + x104*x9;
    double x117 = x101*x16 + x102*x14 + x103*x18 + x104*x20;
    double x118 = x100*(x114*x44 + x115*x47) + x105*(x116*x44 + x117*x47) + x109*x60 + x110*x62 + x111*x70 + x95*(x112*x44 + x113*x47);
    double x119 = rho*x85;
    double x120 = x78 + x79;
    double x121 = x74 + x76;
    double x122 = x68*x73 + x69*x75;
    double x123 = x106 - x120*x59 - x121*x61 - x122*x63;
    double x124 = -x120*x81 - x121*x82 - x122*x83;
    double x125 = x100*(x114*x73 + x115*x75) + x105*(x116*x73 + x117*x75) + x109*x120 + x110*x121 + x111*x122 + x95*(x112*x73 + x113*x75);
    double x126 = x18*x75 + x2*x73;
    double x127 = x18*x47 + x2*x44;
    double x128 = rho*x91;
    double x129 = x20*x75 + x73*x9;
    double x130 = x20*x47 + x44*x9;
    double x131 = rho*x93;
    double x132 = x14*x47 + x44*x7;
    double x133 = x14*x75 + x7*x73;
    double x134 = rho*x89;
    
    res_0[0] = x31*(x107*x72 - x118*x119 + x36*x71);
    res_0[1] = x31*(-x119*x125 + x123*x36 + x124*x72);
    res_0[2] = 0;
    res_0[3] = x31*(x107*x127 - x118*x128 + x126*x71);
    res_0[4] = x31*(x123*x126 + x124*x127 - x125*x128);
    res_0[5] = 0;
    res_0[6] = x31*(x107*x130 - x118*x131 + x129*x71);
    res_0[7] = x31*(x123*x129 + x124*x130 - x125*x131);
    res_0[8] = 0;
    res_0[9] = x31*(x107*x132 - x118*x134 + x133*x71);
    res_0[10] = x31*(x123*x133 + x124*x132 - x125*x134);
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
    double x1 = -x0;
    double x2 = x1 - 1.0/4.0;
    double *x3 = coord;
    double x4 = x3[9];
    double x5 = x0 - 1.0/4.0;
    double x6 = x3[0];
    double x7 = x1 + 1.0/4.0;
    double x8 = x3[3];
    double x9 = x0 + 1.0/4.0;
    double x10 = x3[6];
    double x11 = x10*x9 + x2*x4 + x5*x6 + x7*x8;
    double x12 = (1.0/4.0)*r;
    double x13 = x12 - 1.0/4.0;
    double x14 = -x12;
    double x15 = x14 - 1.0/4.0;
    double x16 = x3[4];
    double x17 = x3[1];
    double x18 = x14 + 1.0/4.0;
    double x19 = x3[10];
    double x20 = x12 + 1.0/4.0;
    double x21 = x3[7];
    double x22 = x13*x17 + x15*x16 + x18*x19 + x20*x21;
    double x23 = x15*x8;
    double x24 = x13*x6;
    double x25 = x18*x4;
    double x26 = x10*x20;
    double x27 = x19*x2;
    double x28 = x17*x5;
    double x29 = x16*x7;
    double x30 = x21*x9;
    double x31 = x11*x22 - (x23 + x24 + x25 + x26)*(x27 + x28 + x29 + x30);
    double x32 = 1.0/x31;
    double x33 = x13*x32;
    double x34 = -x23 - x24 - x25 - x26;
    double x35 = x32*x5;
    double x36 = x11*x33 + x34*x35;
    double *x37 = U;
    double x38 = x37[4];
    double x39 = x37[1];
    double x40 = x37[10];
    double x41 = x37[7];
    double x42 = x13*x39 + x15*x38 + x18*x40 + x20*x41;
    double x43 = -x27 - x28 - x29 - x30;
    double x44 = x32*x43;
    double x45 = x42*x44;
    double x46 = x2*x40 + x38*x7 + x39*x5 + x41*x9;
    double x47 = x22*x32;
    double x48 = x46*x47;
    double x49 = -1.0*x45 - 1.0*x48;
    double x50 = x37[3];
    double x51 = x37[0];
    double x52 = x37[9];
    double x53 = x37[6];
    double x54 = x13*x51 + x15*x50 + x18*x52 + x20*x53;
    double x55 = x44*x54;
    double x56 = x2*x52 + x5*x51 + x50*x7 + x53*x9;
    double x57 = x47*x56;
    double x58 = 1.0*x55 + 1.0*x57 + 1.0;
    double x59 = S11*x49 + S12*x58;
    double x60 = x55 + x57 + 1.0;
    double x61 = S12*x49 + S22*x58;
    double x62 = x45 + x48;
    double x63 = S13*x49;
    double x64 = S23*x58;
    double x65 = x63 + x64;
    double x66 = x37[5];
    double x67 = x37[2];
    double x68 = x37[11];
    double x69 = x37[8];
    double x70 = x13*x67 + x15*x66 + x18*x68 + x20*x69;
    double x71 = x2*x68 + x5*x67 + x66*x7 + x69*x9;
    double x72 = x44*x70 + x47*x71;
    double x73 = -x59*x60 - x61*x62 - x65*x72;
    double x74 = x22*x35 + x33*x43;
    double x75 = x11*x32;
    double x76 = x42*x75;
    double x77 = x32*x34;
    double x78 = x46*x77;
    double x79 = 1.0*x76 + 1.0*x78 + 1.0;
    double x80 = x54*x75;
    double x81 = x56*x77;
    double x82 = -1.0*x80 - 1.0*x81;
    double x83 = S11*x79 + S12*x82;
    double x84 = S12*x79 + S22*x82;
    double x85 = S13*x79;
    double x86 = S23*x82;
    double x87 = x85 + x86;
    double x88 = 1.0*PENER + 1.0*SENER;
    double x89 = -x60*x83 - x62*x84 - x72*x87 + x88;
    double x90 = x80 + x81;
    double x91 = x76 + x78 + 1.0;
    double x92 = x70*x75 + x71*x77;
    double x93 = -x59*x90 - x61*x91 - x65*x92 + x88;
    double x94 = -x83*x90 - x84*x91 - x87*x92;
    double x95 = -1.0*x63 - 1.0*x64;
    double x96 = -1.0*x85 - 1.0*x86;
    double x97 = x15*x75 + x7*x77;
    double x98 = x15*x44 + x47*x7;
    double x99 = x20*x75 + x77*x9;
    double x100 = x20*x44 + x47*x9;
    double x101 = x18*x44 + x2*x47;
    double x102 = x18*x75 + x2*x77;
    
    res_0[0] = x31*(x36*x73 + x74*x89);
    res_0[1] = x31*(x36*x93 + x74*x94);
    res_0[2] = x31*(x36*x95 + x74*x96);
    res_0[3] = x31*(x73*x97 + x89*x98);
    res_0[4] = x31*(x93*x97 + x94*x98);
    res_0[5] = x31*(x95*x97 + x96*x98);
    res_0[6] = x31*(x100*x89 + x73*x99);
    res_0[7] = x31*(x100*x94 + x93*x99);
    res_0[8] = x31*(x100*x96 + x95*x99);
    res_0[9] = x31*(x101*x89 + x102*x73);
    res_0[10] = x31*(x101*x94 + x102*x93);
    res_0[11] = x31*(x101*x96 + x102*x95);
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
    double x1 = -x0;
    double x2 = x1 - 1.0/4.0;
    double *x3 = coord;
    double x4 = x3[9];
    double x5 = x0 - 1.0/4.0;
    double x6 = x3[0];
    double x7 = x1 + 1.0/4.0;
    double x8 = x3[3];
    double x9 = x0 + 1.0/4.0;
    double x10 = x3[6];
    double x11 = x10*x9 + x2*x4 + x5*x6 + x7*x8;
    double x12 = (1.0/4.0)*r;
    double x13 = x12 - 1.0/4.0;
    double x14 = -x12;
    double x15 = x14 - 1.0/4.0;
    double x16 = x3[4];
    double x17 = x3[1];
    double x18 = x14 + 1.0/4.0;
    double x19 = x3[10];
    double x20 = x12 + 1.0/4.0;
    double x21 = x3[7];
    double x22 = x13*x17 + x15*x16 + x18*x19 + x20*x21;
    double x23 = x15*x8;
    double x24 = x13*x6;
    double x25 = x18*x4;
    double x26 = x10*x20;
    double x27 = x19*x2;
    double x28 = x17*x5;
    double x29 = x16*x7;
    double x30 = x21*x9;
    double x31 = x11*x22 - (x23 + x24 + x25 + x26)*(x27 + x28 + x29 + x30);
    double x32 = 1.0/x31;
    double x33 = x13*x32;
    double x34 = -x23 - x24 - x25 - x26;
    double x35 = x32*x5;
    double x36 = x11*x33 + x34*x35;
    double *x37 = U;
    double x38 = x37[4];
    double x39 = x37[1];
    double x40 = x37[10];
    double x41 = x37[7];
    double x42 = x13*x39 + x15*x38 + x18*x40 + x20*x41;
    double x43 = -x27 - x28 - x29 - x30;
    double x44 = x32*x43;
    double x45 = x42*x44;
    double x46 = x2*x40 + x38*x7 + x39*x5 + x41*x9;
    double x47 = x22*x32;
    double x48 = x46*x47;
    double x49 = -1.0*x45 - 1.0*x48;
    double x50 = x37[3];
    double x51 = x37[0];
    double x52 = x37[9];
    double x53 = x37[6];
    double x54 = x13*x51 + x15*x50 + x18*x52 + x20*x53;
    double x55 = x44*x54;
    double x56 = x2*x52 + x5*x51 + x50*x7 + x53*x9;
    double x57 = x47*x56;
    double x58 = 1.0*x55 + 1.0*x57 + 1.0;
    double x59 = S11*x49 + S12*x58;
    double x60 = x55 + x57;
    double x61 = S12*x49 + S22*x58;
    double x62 = x45 + x48;
    double x63 = S13*x49 + S23*x58;
    double x64 = x37[5];
    double x65 = x37[2];
    double x66 = x37[11];
    double x67 = x37[8];
    double x68 = x13*x65 + x15*x64 + x18*x66 + x20*x67;
    double x69 = x2*x66 + x5*x65 + x64*x7 + x67*x9;
    double x70 = x44*x68 + x47*x69;
    double x71 = -x59*x60 - x61*x62 - x63*x70;
    double x72 = x22*x35 + x33*x43;
    double x73 = x11*x32;
    double x74 = x42*x73;
    double x75 = x32*x34;
    double x76 = x46*x75;
    double x77 = 1.0*x74 + 1.0*x76 + 1.0;
    double x78 = x54*x73;
    double x79 = x56*x75;
    double x80 = -1.0*x78 - 1.0*x79;
    double x81 = S11*x77 + S12*x80;
    double x82 = S12*x77 + S22*x80;
    double x83 = S13*x77 + S23*x80;
    double x84 = 1.0*PENER + 1.0*SENER;
    double x85 = -x60*x81 - x62*x82 - x70*x83 + x84;
    double x86 = x78 + x79;
    double x87 = x74 + x76;
    double x88 = x68*x73 + x69*x75;
    double x89 = -x59*x86 - x61*x87 - x63*x88 + x84;
    double x90 = -x81*x86 - x82*x87 - x83*x88;
    double x91 = x15*x73 + x7*x75;
    double x92 = x15*x44 + x47*x7;
    double x93 = x20*x73 + x75*x9;
    double x94 = x20*x44 + x47*x9;
    double x95 = x18*x44 + x2*x47;
    double x96 = x18*x73 + x2*x75;
    
    res_0[0] = x31*(x36*x71 + x72*x85);
    res_0[1] = x31*(x36*x89 + x72*x90);
    res_0[2] = 0;
    res_0[3] = x31*(x71*x91 + x85*x92);
    res_0[4] = x31*(x89*x91 + x90*x92);
    res_0[5] = 0;
    res_0[6] = x31*(x71*x93 + x85*x94);
    res_0[7] = x31*(x89*x93 + x90*x94);
    res_0[8] = 0;
    res_0[9] = x31*(x71*x96 + x85*x95);
    res_0[10] = x31*(x89*x96 + x90*x95);
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
    double x2 = 1.0/2.0 - x1;
    double *x3 = coord;
    double x4 = x3[19];
    double x5 = x1 - 1.0/2.0;
    double x6 = x3[13];
    double x7 = -s;
    double x8 = r*s;
    double x9 = x7 + x8;
    double x10 = x3[22];
    double x11 = -x8;
    double x12 = x11 + x7;
    double x13 = x3[16];
    double x14 = (1.0/4.0)*x0;
    double x15 = -x14;
    double x16 = (1.0/2.0)*x8;
    double x17 = (1.0/2.0)*s;
    double x18 = (1.0/4.0)*r;
    double x19 = x17 - x18;
    double x20 = x15 + x16 + x19;
    double x21 = x3[4];
    double x22 = -x16;
    double x23 = x14 + x19 + x22;
    double x24 = x3[10];
    double x25 = x17 + x18;
    double x26 = x15 + x22 + x25;
    double x27 = x3[1];
    double x28 = x14 + x16 + x25;
    double x29 = x3[7];
    double x30 = x10*x9 + x12*x13 + x2*x4 + x20*x21 + x23*x24 + x26*x27 + x28*x29 + x5*x6;
    double x31 = s*s;
    double x32 = (1.0/2.0)*x31;
    double x33 = 1.0/2.0 - x32;
    double x34 = x3[15];
    double x35 = x32 - 1.0/2.0;
    double x36 = x3[21];
    double x37 = -r;
    double x38 = x37 + x8;
    double x39 = x3[12];
    double x40 = x11 + x37;
    double x41 = x3[18];
    double x42 = (1.0/4.0)*x31;
    double x43 = -x42;
    double x44 = (1.0/2.0)*r;
    double x45 = (1.0/4.0)*s;
    double x46 = x44 - x45;
    double x47 = x16 + x43 + x46;
    double x48 = x3[9];
    double x49 = x22 + x42 + x46;
    double x50 = x3[3];
    double x51 = x44 + x45;
    double x52 = x22 + x43 + x51;
    double x53 = x3[0];
    double x54 = x16 + x42 + x51;
    double x55 = x3[6];
    double x56 = x33*x34 + x35*x36 + x38*x39 + x40*x41 + x47*x48 + x49*x50 + x52*x53 + x54*x55;
    double x57 = x2*x41;
    double x58 = x39*x5;
    double x59 = x36*x9;
    double x60 = x12*x34;
    double x61 = x20*x50;
    double x62 = x23*x48;
    double x63 = x26*x53;
    double x64 = x28*x55;
    double x65 = x13*x33;
    double x66 = x10*x35;
    double x67 = x38*x6;
    double x68 = x4*x40;
    double x69 = x24*x47;
    double x70 = x21*x49;
    double x71 = x27*x52;
    double x72 = x29*x54;
    double x73 = x30*x56 - (x57 + x58 + x59 + x60 + x61 + x62 + x63 + x64)*(x65 + x66 + x67 + x68 + x69 + x70 + x71 + x72);
    double x74 = 1.0/x73;
    double x75 = x30*x74;
    double x76 = x74*(-x65 - x66 - x67 - x68 - x69 - x70 - x71 - x72);
    double x77 = x26*x76 + x52*x75;
    double *x78 = U;
    double x79 = x78[19];
    double x80 = x78[13];
    double x81 = x78[22];
    double x82 = x78[16];
    double x83 = x78[4];
    double x84 = x78[10];
    double x85 = x78[1];
    double x86 = x78[7];
    double x87 = x12*x82 + x2*x79 + x20*x83 + x23*x84 + x26*x85 + x28*x86 + x5*x80 + x81*x9;
    double x88 = x56*x74;
    double x89 = x87*x88;
    double x90 = x33*x82 + x35*x81 + x38*x80 + x40*x79 + x47*x84 + x49*x83 + x52*x85 + x54*x86;
    double x91 = x74*(-x57 - x58 - x59 - x60 - x61 - x62 - x63 - x64);
    double x92 = x90*x91;
    double x93 = 1.0*x89 + 1.0*x92 + 1.0;
    double x94 = x78[18];
    double x95 = x78[12];
    double x96 = x78[21];
    double x97 = x78[15];
    double x98 = x78[3];
    double x99 = x78[9];
    double x100 = x78[0];
    double x101 = x78[6];
    double x102 = x100*x26 + x101*x28 + x12*x97 + x2*x94 + x20*x98 + x23*x99 + x5*x95 + x9*x96;
    double x103 = x102*x88;
    double x104 = x100*x52 + x101*x54 + x33*x97 + x35*x96 + x38*x95 + x40*x94 + x47*x99 + x49*x98;
    double x105 = x104*x91;
    double x106 = -1.0*x103 - 1.0*x105;
    double x107 = S11*x93 + S12*x106;
    double x108 = x102*x76;
    double x109 = x104*x75;
    double x110 = x108 + x109;
    double x111 = S12*x93 + S22*x106;
    double x112 = x76*x87;
    double x113 = x75*x90;
    double x114 = x112 + x113;
    double x115 = S13*x93 + S23*x106;
    double x116 = x78[20];
    double x117 = x78[14];
    double x118 = x78[23];
    double x119 = x78[17];
    double x120 = x78[5];
    double x121 = x78[11];
    double x122 = x78[2];
    double x123 = x78[8];
    double x124 = x116*x2 + x117*x5 + x118*x9 + x119*x12 + x120*x20 + x121*x23 + x122*x26 + x123*x28;
    double x125 = x116*x40 + x117*x38 + x118*x35 + x119*x33 + x120*x49 + x121*x47 + x122*x52 + x123*x54;
    double x126 = x124*x76 + x125*x75;
    double x127 = r*x32;
    double x128 = -x127 + x33 + x44;
    double *x129 = V;
    double x130 = x129[15];
    double x131 = s*x1;
    double x132 = -x131 + x17 + x2;
    double x133 = x129[18];
    double x134 = x127 + x33 - x44;
    double x135 = x129[21];
    double x136 = x131 - x17 + x2;
    double x137 = x129[12];
    double x138 = x18*x31;
    double x139 = -x138;
    double x140 = s*x14;
    double x141 = (1.0/4.0)*x8;
    double x142 = x14 - x141 + x42 - 1.0/4.0;
    double x143 = x139 + x140 + x142;
    double x144 = x129[9];
    double x145 = -x140;
    double x146 = x138 + x142 + x145;
    double x147 = x129[3];
    double x148 = x14 + x141 + x42 - 1.0/4.0;
    double x149 = x139 + x145 + x148;
    double x150 = x129[0];
    double x151 = x138 + x140 + x148;
    double x152 = x129[6];
    double x153 = x128*x130 + x132*x133 + x134*x135 + x136*x137 + x143*x144 + x146*x147 + x149*x150 + x151*x152;
    double x154 = x129[16];
    double x155 = x129[19];
    double x156 = x129[22];
    double x157 = x129[13];
    double x158 = x129[10];
    double x159 = x129[4];
    double x160 = x129[1];
    double x161 = x129[7];
    double x162 = x128*x154 + x132*x155 + x134*x156 + x136*x157 + x143*x158 + x146*x159 + x149*x160 + x151*x161;
    double x163 = x129[17];
    double x164 = x129[20];
    double x165 = x129[23];
    double x166 = x129[14];
    double x167 = x129[11];
    double x168 = x129[5];
    double x169 = x129[2];
    double x170 = x129[8];
    double x171 = x128*x163 + x132*x164 + x134*x165 + x136*x166 + x143*x167 + x146*x168 + x149*x169 + x151*x170;
    double x172 = 1.0*PENER + 1.0*SENER - 0.5*rho*(x153*x153 + x162*x162 + x171*x171);
    double x173 = -x107*x110 - x111*x114 - x115*x126 + x172;
    double x174 = x26*x88 + x52*x91;
    double x175 = -1.0*x112 - 1.0*x113;
    double x176 = 1.0*x108 + 1.0*x109 + 1.0;
    double x177 = S11*x175 + S12*x176;
    double x178 = S12*x175 + S22*x176;
    double x179 = S13*x175 + S23*x176;
    double x180 = -x110*x177 - x114*x178 - x126*x179;
    double *x181 = A;
    double x182 = x128*x181[15] + x132*x181[18] + x134*x181[21] + x136*x181[12] + x143*x181[9] + x146*x181[3] + x149*x181[0] + x151*x181[6];
    double x183 = x128*x181[16] + x132*x181[19] + x134*x181[22] + x136*x181[13] + x143*x181[10] + x146*x181[4] + x149*x181[1] + x151*x181[7];
    double x184 = x128*x181[17] + x132*x181[20] + x134*x181[23] + x136*x181[14] + x143*x181[11] + x146*x181[5] + x149*x181[2] + x151*x181[8];
    double x185 = x12*x130 + x133*x2 + x135*x9 + x137*x5 + x144*x23 + x147*x20 + x150*x26 + x152*x28;
    double x186 = x130*x33 + x133*x40 + x135*x35 + x137*x38 + x144*x47 + x147*x49 + x150*x52 + x152*x54;
    double x187 = x12*x154 + x155*x2 + x156*x9 + x157*x5 + x158*x23 + x159*x20 + x160*x26 + x161*x28;
    double x188 = x154*x33 + x155*x40 + x156*x35 + x157*x38 + x158*x47 + x159*x49 + x160*x52 + x161*x54;
    double x189 = x12*x163 + x164*x2 + x165*x9 + x166*x5 + x167*x23 + x168*x20 + x169*x26 + x170*x28;
    double x190 = x163*x33 + x164*x40 + x165*x35 + x166*x38 + x167*x47 + x168*x49 + x169*x52 + x170*x54;
    double x191 = x110*x182 + x114*x183 + x126*x184 + x153*(x185*x76 + x186*x75) + x162*(x187*x76 + x188*x75) + x171*(x189*x76 + x190*x75);
    double x192 = rho*x149;
    double x193 = x103 + x105;
    double x194 = x89 + x92;
    double x195 = x124*x88 + x125*x91;
    double x196 = -x107*x193 - x111*x194 - x115*x195;
    double x197 = x172 - x177*x193 - x178*x194 - x179*x195;
    double x198 = x153*(x185*x88 + x186*x91) + x162*(x187*x88 + x188*x91) + x171*(x189*x88 + x190*x91) + x182*x193 + x183*x194 + x184*x195;
    double x199 = x20*x76 + x49*x75;
    double x200 = x20*x88 + x49*x91;
    double x201 = rho*x146;
    double x202 = x28*x76 + x54*x75;
    double x203 = x28*x88 + x54*x91;
    double x204 = rho*x151;
    double x205 = x23*x76 + x47*x75;
    double x206 = x23*x88 + x47*x91;
    double x207 = rho*x143;
    double x208 = x38*x91 + x5*x88;
    double x209 = x38*x75 + x5*x76;
    double x210 = rho*x136;
    double x211 = x12*x76 + x33*x75;
    double x212 = x12*x88 + x33*x91;
    double x213 = rho*x128;
    double x214 = x2*x88 + x40*x91;
    double x215 = x2*x76 + x40*x75;
    double x216 = rho*x132;
    double x217 = x35*x75 + x76*x9;
    double x218 = x35*x91 + x88*x9;
    double x219 = rho*x134;
    
    res_0[0] = x73*(x173*x77 + x174*x180 - x191*x192);
    res_0[1] = x73*(x174*x197 - x192*x198 + x196*x77);
    res_0[2] = 0;
    res_0[3] = x73*(x173*x199 + x180*x200 - x191*x201);
    res_0[4] = x73*(x196*x199 + x197*x200 - x198*x201);
    res_0[5] = 0;
    res_0[6] = x73*(x173*x202 + x180*x203 - x191*x204);
    res_0[7] = x73*(x196*x202 + x197*x203 - x198*x204);
    res_0[8] = 0;
    res_0[9] = x73*(x173*x205 + x180*x206 - x191*x207);
    res_0[10] = x73*(x196*x205 + x197*x206 - x198*x207);
    res_0[11] = 0;
    res_0[12] = x73*(x173*x209 + x180*x208 - x191*x210);
    res_0[13] = x73*(x196*x209 + x197*x208 - x198*x210);
    res_0[14] = 0;
    res_0[15] = x73*(x173*x211 + x180*x212 - x191*x213);
    res_0[16] = x73*(x196*x211 + x197*x212 - x198*x213);
    res_0[17] = 0;
    res_0[18] = x73*(x173*x215 + x180*x214 - x191*x216);
    res_0[19] = x73*(x196*x215 + x197*x214 - x198*x216);
    res_0[20] = 0;
    res_0[21] = x73*(x173*x217 + x180*x218 - x191*x219);
    res_0[22] = x73*(x196*x217 + x197*x218 - x198*x219);
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
    
    double x0 = r*s;
    double x1 = (1.0/2.0)*x0;
    double x2 = -x1;
    double x3 = s*s;
    double x4 = (1.0/4.0)*x3;
    double x5 = -x4;
    double x6 = (1.0/2.0)*r;
    double x7 = (1.0/4.0)*s;
    double x8 = x6 + x7;
    double x9 = x2 + x5 + x8;
    double x10 = r*r;
    double x11 = (1.0/2.0)*x10;
    double x12 = 1.0/2.0 - x11;
    double *x13 = coord;
    double x14 = x13[19];
    double x15 = x11 - 1.0/2.0;
    double x16 = x13[13];
    double x17 = -s;
    double x18 = x0 + x17;
    double x19 = x13[22];
    double x20 = -x0;
    double x21 = x17 + x20;
    double x22 = x13[16];
    double x23 = (1.0/4.0)*x10;
    double x24 = -x23;
    double x25 = (1.0/2.0)*s;
    double x26 = (1.0/4.0)*r;
    double x27 = x25 - x26;
    double x28 = x1 + x24 + x27;
    double x29 = x13[4];
    double x30 = x2 + x23 + x27;
    double x31 = x13[10];
    double x32 = x25 + x26;
    double x33 = x2 + x24 + x32;
    double x34 = x13[1];
    double x35 = x1 + x23 + x32;
    double x36 = x13[7];
    double x37 = x12*x14 + x15*x16 + x18*x19 + x21*x22 + x28*x29 + x30*x31 + x33*x34 + x35*x36;
    double x38 = (1.0/2.0)*x3;
    double x39 = 1.0/2.0 - x38;
    double x40 = x13[15];
    double x41 = x38 - 1.0/2.0;
    double x42 = x13[21];
    double x43 = -r;
    double x44 = x0 + x43;
    double x45 = x13[12];
    double x46 = x20 + x43;
    double x47 = x13[18];
    double x48 = x6 - x7;
    double x49 = x1 + x48 + x5;
    double x50 = x13[9];
    double x51 = x2 + x4 + x48;
    double x52 = x13[3];
    double x53 = x13[0];
    double x54 = x1 + x4 + x8;
    double x55 = x13[6];
    double x56 = x39*x40 + x41*x42 + x44*x45 + x46*x47 + x49*x50 + x51*x52 + x53*x9 + x54*x55;
    double x57 = x12*x47;
    double x58 = x15*x45;
    double x59 = x18*x42;
    double x60 = x21*x40;
    double x61 = x28*x52;
    double x62 = x30*x50;
    double x63 = x33*x53;
    double x64 = x35*x55;
    double x65 = x22*x39;
    double x66 = x19*x41;
    double x67 = x16*x44;
    double x68 = x14*x46;
    double x69 = x31*x49;
    double x70 = x29*x51;
    double x71 = x34*x9;
    double x72 = x36*x54;
    double x73 = x37*x56 - (x57 + x58 + x59 + x60 + x61 + x62 + x63 + x64)*(x65 + x66 + x67 + x68 + x69 + x70 + x71 + x72);
    double x74 = 1.0/x73;
    double x75 = x37*x74;
    double x76 = x74*(-x65 - x66 - x67 - x68 - x69 - x70 - x71 - x72);
    double x77 = x33*x76 + x75*x9;
    double *x78 = U;
    double x79 = x78[19];
    double x80 = x78[13];
    double x81 = x78[22];
    double x82 = x78[16];
    double x83 = x78[4];
    double x84 = x78[10];
    double x85 = x78[1];
    double x86 = x78[7];
    double x87 = x12*x79 + x15*x80 + x18*x81 + x21*x82 + x28*x83 + x30*x84 + x33*x85 + x35*x86;
    double x88 = x56*x74;
    double x89 = x87*x88;
    double x90 = x39*x82 + x41*x81 + x44*x80 + x46*x79 + x49*x84 + x51*x83 + x54*x86 + x85*x9;
    double x91 = x74*(-x57 - x58 - x59 - x60 - x61 - x62 - x63 - x64);
    double x92 = x90*x91;
    double x93 = 1.0*x89 + 1.0*x92 + 1.0;
    double x94 = x78[18];
    double x95 = x78[12];
    double x96 = x78[21];
    double x97 = x78[15];
    double x98 = x78[3];
    double x99 = x78[9];
    double x100 = x78[0];
    double x101 = x78[6];
    double x102 = x100*x33 + x101*x35 + x12*x94 + x15*x95 + x18*x96 + x21*x97 + x28*x98 + x30*x99;
    double x103 = x102*x88;
    double x104 = x100*x9 + x101*x54 + x39*x97 + x41*x96 + x44*x95 + x46*x94 + x49*x99 + x51*x98;
    double x105 = x104*x91;
    double x106 = -1.0*x103 - 1.0*x105;
    double x107 = S11*x93 + S12*x106;
    double x108 = x102*x76;
    double x109 = x104*x75;
    double x110 = x108 + x109 + 1.0;
    double x111 = S12*x93 + S22*x106;
    double x112 = x76*x87;
    double x113 = x75*x90;
    double x114 = x112 + x113;
    double x115 = S13*x93;
    double x116 = S23*x106;
    double x117 = x115 + x116;
    double x118 = x78[20];
    double x119 = x78[14];
    double x120 = x78[23];
    double x121 = x78[17];
    double x122 = x78[5];
    double x123 = x78[11];
    double x124 = x78[2];
    double x125 = x78[8];
    double x126 = x118*x12 + x119*x15 + x120*x18 + x121*x21 + x122*x28 + x123*x30 + x124*x33 + x125*x35;
    double x127 = x118*x46 + x119*x44 + x120*x41 + x121*x39 + x122*x51 + x123*x49 + x124*x9 + x125*x54;
    double x128 = x126*x76 + x127*x75;
    double x129 = 1.0*PENER + 1.0*SENER;
    double x130 = -x107*x110 - x111*x114 - x117*x128 + x129;
    double x131 = x33*x88 + x9*x91;
    double x132 = -1.0*x112 - 1.0*x113;
    double x133 = 1.0*x108 + 1.0*x109 + 1.0;
    double x134 = S11*x132 + S12*x133;
    double x135 = S12*x132 + S22*x133;
    double x136 = S13*x132;
    double x137 = S23*x133;
    double x138 = x136 + x137;
    double x139 = -x110*x134 - x114*x135 - x128*x138;
    double x140 = x103 + x105;
    double x141 = x89 + x92 + 1.0;
    double x142 = x126*x88 + x127*x91;
    double x143 = -x107*x140 - x111*x141 - x117*x142;
    double x144 = x129 - x134*x140 - x135*x141 - x138*x142;
    double x145 = -1.0*x136 - 1.0*x137;
    double x146 = -1.0*x115 - 1.0*x116;
    double x147 = x28*x76 + x51*x75;
    double x148 = x28*x88 + x51*x91;
    double x149 = x35*x76 + x54*x75;
    double x150 = x35*x88 + x54*x91;
    double x151 = x30*x76 + x49*x75;
    double x152 = x30*x88 + x49*x91;
    double x153 = x15*x88 + x44*x91;
    double x154 = x15*x76 + x44*x75;
    double x155 = x21*x76 + x39*x75;
    double x156 = x21*x88 + x39*x91;
    double x157 = x12*x88 + x46*x91;
    double x158 = x12*x76 + x46*x75;
    double x159 = x18*x76 + x41*x75;
    double x160 = x18*x88 + x41*x91;
    
    res_0[0] = x73*(x130*x77 + x131*x139);
    res_0[1] = x73*(x131*x144 + x143*x77);
    res_0[2] = x73*(x131*x145 + x146*x77);
    res_0[3] = x73*(x130*x147 + x139*x148);
    res_0[4] = x73*(x143*x147 + x144*x148);
    res_0[5] = x73*(x145*x148 + x146*x147);
    res_0[6] = x73*(x130*x149 + x139*x150);
    res_0[7] = x73*(x143*x149 + x144*x150);
    res_0[8] = x73*(x145*x150 + x146*x149);
    res_0[9] = x73*(x130*x151 + x139*x152);
    res_0[10] = x73*(x143*x151 + x144*x152);
    res_0[11] = x73*(x145*x152 + x146*x151);
    res_0[12] = x73*(x130*x154 + x139*x153);
    res_0[13] = x73*(x143*x154 + x144*x153);
    res_0[14] = x73*(x145*x153 + x146*x154);
    res_0[15] = x73*(x130*x155 + x139*x156);
    res_0[16] = x73*(x143*x155 + x144*x156);
    res_0[17] = x73*(x145*x156 + x146*x155);
    res_0[18] = x73*(x130*x158 + x139*x157);
    res_0[19] = x73*(x143*x158 + x144*x157);
    res_0[20] = x73*(x145*x157 + x146*x158);
    res_0[21] = x73*(x130*x159 + x139*x160);
    res_0[22] = x73*(x143*x159 + x144*x160);
    res_0[23] = x73*(x145*x160 + x146*x159);
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
    
    double x0 = r*s;
    double x1 = (1.0/2.0)*x0;
    double x2 = -x1;
    double x3 = s*s;
    double x4 = (1.0/4.0)*x3;
    double x5 = -x4;
    double x6 = (1.0/2.0)*r;
    double x7 = (1.0/4.0)*s;
    double x8 = x6 + x7;
    double x9 = x2 + x5 + x8;
    double x10 = r*r;
    double x11 = (1.0/2.0)*x10;
    double x12 = 1.0/2.0 - x11;
    double *x13 = coord;
    double x14 = x13[19];
    double x15 = x11 - 1.0/2.0;
    double x16 = x13[13];
    double x17 = -s;
    double x18 = x0 + x17;
    double x19 = x13[22];
    double x20 = -x0;
    double x21 = x17 + x20;
    double x22 = x13[16];
    double x23 = (1.0/4.0)*x10;
    double x24 = -x23;
    double x25 = (1.0/2.0)*s;
    double x26 = (1.0/4.0)*r;
    double x27 = x25 - x26;
    double x28 = x1 + x24 + x27;
    double x29 = x13[4];
    double x30 = x2 + x23 + x27;
    double x31 = x13[10];
    double x32 = x25 + x26;
    double x33 = x2 + x24 + x32;
    double x34 = x13[1];
    double x35 = x1 + x23 + x32;
    double x36 = x13[7];
    double x37 = x12*x14 + x15*x16 + x18*x19 + x21*x22 + x28*x29 + x30*x31 + x33*x34 + x35*x36;
    double x38 = (1.0/2.0)*x3;
    double x39 = 1.0/2.0 - x38;
    double x40 = x13[15];
    double x41 = x38 - 1.0/2.0;
    double x42 = x13[21];
    double x43 = -r;
    double x44 = x0 + x43;
    double x45 = x13[12];
    double x46 = x20 + x43;
    double x47 = x13[18];
    double x48 = x6 - x7;
    double x49 = x1 + x48 + x5;
    double x50 = x13[9];
    double x51 = x2 + x4 + x48;
    double x52 = x13[3];
    double x53 = x13[0];
    double x54 = x1 + x4 + x8;
    double x55 = x13[6];
    double x56 = x39*x40 + x41*x42 + x44*x45 + x46*x47 + x49*x50 + x51*x52 + x53*x9 + x54*x55;
    double x57 = x12*x47;
    double x58 = x15*x45;
    double x59 = x18*x42;
    double x60 = x21*x40;
    double x61 = x28*x52;
    double x62 = x30*x50;
    double x63 = x33*x53;
    double x64 = x35*x55;
    double x65 = x22*x39;
    double x66 = x19*x41;
    double x67 = x16*x44;
    double x68 = x14*x46;
    double x69 = x31*x49;
    double x70 = x29*x51;
    double x71 = x34*x9;
    double x72 = x36*x54;
    double x73 = x37*x56 - (x57 + x58 + x59 + x60 + x61 + x62 + x63 + x64)*(x65 + x66 + x67 + x68 + x69 + x70 + x71 + x72);
    double x74 = 1.0/x73;
    double x75 = x37*x74;
    double x76 = x74*(-x65 - x66 - x67 - x68 - x69 - x70 - x71 - x72);
    double x77 = x33*x76 + x75*x9;
    double *x78 = U;
    double x79 = x78[19];
    double x80 = x78[13];
    double x81 = x78[22];
    double x82 = x78[16];
    double x83 = x78[4];
    double x84 = x78[10];
    double x85 = x78[1];
    double x86 = x78[7];
    double x87 = x12*x79 + x15*x80 + x18*x81 + x21*x82 + x28*x83 + x30*x84 + x33*x85 + x35*x86;
    double x88 = x56*x74;
    double x89 = x87*x88;
    double x90 = x39*x82 + x41*x81 + x44*x80 + x46*x79 + x49*x84 + x51*x83 + x54*x86 + x85*x9;
    double x91 = x74*(-x57 - x58 - x59 - x60 - x61 - x62 - x63 - x64);
    double x92 = x90*x91;
    double x93 = 1.0*x89 + 1.0*x92 + 1.0;
    double x94 = x78[18];
    double x95 = x78[12];
    double x96 = x78[21];
    double x97 = x78[15];
    double x98 = x78[3];
    double x99 = x78[9];
    double x100 = x78[0];
    double x101 = x78[6];
    double x102 = x100*x33 + x101*x35 + x12*x94 + x15*x95 + x18*x96 + x21*x97 + x28*x98 + x30*x99;
    double x103 = x102*x88;
    double x104 = x100*x9 + x101*x54 + x39*x97 + x41*x96 + x44*x95 + x46*x94 + x49*x99 + x51*x98;
    double x105 = x104*x91;
    double x106 = -1.0*x103 - 1.0*x105;
    double x107 = S11*x93 + S12*x106;
    double x108 = x102*x76;
    double x109 = x104*x75;
    double x110 = x108 + x109;
    double x111 = S12*x93 + S22*x106;
    double x112 = x76*x87;
    double x113 = x75*x90;
    double x114 = x112 + x113;
    double x115 = S13*x93 + S23*x106;
    double x116 = x78[20];
    double x117 = x78[14];
    double x118 = x78[23];
    double x119 = x78[17];
    double x120 = x78[5];
    double x121 = x78[11];
    double x122 = x78[2];
    double x123 = x78[8];
    double x124 = x116*x12 + x117*x15 + x118*x18 + x119*x21 + x120*x28 + x121*x30 + x122*x33 + x123*x35;
    double x125 = x116*x46 + x117*x44 + x118*x41 + x119*x39 + x120*x51 + x121*x49 + x122*x9 + x123*x54;
    double x126 = x124*x76 + x125*x75;
    double x127 = 1.0*PENER + 1.0*SENER;
    double x128 = -x107*x110 - x111*x114 - x115*x126 + x127;
    double x129 = x33*x88 + x9*x91;
    double x130 = -1.0*x112 - 1.0*x113;
    double x131 = 1.0*x108 + 1.0*x109 + 1.0;
    double x132 = S11*x130 + S12*x131;
    double x133 = S12*x130 + S22*x131;
    double x134 = S13*x130 + S23*x131;
    double x135 = -x110*x132 - x114*x133 - x126*x134;
    double x136 = x103 + x105;
    double x137 = x89 + x92;
    double x138 = x124*x88 + x125*x91;
    double x139 = -x107*x136 - x111*x137 - x115*x138;
    double x140 = x127 - x132*x136 - x133*x137 - x134*x138;
    double x141 = x28*x76 + x51*x75;
    double x142 = x28*x88 + x51*x91;
    double x143 = x35*x76 + x54*x75;
    double x144 = x35*x88 + x54*x91;
    double x145 = x30*x76 + x49*x75;
    double x146 = x30*x88 + x49*x91;
    double x147 = x15*x88 + x44*x91;
    double x148 = x15*x76 + x44*x75;
    double x149 = x21*x76 + x39*x75;
    double x150 = x21*x88 + x39*x91;
    double x151 = x12*x88 + x46*x91;
    double x152 = x12*x76 + x46*x75;
    double x153 = x18*x76 + x41*x75;
    double x154 = x18*x88 + x41*x91;
    
    res_0[0] = x73*(x128*x77 + x129*x135);
    res_0[1] = x73*(x129*x140 + x139*x77);
    res_0[2] = 0;
    res_0[3] = x73*(x128*x141 + x135*x142);
    res_0[4] = x73*(x139*x141 + x140*x142);
    res_0[5] = 0;
    res_0[6] = x73*(x128*x143 + x135*x144);
    res_0[7] = x73*(x139*x143 + x140*x144);
    res_0[8] = 0;
    res_0[9] = x73*(x128*x145 + x135*x146);
    res_0[10] = x73*(x139*x145 + x140*x146);
    res_0[11] = 0;
    res_0[12] = x73*(x128*x148 + x135*x147);
    res_0[13] = x73*(x139*x148 + x140*x147);
    res_0[14] = 0;
    res_0[15] = x73*(x128*x149 + x135*x150);
    res_0[16] = x73*(x139*x149 + x140*x150);
    res_0[17] = 0;
    res_0[18] = x73*(x128*x152 + x135*x151);
    res_0[19] = x73*(x139*x152 + x140*x151);
    res_0[20] = 0;
    res_0[21] = x73*(x128*x153 + x135*x154);
    res_0[22] = x73*(x139*x153 + x140*x154);
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
    double x2 = 1.0/2.0 - x1;
    double *x3 = coord;
    double x4 = x3[19];
    double x5 = x1 - 1.0/2.0;
    double x6 = x3[13];
    double x7 = -s;
    double x8 = r*s;
    double x9 = x7 + x8;
    double x10 = x3[22];
    double x11 = -x8;
    double x12 = x11 + x7;
    double x13 = x3[16];
    double x14 = (1.0/4.0)*x0;
    double x15 = -x14;
    double x16 = (1.0/2.0)*x8;
    double x17 = (1.0/2.0)*s;
    double x18 = (1.0/4.0)*r;
    double x19 = x17 - x18;
    double x20 = x15 + x16 + x19;
    double x21 = x3[4];
    double x22 = -x16;
    double x23 = x14 + x19 + x22;
    double x24 = x3[10];
    double x25 = x17 + x18;
    double x26 = x15 + x22 + x25;
    double x27 = x3[1];
    double x28 = x14 + x16 + x25;
    double x29 = x3[7];
    double x30 = x10*x9 + x12*x13 + x2*x4 + x20*x21 + x23*x24 + x26*x27 + x28*x29 + x5*x6;
    double x31 = s*s;
    double x32 = (1.0/2.0)*x31;
    double x33 = 1.0/2.0 - x32;
    double x34 = x3[15];
    double x35 = x32 - 1.0/2.0;
    double x36 = x3[21];
    double x37 = -r;
    double x38 = x37 + x8;
    double x39 = x3[12];
    double x40 = x11 + x37;
    double x41 = x3[18];
    double x42 = (1.0/4.0)*x31;
    double x43 = -x42;
    double x44 = (1.0/2.0)*r;
    double x45 = (1.0/4.0)*s;
    double x46 = x44 - x45;
    double x47 = x16 + x43 + x46;
    double x48 = x3[9];
    double x49 = x22 + x42 + x46;
    double x50 = x3[3];
    double x51 = x44 + x45;
    double x52 = x22 + x43 + x51;
    double x53 = x3[0];
    double x54 = x16 + x42 + x51;
    double x55 = x3[6];
    double x56 = x33*x34 + x35*x36 + x38*x39 + x40*x41 + x47*x48 + x49*x50 + x52*x53 + x54*x55;
    double x57 = x2*x41;
    double x58 = x39*x5;
    double x59 = x36*x9;
    double x60 = x12*x34;
    double x61 = x20*x50;
    double x62 = x23*x48;
    double x63 = x26*x53;
    double x64 = x28*x55;
    double x65 = x13*x33;
    double x66 = x10*x35;
    double x67 = x38*x6;
    double x68 = x4*x40;
    double x69 = x24*x47;
    double x70 = x21*x49;
    double x71 = x27*x52;
    double x72 = x29*x54;
    double x73 = x30*x56 - (x57 + x58 + x59 + x60 + x61 + x62 + x63 + x64)*(x65 + x66 + x67 + x68 + x69 + x70 + x71 + x72);
    double x74 = 1.0/x73;
    double x75 = x30*x74;
    double x76 = x74*(-x65 - x66 - x67 - x68 - x69 - x70 - x71 - x72);
    double x77 = x26*x76 + x52*x75;
    double *x78 = U;
    double x79 = x78[19];
    double x80 = x78[13];
    double x81 = x78[22];
    double x82 = x78[16];
    double x83 = x78[4];
    double x84 = x78[10];
    double x85 = x78[1];
    double x86 = x78[7];
    double x87 = x12*x82 + x2*x79 + x20*x83 + x23*x84 + x26*x85 + x28*x86 + x5*x80 + x81*x9;
    double x88 = x56*x74;
    double x89 = x87*x88;
    double x90 = x33*x82 + x35*x81 + x38*x80 + x40*x79 + x47*x84 + x49*x83 + x52*x85 + x54*x86;
    double x91 = x74*(-x57 - x58 - x59 - x60 - x61 - x62 - x63 - x64);
    double x92 = x90*x91;
    double x93 = 1.0*x89 + 1.0*x92 + 1.0;
    double x94 = x78[18];
    double x95 = x78[12];
    double x96 = x78[21];
    double x97 = x78[15];
    double x98 = x78[3];
    double x99 = x78[9];
    double x100 = x78[0];
    double x101 = x78[6];
    double x102 = x100*x26 + x101*x28 + x12*x97 + x2*x94 + x20*x98 + x23*x99 + x5*x95 + x9*x96;
    double x103 = x102*x88;
    double x104 = x100*x52 + x101*x54 + x33*x97 + x35*x96 + x38*x95 + x40*x94 + x47*x99 + x49*x98;
    double x105 = x104*x91;
    double x106 = -1.0*x103 - 1.0*x105;
    double x107 = S11*x93 + S12*x106;
    double x108 = x102*x76;
    double x109 = x104*x75;
    double x110 = x108 + x109;
    double x111 = S12*x93 + S22*x106;
    double x112 = x76*x87;
    double x113 = x75*x90;
    double x114 = x112 + x113;
    double x115 = S13*x93 + S23*x106;
    double x116 = x78[20];
    double x117 = x78[14];
    double x118 = x78[23];
    double x119 = x78[17];
    double x120 = x78[5];
    double x121 = x78[11];
    double x122 = x78[2];
    double x123 = x78[8];
    double x124 = x116*x2 + x117*x5 + x118*x9 + x119*x12 + x120*x20 + x121*x23 + x122*x26 + x123*x28;
    double x125 = x116*x40 + x117*x38 + x118*x35 + x119*x33 + x120*x49 + x121*x47 + x122*x52 + x123*x54;
    double x126 = x124*x76 + x125*x75;
    double x127 = r*x32;
    double x128 = -x127 + x33 + x44;
    double *x129 = V;
    double x130 = x129[15];
    double x131 = s*x1;
    double x132 = -x131 + x17 + x2;
    double x133 = x129[18];
    double x134 = x127 + x33 - x44;
    double x135 = x129[21];
    double x136 = x131 - x17 + x2;
    double x137 = x129[12];
    double x138 = x18*x31;
    double x139 = -x138;
    double x140 = s*x14;
    double x141 = (1.0/4.0)*x8;
    double x142 = x14 - x141 + x42 - 1.0/4.0;
    double x143 = x139 + x140 + x142;
    double x144 = x129[9];
    double x145 = -x140;
    double x146 = x138 + x142 + x145;
    double x147 = x129[3];
    double x148 = x14 + x141 + x42 - 1.0/4.0;
    double x149 = x139 + x145 + x148;
    double x150 = x129[0];
    double x151 = x138 + x140 + x148;
    double x152 = x129[6];
    double x153 = x128*x130 + x132*x133 + x134*x135 + x136*x137 + x143*x144 + x146*x147 + x149*x150 + x151*x152;
    double x154 = x129[16];
    double x155 = x129[19];
    double x156 = x129[22];
    double x157 = x129[13];
    double x158 = x129[10];
    double x159 = x129[4];
    double x160 = x129[1];
    double x161 = x129[7];
    double x162 = x128*x154 + x132*x155 + x134*x156 + x136*x157 + x143*x158 + x146*x159 + x149*x160 + x151*x161;
    double x163 = x129[17];
    double x164 = x129[20];
    double x165 = x129[23];
    double x166 = x129[14];
    double x167 = x129[11];
    double x168 = x129[5];
    double x169 = x129[2];
    double x170 = x129[8];
    double x171 = x128*x163 + x132*x164 + x134*x165 + x136*x166 + x143*x167 + x146*x168 + x149*x169 + x151*x170;
    double x172 = 1.0*PENER + 1.0*SENER - 0.5*rho*(x153*x153 + x162*x162 + x171*x171);
    double x173 = -x107*x110 - x111*x114 - x115*x126 + x172;
    double x174 = x26*x88 + x52*x91;
    double x175 = -1.0*x112 - 1.0*x113;
    double x176 = 1.0*x108 + 1.0*x109 + 1.0;
    double x177 = S11*x175 + S12*x176;
    double x178 = S12*x175 + S22*x176;
    double x179 = S13*x175 + S23*x176;
    double x180 = -x110*x177 - x114*x178 - x126*x179;
    double *x181 = A;
    double x182 = x128*x181[15] + x132*x181[18] + x134*x181[21] + x136*x181[12] + x143*x181[9] + x146*x181[3] + x149*x181[0] + x151*x181[6];
    double x183 = x128*x181[16] + x132*x181[19] + x134*x181[22] + x136*x181[13] + x143*x181[10] + x146*x181[4] + x149*x181[1] + x151*x181[7];
    double x184 = x128*x181[17] + x132*x181[20] + x134*x181[23] + x136*x181[14] + x143*x181[11] + x146*x181[5] + x149*x181[2] + x151*x181[8];
    double x185 = x12*x130 + x133*x2 + x135*x9 + x137*x5 + x144*x23 + x147*x20 + x150*x26 + x152*x28;
    double x186 = x130*x33 + x133*x40 + x135*x35 + x137*x38 + x144*x47 + x147*x49 + x150*x52 + x152*x54;
    double x187 = x12*x154 + x155*x2 + x156*x9 + x157*x5 + x158*x23 + x159*x20 + x160*x26 + x161*x28;
    double x188 = x154*x33 + x155*x40 + x156*x35 + x157*x38 + x158*x47 + x159*x49 + x160*x52 + x161*x54;
    double x189 = x12*x163 + x164*x2 + x165*x9 + x166*x5 + x167*x23 + x168*x20 + x169*x26 + x170*x28;
    double x190 = x163*x33 + x164*x40 + x165*x35 + x166*x38 + x167*x47 + x168*x49 + x169*x52 + x170*x54;
    double x191 = x110*x182 + x114*x183 + x126*x184 + x153*(x185*x76 + x186*x75) + x162*(x187*x76 + x188*x75) + x171*(x189*x76 + x190*x75);
    double x192 = rho*x149;
    double x193 = x103 + x105;
    double x194 = x89 + x92;
    double x195 = x124*x88 + x125*x91;
    double x196 = -x107*x193 - x111*x194 - x115*x195;
    double x197 = x172 - x177*x193 - x178*x194 - x179*x195;
    double x198 = x153*(x185*x88 + x186*x91) + x162*(x187*x88 + x188*x91) + x171*(x189*x88 + x190*x91) + x182*x193 + x183*x194 + x184*x195;
    double x199 = x20*x76 + x49*x75;
    double x200 = x20*x88 + x49*x91;
    double x201 = rho*x146;
    double x202 = x28*x76 + x54*x75;
    double x203 = x28*x88 + x54*x91;
    double x204 = rho*x151;
    double x205 = x23*x76 + x47*x75;
    double x206 = x23*x88 + x47*x91;
    double x207 = rho*x143;
    double x208 = x38*x91 + x5*x88;
    double x209 = x38*x75 + x5*x76;
    double x210 = rho*x136;
    double x211 = x12*x76 + x33*x75;
    double x212 = x12*x88 + x33*x91;
    double x213 = rho*x128;
    double x214 = x2*x88 + x40*x91;
    double x215 = x2*x76 + x40*x75;
    double x216 = rho*x132;
    double x217 = x35*x75 + x76*x9;
    double x218 = x35*x91 + x88*x9;
    double x219 = rho*x134;
    
    res_0[0] = x73*(x173*x77 + x174*x180 - x191*x192);
    res_0[1] = x73*(x174*x197 - x192*x198 + x196*x77);
    res_0[2] = 0;
    res_0[3] = x73*(x173*x199 + x180*x200 - x191*x201);
    res_0[4] = x73*(x196*x199 + x197*x200 - x198*x201);
    res_0[5] = 0;
    res_0[6] = x73*(x173*x202 + x180*x203 - x191*x204);
    res_0[7] = x73*(x196*x202 + x197*x203 - x198*x204);
    res_0[8] = 0;
    res_0[9] = x73*(x173*x205 + x180*x206 - x191*x207);
    res_0[10] = x73*(x196*x205 + x197*x206 - x198*x207);
    res_0[11] = 0;
    res_0[12] = x73*(x173*x209 + x180*x208 - x191*x210);
    res_0[13] = x73*(x196*x209 + x197*x208 - x198*x210);
    res_0[14] = 0;
    res_0[15] = x73*(x173*x211 + x180*x212 - x191*x213);
    res_0[16] = x73*(x196*x211 + x197*x212 - x198*x213);
    res_0[17] = 0;
    res_0[18] = x73*(x173*x215 + x180*x214 - x191*x216);
    res_0[19] = x73*(x196*x215 + x197*x214 - x198*x216);
    res_0[20] = 0;
    res_0[21] = x73*(x173*x217 + x180*x218 - x191*x219);
    res_0[22] = x73*(x196*x217 + x197*x218 - x198*x219);
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
    
    double x0 = r*s;
    double x1 = (1.0/2.0)*x0;
    double x2 = -x1;
    double x3 = s*s;
    double x4 = (1.0/4.0)*x3;
    double x5 = -x4;
    double x6 = (1.0/2.0)*r;
    double x7 = (1.0/4.0)*s;
    double x8 = x6 + x7;
    double x9 = x2 + x5 + x8;
    double x10 = r*r;
    double x11 = (1.0/2.0)*x10;
    double x12 = 1.0/2.0 - x11;
    double *x13 = coord;
    double x14 = x13[19];
    double x15 = x11 - 1.0/2.0;
    double x16 = x13[13];
    double x17 = -s;
    double x18 = x0 + x17;
    double x19 = x13[22];
    double x20 = -x0;
    double x21 = x17 + x20;
    double x22 = x13[16];
    double x23 = (1.0/4.0)*x10;
    double x24 = -x23;
    double x25 = (1.0/2.0)*s;
    double x26 = (1.0/4.0)*r;
    double x27 = x25 - x26;
    double x28 = x1 + x24 + x27;
    double x29 = x13[4];
    double x30 = x2 + x23 + x27;
    double x31 = x13[10];
    double x32 = x25 + x26;
    double x33 = x2 + x24 + x32;
    double x34 = x13[1];
    double x35 = x1 + x23 + x32;
    double x36 = x13[7];
    double x37 = x12*x14 + x15*x16 + x18*x19 + x21*x22 + x28*x29 + x30*x31 + x33*x34 + x35*x36;
    double x38 = (1.0/2.0)*x3;
    double x39 = 1.0/2.0 - x38;
    double x40 = x13[15];
    double x41 = x38 - 1.0/2.0;
    double x42 = x13[21];
    double x43 = -r;
    double x44 = x0 + x43;
    double x45 = x13[12];
    double x46 = x20 + x43;
    double x47 = x13[18];
    double x48 = x6 - x7;
    double x49 = x1 + x48 + x5;
    double x50 = x13[9];
    double x51 = x2 + x4 + x48;
    double x52 = x13[3];
    double x53 = x13[0];
    double x54 = x1 + x4 + x8;
    double x55 = x13[6];
    double x56 = x39*x40 + x41*x42 + x44*x45 + x46*x47 + x49*x50 + x51*x52 + x53*x9 + x54*x55;
    double x57 = x12*x47;
    double x58 = x15*x45;
    double x59 = x18*x42;
    double x60 = x21*x40;
    double x61 = x28*x52;
    double x62 = x30*x50;
    double x63 = x33*x53;
    double x64 = x35*x55;
    double x65 = x22*x39;
    double x66 = x19*x41;
    double x67 = x16*x44;
    double x68 = x14*x46;
    double x69 = x31*x49;
    double x70 = x29*x51;
    double x71 = x34*x9;
    double x72 = x36*x54;
    double x73 = x37*x56 - (x57 + x58 + x59 + x60 + x61 + x62 + x63 + x64)*(x65 + x66 + x67 + x68 + x69 + x70 + x71 + x72);
    double x74 = 1.0/x73;
    double x75 = x37*x74;
    double x76 = x74*(-x65 - x66 - x67 - x68 - x69 - x70 - x71 - x72);
    double x77 = x33*x76 + x75*x9;
    double *x78 = U;
    double x79 = x78[19];
    double x80 = x78[13];
    double x81 = x78[22];
    double x82 = x78[16];
    double x83 = x78[4];
    double x84 = x78[10];
    double x85 = x78[1];
    double x86 = x78[7];
    double x87 = x12*x79 + x15*x80 + x18*x81 + x21*x82 + x28*x83 + x30*x84 + x33*x85 + x35*x86;
    double x88 = x56*x74;
    double x89 = x87*x88;
    double x90 = x39*x82 + x41*x81 + x44*x80 + x46*x79 + x49*x84 + x51*x83 + x54*x86 + x85*x9;
    double x91 = x74*(-x57 - x58 - x59 - x60 - x61 - x62 - x63 - x64);
    double x92 = x90*x91;
    double x93 = 1.0*x89 + 1.0*x92 + 1.0;
    double x94 = x78[18];
    double x95 = x78[12];
    double x96 = x78[21];
    double x97 = x78[15];
    double x98 = x78[3];
    double x99 = x78[9];
    double x100 = x78[0];
    double x101 = x78[6];
    double x102 = x100*x33 + x101*x35 + x12*x94 + x15*x95 + x18*x96 + x21*x97 + x28*x98 + x30*x99;
    double x103 = x102*x88;
    double x104 = x100*x9 + x101*x54 + x39*x97 + x41*x96 + x44*x95 + x46*x94 + x49*x99 + x51*x98;
    double x105 = x104*x91;
    double x106 = -1.0*x103 - 1.0*x105;
    double x107 = S11*x93 + S12*x106;
    double x108 = x102*x76;
    double x109 = x104*x75;
    double x110 = x108 + x109 + 1.0;
    double x111 = S12*x93 + S22*x106;
    double x112 = x76*x87;
    double x113 = x75*x90;
    double x114 = x112 + x113;
    double x115 = S13*x93;
    double x116 = S23*x106;
    double x117 = x115 + x116;
    double x118 = x78[20];
    double x119 = x78[14];
    double x120 = x78[23];
    double x121 = x78[17];
    double x122 = x78[5];
    double x123 = x78[11];
    double x124 = x78[2];
    double x125 = x78[8];
    double x126 = x118*x12 + x119*x15 + x120*x18 + x121*x21 + x122*x28 + x123*x30 + x124*x33 + x125*x35;
    double x127 = x118*x46 + x119*x44 + x120*x41 + x121*x39 + x122*x51 + x123*x49 + x124*x9 + x125*x54;
    double x128 = x126*x76 + x127*x75;
    double x129 = 1.0*PENER + 1.0*SENER;
    double x130 = -x107*x110 - x111*x114 - x117*x128 + x129;
    double x131 = x33*x88 + x9*x91;
    double x132 = -1.0*x112 - 1.0*x113;
    double x133 = 1.0*x108 + 1.0*x109 + 1.0;
    double x134 = S11*x132 + S12*x133;
    double x135 = S12*x132 + S22*x133;
    double x136 = S13*x132;
    double x137 = S23*x133;
    double x138 = x136 + x137;
    double x139 = -x110*x134 - x114*x135 - x128*x138;
    double x140 = x103 + x105;
    double x141 = x89 + x92 + 1.0;
    double x142 = x126*x88 + x127*x91;
    double x143 = -x107*x140 - x111*x141 - x117*x142;
    double x144 = x129 - x134*x140 - x135*x141 - x138*x142;
    double x145 = -1.0*x136 - 1.0*x137;
    double x146 = -1.0*x115 - 1.0*x116;
    double x147 = x28*x76 + x51*x75;
    double x148 = x28*x88 + x51*x91;
    double x149 = x35*x76 + x54*x75;
    double x150 = x35*x88 + x54*x91;
    double x151 = x30*x76 + x49*x75;
    double x152 = x30*x88 + x49*x91;
    double x153 = x15*x88 + x44*x91;
    double x154 = x15*x76 + x44*x75;
    double x155 = x21*x76 + x39*x75;
    double x156 = x21*x88 + x39*x91;
    double x157 = x12*x88 + x46*x91;
    double x158 = x12*x76 + x46*x75;
    double x159 = x18*x76 + x41*x75;
    double x160 = x18*x88 + x41*x91;
    
    res_0[0] = x73*(x130*x77 + x131*x139);
    res_0[1] = x73*(x131*x144 + x143*x77);
    res_0[2] = x73*(x131*x145 + x146*x77);
    res_0[3] = x73*(x130*x147 + x139*x148);
    res_0[4] = x73*(x143*x147 + x144*x148);
    res_0[5] = x73*(x145*x148 + x146*x147);
    res_0[6] = x73*(x130*x149 + x139*x150);
    res_0[7] = x73*(x143*x149 + x144*x150);
    res_0[8] = x73*(x145*x150 + x146*x149);
    res_0[9] = x73*(x130*x151 + x139*x152);
    res_0[10] = x73*(x143*x151 + x144*x152);
    res_0[11] = x73*(x145*x152 + x146*x151);
    res_0[12] = x73*(x130*x154 + x139*x153);
    res_0[13] = x73*(x143*x154 + x144*x153);
    res_0[14] = x73*(x145*x153 + x146*x154);
    res_0[15] = x73*(x130*x155 + x139*x156);
    res_0[16] = x73*(x143*x155 + x144*x156);
    res_0[17] = x73*(x145*x156 + x146*x155);
    res_0[18] = x73*(x130*x158 + x139*x157);
    res_0[19] = x73*(x143*x158 + x144*x157);
    res_0[20] = x73*(x145*x157 + x146*x158);
    res_0[21] = x73*(x130*x159 + x139*x160);
    res_0[22] = x73*(x143*x159 + x144*x160);
    res_0[23] = x73*(x145*x160 + x146*x159);
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
    
    double x0 = r*s;
    double x1 = (1.0/2.0)*x0;
    double x2 = -x1;
    double x3 = s*s;
    double x4 = (1.0/4.0)*x3;
    double x5 = -x4;
    double x6 = (1.0/2.0)*r;
    double x7 = (1.0/4.0)*s;
    double x8 = x6 + x7;
    double x9 = x2 + x5 + x8;
    double x10 = r*r;
    double x11 = (1.0/2.0)*x10;
    double x12 = 1.0/2.0 - x11;
    double *x13 = coord;
    double x14 = x13[19];
    double x15 = x11 - 1.0/2.0;
    double x16 = x13[13];
    double x17 = -s;
    double x18 = x0 + x17;
    double x19 = x13[22];
    double x20 = -x0;
    double x21 = x17 + x20;
    double x22 = x13[16];
    double x23 = (1.0/4.0)*x10;
    double x24 = -x23;
    double x25 = (1.0/2.0)*s;
    double x26 = (1.0/4.0)*r;
    double x27 = x25 - x26;
    double x28 = x1 + x24 + x27;
    double x29 = x13[4];
    double x30 = x2 + x23 + x27;
    double x31 = x13[10];
    double x32 = x25 + x26;
    double x33 = x2 + x24 + x32;
    double x34 = x13[1];
    double x35 = x1 + x23 + x32;
    double x36 = x13[7];
    double x37 = x12*x14 + x15*x16 + x18*x19 + x21*x22 + x28*x29 + x30*x31 + x33*x34 + x35*x36;
    double x38 = (1.0/2.0)*x3;
    double x39 = 1.0/2.0 - x38;
    double x40 = x13[15];
    double x41 = x38 - 1.0/2.0;
    double x42 = x13[21];
    double x43 = -r;
    double x44 = x0 + x43;
    double x45 = x13[12];
    double x46 = x20 + x43;
    double x47 = x13[18];
    double x48 = x6 - x7;
    double x49 = x1 + x48 + x5;
    double x50 = x13[9];
    double x51 = x2 + x4 + x48;
    double x52 = x13[3];
    double x53 = x13[0];
    double x54 = x1 + x4 + x8;
    double x55 = x13[6];
    double x56 = x39*x40 + x41*x42 + x44*x45 + x46*x47 + x49*x50 + x51*x52 + x53*x9 + x54*x55;
    double x57 = x12*x47;
    double x58 = x15*x45;
    double x59 = x18*x42;
    double x60 = x21*x40;
    double x61 = x28*x52;
    double x62 = x30*x50;
    double x63 = x33*x53;
    double x64 = x35*x55;
    double x65 = x22*x39;
    double x66 = x19*x41;
    double x67 = x16*x44;
    double x68 = x14*x46;
    double x69 = x31*x49;
    double x70 = x29*x51;
    double x71 = x34*x9;
    double x72 = x36*x54;
    double x73 = x37*x56 - (x57 + x58 + x59 + x60 + x61 + x62 + x63 + x64)*(x65 + x66 + x67 + x68 + x69 + x70 + x71 + x72);
    double x74 = 1.0/x73;
    double x75 = x37*x74;
    double x76 = x74*(-x65 - x66 - x67 - x68 - x69 - x70 - x71 - x72);
    double x77 = x33*x76 + x75*x9;
    double *x78 = U;
    double x79 = x78[19];
    double x80 = x78[13];
    double x81 = x78[22];
    double x82 = x78[16];
    double x83 = x78[4];
    double x84 = x78[10];
    double x85 = x78[1];
    double x86 = x78[7];
    double x87 = x12*x79 + x15*x80 + x18*x81 + x21*x82 + x28*x83 + x30*x84 + x33*x85 + x35*x86;
    double x88 = x56*x74;
    double x89 = x87*x88;
    double x90 = x39*x82 + x41*x81 + x44*x80 + x46*x79 + x49*x84 + x51*x83 + x54*x86 + x85*x9;
    double x91 = x74*(-x57 - x58 - x59 - x60 - x61 - x62 - x63 - x64);
    double x92 = x90*x91;
    double x93 = 1.0*x89 + 1.0*x92 + 1.0;
    double x94 = x78[18];
    double x95 = x78[12];
    double x96 = x78[21];
    double x97 = x78[15];
    double x98 = x78[3];
    double x99 = x78[9];
    double x100 = x78[0];
    double x101 = x78[6];
    double x102 = x100*x33 + x101*x35 + x12*x94 + x15*x95 + x18*x96 + x21*x97 + x28*x98 + x30*x99;
    double x103 = x102*x88;
    double x104 = x100*x9 + x101*x54 + x39*x97 + x41*x96 + x44*x95 + x46*x94 + x49*x99 + x51*x98;
    double x105 = x104*x91;
    double x106 = -1.0*x103 - 1.0*x105;
    double x107 = S11*x93 + S12*x106;
    double x108 = x102*x76;
    double x109 = x104*x75;
    double x110 = x108 + x109;
    double x111 = S12*x93 + S22*x106;
    double x112 = x76*x87;
    double x113 = x75*x90;
    double x114 = x112 + x113;
    double x115 = S13*x93 + S23*x106;
    double x116 = x78[20];
    double x117 = x78[14];
    double x118 = x78[23];
    double x119 = x78[17];
    double x120 = x78[5];
    double x121 = x78[11];
    double x122 = x78[2];
    double x123 = x78[8];
    double x124 = x116*x12 + x117*x15 + x118*x18 + x119*x21 + x120*x28 + x121*x30 + x122*x33 + x123*x35;
    double x125 = x116*x46 + x117*x44 + x118*x41 + x119*x39 + x120*x51 + x121*x49 + x122*x9 + x123*x54;
    double x126 = x124*x76 + x125*x75;
    double x127 = 1.0*PENER + 1.0*SENER;
    double x128 = -x107*x110 - x111*x114 - x115*x126 + x127;
    double x129 = x33*x88 + x9*x91;
    double x130 = -1.0*x112 - 1.0*x113;
    double x131 = 1.0*x108 + 1.0*x109 + 1.0;
    double x132 = S11*x130 + S12*x131;
    double x133 = S12*x130 + S22*x131;
    double x134 = S13*x130 + S23*x131;
    double x135 = -x110*x132 - x114*x133 - x126*x134;
    double x136 = x103 + x105;
    double x137 = x89 + x92;
    double x138 = x124*x88 + x125*x91;
    double x139 = -x107*x136 - x111*x137 - x115*x138;
    double x140 = x127 - x132*x136 - x133*x137 - x134*x138;
    double x141 = x28*x76 + x51*x75;
    double x142 = x28*x88 + x51*x91;
    double x143 = x35*x76 + x54*x75;
    double x144 = x35*x88 + x54*x91;
    double x145 = x30*x76 + x49*x75;
    double x146 = x30*x88 + x49*x91;
    double x147 = x15*x88 + x44*x91;
    double x148 = x15*x76 + x44*x75;
    double x149 = x21*x76 + x39*x75;
    double x150 = x21*x88 + x39*x91;
    double x151 = x12*x88 + x46*x91;
    double x152 = x12*x76 + x46*x75;
    double x153 = x18*x76 + x41*x75;
    double x154 = x18*x88 + x41*x91;
    
    res_0[0] = x73*(x128*x77 + x129*x135);
    res_0[1] = x73*(x129*x140 + x139*x77);
    res_0[2] = 0;
    res_0[3] = x73*(x128*x141 + x135*x142);
    res_0[4] = x73*(x139*x141 + x140*x142);
    res_0[5] = 0;
    res_0[6] = x73*(x128*x143 + x135*x144);
    res_0[7] = x73*(x139*x143 + x140*x144);
    res_0[8] = 0;
    res_0[9] = x73*(x128*x145 + x135*x146);
    res_0[10] = x73*(x139*x145 + x140*x146);
    res_0[11] = 0;
    res_0[12] = x73*(x128*x148 + x135*x147);
    res_0[13] = x73*(x139*x148 + x140*x147);
    res_0[14] = 0;
    res_0[15] = x73*(x128*x149 + x135*x150);
    res_0[16] = x73*(x139*x149 + x140*x150);
    res_0[17] = 0;
    res_0[18] = x73*(x128*x152 + x135*x151);
    res_0[19] = x73*(x139*x152 + x140*x151);
    res_0[20] = 0;
    res_0[21] = x73*(x128*x153 + x135*x154);
    res_0[22] = x73*(x139*x153 + x140*x154);
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
    
    double *x0 = coord;
    double x1 = x0[0];
    double x2 = -x1;
    double x3 = x2 + x0[3];
    double x4 = x0[1];
    double x5 = -x4;
    double x6 = x5 + x0[7];
    double x7 = x0[6];
    double x8 = x0[4];
    double x9 = x3*x6 - (x2 + x7)*(x5 + x8);
    double x10 = 1.0/x9;
    double x11 = x10*x3;
    double x12 = x10*(x1 - x7);
    double x13 = -x11 - x12;
    double *x14 = U;
    double x15 = -x14[1];
    double x16 = x15 + x14[4];
    double x17 = x10*x6;
    double x18 = x16*x17;
    double x19 = x15 + x14[7];
    double x20 = x10*(x4 - x8);
    double x21 = x19*x20;
    double x22 = -1.0*x18 - 1.0*x21;
    double x23 = -x14[0];
    double x24 = x23 + x14[3];
    double x25 = x17*x24;
    double x26 = x23 + x14[6];
    double x27 = x20*x26;
    double x28 = 1.0*x25 + 1.0*x27 + 1.0;
    double x29 = S11*x22 + S12*x28;
    double x30 = x25 + x27;
    double x31 = S12*x22 + S22*x28;
    double x32 = x18 + x21;
    double x33 = S13*x22 + S23*x28;
    double x34 = -x14[2];
    double x35 = x34 + x14[5];
    double x36 = x34 + x14[8];
    double x37 = x17*x35 + x20*x36;
    double x38 = -x29*x30 - x31*x32 - x33*x37;
    double x39 = -x17 - x20;
    double x40 = x12*x16;
    double x41 = x11*x19;
    double x42 = 1.0*x40 + 1.0*x41 + 1.0;
    double x43 = x12*x24;
    double x44 = x11*x26;
    double x45 = -1.0*x43 - 1.0*x44;
    double x46 = S11*x42 + S12*x45;
    double x47 = S12*x42 + S22*x45;
    double x48 = S13*x42 + S23*x45;
    double *x49 = V;
    double x50 = x49[3];
    double x51 = x49[6];
    double x52 = -r - s + 1;
    double x53 = x49[0];
    double x54 = r*x50 + s*x51 + x52*x53;
    double x55 = x49[4];
    double x56 = x49[7];
    double x57 = x49[1];
    double x58 = r*x55 + s*x56 + x52*x57;
    double x59 = x49[5];
    double x60 = x49[8];
    double x61 = x49[2];
    double x62 = r*x59 + s*x60 + x52*x61;
    double x63 = 1.0*PENER + 1.0*SENER - 0.5*rho*(x54*x54 + x58*x58 + x62*x62);
    double x64 = -x30*x46 - x32*x47 - x37*x48 + x63;
    double *x65 = A;
    double x66 = r*x65[3] + s*x65[6] + x52*x65[0];
    double x67 = r*x65[4] + s*x65[7] + x52*x65[1];
    double x68 = r*x65[5] + s*x65[8] + x52*x65[2];
    double x69 = -x53;
    double x70 = x50 + x69;
    double x71 = x51 + x69;
    double x72 = -x57;
    double x73 = x55 + x72;
    double x74 = x56 + x72;
    double x75 = -x61;
    double x76 = x59 + x75;
    double x77 = x60 + x75;
    double x78 = x30*x66 + x32*x67 + x37*x68 + x54*(x17*x70 + x20*x71) + x58*(x17*x73 + x20*x74) + x62*(x17*x76 + x20*x77);
    double x79 = rho*x52;
    double x80 = x43 + x44;
    double x81 = x40 + x41;
    double x82 = x11*x36 + x12*x35;
    double x83 = -x29*x80 - x31*x81 - x33*x82 + x63;
    double x84 = -x46*x80 - x47*x81 - x48*x82;
    double x85 = x54*(x11*x71 + x12*x70) + x58*(x11*x74 + x12*x73) + x62*(x11*x77 + x12*x76) + x66*x80 + x67*x81 + x68*x82;
    double x86 = r*rho;
    double x87 = rho*s;
    
    res_0[0] = x9*(x13*x38 + x39*x64 - x78*x79);
    res_0[1] = x9*(x13*x83 + x39*x84 - x79*x85);
    res_0[2] = 0;
    res_0[3] = x9*(x12*x38 + x17*x64 - x78*x86);
    res_0[4] = x9*(x12*x83 + x17*x84 - x85*x86);
    res_0[5] = 0;
    res_0[6] = x9*(x11*x38 + x20*x64 - x78*x87);
    res_0[7] = x9*(x11*x83 + x20*x84 - x85*x87);
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
    
    double *x0 = coord;
    double x1 = x0[0];
    double x2 = -x1;
    double x3 = x2 + x0[3];
    double x4 = x0[1];
    double x5 = -x4;
    double x6 = x5 + x0[7];
    double x7 = x0[6];
    double x8 = x0[4];
    double x9 = x3*x6 - (x2 + x7)*(x5 + x8);
    double x10 = 1.0/x9;
    double x11 = x10*x3;
    double x12 = x10*(x1 - x7);
    double x13 = -x11 - x12;
    double *x14 = U;
    double x15 = -x14[1];
    double x16 = x15 + x14[4];
    double x17 = x10*x6;
    double x18 = x16*x17;
    double x19 = x15 + x14[7];
    double x20 = x10*(x4 - x8);
    double x21 = x19*x20;
    double x22 = -1.0*x18 - 1.0*x21;
    double x23 = -x14[0];
    double x24 = x23 + x14[3];
    double x25 = x17*x24;
    double x26 = x23 + x14[6];
    double x27 = x20*x26;
    double x28 = 1.0*x25 + 1.0*x27 + 1.0;
    double x29 = S11*x22 + S12*x28;
    double x30 = x25 + x27 + 1.0;
    double x31 = S12*x22 + S22*x28;
    double x32 = x18 + x21;
    double x33 = S13*x22;
    double x34 = S23*x28;
    double x35 = x33 + x34;
    double x36 = -x14[2];
    double x37 = x36 + x14[5];
    double x38 = x36 + x14[8];
    double x39 = x17*x37 + x20*x38;
    double x40 = -x29*x30 - x31*x32 - x35*x39;
    double x41 = -x17 - x20;
    double x42 = x12*x16;
    double x43 = x11*x19;
    double x44 = 1.0*x42 + 1.0*x43 + 1.0;
    double x45 = x12*x24;
    double x46 = x11*x26;
    double x47 = -1.0*x45 - 1.0*x46;
    double x48 = S11*x44 + S12*x47;
    double x49 = S12*x44 + S22*x47;
    double x50 = S13*x44;
    double x51 = S23*x47;
    double x52 = x50 + x51;
    double x53 = 1.0*PENER + 1.0*SENER;
    double x54 = -x30*x48 - x32*x49 - x39*x52 + x53;
    double x55 = x45 + x46;
    double x56 = x42 + x43 + 1.0;
    double x57 = x11*x38 + x12*x37;
    double x58 = -x29*x55 - x31*x56 - x35*x57 + x53;
    double x59 = -x48*x55 - x49*x56 - x52*x57;
    double x60 = -1.0*x33 - 1.0*x34;
    double x61 = -1.0*x50 - 1.0*x51;
    
    res_0[0] = x9*(x13*x40 + x41*x54);
    res_0[1] = x9*(x13*x58 + x41*x59);
    res_0[2] = x9*(x13*x60 + x41*x61);
    res_0[3] = x9*(x12*x40 + x17*x54);
    res_0[4] = x9*(x12*x58 + x17*x59);
    res_0[5] = x9*(x12*x60 + x17*x61);
    res_0[6] = x9*(x11*x40 + x20*x54);
    res_0[7] = x9*(x11*x58 + x20*x59);
    res_0[8] = x9*(x11*x60 + x20*x61);
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
    
    double *x0 = coord;
    double x1 = x0[0];
    double x2 = -x1;
    double x3 = x2 + x0[3];
    double x4 = x0[1];
    double x5 = -x4;
    double x6 = x5 + x0[7];
    double x7 = x0[6];
    double x8 = x0[4];
    double x9 = x3*x6 - (x2 + x7)*(x5 + x8);
    double x10 = 1.0/x9;
    double x11 = x10*x3;
    double x12 = x10*(x1 - x7);
    double x13 = -x11 - x12;
    double *x14 = U;
    double x15 = -x14[1];
    double x16 = x15 + x14[4];
    double x17 = x10*x6;
    double x18 = x16*x17;
    double x19 = x15 + x14[7];
    double x20 = x10*(x4 - x8);
    double x21 = x19*x20;
    double x22 = -1.0*x18 - 1.0*x21;
    double x23 = -x14[0];
    double x24 = x23 + x14[3];
    double x25 = x17*x24;
    double x26 = x23 + x14[6];
    double x27 = x20*x26;
    double x28 = 1.0*x25 + 1.0*x27 + 1.0;
    double x29 = S11*x22 + S12*x28;
    double x30 = x25 + x27;
    double x31 = S12*x22 + S22*x28;
    double x32 = x18 + x21;
    double x33 = S13*x22 + S23*x28;
    double x34 = -x14[2];
    double x35 = x34 + x14[5];
    double x36 = x34 + x14[8];
    double x37 = x17*x35 + x20*x36;
    double x38 = -x29*x30 - x31*x32 - x33*x37;
    double x39 = -x17 - x20;
    double x40 = x12*x16;
    double x41 = x11*x19;
    double x42 = 1.0*x40 + 1.0*x41 + 1.0;
    double x43 = x12*x24;
    double x44 = x11*x26;
    double x45 = -1.0*x43 - 1.0*x44;
    double x46 = S11*x42 + S12*x45;
    double x47 = S12*x42 + S22*x45;
    double x48 = S13*x42 + S23*x45;
    double x49 = 1.0*PENER + 1.0*SENER;
    double x50 = -x30*x46 - x32*x47 - x37*x48 + x49;
    double x51 = x43 + x44;
    double x52 = x40 + x41;
    double x53 = x11*x36 + x12*x35;
    double x54 = -x29*x51 - x31*x52 - x33*x53 + x49;
    double x55 = -x46*x51 - x47*x52 - x48*x53;
    
    res_0[0] = x9*(x13*x38 + x39*x50);
    res_0[1] = x9*(x13*x54 + x39*x55);
    res_0[2] = 0;
    res_0[3] = x9*(x12*x38 + x17*x50);
    res_0[4] = x9*(x12*x54 + x17*x55);
    res_0[5] = 0;
    res_0[6] = x9*(x11*x38 + x20*x50);
    res_0[7] = x9*(x11*x54 + x20*x55);
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
    double *x2 = coord;
    double x3 = 4*s;
    double x4 = x0 + x3 - 3;
    double x5 = x4*x2[0];
    double x6 = -8*r - x3 + 4;
    double x7 = x2[9];
    double x8 = x2[15];
    double x9 = x2[12];
    double x10 = x1*x2[3] - x3*x8 + x3*x9 + x5 + x6*x7;
    double x11 = x3 - 1;
    double x12 = x4*x2[1];
    double x13 = -8*s - x0 + 4;
    double x14 = x2[16];
    double x15 = x2[10];
    double x16 = x2[13];
    double x17 = -x0*x15 + x0*x16 + x11*x2[7] + x12 + x13*x14;
    double x18 = x1*x2[4];
    double x19 = x15*x6;
    double x20 = x14*x3;
    double x21 = x16*x3;
    double x22 = x11*x2[6];
    double x23 = x13*x8;
    double x24 = x0*x7;
    double x25 = x0*x9;
    double x26 = x10*x17 - (x12 + x18 + x19 - x20 + x21)*(x22 + x23 - x24 + x25 + x5);
    double x27 = 1.0/x26;
    double x28 = x27*x4;
    double x29 = -x22 - x23 + x24 - x25 - x5;
    double x30 = x10*x28 + x28*x29;
    double *x31 = U;
    double x32 = x4*x31[1];
    double x33 = x31[10];
    double x34 = x31[16];
    double x35 = x31[13];
    double x36 = x1*x31[4] - x3*x34 + x3*x35 + x32 + x33*x6;
    double x37 = x17*x27;
    double x38 = x36*x37;
    double x39 = -x0*x33 + x0*x35 + x11*x31[7] + x13*x34 + x32;
    double x40 = -x12 - x18 - x19 + x20 - x21;
    double x41 = x27*x40;
    double x42 = x39*x41;
    double x43 = -1.0*x38 - 1.0*x42;
    double x44 = x4*x31[0];
    double x45 = x31[9];
    double x46 = x31[15];
    double x47 = x31[12];
    double x48 = x1*x31[3] - x3*x46 + x3*x47 + x44 + x45*x6;
    double x49 = x37*x48;
    double x50 = -x0*x45 + x0*x47 + x11*x31[6] + x13*x46 + x44;
    double x51 = x41*x50;
    double x52 = 1.0*x49 + 1.0*x51 + 1.0;
    double x53 = S11*x43 + S12*x52;
    double x54 = x49 + x51;
    double x55 = S12*x43 + S22*x52;
    double x56 = x38 + x42;
    double x57 = S13*x43 + S23*x52;
    double x58 = x4*x31[2];
    double x59 = x31[11];
    double x60 = x31[17];
    double x61 = x31[14];
    double x62 = x1*x31[5] - x3*x60 + x3*x61 + x58 + x59*x6;
    double x63 = -x0*x59 + x0*x61 + x11*x31[8] + x13*x60 + x58;
    double x64 = x37*x62 + x41*x63;
    double x65 = -x53*x54 - x55*x56 - x57*x64;
    double x66 = x17*x28 + x28*x40;
    double x67 = x27*x29;
    double x68 = x36*x67;
    double x69 = x10*x27;
    double x70 = x39*x69;
    double x71 = 1.0*x68 + 1.0*x70 + 1.0;
    double x72 = x48*x67;
    double x73 = x50*x69;
    double x74 = -1.0*x72 - 1.0*x73;
    double x75 = S11*x71 + S12*x74;
    double x76 = S12*x71 + S22*x74;
    double x77 = S13*x71 + S23*x74;
    double x78 = r*r;
    double x79 = 2*x78;
    double x80 = -r + x79;
    double *x81 = V;
    double x82 = x81[3];
    double x83 = s*s;
    double x84 = 2*x83;
    double x85 = -s + x84;
    double x86 = x81[6];
    double x87 = s*x0;
    double x88 = -x87;
    double x89 = x0 - 4*x78 + x88;
    double x90 = x81[9];
    double x91 = x3 - 4*x83 + x88;
    double x92 = x81[15];
    double x93 = -3*r - 3*s + x79 + x84 + x87 + 1;
    double x94 = x81[0];
    double x95 = x81[12];
    double x96 = x80*x82 + x85*x86 + x87*x95 + x89*x90 + x91*x92 + x93*x94;
    double x97 = x81[4];
    double x98 = x81[7];
    double x99 = x81[10];
    double x100 = x81[16];
    double x101 = x81[1];
    double x102 = x81[13];
    double x103 = x100*x91 + x101*x93 + x102*x87 + x80*x97 + x85*x98 + x89*x99;
    double x104 = x81[5];
    double x105 = x81[8];
    double x106 = x81[11];
    double x107 = x81[17];
    double x108 = x81[2];
    double x109 = x81[14];
    double x110 = x104*x80 + x105*x85 + x106*x89 + x107*x91 + x108*x93 + x109*x87;
    double x111 = 1.0*PENER + 1.0*SENER - 0.5*rho*(x103*x103 + x110*x110 + x96*x96);
    double x112 = x111 - x54*x75 - x56*x76 - x64*x77;
    double *x113 = A;
    double x114 = x80*x113[3] + x85*x113[6] + x87*x113[12] + x89*x113[9] + x91*x113[15] + x93*x113[0];
    double x115 = x80*x113[4] + x85*x113[7] + x87*x113[13] + x89*x113[10] + x91*x113[16] + x93*x113[1];
    double x116 = x80*x113[5] + x85*x113[8] + x87*x113[14] + x89*x113[11] + x91*x113[17] + x93*x113[2];
    double x117 = x4*x94;
    double x118 = x1*x82 + x117 - x3*x92 + x3*x95 + x6*x90;
    double x119 = -x0*x90 + x0*x95 + x11*x86 + x117 + x13*x92;
    double x120 = x101*x4;
    double x121 = x1*x97 - x100*x3 + x102*x3 + x120 + x6*x99;
    double x122 = x0*x102 - x0*x99 + x100*x13 + x11*x98 + x120;
    double x123 = x108*x4;
    double x124 = x1*x104 + x106*x6 - x107*x3 + x109*x3 + x123;
    double x125 = -x0*x106 + x0*x109 + x105*x11 + x107*x13 + x123;
    double x126 = x103*(x121*x37 + x122*x41) + x110*(x124*x37 + x125*x41) + x114*x54 + x115*x56 + x116*x64 + x96*(x118*x37 + x119*x41);
    double x127 = rho*x93;
    double x128 = x72 + x73;
    double x129 = x68 + x70;
    double x130 = x62*x67 + x63*x69;
    double x131 = x111 - x128*x53 - x129*x55 - x130*x57;
    double x132 = -x128*x75 - x129*x76 - x130*x77;
    double x133 = x103*(x121*x67 + x122*x69) + x110*(x124*x67 + x125*x69) + x114*x128 + x115*x129 + x116*x130 + x96*(x118*x67 + x119*x69);
    double x134 = rho*x80;
    double x135 = x1*x67;
    double x136 = x1*x37;
    double x137 = rho*x85;
    double x138 = x11*x69;
    double x139 = x11*x41;
    double x140 = x0*x41;
    double x141 = -x140 + x37*x6;
    double x142 = x0*x69;
    double x143 = -x142 + x6*x67;
    double x144 = rho*x89;
    double x145 = x3*x67;
    double x146 = x142 + x145;
    double x147 = x3*x37;
    double x148 = x140 + x147;
    double x149 = rho*x87;
    double x150 = x13*x69 - x145;
    double x151 = x13*x41 - x147;
    double x152 = rho*x91;
    
    res_0[0] = x26*(x112*x66 - x126*x127 + x30*x65);
    res_0[1] = x26*(-x127*x133 + x131*x30 + x132*x66);
    res_0[2] = 0;
    res_0[3] = x26*(x112*x136 - x126*x134 + x135*x65);
    res_0[4] = x26*(x131*x135 + x132*x136 - x133*x134);
    res_0[5] = 0;
    res_0[6] = x26*(x112*x139 - x126*x137 + x138*x65);
    res_0[7] = x26*(x131*x138 + x132*x139 - x133*x137);
    res_0[8] = 0;
    res_0[9] = x26*(x112*x141 - x126*x144 + x143*x65);
    res_0[10] = x26*(x131*x143 + x132*x141 - x133*x144);
    res_0[11] = 0;
    res_0[12] = x26*(x112*x148 - x126*x149 + x146*x65);
    res_0[13] = x26*(x131*x146 + x132*x148 - x133*x149);
    res_0[14] = 0;
    res_0[15] = x26*(x112*x151 - x126*x152 + x150*x65);
    res_0[16] = x26*(x131*x150 + x132*x151 - x133*x152);
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
    double *x2 = coord;
    double x3 = 4*s;
    double x4 = x0 + x3 - 3;
    double x5 = x4*x2[0];
    double x6 = -8*r - x3 + 4;
    double x7 = x2[9];
    double x8 = x2[15];
    double x9 = x2[12];
    double x10 = x1*x2[3] - x3*x8 + x3*x9 + x5 + x6*x7;
    double x11 = x3 - 1;
    double x12 = x4*x2[1];
    double x13 = -8*s - x0 + 4;
    double x14 = x2[16];
    double x15 = x2[10];
    double x16 = x2[13];
    double x17 = -x0*x15 + x0*x16 + x11*x2[7] + x12 + x13*x14;
    double x18 = x1*x2[4];
    double x19 = x15*x6;
    double x20 = x14*x3;
    double x21 = x16*x3;
    double x22 = x11*x2[6];
    double x23 = x13*x8;
    double x24 = x0*x7;
    double x25 = x0*x9;
    double x26 = x10*x17 - (x12 + x18 + x19 - x20 + x21)*(x22 + x23 - x24 + x25 + x5);
    double x27 = 1.0/x26;
    double x28 = x27*x4;
    double x29 = -x22 - x23 + x24 - x25 - x5;
    double x30 = x10*x28 + x28*x29;
    double *x31 = U;
    double x32 = x4*x31[1];
    double x33 = x31[10];
    double x34 = x31[16];
    double x35 = x31[13];
    double x36 = x1*x31[4] - x3*x34 + x3*x35 + x32 + x33*x6;
    double x37 = x17*x27;
    double x38 = x36*x37;
    double x39 = -x0*x33 + x0*x35 + x11*x31[7] + x13*x34 + x32;
    double x40 = -x12 - x18 - x19 + x20 - x21;
    double x41 = x27*x40;
    double x42 = x39*x41;
    double x43 = -1.0*x38 - 1.0*x42;
    double x44 = x4*x31[0];
    double x45 = x31[9];
    double x46 = x31[15];
    double x47 = x31[12];
    double x48 = x1*x31[3] - x3*x46 + x3*x47 + x44 + x45*x6;
    double x49 = x37*x48;
    double x50 = -x0*x45 + x0*x47 + x11*x31[6] + x13*x46 + x44;
    double x51 = x41*x50;
    double x52 = 1.0*x49 + 1.0*x51 + 1.0;
    double x53 = S11*x43 + S12*x52;
    double x54 = x49 + x51 + 1.0;
    double x55 = S12*x43 + S22*x52;
    double x56 = x38 + x42;
    double x57 = S13*x43;
    double x58 = S23*x52;
    double x59 = x57 + x58;
    double x60 = x4*x31[2];
    double x61 = x31[11];
    double x62 = x31[17];
    double x63 = x31[14];
    double x64 = x1*x31[5] - x3*x62 + x3*x63 + x6*x61 + x60;
    double x65 = -x0*x61 + x0*x63 + x11*x31[8] + x13*x62 + x60;
    double x66 = x37*x64 + x41*x65;
    double x67 = -x53*x54 - x55*x56 - x59*x66;
    double x68 = x17*x28 + x28*x40;
    double x69 = x27*x29;
    double x70 = x36*x69;
    double x71 = x10*x27;
    double x72 = x39*x71;
    double x73 = 1.0*x70 + 1.0*x72 + 1.0;
    double x74 = x48*x69;
    double x75 = x50*x71;
    double x76 = -1.0*x74 - 1.0*x75;
    double x77 = S11*x73 + S12*x76;
    double x78 = S12*x73 + S22*x76;
    double x79 = S13*x73;
    double x80 = S23*x76;
    double x81 = x79 + x80;
    double x82 = 1.0*PENER + 1.0*SENER;
    double x83 = -x54*x77 - x56*x78 - x66*x81 + x82;
    double x84 = x74 + x75;
    double x85 = x70 + x72 + 1.0;
    double x86 = x64*x69 + x65*x71;
    double x87 = -x53*x84 - x55*x85 - x59*x86 + x82;
    double x88 = -x77*x84 - x78*x85 - x81*x86;
    double x89 = -1.0*x57 - 1.0*x58;
    double x90 = -1.0*x79 - 1.0*x80;
    double x91 = x1*x69;
    double x92 = x1*x37;
    double x93 = x11*x71;
    double x94 = x11*x41;
    double x95 = x0*x41;
    double x96 = x37*x6 - x95;
    double x97 = x0*x71;
    double x98 = x6*x69 - x97;
    double x99 = x3*x69;
    double x100 = x97 + x99;
    double x101 = x3*x37;
    double x102 = x101 + x95;
    double x103 = x13*x71 - x99;
    double x104 = -x101 + x13*x41;
    
    res_0[0] = x26*(x30*x67 + x68*x83);
    res_0[1] = x26*(x30*x87 + x68*x88);
    res_0[2] = x26*(x30*x89 + x68*x90);
    res_0[3] = x26*(x67*x91 + x83*x92);
    res_0[4] = x26*(x87*x91 + x88*x92);
    res_0[5] = x26*(x89*x91 + x90*x92);
    res_0[6] = x26*(x67*x93 + x83*x94);
    res_0[7] = x26*(x87*x93 + x88*x94);
    res_0[8] = x26*(x89*x93 + x90*x94);
    res_0[9] = x26*(x67*x98 + x83*x96);
    res_0[10] = x26*(x87*x98 + x88*x96);
    res_0[11] = x26*(x89*x98 + x90*x96);
    res_0[12] = x26*(x100*x67 + x102*x83);
    res_0[13] = x26*(x100*x87 + x102*x88);
    res_0[14] = x26*(x100*x89 + x102*x90);
    res_0[15] = x26*(x103*x67 + x104*x83);
    res_0[16] = x26*(x103*x87 + x104*x88);
    res_0[17] = x26*(x103*x89 + x104*x90);
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
    double *x2 = coord;
    double x3 = 4*s;
    double x4 = x0 + x3 - 3;
    double x5 = x4*x2[0];
    double x6 = -8*r - x3 + 4;
    double x7 = x2[9];
    double x8 = x2[15];
    double x9 = x2[12];
    double x10 = x1*x2[3] - x3*x8 + x3*x9 + x5 + x6*x7;
    double x11 = x3 - 1;
    double x12 = x4*x2[1];
    double x13 = -8*s - x0 + 4;
    double x14 = x2[16];
    double x15 = x2[10];
    double x16 = x2[13];
    double x17 = -x0*x15 + x0*x16 + x11*x2[7] + x12 + x13*x14;
    double x18 = x1*x2[4];
    double x19 = x15*x6;
    double x20 = x14*x3;
    double x21 = x16*x3;
    double x22 = x11*x2[6];
    double x23 = x13*x8;
    double x24 = x0*x7;
    double x25 = x0*x9;
    double x26 = x10*x17 - (x12 + x18 + x19 - x20 + x21)*(x22 + x23 - x24 + x25 + x5);
    double x27 = 1.0/x26;
    double x28 = x27*x4;
    double x29 = -x22 - x23 + x24 - x25 - x5;
    double x30 = x10*x28 + x28*x29;
    double *x31 = U;
    double x32 = x4*x31[1];
    double x33 = x31[10];
    double x34 = x31[16];
    double x35 = x31[13];
    double x36 = x1*x31[4] - x3*x34 + x3*x35 + x32 + x33*x6;
    double x37 = x17*x27;
    double x38 = x36*x37;
    double x39 = -x0*x33 + x0*x35 + x11*x31[7] + x13*x34 + x32;
    double x40 = -x12 - x18 - x19 + x20 - x21;
    double x41 = x27*x40;
    double x42 = x39*x41;
    double x43 = -1.0*x38 - 1.0*x42;
    double x44 = x4*x31[0];
    double x45 = x31[9];
    double x46 = x31[15];
    double x47 = x31[12];
    double x48 = x1*x31[3] - x3*x46 + x3*x47 + x44 + x45*x6;
    double x49 = x37*x48;
    double x50 = -x0*x45 + x0*x47 + x11*x31[6] + x13*x46 + x44;
    double x51 = x41*x50;
    double x52 = 1.0*x49 + 1.0*x51 + 1.0;
    double x53 = S11*x43 + S12*x52;
    double x54 = x49 + x51;
    double x55 = S12*x43 + S22*x52;
    double x56 = x38 + x42;
    double x57 = S13*x43 + S23*x52;
    double x58 = x4*x31[2];
    double x59 = x31[11];
    double x60 = x31[17];
    double x61 = x31[14];
    double x62 = x1*x31[5] - x3*x60 + x3*x61 + x58 + x59*x6;
    double x63 = -x0*x59 + x0*x61 + x11*x31[8] + x13*x60 + x58;
    double x64 = x37*x62 + x41*x63;
    double x65 = -x53*x54 - x55*x56 - x57*x64;
    double x66 = x17*x28 + x28*x40;
    double x67 = x27*x29;
    double x68 = x36*x67;
    double x69 = x10*x27;
    double x70 = x39*x69;
    double x71 = 1.0*x68 + 1.0*x70 + 1.0;
    double x72 = x48*x67;
    double x73 = x50*x69;
    double x74 = -1.0*x72 - 1.0*x73;
    double x75 = S11*x71 + S12*x74;
    double x76 = S12*x71 + S22*x74;
    double x77 = S13*x71 + S23*x74;
    double x78 = 1.0*PENER + 1.0*SENER;
    double x79 = -x54*x75 - x56*x76 - x64*x77 + x78;
    double x80 = x72 + x73;
    double x81 = x68 + x70;
    double x82 = x62*x67 + x63*x69;
    double x83 = -x53*x80 - x55*x81 - x57*x82 + x78;
    double x84 = -x75*x80 - x76*x81 - x77*x82;
    double x85 = x1*x67;
    double x86 = x1*x37;
    double x87 = x11*x69;
    double x88 = x11*x41;
    double x89 = x0*x41;
    double x90 = x37*x6 - x89;
    double x91 = x0*x69;
    double x92 = x6*x67 - x91;
    double x93 = x3*x67;
    double x94 = x91 + x93;
    double x95 = x3*x37;
    double x96 = x89 + x95;
    double x97 = x13*x69 - x93;
    double x98 = x13*x41 - x95;
    
    res_0[0] = x26*(x30*x65 + x66*x79);
    res_0[1] = x26*(x30*x83 + x66*x84);
    res_0[2] = 0;
    res_0[3] = x26*(x65*x85 + x79*x86);
    res_0[4] = x26*(x83*x85 + x84*x86);
    res_0[5] = 0;
    res_0[6] = x26*(x65*x87 + x79*x88);
    res_0[7] = x26*(x83*x87 + x84*x88);
    res_0[8] = 0;
    res_0[9] = x26*(x65*x92 + x79*x90);
    res_0[10] = x26*(x83*x92 + x84*x90);
    res_0[11] = 0;
    res_0[12] = x26*(x65*x94 + x79*x96);
    res_0[13] = x26*(x83*x94 + x84*x96);
    res_0[14] = 0;
    res_0[15] = x26*(x65*x97 + x79*x98);
    res_0[16] = x26*(x83*x97 + x84*x98);
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
    double *x2 = coord;
    double x3 = x2[12];
    double x4 = -x0;
    double x5 = x4 + 1.0/2.0;
    double x6 = x2[3];
    double x7 = x0 - 1.0/2.0;
    double x8 = x2[0];
    double x9 = x4 - 1.0/2.0;
    double x10 = x2[9];
    double x11 = x10*x9 + x7*x8;
    double x12 = x1*x3 + x11 + x5*x6;
    double x13 = x2[16];
    double x14 = x2[7];
    double x15 = x2[1];
    double x16 = x2[10];
    double x17 = x15*x7 + x16*x9;
    double x18 = x1*x13 + x14*x5 + x17;
    double x19 = x12*x18;
    double x20 = x2[13];
    double x21 = x2[4];
    double x22 = x1*x20 + x17 + x21*x5;
    double x23 = x2[15];
    double x24 = x2[6];
    double x25 = x1*x23 + x11 + x24*x5;
    double x26 = x22*x25;
    double x27 = x19 - x26;
    double x28 = (1.0/2.0)*r;
    double x29 = (1.0/2.0)*s;
    double x30 = -x28 - x29;
    double x31 = x30 + 1.0/2.0;
    double x32 = x2[11];
    double x33 = x28 + x29 - 1.0/2.0;
    double x34 = x2[2];
    double x35 = x2[14];
    double x36 = x2[17];
    double x37 = x2[5];
    double x38 = x2[8];
    double x39 = x28*x35 - x28*x37 + x29*x36 - x29*x38 + x31*x32 + x33*x34;
    double x40 = x32*x9 + x34*x7;
    double x41 = x1*x36 + x38*x5 + x40;
    double x42 = x10*x31 + x23*x29 - x24*x29 + x28*x3 - x28*x6 + x33*x8;
    double x43 = x22*x42;
    double x44 = x1*x35 + x37*x5 + x40;
    double x45 = x13*x29 - x14*x29 + x15*x33 + x16*x31 + x20*x28 - x21*x28;
    double x46 = x25*x45;
    double x47 = x12*x45;
    double x48 = x18*x42;
    double x49 = x19*x39 - x26*x39 + x41*x43 - x41*x47 + x44*x46 - x44*x48;
    double x50 = 1.0/x49;
    double *x51 = U;
    double x52 = x51[10];
    double x53 = x51[1];
    double x54 = x51[13];
    double x55 = x51[16];
    double x56 = x51[4];
    double x57 = x51[7];
    double x58 = x50*(x28*x54 - x28*x56 + x29*x55 - x29*x57 + x31*x52 + x33*x53);
    double x59 = x43 - x47;
    double x60 = x52*x9 + x53*x7;
    double x61 = x50*(x1*x55 + x5*x57 + x60);
    double x62 = x46 - x48;
    double x63 = x50*(x1*x54 + x5*x56 + x60);
    double x64 = x27*x58 + x59*x61 + x62*x63;
    double x65 = -x18*x44 + x22*x41;
    double x66 = x51[11];
    double x67 = x51[2];
    double x68 = x51[14];
    double x69 = x51[17];
    double x70 = x51[5];
    double x71 = x51[8];
    double x72 = x50*(x28*x68 - x28*x70 + x29*x69 - x29*x71 + x31*x66 + x33*x67);
    double x73 = -x22*x39 + x44*x45;
    double x74 = x66*x9 + x67*x7;
    double x75 = x50*(x1*x69 + x5*x71 + x74);
    double x76 = x18*x39 - x41*x45;
    double x77 = x50*(x1*x68 + x5*x70 + x74);
    double x78 = x65*x72 + x73*x75 + x76*x77;
    double x79 = x58*x65 + x61*x73 + x63*x76;
    double x80 = x27*x72 + x59*x75 + x62*x77;
    double x81 = x80 + 1.0;
    double x82 = x64*x78 - x79*x81;
    double x83 = x51[9];
    double x84 = x51[0];
    double x85 = x51[12];
    double x86 = x51[15];
    double x87 = x51[3];
    double x88 = x51[6];
    double x89 = x50*(x28*x85 - x28*x87 + x29*x86 - x29*x88 + x31*x83 + x33*x84);
    double x90 = x7*x84 + x83*x9;
    double x91 = x50*(x1*x86 + x5*x88 + x90);
    double x92 = x50*(x1*x85 + x5*x87 + x90);
    double x93 = x65*x89 + x73*x91 + x76*x92;
    double x94 = x93 + 1.0;
    double x95 = x27*x89 + x59*x91 + x62*x92;
    double x96 = -x78*x95 + x81*x94;
    double x97 = -x64*x94 + x79*x95;
    double x98 = S11*x82 + S12*x96 + S13*x97;
    double x99 = S12*x82 + S22*x96 + S23*x97;
    double x100 = S13*x82 + S23*x96 + S33*x97;
    double x101 = -x100*x78 - x79*x99 - x93*x98;
    double x102 = x12*x39 - x42*x44;
    double x103 = x50*x7;
    double x104 = -x25*x39 + x41*x42;
    double x105 = -x12*x41 + x25*x44;
    double x106 = x33*x50;
    double x107 = x102*x103 + x103*x104 + x105*x106;
    double x108 = x102*x75 + x104*x77 + x105*x72;
    double x109 = x102*x61 + x104*x63 + x105*x58;
    double x110 = x109 + 1.0;
    double x111 = x108*x79 - x110*x78;
    double x112 = x102*x91 + x104*x92 + x105*x89;
    double x113 = -x108*x94 + x112*x78;
    double x114 = x110*x94 - x112*x79;
    double x115 = S11*x111 + S12*x113 + S13*x114;
    double x116 = S12*x111 + S22*x113 + S23*x114;
    double x117 = S13*x111 + S23*x113 + S33*x114;
    double x118 = -x115*x93 - x116*x79 - x117*x78;
    double x119 = x103*x59 + x103*x62 + x106*x27;
    double x120 = x103*x73 + x103*x76 + x106*x65;
    double x121 = -x108*x64 + x110*x81;
    double x122 = x108*x95 - x112*x81;
    double x123 = -x110*x95 + x112*x64;
    double x124 = S11*x121 + S12*x122 + S13*x123;
    double x125 = S12*x121 + S22*x122 + S23*x123;
    double x126 = S13*x121 + S23*x122 + S33*x123;
    double x127 = r*x0;
    double x128 = x127 + x28;
    double *x129 = V;
    double x130 = x129[12];
    double x131 = -x127;
    double x132 = x131 + x28;
    double x133 = x129[3];
    double x134 = s*x0;
    double x135 = x134 + x29;
    double x136 = x129[15];
    double x137 = -x134;
    double x138 = x137 + x29;
    double x139 = x129[6];
    double x140 = x1 + x131 + x137 + x30;
    double x141 = x129[9];
    double x142 = x127 + x134 + x30 + x5;
    double x143 = x129[0];
    double x144 = x128*x130 + x132*x133 + x135*x136 + x138*x139 + x140*x141 + x142*x143;
    double x145 = x129[13];
    double x146 = x129[4];
    double x147 = x129[16];
    double x148 = x129[7];
    double x149 = x129[10];
    double x150 = x129[1];
    double x151 = x128*x145 + x132*x146 + x135*x147 + x138*x148 + x140*x149 + x142*x150;
    double x152 = x129[14];
    double x153 = x129[5];
    double x154 = x129[17];
    double x155 = x129[8];
    double x156 = x129[11];
    double x157 = x129[2];
    double x158 = x128*x152 + x132*x153 + x135*x154 + x138*x155 + x140*x156 + x142*x157;
    double x159 = 1.0*PENER + 1.0*SENER - 0.5*rho*(x144*x144 + x151*x151 + x158*x158);
    double x160 = -x124*x93 - x125*x79 - x126*x78 + x159;
    double *x161 = A;
    double x162 = x128*x161[12] + x132*x161[3] + x135*x161[15] + x138*x161[6] + x140*x161[9] + x142*x161[0];
    double x163 = x128*x161[13] + x132*x161[4] + x135*x161[16] + x138*x161[7] + x140*x161[10] + x142*x161[1];
    double x164 = x128*x161[14] + x132*x161[5] + x135*x161[17] + x138*x161[8] + x140*x161[11] + x142*x161[2];
    double x165 = x130*x28 - x133*x28 + x136*x29 - x139*x29 + x141*x31 + x143*x33;
    double x166 = x50*x65;
    double x167 = x141*x9 + x143*x7;
    double x168 = x1*x136 + x139*x5 + x167;
    double x169 = x50*x73;
    double x170 = x1*x130 + x133*x5 + x167;
    double x171 = x50*x76;
    double x172 = x145*x28 - x146*x28 + x147*x29 - x148*x29 + x149*x31 + x150*x33;
    double x173 = x149*x9 + x150*x7;
    double x174 = x1*x147 + x148*x5 + x173;
    double x175 = x1*x145 + x146*x5 + x173;
    double x176 = x152*x28 - x153*x28 + x154*x29 - x155*x29 + x156*x31 + x157*x33;
    double x177 = x156*x9 + x157*x7;
    double x178 = x1*x154 + x155*x5 + x177;
    double x179 = x1*x152 + x153*x5 + x177;
    double x180 = x144*(x165*x166 + x168*x169 + x170*x171) + x151*(x166*x172 + x169*x174 + x171*x175) + x158*(x166*x176 + x169*x178 + x171*x179) + x162*x93 + x163*x79 + x164*x78;
    double x181 = rho*x142;
    double x182 = -x108*x117 - x109*x116 - x112*x115;
    double x183 = -x108*x126 - x109*x125 - x112*x124;
    double x184 = -x100*x108 - x109*x99 - x112*x98 + x159;
    double x185 = x102*x50;
    double x186 = x105*x50;
    double x187 = x104*x50;
    double x188 = x108*x164 + x109*x163 + x112*x162 + x144*(x165*x186 + x168*x185 + x170*x187) + x151*(x172*x186 + x174*x185 + x175*x187) + x158*(x176*x186 + x178*x185 + x179*x187);
    double x189 = -x100*x80 - x64*x99 - x95*x98;
    double x190 = -x124*x95 - x125*x64 - x126*x80;
    double x191 = -x115*x95 - x116*x64 - x117*x80 + x159;
    double x192 = x27*x50;
    double x193 = x50*x59;
    double x194 = x50*x62;
    double x195 = x144*(x165*x192 + x168*x193 + x170*x194) + x151*(x172*x192 + x174*x193 + x175*x194) + x158*(x176*x192 + x178*x193 + x179*x194) + x162*x95 + x163*x64 + x164*x80;
    double x196 = x192*x28;
    double x197 = x194*x5 - x196;
    double x198 = x166*x28;
    double x199 = x171*x5 - x198;
    double x200 = x186*x28;
    double x201 = x187*x5 - x200;
    double x202 = rho*x132;
    double x203 = x186*x29;
    double x204 = x185*x5 - x203;
    double x205 = x192*x29;
    double x206 = x193*x5 - x205;
    double x207 = x166*x29;
    double x208 = x169*x5 - x207;
    double x209 = rho*x138;
    double x210 = x185*x9 + x186*x31 + x187*x9;
    double x211 = x192*x31 + x193*x9 + x194*x9;
    double x212 = x166*x31 + x169*x9 + x171*x9;
    double x213 = rho*x140;
    double x214 = x1*x194 + x196;
    double x215 = x1*x171 + x198;
    double x216 = x1*x187 + x200;
    double x217 = rho*x128;
    double x218 = x1*x185 + x203;
    double x219 = x1*x193 + x205;
    double x220 = x1*x169 + x207;
    double x221 = rho*x135;
    
    res_0[0] = x49*(x101*x107 + x118*x119 + x120*x160 - x180*x181);
    res_0[1] = x49*(x107*x184 + x119*x182 + x120*x183 - x181*x188);
    res_0[2] = x49*(x107*x189 + x119*x191 + x120*x190 - x181*x195);
    res_0[3] = x49*(x101*x201 + x118*x197 + x160*x199 - x180*x202);
    res_0[4] = x49*(x182*x197 + x183*x199 + x184*x201 - x188*x202);
    res_0[5] = x49*(x189*x201 + x190*x199 + x191*x197 - x195*x202);
    res_0[6] = x49*(x101*x204 + x118*x206 + x160*x208 - x180*x209);
    res_0[7] = x49*(x182*x206 + x183*x208 + x184*x204 - x188*x209);
    res_0[8] = x49*(x189*x204 + x190*x208 + x191*x206 - x195*x209);
    res_0[9] = x49*(x101*x210 + x118*x211 + x160*x212 - x180*x213);
    res_0[10] = x49*(x182*x211 + x183*x212 + x184*x210 - x188*x213);
    res_0[11] = x49*(x189*x210 + x190*x212 + x191*x211 - x195*x213);
    res_0[12] = x49*(x101*x216 + x118*x214 + x160*x215 - x180*x217);
    res_0[13] = x49*(x182*x214 + x183*x215 + x184*x216 - x188*x217);
    res_0[14] = x49*(x189*x216 + x190*x215 + x191*x214 - x195*x217);
    res_0[15] = x49*(x101*x218 + x118*x219 + x160*x220 - x180*x221);
    res_0[16] = x49*(x182*x219 + x183*x220 + x184*x218 - x188*x221);
    res_0[17] = x49*(x189*x218 + x190*x220 + x191*x219 - x195*x221);
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
    double *x2 = coord;
    double x3 = x2[12];
    double x4 = -x0;
    double x5 = x4 + 1.0/2.0;
    double x6 = x2[3];
    double x7 = x0 - 1.0/2.0;
    double x8 = x2[0];
    double x9 = x4 - 1.0/2.0;
    double x10 = x2[9];
    double x11 = x10*x9 + x7*x8;
    double x12 = x1*x3 + x11 + x5*x6;
    double x13 = x2[16];
    double x14 = x2[7];
    double x15 = x2[1];
    double x16 = x2[10];
    double x17 = x15*x7 + x16*x9;
    double x18 = x1*x13 + x14*x5 + x17;
    double x19 = x12*x18;
    double x20 = x2[13];
    double x21 = x2[4];
    double x22 = x1*x20 + x17 + x21*x5;
    double x23 = x2[15];
    double x24 = x2[6];
    double x25 = x1*x23 + x11 + x24*x5;
    double x26 = x22*x25;
    double x27 = x19 - x26;
    double x28 = (1.0/2.0)*r;
    double x29 = (1.0/2.0)*s;
    double x30 = -x28 - x29 + 1.0/2.0;
    double x31 = x2[11];
    double x32 = x28 + x29 - 1.0/2.0;
    double x33 = x2[2];
    double x34 = x2[14];
    double x35 = x2[17];
    double x36 = x2[5];
    double x37 = x2[8];
    double x38 = x28*x34 - x28*x36 + x29*x35 - x29*x37 + x30*x31 + x32*x33;
    double x39 = x31*x9 + x33*x7;
    double x40 = x1*x35 + x37*x5 + x39;
    double x41 = x10*x30 + x23*x29 - x24*x29 + x28*x3 - x28*x6 + x32*x8;
    double x42 = x22*x41;
    double x43 = x1*x34 + x36*x5 + x39;
    double x44 = x13*x29 - x14*x29 + x15*x32 + x16*x30 + x20*x28 - x21*x28;
    double x45 = x25*x44;
    double x46 = x12*x44;
    double x47 = x18*x41;
    double x48 = x19*x38 - x26*x38 + x40*x42 - x40*x46 + x43*x45 - x43*x47;
    double x49 = 1.0/x48;
    double *x50 = U;
    double x51 = x50[10];
    double x52 = x50[1];
    double x53 = x50[13];
    double x54 = x50[16];
    double x55 = x50[4];
    double x56 = x50[7];
    double x57 = x49*(x28*x53 - x28*x55 + x29*x54 - x29*x56 + x30*x51 + x32*x52);
    double x58 = x42 - x46;
    double x59 = x51*x9 + x52*x7;
    double x60 = x49*(x1*x54 + x5*x56 + x59);
    double x61 = x45 - x47;
    double x62 = x49*(x1*x53 + x5*x55 + x59);
    double x63 = x27*x57 + x58*x60 + x61*x62;
    double x64 = -x18*x43 + x22*x40;
    double x65 = x50[11];
    double x66 = x50[2];
    double x67 = x50[14];
    double x68 = x50[17];
    double x69 = x50[5];
    double x70 = x50[8];
    double x71 = x49*(x28*x67 - x28*x69 + x29*x68 - x29*x70 + x30*x65 + x32*x66);
    double x72 = -x22*x38 + x43*x44;
    double x73 = x65*x9 + x66*x7;
    double x74 = x49*(x1*x68 + x5*x70 + x73);
    double x75 = x18*x38 - x40*x44;
    double x76 = x49*(x1*x67 + x5*x69 + x73);
    double x77 = x64*x71 + x72*x74 + x75*x76;
    double x78 = x57*x64 + x60*x72 + x62*x75;
    double x79 = x27*x71 + x58*x74 + x61*x76 + 1.0;
    double x80 = x63*x77 - x78*x79;
    double x81 = x50[9];
    double x82 = x50[0];
    double x83 = x50[12];
    double x84 = x50[15];
    double x85 = x50[3];
    double x86 = x50[6];
    double x87 = x49*(x28*x83 - x28*x85 + x29*x84 - x29*x86 + x30*x81 + x32*x82);
    double x88 = x7*x82 + x81*x9;
    double x89 = x49*(x1*x84 + x5*x86 + x88);
    double x90 = x49*(x1*x83 + x5*x85 + x88);
    double x91 = x64*x87 + x72*x89 + x75*x90 + 1.0;
    double x92 = x27*x87 + x58*x89 + x61*x90;
    double x93 = -x77*x92 + x79*x91;
    double x94 = -x63*x91 + x78*x92;
    double x95 = S11*x80 + S12*x93 + S13*x94;
    double x96 = S12*x80 + S22*x93 + S23*x94;
    double x97 = S13*x80 + S23*x93 + S33*x94;
    double x98 = -x77*x97 - x78*x96 - x91*x95;
    double x99 = x12*x38 - x41*x43;
    double x100 = x49*x7;
    double x101 = -x25*x38 + x40*x41;
    double x102 = -x12*x40 + x25*x43;
    double x103 = x32*x49;
    double x104 = x100*x101 + x100*x99 + x102*x103;
    double x105 = x101*x76 + x102*x71 + x74*x99;
    double x106 = x101*x62 + x102*x57 + x60*x99 + 1.0;
    double x107 = x105*x78 - x106*x77;
    double x108 = x101*x90 + x102*x87 + x89*x99;
    double x109 = -x105*x91 + x108*x77;
    double x110 = x106*x91 - x108*x78;
    double x111 = S11*x107 + S12*x109 + S13*x110;
    double x112 = S12*x107 + S22*x109 + S23*x110;
    double x113 = S13*x107 + S23*x109 + S33*x110;
    double x114 = -x111*x91 - x112*x78 - x113*x77;
    double x115 = x100*x58 + x100*x61 + x103*x27;
    double x116 = x100*x72 + x100*x75 + x103*x64;
    double x117 = -x105*x63 + x106*x79;
    double x118 = x105*x92 - x108*x79;
    double x119 = -x106*x92 + x108*x63;
    double x120 = S11*x117 + S12*x118 + S13*x119;
    double x121 = S12*x117 + S22*x118 + S23*x119;
    double x122 = S13*x117 + S23*x118 + S33*x119;
    double x123 = 1.0*PENER + 1.0*SENER;
    double x124 = -x120*x91 - x121*x78 - x122*x77 + x123;
    double x125 = -x105*x113 - x106*x112 - x108*x111;
    double x126 = -x105*x122 - x106*x121 - x108*x120;
    double x127 = -x105*x97 - x106*x96 - x108*x95 + x123;
    double x128 = -x63*x96 - x79*x97 - x92*x95;
    double x129 = -x120*x92 - x121*x63 - x122*x79;
    double x130 = -x111*x92 - x112*x63 - x113*x79 + x123;
    double x131 = x49*x5;
    double x132 = x28*x49;
    double x133 = x132*x27;
    double x134 = x131*x61 - x133;
    double x135 = x132*x64;
    double x136 = x131*x75 - x135;
    double x137 = x102*x132;
    double x138 = x101*x131 - x137;
    double x139 = x29*x49;
    double x140 = x102*x139;
    double x141 = x131*x99 - x140;
    double x142 = x139*x27;
    double x143 = x131*x58 - x142;
    double x144 = x139*x64;
    double x145 = x131*x72 - x144;
    double x146 = x49*x9;
    double x147 = x30*x49;
    double x148 = x101*x146 + x102*x147 + x146*x99;
    double x149 = x146*x58 + x146*x61 + x147*x27;
    double x150 = x146*x72 + x146*x75 + x147*x64;
    double x151 = x1*x49;
    double x152 = x133 + x151*x61;
    double x153 = x135 + x151*x75;
    double x154 = x101*x151 + x137;
    double x155 = x140 + x151*x99;
    double x156 = x142 + x151*x58;
    double x157 = x144 + x151*x72;
    
    res_0[0] = x48*(x104*x98 + x114*x115 + x116*x124);
    res_0[1] = x48*(x104*x127 + x115*x125 + x116*x126);
    res_0[2] = x48*(x104*x128 + x115*x130 + x116*x129);
    res_0[3] = x48*(x114*x134 + x124*x136 + x138*x98);
    res_0[4] = x48*(x125*x134 + x126*x136 + x127*x138);
    res_0[5] = x48*(x128*x138 + x129*x136 + x130*x134);
    res_0[6] = x48*(x114*x143 + x124*x145 + x141*x98);
    res_0[7] = x48*(x125*x143 + x126*x145 + x127*x141);
    res_0[8] = x48*(x128*x141 + x129*x145 + x130*x143);
    res_0[9] = x48*(x114*x149 + x124*x150 + x148*x98);
    res_0[10] = x48*(x125*x149 + x126*x150 + x127*x148);
    res_0[11] = x48*(x128*x148 + x129*x150 + x130*x149);
    res_0[12] = x48*(x114*x152 + x124*x153 + x154*x98);
    res_0[13] = x48*(x125*x152 + x126*x153 + x127*x154);
    res_0[14] = x48*(x128*x154 + x129*x153 + x130*x152);
    res_0[15] = x48*(x114*x156 + x124*x157 + x155*x98);
    res_0[16] = x48*(x125*x156 + x126*x157 + x127*x155);
    res_0[17] = x48*(x128*x155 + x129*x157 + x130*x156);
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
    double *x2 = coord;
    double x3 = x2[12];
    double x4 = -x0;
    double x5 = x4 + 1.0/2.0;
    double x6 = x2[3];
    double x7 = x0 - 1.0/2.0;
    double x8 = x2[0];
    double x9 = x4 - 1.0/2.0;
    double x10 = x2[9];
    double x11 = x10*x9 + x7*x8;
    double x12 = x1*x3 + x11 + x5*x6;
    double x13 = x2[16];
    double x14 = x2[7];
    double x15 = x2[1];
    double x16 = x2[10];
    double x17 = x15*x7 + x16*x9;
    double x18 = x1*x13 + x14*x5 + x17;
    double x19 = x12*x18;
    double x20 = x2[13];
    double x21 = x2[4];
    double x22 = x1*x20 + x17 + x21*x5;
    double x23 = x2[15];
    double x24 = x2[6];
    double x25 = x1*x23 + x11 + x24*x5;
    double x26 = x22*x25;
    double x27 = x19 - x26;
    double x28 = (1.0/2.0)*r;
    double x29 = (1.0/2.0)*s;
    double x30 = -x28 - x29 + 1.0/2.0;
    double x31 = x2[11];
    double x32 = x28 + x29 - 1.0/2.0;
    double x33 = x2[2];
    double x34 = x2[14];
    double x35 = x2[17];
    double x36 = x2[5];
    double x37 = x2[8];
    double x38 = x28*x34 - x28*x36 + x29*x35 - x29*x37 + x30*x31 + x32*x33;
    double x39 = x31*x9 + x33*x7;
    double x40 = x1*x35 + x37*x5 + x39;
    double x41 = x10*x30 + x23*x29 - x24*x29 + x28*x3 - x28*x6 + x32*x8;
    double x42 = x22*x41;
    double x43 = x1*x34 + x36*x5 + x39;
    double x44 = x13*x29 - x14*x29 + x15*x32 + x16*x30 + x20*x28 - x21*x28;
    double x45 = x25*x44;
    double x46 = x12*x44;
    double x47 = x18*x41;
    double x48 = x19*x38 - x26*x38 + x40*x42 - x40*x46 + x43*x45 - x43*x47;
    double x49 = 1.0/x48;
    double *x50 = U;
    double x51 = x50[10];
    double x52 = x50[1];
    double x53 = x50[13];
    double x54 = x50[16];
    double x55 = x50[4];
    double x56 = x50[7];
    double x57 = x49*(x28*x53 - x28*x55 + x29*x54 - x29*x56 + x30*x51 + x32*x52);
    double x58 = x42 - x46;
    double x59 = x51*x9 + x52*x7;
    double x60 = x49*(x1*x54 + x5*x56 + x59);
    double x61 = x45 - x47;
    double x62 = x49*(x1*x53 + x5*x55 + x59);
    double x63 = x27*x57 + x58*x60 + x61*x62;
    double x64 = -x18*x43 + x22*x40;
    double x65 = x50[11];
    double x66 = x50[2];
    double x67 = x50[14];
    double x68 = x50[17];
    double x69 = x50[5];
    double x70 = x50[8];
    double x71 = x49*(x28*x67 - x28*x69 + x29*x68 - x29*x70 + x30*x65 + x32*x66);
    double x72 = -x22*x38 + x43*x44;
    double x73 = x65*x9 + x66*x7;
    double x74 = x49*(x1*x68 + x5*x70 + x73);
    double x75 = x18*x38 - x40*x44;
    double x76 = x49*(x1*x67 + x5*x69 + x73);
    double x77 = x64*x71 + x72*x74 + x75*x76;
    double x78 = x57*x64 + x60*x72 + x62*x75;
    double x79 = x27*x71 + x58*x74 + x61*x76;
    double x80 = x79 + 1.0;
    double x81 = x63*x77 - x78*x80;
    double x82 = x50[9];
    double x83 = x50[0];
    double x84 = x50[12];
    double x85 = x50[15];
    double x86 = x50[3];
    double x87 = x50[6];
    double x88 = x49*(x28*x84 - x28*x86 + x29*x85 - x29*x87 + x30*x82 + x32*x83);
    double x89 = x7*x83 + x82*x9;
    double x90 = x49*(x1*x85 + x5*x87 + x89);
    double x91 = x49*(x1*x84 + x5*x86 + x89);
    double x92 = x64*x88 + x72*x90 + x75*x91;
    double x93 = x92 + 1.0;
    double x94 = x27*x88 + x58*x90 + x61*x91;
    double x95 = -x77*x94 + x80*x93;
    double x96 = -x63*x93 + x78*x94;
    double x97 = S11*x81 + S12*x95 + S13*x96;
    double x98 = S12*x81 + S22*x95 + S23*x96;
    double x99 = S13*x81 + S23*x95 + S33*x96;
    double x100 = -x77*x99 - x78*x98 - x92*x97;
    double x101 = x12*x38 - x41*x43;
    double x102 = x49*x7;
    double x103 = -x25*x38 + x40*x41;
    double x104 = -x12*x40 + x25*x43;
    double x105 = x32*x49;
    double x106 = x101*x102 + x102*x103 + x104*x105;
    double x107 = x101*x74 + x103*x76 + x104*x71;
    double x108 = x101*x60 + x103*x62 + x104*x57;
    double x109 = x108 + 1.0;
    double x110 = x107*x78 - x109*x77;
    double x111 = x101*x90 + x103*x91 + x104*x88;
    double x112 = -x107*x93 + x111*x77;
    double x113 = x109*x93 - x111*x78;
    double x114 = S11*x110 + S12*x112 + S13*x113;
    double x115 = S12*x110 + S22*x112 + S23*x113;
    double x116 = S13*x110 + S23*x112 + S33*x113;
    double x117 = -x114*x92 - x115*x78 - x116*x77;
    double x118 = x102*x58 + x102*x61 + x105*x27;
    double x119 = x102*x72 + x102*x75 + x105*x64;
    double x120 = -x107*x63 + x109*x80;
    double x121 = x107*x94 - x111*x80;
    double x122 = -x109*x94 + x111*x63;
    double x123 = S11*x120 + S12*x121 + S13*x122;
    double x124 = S12*x120 + S22*x121 + S23*x122;
    double x125 = S13*x120 + S23*x121 + S33*x122;
    double x126 = 1.0*PENER + 1.0*SENER;
    double x127 = -x123*x92 - x124*x78 - x125*x77 + x126;
    double x128 = -x107*x116 - x108*x115 - x111*x114;
    double x129 = -x107*x125 - x108*x124 - x111*x123;
    double x130 = -x107*x99 - x108*x98 - x111*x97 + x126;
    double x131 = -x63*x98 - x79*x99 - x94*x97;
    double x132 = -x123*x94 - x124*x63 - x125*x79;
    double x133 = -x114*x94 - x115*x63 - x116*x79 + x126;
    double x134 = x49*x5;
    double x135 = x28*x49;
    double x136 = x135*x27;
    double x137 = x134*x61 - x136;
    double x138 = x135*x64;
    double x139 = x134*x75 - x138;
    double x140 = x104*x135;
    double x141 = x103*x134 - x140;
    double x142 = x29*x49;
    double x143 = x104*x142;
    double x144 = x101*x134 - x143;
    double x145 = x142*x27;
    double x146 = x134*x58 - x145;
    double x147 = x142*x64;
    double x148 = x134*x72 - x147;
    double x149 = x49*x9;
    double x150 = x30*x49;
    double x151 = x101*x149 + x103*x149 + x104*x150;
    double x152 = x149*x58 + x149*x61 + x150*x27;
    double x153 = x149*x72 + x149*x75 + x150*x64;
    double x154 = x1*x49;
    double x155 = x136 + x154*x61;
    double x156 = x138 + x154*x75;
    double x157 = x103*x154 + x140;
    double x158 = x101*x154 + x143;
    double x159 = x145 + x154*x58;
    double x160 = x147 + x154*x72;
    
    res_0[0] = x48*(x100*x106 + x117*x118 + x119*x127);
    res_0[1] = x48*(x106*x130 + x118*x128 + x119*x129);
    res_0[2] = x48*(x106*x131 + x118*x133 + x119*x132);
    res_0[3] = x48*(x100*x141 + x117*x137 + x127*x139);
    res_0[4] = x48*(x128*x137 + x129*x139 + x130*x141);
    res_0[5] = x48*(x131*x141 + x132*x139 + x133*x137);
    res_0[6] = x48*(x100*x144 + x117*x146 + x127*x148);
    res_0[7] = x48*(x128*x146 + x129*x148 + x130*x144);
    res_0[8] = x48*(x131*x144 + x132*x148 + x133*x146);
    res_0[9] = x48*(x100*x151 + x117*x152 + x127*x153);
    res_0[10] = x48*(x128*x152 + x129*x153 + x130*x151);
    res_0[11] = x48*(x131*x151 + x132*x153 + x133*x152);
    res_0[12] = x48*(x100*x157 + x117*x155 + x127*x156);
    res_0[13] = x48*(x128*x155 + x129*x156 + x130*x157);
    res_0[14] = x48*(x131*x157 + x132*x156 + x133*x155);
    res_0[15] = x48*(x100*x158 + x117*x159 + x127*x160);
    res_0[16] = x48*(x128*x159 + x129*x160 + x130*x158);
    res_0[17] = x48*(x131*x158 + x132*x160 + x133*x159);
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
    
    double x0 = (1.0/8.0)*r;
    double x1 = -x0;
    double x2 = x1 - 1.0/8.0;
    double x3 = (1.0/8.0)*s;
    double x4 = -x3;
    double x5 = s*x0;
    double x6 = -x5;
    double x7 = x4 + x6;
    double x8 = x2 + x7;
    double *x9 = coord;
    double x10 = x9[6];
    double x11 = x3 + x5;
    double x12 = x11 + x2;
    double x13 = x9[3];
    double x14 = x0 - 1.0/8.0;
    double x15 = x4 + x5;
    double x16 = x14 + x15;
    double x17 = x9[9];
    double x18 = x3 + x6;
    double x19 = x14 + x18;
    double x20 = x9[0];
    double x21 = x1 + 1.0/8.0;
    double x22 = x15 + x21;
    double x23 = x9[12];
    double x24 = x18 + x21;
    double x25 = x9[21];
    double x26 = x0 + 1.0/8.0;
    double x27 = x26 + x7;
    double x28 = x9[15];
    double x29 = x11 + x26;
    double x30 = x9[18];
    double x31 = x10*x8 + x12*x13 + x16*x17 + x19*x20 + x22*x23 + x24*x25 + x27*x28 + x29*x30;
    double x32 = (1.0/8.0)*t;
    double x33 = -x32;
    double x34 = t*x0;
    double x35 = -x34;
    double x36 = x33 + x35;
    double x37 = x2 + x36;
    double x38 = x9[17];
    double x39 = x32 + x34;
    double x40 = x2 + x39;
    double x41 = x9[5];
    double x42 = x33 + x34;
    double x43 = x14 + x42;
    double x44 = x9[14];
    double x45 = x32 + x35;
    double x46 = x14 + x45;
    double x47 = x9[2];
    double x48 = x21 + x42;
    double x49 = x9[11];
    double x50 = x21 + x45;
    double x51 = x9[23];
    double x52 = x26 + x36;
    double x53 = x9[8];
    double x54 = x26 + x39;
    double x55 = x9[20];
    double x56 = x37*x38 + x40*x41 + x43*x44 + x46*x47 + x48*x49 + x50*x51 + x52*x53 + x54*x55;
    double x57 = x31*x56;
    double x58 = x12*x41 + x16*x49 + x19*x47 + x22*x44 + x24*x51 + x27*x38 + x29*x55 + x53*x8;
    double x59 = x10*x52 + x13*x40 + x17*x48 + x20*x46 + x23*x43 + x25*x50 + x28*x37 + x30*x54;
    double x60 = x58*x59;
    double x61 = x57 - x60;
    double x62 = x4 - 1.0/8.0;
    double x63 = t*x3;
    double x64 = -x63;
    double x65 = x33 + x64;
    double x66 = x62 + x65;
    double x67 = x9[22];
    double x68 = x32 + x63;
    double x69 = x62 + x68;
    double x70 = x9[10];
    double x71 = x3 - 1.0/8.0;
    double x72 = x33 + x63;
    double x73 = x71 + x72;
    double x74 = x9[13];
    double x75 = x32 + x64;
    double x76 = x71 + x75;
    double x77 = x9[1];
    double x78 = x4 + 1.0/8.0;
    double x79 = x72 + x78;
    double x80 = x9[4];
    double x81 = x75 + x78;
    double x82 = x9[16];
    double x83 = x3 + 1.0/8.0;
    double x84 = x65 + x83;
    double x85 = x9[7];
    double x86 = x68 + x83;
    double x87 = x9[19];
    double x88 = x66*x67 + x69*x70 + x73*x74 + x76*x77 + x79*x80 + x81*x82 + x84*x85 + x86*x87;
    double x89 = x12*x80 + x16*x70 + x19*x77 + x22*x74 + x24*x67 + x27*x82 + x29*x87 + x8*x85;
    double x90 = x38*x81 + x41*x79 + x44*x73 + x47*x76 + x49*x69 + x51*x66 + x53*x84 + x55*x86;
    double x91 = x59*x90;
    double x92 = x37*x82 + x40*x80 + x43*x74 + x46*x77 + x48*x70 + x50*x67 + x52*x85 + x54*x87;
    double x93 = x10*x84 + x13*x79 + x17*x69 + x20*x76 + x23*x73 + x25*x66 + x28*x81 + x30*x86;
    double x94 = x58*x93;
    double x95 = x31*x90;
    double x96 = x56*x93;
    double x97 = x57*x88 - x60*x88 + x89*x91 - x89*x96 + x92*x94 - x92*x95;
    double x98 = 1.0/x97;
    double *x99 = U;
    double x100 = x99[23];
    double x101 = x99[11];
    double x102 = x99[14];
    double x103 = x99[2];
    double x104 = x99[5];
    double x105 = x99[17];
    double x106 = x99[8];
    double x107 = x99[20];
    double x108 = x98*(x100*x66 + x101*x69 + x102*x73 + x103*x76 + x104*x79 + x105*x81 + x106*x84 + x107*x86);
    double x109 = x94 - x95;
    double x110 = x98*(x100*x50 + x101*x48 + x102*x43 + x103*x46 + x104*x40 + x105*x37 + x106*x52 + x107*x54);
    double x111 = x91 - x96;
    double x112 = x98*(x100*x24 + x101*x16 + x102*x22 + x103*x19 + x104*x12 + x105*x27 + x106*x8 + x107*x29);
    double x113 = x108*x61 + x109*x110 + x111*x112;
    double x114 = -x58*x88 + x89*x90;
    double x115 = x99[16];
    double x116 = x99[4];
    double x117 = x99[13];
    double x118 = x99[1];
    double x119 = x99[10];
    double x120 = x99[22];
    double x121 = x99[7];
    double x122 = x99[19];
    double x123 = x98*(x115*x37 + x116*x40 + x117*x43 + x118*x46 + x119*x48 + x120*x50 + x121*x52 + x122*x54);
    double x124 = -x56*x89 + x58*x92;
    double x125 = x98*(x115*x81 + x116*x79 + x117*x73 + x118*x76 + x119*x69 + x120*x66 + x121*x84 + x122*x86);
    double x126 = x56*x88 - x90*x92;
    double x127 = x98*(x115*x27 + x116*x12 + x117*x22 + x118*x19 + x119*x16 + x120*x24 + x121*x8 + x122*x29);
    double x128 = x114*x123 + x124*x125 + x126*x127;
    double x129 = x108*x124 + x110*x114 + x112*x126;
    double x130 = x109*x123 + x111*x127 + x125*x61;
    double x131 = x130 + 1.0;
    double x132 = x113*x128 - x129*x131;
    double x133 = x99[21];
    double x134 = x99[9];
    double x135 = x99[12];
    double x136 = x99[0];
    double x137 = x99[3];
    double x138 = x99[15];
    double x139 = x99[6];
    double x140 = x99[18];
    double x141 = x98*(x133*x66 + x134*x69 + x135*x73 + x136*x76 + x137*x79 + x138*x81 + x139*x84 + x140*x86);
    double x142 = x98*(x133*x50 + x134*x48 + x135*x43 + x136*x46 + x137*x40 + x138*x37 + x139*x52 + x140*x54);
    double x143 = x98*(x12*x137 + x133*x24 + x134*x16 + x135*x22 + x136*x19 + x138*x27 + x139*x8 + x140*x29);
    double x144 = x109*x142 + x111*x143 + x141*x61;
    double x145 = x114*x142 + x124*x141 + x126*x143;
    double x146 = x145 + 1.0;
    double x147 = -x113*x146 + x129*x144;
    double x148 = -x128*x144 + x131*x146;
    double x149 = S11*x132 + S12*x147 + S13*x148;
    double x150 = S12*x132 + S22*x147 + S23*x148;
    double x151 = S13*x132 + S23*x147 + S33*x148;
    double x152 = -x128*x150 - x129*x151 - x145*x149;
    double x153 = x31*x88 - x89*x93;
    double x154 = x46*x98;
    double x155 = -x31*x92 + x59*x89;
    double x156 = x76*x98;
    double x157 = -x59*x88 + x92*x93;
    double x158 = x19*x98;
    double x159 = x153*x154 + x155*x156 + x157*x158;
    double x160 = x123*x153 + x125*x155 + x127*x157;
    double x161 = x108*x155 + x110*x153 + x112*x157;
    double x162 = x161 + 1.0;
    double x163 = -x128*x162 + x129*x160;
    double x164 = x141*x155 + x142*x153 + x143*x157;
    double x165 = -x129*x164 + x146*x162;
    double x166 = x128*x164 - x146*x160;
    double x167 = S11*x163 + S12*x165 + S13*x166;
    double x168 = S12*x163 + S22*x165 + S23*x166;
    double x169 = S13*x163 + S23*x165 + S33*x166;
    double x170 = -x128*x168 - x129*x169 - x145*x167;
    double x171 = x109*x154 + x111*x158 + x156*x61;
    double x172 = x114*x154 + x124*x156 + x126*x158;
    double x173 = -x113*x160 + x131*x162;
    double x174 = x113*x164 - x144*x162;
    double x175 = -x131*x164 + x144*x160;
    double x176 = S11*x173 + S12*x174 + S13*x175;
    double x177 = S12*x173 + S22*x174 + S23*x175;
    double x178 = S13*x173 + S23*x174 + S33*x175;
    double x179 = t*x5;
    double x180 = -x179;
    double x181 = x180 + x63;
    double x182 = x181 + x22 + x42;
    double *x183 = V;
    double x184 = x183[0];
    double x185 = x179 + x64;
    double x186 = x185 + x22 + x45;
    double x187 = x183[12];
    double x188 = x185 + x24 + x42;
    double x189 = x183[9];
    double x190 = x181 + x24 + x45;
    double x191 = x183[21];
    double x192 = x179 + x63;
    double x193 = x192 + x27 + x36;
    double x194 = x183[3];
    double x195 = x180 + x64;
    double x196 = x195 + x27 + x39;
    double x197 = x183[15];
    double x198 = x195 + x29 + x36;
    double x199 = x183[6];
    double x200 = x192 + x29 + x39;
    double x201 = x183[18];
    double x202 = x182*x184 + x186*x187 + x188*x189 + x190*x191 + x193*x194 + x196*x197 + x198*x199 + x200*x201;
    double x203 = x183[1];
    double x204 = x183[13];
    double x205 = x183[10];
    double x206 = x183[22];
    double x207 = x183[4];
    double x208 = x183[16];
    double x209 = x183[7];
    double x210 = x183[19];
    double x211 = x182*x203 + x186*x204 + x188*x205 + x190*x206 + x193*x207 + x196*x208 + x198*x209 + x200*x210;
    double x212 = x183[2];
    double x213 = x183[14];
    double x214 = x183[11];
    double x215 = x183[23];
    double x216 = x183[5];
    double x217 = x183[17];
    double x218 = x183[8];
    double x219 = x183[20];
    double x220 = x182*x212 + x186*x213 + x188*x214 + x190*x215 + x193*x216 + x196*x217 + x198*x218 + x200*x219;
    double x221 = 1.0*PENER + 1.0*SENER - 0.5*rho*(x202*x202 + x211*x211 + x220*x220);
    double x222 = -x128*x177 - x129*x178 - x145*x176 + x221;
    double *x223 = A;
    double x224 = x182*x223[0] + x186*x223[12] + x188*x223[9] + x190*x223[21] + x193*x223[3] + x196*x223[15] + x198*x223[6] + x200*x223[18];
    double x225 = x182*x223[1] + x186*x223[13] + x188*x223[10] + x190*x223[22] + x193*x223[4] + x196*x223[16] + x198*x223[7] + x200*x223[19];
    double x226 = x182*x223[2] + x186*x223[14] + x188*x223[11] + x190*x223[23] + x193*x223[5] + x196*x223[17] + x198*x223[8] + x200*x223[20];
    double x227 = x184*x46 + x187*x43 + x189*x48 + x191*x50 + x194*x40 + x197*x37 + x199*x52 + x201*x54;
    double x228 = x114*x98;
    double x229 = x184*x76 + x187*x73 + x189*x69 + x191*x66 + x194*x79 + x197*x81 + x199*x84 + x201*x86;
    double x230 = x124*x98;
    double x231 = x12*x194 + x16*x189 + x184*x19 + x187*x22 + x191*x24 + x197*x27 + x199*x8 + x201*x29;
    double x232 = x126*x98;
    double x233 = x203*x46 + x204*x43 + x205*x48 + x206*x50 + x207*x40 + x208*x37 + x209*x52 + x210*x54;
    double x234 = x203*x76 + x204*x73 + x205*x69 + x206*x66 + x207*x79 + x208*x81 + x209*x84 + x210*x86;
    double x235 = x12*x207 + x16*x205 + x19*x203 + x204*x22 + x206*x24 + x208*x27 + x209*x8 + x210*x29;
    double x236 = x212*x46 + x213*x43 + x214*x48 + x215*x50 + x216*x40 + x217*x37 + x218*x52 + x219*x54;
    double x237 = x212*x76 + x213*x73 + x214*x69 + x215*x66 + x216*x79 + x217*x81 + x218*x84 + x219*x86;
    double x238 = x12*x216 + x16*x214 + x19*x212 + x213*x22 + x215*x24 + x217*x27 + x218*x8 + x219*x29;
    double x239 = x128*x225 + x129*x226 + x145*x224 + x202*(x227*x228 + x229*x230 + x231*x232) + x211*(x228*x233 + x230*x234 + x232*x235) + x220*(x228*x236 + x230*x237 + x232*x238);
    double x240 = rho*x182;
    double x241 = -x113*x151 - x130*x150 - x144*x149;
    double x242 = -x113*x178 - x130*x177 - x144*x176;
    double x243 = -x113*x169 - x130*x168 - x144*x167 + x221;
    double x244 = x61*x98;
    double x245 = x109*x98;
    double x246 = x111*x98;
    double x247 = x113*x226 + x130*x225 + x144*x224 + x202*(x227*x245 + x229*x244 + x231*x246) + x211*(x233*x245 + x234*x244 + x235*x246) + x220*(x236*x245 + x237*x244 + x238*x246);
    double x248 = -x160*x168 - x161*x169 - x164*x167;
    double x249 = -x160*x177 - x161*x178 - x164*x176;
    double x250 = -x149*x164 - x150*x160 - x151*x161 + x221;
    double x251 = x153*x98;
    double x252 = x155*x98;
    double x253 = x157*x98;
    double x254 = x160*x225 + x161*x226 + x164*x224 + x202*(x227*x251 + x229*x252 + x231*x253) + x211*(x233*x251 + x234*x252 + x235*x253) + x220*(x236*x251 + x237*x252 + x238*x253);
    double x255 = x12*x253 + x251*x40 + x252*x79;
    double x256 = x12*x246 + x244*x79 + x245*x40;
    double x257 = x12*x232 + x228*x40 + x230*x79;
    double x258 = rho*x193;
    double x259 = x251*x52 + x252*x84 + x253*x8;
    double x260 = x244*x84 + x245*x52 + x246*x8;
    double x261 = x228*x52 + x230*x84 + x232*x8;
    double x262 = rho*x198;
    double x263 = x16*x253 + x251*x48 + x252*x69;
    double x264 = x16*x246 + x244*x69 + x245*x48;
    double x265 = x16*x232 + x228*x48 + x230*x69;
    double x266 = rho*x188;
    double x267 = x22*x253 + x251*x43 + x252*x73;
    double x268 = x22*x246 + x244*x73 + x245*x43;
    double x269 = x22*x232 + x228*x43 + x230*x73;
    double x270 = rho*x186;
    double x271 = x251*x37 + x252*x81 + x253*x27;
    double x272 = x244*x81 + x245*x37 + x246*x27;
    double x273 = x228*x37 + x230*x81 + x232*x27;
    double x274 = rho*x196;
    double x275 = x251*x54 + x252*x86 + x253*x29;
    double x276 = x244*x86 + x245*x54 + x246*x29;
    double x277 = x228*x54 + x230*x86 + x232*x29;
    double x278 = rho*x200;
    double x279 = x24*x253 + x251*x50 + x252*x66;
    double x280 = x24*x246 + x244*x66 + x245*x50;
    double x281 = x228*x50 + x230*x66 + x232*x24;
    double x282 = rho*x190;
    
    res_0[0] = x97*(x152*x159 + x170*x171 + x172*x222 - x239*x240);
    res_0[1] = x97*(x159*x241 + x171*x243 + x172*x242 - x240*x247);
    res_0[2] = x97*(x159*x250 + x171*x248 + x172*x249 - x240*x254);
    res_0[3] = x97*(x152*x255 + x170*x256 + x222*x257 - x239*x258);
    res_0[4] = x97*(x241*x255 + x242*x257 + x243*x256 - x247*x258);
    res_0[5] = x97*(x248*x256 + x249*x257 + x250*x255 - x254*x258);
    res_0[6] = x97*(x152*x259 + x170*x260 + x222*x261 - x239*x262);
    res_0[7] = x97*(x241*x259 + x242*x261 + x243*x260 - x247*x262);
    res_0[8] = x97*(x248*x260 + x249*x261 + x250*x259 - x254*x262);
    res_0[9] = x97*(x152*x263 + x170*x264 + x222*x265 - x239*x266);
    res_0[10] = x97*(x241*x263 + x242*x265 + x243*x264 - x247*x266);
    res_0[11] = x97*(x248*x264 + x249*x265 + x250*x263 - x254*x266);
    res_0[12] = x97*(x152*x267 + x170*x268 + x222*x269 - x239*x270);
    res_0[13] = x97*(x241*x267 + x242*x269 + x243*x268 - x247*x270);
    res_0[14] = x97*(x248*x268 + x249*x269 + x250*x267 - x254*x270);
    res_0[15] = x97*(x152*x271 + x170*x272 + x222*x273 - x239*x274);
    res_0[16] = x97*(x241*x271 + x242*x273 + x243*x272 - x247*x274);
    res_0[17] = x97*(x248*x272 + x249*x273 + x250*x271 - x254*x274);
    res_0[18] = x97*(x152*x275 + x170*x276 + x222*x277 - x239*x278);
    res_0[19] = x97*(x241*x275 + x242*x277 + x243*x276 - x247*x278);
    res_0[20] = x97*(x248*x276 + x249*x277 + x250*x275 - x254*x278);
    res_0[21] = x97*(x152*x279 + x170*x280 + x222*x281 - x239*x282);
    res_0[22] = x97*(x241*x279 + x242*x281 + x243*x280 - x247*x282);
    res_0[23] = x97*(x248*x280 + x249*x281 + x250*x279 - x254*x282);
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
    
    double x0 = (1.0/8.0)*r;
    double x1 = -x0;
    double x2 = x1 - 1.0/8.0;
    double x3 = (1.0/8.0)*s;
    double x4 = -x3;
    double x5 = s*x0;
    double x6 = -x5;
    double x7 = x4 + x6;
    double x8 = x2 + x7;
    double *x9 = coord;
    double x10 = x9[6];
    double x11 = x3 + x5;
    double x12 = x11 + x2;
    double x13 = x9[3];
    double x14 = x0 - 1.0/8.0;
    double x15 = x4 + x5;
    double x16 = x14 + x15;
    double x17 = x9[9];
    double x18 = x3 + x6;
    double x19 = x14 + x18;
    double x20 = x9[0];
    double x21 = x1 + 1.0/8.0;
    double x22 = x15 + x21;
    double x23 = x9[12];
    double x24 = x18 + x21;
    double x25 = x9[21];
    double x26 = x0 + 1.0/8.0;
    double x27 = x26 + x7;
    double x28 = x9[15];
    double x29 = x11 + x26;
    double x30 = x9[18];
    double x31 = x10*x8 + x12*x13 + x16*x17 + x19*x20 + x22*x23 + x24*x25 + x27*x28 + x29*x30;
    double x32 = (1.0/8.0)*t;
    double x33 = -x32;
    double x34 = t*x0;
    double x35 = -x34;
    double x36 = x33 + x35;
    double x37 = x2 + x36;
    double x38 = x9[17];
    double x39 = x32 + x34;
    double x40 = x2 + x39;
    double x41 = x9[5];
    double x42 = x33 + x34;
    double x43 = x14 + x42;
    double x44 = x9[14];
    double x45 = x32 + x35;
    double x46 = x14 + x45;
    double x47 = x9[2];
    double x48 = x21 + x42;
    double x49 = x9[11];
    double x50 = x21 + x45;
    double x51 = x9[23];
    double x52 = x26 + x36;
    double x53 = x9[8];
    double x54 = x26 + x39;
    double x55 = x9[20];
    double x56 = x37*x38 + x40*x41 + x43*x44 + x46*x47 + x48*x49 + x50*x51 + x52*x53 + x54*x55;
    double x57 = x31*x56;
    double x58 = x12*x41 + x16*x49 + x19*x47 + x22*x44 + x24*x51 + x27*x38 + x29*x55 + x53*x8;
    double x59 = x10*x52 + x13*x40 + x17*x48 + x20*x46 + x23*x43 + x25*x50 + x28*x37 + x30*x54;
    double x60 = x58*x59;
    double x61 = x57 - x60;
    double x62 = x4 - 1.0/8.0;
    double x63 = t*x3;
    double x64 = -x63;
    double x65 = x33 + x64;
    double x66 = x62 + x65;
    double x67 = x9[22];
    double x68 = x32 + x63;
    double x69 = x62 + x68;
    double x70 = x9[10];
    double x71 = x3 - 1.0/8.0;
    double x72 = x33 + x63;
    double x73 = x71 + x72;
    double x74 = x9[13];
    double x75 = x32 + x64;
    double x76 = x71 + x75;
    double x77 = x9[1];
    double x78 = x4 + 1.0/8.0;
    double x79 = x72 + x78;
    double x80 = x9[4];
    double x81 = x75 + x78;
    double x82 = x9[16];
    double x83 = x3 + 1.0/8.0;
    double x84 = x65 + x83;
    double x85 = x9[7];
    double x86 = x68 + x83;
    double x87 = x9[19];
    double x88 = x66*x67 + x69*x70 + x73*x74 + x76*x77 + x79*x80 + x81*x82 + x84*x85 + x86*x87;
    double x89 = x12*x80 + x16*x70 + x19*x77 + x22*x74 + x24*x67 + x27*x82 + x29*x87 + x8*x85;
    double x90 = x38*x81 + x41*x79 + x44*x73 + x47*x76 + x49*x69 + x51*x66 + x53*x84 + x55*x86;
    double x91 = x59*x90;
    double x92 = x37*x82 + x40*x80 + x43*x74 + x46*x77 + x48*x70 + x50*x67 + x52*x85 + x54*x87;
    double x93 = x10*x84 + x13*x79 + x17*x69 + x20*x76 + x23*x73 + x25*x66 + x28*x81 + x30*x86;
    double x94 = x58*x93;
    double x95 = x31*x90;
    double x96 = x56*x93;
    double x97 = x57*x88 - x60*x88 + x89*x91 - x89*x96 + x92*x94 - x92*x95;
    double x98 = 1.0/x97;
    double *x99 = U;
    double x100 = x99[23];
    double x101 = x99[11];
    double x102 = x99[14];
    double x103 = x99[2];
    double x104 = x99[5];
    double x105 = x99[17];
    double x106 = x99[8];
    double x107 = x99[20];
    double x108 = x98*(x100*x66 + x101*x69 + x102*x73 + x103*x76 + x104*x79 + x105*x81 + x106*x84 + x107*x86);
    double x109 = x94 - x95;
    double x110 = x98*(x100*x50 + x101*x48 + x102*x43 + x103*x46 + x104*x40 + x105*x37 + x106*x52 + x107*x54);
    double x111 = x91 - x96;
    double x112 = x98*(x100*x24 + x101*x16 + x102*x22 + x103*x19 + x104*x12 + x105*x27 + x106*x8 + x107*x29);
    double x113 = x108*x61 + x109*x110 + x111*x112;
    double x114 = -x58*x88 + x89*x90;
    double x115 = x99[16];
    double x116 = x99[4];
    double x117 = x99[13];
    double x118 = x99[1];
    double x119 = x99[10];
    double x120 = x99[22];
    double x121 = x99[7];
    double x122 = x99[19];
    double x123 = x98*(x115*x37 + x116*x40 + x117*x43 + x118*x46 + x119*x48 + x120*x50 + x121*x52 + x122*x54);
    double x124 = -x56*x89 + x58*x92;
    double x125 = x98*(x115*x81 + x116*x79 + x117*x73 + x118*x76 + x119*x69 + x120*x66 + x121*x84 + x122*x86);
    double x126 = x56*x88 - x90*x92;
    double x127 = x98*(x115*x27 + x116*x12 + x117*x22 + x118*x19 + x119*x16 + x120*x24 + x121*x8 + x122*x29);
    double x128 = x114*x123 + x124*x125 + x126*x127;
    double x129 = x108*x124 + x110*x114 + x112*x126;
    double x130 = x109*x123 + x111*x127 + x125*x61 + 1.0;
    double x131 = x113*x128 - x129*x130;
    double x132 = x99[21];
    double x133 = x99[9];
    double x134 = x99[12];
    double x135 = x99[0];
    double x136 = x99[3];
    double x137 = x99[15];
    double x138 = x99[6];
    double x139 = x99[18];
    double x140 = x98*(x132*x66 + x133*x69 + x134*x73 + x135*x76 + x136*x79 + x137*x81 + x138*x84 + x139*x86);
    double x141 = x98*(x132*x50 + x133*x48 + x134*x43 + x135*x46 + x136*x40 + x137*x37 + x138*x52 + x139*x54);
    double x142 = x98*(x12*x136 + x132*x24 + x133*x16 + x134*x22 + x135*x19 + x137*x27 + x138*x8 + x139*x29);
    double x143 = x109*x141 + x111*x142 + x140*x61;
    double x144 = x114*x141 + x124*x140 + x126*x142 + 1.0;
    double x145 = -x113*x144 + x129*x143;
    double x146 = -x128*x143 + x130*x144;
    double x147 = S11*x131 + S12*x145 + S13*x146;
    double x148 = S12*x131 + S22*x145 + S23*x146;
    double x149 = S13*x131 + S23*x145 + S33*x146;
    double x150 = -x128*x148 - x129*x149 - x144*x147;
    double x151 = x31*x88 - x89*x93;
    double x152 = x46*x98;
    double x153 = -x31*x92 + x59*x89;
    double x154 = x76*x98;
    double x155 = -x59*x88 + x92*x93;
    double x156 = x19*x98;
    double x157 = x151*x152 + x153*x154 + x155*x156;
    double x158 = x123*x151 + x125*x153 + x127*x155;
    double x159 = x108*x153 + x110*x151 + x112*x155 + 1.0;
    double x160 = -x128*x159 + x129*x158;
    double x161 = x140*x153 + x141*x151 + x142*x155;
    double x162 = -x129*x161 + x144*x159;
    double x163 = x128*x161 - x144*x158;
    double x164 = S11*x160 + S12*x162 + S13*x163;
    double x165 = S12*x160 + S22*x162 + S23*x163;
    double x166 = S13*x160 + S23*x162 + S33*x163;
    double x167 = -x128*x165 - x129*x166 - x144*x164;
    double x168 = x109*x152 + x111*x156 + x154*x61;
    double x169 = x114*x152 + x124*x154 + x126*x156;
    double x170 = -x113*x158 + x130*x159;
    double x171 = x113*x161 - x143*x159;
    double x172 = -x130*x161 + x143*x158;
    double x173 = S11*x170 + S12*x171 + S13*x172;
    double x174 = S12*x170 + S22*x171 + S23*x172;
    double x175 = S13*x170 + S23*x171 + S33*x172;
    double x176 = 1.0*PENER + 1.0*SENER;
    double x177 = -x128*x174 - x129*x175 - x144*x173 + x176;
    double x178 = -x113*x149 - x130*x148 - x143*x147;
    double x179 = -x113*x175 - x130*x174 - x143*x173;
    double x180 = -x113*x166 - x130*x165 - x143*x164 + x176;
    double x181 = -x158*x165 - x159*x166 - x161*x164;
    double x182 = -x158*x174 - x159*x175 - x161*x173;
    double x183 = -x147*x161 - x148*x158 - x149*x159 + x176;
    double x184 = x40*x98;
    double x185 = x79*x98;
    double x186 = x12*x98;
    double x187 = x151*x184 + x153*x185 + x155*x186;
    double x188 = x109*x184 + x111*x186 + x185*x61;
    double x189 = x114*x184 + x124*x185 + x126*x186;
    double x190 = x52*x98;
    double x191 = x84*x98;
    double x192 = x8*x98;
    double x193 = x151*x190 + x153*x191 + x155*x192;
    double x194 = x109*x190 + x111*x192 + x191*x61;
    double x195 = x114*x190 + x124*x191 + x126*x192;
    double x196 = x48*x98;
    double x197 = x69*x98;
    double x198 = x16*x98;
    double x199 = x151*x196 + x153*x197 + x155*x198;
    double x200 = x109*x196 + x111*x198 + x197*x61;
    double x201 = x114*x196 + x124*x197 + x126*x198;
    double x202 = x43*x98;
    double x203 = x73*x98;
    double x204 = x22*x98;
    double x205 = x151*x202 + x153*x203 + x155*x204;
    double x206 = x109*x202 + x111*x204 + x203*x61;
    double x207 = x114*x202 + x124*x203 + x126*x204;
    double x208 = x37*x98;
    double x209 = x81*x98;
    double x210 = x27*x98;
    double x211 = x151*x208 + x153*x209 + x155*x210;
    double x212 = x109*x208 + x111*x210 + x209*x61;
    double x213 = x114*x208 + x124*x209 + x126*x210;
    double x214 = x54*x98;
    double x215 = x86*x98;
    double x216 = x29*x98;
    double x217 = x151*x214 + x153*x215 + x155*x216;
    double x218 = x109*x214 + x111*x216 + x215*x61;
    double x219 = x114*x214 + x124*x215 + x126*x216;
    double x220 = x50*x98;
    double x221 = x66*x98;
    double x222 = x24*x98;
    double x223 = x151*x220 + x153*x221 + x155*x222;
    double x224 = x109*x220 + x111*x222 + x221*x61;
    double x225 = x114*x220 + x124*x221 + x126*x222;
    
    res_0[0] = x97*(x150*x157 + x167*x168 + x169*x177);
    res_0[1] = x97*(x157*x178 + x168*x180 + x169*x179);
    res_0[2] = x97*(x157*x183 + x168*x181 + x169*x182);
    res_0[3] = x97*(x150*x187 + x167*x188 + x177*x189);
    res_0[4] = x97*(x178*x187 + x179*x189 + x180*x188);
    res_0[5] = x97*(x181*x188 + x182*x189 + x183*x187);
    res_0[6] = x97*(x150*x193 + x167*x194 + x177*x195);
    res_0[7] = x97*(x178*x193 + x179*x195 + x180*x194);
    res_0[8] = x97*(x181*x194 + x182*x195 + x183*x193);
    res_0[9] = x97*(x150*x199 + x167*x200 + x177*x201);
    res_0[10] = x97*(x178*x199 + x179*x201 + x180*x200);
    res_0[11] = x97*(x181*x200 + x182*x201 + x183*x199);
    res_0[12] = x97*(x150*x205 + x167*x206 + x177*x207);
    res_0[13] = x97*(x178*x205 + x179*x207 + x180*x206);
    res_0[14] = x97*(x181*x206 + x182*x207 + x183*x205);
    res_0[15] = x97*(x150*x211 + x167*x212 + x177*x213);
    res_0[16] = x97*(x178*x211 + x179*x213 + x180*x212);
    res_0[17] = x97*(x181*x212 + x182*x213 + x183*x211);
    res_0[18] = x97*(x150*x217 + x167*x218 + x177*x219);
    res_0[19] = x97*(x178*x217 + x179*x219 + x180*x218);
    res_0[20] = x97*(x181*x218 + x182*x219 + x183*x217);
    res_0[21] = x97*(x150*x223 + x167*x224 + x177*x225);
    res_0[22] = x97*(x178*x223 + x179*x225 + x180*x224);
    res_0[23] = x97*(x181*x224 + x182*x225 + x183*x223);
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
    
    double x0 = (1.0/8.0)*r;
    double x1 = -x0;
    double x2 = x1 - 1.0/8.0;
    double x3 = (1.0/8.0)*s;
    double x4 = -x3;
    double x5 = s*x0;
    double x6 = -x5;
    double x7 = x4 + x6;
    double x8 = x2 + x7;
    double *x9 = coord;
    double x10 = x9[6];
    double x11 = x3 + x5;
    double x12 = x11 + x2;
    double x13 = x9[3];
    double x14 = x0 - 1.0/8.0;
    double x15 = x4 + x5;
    double x16 = x14 + x15;
    double x17 = x9[9];
    double x18 = x3 + x6;
    double x19 = x14 + x18;
    double x20 = x9[0];
    double x21 = x1 + 1.0/8.0;
    double x22 = x15 + x21;
    double x23 = x9[12];
    double x24 = x18 + x21;
    double x25 = x9[21];
    double x26 = x0 + 1.0/8.0;
    double x27 = x26 + x7;
    double x28 = x9[15];
    double x29 = x11 + x26;
    double x30 = x9[18];
    double x31 = x10*x8 + x12*x13 + x16*x17 + x19*x20 + x22*x23 + x24*x25 + x27*x28 + x29*x30;
    double x32 = (1.0/8.0)*t;
    double x33 = -x32;
    double x34 = t*x0;
    double x35 = -x34;
    double x36 = x33 + x35;
    double x37 = x2 + x36;
    double x38 = x9[17];
    double x39 = x32 + x34;
    double x40 = x2 + x39;
    double x41 = x9[5];
    double x42 = x33 + x34;
    double x43 = x14 + x42;
    double x44 = x9[14];
    double x45 = x32 + x35;
    double x46 = x14 + x45;
    double x47 = x9[2];
    double x48 = x21 + x42;
    double x49 = x9[11];
    double x50 = x21 + x45;
    double x51 = x9[23];
    double x52 = x26 + x36;
    double x53 = x9[8];
    double x54 = x26 + x39;
    double x55 = x9[20];
    double x56 = x37*x38 + x40*x41 + x43*x44 + x46*x47 + x48*x49 + x50*x51 + x52*x53 + x54*x55;
    double x57 = x31*x56;
    double x58 = x12*x41 + x16*x49 + x19*x47 + x22*x44 + x24*x51 + x27*x38 + x29*x55 + x53*x8;
    double x59 = x10*x52 + x13*x40 + x17*x48 + x20*x46 + x23*x43 + x25*x50 + x28*x37 + x30*x54;
    double x60 = x58*x59;
    double x61 = x57 - x60;
    double x62 = x4 - 1.0/8.0;
    double x63 = t*x3;
    double x64 = -x63;
    double x65 = x33 + x64;
    double x66 = x62 + x65;
    double x67 = x9[22];
    double x68 = x32 + x63;
    double x69 = x62 + x68;
    double x70 = x9[10];
    double x71 = x3 - 1.0/8.0;
    double x72 = x33 + x63;
    double x73 = x71 + x72;
    double x74 = x9[13];
    double x75 = x32 + x64;
    double x76 = x71 + x75;
    double x77 = x9[1];
    double x78 = x4 + 1.0/8.0;
    double x79 = x72 + x78;
    double x80 = x9[4];
    double x81 = x75 + x78;
    double x82 = x9[16];
    double x83 = x3 + 1.0/8.0;
    double x84 = x65 + x83;
    double x85 = x9[7];
    double x86 = x68 + x83;
    double x87 = x9[19];
    double x88 = x66*x67 + x69*x70 + x73*x74 + x76*x77 + x79*x80 + x81*x82 + x84*x85 + x86*x87;
    double x89 = x12*x80 + x16*x70 + x19*x77 + x22*x74 + x24*x67 + x27*x82 + x29*x87 + x8*x85;
    double x90 = x38*x81 + x41*x79 + x44*x73 + x47*x76 + x49*x69 + x51*x66 + x53*x84 + x55*x86;
    double x91 = x59*x90;
    double x92 = x37*x82 + x40*x80 + x43*x74 + x46*x77 + x48*x70 + x50*x67 + x52*x85 + x54*x87;
    double x93 = x10*x84 + x13*x79 + x17*x69 + x20*x76 + x23*x73 + x25*x66 + x28*x81 + x30*x86;
    double x94 = x58*x93;
    double x95 = x31*x90;
    double x96 = x56*x93;
    double x97 = x57*x88 - x60*x88 + x89*x91 - x89*x96 + x92*x94 - x92*x95;
    double x98 = 1.0/x97;
    double *x99 = U;
    double x100 = x99[23];
    double x101 = x99[11];
    double x102 = x99[14];
    double x103 = x99[2];
    double x104 = x99[5];
    double x105 = x99[17];
    double x106 = x99[8];
    double x107 = x99[20];
    double x108 = x98*(x100*x66 + x101*x69 + x102*x73 + x103*x76 + x104*x79 + x105*x81 + x106*x84 + x107*x86);
    double x109 = x94 - x95;
    double x110 = x98*(x100*x50 + x101*x48 + x102*x43 + x103*x46 + x104*x40 + x105*x37 + x106*x52 + x107*x54);
    double x111 = x91 - x96;
    double x112 = x98*(x100*x24 + x101*x16 + x102*x22 + x103*x19 + x104*x12 + x105*x27 + x106*x8 + x107*x29);
    double x113 = x108*x61 + x109*x110 + x111*x112;
    double x114 = -x58*x88 + x89*x90;
    double x115 = x99[16];
    double x116 = x99[4];
    double x117 = x99[13];
    double x118 = x99[1];
    double x119 = x99[10];
    double x120 = x99[22];
    double x121 = x99[7];
    double x122 = x99[19];
    double x123 = x98*(x115*x37 + x116*x40 + x117*x43 + x118*x46 + x119*x48 + x120*x50 + x121*x52 + x122*x54);
    double x124 = -x56*x89 + x58*x92;
    double x125 = x98*(x115*x81 + x116*x79 + x117*x73 + x118*x76 + x119*x69 + x120*x66 + x121*x84 + x122*x86);
    double x126 = x56*x88 - x90*x92;
    double x127 = x98*(x115*x27 + x116*x12 + x117*x22 + x118*x19 + x119*x16 + x120*x24 + x121*x8 + x122*x29);
    double x128 = x114*x123 + x124*x125 + x126*x127;
    double x129 = x108*x124 + x110*x114 + x112*x126;
    double x130 = x109*x123 + x111*x127 + x125*x61;
    double x131 = x130 + 1.0;
    double x132 = x113*x128 - x129*x131;
    double x133 = x99[21];
    double x134 = x99[9];
    double x135 = x99[12];
    double x136 = x99[0];
    double x137 = x99[3];
    double x138 = x99[15];
    double x139 = x99[6];
    double x140 = x99[18];
    double x141 = x98*(x133*x66 + x134*x69 + x135*x73 + x136*x76 + x137*x79 + x138*x81 + x139*x84 + x140*x86);
    double x142 = x98*(x133*x50 + x134*x48 + x135*x43 + x136*x46 + x137*x40 + x138*x37 + x139*x52 + x140*x54);
    double x143 = x98*(x12*x137 + x133*x24 + x134*x16 + x135*x22 + x136*x19 + x138*x27 + x139*x8 + x140*x29);
    double x144 = x109*x142 + x111*x143 + x141*x61;
    double x145 = x114*x142 + x124*x141 + x126*x143;
    double x146 = x145 + 1.0;
    double x147 = -x113*x146 + x129*x144;
    double x148 = -x128*x144 + x131*x146;
    double x149 = S11*x132 + S12*x147 + S13*x148;
    double x150 = S12*x132 + S22*x147 + S23*x148;
    double x151 = S13*x132 + S23*x147 + S33*x148;
    double x152 = -x128*x150 - x129*x151 - x145*x149;
    double x153 = x31*x88 - x89*x93;
    double x154 = x46*x98;
    double x155 = -x31*x92 + x59*x89;
    double x156 = x76*x98;
    double x157 = -x59*x88 + x92*x93;
    double x158 = x19*x98;
    double x159 = x153*x154 + x155*x156 + x157*x158;
    double x160 = x123*x153 + x125*x155 + x127*x157;
    double x161 = x108*x155 + x110*x153 + x112*x157;
    double x162 = x161 + 1.0;
    double x163 = -x128*x162 + x129*x160;
    double x164 = x141*x155 + x142*x153 + x143*x157;
    double x165 = -x129*x164 + x146*x162;
    double x166 = x128*x164 - x146*x160;
    double x167 = S11*x163 + S12*x165 + S13*x166;
    double x168 = S12*x163 + S22*x165 + S23*x166;
    double x169 = S13*x163 + S23*x165 + S33*x166;
    double x170 = -x128*x168 - x129*x169 - x145*x167;
    double x171 = x109*x154 + x111*x158 + x156*x61;
    double x172 = x114*x154 + x124*x156 + x126*x158;
    double x173 = -x113*x160 + x131*x162;
    double x174 = x113*x164 - x144*x162;
    double x175 = -x131*x164 + x144*x160;
    double x176 = S11*x173 + S12*x174 + S13*x175;
    double x177 = S12*x173 + S22*x174 + S23*x175;
    double x178 = S13*x173 + S23*x174 + S33*x175;
    double x179 = 1.0*PENER + 1.0*SENER;
    double x180 = -x128*x177 - x129*x178 - x145*x176 + x179;
    double x181 = -x113*x151 - x130*x150 - x144*x149;
    double x182 = -x113*x178 - x130*x177 - x144*x176;
    double x183 = -x113*x169 - x130*x168 - x144*x167 + x179;
    double x184 = -x160*x168 - x161*x169 - x164*x167;
    double x185 = -x160*x177 - x161*x178 - x164*x176;
    double x186 = -x149*x164 - x150*x160 - x151*x161 + x179;
    double x187 = x40*x98;
    double x188 = x79*x98;
    double x189 = x12*x98;
    double x190 = x153*x187 + x155*x188 + x157*x189;
    double x191 = x109*x187 + x111*x189 + x188*x61;
    double x192 = x114*x187 + x124*x188 + x126*x189;
    double x193 = x52*x98;
    double x194 = x84*x98;
    double x195 = x8*x98;
    double x196 = x153*x193 + x155*x194 + x157*x195;
    double x197 = x109*x193 + x111*x195 + x194*x61;
    double x198 = x114*x193 + x124*x194 + x126*x195;
    double x199 = x48*x98;
    double x200 = x69*x98;
    double x201 = x16*x98;
    double x202 = x153*x199 + x155*x200 + x157*x201;
    double x203 = x109*x199 + x111*x201 + x200*x61;
    double x204 = x114*x199 + x124*x200 + x126*x201;
    double x205 = x43*x98;
    double x206 = x73*x98;
    double x207 = x22*x98;
    double x208 = x153*x205 + x155*x206 + x157*x207;
    double x209 = x109*x205 + x111*x207 + x206*x61;
    double x210 = x114*x205 + x124*x206 + x126*x207;
    double x211 = x37*x98;
    double x212 = x81*x98;
    double x213 = x27*x98;
    double x214 = x153*x211 + x155*x212 + x157*x213;
    double x215 = x109*x211 + x111*x213 + x212*x61;
    double x216 = x114*x211 + x124*x212 + x126*x213;
    double x217 = x54*x98;
    double x218 = x86*x98;
    double x219 = x29*x98;
    double x220 = x153*x217 + x155*x218 + x157*x219;
    double x221 = x109*x217 + x111*x219 + x218*x61;
    double x222 = x114*x217 + x124*x218 + x126*x219;
    double x223 = x50*x98;
    double x224 = x66*x98;
    double x225 = x24*x98;
    double x226 = x153*x223 + x155*x224 + x157*x225;
    double x227 = x109*x223 + x111*x225 + x224*x61;
    double x228 = x114*x223 + x124*x224 + x126*x225;
    
    res_0[0] = x97*(x152*x159 + x170*x171 + x172*x180);
    res_0[1] = x97*(x159*x181 + x171*x183 + x172*x182);
    res_0[2] = x97*(x159*x186 + x171*x184 + x172*x185);
    res_0[3] = x97*(x152*x190 + x170*x191 + x180*x192);
    res_0[4] = x97*(x181*x190 + x182*x192 + x183*x191);
    res_0[5] = x97*(x184*x191 + x185*x192 + x186*x190);
    res_0[6] = x97*(x152*x196 + x170*x197 + x180*x198);
    res_0[7] = x97*(x181*x196 + x182*x198 + x183*x197);
    res_0[8] = x97*(x184*x197 + x185*x198 + x186*x196);
    res_0[9] = x97*(x152*x202 + x170*x203 + x180*x204);
    res_0[10] = x97*(x181*x202 + x182*x204 + x183*x203);
    res_0[11] = x97*(x184*x203 + x185*x204 + x186*x202);
    res_0[12] = x97*(x152*x208 + x170*x209 + x180*x210);
    res_0[13] = x97*(x181*x208 + x182*x210 + x183*x209);
    res_0[14] = x97*(x184*x209 + x185*x210 + x186*x208);
    res_0[15] = x97*(x152*x214 + x170*x215 + x180*x216);
    res_0[16] = x97*(x181*x214 + x182*x216 + x183*x215);
    res_0[17] = x97*(x184*x215 + x185*x216 + x186*x214);
    res_0[18] = x97*(x152*x220 + x170*x221 + x180*x222);
    res_0[19] = x97*(x181*x220 + x182*x222 + x183*x221);
    res_0[20] = x97*(x184*x221 + x185*x222 + x186*x220);
    res_0[21] = x97*(x152*x226 + x170*x227 + x180*x228);
    res_0[22] = x97*(x181*x226 + x182*x228 + x183*x227);
    res_0[23] = x97*(x184*x227 + x185*x228 + x186*x226);
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
    
    double x0 = (1.0/8.0)*r;
    double x1 = -x0;
    double x2 = x1 - 1.0/8.0;
    double x3 = (1.0/8.0)*s;
    double x4 = -x3;
    double x5 = s*x0;
    double x6 = -x5;
    double x7 = x4 + x6;
    double x8 = x2 + x7;
    double *x9 = coord;
    double x10 = x9[6];
    double x11 = x3 + x5;
    double x12 = x11 + x2;
    double x13 = x9[3];
    double x14 = x0 - 1.0/8.0;
    double x15 = x4 + x5;
    double x16 = x14 + x15;
    double x17 = x9[9];
    double x18 = x3 + x6;
    double x19 = x14 + x18;
    double x20 = x9[0];
    double x21 = x1 + 1.0/8.0;
    double x22 = x15 + x21;
    double x23 = x9[12];
    double x24 = x18 + x21;
    double x25 = x9[21];
    double x26 = x0 + 1.0/8.0;
    double x27 = x26 + x7;
    double x28 = x9[15];
    double x29 = x11 + x26;
    double x30 = x9[18];
    double x31 = x10*x8 + x12*x13 + x16*x17 + x19*x20 + x22*x23 + x24*x25 + x27*x28 + x29*x30;
    double x32 = (1.0/8.0)*t;
    double x33 = -x32;
    double x34 = t*x0;
    double x35 = -x34;
    double x36 = x33 + x35;
    double x37 = x2 + x36;
    double x38 = x9[17];
    double x39 = x32 + x34;
    double x40 = x2 + x39;
    double x41 = x9[5];
    double x42 = x33 + x34;
    double x43 = x14 + x42;
    double x44 = x9[14];
    double x45 = x32 + x35;
    double x46 = x14 + x45;
    double x47 = x9[2];
    double x48 = x21 + x42;
    double x49 = x9[11];
    double x50 = x21 + x45;
    double x51 = x9[23];
    double x52 = x26 + x36;
    double x53 = x9[8];
    double x54 = x26 + x39;
    double x55 = x9[20];
    double x56 = x37*x38 + x40*x41 + x43*x44 + x46*x47 + x48*x49 + x50*x51 + x52*x53 + x54*x55;
    double x57 = x31*x56;
    double x58 = x12*x41 + x16*x49 + x19*x47 + x22*x44 + x24*x51 + x27*x38 + x29*x55 + x53*x8;
    double x59 = x10*x52 + x13*x40 + x17*x48 + x20*x46 + x23*x43 + x25*x50 + x28*x37 + x30*x54;
    double x60 = x58*x59;
    double x61 = x57 - x60;
    double x62 = x4 - 1.0/8.0;
    double x63 = t*x3;
    double x64 = -x63;
    double x65 = x33 + x64;
    double x66 = x62 + x65;
    double x67 = x9[22];
    double x68 = x32 + x63;
    double x69 = x62 + x68;
    double x70 = x9[10];
    double x71 = x3 - 1.0/8.0;
    double x72 = x33 + x63;
    double x73 = x71 + x72;
    double x74 = x9[13];
    double x75 = x32 + x64;
    double x76 = x71 + x75;
    double x77 = x9[1];
    double x78 = x4 + 1.0/8.0;
    double x79 = x72 + x78;
    double x80 = x9[4];
    double x81 = x75 + x78;
    double x82 = x9[16];
    double x83 = x3 + 1.0/8.0;
    double x84 = x65 + x83;
    double x85 = x9[7];
    double x86 = x68 + x83;
    double x87 = x9[19];
    double x88 = x66*x67 + x69*x70 + x73*x74 + x76*x77 + x79*x80 + x81*x82 + x84*x85 + x86*x87;
    double x89 = x12*x80 + x16*x70 + x19*x77 + x22*x74 + x24*x67 + x27*x82 + x29*x87 + x8*x85;
    double x90 = x38*x81 + x41*x79 + x44*x73 + x47*x76 + x49*x69 + x51*x66 + x53*x84 + x55*x86;
    double x91 = x59*x90;
    double x92 = x37*x82 + x40*x80 + x43*x74 + x46*x77 + x48*x70 + x50*x67 + x52*x85 + x54*x87;
    double x93 = x10*x84 + x13*x79 + x17*x69 + x20*x76 + x23*x73 + x25*x66 + x28*x81 + x30*x86;
    double x94 = x58*x93;
    double x95 = x31*x90;
    double x96 = x56*x93;
    double x97 = x57*x88 - x60*x88 + x89*x91 - x89*x96 + x92*x94 - x92*x95;
    double x98 = 1.0/x97;
    double *x99 = U;
    double x100 = x99[23];
    double x101 = x99[11];
    double x102 = x99[14];
    double x103 = x99[2];
    double x104 = x99[5];
    double x105 = x99[17];
    double x106 = x99[8];
    double x107 = x99[20];
    double x108 = x98*(x100*x66 + x101*x69 + x102*x73 + x103*x76 + x104*x79 + x105*x81 + x106*x84 + x107*x86);
    double x109 = x94 - x95;
    double x110 = x98*(x100*x50 + x101*x48 + x102*x43 + x103*x46 + x104*x40 + x105*x37 + x106*x52 + x107*x54);
    double x111 = x91 - x96;
    double x112 = x98*(x100*x24 + x101*x16 + x102*x22 + x103*x19 + x104*x12 + x105*x27 + x106*x8 + x107*x29);
    double x113 = x108*x61 + x109*x110 + x111*x112;
    double x114 = -x58*x88 + x89*x90;
    double x115 = x99[16];
    double x116 = x99[4];
    double x117 = x99[13];
    double x118 = x99[1];
    double x119 = x99[10];
    double x120 = x99[22];
    double x121 = x99[7];
    double x122 = x99[19];
    double x123 = x98*(x115*x37 + x116*x40 + x117*x43 + x118*x46 + x119*x48 + x120*x50 + x121*x52 + x122*x54);
    double x124 = -x56*x89 + x58*x92;
    double x125 = x98*(x115*x81 + x116*x79 + x117*x73 + x118*x76 + x119*x69 + x120*x66 + x121*x84 + x122*x86);
    double x126 = x56*x88 - x90*x92;
    double x127 = x98*(x115*x27 + x116*x12 + x117*x22 + x118*x19 + x119*x16 + x120*x24 + x121*x8 + x122*x29);
    double x128 = x114*x123 + x124*x125 + x126*x127;
    double x129 = x108*x124 + x110*x114 + x112*x126;
    double x130 = x109*x123 + x111*x127 + x125*x61;
    double x131 = x130 + 1.0;
    double x132 = x113*x128 - x129*x131;
    double x133 = x99[21];
    double x134 = x99[9];
    double x135 = x99[12];
    double x136 = x99[0];
    double x137 = x99[3];
    double x138 = x99[15];
    double x139 = x99[6];
    double x140 = x99[18];
    double x141 = x98*(x133*x66 + x134*x69 + x135*x73 + x136*x76 + x137*x79 + x138*x81 + x139*x84 + x140*x86);
    double x142 = x98*(x133*x50 + x134*x48 + x135*x43 + x136*x46 + x137*x40 + x138*x37 + x139*x52 + x140*x54);
    double x143 = x98*(x12*x137 + x133*x24 + x134*x16 + x135*x22 + x136*x19 + x138*x27 + x139*x8 + x140*x29);
    double x144 = x109*x142 + x111*x143 + x141*x61;
    double x145 = x114*x142 + x124*x141 + x126*x143;
    double x146 = x145 + 1.0;
    double x147 = -x113*x146 + x129*x144;
    double x148 = -x128*x144 + x131*x146;
    double x149 = S11*x132 + S12*x147 + S13*x148;
    double x150 = S12*x132 + S22*x147 + S23*x148;
    double x151 = S13*x132 + S23*x147 + S33*x148;
    double x152 = -x128*x150 - x129*x151 - x145*x149;
    double x153 = x31*x88 - x89*x93;
    double x154 = x46*x98;
    double x155 = -x31*x92 + x59*x89;
    double x156 = x76*x98;
    double x157 = -x59*x88 + x92*x93;
    double x158 = x19*x98;
    double x159 = x153*x154 + x155*x156 + x157*x158;
    double x160 = x123*x153 + x125*x155 + x127*x157;
    double x161 = x108*x155 + x110*x153 + x112*x157;
    double x162 = x161 + 1.0;
    double x163 = -x128*x162 + x129*x160;
    double x164 = x141*x155 + x142*x153 + x143*x157;
    double x165 = -x129*x164 + x146*x162;
    double x166 = x128*x164 - x146*x160;
    double x167 = S11*x163 + S12*x165 + S13*x166;
    double x168 = S12*x163 + S22*x165 + S23*x166;
    double x169 = S13*x163 + S23*x165 + S33*x166;
    double x170 = -x128*x168 - x129*x169 - x145*x167;
    double x171 = x109*x154 + x111*x158 + x156*x61;
    double x172 = x114*x154 + x124*x156 + x126*x158;
    double x173 = -x113*x160 + x131*x162;
    double x174 = x113*x164 - x144*x162;
    double x175 = -x131*x164 + x144*x160;
    double x176 = S11*x173 + S12*x174 + S13*x175;
    double x177 = S12*x173 + S22*x174 + S23*x175;
    double x178 = S13*x173 + S23*x174 + S33*x175;
    double x179 = t*x5;
    double x180 = -x179;
    double x181 = x180 + x63;
    double x182 = x181 + x22 + x42;
    double *x183 = V;
    double x184 = x183[0];
    double x185 = x179 + x64;
    double x186 = x185 + x22 + x45;
    double x187 = x183[12];
    double x188 = x185 + x24 + x42;
    double x189 = x183[9];
    double x190 = x181 + x24 + x45;
    double x191 = x183[21];
    double x192 = x179 + x63;
    double x193 = x192 + x27 + x36;
    double x194 = x183[3];
    double x195 = x180 + x64;
    double x196 = x195 + x27 + x39;
    double x197 = x183[15];
    double x198 = x195 + x29 + x36;
    double x199 = x183[6];
    double x200 = x192 + x29 + x39;
    double x201 = x183[18];
    double x202 = x182*x184 + x186*x187 + x188*x189 + x190*x191 + x193*x194 + x196*x197 + x198*x199 + x200*x201;
    double x203 = x183[1];
    double x204 = x183[13];
    double x205 = x183[10];
    double x206 = x183[22];
    double x207 = x183[4];
    double x208 = x183[16];
    double x209 = x183[7];
    double x210 = x183[19];
    double x211 = x182*x203 + x186*x204 + x188*x205 + x190*x206 + x193*x207 + x196*x208 + x198*x209 + x200*x210;
    double x212 = x183[2];
    double x213 = x183[14];
    double x214 = x183[11];
    double x215 = x183[23];
    double x216 = x183[5];
    double x217 = x183[17];
    double x218 = x183[8];
    double x219 = x183[20];
    double x220 = x182*x212 + x186*x213 + x188*x214 + x190*x215 + x193*x216 + x196*x217 + x198*x218 + x200*x219;
    double x221 = 1.0*PENER + 1.0*SENER - 0.5*rho*(x202*x202 + x211*x211 + x220*x220);
    double x222 = -x128*x177 - x129*x178 - x145*x176 + x221;
    double *x223 = A;
    double x224 = x182*x223[0] + x186*x223[12] + x188*x223[9] + x190*x223[21] + x193*x223[3] + x196*x223[15] + x198*x223[6] + x200*x223[18];
    double x225 = x182*x223[1] + x186*x223[13] + x188*x223[10] + x190*x223[22] + x193*x223[4] + x196*x223[16] + x198*x223[7] + x200*x223[19];
    double x226 = x182*x223[2] + x186*x223[14] + x188*x223[11] + x190*x223[23] + x193*x223[5] + x196*x223[17] + x198*x223[8] + x200*x223[20];
    double x227 = x184*x46 + x187*x43 + x189*x48 + x191*x50 + x194*x40 + x197*x37 + x199*x52 + x201*x54;
    double x228 = x114*x98;
    double x229 = x184*x76 + x187*x73 + x189*x69 + x191*x66 + x194*x79 + x197*x81 + x199*x84 + x201*x86;
    double x230 = x124*x98;
    double x231 = x12*x194 + x16*x189 + x184*x19 + x187*x22 + x191*x24 + x197*x27 + x199*x8 + x201*x29;
    double x232 = x126*x98;
    double x233 = x203*x46 + x204*x43 + x205*x48 + x206*x50 + x207*x40 + x208*x37 + x209*x52 + x210*x54;
    double x234 = x203*x76 + x204*x73 + x205*x69 + x206*x66 + x207*x79 + x208*x81 + x209*x84 + x210*x86;
    double x235 = x12*x207 + x16*x205 + x19*x203 + x204*x22 + x206*x24 + x208*x27 + x209*x8 + x210*x29;
    double x236 = x212*x46 + x213*x43 + x214*x48 + x215*x50 + x216*x40 + x217*x37 + x218*x52 + x219*x54;
    double x237 = x212*x76 + x213*x73 + x214*x69 + x215*x66 + x216*x79 + x217*x81 + x218*x84 + x219*x86;
    double x238 = x12*x216 + x16*x214 + x19*x212 + x213*x22 + x215*x24 + x217*x27 + x218*x8 + x219*x29;
    double x239 = x128*x225 + x129*x226 + x145*x224 + x202*(x227*x228 + x229*x230 + x231*x232) + x211*(x228*x233 + x230*x234 + x232*x235) + x220*(x228*x236 + x230*x237 + x232*x238);
    double x240 = rho*x182;
    double x241 = -x113*x151 - x130*x150 - x144*x149;
    double x242 = -x113*x178 - x130*x177 - x144*x176;
    double x243 = -x113*x169 - x130*x168 - x144*x167 + x221;
    double x244 = x61*x98;
    double x245 = x109*x98;
    double x246 = x111*x98;
    double x247 = x113*x226 + x130*x225 + x144*x224 + x202*(x227*x245 + x229*x244 + x231*x246) + x211*(x233*x245 + x234*x244 + x235*x246) + x220*(x236*x245 + x237*x244 + x238*x246);
    double x248 = -x160*x168 - x161*x169 - x164*x167;
    double x249 = -x160*x177 - x161*x178 - x164*x176;
    double x250 = -x149*x164 - x150*x160 - x151*x161 + x221;
    double x251 = x153*x98;
    double x252 = x155*x98;
    double x253 = x157*x98;
    double x254 = x160*x225 + x161*x226 + x164*x224 + x202*(x227*x251 + x229*x252 + x231*x253) + x211*(x233*x251 + x234*x252 + x235*x253) + x220*(x236*x251 + x237*x252 + x238*x253);
    double x255 = x12*x253 + x251*x40 + x252*x79;
    double x256 = x12*x246 + x244*x79 + x245*x40;
    double x257 = x12*x232 + x228*x40 + x230*x79;
    double x258 = rho*x193;
    double x259 = x251*x52 + x252*x84 + x253*x8;
    double x260 = x244*x84 + x245*x52 + x246*x8;
    double x261 = x228*x52 + x230*x84 + x232*x8;
    double x262 = rho*x198;
    double x263 = x16*x253 + x251*x48 + x252*x69;
    double x264 = x16*x246 + x244*x69 + x245*x48;
    double x265 = x16*x232 + x228*x48 + x230*x69;
    double x266 = rho*x188;
    double x267 = x22*x253 + x251*x43 + x252*x73;
    double x268 = x22*x246 + x244*x73 + x245*x43;
    double x269 = x22*x232 + x228*x43 + x230*x73;
    double x270 = rho*x186;
    double x271 = x251*x37 + x252*x81 + x253*x27;
    double x272 = x244*x81 + x245*x37 + x246*x27;
    double x273 = x228*x37 + x230*x81 + x232*x27;
    double x274 = rho*x196;
    double x275 = x251*x54 + x252*x86 + x253*x29;
    double x276 = x244*x86 + x245*x54 + x246*x29;
    double x277 = x228*x54 + x230*x86 + x232*x29;
    double x278 = rho*x200;
    double x279 = x24*x253 + x251*x50 + x252*x66;
    double x280 = x24*x246 + x244*x66 + x245*x50;
    double x281 = x228*x50 + x230*x66 + x232*x24;
    double x282 = rho*x190;
    
    res_0[0] = x97*(x152*x159 + x170*x171 + x172*x222 - x239*x240);
    res_0[1] = x97*(x159*x241 + x171*x243 + x172*x242 - x240*x247);
    res_0[2] = x97*(x159*x250 + x171*x248 + x172*x249 - x240*x254);
    res_0[3] = x97*(x152*x255 + x170*x256 + x222*x257 - x239*x258);
    res_0[4] = x97*(x241*x255 + x242*x257 + x243*x256 - x247*x258);
    res_0[5] = x97*(x248*x256 + x249*x257 + x250*x255 - x254*x258);
    res_0[6] = x97*(x152*x259 + x170*x260 + x222*x261 - x239*x262);
    res_0[7] = x97*(x241*x259 + x242*x261 + x243*x260 - x247*x262);
    res_0[8] = x97*(x248*x260 + x249*x261 + x250*x259 - x254*x262);
    res_0[9] = x97*(x152*x263 + x170*x264 + x222*x265 - x239*x266);
    res_0[10] = x97*(x241*x263 + x242*x265 + x243*x264 - x247*x266);
    res_0[11] = x97*(x248*x264 + x249*x265 + x250*x263 - x254*x266);
    res_0[12] = x97*(x152*x267 + x170*x268 + x222*x269 - x239*x270);
    res_0[13] = x97*(x241*x267 + x242*x269 + x243*x268 - x247*x270);
    res_0[14] = x97*(x248*x268 + x249*x269 + x250*x267 - x254*x270);
    res_0[15] = x97*(x152*x271 + x170*x272 + x222*x273 - x239*x274);
    res_0[16] = x97*(x241*x271 + x242*x273 + x243*x272 - x247*x274);
    res_0[17] = x97*(x248*x272 + x249*x273 + x250*x271 - x254*x274);
    res_0[18] = x97*(x152*x275 + x170*x276 + x222*x277 - x239*x278);
    res_0[19] = x97*(x241*x275 + x242*x277 + x243*x276 - x247*x278);
    res_0[20] = x97*(x248*x276 + x249*x277 + x250*x275 - x254*x278);
    res_0[21] = x97*(x152*x279 + x170*x280 + x222*x281 - x239*x282);
    res_0[22] = x97*(x241*x279 + x242*x281 + x243*x280 - x247*x282);
    res_0[23] = x97*(x248*x280 + x249*x281 + x250*x279 - x254*x282);
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
    
    double x0 = (1.0/8.0)*r;
    double x1 = -x0;
    double x2 = x1 - 1.0/8.0;
    double x3 = (1.0/8.0)*s;
    double x4 = -x3;
    double x5 = s*x0;
    double x6 = -x5;
    double x7 = x4 + x6;
    double x8 = x2 + x7;
    double *x9 = coord;
    double x10 = x9[6];
    double x11 = x3 + x5;
    double x12 = x11 + x2;
    double x13 = x9[3];
    double x14 = x0 - 1.0/8.0;
    double x15 = x4 + x5;
    double x16 = x14 + x15;
    double x17 = x9[9];
    double x18 = x3 + x6;
    double x19 = x14 + x18;
    double x20 = x9[0];
    double x21 = x1 + 1.0/8.0;
    double x22 = x15 + x21;
    double x23 = x9[12];
    double x24 = x18 + x21;
    double x25 = x9[21];
    double x26 = x0 + 1.0/8.0;
    double x27 = x26 + x7;
    double x28 = x9[15];
    double x29 = x11 + x26;
    double x30 = x9[18];
    double x31 = x10*x8 + x12*x13 + x16*x17 + x19*x20 + x22*x23 + x24*x25 + x27*x28 + x29*x30;
    double x32 = (1.0/8.0)*t;
    double x33 = -x32;
    double x34 = t*x0;
    double x35 = -x34;
    double x36 = x33 + x35;
    double x37 = x2 + x36;
    double x38 = x9[17];
    double x39 = x32 + x34;
    double x40 = x2 + x39;
    double x41 = x9[5];
    double x42 = x33 + x34;
    double x43 = x14 + x42;
    double x44 = x9[14];
    double x45 = x32 + x35;
    double x46 = x14 + x45;
    double x47 = x9[2];
    double x48 = x21 + x42;
    double x49 = x9[11];
    double x50 = x21 + x45;
    double x51 = x9[23];
    double x52 = x26 + x36;
    double x53 = x9[8];
    double x54 = x26 + x39;
    double x55 = x9[20];
    double x56 = x37*x38 + x40*x41 + x43*x44 + x46*x47 + x48*x49 + x50*x51 + x52*x53 + x54*x55;
    double x57 = x31*x56;
    double x58 = x12*x41 + x16*x49 + x19*x47 + x22*x44 + x24*x51 + x27*x38 + x29*x55 + x53*x8;
    double x59 = x10*x52 + x13*x40 + x17*x48 + x20*x46 + x23*x43 + x25*x50 + x28*x37 + x30*x54;
    double x60 = x58*x59;
    double x61 = x57 - x60;
    double x62 = x4 - 1.0/8.0;
    double x63 = t*x3;
    double x64 = -x63;
    double x65 = x33 + x64;
    double x66 = x62 + x65;
    double x67 = x9[22];
    double x68 = x32 + x63;
    double x69 = x62 + x68;
    double x70 = x9[10];
    double x71 = x3 - 1.0/8.0;
    double x72 = x33 + x63;
    double x73 = x71 + x72;
    double x74 = x9[13];
    double x75 = x32 + x64;
    double x76 = x71 + x75;
    double x77 = x9[1];
    double x78 = x4 + 1.0/8.0;
    double x79 = x72 + x78;
    double x80 = x9[4];
    double x81 = x75 + x78;
    double x82 = x9[16];
    double x83 = x3 + 1.0/8.0;
    double x84 = x65 + x83;
    double x85 = x9[7];
    double x86 = x68 + x83;
    double x87 = x9[19];
    double x88 = x66*x67 + x69*x70 + x73*x74 + x76*x77 + x79*x80 + x81*x82 + x84*x85 + x86*x87;
    double x89 = x12*x80 + x16*x70 + x19*x77 + x22*x74 + x24*x67 + x27*x82 + x29*x87 + x8*x85;
    double x90 = x38*x81 + x41*x79 + x44*x73 + x47*x76 + x49*x69 + x51*x66 + x53*x84 + x55*x86;
    double x91 = x59*x90;
    double x92 = x37*x82 + x40*x80 + x43*x74 + x46*x77 + x48*x70 + x50*x67 + x52*x85 + x54*x87;
    double x93 = x10*x84 + x13*x79 + x17*x69 + x20*x76 + x23*x73 + x25*x66 + x28*x81 + x30*x86;
    double x94 = x58*x93;
    double x95 = x31*x90;
    double x96 = x56*x93;
    double x97 = x57*x88 - x60*x88 + x89*x91 - x89*x96 + x92*x94 - x92*x95;
    double x98 = 1.0/x97;
    double *x99 = U;
    double x100 = x99[23];
    double x101 = x99[11];
    double x102 = x99[14];
    double x103 = x99[2];
    double x104 = x99[5];
    double x105 = x99[17];
    double x106 = x99[8];
    double x107 = x99[20];
    double x108 = x98*(x100*x66 + x101*x69 + x102*x73 + x103*x76 + x104*x79 + x105*x81 + x106*x84 + x107*x86);
    double x109 = x94 - x95;
    double x110 = x98*(x100*x50 + x101*x48 + x102*x43 + x103*x46 + x104*x40 + x105*x37 + x106*x52 + x107*x54);
    double x111 = x91 - x96;
    double x112 = x98*(x100*x24 + x101*x16 + x102*x22 + x103*x19 + x104*x12 + x105*x27 + x106*x8 + x107*x29);
    double x113 = x108*x61 + x109*x110 + x111*x112;
    double x114 = -x58*x88 + x89*x90;
    double x115 = x99[16];
    double x116 = x99[4];
    double x117 = x99[13];
    double x118 = x99[1];
    double x119 = x99[10];
    double x120 = x99[22];
    double x121 = x99[7];
    double x122 = x99[19];
    double x123 = x98*(x115*x37 + x116*x40 + x117*x43 + x118*x46 + x119*x48 + x120*x50 + x121*x52 + x122*x54);
    double x124 = -x56*x89 + x58*x92;
    double x125 = x98*(x115*x81 + x116*x79 + x117*x73 + x118*x76 + x119*x69 + x120*x66 + x121*x84 + x122*x86);
    double x126 = x56*x88 - x90*x92;
    double x127 = x98*(x115*x27 + x116*x12 + x117*x22 + x118*x19 + x119*x16 + x120*x24 + x121*x8 + x122*x29);
    double x128 = x114*x123 + x124*x125 + x126*x127;
    double x129 = x108*x124 + x110*x114 + x112*x126;
    double x130 = x109*x123 + x111*x127 + x125*x61 + 1.0;
    double x131 = x113*x128 - x129*x130;
    double x132 = x99[21];
    double x133 = x99[9];
    double x134 = x99[12];
    double x135 = x99[0];
    double x136 = x99[3];
    double x137 = x99[15];
    double x138 = x99[6];
    double x139 = x99[18];
    double x140 = x98*(x132*x66 + x133*x69 + x134*x73 + x135*x76 + x136*x79 + x137*x81 + x138*x84 + x139*x86);
    double x141 = x98*(x132*x50 + x133*x48 + x134*x43 + x135*x46 + x136*x40 + x137*x37 + x138*x52 + x139*x54);
    double x142 = x98*(x12*x136 + x132*x24 + x133*x16 + x134*x22 + x135*x19 + x137*x27 + x138*x8 + x139*x29);
    double x143 = x109*x141 + x111*x142 + x140*x61;
    double x144 = x114*x141 + x124*x140 + x126*x142 + 1.0;
    double x145 = -x113*x144 + x129*x143;
    double x146 = -x128*x143 + x130*x144;
    double x147 = S11*x131 + S12*x145 + S13*x146;
    double x148 = S12*x131 + S22*x145 + S23*x146;
    double x149 = S13*x131 + S23*x145 + S33*x146;
    double x150 = -x128*x148 - x129*x149 - x144*x147;
    double x151 = x31*x88 - x89*x93;
    double x152 = x46*x98;
    double x153 = -x31*x92 + x59*x89;
    double x154 = x76*x98;
    double x155 = -x59*x88 + x92*x93;
    double x156 = x19*x98;
    double x157 = x151*x152 + x153*x154 + x155*x156;
    double x158 = x123*x151 + x125*x153 + x127*x155;
    double x159 = x108*x153 + x110*x151 + x112*x155 + 1.0;
    double x160 = -x128*x159 + x129*x158;
    double x161 = x140*x153 + x141*x151 + x142*x155;
    double x162 = -x129*x161 + x144*x159;
    double x163 = x128*x161 - x144*x158;
    double x164 = S11*x160 + S12*x162 + S13*x163;
    double x165 = S12*x160 + S22*x162 + S23*x163;
    double x166 = S13*x160 + S23*x162 + S33*x163;
    double x167 = -x128*x165 - x129*x166 - x144*x164;
    double x168 = x109*x152 + x111*x156 + x154*x61;
    double x169 = x114*x152 + x124*x154 + x126*x156;
    double x170 = -x113*x158 + x130*x159;
    double x171 = x113*x161 - x143*x159;
    double x172 = -x130*x161 + x143*x158;
    double x173 = S11*x170 + S12*x171 + S13*x172;
    double x174 = S12*x170 + S22*x171 + S23*x172;
    double x175 = S13*x170 + S23*x171 + S33*x172;
    double x176 = 1.0*PENER + 1.0*SENER;
    double x177 = -x128*x174 - x129*x175 - x144*x173 + x176;
    double x178 = -x113*x149 - x130*x148 - x143*x147;
    double x179 = -x113*x175 - x130*x174 - x143*x173;
    double x180 = -x113*x166 - x130*x165 - x143*x164 + x176;
    double x181 = -x158*x165 - x159*x166 - x161*x164;
    double x182 = -x158*x174 - x159*x175 - x161*x173;
    double x183 = -x147*x161 - x148*x158 - x149*x159 + x176;
    double x184 = x40*x98;
    double x185 = x79*x98;
    double x186 = x12*x98;
    double x187 = x151*x184 + x153*x185 + x155*x186;
    double x188 = x109*x184 + x111*x186 + x185*x61;
    double x189 = x114*x184 + x124*x185 + x126*x186;
    double x190 = x52*x98;
    double x191 = x84*x98;
    double x192 = x8*x98;
    double x193 = x151*x190 + x153*x191 + x155*x192;
    double x194 = x109*x190 + x111*x192 + x191*x61;
    double x195 = x114*x190 + x124*x191 + x126*x192;
    double x196 = x48*x98;
    double x197 = x69*x98;
    double x198 = x16*x98;
    double x199 = x151*x196 + x153*x197 + x155*x198;
    double x200 = x109*x196 + x111*x198 + x197*x61;
    double x201 = x114*x196 + x124*x197 + x126*x198;
    double x202 = x43*x98;
    double x203 = x73*x98;
    double x204 = x22*x98;
    double x205 = x151*x202 + x153*x203 + x155*x204;
    double x206 = x109*x202 + x111*x204 + x203*x61;
    double x207 = x114*x202 + x124*x203 + x126*x204;
    double x208 = x37*x98;
    double x209 = x81*x98;
    double x210 = x27*x98;
    double x211 = x151*x208 + x153*x209 + x155*x210;
    double x212 = x109*x208 + x111*x210 + x209*x61;
    double x213 = x114*x208 + x124*x209 + x126*x210;
    double x214 = x54*x98;
    double x215 = x86*x98;
    double x216 = x29*x98;
    double x217 = x151*x214 + x153*x215 + x155*x216;
    double x218 = x109*x214 + x111*x216 + x215*x61;
    double x219 = x114*x214 + x124*x215 + x126*x216;
    double x220 = x50*x98;
    double x221 = x66*x98;
    double x222 = x24*x98;
    double x223 = x151*x220 + x153*x221 + x155*x222;
    double x224 = x109*x220 + x111*x222 + x221*x61;
    double x225 = x114*x220 + x124*x221 + x126*x222;
    
    res_0[0] = x97*(x150*x157 + x167*x168 + x169*x177);
    res_0[1] = x97*(x157*x178 + x168*x180 + x169*x179);
    res_0[2] = x97*(x157*x183 + x168*x181 + x169*x182);
    res_0[3] = x97*(x150*x187 + x167*x188 + x177*x189);
    res_0[4] = x97*(x178*x187 + x179*x189 + x180*x188);
    res_0[5] = x97*(x181*x188 + x182*x189 + x183*x187);
    res_0[6] = x97*(x150*x193 + x167*x194 + x177*x195);
    res_0[7] = x97*(x178*x193 + x179*x195 + x180*x194);
    res_0[8] = x97*(x181*x194 + x182*x195 + x183*x193);
    res_0[9] = x97*(x150*x199 + x167*x200 + x177*x201);
    res_0[10] = x97*(x178*x199 + x179*x201 + x180*x200);
    res_0[11] = x97*(x181*x200 + x182*x201 + x183*x199);
    res_0[12] = x97*(x150*x205 + x167*x206 + x177*x207);
    res_0[13] = x97*(x178*x205 + x179*x207 + x180*x206);
    res_0[14] = x97*(x181*x206 + x182*x207 + x183*x205);
    res_0[15] = x97*(x150*x211 + x167*x212 + x177*x213);
    res_0[16] = x97*(x178*x211 + x179*x213 + x180*x212);
    res_0[17] = x97*(x181*x212 + x182*x213 + x183*x211);
    res_0[18] = x97*(x150*x217 + x167*x218 + x177*x219);
    res_0[19] = x97*(x178*x217 + x179*x219 + x180*x218);
    res_0[20] = x97*(x181*x218 + x182*x219 + x183*x217);
    res_0[21] = x97*(x150*x223 + x167*x224 + x177*x225);
    res_0[22] = x97*(x178*x223 + x179*x225 + x180*x224);
    res_0[23] = x97*(x181*x224 + x182*x225 + x183*x223);
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
    
    double x0 = (1.0/8.0)*r;
    double x1 = -x0;
    double x2 = x1 - 1.0/8.0;
    double x3 = (1.0/8.0)*s;
    double x4 = -x3;
    double x5 = s*x0;
    double x6 = -x5;
    double x7 = x4 + x6;
    double x8 = x2 + x7;
    double *x9 = coord;
    double x10 = x9[6];
    double x11 = x3 + x5;
    double x12 = x11 + x2;
    double x13 = x9[3];
    double x14 = x0 - 1.0/8.0;
    double x15 = x4 + x5;
    double x16 = x14 + x15;
    double x17 = x9[9];
    double x18 = x3 + x6;
    double x19 = x14 + x18;
    double x20 = x9[0];
    double x21 = x1 + 1.0/8.0;
    double x22 = x15 + x21;
    double x23 = x9[12];
    double x24 = x18 + x21;
    double x25 = x9[21];
    double x26 = x0 + 1.0/8.0;
    double x27 = x26 + x7;
    double x28 = x9[15];
    double x29 = x11 + x26;
    double x30 = x9[18];
    double x31 = x10*x8 + x12*x13 + x16*x17 + x19*x20 + x22*x23 + x24*x25 + x27*x28 + x29*x30;
    double x32 = (1.0/8.0)*t;
    double x33 = -x32;
    double x34 = t*x0;
    double x35 = -x34;
    double x36 = x33 + x35;
    double x37 = x2 + x36;
    double x38 = x9[17];
    double x39 = x32 + x34;
    double x40 = x2 + x39;
    double x41 = x9[5];
    double x42 = x33 + x34;
    double x43 = x14 + x42;
    double x44 = x9[14];
    double x45 = x32 + x35;
    double x46 = x14 + x45;
    double x47 = x9[2];
    double x48 = x21 + x42;
    double x49 = x9[11];
    double x50 = x21 + x45;
    double x51 = x9[23];
    double x52 = x26 + x36;
    double x53 = x9[8];
    double x54 = x26 + x39;
    double x55 = x9[20];
    double x56 = x37*x38 + x40*x41 + x43*x44 + x46*x47 + x48*x49 + x50*x51 + x52*x53 + x54*x55;
    double x57 = x31*x56;
    double x58 = x12*x41 + x16*x49 + x19*x47 + x22*x44 + x24*x51 + x27*x38 + x29*x55 + x53*x8;
    double x59 = x10*x52 + x13*x40 + x17*x48 + x20*x46 + x23*x43 + x25*x50 + x28*x37 + x30*x54;
    double x60 = x58*x59;
    double x61 = x57 - x60;
    double x62 = x4 - 1.0/8.0;
    double x63 = t*x3;
    double x64 = -x63;
    double x65 = x33 + x64;
    double x66 = x62 + x65;
    double x67 = x9[22];
    double x68 = x32 + x63;
    double x69 = x62 + x68;
    double x70 = x9[10];
    double x71 = x3 - 1.0/8.0;
    double x72 = x33 + x63;
    double x73 = x71 + x72;
    double x74 = x9[13];
    double x75 = x32 + x64;
    double x76 = x71 + x75;
    double x77 = x9[1];
    double x78 = x4 + 1.0/8.0;
    double x79 = x72 + x78;
    double x80 = x9[4];
    double x81 = x75 + x78;
    double x82 = x9[16];
    double x83 = x3 + 1.0/8.0;
    double x84 = x65 + x83;
    double x85 = x9[7];
    double x86 = x68 + x83;
    double x87 = x9[19];
    double x88 = x66*x67 + x69*x70 + x73*x74 + x76*x77 + x79*x80 + x81*x82 + x84*x85 + x86*x87;
    double x89 = x12*x80 + x16*x70 + x19*x77 + x22*x74 + x24*x67 + x27*x82 + x29*x87 + x8*x85;
    double x90 = x38*x81 + x41*x79 + x44*x73 + x47*x76 + x49*x69 + x51*x66 + x53*x84 + x55*x86;
    double x91 = x59*x90;
    double x92 = x37*x82 + x40*x80 + x43*x74 + x46*x77 + x48*x70 + x50*x67 + x52*x85 + x54*x87;
    double x93 = x10*x84 + x13*x79 + x17*x69 + x20*x76 + x23*x73 + x25*x66 + x28*x81 + x30*x86;
    double x94 = x58*x93;
    double x95 = x31*x90;
    double x96 = x56*x93;
    double x97 = x57*x88 - x60*x88 + x89*x91 - x89*x96 + x92*x94 - x92*x95;
    double x98 = 1.0/x97;
    double *x99 = U;
    double x100 = x99[23];
    double x101 = x99[11];
    double x102 = x99[14];
    double x103 = x99[2];
    double x104 = x99[5];
    double x105 = x99[17];
    double x106 = x99[8];
    double x107 = x99[20];
    double x108 = x98*(x100*x66 + x101*x69 + x102*x73 + x103*x76 + x104*x79 + x105*x81 + x106*x84 + x107*x86);
    double x109 = x94 - x95;
    double x110 = x98*(x100*x50 + x101*x48 + x102*x43 + x103*x46 + x104*x40 + x105*x37 + x106*x52 + x107*x54);
    double x111 = x91 - x96;
    double x112 = x98*(x100*x24 + x101*x16 + x102*x22 + x103*x19 + x104*x12 + x105*x27 + x106*x8 + x107*x29);
    double x113 = x108*x61 + x109*x110 + x111*x112;
    double x114 = -x58*x88 + x89*x90;
    double x115 = x99[16];
    double x116 = x99[4];
    double x117 = x99[13];
    double x118 = x99[1];
    double x119 = x99[10];
    double x120 = x99[22];
    double x121 = x99[7];
    double x122 = x99[19];
    double x123 = x98*(x115*x37 + x116*x40 + x117*x43 + x118*x46 + x119*x48 + x120*x50 + x121*x52 + x122*x54);
    double x124 = -x56*x89 + x58*x92;
    double x125 = x98*(x115*x81 + x116*x79 + x117*x73 + x118*x76 + x119*x69 + x120*x66 + x121*x84 + x122*x86);
    double x126 = x56*x88 - x90*x92;
    double x127 = x98*(x115*x27 + x116*x12 + x117*x22 + x118*x19 + x119*x16 + x120*x24 + x121*x8 + x122*x29);
    double x128 = x114*x123 + x124*x125 + x126*x127;
    double x129 = x108*x124 + x110*x114 + x112*x126;
    double x130 = x109*x123 + x111*x127 + x125*x61;
    double x131 = x130 + 1.0;
    double x132 = x113*x128 - x129*x131;
    double x133 = x99[21];
    double x134 = x99[9];
    double x135 = x99[12];
    double x136 = x99[0];
    double x137 = x99[3];
    double x138 = x99[15];
    double x139 = x99[6];
    double x140 = x99[18];
    double x141 = x98*(x133*x66 + x134*x69 + x135*x73 + x136*x76 + x137*x79 + x138*x81 + x139*x84 + x140*x86);
    double x142 = x98*(x133*x50 + x134*x48 + x135*x43 + x136*x46 + x137*x40 + x138*x37 + x139*x52 + x140*x54);
    double x143 = x98*(x12*x137 + x133*x24 + x134*x16 + x135*x22 + x136*x19 + x138*x27 + x139*x8 + x140*x29);
    double x144 = x109*x142 + x111*x143 + x141*x61;
    double x145 = x114*x142 + x124*x141 + x126*x143;
    double x146 = x145 + 1.0;
    double x147 = -x113*x146 + x129*x144;
    double x148 = -x128*x144 + x131*x146;
    double x149 = S11*x132 + S12*x147 + S13*x148;
    double x150 = S12*x132 + S22*x147 + S23*x148;
    double x151 = S13*x132 + S23*x147 + S33*x148;
    double x152 = -x128*x150 - x129*x151 - x145*x149;
    double x153 = x31*x88 - x89*x93;
    double x154 = x46*x98;
    double x155 = -x31*x92 + x59*x89;
    double x156 = x76*x98;
    double x157 = -x59*x88 + x92*x93;
    double x158 = x19*x98;
    double x159 = x153*x154 + x155*x156 + x157*x158;
    double x160 = x123*x153 + x125*x155 + x127*x157;
    double x161 = x108*x155 + x110*x153 + x112*x157;
    double x162 = x161 + 1.0;
    double x163 = -x128*x162 + x129*x160;
    double x164 = x141*x155 + x142*x153 + x143*x157;
    double x165 = -x129*x164 + x146*x162;
    double x166 = x128*x164 - x146*x160;
    double x167 = S11*x163 + S12*x165 + S13*x166;
    double x168 = S12*x163 + S22*x165 + S23*x166;
    double x169 = S13*x163 + S23*x165 + S33*x166;
    double x170 = -x128*x168 - x129*x169 - x145*x167;
    double x171 = x109*x154 + x111*x158 + x156*x61;
    double x172 = x114*x154 + x124*x156 + x126*x158;
    double x173 = -x113*x160 + x131*x162;
    double x174 = x113*x164 - x144*x162;
    double x175 = -x131*x164 + x144*x160;
    double x176 = S11*x173 + S12*x174 + S13*x175;
    double x177 = S12*x173 + S22*x174 + S23*x175;
    double x178 = S13*x173 + S23*x174 + S33*x175;
    double x179 = 1.0*PENER + 1.0*SENER;
    double x180 = -x128*x177 - x129*x178 - x145*x176 + x179;
    double x181 = -x113*x151 - x130*x150 - x144*x149;
    double x182 = -x113*x178 - x130*x177 - x144*x176;
    double x183 = -x113*x169 - x130*x168 - x144*x167 + x179;
    double x184 = -x160*x168 - x161*x169 - x164*x167;
    double x185 = -x160*x177 - x161*x178 - x164*x176;
    double x186 = -x149*x164 - x150*x160 - x151*x161 + x179;
    double x187 = x40*x98;
    double x188 = x79*x98;
    double x189 = x12*x98;
    double x190 = x153*x187 + x155*x188 + x157*x189;
    double x191 = x109*x187 + x111*x189 + x188*x61;
    double x192 = x114*x187 + x124*x188 + x126*x189;
    double x193 = x52*x98;
    double x194 = x84*x98;
    double x195 = x8*x98;
    double x196 = x153*x193 + x155*x194 + x157*x195;
    double x197 = x109*x193 + x111*x195 + x194*x61;
    double x198 = x114*x193 + x124*x194 + x126*x195;
    double x199 = x48*x98;
    double x200 = x69*x98;
    double x201 = x16*x98;
    double x202 = x153*x199 + x155*x200 + x157*x201;
    double x203 = x109*x199 + x111*x201 + x200*x61;
    double x204 = x114*x199 + x124*x200 + x126*x201;
    double x205 = x43*x98;
    double x206 = x73*x98;
    double x207 = x22*x98;
    double x208 = x153*x205 + x155*x206 + x157*x207;
    double x209 = x109*x205 + x111*x207 + x206*x61;
    double x210 = x114*x205 + x124*x206 + x126*x207;
    double x211 = x37*x98;
    double x212 = x81*x98;
    double x213 = x27*x98;
    double x214 = x153*x211 + x155*x212 + x157*x213;
    double x215 = x109*x211 + x111*x213 + x212*x61;
    double x216 = x114*x211 + x124*x212 + x126*x213;
    double x217 = x54*x98;
    double x218 = x86*x98;
    double x219 = x29*x98;
    double x220 = x153*x217 + x155*x218 + x157*x219;
    double x221 = x109*x217 + x111*x219 + x218*x61;
    double x222 = x114*x217 + x124*x218 + x126*x219;
    double x223 = x50*x98;
    double x224 = x66*x98;
    double x225 = x24*x98;
    double x226 = x153*x223 + x155*x224 + x157*x225;
    double x227 = x109*x223 + x111*x225 + x224*x61;
    double x228 = x114*x223 + x124*x224 + x126*x225;
    
    res_0[0] = x97*(x152*x159 + x170*x171 + x172*x180);
    res_0[1] = x97*(x159*x181 + x171*x183 + x172*x182);
    res_0[2] = x97*(x159*x186 + x171*x184 + x172*x185);
    res_0[3] = x97*(x152*x190 + x170*x191 + x180*x192);
    res_0[4] = x97*(x181*x190 + x182*x192 + x183*x191);
    res_0[5] = x97*(x184*x191 + x185*x192 + x186*x190);
    res_0[6] = x97*(x152*x196 + x170*x197 + x180*x198);
    res_0[7] = x97*(x181*x196 + x182*x198 + x183*x197);
    res_0[8] = x97*(x184*x197 + x185*x198 + x186*x196);
    res_0[9] = x97*(x152*x202 + x170*x203 + x180*x204);
    res_0[10] = x97*(x181*x202 + x182*x204 + x183*x203);
    res_0[11] = x97*(x184*x203 + x185*x204 + x186*x202);
    res_0[12] = x97*(x152*x208 + x170*x209 + x180*x210);
    res_0[13] = x97*(x181*x208 + x182*x210 + x183*x209);
    res_0[14] = x97*(x184*x209 + x185*x210 + x186*x208);
    res_0[15] = x97*(x152*x214 + x170*x215 + x180*x216);
    res_0[16] = x97*(x181*x214 + x182*x216 + x183*x215);
    res_0[17] = x97*(x184*x215 + x185*x216 + x186*x214);
    res_0[18] = x97*(x152*x220 + x170*x221 + x180*x222);
    res_0[19] = x97*(x181*x220 + x182*x222 + x183*x221);
    res_0[20] = x97*(x184*x221 + x185*x222 + x186*x220);
    res_0[21] = x97*(x152*x226 + x170*x227 + x180*x228);
    res_0[22] = x97*(x181*x226 + x182*x228 + x183*x227);
    res_0[23] = x97*(x184*x227 + x185*x228 + x186*x226);
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
    
    double x0 = s*s;
    double x1 = (1.0/4.0)*x0;
    double x2 = x1 - 1.0/4.0;
    double x3 = (1.0/4.0)*r;
    double x4 = -x3;
    double x5 = x0*x3;
    double x6 = x4 + x5;
    double x7 = x2 + x6;
    double *x8 = coord;
    double x9 = x8[27];
    double x10 = r*r;
    double x11 = (1.0/4.0)*x10;
    double x12 = x11 - 1.0/4.0;
    double x13 = (1.0/4.0)*s;
    double x14 = -x13;
    double x15 = x10*x13;
    double x16 = x14 + x15;
    double x17 = x12 + x16;
    double x18 = x8[30];
    double x19 = x3 - x5;
    double x20 = x19 + x2;
    double x21 = x8[33];
    double x22 = x13 - x15;
    double x23 = x12 + x22;
    double x24 = x8[24];
    double x25 = 1.0/4.0 - x1;
    double x26 = x25 + x6;
    double x27 = x8[45];
    double x28 = 1.0/4.0 - x11;
    double x29 = x16 + x28;
    double x30 = x8[36];
    double x31 = x22 + x28;
    double x32 = x8[42];
    double x33 = x19 + x25;
    double x34 = x8[39];
    double x35 = (1.0/2.0)*t;
    double x36 = -x35;
    double x37 = r*x35;
    double x38 = x36 + x37;
    double x39 = s*x35;
    double x40 = s*x37;
    double x41 = -x40;
    double x42 = x39 + x41;
    double x43 = x38 + x42;
    double x44 = x8[48];
    double x45 = -x39;
    double x46 = x40 + x45;
    double x47 = x38 + x46;
    double x48 = x8[57];
    double x49 = -x37;
    double x50 = x36 + x49;
    double x51 = x39 + x40;
    double x52 = x50 + x51;
    double x53 = x8[51];
    double x54 = x41 + x45;
    double x55 = x50 + x54;
    double x56 = x8[54];
    double x57 = (1.0/4.0)*t;
    double x58 = t*x13;
    double x59 = -x58;
    double x60 = t*x3;
    double x61 = s*x60;
    double x62 = x59 + x61;
    double x63 = x57 + x62;
    double x64 = (1.0/8.0)*x10;
    double x65 = s*x64;
    double x66 = -x65;
    double x67 = (1.0/8.0)*r;
    double x68 = s*x67;
    double x69 = x66 + x68;
    double x70 = (1.0/8.0)*x0;
    double x71 = r*x70;
    double x72 = -x71;
    double x73 = x64 - 1.0/8.0;
    double x74 = -x60;
    double x75 = x70 + x74;
    double x76 = x72 + x73 + x75;
    double x77 = x63 + x69 + x76;
    double x78 = x8[12];
    double x79 = -x68;
    double x80 = x65 + x79;
    double x81 = -x61;
    double x82 = x58 + x81;
    double x83 = x57 + x82;
    double x84 = x76 + x80 + x83;
    double x85 = x8[21];
    double x86 = x60 + x70;
    double x87 = x59 + x81;
    double x88 = x73 + x87;
    double x89 = x66 + x71;
    double x90 = x57 + x79;
    double x91 = x86 + x88 + x89 + x90;
    double x92 = x8[15];
    double x93 = x65 + x71;
    double x94 = x57 + x68;
    double x95 = x58 + x61;
    double x96 = x73 + x95;
    double x97 = x86 + x93 + x94 + x96;
    double x98 = x8[18];
    double x99 = -x70;
    double x100 = 1.0/8.0 - x64;
    double x101 = x100 + x99;
    double x102 = x101 + x71 + x74;
    double x103 = x102 + x63 + x80;
    double x104 = x8[0];
    double x105 = x102 + x69 + x83;
    double x106 = x8[9];
    double x107 = x65 + x72;
    double x108 = x101 + x60;
    double x109 = x107 + x108 + x87 + x94;
    double x110 = x8[3];
    double x111 = x66 + x72;
    double x112 = x108 + x111 + x90 + x95;
    double x113 = x8[6];
    double x114 = x103*x104 + x105*x106 + x109*x110 + x112*x113 + x17*x18 + x20*x21 + x23*x24 + x26*x27 + x29*x30 + x31*x32 + x33*x34 + x43*x44 + x47*x48 + x52*x53 + x55*x56 + x7*x9 + x77*x78 + x84*x85 + x91*x92 + x97*x98;
    double x115 = t*t;
    double x116 = (1.0/4.0)*x115;
    double x117 = x116 - 1.0/4.0;
    double x118 = x115*x3;
    double x119 = x118 + x4;
    double x120 = x117 + x119;
    double x121 = x8[53];
    double x122 = -x57;
    double x123 = t*x11;
    double x124 = x122 + x123;
    double x125 = x12 + x124;
    double x126 = x8[38];
    double x127 = -x118 + x3;
    double x128 = x117 + x127;
    double x129 = x8[50];
    double x130 = -x123 + x57;
    double x131 = x12 + x130;
    double x132 = x8[26];
    double x133 = 1.0/4.0 - x116;
    double x134 = x119 + x133;
    double x135 = x8[59];
    double x136 = x124 + x28;
    double x137 = x8[32];
    double x138 = x130 + x28;
    double x139 = x8[44];
    double x140 = x127 + x133;
    double x141 = x8[56];
    double x142 = (1.0/2.0)*s;
    double x143 = -x142;
    double x144 = r*x142;
    double x145 = x143 + x144;
    double x146 = x145 + x42;
    double x147 = x8[35];
    double x148 = x145 + x46;
    double x149 = x8[47];
    double x150 = -x144;
    double x151 = x143 + x150;
    double x152 = x151 + x51;
    double x153 = x8[29];
    double x154 = x151 + x54;
    double x155 = x8[41];
    double x156 = t*x64;
    double x157 = -x156;
    double x158 = x13 + x157;
    double x159 = (1.0/8.0)*x115;
    double x160 = r*x159;
    double x161 = -x160;
    double x162 = t*x67;
    double x163 = x161 + x162;
    double x164 = s*x3;
    double x165 = -x164;
    double x166 = x159 + x165;
    double x167 = x166 + x73;
    double x168 = x158 + x163 + x167 + x62;
    double x169 = x8[11];
    double x170 = x13 + x156;
    double x171 = -x162;
    double x172 = x161 + x171;
    double x173 = x167 + x170 + x172 + x82;
    double x174 = x8[23];
    double x175 = x160 + x171;
    double x176 = x158 + x164;
    double x177 = x159 + x175 + x176 + x88;
    double x178 = x8[8];
    double x179 = x160 + x162;
    double x180 = x164 + x170;
    double x181 = x159 + x179 + x180 + x96;
    double x182 = x8[20];
    double x183 = -x159;
    double x184 = x100 + x183;
    double x185 = x165 + x184;
    double x186 = x170 + x175 + x185 + x62;
    double x187 = x8[2];
    double x188 = x158 + x179 + x185 + x82;
    double x189 = x8[14];
    double x190 = x163 + x180 + x184 + x87;
    double x191 = x8[5];
    double x192 = x172 + x176 + x184 + x95;
    double x193 = x8[17];
    double x194 = x120*x121 + x125*x126 + x128*x129 + x131*x132 + x134*x135 + x136*x137 + x138*x139 + x140*x141 + x146*x147 + x148*x149 + x152*x153 + x154*x155 + x168*x169 + x173*x174 + x177*x178 + x181*x182 + x186*x187 + x188*x189 + x190*x191 + x192*x193;
    double x195 = x114*x194;
    double x196 = x103*x187 + x105*x169 + x109*x191 + x112*x178 + x121*x52 + x126*x29 + x129*x43 + x132*x23 + x135*x47 + x137*x17 + x139*x31 + x141*x55 + x147*x20 + x149*x26 + x153*x7 + x155*x33 + x174*x84 + x182*x97 + x189*x77 + x193*x91;
    double x197 = x104*x186 + x106*x168 + x110*x190 + x113*x177 + x120*x53 + x125*x30 + x128*x44 + x131*x24 + x134*x48 + x136*x18 + x138*x32 + x140*x56 + x146*x21 + x148*x27 + x152*x9 + x154*x34 + x173*x85 + x181*x98 + x188*x78 + x192*x92;
    double x198 = x196*x197;
    double x199 = x195 - x198;
    double x200 = x115*x13;
    double x201 = x14 + x200;
    double x202 = x117 + x201;
    double x203 = x8[58];
    double x204 = t*x1;
    double x205 = x122 + x204;
    double x206 = x2 + x205;
    double x207 = x8[46];
    double x208 = x13 - x200;
    double x209 = x117 + x208;
    double x210 = x8[49];
    double x211 = -x204 + x57;
    double x212 = x2 + x211;
    double x213 = x8[34];
    double x214 = x133 + x201;
    double x215 = x8[52];
    double x216 = x205 + x25;
    double x217 = x8[28];
    double x218 = x211 + x25;
    double x219 = x8[40];
    double x220 = x133 + x208;
    double x221 = x8[55];
    double x222 = -1.0/2.0*r;
    double x223 = x144 + x222;
    double x224 = x223 + x37 + x41;
    double x225 = x8[25];
    double x226 = x223 + x40 + x49;
    double x227 = x8[37];
    double x228 = x150 + x222;
    double x229 = x228 + x37 + x40;
    double x230 = x8[31];
    double x231 = x228 + x41 + x49;
    double x232 = x8[43];
    double x233 = s*x159;
    double x234 = -x233;
    double x235 = t*x70;
    double x236 = -x235;
    double x237 = x234 + x236;
    double x238 = x166 - 1.0/8.0;
    double x239 = (1.0/8.0)*s;
    double x240 = t*x239;
    double x241 = x3 + x61;
    double x242 = x240 + x241;
    double x243 = x237 + x238 + x242 + x75;
    double x244 = x8[4];
    double x245 = x235 + x86;
    double x246 = -x240;
    double x247 = x3 + x81;
    double x248 = x246 + x247;
    double x249 = x234 + x238 + x245 + x248;
    double x250 = x8[16];
    double x251 = x233 + x236;
    double x252 = x159 + x251;
    double x253 = x164 - 1.0/8.0;
    double x254 = x248 + x252 + x253 + x75;
    double x255 = x8[7];
    double x256 = x159 + x233;
    double x257 = x242 + x245 + x253 + x256;
    double x258 = x8[19];
    double x259 = x165 + x183 + x99 + 1.0/8.0;
    double x260 = x235 + x74;
    double x261 = x241 + x246;
    double x262 = x233 + x259 + x260 + x261;
    double x263 = x8[1];
    double x264 = x240 + x247;
    double x265 = x251 + x259 + x264 + x60;
    double x266 = x8[13];
    double x267 = x164 + x183 + x99 + 1.0/8.0;
    double x268 = x234 + x260 + x264 + x267;
    double x269 = x8[10];
    double x270 = x237 + x261 + x267 + x60;
    double x271 = x8[22];
    double x272 = x202*x203 + x206*x207 + x209*x210 + x212*x213 + x214*x215 + x216*x217 + x218*x219 + x220*x221 + x224*x225 + x226*x227 + x229*x230 + x231*x232 + x243*x244 + x249*x250 + x254*x255 + x257*x258 + x262*x263 + x265*x266 + x268*x269 + x270*x271;
    double x273 = x103*x263 + x105*x269 + x109*x244 + x112*x255 + x17*x230 + x20*x213 + x203*x47 + x207*x26 + x210*x43 + x215*x52 + x217*x7 + x219*x33 + x221*x55 + x225*x23 + x227*x29 + x232*x31 + x250*x91 + x258*x97 + x266*x77 + x271*x84;
    double x274 = x121*x214 + x126*x226 + x129*x209 + x132*x224 + x135*x202 + x137*x229 + x139*x231 + x141*x220 + x147*x212 + x149*x206 + x153*x216 + x155*x218 + x169*x268 + x174*x270 + x178*x254 + x182*x257 + x187*x262 + x189*x265 + x191*x243 + x193*x249;
    double x275 = x197*x274;
    double x276 = x120*x215 + x125*x227 + x128*x210 + x131*x225 + x134*x203 + x136*x230 + x138*x232 + x140*x221 + x146*x213 + x148*x207 + x152*x217 + x154*x219 + x168*x269 + x173*x271 + x177*x255 + x181*x258 + x186*x263 + x188*x266 + x190*x244 + x192*x250;
    double x277 = x104*x262 + x106*x268 + x110*x243 + x113*x254 + x18*x229 + x202*x48 + x206*x27 + x209*x44 + x21*x212 + x214*x53 + x216*x9 + x218*x34 + x220*x56 + x224*x24 + x226*x30 + x231*x32 + x249*x92 + x257*x98 + x265*x78 + x270*x85;
    double x278 = x196*x277;
    double x279 = x114*x274;
    double x280 = x194*x277;
    double x281 = x195*x272 - x198*x272 + x273*x275 - x273*x280 + x276*x278 - x276*x279;
    double x282 = 1.0/x281;
    double *x283 = U;
    double x284 = x283[59];
    double x285 = x283[47];
    double x286 = x283[50];
    double x287 = x283[35];
    double x288 = x283[53];
    double x289 = x283[29];
    double x290 = x283[41];
    double x291 = x283[56];
    double x292 = x283[26];
    double x293 = x283[38];
    double x294 = x283[32];
    double x295 = x283[44];
    double x296 = x283[5];
    double x297 = x283[17];
    double x298 = x283[8];
    double x299 = x283[20];
    double x300 = x283[2];
    double x301 = x283[14];
    double x302 = x283[11];
    double x303 = x283[23];
    double x304 = x282*(x202*x284 + x206*x285 + x209*x286 + x212*x287 + x214*x288 + x216*x289 + x218*x290 + x220*x291 + x224*x292 + x226*x293 + x229*x294 + x231*x295 + x243*x296 + x249*x297 + x254*x298 + x257*x299 + x262*x300 + x265*x301 + x268*x302 + x270*x303);
    double x305 = x278 - x279;
    double x306 = x282*(x120*x288 + x125*x293 + x128*x286 + x131*x292 + x134*x284 + x136*x294 + x138*x295 + x140*x291 + x146*x287 + x148*x285 + x152*x289 + x154*x290 + x168*x302 + x173*x303 + x177*x298 + x181*x299 + x186*x300 + x188*x301 + x190*x296 + x192*x297);
    double x307 = x275 - x280;
    double x308 = x282*(x103*x300 + x105*x302 + x109*x296 + x112*x298 + x17*x294 + x20*x287 + x23*x292 + x26*x285 + x284*x47 + x286*x43 + x288*x52 + x289*x7 + x29*x293 + x290*x33 + x291*x55 + x295*x31 + x297*x91 + x299*x97 + x301*x77 + x303*x84);
    double x309 = x199*x304 + x305*x306 + x307*x308;
    double x310 = -x196*x272 + x273*x274;
    double x311 = x283[52];
    double x312 = x283[37];
    double x313 = x283[49];
    double x314 = x283[25];
    double x315 = x283[58];
    double x316 = x283[31];
    double x317 = x283[43];
    double x318 = x283[55];
    double x319 = x283[34];
    double x320 = x283[46];
    double x321 = x283[28];
    double x322 = x283[40];
    double x323 = x283[10];
    double x324 = x283[22];
    double x325 = x283[7];
    double x326 = x283[19];
    double x327 = x283[1];
    double x328 = x283[13];
    double x329 = x283[4];
    double x330 = x283[16];
    double x331 = x282*(x120*x311 + x125*x312 + x128*x313 + x131*x314 + x134*x315 + x136*x316 + x138*x317 + x140*x318 + x146*x319 + x148*x320 + x152*x321 + x154*x322 + x168*x323 + x173*x324 + x177*x325 + x181*x326 + x186*x327 + x188*x328 + x190*x329 + x192*x330);
    double x332 = -x194*x273 + x196*x276;
    double x333 = x282*(x202*x315 + x206*x320 + x209*x313 + x212*x319 + x214*x311 + x216*x321 + x218*x322 + x220*x318 + x224*x314 + x226*x312 + x229*x316 + x231*x317 + x243*x329 + x249*x330 + x254*x325 + x257*x326 + x262*x327 + x265*x328 + x268*x323 + x270*x324);
    double x334 = x194*x272 - x274*x276;
    double x335 = x282*(x103*x327 + x105*x323 + x109*x329 + x112*x325 + x17*x316 + x20*x319 + x23*x314 + x26*x320 + x29*x312 + x31*x317 + x311*x52 + x313*x43 + x315*x47 + x318*x55 + x321*x7 + x322*x33 + x324*x84 + x326*x97 + x328*x77 + x330*x91);
    double x336 = x310*x331 + x332*x333 + x334*x335;
    double x337 = x304*x332 + x306*x310 + x308*x334;
    double x338 = x199*x333 + x305*x331 + x307*x335;
    double x339 = x338 + 1.0;
    double x340 = x309*x336 - x337*x339;
    double x341 = x283[57];
    double x342 = x283[45];
    double x343 = x283[48];
    double x344 = x283[33];
    double x345 = x283[51];
    double x346 = x283[27];
    double x347 = x283[39];
    double x348 = x283[54];
    double x349 = x283[24];
    double x350 = x283[36];
    double x351 = x283[30];
    double x352 = x283[42];
    double x353 = x283[3];
    double x354 = x283[15];
    double x355 = x283[6];
    double x356 = x283[18];
    double x357 = x283[0];
    double x358 = x283[12];
    double x359 = x283[9];
    double x360 = x283[21];
    double x361 = x282*(x202*x341 + x206*x342 + x209*x343 + x212*x344 + x214*x345 + x216*x346 + x218*x347 + x220*x348 + x224*x349 + x226*x350 + x229*x351 + x231*x352 + x243*x353 + x249*x354 + x254*x355 + x257*x356 + x262*x357 + x265*x358 + x268*x359 + x270*x360);
    double x362 = x282*(x120*x345 + x125*x350 + x128*x343 + x131*x349 + x134*x341 + x136*x351 + x138*x352 + x140*x348 + x146*x344 + x148*x342 + x152*x346 + x154*x347 + x168*x359 + x173*x360 + x177*x355 + x181*x356 + x186*x357 + x188*x358 + x190*x353 + x192*x354);
    double x363 = x282*(x103*x357 + x105*x359 + x109*x353 + x112*x355 + x17*x351 + x20*x344 + x23*x349 + x26*x342 + x29*x350 + x31*x352 + x33*x347 + x341*x47 + x343*x43 + x345*x52 + x346*x7 + x348*x55 + x354*x91 + x356*x97 + x358*x77 + x360*x84);
    double x364 = x199*x361 + x305*x362 + x307*x363;
    double x365 = x310*x362 + x332*x361 + x334*x363;
    double x366 = x365 + 1.0;
    double x367 = -x309*x366 + x337*x364;
    double x368 = -x336*x364 + x339*x366;
    double x369 = S11*x340 + S12*x367 + S13*x368;
    double x370 = S12*x340 + S22*x367 + S23*x368;
    double x371 = S13*x340 + S23*x367 + S33*x368;
    double x372 = -x336*x370 - x337*x371 - x365*x369;
    double x373 = x114*x272 - x273*x277;
    double x374 = x186*x282;
    double x375 = -x114*x276 + x197*x273;
    double x376 = x262*x282;
    double x377 = -x197*x272 + x276*x277;
    double x378 = x103*x282;
    double x379 = x373*x374 + x375*x376 + x377*x378;
    double x380 = x331*x373 + x333*x375 + x335*x377;
    double x381 = x304*x375 + x306*x373 + x308*x377;
    double x382 = x381 + 1.0;
    double x383 = -x336*x382 + x337*x380;
    double x384 = x361*x375 + x362*x373 + x363*x377;
    double x385 = -x337*x384 + x366*x382;
    double x386 = x336*x384 - x366*x380;
    double x387 = S11*x383 + S12*x385 + S13*x386;
    double x388 = S12*x383 + S22*x385 + S23*x386;
    double x389 = S13*x383 + S23*x385 + S33*x386;
    double x390 = -x336*x388 - x337*x389 - x365*x387;
    double x391 = x199*x376 + x305*x374 + x307*x378;
    double x392 = x310*x374 + x332*x376 + x334*x378;
    double x393 = -x309*x380 + x339*x382;
    double x394 = x309*x384 - x364*x382;
    double x395 = -x339*x384 + x364*x380;
    double x396 = S11*x393 + S12*x394 + S13*x395;
    double x397 = S12*x393 + S22*x394 + S23*x395;
    double x398 = S13*x393 + S23*x394 + S33*x395;
    double x399 = s*x118;
    double x400 = x164 - x399;
    double x401 = x134 + x201 + x400;
    double *x402 = V;
    double x403 = x402[48];
    double x404 = t*x5;
    double x405 = -x404 + x60;
    double x406 = x205 + x26 + x405;
    double x407 = x402[33];
    double x408 = x404 + x74;
    double x409 = x211 + x26 + x408;
    double x410 = x402[45];
    double x411 = x165 + x399;
    double x412 = x134 + x208 + x411;
    double x413 = x402[57];
    double x414 = t*x15;
    double x415 = -x414 + x58;
    double x416 = x124 + x29 + x415;
    double x417 = x402[24];
    double x418 = x414 + x59;
    double x419 = x130 + x29 + x418;
    double x420 = x402[36];
    double x421 = x140 + x201 + x411;
    double x422 = x402[51];
    double x423 = x124 + x31 + x418;
    double x424 = x402[30];
    double x425 = x205 + x33 + x408;
    double x426 = x402[27];
    double x427 = x130 + x31 + x415;
    double x428 = x402[42];
    double x429 = x211 + x33 + x405;
    double x430 = x402[39];
    double x431 = x140 + x208 + x400;
    double x432 = x402[54];
    double x433 = -x239;
    double x434 = t*x68;
    double x435 = t*x71;
    double x436 = x433 + x434 + x435 + x64 + x70 - 1.0/4.0;
    double x437 = s*x160;
    double x438 = t*x65;
    double x439 = x437 + x438;
    double x440 = x160 - x67;
    double x441 = x440 + x93;
    double x442 = (1.0/8.0)*t;
    double x443 = x156 + x235 - x442;
    double x444 = x256 + x443;
    double x445 = x436 + x439 + x441 + x444;
    double x446 = x402[18];
    double x447 = -x438;
    double x448 = x157 + x252 + x442 + x447;
    double x449 = -x435;
    double x450 = -x434 + x64 + x70 - 1.0/4.0;
    double x451 = x433 + x449 + x450;
    double x452 = x437 + x441 + x448 + x451;
    double x453 = x402[6];
    double x454 = x159 + x239 + x435 + x450;
    double x455 = -x437;
    double x456 = x440 + x455 + x89;
    double x457 = x234 + x443 + x447;
    double x458 = x454 + x456 + x457;
    double x459 = x402[15];
    double x460 = x157 + x237 + x442;
    double x461 = x159 + x239 + x434 + x449 + x64 + x70 - 1.0/4.0;
    double x462 = x438 + x456 + x460 + x461;
    double x463 = x402[3];
    double x464 = x161 + x67;
    double x465 = x107 + x455 + x464;
    double x466 = x438 + x444 + x451 + x465;
    double x467 = x402[21];
    double x468 = x436 + x448 + x465;
    double x469 = x402[9];
    double x470 = x111 + x464;
    double x471 = x437 + x457 + x461 + x470;
    double x472 = x402[12];
    double x473 = x439 + x454 + x460 + x470;
    double x474 = x402[0];
    double x475 = x401*x403 + x406*x407 + x409*x410 + x412*x413 + x416*x417 + x419*x420 + x421*x422 + x423*x424 + x425*x426 + x427*x428 + x429*x430 + x431*x432 + x445*x446 + x452*x453 + x458*x459 + x462*x463 + x466*x467 + x468*x469 + x471*x472 + x473*x474;
    double x476 = x402[49];
    double x477 = x402[34];
    double x478 = x402[46];
    double x479 = x402[58];
    double x480 = x402[25];
    double x481 = x402[37];
    double x482 = x402[52];
    double x483 = x402[31];
    double x484 = x402[28];
    double x485 = x402[43];
    double x486 = x402[40];
    double x487 = x402[55];
    double x488 = x402[19];
    double x489 = x402[7];
    double x490 = x402[16];
    double x491 = x402[4];
    double x492 = x402[22];
    double x493 = x402[10];
    double x494 = x402[13];
    double x495 = x402[1];
    double x496 = x401*x476 + x406*x477 + x409*x478 + x412*x479 + x416*x480 + x419*x481 + x421*x482 + x423*x483 + x425*x484 + x427*x485 + x429*x486 + x431*x487 + x445*x488 + x452*x489 + x458*x490 + x462*x491 + x466*x492 + x468*x493 + x471*x494 + x473*x495;
    double x497 = x402[50];
    double x498 = x402[35];
    double x499 = x402[47];
    double x500 = x402[59];
    double x501 = x402[26];
    double x502 = x402[38];
    double x503 = x402[53];
    double x504 = x402[32];
    double x505 = x402[29];
    double x506 = x402[44];
    double x507 = x402[41];
    double x508 = x402[56];
    double x509 = x402[20];
    double x510 = x402[8];
    double x511 = x402[17];
    double x512 = x402[5];
    double x513 = x402[23];
    double x514 = x402[11];
    double x515 = x402[14];
    double x516 = x402[2];
    double x517 = x401*x497 + x406*x498 + x409*x499 + x412*x500 + x416*x501 + x419*x502 + x421*x503 + x423*x504 + x425*x505 + x427*x506 + x429*x507 + x431*x508 + x445*x509 + x452*x510 + x458*x511 + x462*x512 + x466*x513 + x468*x514 + x471*x515 + x473*x516;
    double x518 = 1.0*PENER + 1.0*SENER - 0.5*rho*(x475*x475 + x496*x496 + x517*x517);
    double x519 = -x336*x397 - x337*x398 - x365*x396 + x518;
    double *x520 = A;
    double x521 = x401*x520[48] + x406*x520[33] + x409*x520[45] + x412*x520[57] + x416*x520[24] + x419*x520[36] + x421*x520[51] + x423*x520[30] + x425*x520[27] + x427*x520[42] + x429*x520[39] + x431*x520[54] + x445*x520[18] + x452*x520[6] + x458*x520[15] + x462*x520[3] + x466*x520[21] + x468*x520[9] + x471*x520[12] + x473*x520[0];
    double x522 = x401*x520[49] + x406*x520[34] + x409*x520[46] + x412*x520[58] + x416*x520[25] + x419*x520[37] + x421*x520[52] + x423*x520[31] + x425*x520[28] + x427*x520[43] + x429*x520[40] + x431*x520[55] + x445*x520[19] + x452*x520[7] + x458*x520[16] + x462*x520[4] + x466*x520[22] + x468*x520[10] + x471*x520[13] + x473*x520[1];
    double x523 = x401*x520[50] + x406*x520[35] + x409*x520[47] + x412*x520[59] + x416*x520[26] + x419*x520[38] + x421*x520[53] + x423*x520[32] + x425*x520[29] + x427*x520[44] + x429*x520[41] + x431*x520[56] + x445*x520[20] + x452*x520[8] + x458*x520[17] + x462*x520[5] + x466*x520[23] + x468*x520[11] + x471*x520[14] + x473*x520[2];
    double x524 = x120*x422 + x125*x420 + x128*x403 + x131*x417 + x134*x413 + x136*x424 + x138*x428 + x140*x432 + x146*x407 + x148*x410 + x152*x426 + x154*x430 + x168*x469 + x173*x467 + x177*x453 + x181*x446 + x186*x474 + x188*x472 + x190*x463 + x192*x459;
    double x525 = x282*x310;
    double x526 = x202*x413 + x206*x410 + x209*x403 + x212*x407 + x214*x422 + x216*x426 + x218*x430 + x220*x432 + x224*x417 + x226*x420 + x229*x424 + x231*x428 + x243*x463 + x249*x459 + x254*x453 + x257*x446 + x262*x474 + x265*x472 + x268*x469 + x270*x467;
    double x527 = x282*x332;
    double x528 = x103*x474 + x105*x469 + x109*x463 + x112*x453 + x17*x424 + x20*x407 + x23*x417 + x26*x410 + x29*x420 + x31*x428 + x33*x430 + x403*x43 + x413*x47 + x422*x52 + x426*x7 + x432*x55 + x446*x97 + x459*x91 + x467*x84 + x472*x77;
    double x529 = x282*x334;
    double x530 = x120*x482 + x125*x481 + x128*x476 + x131*x480 + x134*x479 + x136*x483 + x138*x485 + x140*x487 + x146*x477 + x148*x478 + x152*x484 + x154*x486 + x168*x493 + x173*x492 + x177*x489 + x181*x488 + x186*x495 + x188*x494 + x190*x491 + x192*x490;
    double x531 = x202*x479 + x206*x478 + x209*x476 + x212*x477 + x214*x482 + x216*x484 + x218*x486 + x220*x487 + x224*x480 + x226*x481 + x229*x483 + x231*x485 + x243*x491 + x249*x490 + x254*x489 + x257*x488 + x262*x495 + x265*x494 + x268*x493 + x270*x492;
    double x532 = x103*x495 + x105*x493 + x109*x491 + x112*x489 + x17*x483 + x20*x477 + x23*x480 + x26*x478 + x29*x481 + x31*x485 + x33*x486 + x43*x476 + x47*x479 + x482*x52 + x484*x7 + x487*x55 + x488*x97 + x490*x91 + x492*x84 + x494*x77;
    double x533 = x120*x503 + x125*x502 + x128*x497 + x131*x501 + x134*x500 + x136*x504 + x138*x506 + x140*x508 + x146*x498 + x148*x499 + x152*x505 + x154*x507 + x168*x514 + x173*x513 + x177*x510 + x181*x509 + x186*x516 + x188*x515 + x190*x512 + x192*x511;
    double x534 = x202*x500 + x206*x499 + x209*x497 + x212*x498 + x214*x503 + x216*x505 + x218*x507 + x220*x508 + x224*x501 + x226*x502 + x229*x504 + x231*x506 + x243*x512 + x249*x511 + x254*x510 + x257*x509 + x262*x516 + x265*x515 + x268*x514 + x270*x513;
    double x535 = x103*x516 + x105*x514 + x109*x512 + x112*x510 + x17*x504 + x20*x498 + x23*x501 + x26*x499 + x29*x502 + x31*x506 + x33*x507 + x43*x497 + x47*x500 + x503*x52 + x505*x7 + x508*x55 + x509*x97 + x511*x91 + x513*x84 + x515*x77;
    double x536 = x336*x522 + x337*x523 + x365*x521 + x475*(x524*x525 + x526*x527 + x528*x529) + x496*(x525*x530 + x527*x531 + x529*x532) + x517*(x525*x533 + x527*x534 + x529*x535);
    double x537 = rho*x473;
    double x538 = -x309*x371 - x338*x370 - x364*x369;
    double x539 = -x309*x398 - x338*x397 - x364*x396;
    double x540 = -x309*x389 - x338*x388 - x364*x387 + x518;
    double x541 = x199*x282;
    double x542 = x282*x305;
    double x543 = x282*x307;
    double x544 = x309*x523 + x338*x522 + x364*x521 + x475*(x524*x542 + x526*x541 + x528*x543) + x496*(x530*x542 + x531*x541 + x532*x543) + x517*(x533*x542 + x534*x541 + x535*x543);
    double x545 = -x380*x388 - x381*x389 - x384*x387;
    double x546 = -x380*x397 - x381*x398 - x384*x396;
    double x547 = -x369*x384 - x370*x380 - x371*x381 + x518;
    double x548 = x282*x373;
    double x549 = x282*x375;
    double x550 = x282*x377;
    double x551 = x380*x522 + x381*x523 + x384*x521 + x475*(x524*x548 + x526*x549 + x528*x550) + x496*(x530*x548 + x531*x549 + x532*x550) + x517*(x533*x548 + x534*x549 + x535*x550);
    double x552 = x109*x550 + x190*x548 + x243*x549;
    double x553 = x109*x543 + x190*x542 + x243*x541;
    double x554 = x109*x529 + x190*x525 + x243*x527;
    double x555 = rho*x462;
    double x556 = x112*x550 + x177*x548 + x254*x549;
    double x557 = x112*x543 + x177*x542 + x254*x541;
    double x558 = x112*x529 + x177*x525 + x254*x527;
    double x559 = rho*x452;
    double x560 = x105*x550 + x168*x548 + x268*x549;
    double x561 = x105*x543 + x168*x542 + x268*x541;
    double x562 = x105*x529 + x168*x525 + x268*x527;
    double x563 = rho*x468;
    double x564 = x188*x548 + x265*x549 + x550*x77;
    double x565 = x188*x542 + x265*x541 + x543*x77;
    double x566 = x188*x525 + x265*x527 + x529*x77;
    double x567 = rho*x471;
    double x568 = x192*x548 + x249*x549 + x550*x91;
    double x569 = x192*x542 + x249*x541 + x543*x91;
    double x570 = x192*x525 + x249*x527 + x529*x91;
    double x571 = rho*x458;
    double x572 = x181*x548 + x257*x549 + x550*x97;
    double x573 = x181*x542 + x257*x541 + x543*x97;
    double x574 = x181*x525 + x257*x527 + x529*x97;
    double x575 = rho*x445;
    double x576 = x173*x548 + x270*x549 + x550*x84;
    double x577 = x173*x542 + x270*x541 + x543*x84;
    double x578 = x173*x525 + x270*x527 + x529*x84;
    double x579 = rho*x466;
    double x580 = x131*x548 + x224*x549 + x23*x550;
    double x581 = x131*x542 + x224*x541 + x23*x543;
    double x582 = x131*x525 + x224*x527 + x23*x529;
    double x583 = rho*x416;
    double x584 = x152*x548 + x216*x549 + x550*x7;
    double x585 = x152*x542 + x216*x541 + x543*x7;
    double x586 = x152*x525 + x216*x527 + x529*x7;
    double x587 = rho*x425;
    double x588 = x136*x548 + x17*x550 + x229*x549;
    double x589 = x136*x542 + x17*x543 + x229*x541;
    double x590 = x136*x525 + x17*x529 + x229*x527;
    double x591 = rho*x423;
    double x592 = x146*x548 + x20*x550 + x212*x549;
    double x593 = x146*x542 + x20*x543 + x212*x541;
    double x594 = x146*x525 + x20*x529 + x212*x527;
    double x595 = rho*x406;
    double x596 = x125*x548 + x226*x549 + x29*x550;
    double x597 = x125*x542 + x226*x541 + x29*x543;
    double x598 = x125*x525 + x226*x527 + x29*x529;
    double x599 = rho*x419;
    double x600 = x154*x548 + x218*x549 + x33*x550;
    double x601 = x154*x542 + x218*x541 + x33*x543;
    double x602 = x154*x525 + x218*x527 + x33*x529;
    double x603 = rho*x429;
    double x604 = x138*x548 + x231*x549 + x31*x550;
    double x605 = x138*x542 + x231*x541 + x31*x543;
    double x606 = x138*x525 + x231*x527 + x31*x529;
    double x607 = rho*x427;
    double x608 = x148*x548 + x206*x549 + x26*x550;
    double x609 = x148*x542 + x206*x541 + x26*x543;
    double x610 = x148*x525 + x206*x527 + x26*x529;
    double x611 = rho*x409;
    double x612 = x128*x548 + x209*x549 + x43*x550;
    double x613 = x128*x542 + x209*x541 + x43*x543;
    double x614 = x128*x525 + x209*x527 + x43*x529;
    double x615 = rho*x401;
    double x616 = x120*x548 + x214*x549 + x52*x550;
    double x617 = x120*x542 + x214*x541 + x52*x543;
    double x618 = x120*x525 + x214*x527 + x52*x529;
    double x619 = rho*x421;
    double x620 = x140*x548 + x220*x549 + x55*x550;
    double x621 = x140*x542 + x220*x541 + x543*x55;
    double x622 = x140*x525 + x220*x527 + x529*x55;
    double x623 = rho*x431;
    double x624 = x134*x548 + x202*x549 + x47*x550;
    double x625 = x134*x542 + x202*x541 + x47*x543;
    double x626 = x134*x525 + x202*x527 + x47*x529;
    double x627 = rho*x412;
    
    res_0[0] = x281*(x372*x379 + x390*x391 + x392*x519 - x536*x537);
    res_0[1] = x281*(x379*x538 + x391*x540 + x392*x539 - x537*x544);
    res_0[2] = x281*(x379*x547 + x391*x545 + x392*x546 - x537*x551);
    res_0[3] = x281*(x372*x552 + x390*x553 + x519*x554 - x536*x555);
    res_0[4] = x281*(x538*x552 + x539*x554 + x540*x553 - x544*x555);
    res_0[5] = x281*(x545*x553 + x546*x554 + x547*x552 - x551*x555);
    res_0[6] = x281*(x372*x556 + x390*x557 + x519*x558 - x536*x559);
    res_0[7] = x281*(x538*x556 + x539*x558 + x540*x557 - x544*x559);
    res_0[8] = x281*(x545*x557 + x546*x558 + x547*x556 - x551*x559);
    res_0[9] = x281*(x372*x560 + x390*x561 + x519*x562 - x536*x563);
    res_0[10] = x281*(x538*x560 + x539*x562 + x540*x561 - x544*x563);
    res_0[11] = x281*(x545*x561 + x546*x562 + x547*x560 - x551*x563);
    res_0[12] = x281*(x372*x564 + x390*x565 + x519*x566 - x536*x567);
    res_0[13] = x281*(x538*x564 + x539*x566 + x540*x565 - x544*x567);
    res_0[14] = x281*(x545*x565 + x546*x566 + x547*x564 - x551*x567);
    res_0[15] = x281*(x372*x568 + x390*x569 + x519*x570 - x536*x571);
    res_0[16] = x281*(x538*x568 + x539*x570 + x540*x569 - x544*x571);
    res_0[17] = x281*(x545*x569 + x546*x570 + x547*x568 - x551*x571);
    res_0[18] = x281*(x372*x572 + x390*x573 + x519*x574 - x536*x575);
    res_0[19] = x281*(x538*x572 + x539*x574 + x540*x573 - x544*x575);
    res_0[20] = x281*(x545*x573 + x546*x574 + x547*x572 - x551*x575);
    res_0[21] = x281*(x372*x576 + x390*x577 + x519*x578 - x536*x579);
    res_0[22] = x281*(x538*x576 + x539*x578 + x540*x577 - x544*x579);
    res_0[23] = x281*(x545*x577 + x546*x578 + x547*x576 - x551*x579);
    res_0[24] = x281*(x372*x580 + x390*x581 + x519*x582 - x536*x583);
    res_0[25] = x281*(x538*x580 + x539*x582 + x540*x581 - x544*x583);
    res_0[26] = x281*(x545*x581 + x546*x582 + x547*x580 - x551*x583);
    res_0[27] = x281*(x372*x584 + x390*x585 + x519*x586 - x536*x587);
    res_0[28] = x281*(x538*x584 + x539*x586 + x540*x585 - x544*x587);
    res_0[29] = x281*(x545*x585 + x546*x586 + x547*x584 - x551*x587);
    res_0[30] = x281*(x372*x588 + x390*x589 + x519*x590 - x536*x591);
    res_0[31] = x281*(x538*x588 + x539*x590 + x540*x589 - x544*x591);
    res_0[32] = x281*(x545*x589 + x546*x590 + x547*x588 - x551*x591);
    res_0[33] = x281*(x372*x592 + x390*x593 + x519*x594 - x536*x595);
    res_0[34] = x281*(x538*x592 + x539*x594 + x540*x593 - x544*x595);
    res_0[35] = x281*(x545*x593 + x546*x594 + x547*x592 - x551*x595);
    res_0[36] = x281*(x372*x596 + x390*x597 + x519*x598 - x536*x599);
    res_0[37] = x281*(x538*x596 + x539*x598 + x540*x597 - x544*x599);
    res_0[38] = x281*(x545*x597 + x546*x598 + x547*x596 - x551*x599);
    res_0[39] = x281*(x372*x600 + x390*x601 + x519*x602 - x536*x603);
    res_0[40] = x281*(x538*x600 + x539*x602 + x540*x601 - x544*x603);
    res_0[41] = x281*(x545*x601 + x546*x602 + x547*x600 - x551*x603);
    res_0[42] = x281*(x372*x604 + x390*x605 + x519*x606 - x536*x607);
    res_0[43] = x281*(x538*x604 + x539*x606 + x540*x605 - x544*x607);
    res_0[44] = x281*(x545*x605 + x546*x606 + x547*x604 - x551*x607);
    res_0[45] = x281*(x372*x608 + x390*x609 + x519*x610 - x536*x611);
    res_0[46] = x281*(x538*x608 + x539*x610 + x540*x609 - x544*x611);
    res_0[47] = x281*(x545*x609 + x546*x610 + x547*x608 - x551*x611);
    res_0[48] = x281*(x372*x612 + x390*x613 + x519*x614 - x536*x615);
    res_0[49] = x281*(x538*x612 + x539*x614 + x540*x613 - x544*x615);
    res_0[50] = x281*(x545*x613 + x546*x614 + x547*x612 - x551*x615);
    res_0[51] = x281*(x372*x616 + x390*x617 + x519*x618 - x536*x619);
    res_0[52] = x281*(x538*x616 + x539*x618 + x540*x617 - x544*x619);
    res_0[53] = x281*(x545*x617 + x546*x618 + x547*x616 - x551*x619);
    res_0[54] = x281*(x372*x620 + x390*x621 + x519*x622 - x536*x623);
    res_0[55] = x281*(x538*x620 + x539*x622 + x540*x621 - x544*x623);
    res_0[56] = x281*(x545*x621 + x546*x622 + x547*x620 - x551*x623);
    res_0[57] = x281*(x372*x624 + x390*x625 + x519*x626 - x536*x627);
    res_0[58] = x281*(x538*x624 + x539*x626 + x540*x625 - x544*x627);
    res_0[59] = x281*(x545*x625 + x546*x626 + x547*x624 - x551*x627);
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
    
    double x0 = s*s;
    double x1 = (1.0/4.0)*x0;
    double x2 = x1 - 1.0/4.0;
    double x3 = (1.0/4.0)*r;
    double x4 = -x3;
    double x5 = x0*x3;
    double x6 = x4 + x5;
    double x7 = x2 + x6;
    double *x8 = coord;
    double x9 = x8[27];
    double x10 = r*r;
    double x11 = (1.0/4.0)*x10;
    double x12 = x11 - 1.0/4.0;
    double x13 = (1.0/4.0)*s;
    double x14 = -x13;
    double x15 = x10*x13;
    double x16 = x14 + x15;
    double x17 = x12 + x16;
    double x18 = x8[30];
    double x19 = x3 - x5;
    double x20 = x19 + x2;
    double x21 = x8[33];
    double x22 = x13 - x15;
    double x23 = x12 + x22;
    double x24 = x8[24];
    double x25 = 1.0/4.0 - x1;
    double x26 = x25 + x6;
    double x27 = x8[45];
    double x28 = 1.0/4.0 - x11;
    double x29 = x16 + x28;
    double x30 = x8[36];
    double x31 = x22 + x28;
    double x32 = x8[42];
    double x33 = x19 + x25;
    double x34 = x8[39];
    double x35 = (1.0/2.0)*t;
    double x36 = -x35;
    double x37 = r*x35;
    double x38 = x36 + x37;
    double x39 = s*x35;
    double x40 = s*x37;
    double x41 = -x40;
    double x42 = x39 + x41;
    double x43 = x38 + x42;
    double x44 = x8[48];
    double x45 = -x39;
    double x46 = x40 + x45;
    double x47 = x38 + x46;
    double x48 = x8[57];
    double x49 = -x37;
    double x50 = x36 + x49;
    double x51 = x39 + x40;
    double x52 = x50 + x51;
    double x53 = x8[51];
    double x54 = x41 + x45;
    double x55 = x50 + x54;
    double x56 = x8[54];
    double x57 = (1.0/4.0)*t;
    double x58 = t*x13;
    double x59 = -x58;
    double x60 = t*x3;
    double x61 = s*x60;
    double x62 = x59 + x61;
    double x63 = x57 + x62;
    double x64 = (1.0/8.0)*x10;
    double x65 = s*x64;
    double x66 = -x65;
    double x67 = (1.0/8.0)*r;
    double x68 = s*x67;
    double x69 = x66 + x68;
    double x70 = (1.0/8.0)*x0;
    double x71 = r*x70;
    double x72 = -x71;
    double x73 = x64 - 1.0/8.0;
    double x74 = -x60;
    double x75 = x70 + x74;
    double x76 = x72 + x73 + x75;
    double x77 = x63 + x69 + x76;
    double x78 = x8[12];
    double x79 = -x68;
    double x80 = x65 + x79;
    double x81 = -x61;
    double x82 = x58 + x81;
    double x83 = x57 + x82;
    double x84 = x76 + x80 + x83;
    double x85 = x8[21];
    double x86 = x59 + x81;
    double x87 = x73 + x86;
    double x88 = x57 + x66 + x79;
    double x89 = x60 + x70;
    double x90 = x71 + x89;
    double x91 = x87 + x88 + x90;
    double x92 = x8[15];
    double x93 = x57 + x65 + x68;
    double x94 = x58 + x61;
    double x95 = x73 + x94;
    double x96 = x90 + x93 + x95;
    double x97 = x8[18];
    double x98 = -x70;
    double x99 = 1.0/8.0 - x64;
    double x100 = x98 + x99;
    double x101 = x100 + x71 + x74;
    double x102 = x101 + x63 + x80;
    double x103 = x8[0];
    double x104 = x101 + x69 + x83;
    double x105 = x8[9];
    double x106 = x100 + x60 + x72;
    double x107 = x106 + x86 + x93;
    double x108 = x8[3];
    double x109 = x106 + x88 + x94;
    double x110 = x8[6];
    double x111 = x102*x103 + x104*x105 + x107*x108 + x109*x110 + x17*x18 + x20*x21 + x23*x24 + x26*x27 + x29*x30 + x31*x32 + x33*x34 + x43*x44 + x47*x48 + x52*x53 + x55*x56 + x7*x9 + x77*x78 + x84*x85 + x91*x92 + x96*x97;
    double x112 = t*t;
    double x113 = (1.0/4.0)*x112;
    double x114 = x113 - 1.0/4.0;
    double x115 = x112*x3;
    double x116 = x115 + x4;
    double x117 = x114 + x116;
    double x118 = x8[53];
    double x119 = -x57;
    double x120 = t*x11;
    double x121 = x119 + x120;
    double x122 = x12 + x121;
    double x123 = x8[38];
    double x124 = -x115 + x3;
    double x125 = x114 + x124;
    double x126 = x8[50];
    double x127 = -x120 + x57;
    double x128 = x12 + x127;
    double x129 = x8[26];
    double x130 = 1.0/4.0 - x113;
    double x131 = x116 + x130;
    double x132 = x8[59];
    double x133 = x121 + x28;
    double x134 = x8[32];
    double x135 = x127 + x28;
    double x136 = x8[44];
    double x137 = x124 + x130;
    double x138 = x8[56];
    double x139 = (1.0/2.0)*s;
    double x140 = -x139;
    double x141 = r*x139;
    double x142 = x140 + x141;
    double x143 = x142 + x42;
    double x144 = x8[35];
    double x145 = x142 + x46;
    double x146 = x8[47];
    double x147 = -x141;
    double x148 = x140 + x147;
    double x149 = x148 + x51;
    double x150 = x8[29];
    double x151 = x148 + x54;
    double x152 = x8[41];
    double x153 = t*x64;
    double x154 = x13 - x153;
    double x155 = (1.0/8.0)*x112;
    double x156 = r*x155;
    double x157 = -x156;
    double x158 = t*x67;
    double x159 = x157 + x158;
    double x160 = s*x3;
    double x161 = -x160;
    double x162 = x155 + x161;
    double x163 = x162 + x73;
    double x164 = x154 + x159 + x163 + x62;
    double x165 = x8[11];
    double x166 = x13 + x153;
    double x167 = -x158;
    double x168 = x157 + x167;
    double x169 = x163 + x166 + x168 + x82;
    double x170 = x8[23];
    double x171 = x156 + x167;
    double x172 = x154 + x160;
    double x173 = x155 + x171 + x172 + x87;
    double x174 = x8[8];
    double x175 = x156 + x158;
    double x176 = x160 + x166;
    double x177 = x155 + x175 + x176 + x95;
    double x178 = x8[20];
    double x179 = -x155;
    double x180 = x179 + x99;
    double x181 = x161 + x180;
    double x182 = x166 + x171 + x181 + x62;
    double x183 = x8[2];
    double x184 = x154 + x175 + x181 + x82;
    double x185 = x8[14];
    double x186 = x159 + x176 + x180 + x86;
    double x187 = x8[5];
    double x188 = x168 + x172 + x180 + x94;
    double x189 = x8[17];
    double x190 = x117*x118 + x122*x123 + x125*x126 + x128*x129 + x131*x132 + x133*x134 + x135*x136 + x137*x138 + x143*x144 + x145*x146 + x149*x150 + x151*x152 + x164*x165 + x169*x170 + x173*x174 + x177*x178 + x182*x183 + x184*x185 + x186*x187 + x188*x189;
    double x191 = x111*x190;
    double x192 = x102*x183 + x104*x165 + x107*x187 + x109*x174 + x118*x52 + x123*x29 + x126*x43 + x129*x23 + x132*x47 + x134*x17 + x136*x31 + x138*x55 + x144*x20 + x146*x26 + x150*x7 + x152*x33 + x170*x84 + x178*x96 + x185*x77 + x189*x91;
    double x193 = x103*x182 + x105*x164 + x108*x186 + x110*x173 + x117*x53 + x122*x30 + x125*x44 + x128*x24 + x131*x48 + x133*x18 + x135*x32 + x137*x56 + x143*x21 + x145*x27 + x149*x9 + x151*x34 + x169*x85 + x177*x97 + x184*x78 + x188*x92;
    double x194 = x192*x193;
    double x195 = x191 - x194;
    double x196 = x112*x13;
    double x197 = x14 + x196;
    double x198 = x114 + x197;
    double x199 = x8[58];
    double x200 = t*x1;
    double x201 = x119 + x200;
    double x202 = x2 + x201;
    double x203 = x8[46];
    double x204 = x13 - x196;
    double x205 = x114 + x204;
    double x206 = x8[49];
    double x207 = -x200 + x57;
    double x208 = x2 + x207;
    double x209 = x8[34];
    double x210 = x130 + x197;
    double x211 = x8[52];
    double x212 = x201 + x25;
    double x213 = x8[28];
    double x214 = x207 + x25;
    double x215 = x8[40];
    double x216 = x130 + x204;
    double x217 = x8[55];
    double x218 = -1.0/2.0*r;
    double x219 = x141 + x218;
    double x220 = x219 + x37 + x41;
    double x221 = x8[25];
    double x222 = x219 + x40 + x49;
    double x223 = x8[37];
    double x224 = x147 + x218;
    double x225 = x224 + x37 + x40;
    double x226 = x8[31];
    double x227 = x224 + x41 + x49;
    double x228 = x8[43];
    double x229 = x3 + x61;
    double x230 = t*x70;
    double x231 = -x230;
    double x232 = (1.0/8.0)*s*t;
    double x233 = x231 + x232;
    double x234 = s*x155;
    double x235 = -x234;
    double x236 = x162 + x235 - 1.0/8.0;
    double x237 = x229 + x233 + x236 + x75;
    double x238 = x8[4];
    double x239 = x230 + x89;
    double x240 = -x232;
    double x241 = x3 + x81;
    double x242 = x240 + x241;
    double x243 = x236 + x239 + x242;
    double x244 = x8[16];
    double x245 = x160 + x231;
    double x246 = x155 + x234 - 1.0/8.0;
    double x247 = x242 + x245 + x246 + x75;
    double x248 = x8[7];
    double x249 = x160 + x232;
    double x250 = x229 + x239 + x246 + x249;
    double x251 = x8[19];
    double x252 = x161 + x179 + x234 + x98 + 1.0/8.0;
    double x253 = x230 + x74;
    double x254 = x229 + x240;
    double x255 = x252 + x253 + x254;
    double x256 = x8[1];
    double x257 = x233 + x241 + x252 + x60;
    double x258 = x8[13];
    double x259 = x179 + x235 + x98 + 1.0/8.0;
    double x260 = x241 + x249 + x253 + x259;
    double x261 = x8[10];
    double x262 = x245 + x254 + x259 + x60;
    double x263 = x8[22];
    double x264 = x198*x199 + x202*x203 + x205*x206 + x208*x209 + x210*x211 + x212*x213 + x214*x215 + x216*x217 + x220*x221 + x222*x223 + x225*x226 + x227*x228 + x237*x238 + x243*x244 + x247*x248 + x250*x251 + x255*x256 + x257*x258 + x260*x261 + x262*x263;
    double x265 = x102*x256 + x104*x261 + x107*x238 + x109*x248 + x17*x226 + x199*x47 + x20*x209 + x203*x26 + x206*x43 + x211*x52 + x213*x7 + x215*x33 + x217*x55 + x221*x23 + x223*x29 + x228*x31 + x244*x91 + x251*x96 + x258*x77 + x263*x84;
    double x266 = x118*x210 + x123*x222 + x126*x205 + x129*x220 + x132*x198 + x134*x225 + x136*x227 + x138*x216 + x144*x208 + x146*x202 + x150*x212 + x152*x214 + x165*x260 + x170*x262 + x174*x247 + x178*x250 + x183*x255 + x185*x257 + x187*x237 + x189*x243;
    double x267 = x193*x266;
    double x268 = x117*x211 + x122*x223 + x125*x206 + x128*x221 + x131*x199 + x133*x226 + x135*x228 + x137*x217 + x143*x209 + x145*x203 + x149*x213 + x151*x215 + x164*x261 + x169*x263 + x173*x248 + x177*x251 + x182*x256 + x184*x258 + x186*x238 + x188*x244;
    double x269 = x103*x255 + x105*x260 + x108*x237 + x110*x247 + x18*x225 + x198*x48 + x202*x27 + x205*x44 + x208*x21 + x210*x53 + x212*x9 + x214*x34 + x216*x56 + x220*x24 + x222*x30 + x227*x32 + x243*x92 + x250*x97 + x257*x78 + x262*x85;
    double x270 = x192*x269;
    double x271 = x111*x266;
    double x272 = x190*x269;
    double x273 = x191*x264 - x194*x264 + x265*x267 - x265*x272 + x268*x270 - x268*x271;
    double x274 = 1.0/x273;
    double *x275 = U;
    double x276 = x275[59];
    double x277 = x275[47];
    double x278 = x275[50];
    double x279 = x275[35];
    double x280 = x275[53];
    double x281 = x275[29];
    double x282 = x275[41];
    double x283 = x275[56];
    double x284 = x275[26];
    double x285 = x275[38];
    double x286 = x275[32];
    double x287 = x275[44];
    double x288 = x275[5];
    double x289 = x275[17];
    double x290 = x275[8];
    double x291 = x275[20];
    double x292 = x275[2];
    double x293 = x275[14];
    double x294 = x275[11];
    double x295 = x275[23];
    double x296 = x274*(x198*x276 + x202*x277 + x205*x278 + x208*x279 + x210*x280 + x212*x281 + x214*x282 + x216*x283 + x220*x284 + x222*x285 + x225*x286 + x227*x287 + x237*x288 + x243*x289 + x247*x290 + x250*x291 + x255*x292 + x257*x293 + x260*x294 + x262*x295);
    double x297 = x270 - x271;
    double x298 = x274*(x117*x280 + x122*x285 + x125*x278 + x128*x284 + x131*x276 + x133*x286 + x135*x287 + x137*x283 + x143*x279 + x145*x277 + x149*x281 + x151*x282 + x164*x294 + x169*x295 + x173*x290 + x177*x291 + x182*x292 + x184*x293 + x186*x288 + x188*x289);
    double x299 = x267 - x272;
    double x300 = x274*(x102*x292 + x104*x294 + x107*x288 + x109*x290 + x17*x286 + x20*x279 + x23*x284 + x26*x277 + x276*x47 + x278*x43 + x280*x52 + x281*x7 + x282*x33 + x283*x55 + x285*x29 + x287*x31 + x289*x91 + x291*x96 + x293*x77 + x295*x84);
    double x301 = x195*x296 + x297*x298 + x299*x300;
    double x302 = -x192*x264 + x265*x266;
    double x303 = x275[52];
    double x304 = x275[37];
    double x305 = x275[49];
    double x306 = x275[25];
    double x307 = x275[58];
    double x308 = x275[31];
    double x309 = x275[43];
    double x310 = x275[55];
    double x311 = x275[34];
    double x312 = x275[46];
    double x313 = x275[28];
    double x314 = x275[40];
    double x315 = x275[10];
    double x316 = x275[22];
    double x317 = x275[7];
    double x318 = x275[19];
    double x319 = x275[1];
    double x320 = x275[13];
    double x321 = x275[4];
    double x322 = x275[16];
    double x323 = x274*(x117*x303 + x122*x304 + x125*x305 + x128*x306 + x131*x307 + x133*x308 + x135*x309 + x137*x310 + x143*x311 + x145*x312 + x149*x313 + x151*x314 + x164*x315 + x169*x316 + x173*x317 + x177*x318 + x182*x319 + x184*x320 + x186*x321 + x188*x322);
    double x324 = -x190*x265 + x192*x268;
    double x325 = x274*(x198*x307 + x202*x312 + x205*x305 + x208*x311 + x210*x303 + x212*x313 + x214*x314 + x216*x310 + x220*x306 + x222*x304 + x225*x308 + x227*x309 + x237*x321 + x243*x322 + x247*x317 + x250*x318 + x255*x319 + x257*x320 + x260*x315 + x262*x316);
    double x326 = x190*x264 - x266*x268;
    double x327 = x274*(x102*x319 + x104*x315 + x107*x321 + x109*x317 + x17*x308 + x20*x311 + x23*x306 + x26*x312 + x29*x304 + x303*x52 + x305*x43 + x307*x47 + x309*x31 + x310*x55 + x313*x7 + x314*x33 + x316*x84 + x318*x96 + x320*x77 + x322*x91);
    double x328 = x302*x323 + x324*x325 + x326*x327;
    double x329 = x296*x324 + x298*x302 + x300*x326;
    double x330 = x195*x325 + x297*x323 + x299*x327 + 1.0;
    double x331 = x301*x328 - x329*x330;
    double x332 = x275[57];
    double x333 = x275[45];
    double x334 = x275[48];
    double x335 = x275[33];
    double x336 = x275[51];
    double x337 = x275[27];
    double x338 = x275[39];
    double x339 = x275[54];
    double x340 = x275[24];
    double x341 = x275[36];
    double x342 = x275[30];
    double x343 = x275[42];
    double x344 = x275[3];
    double x345 = x275[15];
    double x346 = x275[6];
    double x347 = x275[18];
    double x348 = x275[0];
    double x349 = x275[12];
    double x350 = x275[9];
    double x351 = x275[21];
    double x352 = x274*(x198*x332 + x202*x333 + x205*x334 + x208*x335 + x210*x336 + x212*x337 + x214*x338 + x216*x339 + x220*x340 + x222*x341 + x225*x342 + x227*x343 + x237*x344 + x243*x345 + x247*x346 + x250*x347 + x255*x348 + x257*x349 + x260*x350 + x262*x351);
    double x353 = x274*(x117*x336 + x122*x341 + x125*x334 + x128*x340 + x131*x332 + x133*x342 + x135*x343 + x137*x339 + x143*x335 + x145*x333 + x149*x337 + x151*x338 + x164*x350 + x169*x351 + x173*x346 + x177*x347 + x182*x348 + x184*x349 + x186*x344 + x188*x345);
    double x354 = x274*(x102*x348 + x104*x350 + x107*x344 + x109*x346 + x17*x342 + x20*x335 + x23*x340 + x26*x333 + x29*x341 + x31*x343 + x33*x338 + x332*x47 + x334*x43 + x336*x52 + x337*x7 + x339*x55 + x345*x91 + x347*x96 + x349*x77 + x351*x84);
    double x355 = x195*x352 + x297*x353 + x299*x354;
    double x356 = x302*x353 + x324*x352 + x326*x354 + 1.0;
    double x357 = -x301*x356 + x329*x355;
    double x358 = -x328*x355 + x330*x356;
    double x359 = S11*x331 + S12*x357 + S13*x358;
    double x360 = S12*x331 + S22*x357 + S23*x358;
    double x361 = S13*x331 + S23*x357 + S33*x358;
    double x362 = -x328*x360 - x329*x361 - x356*x359;
    double x363 = x111*x264 - x265*x269;
    double x364 = x182*x274;
    double x365 = -x111*x268 + x193*x265;
    double x366 = x255*x274;
    double x367 = -x193*x264 + x268*x269;
    double x368 = x102*x274;
    double x369 = x363*x364 + x365*x366 + x367*x368;
    double x370 = x323*x363 + x325*x365 + x327*x367;
    double x371 = x296*x365 + x298*x363 + x300*x367 + 1.0;
    double x372 = -x328*x371 + x329*x370;
    double x373 = x352*x365 + x353*x363 + x354*x367;
    double x374 = -x329*x373 + x356*x371;
    double x375 = x328*x373 - x356*x370;
    double x376 = S11*x372 + S12*x374 + S13*x375;
    double x377 = S12*x372 + S22*x374 + S23*x375;
    double x378 = S13*x372 + S23*x374 + S33*x375;
    double x379 = -x328*x377 - x329*x378 - x356*x376;
    double x380 = x195*x366 + x297*x364 + x299*x368;
    double x381 = x302*x364 + x324*x366 + x326*x368;
    double x382 = -x301*x370 + x330*x371;
    double x383 = x301*x373 - x355*x371;
    double x384 = -x330*x373 + x355*x370;
    double x385 = S11*x382 + S12*x383 + S13*x384;
    double x386 = S12*x382 + S22*x383 + S23*x384;
    double x387 = S13*x382 + S23*x383 + S33*x384;
    double x388 = 1.0*PENER + 1.0*SENER;
    double x389 = -x328*x386 - x329*x387 - x356*x385 + x388;
    double x390 = -x301*x361 - x330*x360 - x355*x359;
    double x391 = -x301*x387 - x330*x386 - x355*x385;
    double x392 = -x301*x378 - x330*x377 - x355*x376 + x388;
    double x393 = -x370*x377 - x371*x378 - x373*x376;
    double x394 = -x370*x386 - x371*x387 - x373*x385;
    double x395 = -x359*x373 - x360*x370 - x361*x371 + x388;
    double x396 = x186*x274;
    double x397 = x237*x274;
    double x398 = x107*x274;
    double x399 = x363*x396 + x365*x397 + x367*x398;
    double x400 = x195*x397 + x297*x396 + x299*x398;
    double x401 = x302*x396 + x324*x397 + x326*x398;
    double x402 = x173*x274;
    double x403 = x247*x274;
    double x404 = x109*x274;
    double x405 = x363*x402 + x365*x403 + x367*x404;
    double x406 = x195*x403 + x297*x402 + x299*x404;
    double x407 = x302*x402 + x324*x403 + x326*x404;
    double x408 = x164*x274;
    double x409 = x260*x274;
    double x410 = x104*x274;
    double x411 = x363*x408 + x365*x409 + x367*x410;
    double x412 = x195*x409 + x297*x408 + x299*x410;
    double x413 = x302*x408 + x324*x409 + x326*x410;
    double x414 = x184*x274;
    double x415 = x257*x274;
    double x416 = x274*x77;
    double x417 = x363*x414 + x365*x415 + x367*x416;
    double x418 = x195*x415 + x297*x414 + x299*x416;
    double x419 = x302*x414 + x324*x415 + x326*x416;
    double x420 = x188*x274;
    double x421 = x243*x274;
    double x422 = x274*x91;
    double x423 = x363*x420 + x365*x421 + x367*x422;
    double x424 = x195*x421 + x297*x420 + x299*x422;
    double x425 = x302*x420 + x324*x421 + x326*x422;
    double x426 = x177*x274;
    double x427 = x250*x274;
    double x428 = x274*x96;
    double x429 = x363*x426 + x365*x427 + x367*x428;
    double x430 = x195*x427 + x297*x426 + x299*x428;
    double x431 = x302*x426 + x324*x427 + x326*x428;
    double x432 = x169*x274;
    double x433 = x262*x274;
    double x434 = x274*x84;
    double x435 = x363*x432 + x365*x433 + x367*x434;
    double x436 = x195*x433 + x297*x432 + x299*x434;
    double x437 = x302*x432 + x324*x433 + x326*x434;
    double x438 = x128*x274;
    double x439 = x220*x274;
    double x440 = x23*x274;
    double x441 = x363*x438 + x365*x439 + x367*x440;
    double x442 = x195*x439 + x297*x438 + x299*x440;
    double x443 = x302*x438 + x324*x439 + x326*x440;
    double x444 = x149*x274;
    double x445 = x212*x274;
    double x446 = x274*x7;
    double x447 = x363*x444 + x365*x445 + x367*x446;
    double x448 = x195*x445 + x297*x444 + x299*x446;
    double x449 = x302*x444 + x324*x445 + x326*x446;
    double x450 = x133*x274;
    double x451 = x225*x274;
    double x452 = x17*x274;
    double x453 = x363*x450 + x365*x451 + x367*x452;
    double x454 = x195*x451 + x297*x450 + x299*x452;
    double x455 = x302*x450 + x324*x451 + x326*x452;
    double x456 = x143*x274;
    double x457 = x208*x274;
    double x458 = x20*x274;
    double x459 = x363*x456 + x365*x457 + x367*x458;
    double x460 = x195*x457 + x297*x456 + x299*x458;
    double x461 = x302*x456 + x324*x457 + x326*x458;
    double x462 = x122*x274;
    double x463 = x222*x274;
    double x464 = x274*x29;
    double x465 = x363*x462 + x365*x463 + x367*x464;
    double x466 = x195*x463 + x297*x462 + x299*x464;
    double x467 = x302*x462 + x324*x463 + x326*x464;
    double x468 = x151*x274;
    double x469 = x214*x274;
    double x470 = x274*x33;
    double x471 = x363*x468 + x365*x469 + x367*x470;
    double x472 = x195*x469 + x297*x468 + x299*x470;
    double x473 = x302*x468 + x324*x469 + x326*x470;
    double x474 = x135*x274;
    double x475 = x227*x274;
    double x476 = x274*x31;
    double x477 = x363*x474 + x365*x475 + x367*x476;
    double x478 = x195*x475 + x297*x474 + x299*x476;
    double x479 = x302*x474 + x324*x475 + x326*x476;
    double x480 = x145*x274;
    double x481 = x202*x274;
    double x482 = x26*x274;
    double x483 = x363*x480 + x365*x481 + x367*x482;
    double x484 = x195*x481 + x297*x480 + x299*x482;
    double x485 = x302*x480 + x324*x481 + x326*x482;
    double x486 = x125*x274;
    double x487 = x205*x274;
    double x488 = x274*x43;
    double x489 = x363*x486 + x365*x487 + x367*x488;
    double x490 = x195*x487 + x297*x486 + x299*x488;
    double x491 = x302*x486 + x324*x487 + x326*x488;
    double x492 = x117*x274;
    double x493 = x210*x274;
    double x494 = x274*x52;
    double x495 = x363*x492 + x365*x493 + x367*x494;
    double x496 = x195*x493 + x297*x492 + x299*x494;
    double x497 = x302*x492 + x324*x493 + x326*x494;
    double x498 = x137*x274;
    double x499 = x216*x274;
    double x500 = x274*x55;
    double x501 = x363*x498 + x365*x499 + x367*x500;
    double x502 = x195*x499 + x297*x498 + x299*x500;
    double x503 = x302*x498 + x324*x499 + x326*x500;
    double x504 = x131*x274;
    double x505 = x198*x274;
    double x506 = x274*x47;
    double x507 = x363*x504 + x365*x505 + x367*x506;
    double x508 = x195*x505 + x297*x504 + x299*x506;
    double x509 = x302*x504 + x324*x505 + x326*x506;
    
    res_0[0] = x273*(x362*x369 + x379*x380 + x381*x389);
    res_0[1] = x273*(x369*x390 + x380*x392 + x381*x391);
    res_0[2] = x273*(x369*x395 + x380*x393 + x381*x394);
    res_0[3] = x273*(x362*x399 + x379*x400 + x389*x401);
    res_0[4] = x273*(x390*x399 + x391*x401 + x392*x400);
    res_0[5] = x273*(x393*x400 + x394*x401 + x395*x399);
    res_0[6] = x273*(x362*x405 + x379*x406 + x389*x407);
    res_0[7] = x273*(x390*x405 + x391*x407 + x392*x406);
    res_0[8] = x273*(x393*x406 + x394*x407 + x395*x405);
    res_0[9] = x273*(x362*x411 + x379*x412 + x389*x413);
    res_0[10] = x273*(x390*x411 + x391*x413 + x392*x412);
    res_0[11] = x273*(x393*x412 + x394*x413 + x395*x411);
    res_0[12] = x273*(x362*x417 + x379*x418 + x389*x419);
    res_0[13] = x273*(x390*x417 + x391*x419 + x392*x418);
    res_0[14] = x273*(x393*x418 + x394*x419 + x395*x417);
    res_0[15] = x273*(x362*x423 + x379*x424 + x389*x425);
    res_0[16] = x273*(x390*x423 + x391*x425 + x392*x424);
    res_0[17] = x273*(x393*x424 + x394*x425 + x395*x423);
    res_0[18] = x273*(x362*x429 + x379*x430 + x389*x431);
    res_0[19] = x273*(x390*x429 + x391*x431 + x392*x430);
    res_0[20] = x273*(x393*x430 + x394*x431 + x395*x429);
    res_0[21] = x273*(x362*x435 + x379*x436 + x389*x437);
    res_0[22] = x273*(x390*x435 + x391*x437 + x392*x436);
    res_0[23] = x273*(x393*x436 + x394*x437 + x395*x435);
    res_0[24] = x273*(x362*x441 + x379*x442 + x389*x443);
    res_0[25] = x273*(x390*x441 + x391*x443 + x392*x442);
    res_0[26] = x273*(x393*x442 + x394*x443 + x395*x441);
    res_0[27] = x273*(x362*x447 + x379*x448 + x389*x449);
    res_0[28] = x273*(x390*x447 + x391*x449 + x392*x448);
    res_0[29] = x273*(x393*x448 + x394*x449 + x395*x447);
    res_0[30] = x273*(x362*x453 + x379*x454 + x389*x455);
    res_0[31] = x273*(x390*x453 + x391*x455 + x392*x454);
    res_0[32] = x273*(x393*x454 + x394*x455 + x395*x453);
    res_0[33] = x273*(x362*x459 + x379*x460 + x389*x461);
    res_0[34] = x273*(x390*x459 + x391*x461 + x392*x460);
    res_0[35] = x273*(x393*x460 + x394*x461 + x395*x459);
    res_0[36] = x273*(x362*x465 + x379*x466 + x389*x467);
    res_0[37] = x273*(x390*x465 + x391*x467 + x392*x466);
    res_0[38] = x273*(x393*x466 + x394*x467 + x395*x465);
    res_0[39] = x273*(x362*x471 + x379*x472 + x389*x473);
    res_0[40] = x273*(x390*x471 + x391*x473 + x392*x472);
    res_0[41] = x273*(x393*x472 + x394*x473 + x395*x471);
    res_0[42] = x273*(x362*x477 + x379*x478 + x389*x479);
    res_0[43] = x273*(x390*x477 + x391*x479 + x392*x478);
    res_0[44] = x273*(x393*x478 + x394*x479 + x395*x477);
    res_0[45] = x273*(x362*x483 + x379*x484 + x389*x485);
    res_0[46] = x273*(x390*x483 + x391*x485 + x392*x484);
    res_0[47] = x273*(x393*x484 + x394*x485 + x395*x483);
    res_0[48] = x273*(x362*x489 + x379*x490 + x389*x491);
    res_0[49] = x273*(x390*x489 + x391*x491 + x392*x490);
    res_0[50] = x273*(x393*x490 + x394*x491 + x395*x489);
    res_0[51] = x273*(x362*x495 + x379*x496 + x389*x497);
    res_0[52] = x273*(x390*x495 + x391*x497 + x392*x496);
    res_0[53] = x273*(x393*x496 + x394*x497 + x395*x495);
    res_0[54] = x273*(x362*x501 + x379*x502 + x389*x503);
    res_0[55] = x273*(x390*x501 + x391*x503 + x392*x502);
    res_0[56] = x273*(x393*x502 + x394*x503 + x395*x501);
    res_0[57] = x273*(x362*x507 + x379*x508 + x389*x509);
    res_0[58] = x273*(x390*x507 + x391*x509 + x392*x508);
    res_0[59] = x273*(x393*x508 + x394*x509 + x395*x507);
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
    
    double x0 = s*s;
    double x1 = (1.0/4.0)*x0;
    double x2 = x1 - 1.0/4.0;
    double x3 = (1.0/4.0)*r;
    double x4 = -x3;
    double x5 = x0*x3;
    double x6 = x4 + x5;
    double x7 = x2 + x6;
    double *x8 = coord;
    double x9 = x8[27];
    double x10 = r*r;
    double x11 = (1.0/4.0)*x10;
    double x12 = x11 - 1.0/4.0;
    double x13 = (1.0/4.0)*s;
    double x14 = -x13;
    double x15 = x10*x13;
    double x16 = x14 + x15;
    double x17 = x12 + x16;
    double x18 = x8[30];
    double x19 = x3 - x5;
    double x20 = x19 + x2;
    double x21 = x8[33];
    double x22 = x13 - x15;
    double x23 = x12 + x22;
    double x24 = x8[24];
    double x25 = 1.0/4.0 - x1;
    double x26 = x25 + x6;
    double x27 = x8[45];
    double x28 = 1.0/4.0 - x11;
    double x29 = x16 + x28;
    double x30 = x8[36];
    double x31 = x22 + x28;
    double x32 = x8[42];
    double x33 = x19 + x25;
    double x34 = x8[39];
    double x35 = (1.0/2.0)*t;
    double x36 = -x35;
    double x37 = r*x35;
    double x38 = x36 + x37;
    double x39 = s*x35;
    double x40 = s*x37;
    double x41 = -x40;
    double x42 = x39 + x41;
    double x43 = x38 + x42;
    double x44 = x8[48];
    double x45 = -x39;
    double x46 = x40 + x45;
    double x47 = x38 + x46;
    double x48 = x8[57];
    double x49 = -x37;
    double x50 = x36 + x49;
    double x51 = x39 + x40;
    double x52 = x50 + x51;
    double x53 = x8[51];
    double x54 = x41 + x45;
    double x55 = x50 + x54;
    double x56 = x8[54];
    double x57 = (1.0/4.0)*t;
    double x58 = t*x13;
    double x59 = -x58;
    double x60 = t*x3;
    double x61 = s*x60;
    double x62 = x59 + x61;
    double x63 = x57 + x62;
    double x64 = (1.0/8.0)*x10;
    double x65 = s*x64;
    double x66 = -x65;
    double x67 = (1.0/8.0)*r;
    double x68 = s*x67;
    double x69 = x66 + x68;
    double x70 = (1.0/8.0)*x0;
    double x71 = r*x70;
    double x72 = -x71;
    double x73 = x64 - 1.0/8.0;
    double x74 = -x60;
    double x75 = x70 + x74;
    double x76 = x72 + x73 + x75;
    double x77 = x63 + x69 + x76;
    double x78 = x8[12];
    double x79 = -x68;
    double x80 = x65 + x79;
    double x81 = -x61;
    double x82 = x58 + x81;
    double x83 = x57 + x82;
    double x84 = x76 + x80 + x83;
    double x85 = x8[21];
    double x86 = x59 + x81;
    double x87 = x73 + x86;
    double x88 = x57 + x66 + x79;
    double x89 = x60 + x70;
    double x90 = x71 + x89;
    double x91 = x87 + x88 + x90;
    double x92 = x8[15];
    double x93 = x57 + x65 + x68;
    double x94 = x58 + x61;
    double x95 = x73 + x94;
    double x96 = x90 + x93 + x95;
    double x97 = x8[18];
    double x98 = -x70;
    double x99 = 1.0/8.0 - x64;
    double x100 = x98 + x99;
    double x101 = x100 + x71 + x74;
    double x102 = x101 + x63 + x80;
    double x103 = x8[0];
    double x104 = x101 + x69 + x83;
    double x105 = x8[9];
    double x106 = x100 + x60 + x72;
    double x107 = x106 + x86 + x93;
    double x108 = x8[3];
    double x109 = x106 + x88 + x94;
    double x110 = x8[6];
    double x111 = x102*x103 + x104*x105 + x107*x108 + x109*x110 + x17*x18 + x20*x21 + x23*x24 + x26*x27 + x29*x30 + x31*x32 + x33*x34 + x43*x44 + x47*x48 + x52*x53 + x55*x56 + x7*x9 + x77*x78 + x84*x85 + x91*x92 + x96*x97;
    double x112 = t*t;
    double x113 = (1.0/4.0)*x112;
    double x114 = x113 - 1.0/4.0;
    double x115 = x112*x3;
    double x116 = x115 + x4;
    double x117 = x114 + x116;
    double x118 = x8[53];
    double x119 = -x57;
    double x120 = t*x11;
    double x121 = x119 + x120;
    double x122 = x12 + x121;
    double x123 = x8[38];
    double x124 = -x115 + x3;
    double x125 = x114 + x124;
    double x126 = x8[50];
    double x127 = -x120 + x57;
    double x128 = x12 + x127;
    double x129 = x8[26];
    double x130 = 1.0/4.0 - x113;
    double x131 = x116 + x130;
    double x132 = x8[59];
    double x133 = x121 + x28;
    double x134 = x8[32];
    double x135 = x127 + x28;
    double x136 = x8[44];
    double x137 = x124 + x130;
    double x138 = x8[56];
    double x139 = (1.0/2.0)*s;
    double x140 = -x139;
    double x141 = r*x139;
    double x142 = x140 + x141;
    double x143 = x142 + x42;
    double x144 = x8[35];
    double x145 = x142 + x46;
    double x146 = x8[47];
    double x147 = -x141;
    double x148 = x140 + x147;
    double x149 = x148 + x51;
    double x150 = x8[29];
    double x151 = x148 + x54;
    double x152 = x8[41];
    double x153 = t*x64;
    double x154 = x13 - x153;
    double x155 = (1.0/8.0)*x112;
    double x156 = r*x155;
    double x157 = -x156;
    double x158 = t*x67;
    double x159 = x157 + x158;
    double x160 = s*x3;
    double x161 = -x160;
    double x162 = x155 + x161;
    double x163 = x162 + x73;
    double x164 = x154 + x159 + x163 + x62;
    double x165 = x8[11];
    double x166 = x13 + x153;
    double x167 = -x158;
    double x168 = x157 + x167;
    double x169 = x163 + x166 + x168 + x82;
    double x170 = x8[23];
    double x171 = x156 + x167;
    double x172 = x154 + x160;
    double x173 = x155 + x171 + x172 + x87;
    double x174 = x8[8];
    double x175 = x156 + x158;
    double x176 = x160 + x166;
    double x177 = x155 + x175 + x176 + x95;
    double x178 = x8[20];
    double x179 = -x155;
    double x180 = x179 + x99;
    double x181 = x161 + x180;
    double x182 = x166 + x171 + x181 + x62;
    double x183 = x8[2];
    double x184 = x154 + x175 + x181 + x82;
    double x185 = x8[14];
    double x186 = x159 + x176 + x180 + x86;
    double x187 = x8[5];
    double x188 = x168 + x172 + x180 + x94;
    double x189 = x8[17];
    double x190 = x117*x118 + x122*x123 + x125*x126 + x128*x129 + x131*x132 + x133*x134 + x135*x136 + x137*x138 + x143*x144 + x145*x146 + x149*x150 + x151*x152 + x164*x165 + x169*x170 + x173*x174 + x177*x178 + x182*x183 + x184*x185 + x186*x187 + x188*x189;
    double x191 = x111*x190;
    double x192 = x102*x183 + x104*x165 + x107*x187 + x109*x174 + x118*x52 + x123*x29 + x126*x43 + x129*x23 + x132*x47 + x134*x17 + x136*x31 + x138*x55 + x144*x20 + x146*x26 + x150*x7 + x152*x33 + x170*x84 + x178*x96 + x185*x77 + x189*x91;
    double x193 = x103*x182 + x105*x164 + x108*x186 + x110*x173 + x117*x53 + x122*x30 + x125*x44 + x128*x24 + x131*x48 + x133*x18 + x135*x32 + x137*x56 + x143*x21 + x145*x27 + x149*x9 + x151*x34 + x169*x85 + x177*x97 + x184*x78 + x188*x92;
    double x194 = x192*x193;
    double x195 = x191 - x194;
    double x196 = x112*x13;
    double x197 = x14 + x196;
    double x198 = x114 + x197;
    double x199 = x8[58];
    double x200 = t*x1;
    double x201 = x119 + x200;
    double x202 = x2 + x201;
    double x203 = x8[46];
    double x204 = x13 - x196;
    double x205 = x114 + x204;
    double x206 = x8[49];
    double x207 = -x200 + x57;
    double x208 = x2 + x207;
    double x209 = x8[34];
    double x210 = x130 + x197;
    double x211 = x8[52];
    double x212 = x201 + x25;
    double x213 = x8[28];
    double x214 = x207 + x25;
    double x215 = x8[40];
    double x216 = x130 + x204;
    double x217 = x8[55];
    double x218 = -1.0/2.0*r;
    double x219 = x141 + x218;
    double x220 = x219 + x37 + x41;
    double x221 = x8[25];
    double x222 = x219 + x40 + x49;
    double x223 = x8[37];
    double x224 = x147 + x218;
    double x225 = x224 + x37 + x40;
    double x226 = x8[31];
    double x227 = x224 + x41 + x49;
    double x228 = x8[43];
    double x229 = x3 + x61;
    double x230 = t*x70;
    double x231 = -x230;
    double x232 = (1.0/8.0)*s*t;
    double x233 = x231 + x232;
    double x234 = s*x155;
    double x235 = -x234;
    double x236 = x162 + x235 - 1.0/8.0;
    double x237 = x229 + x233 + x236 + x75;
    double x238 = x8[4];
    double x239 = x230 + x89;
    double x240 = -x232;
    double x241 = x3 + x81;
    double x242 = x240 + x241;
    double x243 = x236 + x239 + x242;
    double x244 = x8[16];
    double x245 = x160 + x231;
    double x246 = x155 + x234 - 1.0/8.0;
    double x247 = x242 + x245 + x246 + x75;
    double x248 = x8[7];
    double x249 = x160 + x232;
    double x250 = x229 + x239 + x246 + x249;
    double x251 = x8[19];
    double x252 = x161 + x179 + x234 + x98 + 1.0/8.0;
    double x253 = x230 + x74;
    double x254 = x229 + x240;
    double x255 = x252 + x253 + x254;
    double x256 = x8[1];
    double x257 = x233 + x241 + x252 + x60;
    double x258 = x8[13];
    double x259 = x179 + x235 + x98 + 1.0/8.0;
    double x260 = x241 + x249 + x253 + x259;
    double x261 = x8[10];
    double x262 = x245 + x254 + x259 + x60;
    double x263 = x8[22];
    double x264 = x198*x199 + x202*x203 + x205*x206 + x208*x209 + x210*x211 + x212*x213 + x214*x215 + x216*x217 + x220*x221 + x222*x223 + x225*x226 + x227*x228 + x237*x238 + x243*x244 + x247*x248 + x250*x251 + x255*x256 + x257*x258 + x260*x261 + x262*x263;
    double x265 = x102*x256 + x104*x261 + x107*x238 + x109*x248 + x17*x226 + x199*x47 + x20*x209 + x203*x26 + x206*x43 + x211*x52 + x213*x7 + x215*x33 + x217*x55 + x221*x23 + x223*x29 + x228*x31 + x244*x91 + x251*x96 + x258*x77 + x263*x84;
    double x266 = x118*x210 + x123*x222 + x126*x205 + x129*x220 + x132*x198 + x134*x225 + x136*x227 + x138*x216 + x144*x208 + x146*x202 + x150*x212 + x152*x214 + x165*x260 + x170*x262 + x174*x247 + x178*x250 + x183*x255 + x185*x257 + x187*x237 + x189*x243;
    double x267 = x193*x266;
    double x268 = x117*x211 + x122*x223 + x125*x206 + x128*x221 + x131*x199 + x133*x226 + x135*x228 + x137*x217 + x143*x209 + x145*x203 + x149*x213 + x151*x215 + x164*x261 + x169*x263 + x173*x248 + x177*x251 + x182*x256 + x184*x258 + x186*x238 + x188*x244;
    double x269 = x103*x255 + x105*x260 + x108*x237 + x110*x247 + x18*x225 + x198*x48 + x202*x27 + x205*x44 + x208*x21 + x210*x53 + x212*x9 + x214*x34 + x216*x56 + x220*x24 + x222*x30 + x227*x32 + x243*x92 + x250*x97 + x257*x78 + x262*x85;
    double x270 = x192*x269;
    double x271 = x111*x266;
    double x272 = x190*x269;
    double x273 = x191*x264 - x194*x264 + x265*x267 - x265*x272 + x268*x270 - x268*x271;
    double x274 = 1.0/x273;
    double *x275 = U;
    double x276 = x275[59];
    double x277 = x275[47];
    double x278 = x275[50];
    double x279 = x275[35];
    double x280 = x275[53];
    double x281 = x275[29];
    double x282 = x275[41];
    double x283 = x275[56];
    double x284 = x275[26];
    double x285 = x275[38];
    double x286 = x275[32];
    double x287 = x275[44];
    double x288 = x275[5];
    double x289 = x275[17];
    double x290 = x275[8];
    double x291 = x275[20];
    double x292 = x275[2];
    double x293 = x275[14];
    double x294 = x275[11];
    double x295 = x275[23];
    double x296 = x274*(x198*x276 + x202*x277 + x205*x278 + x208*x279 + x210*x280 + x212*x281 + x214*x282 + x216*x283 + x220*x284 + x222*x285 + x225*x286 + x227*x287 + x237*x288 + x243*x289 + x247*x290 + x250*x291 + x255*x292 + x257*x293 + x260*x294 + x262*x295);
    double x297 = x270 - x271;
    double x298 = x274*(x117*x280 + x122*x285 + x125*x278 + x128*x284 + x131*x276 + x133*x286 + x135*x287 + x137*x283 + x143*x279 + x145*x277 + x149*x281 + x151*x282 + x164*x294 + x169*x295 + x173*x290 + x177*x291 + x182*x292 + x184*x293 + x186*x288 + x188*x289);
    double x299 = x267 - x272;
    double x300 = x274*(x102*x292 + x104*x294 + x107*x288 + x109*x290 + x17*x286 + x20*x279 + x23*x284 + x26*x277 + x276*x47 + x278*x43 + x280*x52 + x281*x7 + x282*x33 + x283*x55 + x285*x29 + x287*x31 + x289*x91 + x291*x96 + x293*x77 + x295*x84);
    double x301 = x195*x296 + x297*x298 + x299*x300;
    double x302 = -x192*x264 + x265*x266;
    double x303 = x275[52];
    double x304 = x275[37];
    double x305 = x275[49];
    double x306 = x275[25];
    double x307 = x275[58];
    double x308 = x275[31];
    double x309 = x275[43];
    double x310 = x275[55];
    double x311 = x275[34];
    double x312 = x275[46];
    double x313 = x275[28];
    double x314 = x275[40];
    double x315 = x275[10];
    double x316 = x275[22];
    double x317 = x275[7];
    double x318 = x275[19];
    double x319 = x275[1];
    double x320 = x275[13];
    double x321 = x275[4];
    double x322 = x275[16];
    double x323 = x274*(x117*x303 + x122*x304 + x125*x305 + x128*x306 + x131*x307 + x133*x308 + x135*x309 + x137*x310 + x143*x311 + x145*x312 + x149*x313 + x151*x314 + x164*x315 + x169*x316 + x173*x317 + x177*x318 + x182*x319 + x184*x320 + x186*x321 + x188*x322);
    double x324 = -x190*x265 + x192*x268;
    double x325 = x274*(x198*x307 + x202*x312 + x205*x305 + x208*x311 + x210*x303 + x212*x313 + x214*x314 + x216*x310 + x220*x306 + x222*x304 + x225*x308 + x227*x309 + x237*x321 + x243*x322 + x247*x317 + x250*x318 + x255*x319 + x257*x320 + x260*x315 + x262*x316);
    double x326 = x190*x264 - x266*x268;
    double x327 = x274*(x102*x319 + x104*x315 + x107*x321 + x109*x317 + x17*x308 + x20*x311 + x23*x306 + x26*x312 + x29*x304 + x303*x52 + x305*x43 + x307*x47 + x309*x31 + x310*x55 + x313*x7 + x314*x33 + x316*x84 + x318*x96 + x320*x77 + x322*x91);
    double x328 = x302*x323 + x324*x325 + x326*x327;
    double x329 = x296*x324 + x298*x302 + x300*x326;
    double x330 = x195*x325 + x297*x323 + x299*x327;
    double x331 = x330 + 1.0;
    double x332 = x301*x328 - x329*x331;
    double x333 = x275[57];
    double x334 = x275[45];
    double x335 = x275[48];
    double x336 = x275[33];
    double x337 = x275[51];
    double x338 = x275[27];
    double x339 = x275[39];
    double x340 = x275[54];
    double x341 = x275[24];
    double x342 = x275[36];
    double x343 = x275[30];
    double x344 = x275[42];
    double x345 = x275[3];
    double x346 = x275[15];
    double x347 = x275[6];
    double x348 = x275[18];
    double x349 = x275[0];
    double x350 = x275[12];
    double x351 = x275[9];
    double x352 = x275[21];
    double x353 = x274*(x198*x333 + x202*x334 + x205*x335 + x208*x336 + x210*x337 + x212*x338 + x214*x339 + x216*x340 + x220*x341 + x222*x342 + x225*x343 + x227*x344 + x237*x345 + x243*x346 + x247*x347 + x250*x348 + x255*x349 + x257*x350 + x260*x351 + x262*x352);
    double x354 = x274*(x117*x337 + x122*x342 + x125*x335 + x128*x341 + x131*x333 + x133*x343 + x135*x344 + x137*x340 + x143*x336 + x145*x334 + x149*x338 + x151*x339 + x164*x351 + x169*x352 + x173*x347 + x177*x348 + x182*x349 + x184*x350 + x186*x345 + x188*x346);
    double x355 = x274*(x102*x349 + x104*x351 + x107*x345 + x109*x347 + x17*x343 + x20*x336 + x23*x341 + x26*x334 + x29*x342 + x31*x344 + x33*x339 + x333*x47 + x335*x43 + x337*x52 + x338*x7 + x340*x55 + x346*x91 + x348*x96 + x350*x77 + x352*x84);
    double x356 = x195*x353 + x297*x354 + x299*x355;
    double x357 = x302*x354 + x324*x353 + x326*x355;
    double x358 = x357 + 1.0;
    double x359 = -x301*x358 + x329*x356;
    double x360 = -x328*x356 + x331*x358;
    double x361 = S11*x332 + S12*x359 + S13*x360;
    double x362 = S12*x332 + S22*x359 + S23*x360;
    double x363 = S13*x332 + S23*x359 + S33*x360;
    double x364 = -x328*x362 - x329*x363 - x357*x361;
    double x365 = x111*x264 - x265*x269;
    double x366 = x182*x274;
    double x367 = -x111*x268 + x193*x265;
    double x368 = x255*x274;
    double x369 = -x193*x264 + x268*x269;
    double x370 = x102*x274;
    double x371 = x365*x366 + x367*x368 + x369*x370;
    double x372 = x323*x365 + x325*x367 + x327*x369;
    double x373 = x296*x367 + x298*x365 + x300*x369;
    double x374 = x373 + 1.0;
    double x375 = -x328*x374 + x329*x372;
    double x376 = x353*x367 + x354*x365 + x355*x369;
    double x377 = -x329*x376 + x358*x374;
    double x378 = x328*x376 - x358*x372;
    double x379 = S11*x375 + S12*x377 + S13*x378;
    double x380 = S12*x375 + S22*x377 + S23*x378;
    double x381 = S13*x375 + S23*x377 + S33*x378;
    double x382 = -x328*x380 - x329*x381 - x357*x379;
    double x383 = x195*x368 + x297*x366 + x299*x370;
    double x384 = x302*x366 + x324*x368 + x326*x370;
    double x385 = -x301*x372 + x331*x374;
    double x386 = x301*x376 - x356*x374;
    double x387 = -x331*x376 + x356*x372;
    double x388 = S11*x385 + S12*x386 + S13*x387;
    double x389 = S12*x385 + S22*x386 + S23*x387;
    double x390 = S13*x385 + S23*x386 + S33*x387;
    double x391 = 1.0*PENER + 1.0*SENER;
    double x392 = -x328*x389 - x329*x390 - x357*x388 + x391;
    double x393 = -x301*x363 - x330*x362 - x356*x361;
    double x394 = -x301*x390 - x330*x389 - x356*x388;
    double x395 = -x301*x381 - x330*x380 - x356*x379 + x391;
    double x396 = -x372*x380 - x373*x381 - x376*x379;
    double x397 = -x372*x389 - x373*x390 - x376*x388;
    double x398 = -x361*x376 - x362*x372 - x363*x373 + x391;
    double x399 = x186*x274;
    double x400 = x237*x274;
    double x401 = x107*x274;
    double x402 = x365*x399 + x367*x400 + x369*x401;
    double x403 = x195*x400 + x297*x399 + x299*x401;
    double x404 = x302*x399 + x324*x400 + x326*x401;
    double x405 = x173*x274;
    double x406 = x247*x274;
    double x407 = x109*x274;
    double x408 = x365*x405 + x367*x406 + x369*x407;
    double x409 = x195*x406 + x297*x405 + x299*x407;
    double x410 = x302*x405 + x324*x406 + x326*x407;
    double x411 = x164*x274;
    double x412 = x260*x274;
    double x413 = x104*x274;
    double x414 = x365*x411 + x367*x412 + x369*x413;
    double x415 = x195*x412 + x297*x411 + x299*x413;
    double x416 = x302*x411 + x324*x412 + x326*x413;
    double x417 = x184*x274;
    double x418 = x257*x274;
    double x419 = x274*x77;
    double x420 = x365*x417 + x367*x418 + x369*x419;
    double x421 = x195*x418 + x297*x417 + x299*x419;
    double x422 = x302*x417 + x324*x418 + x326*x419;
    double x423 = x188*x274;
    double x424 = x243*x274;
    double x425 = x274*x91;
    double x426 = x365*x423 + x367*x424 + x369*x425;
    double x427 = x195*x424 + x297*x423 + x299*x425;
    double x428 = x302*x423 + x324*x424 + x326*x425;
    double x429 = x177*x274;
    double x430 = x250*x274;
    double x431 = x274*x96;
    double x432 = x365*x429 + x367*x430 + x369*x431;
    double x433 = x195*x430 + x297*x429 + x299*x431;
    double x434 = x302*x429 + x324*x430 + x326*x431;
    double x435 = x169*x274;
    double x436 = x262*x274;
    double x437 = x274*x84;
    double x438 = x365*x435 + x367*x436 + x369*x437;
    double x439 = x195*x436 + x297*x435 + x299*x437;
    double x440 = x302*x435 + x324*x436 + x326*x437;
    double x441 = x128*x274;
    double x442 = x220*x274;
    double x443 = x23*x274;
    double x444 = x365*x441 + x367*x442 + x369*x443;
    double x445 = x195*x442 + x297*x441 + x299*x443;
    double x446 = x302*x441 + x324*x442 + x326*x443;
    double x447 = x149*x274;
    double x448 = x212*x274;
    double x449 = x274*x7;
    double x450 = x365*x447 + x367*x448 + x369*x449;
    double x451 = x195*x448 + x297*x447 + x299*x449;
    double x452 = x302*x447 + x324*x448 + x326*x449;
    double x453 = x133*x274;
    double x454 = x225*x274;
    double x455 = x17*x274;
    double x456 = x365*x453 + x367*x454 + x369*x455;
    double x457 = x195*x454 + x297*x453 + x299*x455;
    double x458 = x302*x453 + x324*x454 + x326*x455;
    double x459 = x143*x274;
    double x460 = x208*x274;
    double x461 = x20*x274;
    double x462 = x365*x459 + x367*x460 + x369*x461;
    double x463 = x195*x460 + x297*x459 + x299*x461;
    double x464 = x302*x459 + x324*x460 + x326*x461;
    double x465 = x122*x274;
    double x466 = x222*x274;
    double x467 = x274*x29;
    double x468 = x365*x465 + x367*x466 + x369*x467;
    double x469 = x195*x466 + x297*x465 + x299*x467;
    double x470 = x302*x465 + x324*x466 + x326*x467;
    double x471 = x151*x274;
    double x472 = x214*x274;
    double x473 = x274*x33;
    double x474 = x365*x471 + x367*x472 + x369*x473;
    double x475 = x195*x472 + x297*x471 + x299*x473;
    double x476 = x302*x471 + x324*x472 + x326*x473;
    double x477 = x135*x274;
    double x478 = x227*x274;
    double x479 = x274*x31;
    double x480 = x365*x477 + x367*x478 + x369*x479;
    double x481 = x195*x478 + x297*x477 + x299*x479;
    double x482 = x302*x477 + x324*x478 + x326*x479;
    double x483 = x145*x274;
    double x484 = x202*x274;
    double x485 = x26*x274;
    double x486 = x365*x483 + x367*x484 + x369*x485;
    double x487 = x195*x484 + x297*x483 + x299*x485;
    double x488 = x302*x483 + x324*x484 + x326*x485;
    double x489 = x125*x274;
    double x490 = x205*x274;
    double x491 = x274*x43;
    double x492 = x365*x489 + x367*x490 + x369*x491;
    double x493 = x195*x490 + x297*x489 + x299*x491;
    double x494 = x302*x489 + x324*x490 + x326*x491;
    double x495 = x117*x274;
    double x496 = x210*x274;
    double x497 = x274*x52;
    double x498 = x365*x495 + x367*x496 + x369*x497;
    double x499 = x195*x496 + x297*x495 + x299*x497;
    double x500 = x302*x495 + x324*x496 + x326*x497;
    double x501 = x137*x274;
    double x502 = x216*x274;
    double x503 = x274*x55;
    double x504 = x365*x501 + x367*x502 + x369*x503;
    double x505 = x195*x502 + x297*x501 + x299*x503;
    double x506 = x302*x501 + x324*x502 + x326*x503;
    double x507 = x131*x274;
    double x508 = x198*x274;
    double x509 = x274*x47;
    double x510 = x365*x507 + x367*x508 + x369*x509;
    double x511 = x195*x508 + x297*x507 + x299*x509;
    double x512 = x302*x507 + x324*x508 + x326*x509;
    
    res_0[0] = x273*(x364*x371 + x382*x383 + x384*x392);
    res_0[1] = x273*(x371*x393 + x383*x395 + x384*x394);
    res_0[2] = x273*(x371*x398 + x383*x396 + x384*x397);
    res_0[3] = x273*(x364*x402 + x382*x403 + x392*x404);
    res_0[4] = x273*(x393*x402 + x394*x404 + x395*x403);
    res_0[5] = x273*(x396*x403 + x397*x404 + x398*x402);
    res_0[6] = x273*(x364*x408 + x382*x409 + x392*x410);
    res_0[7] = x273*(x393*x408 + x394*x410 + x395*x409);
    res_0[8] = x273*(x396*x409 + x397*x410 + x398*x408);
    res_0[9] = x273*(x364*x414 + x382*x415 + x392*x416);
    res_0[10] = x273*(x393*x414 + x394*x416 + x395*x415);
    res_0[11] = x273*(x396*x415 + x397*x416 + x398*x414);
    res_0[12] = x273*(x364*x420 + x382*x421 + x392*x422);
    res_0[13] = x273*(x393*x420 + x394*x422 + x395*x421);
    res_0[14] = x273*(x396*x421 + x397*x422 + x398*x420);
    res_0[15] = x273*(x364*x426 + x382*x427 + x392*x428);
    res_0[16] = x273*(x393*x426 + x394*x428 + x395*x427);
    res_0[17] = x273*(x396*x427 + x397*x428 + x398*x426);
    res_0[18] = x273*(x364*x432 + x382*x433 + x392*x434);
    res_0[19] = x273*(x393*x432 + x394*x434 + x395*x433);
    res_0[20] = x273*(x396*x433 + x397*x434 + x398*x432);
    res_0[21] = x273*(x364*x438 + x382*x439 + x392*x440);
    res_0[22] = x273*(x393*x438 + x394*x440 + x395*x439);
    res_0[23] = x273*(x396*x439 + x397*x440 + x398*x438);
    res_0[24] = x273*(x364*x444 + x382*x445 + x392*x446);
    res_0[25] = x273*(x393*x444 + x394*x446 + x395*x445);
    res_0[26] = x273*(x396*x445 + x397*x446 + x398*x444);
    res_0[27] = x273*(x364*x450 + x382*x451 + x392*x452);
    res_0[28] = x273*(x393*x450 + x394*x452 + x395*x451);
    res_0[29] = x273*(x396*x451 + x397*x452 + x398*x450);
    res_0[30] = x273*(x364*x456 + x382*x457 + x392*x458);
    res_0[31] = x273*(x393*x456 + x394*x458 + x395*x457);
    res_0[32] = x273*(x396*x457 + x397*x458 + x398*x456);
    res_0[33] = x273*(x364*x462 + x382*x463 + x392*x464);
    res_0[34] = x273*(x393*x462 + x394*x464 + x395*x463);
    res_0[35] = x273*(x396*x463 + x397*x464 + x398*x462);
    res_0[36] = x273*(x364*x468 + x382*x469 + x392*x470);
    res_0[37] = x273*(x393*x468 + x394*x470 + x395*x469);
    res_0[38] = x273*(x396*x469 + x397*x470 + x398*x468);
    res_0[39] = x273*(x364*x474 + x382*x475 + x392*x476);
    res_0[40] = x273*(x393*x474 + x394*x476 + x395*x475);
    res_0[41] = x273*(x396*x475 + x397*x476 + x398*x474);
    res_0[42] = x273*(x364*x480 + x382*x481 + x392*x482);
    res_0[43] = x273*(x393*x480 + x394*x482 + x395*x481);
    res_0[44] = x273*(x396*x481 + x397*x482 + x398*x480);
    res_0[45] = x273*(x364*x486 + x382*x487 + x392*x488);
    res_0[46] = x273*(x393*x486 + x394*x488 + x395*x487);
    res_0[47] = x273*(x396*x487 + x397*x488 + x398*x486);
    res_0[48] = x273*(x364*x492 + x382*x493 + x392*x494);
    res_0[49] = x273*(x393*x492 + x394*x494 + x395*x493);
    res_0[50] = x273*(x396*x493 + x397*x494 + x398*x492);
    res_0[51] = x273*(x364*x498 + x382*x499 + x392*x500);
    res_0[52] = x273*(x393*x498 + x394*x500 + x395*x499);
    res_0[53] = x273*(x396*x499 + x397*x500 + x398*x498);
    res_0[54] = x273*(x364*x504 + x382*x505 + x392*x506);
    res_0[55] = x273*(x393*x504 + x394*x506 + x395*x505);
    res_0[56] = x273*(x396*x505 + x397*x506 + x398*x504);
    res_0[57] = x273*(x364*x510 + x382*x511 + x392*x512);
    res_0[58] = x273*(x393*x510 + x394*x512 + x395*x511);
    res_0[59] = x273*(x396*x511 + x397*x512 + x398*x510);
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
    
    double x0 = s*s;
    double x1 = (1.0/4.0)*x0;
    double x2 = x1 - 1.0/4.0;
    double x3 = (1.0/4.0)*r;
    double x4 = -x3;
    double x5 = x0*x3;
    double x6 = x4 + x5;
    double x7 = x2 + x6;
    double *x8 = coord;
    double x9 = x8[27];
    double x10 = r*r;
    double x11 = (1.0/4.0)*x10;
    double x12 = x11 - 1.0/4.0;
    double x13 = (1.0/4.0)*s;
    double x14 = -x13;
    double x15 = x10*x13;
    double x16 = x14 + x15;
    double x17 = x12 + x16;
    double x18 = x8[30];
    double x19 = x3 - x5;
    double x20 = x19 + x2;
    double x21 = x8[33];
    double x22 = x13 - x15;
    double x23 = x12 + x22;
    double x24 = x8[24];
    double x25 = 1.0/4.0 - x1;
    double x26 = x25 + x6;
    double x27 = x8[45];
    double x28 = 1.0/4.0 - x11;
    double x29 = x16 + x28;
    double x30 = x8[36];
    double x31 = x22 + x28;
    double x32 = x8[42];
    double x33 = x19 + x25;
    double x34 = x8[39];
    double x35 = (1.0/2.0)*t;
    double x36 = -x35;
    double x37 = r*x35;
    double x38 = x36 + x37;
    double x39 = s*x35;
    double x40 = s*x37;
    double x41 = -x40;
    double x42 = x39 + x41;
    double x43 = x38 + x42;
    double x44 = x8[48];
    double x45 = -x39;
    double x46 = x40 + x45;
    double x47 = x38 + x46;
    double x48 = x8[57];
    double x49 = -x37;
    double x50 = x36 + x49;
    double x51 = x39 + x40;
    double x52 = x50 + x51;
    double x53 = x8[51];
    double x54 = x41 + x45;
    double x55 = x50 + x54;
    double x56 = x8[54];
    double x57 = (1.0/4.0)*t;
    double x58 = t*x13;
    double x59 = -x58;
    double x60 = t*x3;
    double x61 = s*x60;
    double x62 = x59 + x61;
    double x63 = x57 + x62;
    double x64 = (1.0/8.0)*x10;
    double x65 = s*x64;
    double x66 = -x65;
    double x67 = (1.0/8.0)*r;
    double x68 = s*x67;
    double x69 = x66 + x68;
    double x70 = (1.0/8.0)*x0;
    double x71 = r*x70;
    double x72 = -x71;
    double x73 = x64 - 1.0/8.0;
    double x74 = -x60;
    double x75 = x70 + x74;
    double x76 = x72 + x73 + x75;
    double x77 = x63 + x69 + x76;
    double x78 = x8[12];
    double x79 = -x68;
    double x80 = x65 + x79;
    double x81 = -x61;
    double x82 = x58 + x81;
    double x83 = x57 + x82;
    double x84 = x76 + x80 + x83;
    double x85 = x8[21];
    double x86 = x60 + x70;
    double x87 = x59 + x81;
    double x88 = x73 + x87;
    double x89 = x66 + x71;
    double x90 = x57 + x79;
    double x91 = x86 + x88 + x89 + x90;
    double x92 = x8[15];
    double x93 = x65 + x71;
    double x94 = x57 + x68;
    double x95 = x58 + x61;
    double x96 = x73 + x95;
    double x97 = x86 + x93 + x94 + x96;
    double x98 = x8[18];
    double x99 = -x70;
    double x100 = 1.0/8.0 - x64;
    double x101 = x100 + x99;
    double x102 = x101 + x71 + x74;
    double x103 = x102 + x63 + x80;
    double x104 = x8[0];
    double x105 = x102 + x69 + x83;
    double x106 = x8[9];
    double x107 = x65 + x72;
    double x108 = x101 + x60;
    double x109 = x107 + x108 + x87 + x94;
    double x110 = x8[3];
    double x111 = x66 + x72;
    double x112 = x108 + x111 + x90 + x95;
    double x113 = x8[6];
    double x114 = x103*x104 + x105*x106 + x109*x110 + x112*x113 + x17*x18 + x20*x21 + x23*x24 + x26*x27 + x29*x30 + x31*x32 + x33*x34 + x43*x44 + x47*x48 + x52*x53 + x55*x56 + x7*x9 + x77*x78 + x84*x85 + x91*x92 + x97*x98;
    double x115 = t*t;
    double x116 = (1.0/4.0)*x115;
    double x117 = x116 - 1.0/4.0;
    double x118 = x115*x3;
    double x119 = x118 + x4;
    double x120 = x117 + x119;
    double x121 = x8[53];
    double x122 = -x57;
    double x123 = t*x11;
    double x124 = x122 + x123;
    double x125 = x12 + x124;
    double x126 = x8[38];
    double x127 = -x118 + x3;
    double x128 = x117 + x127;
    double x129 = x8[50];
    double x130 = -x123 + x57;
    double x131 = x12 + x130;
    double x132 = x8[26];
    double x133 = 1.0/4.0 - x116;
    double x134 = x119 + x133;
    double x135 = x8[59];
    double x136 = x124 + x28;
    double x137 = x8[32];
    double x138 = x130 + x28;
    double x139 = x8[44];
    double x140 = x127 + x133;
    double x141 = x8[56];
    double x142 = (1.0/2.0)*s;
    double x143 = -x142;
    double x144 = r*x142;
    double x145 = x143 + x144;
    double x146 = x145 + x42;
    double x147 = x8[35];
    double x148 = x145 + x46;
    double x149 = x8[47];
    double x150 = -x144;
    double x151 = x143 + x150;
    double x152 = x151 + x51;
    double x153 = x8[29];
    double x154 = x151 + x54;
    double x155 = x8[41];
    double x156 = t*x64;
    double x157 = -x156;
    double x158 = x13 + x157;
    double x159 = (1.0/8.0)*x115;
    double x160 = r*x159;
    double x161 = -x160;
    double x162 = t*x67;
    double x163 = x161 + x162;
    double x164 = s*x3;
    double x165 = -x164;
    double x166 = x159 + x165;
    double x167 = x166 + x73;
    double x168 = x158 + x163 + x167 + x62;
    double x169 = x8[11];
    double x170 = x13 + x156;
    double x171 = -x162;
    double x172 = x161 + x171;
    double x173 = x167 + x170 + x172 + x82;
    double x174 = x8[23];
    double x175 = x160 + x171;
    double x176 = x158 + x164;
    double x177 = x159 + x175 + x176 + x88;
    double x178 = x8[8];
    double x179 = x160 + x162;
    double x180 = x164 + x170;
    double x181 = x159 + x179 + x180 + x96;
    double x182 = x8[20];
    double x183 = -x159;
    double x184 = x100 + x183;
    double x185 = x165 + x184;
    double x186 = x170 + x175 + x185 + x62;
    double x187 = x8[2];
    double x188 = x158 + x179 + x185 + x82;
    double x189 = x8[14];
    double x190 = x163 + x180 + x184 + x87;
    double x191 = x8[5];
    double x192 = x172 + x176 + x184 + x95;
    double x193 = x8[17];
    double x194 = x120*x121 + x125*x126 + x128*x129 + x131*x132 + x134*x135 + x136*x137 + x138*x139 + x140*x141 + x146*x147 + x148*x149 + x152*x153 + x154*x155 + x168*x169 + x173*x174 + x177*x178 + x181*x182 + x186*x187 + x188*x189 + x190*x191 + x192*x193;
    double x195 = x114*x194;
    double x196 = x103*x187 + x105*x169 + x109*x191 + x112*x178 + x121*x52 + x126*x29 + x129*x43 + x132*x23 + x135*x47 + x137*x17 + x139*x31 + x141*x55 + x147*x20 + x149*x26 + x153*x7 + x155*x33 + x174*x84 + x182*x97 + x189*x77 + x193*x91;
    double x197 = x104*x186 + x106*x168 + x110*x190 + x113*x177 + x120*x53 + x125*x30 + x128*x44 + x131*x24 + x134*x48 + x136*x18 + x138*x32 + x140*x56 + x146*x21 + x148*x27 + x152*x9 + x154*x34 + x173*x85 + x181*x98 + x188*x78 + x192*x92;
    double x198 = x196*x197;
    double x199 = x195 - x198;
    double x200 = x115*x13;
    double x201 = x14 + x200;
    double x202 = x117 + x201;
    double x203 = x8[58];
    double x204 = t*x1;
    double x205 = x122 + x204;
    double x206 = x2 + x205;
    double x207 = x8[46];
    double x208 = x13 - x200;
    double x209 = x117 + x208;
    double x210 = x8[49];
    double x211 = -x204 + x57;
    double x212 = x2 + x211;
    double x213 = x8[34];
    double x214 = x133 + x201;
    double x215 = x8[52];
    double x216 = x205 + x25;
    double x217 = x8[28];
    double x218 = x211 + x25;
    double x219 = x8[40];
    double x220 = x133 + x208;
    double x221 = x8[55];
    double x222 = -1.0/2.0*r;
    double x223 = x144 + x222;
    double x224 = x223 + x37 + x41;
    double x225 = x8[25];
    double x226 = x223 + x40 + x49;
    double x227 = x8[37];
    double x228 = x150 + x222;
    double x229 = x228 + x37 + x40;
    double x230 = x8[31];
    double x231 = x228 + x41 + x49;
    double x232 = x8[43];
    double x233 = s*x159;
    double x234 = -x233;
    double x235 = t*x70;
    double x236 = -x235;
    double x237 = x234 + x236;
    double x238 = x166 - 1.0/8.0;
    double x239 = (1.0/8.0)*s;
    double x240 = t*x239;
    double x241 = x3 + x61;
    double x242 = x240 + x241;
    double x243 = x237 + x238 + x242 + x75;
    double x244 = x8[4];
    double x245 = x235 + x86;
    double x246 = -x240;
    double x247 = x3 + x81;
    double x248 = x246 + x247;
    double x249 = x234 + x238 + x245 + x248;
    double x250 = x8[16];
    double x251 = x233 + x236;
    double x252 = x159 + x251;
    double x253 = x164 - 1.0/8.0;
    double x254 = x248 + x252 + x253 + x75;
    double x255 = x8[7];
    double x256 = x159 + x233;
    double x257 = x242 + x245 + x253 + x256;
    double x258 = x8[19];
    double x259 = x165 + x183 + x99 + 1.0/8.0;
    double x260 = x235 + x74;
    double x261 = x241 + x246;
    double x262 = x233 + x259 + x260 + x261;
    double x263 = x8[1];
    double x264 = x240 + x247;
    double x265 = x251 + x259 + x264 + x60;
    double x266 = x8[13];
    double x267 = x164 + x183 + x99 + 1.0/8.0;
    double x268 = x234 + x260 + x264 + x267;
    double x269 = x8[10];
    double x270 = x237 + x261 + x267 + x60;
    double x271 = x8[22];
    double x272 = x202*x203 + x206*x207 + x209*x210 + x212*x213 + x214*x215 + x216*x217 + x218*x219 + x220*x221 + x224*x225 + x226*x227 + x229*x230 + x231*x232 + x243*x244 + x249*x250 + x254*x255 + x257*x258 + x262*x263 + x265*x266 + x268*x269 + x270*x271;
    double x273 = x103*x263 + x105*x269 + x109*x244 + x112*x255 + x17*x230 + x20*x213 + x203*x47 + x207*x26 + x210*x43 + x215*x52 + x217*x7 + x219*x33 + x221*x55 + x225*x23 + x227*x29 + x232*x31 + x250*x91 + x258*x97 + x266*x77 + x271*x84;
    double x274 = x121*x214 + x126*x226 + x129*x209 + x132*x224 + x135*x202 + x137*x229 + x139*x231 + x141*x220 + x147*x212 + x149*x206 + x153*x216 + x155*x218 + x169*x268 + x174*x270 + x178*x254 + x182*x257 + x187*x262 + x189*x265 + x191*x243 + x193*x249;
    double x275 = x197*x274;
    double x276 = x120*x215 + x125*x227 + x128*x210 + x131*x225 + x134*x203 + x136*x230 + x138*x232 + x140*x221 + x146*x213 + x148*x207 + x152*x217 + x154*x219 + x168*x269 + x173*x271 + x177*x255 + x181*x258 + x186*x263 + x188*x266 + x190*x244 + x192*x250;
    double x277 = x104*x262 + x106*x268 + x110*x243 + x113*x254 + x18*x229 + x202*x48 + x206*x27 + x209*x44 + x21*x212 + x214*x53 + x216*x9 + x218*x34 + x220*x56 + x224*x24 + x226*x30 + x231*x32 + x249*x92 + x257*x98 + x265*x78 + x270*x85;
    double x278 = x196*x277;
    double x279 = x114*x274;
    double x280 = x194*x277;
    double x281 = x195*x272 - x198*x272 + x273*x275 - x273*x280 + x276*x278 - x276*x279;
    double x282 = 1.0/x281;
    double *x283 = U;
    double x284 = x283[59];
    double x285 = x283[47];
    double x286 = x283[50];
    double x287 = x283[35];
    double x288 = x283[53];
    double x289 = x283[29];
    double x290 = x283[41];
    double x291 = x283[56];
    double x292 = x283[26];
    double x293 = x283[38];
    double x294 = x283[32];
    double x295 = x283[44];
    double x296 = x283[5];
    double x297 = x283[17];
    double x298 = x283[8];
    double x299 = x283[20];
    double x300 = x283[2];
    double x301 = x283[14];
    double x302 = x283[11];
    double x303 = x283[23];
    double x304 = x282*(x202*x284 + x206*x285 + x209*x286 + x212*x287 + x214*x288 + x216*x289 + x218*x290 + x220*x291 + x224*x292 + x226*x293 + x229*x294 + x231*x295 + x243*x296 + x249*x297 + x254*x298 + x257*x299 + x262*x300 + x265*x301 + x268*x302 + x270*x303);
    double x305 = x278 - x279;
    double x306 = x282*(x120*x288 + x125*x293 + x128*x286 + x131*x292 + x134*x284 + x136*x294 + x138*x295 + x140*x291 + x146*x287 + x148*x285 + x152*x289 + x154*x290 + x168*x302 + x173*x303 + x177*x298 + x181*x299 + x186*x300 + x188*x301 + x190*x296 + x192*x297);
    double x307 = x275 - x280;
    double x308 = x282*(x103*x300 + x105*x302 + x109*x296 + x112*x298 + x17*x294 + x20*x287 + x23*x292 + x26*x285 + x284*x47 + x286*x43 + x288*x52 + x289*x7 + x29*x293 + x290*x33 + x291*x55 + x295*x31 + x297*x91 + x299*x97 + x301*x77 + x303*x84);
    double x309 = x199*x304 + x305*x306 + x307*x308;
    double x310 = -x196*x272 + x273*x274;
    double x311 = x283[52];
    double x312 = x283[37];
    double x313 = x283[49];
    double x314 = x283[25];
    double x315 = x283[58];
    double x316 = x283[31];
    double x317 = x283[43];
    double x318 = x283[55];
    double x319 = x283[34];
    double x320 = x283[46];
    double x321 = x283[28];
    double x322 = x283[40];
    double x323 = x283[10];
    double x324 = x283[22];
    double x325 = x283[7];
    double x326 = x283[19];
    double x327 = x283[1];
    double x328 = x283[13];
    double x329 = x283[4];
    double x330 = x283[16];
    double x331 = x282*(x120*x311 + x125*x312 + x128*x313 + x131*x314 + x134*x315 + x136*x316 + x138*x317 + x140*x318 + x146*x319 + x148*x320 + x152*x321 + x154*x322 + x168*x323 + x173*x324 + x177*x325 + x181*x326 + x186*x327 + x188*x328 + x190*x329 + x192*x330);
    double x332 = -x194*x273 + x196*x276;
    double x333 = x282*(x202*x315 + x206*x320 + x209*x313 + x212*x319 + x214*x311 + x216*x321 + x218*x322 + x220*x318 + x224*x314 + x226*x312 + x229*x316 + x231*x317 + x243*x329 + x249*x330 + x254*x325 + x257*x326 + x262*x327 + x265*x328 + x268*x323 + x270*x324);
    double x334 = x194*x272 - x274*x276;
    double x335 = x282*(x103*x327 + x105*x323 + x109*x329 + x112*x325 + x17*x316 + x20*x319 + x23*x314 + x26*x320 + x29*x312 + x31*x317 + x311*x52 + x313*x43 + x315*x47 + x318*x55 + x321*x7 + x322*x33 + x324*x84 + x326*x97 + x328*x77 + x330*x91);
    double x336 = x310*x331 + x332*x333 + x334*x335;
    double x337 = x304*x332 + x306*x310 + x308*x334;
    double x338 = x199*x333 + x305*x331 + x307*x335;
    double x339 = x338 + 1.0;
    double x340 = x309*x336 - x337*x339;
    double x341 = x283[57];
    double x342 = x283[45];
    double x343 = x283[48];
    double x344 = x283[33];
    double x345 = x283[51];
    double x346 = x283[27];
    double x347 = x283[39];
    double x348 = x283[54];
    double x349 = x283[24];
    double x350 = x283[36];
    double x351 = x283[30];
    double x352 = x283[42];
    double x353 = x283[3];
    double x354 = x283[15];
    double x355 = x283[6];
    double x356 = x283[18];
    double x357 = x283[0];
    double x358 = x283[12];
    double x359 = x283[9];
    double x360 = x283[21];
    double x361 = x282*(x202*x341 + x206*x342 + x209*x343 + x212*x344 + x214*x345 + x216*x346 + x218*x347 + x220*x348 + x224*x349 + x226*x350 + x229*x351 + x231*x352 + x243*x353 + x249*x354 + x254*x355 + x257*x356 + x262*x357 + x265*x358 + x268*x359 + x270*x360);
    double x362 = x282*(x120*x345 + x125*x350 + x128*x343 + x131*x349 + x134*x341 + x136*x351 + x138*x352 + x140*x348 + x146*x344 + x148*x342 + x152*x346 + x154*x347 + x168*x359 + x173*x360 + x177*x355 + x181*x356 + x186*x357 + x188*x358 + x190*x353 + x192*x354);
    double x363 = x282*(x103*x357 + x105*x359 + x109*x353 + x112*x355 + x17*x351 + x20*x344 + x23*x349 + x26*x342 + x29*x350 + x31*x352 + x33*x347 + x341*x47 + x343*x43 + x345*x52 + x346*x7 + x348*x55 + x354*x91 + x356*x97 + x358*x77 + x360*x84);
    double x364 = x199*x361 + x305*x362 + x307*x363;
    double x365 = x310*x362 + x332*x361 + x334*x363;
    double x366 = x365 + 1.0;
    double x367 = -x309*x366 + x337*x364;
    double x368 = -x336*x364 + x339*x366;
    double x369 = S11*x340 + S12*x367 + S13*x368;
    double x370 = S12*x340 + S22*x367 + S23*x368;
    double x371 = S13*x340 + S23*x367 + S33*x368;
    double x372 = -x336*x370 - x337*x371 - x365*x369;
    double x373 = x114*x272 - x273*x277;
    double x374 = x186*x282;
    double x375 = -x114*x276 + x197*x273;
    double x376 = x262*x282;
    double x377 = -x197*x272 + x276*x277;
    double x378 = x103*x282;
    double x379 = x373*x374 + x375*x376 + x377*x378;
    double x380 = x331*x373 + x333*x375 + x335*x377;
    double x381 = x304*x375 + x306*x373 + x308*x377;
    double x382 = x381 + 1.0;
    double x383 = -x336*x382 + x337*x380;
    double x384 = x361*x375 + x362*x373 + x363*x377;
    double x385 = -x337*x384 + x366*x382;
    double x386 = x336*x384 - x366*x380;
    double x387 = S11*x383 + S12*x385 + S13*x386;
    double x388 = S12*x383 + S22*x385 + S23*x386;
    double x389 = S13*x383 + S23*x385 + S33*x386;
    double x390 = -x336*x388 - x337*x389 - x365*x387;
    double x391 = x199*x376 + x305*x374 + x307*x378;
    double x392 = x310*x374 + x332*x376 + x334*x378;
    double x393 = -x309*x380 + x339*x382;
    double x394 = x309*x384 - x364*x382;
    double x395 = -x339*x384 + x364*x380;
    double x396 = S11*x393 + S12*x394 + S13*x395;
    double x397 = S12*x393 + S22*x394 + S23*x395;
    double x398 = S13*x393 + S23*x394 + S33*x395;
    double x399 = s*x118;
    double x400 = x164 - x399;
    double x401 = x134 + x201 + x400;
    double *x402 = V;
    double x403 = x402[48];
    double x404 = t*x5;
    double x405 = -x404 + x60;
    double x406 = x205 + x26 + x405;
    double x407 = x402[33];
    double x408 = x404 + x74;
    double x409 = x211 + x26 + x408;
    double x410 = x402[45];
    double x411 = x165 + x399;
    double x412 = x134 + x208 + x411;
    double x413 = x402[57];
    double x414 = t*x15;
    double x415 = -x414 + x58;
    double x416 = x124 + x29 + x415;
    double x417 = x402[24];
    double x418 = x414 + x59;
    double x419 = x130 + x29 + x418;
    double x420 = x402[36];
    double x421 = x140 + x201 + x411;
    double x422 = x402[51];
    double x423 = x124 + x31 + x418;
    double x424 = x402[30];
    double x425 = x205 + x33 + x408;
    double x426 = x402[27];
    double x427 = x130 + x31 + x415;
    double x428 = x402[42];
    double x429 = x211 + x33 + x405;
    double x430 = x402[39];
    double x431 = x140 + x208 + x400;
    double x432 = x402[54];
    double x433 = -x239;
    double x434 = t*x68;
    double x435 = t*x71;
    double x436 = x433 + x434 + x435 + x64 + x70 - 1.0/4.0;
    double x437 = s*x160;
    double x438 = t*x65;
    double x439 = x437 + x438;
    double x440 = x160 - x67;
    double x441 = x440 + x93;
    double x442 = (1.0/8.0)*t;
    double x443 = x156 + x235 - x442;
    double x444 = x256 + x443;
    double x445 = x436 + x439 + x441 + x444;
    double x446 = x402[18];
    double x447 = -x438;
    double x448 = x157 + x252 + x442 + x447;
    double x449 = -x435;
    double x450 = -x434 + x64 + x70 - 1.0/4.0;
    double x451 = x433 + x449 + x450;
    double x452 = x437 + x441 + x448 + x451;
    double x453 = x402[6];
    double x454 = x159 + x239 + x435 + x450;
    double x455 = -x437;
    double x456 = x440 + x455 + x89;
    double x457 = x234 + x443 + x447;
    double x458 = x454 + x456 + x457;
    double x459 = x402[15];
    double x460 = x157 + x237 + x442;
    double x461 = x159 + x239 + x434 + x449 + x64 + x70 - 1.0/4.0;
    double x462 = x438 + x456 + x460 + x461;
    double x463 = x402[3];
    double x464 = x161 + x67;
    double x465 = x107 + x455 + x464;
    double x466 = x438 + x444 + x451 + x465;
    double x467 = x402[21];
    double x468 = x436 + x448 + x465;
    double x469 = x402[9];
    double x470 = x111 + x464;
    double x471 = x437 + x457 + x461 + x470;
    double x472 = x402[12];
    double x473 = x439 + x454 + x460 + x470;
    double x474 = x402[0];
    double x475 = x401*x403 + x406*x407 + x409*x410 + x412*x413 + x416*x417 + x419*x420 + x421*x422 + x423*x424 + x425*x426 + x427*x428 + x429*x430 + x431*x432 + x445*x446 + x452*x453 + x458*x459 + x462*x463 + x466*x467 + x468*x469 + x471*x472 + x473*x474;
    double x476 = x402[49];
    double x477 = x402[34];
    double x478 = x402[46];
    double x479 = x402[58];
    double x480 = x402[25];
    double x481 = x402[37];
    double x482 = x402[52];
    double x483 = x402[31];
    double x484 = x402[28];
    double x485 = x402[43];
    double x486 = x402[40];
    double x487 = x402[55];
    double x488 = x402[19];
    double x489 = x402[7];
    double x490 = x402[16];
    double x491 = x402[4];
    double x492 = x402[22];
    double x493 = x402[10];
    double x494 = x402[13];
    double x495 = x402[1];
    double x496 = x401*x476 + x406*x477 + x409*x478 + x412*x479 + x416*x480 + x419*x481 + x421*x482 + x423*x483 + x425*x484 + x427*x485 + x429*x486 + x431*x487 + x445*x488 + x452*x489 + x458*x490 + x462*x491 + x466*x492 + x468*x493 + x471*x494 + x473*x495;
    double x497 = x402[50];
    double x498 = x402[35];
    double x499 = x402[47];
    double x500 = x402[59];
    double x501 = x402[26];
    double x502 = x402[38];
    double x503 = x402[53];
    double x504 = x402[32];
    double x505 = x402[29];
    double x506 = x402[44];
    double x507 = x402[41];
    double x508 = x402[56];
    double x509 = x402[20];
    double x510 = x402[8];
    double x511 = x402[17];
    double x512 = x402[5];
    double x513 = x402[23];
    double x514 = x402[11];
    double x515 = x402[14];
    double x516 = x402[2];
    double x517 = x401*x497 + x406*x498 + x409*x499 + x412*x500 + x416*x501 + x419*x502 + x421*x503 + x423*x504 + x425*x505 + x427*x506 + x429*x507 + x431*x508 + x445*x509 + x452*x510 + x458*x511 + x462*x512 + x466*x513 + x468*x514 + x471*x515 + x473*x516;
    double x518 = 1.0*PENER + 1.0*SENER - 0.5*rho*(x475*x475 + x496*x496 + x517*x517);
    double x519 = -x336*x397 - x337*x398 - x365*x396 + x518;
    double *x520 = A;
    double x521 = x401*x520[48] + x406*x520[33] + x409*x520[45] + x412*x520[57] + x416*x520[24] + x419*x520[36] + x421*x520[51] + x423*x520[30] + x425*x520[27] + x427*x520[42] + x429*x520[39] + x431*x520[54] + x445*x520[18] + x452*x520[6] + x458*x520[15] + x462*x520[3] + x466*x520[21] + x468*x520[9] + x471*x520[12] + x473*x520[0];
    double x522 = x401*x520[49] + x406*x520[34] + x409*x520[46] + x412*x520[58] + x416*x520[25] + x419*x520[37] + x421*x520[52] + x423*x520[31] + x425*x520[28] + x427*x520[43] + x429*x520[40] + x431*x520[55] + x445*x520[19] + x452*x520[7] + x458*x520[16] + x462*x520[4] + x466*x520[22] + x468*x520[10] + x471*x520[13] + x473*x520[1];
    double x523 = x401*x520[50] + x406*x520[35] + x409*x520[47] + x412*x520[59] + x416*x520[26] + x419*x520[38] + x421*x520[53] + x423*x520[32] + x425*x520[29] + x427*x520[44] + x429*x520[41] + x431*x520[56] + x445*x520[20] + x452*x520[8] + x458*x520[17] + x462*x520[5] + x466*x520[23] + x468*x520[11] + x471*x520[14] + x473*x520[2];
    double x524 = x120*x422 + x125*x420 + x128*x403 + x131*x417 + x134*x413 + x136*x424 + x138*x428 + x140*x432 + x146*x407 + x148*x410 + x152*x426 + x154*x430 + x168*x469 + x173*x467 + x177*x453 + x181*x446 + x186*x474 + x188*x472 + x190*x463 + x192*x459;
    double x525 = x282*x310;
    double x526 = x202*x413 + x206*x410 + x209*x403 + x212*x407 + x214*x422 + x216*x426 + x218*x430 + x220*x432 + x224*x417 + x226*x420 + x229*x424 + x231*x428 + x243*x463 + x249*x459 + x254*x453 + x257*x446 + x262*x474 + x265*x472 + x268*x469 + x270*x467;
    double x527 = x282*x332;
    double x528 = x103*x474 + x105*x469 + x109*x463 + x112*x453 + x17*x424 + x20*x407 + x23*x417 + x26*x410 + x29*x420 + x31*x428 + x33*x430 + x403*x43 + x413*x47 + x422*x52 + x426*x7 + x432*x55 + x446*x97 + x459*x91 + x467*x84 + x472*x77;
    double x529 = x282*x334;
    double x530 = x120*x482 + x125*x481 + x128*x476 + x131*x480 + x134*x479 + x136*x483 + x138*x485 + x140*x487 + x146*x477 + x148*x478 + x152*x484 + x154*x486 + x168*x493 + x173*x492 + x177*x489 + x181*x488 + x186*x495 + x188*x494 + x190*x491 + x192*x490;
    double x531 = x202*x479 + x206*x478 + x209*x476 + x212*x477 + x214*x482 + x216*x484 + x218*x486 + x220*x487 + x224*x480 + x226*x481 + x229*x483 + x231*x485 + x243*x491 + x249*x490 + x254*x489 + x257*x488 + x262*x495 + x265*x494 + x268*x493 + x270*x492;
    double x532 = x103*x495 + x105*x493 + x109*x491 + x112*x489 + x17*x483 + x20*x477 + x23*x480 + x26*x478 + x29*x481 + x31*x485 + x33*x486 + x43*x476 + x47*x479 + x482*x52 + x484*x7 + x487*x55 + x488*x97 + x490*x91 + x492*x84 + x494*x77;
    double x533 = x120*x503 + x125*x502 + x128*x497 + x131*x501 + x134*x500 + x136*x504 + x138*x506 + x140*x508 + x146*x498 + x148*x499 + x152*x505 + x154*x507 + x168*x514 + x173*x513 + x177*x510 + x181*x509 + x186*x516 + x188*x515 + x190*x512 + x192*x511;
    double x534 = x202*x500 + x206*x499 + x209*x497 + x212*x498 + x214*x503 + x216*x505 + x218*x507 + x220*x508 + x224*x501 + x226*x502 + x229*x504 + x231*x506 + x243*x512 + x249*x511 + x254*x510 + x257*x509 + x262*x516 + x265*x515 + x268*x514 + x270*x513;
    double x535 = x103*x516 + x105*x514 + x109*x512 + x112*x510 + x17*x504 + x20*x498 + x23*x501 + x26*x499 + x29*x502 + x31*x506 + x33*x507 + x43*x497 + x47*x500 + x503*x52 + x505*x7 + x508*x55 + x509*x97 + x511*x91 + x513*x84 + x515*x77;
    double x536 = x336*x522 + x337*x523 + x365*x521 + x475*(x524*x525 + x526*x527 + x528*x529) + x496*(x525*x530 + x527*x531 + x529*x532) + x517*(x525*x533 + x527*x534 + x529*x535);
    double x537 = rho*x473;
    double x538 = -x309*x371 - x338*x370 - x364*x369;
    double x539 = -x309*x398 - x338*x397 - x364*x396;
    double x540 = -x309*x389 - x338*x388 - x364*x387 + x518;
    double x541 = x199*x282;
    double x542 = x282*x305;
    double x543 = x282*x307;
    double x544 = x309*x523 + x338*x522 + x364*x521 + x475*(x524*x542 + x526*x541 + x528*x543) + x496*(x530*x542 + x531*x541 + x532*x543) + x517*(x533*x542 + x534*x541 + x535*x543);
    double x545 = -x380*x388 - x381*x389 - x384*x387;
    double x546 = -x380*x397 - x381*x398 - x384*x396;
    double x547 = -x369*x384 - x370*x380 - x371*x381 + x518;
    double x548 = x282*x373;
    double x549 = x282*x375;
    double x550 = x282*x377;
    double x551 = x380*x522 + x381*x523 + x384*x521 + x475*(x524*x548 + x526*x549 + x528*x550) + x496*(x530*x548 + x531*x549 + x532*x550) + x517*(x533*x548 + x534*x549 + x535*x550);
    double x552 = x109*x550 + x190*x548 + x243*x549;
    double x553 = x109*x543 + x190*x542 + x243*x541;
    double x554 = x109*x529 + x190*x525 + x243*x527;
    double x555 = rho*x462;
    double x556 = x112*x550 + x177*x548 + x254*x549;
    double x557 = x112*x543 + x177*x542 + x254*x541;
    double x558 = x112*x529 + x177*x525 + x254*x527;
    double x559 = rho*x452;
    double x560 = x105*x550 + x168*x548 + x268*x549;
    double x561 = x105*x543 + x168*x542 + x268*x541;
    double x562 = x105*x529 + x168*x525 + x268*x527;
    double x563 = rho*x468;
    double x564 = x188*x548 + x265*x549 + x550*x77;
    double x565 = x188*x542 + x265*x541 + x543*x77;
    double x566 = x188*x525 + x265*x527 + x529*x77;
    double x567 = rho*x471;
    double x568 = x192*x548 + x249*x549 + x550*x91;
    double x569 = x192*x542 + x249*x541 + x543*x91;
    double x570 = x192*x525 + x249*x527 + x529*x91;
    double x571 = rho*x458;
    double x572 = x181*x548 + x257*x549 + x550*x97;
    double x573 = x181*x542 + x257*x541 + x543*x97;
    double x574 = x181*x525 + x257*x527 + x529*x97;
    double x575 = rho*x445;
    double x576 = x173*x548 + x270*x549 + x550*x84;
    double x577 = x173*x542 + x270*x541 + x543*x84;
    double x578 = x173*x525 + x270*x527 + x529*x84;
    double x579 = rho*x466;
    double x580 = x131*x548 + x224*x549 + x23*x550;
    double x581 = x131*x542 + x224*x541 + x23*x543;
    double x582 = x131*x525 + x224*x527 + x23*x529;
    double x583 = rho*x416;
    double x584 = x152*x548 + x216*x549 + x550*x7;
    double x585 = x152*x542 + x216*x541 + x543*x7;
    double x586 = x152*x525 + x216*x527 + x529*x7;
    double x587 = rho*x425;
    double x588 = x136*x548 + x17*x550 + x229*x549;
    double x589 = x136*x542 + x17*x543 + x229*x541;
    double x590 = x136*x525 + x17*x529 + x229*x527;
    double x591 = rho*x423;
    double x592 = x146*x548 + x20*x550 + x212*x549;
    double x593 = x146*x542 + x20*x543 + x212*x541;
    double x594 = x146*x525 + x20*x529 + x212*x527;
    double x595 = rho*x406;
    double x596 = x125*x548 + x226*x549 + x29*x550;
    double x597 = x125*x542 + x226*x541 + x29*x543;
    double x598 = x125*x525 + x226*x527 + x29*x529;
    double x599 = rho*x419;
    double x600 = x154*x548 + x218*x549 + x33*x550;
    double x601 = x154*x542 + x218*x541 + x33*x543;
    double x602 = x154*x525 + x218*x527 + x33*x529;
    double x603 = rho*x429;
    double x604 = x138*x548 + x231*x549 + x31*x550;
    double x605 = x138*x542 + x231*x541 + x31*x543;
    double x606 = x138*x525 + x231*x527 + x31*x529;
    double x607 = rho*x427;
    double x608 = x148*x548 + x206*x549 + x26*x550;
    double x609 = x148*x542 + x206*x541 + x26*x543;
    double x610 = x148*x525 + x206*x527 + x26*x529;
    double x611 = rho*x409;
    double x612 = x128*x548 + x209*x549 + x43*x550;
    double x613 = x128*x542 + x209*x541 + x43*x543;
    double x614 = x128*x525 + x209*x527 + x43*x529;
    double x615 = rho*x401;
    double x616 = x120*x548 + x214*x549 + x52*x550;
    double x617 = x120*x542 + x214*x541 + x52*x543;
    double x618 = x120*x525 + x214*x527 + x52*x529;
    double x619 = rho*x421;
    double x620 = x140*x548 + x220*x549 + x55*x550;
    double x621 = x140*x542 + x220*x541 + x543*x55;
    double x622 = x140*x525 + x220*x527 + x529*x55;
    double x623 = rho*x431;
    double x624 = x134*x548 + x202*x549 + x47*x550;
    double x625 = x134*x542 + x202*x541 + x47*x543;
    double x626 = x134*x525 + x202*x527 + x47*x529;
    double x627 = rho*x412;
    
    res_0[0] = x281*(x372*x379 + x390*x391 + x392*x519 - x536*x537);
    res_0[1] = x281*(x379*x538 + x391*x540 + x392*x539 - x537*x544);
    res_0[2] = x281*(x379*x547 + x391*x545 + x392*x546 - x537*x551);
    res_0[3] = x281*(x372*x552 + x390*x553 + x519*x554 - x536*x555);
    res_0[4] = x281*(x538*x552 + x539*x554 + x540*x553 - x544*x555);
    res_0[5] = x281*(x545*x553 + x546*x554 + x547*x552 - x551*x555);
    res_0[6] = x281*(x372*x556 + x390*x557 + x519*x558 - x536*x559);
    res_0[7] = x281*(x538*x556 + x539*x558 + x540*x557 - x544*x559);
    res_0[8] = x281*(x545*x557 + x546*x558 + x547*x556 - x551*x559);
    res_0[9] = x281*(x372*x560 + x390*x561 + x519*x562 - x536*x563);
    res_0[10] = x281*(x538*x560 + x539*x562 + x540*x561 - x544*x563);
    res_0[11] = x281*(x545*x561 + x546*x562 + x547*x560 - x551*x563);
    res_0[12] = x281*(x372*x564 + x390*x565 + x519*x566 - x536*x567);
    res_0[13] = x281*(x538*x564 + x539*x566 + x540*x565 - x544*x567);
    res_0[14] = x281*(x545*x565 + x546*x566 + x547*x564 - x551*x567);
    res_0[15] = x281*(x372*x568 + x390*x569 + x519*x570 - x536*x571);
    res_0[16] = x281*(x538*x568 + x539*x570 + x540*x569 - x544*x571);
    res_0[17] = x281*(x545*x569 + x546*x570 + x547*x568 - x551*x571);
    res_0[18] = x281*(x372*x572 + x390*x573 + x519*x574 - x536*x575);
    res_0[19] = x281*(x538*x572 + x539*x574 + x540*x573 - x544*x575);
    res_0[20] = x281*(x545*x573 + x546*x574 + x547*x572 - x551*x575);
    res_0[21] = x281*(x372*x576 + x390*x577 + x519*x578 - x536*x579);
    res_0[22] = x281*(x538*x576 + x539*x578 + x540*x577 - x544*x579);
    res_0[23] = x281*(x545*x577 + x546*x578 + x547*x576 - x551*x579);
    res_0[24] = x281*(x372*x580 + x390*x581 + x519*x582 - x536*x583);
    res_0[25] = x281*(x538*x580 + x539*x582 + x540*x581 - x544*x583);
    res_0[26] = x281*(x545*x581 + x546*x582 + x547*x580 - x551*x583);
    res_0[27] = x281*(x372*x584 + x390*x585 + x519*x586 - x536*x587);
    res_0[28] = x281*(x538*x584 + x539*x586 + x540*x585 - x544*x587);
    res_0[29] = x281*(x545*x585 + x546*x586 + x547*x584 - x551*x587);
    res_0[30] = x281*(x372*x588 + x390*x589 + x519*x590 - x536*x591);
    res_0[31] = x281*(x538*x588 + x539*x590 + x540*x589 - x544*x591);
    res_0[32] = x281*(x545*x589 + x546*x590 + x547*x588 - x551*x591);
    res_0[33] = x281*(x372*x592 + x390*x593 + x519*x594 - x536*x595);
    res_0[34] = x281*(x538*x592 + x539*x594 + x540*x593 - x544*x595);
    res_0[35] = x281*(x545*x593 + x546*x594 + x547*x592 - x551*x595);
    res_0[36] = x281*(x372*x596 + x390*x597 + x519*x598 - x536*x599);
    res_0[37] = x281*(x538*x596 + x539*x598 + x540*x597 - x544*x599);
    res_0[38] = x281*(x545*x597 + x546*x598 + x547*x596 - x551*x599);
    res_0[39] = x281*(x372*x600 + x390*x601 + x519*x602 - x536*x603);
    res_0[40] = x281*(x538*x600 + x539*x602 + x540*x601 - x544*x603);
    res_0[41] = x281*(x545*x601 + x546*x602 + x547*x600 - x551*x603);
    res_0[42] = x281*(x372*x604 + x390*x605 + x519*x606 - x536*x607);
    res_0[43] = x281*(x538*x604 + x539*x606 + x540*x605 - x544*x607);
    res_0[44] = x281*(x545*x605 + x546*x606 + x547*x604 - x551*x607);
    res_0[45] = x281*(x372*x608 + x390*x609 + x519*x610 - x536*x611);
    res_0[46] = x281*(x538*x608 + x539*x610 + x540*x609 - x544*x611);
    res_0[47] = x281*(x545*x609 + x546*x610 + x547*x608 - x551*x611);
    res_0[48] = x281*(x372*x612 + x390*x613 + x519*x614 - x536*x615);
    res_0[49] = x281*(x538*x612 + x539*x614 + x540*x613 - x544*x615);
    res_0[50] = x281*(x545*x613 + x546*x614 + x547*x612 - x551*x615);
    res_0[51] = x281*(x372*x616 + x390*x617 + x519*x618 - x536*x619);
    res_0[52] = x281*(x538*x616 + x539*x618 + x540*x617 - x544*x619);
    res_0[53] = x281*(x545*x617 + x546*x618 + x547*x616 - x551*x619);
    res_0[54] = x281*(x372*x620 + x390*x621 + x519*x622 - x536*x623);
    res_0[55] = x281*(x538*x620 + x539*x622 + x540*x621 - x544*x623);
    res_0[56] = x281*(x545*x621 + x546*x622 + x547*x620 - x551*x623);
    res_0[57] = x281*(x372*x624 + x390*x625 + x519*x626 - x536*x627);
    res_0[58] = x281*(x538*x624 + x539*x626 + x540*x625 - x544*x627);
    res_0[59] = x281*(x545*x625 + x546*x626 + x547*x624 - x551*x627);
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
    
    double x0 = s*s;
    double x1 = (1.0/4.0)*x0;
    double x2 = x1 - 1.0/4.0;
    double x3 = (1.0/4.0)*r;
    double x4 = -x3;
    double x5 = x0*x3;
    double x6 = x4 + x5;
    double x7 = x2 + x6;
    double *x8 = coord;
    double x9 = x8[27];
    double x10 = r*r;
    double x11 = (1.0/4.0)*x10;
    double x12 = x11 - 1.0/4.0;
    double x13 = (1.0/4.0)*s;
    double x14 = -x13;
    double x15 = x10*x13;
    double x16 = x14 + x15;
    double x17 = x12 + x16;
    double x18 = x8[30];
    double x19 = x3 - x5;
    double x20 = x19 + x2;
    double x21 = x8[33];
    double x22 = x13 - x15;
    double x23 = x12 + x22;
    double x24 = x8[24];
    double x25 = 1.0/4.0 - x1;
    double x26 = x25 + x6;
    double x27 = x8[45];
    double x28 = 1.0/4.0 - x11;
    double x29 = x16 + x28;
    double x30 = x8[36];
    double x31 = x22 + x28;
    double x32 = x8[42];
    double x33 = x19 + x25;
    double x34 = x8[39];
    double x35 = (1.0/2.0)*t;
    double x36 = -x35;
    double x37 = r*x35;
    double x38 = x36 + x37;
    double x39 = s*x35;
    double x40 = s*x37;
    double x41 = -x40;
    double x42 = x39 + x41;
    double x43 = x38 + x42;
    double x44 = x8[48];
    double x45 = -x39;
    double x46 = x40 + x45;
    double x47 = x38 + x46;
    double x48 = x8[57];
    double x49 = -x37;
    double x50 = x36 + x49;
    double x51 = x39 + x40;
    double x52 = x50 + x51;
    double x53 = x8[51];
    double x54 = x41 + x45;
    double x55 = x50 + x54;
    double x56 = x8[54];
    double x57 = (1.0/4.0)*t;
    double x58 = t*x13;
    double x59 = -x58;
    double x60 = t*x3;
    double x61 = s*x60;
    double x62 = x59 + x61;
    double x63 = x57 + x62;
    double x64 = (1.0/8.0)*x10;
    double x65 = s*x64;
    double x66 = -x65;
    double x67 = (1.0/8.0)*r;
    double x68 = s*x67;
    double x69 = x66 + x68;
    double x70 = (1.0/8.0)*x0;
    double x71 = r*x70;
    double x72 = -x71;
    double x73 = x64 - 1.0/8.0;
    double x74 = -x60;
    double x75 = x70 + x74;
    double x76 = x72 + x73 + x75;
    double x77 = x63 + x69 + x76;
    double x78 = x8[12];
    double x79 = -x68;
    double x80 = x65 + x79;
    double x81 = -x61;
    double x82 = x58 + x81;
    double x83 = x57 + x82;
    double x84 = x76 + x80 + x83;
    double x85 = x8[21];
    double x86 = x59 + x81;
    double x87 = x73 + x86;
    double x88 = x57 + x66 + x79;
    double x89 = x60 + x70;
    double x90 = x71 + x89;
    double x91 = x87 + x88 + x90;
    double x92 = x8[15];
    double x93 = x57 + x65 + x68;
    double x94 = x58 + x61;
    double x95 = x73 + x94;
    double x96 = x90 + x93 + x95;
    double x97 = x8[18];
    double x98 = -x70;
    double x99 = 1.0/8.0 - x64;
    double x100 = x98 + x99;
    double x101 = x100 + x71 + x74;
    double x102 = x101 + x63 + x80;
    double x103 = x8[0];
    double x104 = x101 + x69 + x83;
    double x105 = x8[9];
    double x106 = x100 + x60 + x72;
    double x107 = x106 + x86 + x93;
    double x108 = x8[3];
    double x109 = x106 + x88 + x94;
    double x110 = x8[6];
    double x111 = x102*x103 + x104*x105 + x107*x108 + x109*x110 + x17*x18 + x20*x21 + x23*x24 + x26*x27 + x29*x30 + x31*x32 + x33*x34 + x43*x44 + x47*x48 + x52*x53 + x55*x56 + x7*x9 + x77*x78 + x84*x85 + x91*x92 + x96*x97;
    double x112 = t*t;
    double x113 = (1.0/4.0)*x112;
    double x114 = x113 - 1.0/4.0;
    double x115 = x112*x3;
    double x116 = x115 + x4;
    double x117 = x114 + x116;
    double x118 = x8[53];
    double x119 = -x57;
    double x120 = t*x11;
    double x121 = x119 + x120;
    double x122 = x12 + x121;
    double x123 = x8[38];
    double x124 = -x115 + x3;
    double x125 = x114 + x124;
    double x126 = x8[50];
    double x127 = -x120 + x57;
    double x128 = x12 + x127;
    double x129 = x8[26];
    double x130 = 1.0/4.0 - x113;
    double x131 = x116 + x130;
    double x132 = x8[59];
    double x133 = x121 + x28;
    double x134 = x8[32];
    double x135 = x127 + x28;
    double x136 = x8[44];
    double x137 = x124 + x130;
    double x138 = x8[56];
    double x139 = (1.0/2.0)*s;
    double x140 = -x139;
    double x141 = r*x139;
    double x142 = x140 + x141;
    double x143 = x142 + x42;
    double x144 = x8[35];
    double x145 = x142 + x46;
    double x146 = x8[47];
    double x147 = -x141;
    double x148 = x140 + x147;
    double x149 = x148 + x51;
    double x150 = x8[29];
    double x151 = x148 + x54;
    double x152 = x8[41];
    double x153 = t*x64;
    double x154 = x13 - x153;
    double x155 = (1.0/8.0)*x112;
    double x156 = r*x155;
    double x157 = -x156;
    double x158 = t*x67;
    double x159 = x157 + x158;
    double x160 = s*x3;
    double x161 = -x160;
    double x162 = x155 + x161;
    double x163 = x162 + x73;
    double x164 = x154 + x159 + x163 + x62;
    double x165 = x8[11];
    double x166 = x13 + x153;
    double x167 = -x158;
    double x168 = x157 + x167;
    double x169 = x163 + x166 + x168 + x82;
    double x170 = x8[23];
    double x171 = x156 + x167;
    double x172 = x154 + x160;
    double x173 = x155 + x171 + x172 + x87;
    double x174 = x8[8];
    double x175 = x156 + x158;
    double x176 = x160 + x166;
    double x177 = x155 + x175 + x176 + x95;
    double x178 = x8[20];
    double x179 = -x155;
    double x180 = x179 + x99;
    double x181 = x161 + x180;
    double x182 = x166 + x171 + x181 + x62;
    double x183 = x8[2];
    double x184 = x154 + x175 + x181 + x82;
    double x185 = x8[14];
    double x186 = x159 + x176 + x180 + x86;
    double x187 = x8[5];
    double x188 = x168 + x172 + x180 + x94;
    double x189 = x8[17];
    double x190 = x117*x118 + x122*x123 + x125*x126 + x128*x129 + x131*x132 + x133*x134 + x135*x136 + x137*x138 + x143*x144 + x145*x146 + x149*x150 + x151*x152 + x164*x165 + x169*x170 + x173*x174 + x177*x178 + x182*x183 + x184*x185 + x186*x187 + x188*x189;
    double x191 = x111*x190;
    double x192 = x102*x183 + x104*x165 + x107*x187 + x109*x174 + x118*x52 + x123*x29 + x126*x43 + x129*x23 + x132*x47 + x134*x17 + x136*x31 + x138*x55 + x144*x20 + x146*x26 + x150*x7 + x152*x33 + x170*x84 + x178*x96 + x185*x77 + x189*x91;
    double x193 = x103*x182 + x105*x164 + x108*x186 + x110*x173 + x117*x53 + x122*x30 + x125*x44 + x128*x24 + x131*x48 + x133*x18 + x135*x32 + x137*x56 + x143*x21 + x145*x27 + x149*x9 + x151*x34 + x169*x85 + x177*x97 + x184*x78 + x188*x92;
    double x194 = x192*x193;
    double x195 = x191 - x194;
    double x196 = x112*x13;
    double x197 = x14 + x196;
    double x198 = x114 + x197;
    double x199 = x8[58];
    double x200 = t*x1;
    double x201 = x119 + x200;
    double x202 = x2 + x201;
    double x203 = x8[46];
    double x204 = x13 - x196;
    double x205 = x114 + x204;
    double x206 = x8[49];
    double x207 = -x200 + x57;
    double x208 = x2 + x207;
    double x209 = x8[34];
    double x210 = x130 + x197;
    double x211 = x8[52];
    double x212 = x201 + x25;
    double x213 = x8[28];
    double x214 = x207 + x25;
    double x215 = x8[40];
    double x216 = x130 + x204;
    double x217 = x8[55];
    double x218 = -1.0/2.0*r;
    double x219 = x141 + x218;
    double x220 = x219 + x37 + x41;
    double x221 = x8[25];
    double x222 = x219 + x40 + x49;
    double x223 = x8[37];
    double x224 = x147 + x218;
    double x225 = x224 + x37 + x40;
    double x226 = x8[31];
    double x227 = x224 + x41 + x49;
    double x228 = x8[43];
    double x229 = x3 + x61;
    double x230 = t*x70;
    double x231 = -x230;
    double x232 = (1.0/8.0)*s*t;
    double x233 = x231 + x232;
    double x234 = s*x155;
    double x235 = -x234;
    double x236 = x162 + x235 - 1.0/8.0;
    double x237 = x229 + x233 + x236 + x75;
    double x238 = x8[4];
    double x239 = x230 + x89;
    double x240 = -x232;
    double x241 = x3 + x81;
    double x242 = x240 + x241;
    double x243 = x236 + x239 + x242;
    double x244 = x8[16];
    double x245 = x160 + x231;
    double x246 = x155 + x234 - 1.0/8.0;
    double x247 = x242 + x245 + x246 + x75;
    double x248 = x8[7];
    double x249 = x160 + x232;
    double x250 = x229 + x239 + x246 + x249;
    double x251 = x8[19];
    double x252 = x161 + x179 + x234 + x98 + 1.0/8.0;
    double x253 = x230 + x74;
    double x254 = x229 + x240;
    double x255 = x252 + x253 + x254;
    double x256 = x8[1];
    double x257 = x233 + x241 + x252 + x60;
    double x258 = x8[13];
    double x259 = x179 + x235 + x98 + 1.0/8.0;
    double x260 = x241 + x249 + x253 + x259;
    double x261 = x8[10];
    double x262 = x245 + x254 + x259 + x60;
    double x263 = x8[22];
    double x264 = x198*x199 + x202*x203 + x205*x206 + x208*x209 + x210*x211 + x212*x213 + x214*x215 + x216*x217 + x220*x221 + x222*x223 + x225*x226 + x227*x228 + x237*x238 + x243*x244 + x247*x248 + x250*x251 + x255*x256 + x257*x258 + x260*x261 + x262*x263;
    double x265 = x102*x256 + x104*x261 + x107*x238 + x109*x248 + x17*x226 + x199*x47 + x20*x209 + x203*x26 + x206*x43 + x211*x52 + x213*x7 + x215*x33 + x217*x55 + x221*x23 + x223*x29 + x228*x31 + x244*x91 + x251*x96 + x258*x77 + x263*x84;
    double x266 = x118*x210 + x123*x222 + x126*x205 + x129*x220 + x132*x198 + x134*x225 + x136*x227 + x138*x216 + x144*x208 + x146*x202 + x150*x212 + x152*x214 + x165*x260 + x170*x262 + x174*x247 + x178*x250 + x183*x255 + x185*x257 + x187*x237 + x189*x243;
    double x267 = x193*x266;
    double x268 = x117*x211 + x122*x223 + x125*x206 + x128*x221 + x131*x199 + x133*x226 + x135*x228 + x137*x217 + x143*x209 + x145*x203 + x149*x213 + x151*x215 + x164*x261 + x169*x263 + x173*x248 + x177*x251 + x182*x256 + x184*x258 + x186*x238 + x188*x244;
    double x269 = x103*x255 + x105*x260 + x108*x237 + x110*x247 + x18*x225 + x198*x48 + x202*x27 + x205*x44 + x208*x21 + x210*x53 + x212*x9 + x214*x34 + x216*x56 + x220*x24 + x222*x30 + x227*x32 + x243*x92 + x250*x97 + x257*x78 + x262*x85;
    double x270 = x192*x269;
    double x271 = x111*x266;
    double x272 = x190*x269;
    double x273 = x191*x264 - x194*x264 + x265*x267 - x265*x272 + x268*x270 - x268*x271;
    double x274 = 1.0/x273;
    double *x275 = U;
    double x276 = x275[59];
    double x277 = x275[47];
    double x278 = x275[50];
    double x279 = x275[35];
    double x280 = x275[53];
    double x281 = x275[29];
    double x282 = x275[41];
    double x283 = x275[56];
    double x284 = x275[26];
    double x285 = x275[38];
    double x286 = x275[32];
    double x287 = x275[44];
    double x288 = x275[5];
    double x289 = x275[17];
    double x290 = x275[8];
    double x291 = x275[20];
    double x292 = x275[2];
    double x293 = x275[14];
    double x294 = x275[11];
    double x295 = x275[23];
    double x296 = x274*(x198*x276 + x202*x277 + x205*x278 + x208*x279 + x210*x280 + x212*x281 + x214*x282 + x216*x283 + x220*x284 + x222*x285 + x225*x286 + x227*x287 + x237*x288 + x243*x289 + x247*x290 + x250*x291 + x255*x292 + x257*x293 + x260*x294 + x262*x295);
    double x297 = x270 - x271;
    double x298 = x274*(x117*x280 + x122*x285 + x125*x278 + x128*x284 + x131*x276 + x133*x286 + x135*x287 + x137*x283 + x143*x279 + x145*x277 + x149*x281 + x151*x282 + x164*x294 + x169*x295 + x173*x290 + x177*x291 + x182*x292 + x184*x293 + x186*x288 + x188*x289);
    double x299 = x267 - x272;
    double x300 = x274*(x102*x292 + x104*x294 + x107*x288 + x109*x290 + x17*x286 + x20*x279 + x23*x284 + x26*x277 + x276*x47 + x278*x43 + x280*x52 + x281*x7 + x282*x33 + x283*x55 + x285*x29 + x287*x31 + x289*x91 + x291*x96 + x293*x77 + x295*x84);
    double x301 = x195*x296 + x297*x298 + x299*x300;
    double x302 = -x192*x264 + x265*x266;
    double x303 = x275[52];
    double x304 = x275[37];
    double x305 = x275[49];
    double x306 = x275[25];
    double x307 = x275[58];
    double x308 = x275[31];
    double x309 = x275[43];
    double x310 = x275[55];
    double x311 = x275[34];
    double x312 = x275[46];
    double x313 = x275[28];
    double x314 = x275[40];
    double x315 = x275[10];
    double x316 = x275[22];
    double x317 = x275[7];
    double x318 = x275[19];
    double x319 = x275[1];
    double x320 = x275[13];
    double x321 = x275[4];
    double x322 = x275[16];
    double x323 = x274*(x117*x303 + x122*x304 + x125*x305 + x128*x306 + x131*x307 + x133*x308 + x135*x309 + x137*x310 + x143*x311 + x145*x312 + x149*x313 + x151*x314 + x164*x315 + x169*x316 + x173*x317 + x177*x318 + x182*x319 + x184*x320 + x186*x321 + x188*x322);
    double x324 = -x190*x265 + x192*x268;
    double x325 = x274*(x198*x307 + x202*x312 + x205*x305 + x208*x311 + x210*x303 + x212*x313 + x214*x314 + x216*x310 + x220*x306 + x222*x304 + x225*x308 + x227*x309 + x237*x321 + x243*x322 + x247*x317 + x250*x318 + x255*x319 + x257*x320 + x260*x315 + x262*x316);
    double x326 = x190*x264 - x266*x268;
    double x327 = x274*(x102*x319 + x104*x315 + x107*x321 + x109*x317 + x17*x308 + x20*x311 + x23*x306 + x26*x312 + x29*x304 + x303*x52 + x305*x43 + x307*x47 + x309*x31 + x310*x55 + x313*x7 + x314*x33 + x316*x84 + x318*x96 + x320*x77 + x322*x91);
    double x328 = x302*x323 + x324*x325 + x326*x327;
    double x329 = x296*x324 + x298*x302 + x300*x326;
    double x330 = x195*x325 + x297*x323 + x299*x327 + 1.0;
    double x331 = x301*x328 - x329*x330;
    double x332 = x275[57];
    double x333 = x275[45];
    double x334 = x275[48];
    double x335 = x275[33];
    double x336 = x275[51];
    double x337 = x275[27];
    double x338 = x275[39];
    double x339 = x275[54];
    double x340 = x275[24];
    double x341 = x275[36];
    double x342 = x275[30];
    double x343 = x275[42];
    double x344 = x275[3];
    double x345 = x275[15];
    double x346 = x275[6];
    double x347 = x275[18];
    double x348 = x275[0];
    double x349 = x275[12];
    double x350 = x275[9];
    double x351 = x275[21];
    double x352 = x274*(x198*x332 + x202*x333 + x205*x334 + x208*x335 + x210*x336 + x212*x337 + x214*x338 + x216*x339 + x220*x340 + x222*x341 + x225*x342 + x227*x343 + x237*x344 + x243*x345 + x247*x346 + x250*x347 + x255*x348 + x257*x349 + x260*x350 + x262*x351);
    double x353 = x274*(x117*x336 + x122*x341 + x125*x334 + x128*x340 + x131*x332 + x133*x342 + x135*x343 + x137*x339 + x143*x335 + x145*x333 + x149*x337 + x151*x338 + x164*x350 + x169*x351 + x173*x346 + x177*x347 + x182*x348 + x184*x349 + x186*x344 + x188*x345);
    double x354 = x274*(x102*x348 + x104*x350 + x107*x344 + x109*x346 + x17*x342 + x20*x335 + x23*x340 + x26*x333 + x29*x341 + x31*x343 + x33*x338 + x332*x47 + x334*x43 + x336*x52 + x337*x7 + x339*x55 + x345*x91 + x347*x96 + x349*x77 + x351*x84);
    double x355 = x195*x352 + x297*x353 + x299*x354;
    double x356 = x302*x353 + x324*x352 + x326*x354 + 1.0;
    double x357 = -x301*x356 + x329*x355;
    double x358 = -x328*x355 + x330*x356;
    double x359 = S11*x331 + S12*x357 + S13*x358;
    double x360 = S12*x331 + S22*x357 + S23*x358;
    double x361 = S13*x331 + S23*x357 + S33*x358;
    double x362 = -x328*x360 - x329*x361 - x356*x359;
    double x363 = x111*x264 - x265*x269;
    double x364 = x182*x274;
    double x365 = -x111*x268 + x193*x265;
    double x366 = x255*x274;
    double x367 = -x193*x264 + x268*x269;
    double x368 = x102*x274;
    double x369 = x363*x364 + x365*x366 + x367*x368;
    double x370 = x323*x363 + x325*x365 + x327*x367;
    double x371 = x296*x365 + x298*x363 + x300*x367 + 1.0;
    double x372 = -x328*x371 + x329*x370;
    double x373 = x352*x365 + x353*x363 + x354*x367;
    double x374 = -x329*x373 + x356*x371;
    double x375 = x328*x373 - x356*x370;
    double x376 = S11*x372 + S12*x374 + S13*x375;
    double x377 = S12*x372 + S22*x374 + S23*x375;
    double x378 = S13*x372 + S23*x374 + S33*x375;
    double x379 = -x328*x377 - x329*x378 - x356*x376;
    double x380 = x195*x366 + x297*x364 + x299*x368;
    double x381 = x302*x364 + x324*x366 + x326*x368;
    double x382 = -x301*x370 + x330*x371;
    double x383 = x301*x373 - x355*x371;
    double x384 = -x330*x373 + x355*x370;
    double x385 = S11*x382 + S12*x383 + S13*x384;
    double x386 = S12*x382 + S22*x383 + S23*x384;
    double x387 = S13*x382 + S23*x383 + S33*x384;
    double x388 = 1.0*PENER + 1.0*SENER;
    double x389 = -x328*x386 - x329*x387 - x356*x385 + x388;
    double x390 = -x301*x361 - x330*x360 - x355*x359;
    double x391 = -x301*x387 - x330*x386 - x355*x385;
    double x392 = -x301*x378 - x330*x377 - x355*x376 + x388;
    double x393 = -x370*x377 - x371*x378 - x373*x376;
    double x394 = -x370*x386 - x371*x387 - x373*x385;
    double x395 = -x359*x373 - x360*x370 - x361*x371 + x388;
    double x396 = x186*x274;
    double x397 = x237*x274;
    double x398 = x107*x274;
    double x399 = x363*x396 + x365*x397 + x367*x398;
    double x400 = x195*x397 + x297*x396 + x299*x398;
    double x401 = x302*x396 + x324*x397 + x326*x398;
    double x402 = x173*x274;
    double x403 = x247*x274;
    double x404 = x109*x274;
    double x405 = x363*x402 + x365*x403 + x367*x404;
    double x406 = x195*x403 + x297*x402 + x299*x404;
    double x407 = x302*x402 + x324*x403 + x326*x404;
    double x408 = x164*x274;
    double x409 = x260*x274;
    double x410 = x104*x274;
    double x411 = x363*x408 + x365*x409 + x367*x410;
    double x412 = x195*x409 + x297*x408 + x299*x410;
    double x413 = x302*x408 + x324*x409 + x326*x410;
    double x414 = x184*x274;
    double x415 = x257*x274;
    double x416 = x274*x77;
    double x417 = x363*x414 + x365*x415 + x367*x416;
    double x418 = x195*x415 + x297*x414 + x299*x416;
    double x419 = x302*x414 + x324*x415 + x326*x416;
    double x420 = x188*x274;
    double x421 = x243*x274;
    double x422 = x274*x91;
    double x423 = x363*x420 + x365*x421 + x367*x422;
    double x424 = x195*x421 + x297*x420 + x299*x422;
    double x425 = x302*x420 + x324*x421 + x326*x422;
    double x426 = x177*x274;
    double x427 = x250*x274;
    double x428 = x274*x96;
    double x429 = x363*x426 + x365*x427 + x367*x428;
    double x430 = x195*x427 + x297*x426 + x299*x428;
    double x431 = x302*x426 + x324*x427 + x326*x428;
    double x432 = x169*x274;
    double x433 = x262*x274;
    double x434 = x274*x84;
    double x435 = x363*x432 + x365*x433 + x367*x434;
    double x436 = x195*x433 + x297*x432 + x299*x434;
    double x437 = x302*x432 + x324*x433 + x326*x434;
    double x438 = x128*x274;
    double x439 = x220*x274;
    double x440 = x23*x274;
    double x441 = x363*x438 + x365*x439 + x367*x440;
    double x442 = x195*x439 + x297*x438 + x299*x440;
    double x443 = x302*x438 + x324*x439 + x326*x440;
    double x444 = x149*x274;
    double x445 = x212*x274;
    double x446 = x274*x7;
    double x447 = x363*x444 + x365*x445 + x367*x446;
    double x448 = x195*x445 + x297*x444 + x299*x446;
    double x449 = x302*x444 + x324*x445 + x326*x446;
    double x450 = x133*x274;
    double x451 = x225*x274;
    double x452 = x17*x274;
    double x453 = x363*x450 + x365*x451 + x367*x452;
    double x454 = x195*x451 + x297*x450 + x299*x452;
    double x455 = x302*x450 + x324*x451 + x326*x452;
    double x456 = x143*x274;
    double x457 = x208*x274;
    double x458 = x20*x274;
    double x459 = x363*x456 + x365*x457 + x367*x458;
    double x460 = x195*x457 + x297*x456 + x299*x458;
    double x461 = x302*x456 + x324*x457 + x326*x458;
    double x462 = x122*x274;
    double x463 = x222*x274;
    double x464 = x274*x29;
    double x465 = x363*x462 + x365*x463 + x367*x464;
    double x466 = x195*x463 + x297*x462 + x299*x464;
    double x467 = x302*x462 + x324*x463 + x326*x464;
    double x468 = x151*x274;
    double x469 = x214*x274;
    double x470 = x274*x33;
    double x471 = x363*x468 + x365*x469 + x367*x470;
    double x472 = x195*x469 + x297*x468 + x299*x470;
    double x473 = x302*x468 + x324*x469 + x326*x470;
    double x474 = x135*x274;
    double x475 = x227*x274;
    double x476 = x274*x31;
    double x477 = x363*x474 + x365*x475 + x367*x476;
    double x478 = x195*x475 + x297*x474 + x299*x476;
    double x479 = x302*x474 + x324*x475 + x326*x476;
    double x480 = x145*x274;
    double x481 = x202*x274;
    double x482 = x26*x274;
    double x483 = x363*x480 + x365*x481 + x367*x482;
    double x484 = x195*x481 + x297*x480 + x299*x482;
    double x485 = x302*x480 + x324*x481 + x326*x482;
    double x486 = x125*x274;
    double x487 = x205*x274;
    double x488 = x274*x43;
    double x489 = x363*x486 + x365*x487 + x367*x488;
    double x490 = x195*x487 + x297*x486 + x299*x488;
    double x491 = x302*x486 + x324*x487 + x326*x488;
    double x492 = x117*x274;
    double x493 = x210*x274;
    double x494 = x274*x52;
    double x495 = x363*x492 + x365*x493 + x367*x494;
    double x496 = x195*x493 + x297*x492 + x299*x494;
    double x497 = x302*x492 + x324*x493 + x326*x494;
    double x498 = x137*x274;
    double x499 = x216*x274;
    double x500 = x274*x55;
    double x501 = x363*x498 + x365*x499 + x367*x500;
    double x502 = x195*x499 + x297*x498 + x299*x500;
    double x503 = x302*x498 + x324*x499 + x326*x500;
    double x504 = x131*x274;
    double x505 = x198*x274;
    double x506 = x274*x47;
    double x507 = x363*x504 + x365*x505 + x367*x506;
    double x508 = x195*x505 + x297*x504 + x299*x506;
    double x509 = x302*x504 + x324*x505 + x326*x506;
    
    res_0[0] = x273*(x362*x369 + x379*x380 + x381*x389);
    res_0[1] = x273*(x369*x390 + x380*x392 + x381*x391);
    res_0[2] = x273*(x369*x395 + x380*x393 + x381*x394);
    res_0[3] = x273*(x362*x399 + x379*x400 + x389*x401);
    res_0[4] = x273*(x390*x399 + x391*x401 + x392*x400);
    res_0[5] = x273*(x393*x400 + x394*x401 + x395*x399);
    res_0[6] = x273*(x362*x405 + x379*x406 + x389*x407);
    res_0[7] = x273*(x390*x405 + x391*x407 + x392*x406);
    res_0[8] = x273*(x393*x406 + x394*x407 + x395*x405);
    res_0[9] = x273*(x362*x411 + x379*x412 + x389*x413);
    res_0[10] = x273*(x390*x411 + x391*x413 + x392*x412);
    res_0[11] = x273*(x393*x412 + x394*x413 + x395*x411);
    res_0[12] = x273*(x362*x417 + x379*x418 + x389*x419);
    res_0[13] = x273*(x390*x417 + x391*x419 + x392*x418);
    res_0[14] = x273*(x393*x418 + x394*x419 + x395*x417);
    res_0[15] = x273*(x362*x423 + x379*x424 + x389*x425);
    res_0[16] = x273*(x390*x423 + x391*x425 + x392*x424);
    res_0[17] = x273*(x393*x424 + x394*x425 + x395*x423);
    res_0[18] = x273*(x362*x429 + x379*x430 + x389*x431);
    res_0[19] = x273*(x390*x429 + x391*x431 + x392*x430);
    res_0[20] = x273*(x393*x430 + x394*x431 + x395*x429);
    res_0[21] = x273*(x362*x435 + x379*x436 + x389*x437);
    res_0[22] = x273*(x390*x435 + x391*x437 + x392*x436);
    res_0[23] = x273*(x393*x436 + x394*x437 + x395*x435);
    res_0[24] = x273*(x362*x441 + x379*x442 + x389*x443);
    res_0[25] = x273*(x390*x441 + x391*x443 + x392*x442);
    res_0[26] = x273*(x393*x442 + x394*x443 + x395*x441);
    res_0[27] = x273*(x362*x447 + x379*x448 + x389*x449);
    res_0[28] = x273*(x390*x447 + x391*x449 + x392*x448);
    res_0[29] = x273*(x393*x448 + x394*x449 + x395*x447);
    res_0[30] = x273*(x362*x453 + x379*x454 + x389*x455);
    res_0[31] = x273*(x390*x453 + x391*x455 + x392*x454);
    res_0[32] = x273*(x393*x454 + x394*x455 + x395*x453);
    res_0[33] = x273*(x362*x459 + x379*x460 + x389*x461);
    res_0[34] = x273*(x390*x459 + x391*x461 + x392*x460);
    res_0[35] = x273*(x393*x460 + x394*x461 + x395*x459);
    res_0[36] = x273*(x362*x465 + x379*x466 + x389*x467);
    res_0[37] = x273*(x390*x465 + x391*x467 + x392*x466);
    res_0[38] = x273*(x393*x466 + x394*x467 + x395*x465);
    res_0[39] = x273*(x362*x471 + x379*x472 + x389*x473);
    res_0[40] = x273*(x390*x471 + x391*x473 + x392*x472);
    res_0[41] = x273*(x393*x472 + x394*x473 + x395*x471);
    res_0[42] = x273*(x362*x477 + x379*x478 + x389*x479);
    res_0[43] = x273*(x390*x477 + x391*x479 + x392*x478);
    res_0[44] = x273*(x393*x478 + x394*x479 + x395*x477);
    res_0[45] = x273*(x362*x483 + x379*x484 + x389*x485);
    res_0[46] = x273*(x390*x483 + x391*x485 + x392*x484);
    res_0[47] = x273*(x393*x484 + x394*x485 + x395*x483);
    res_0[48] = x273*(x362*x489 + x379*x490 + x389*x491);
    res_0[49] = x273*(x390*x489 + x391*x491 + x392*x490);
    res_0[50] = x273*(x393*x490 + x394*x491 + x395*x489);
    res_0[51] = x273*(x362*x495 + x379*x496 + x389*x497);
    res_0[52] = x273*(x390*x495 + x391*x497 + x392*x496);
    res_0[53] = x273*(x393*x496 + x394*x497 + x395*x495);
    res_0[54] = x273*(x362*x501 + x379*x502 + x389*x503);
    res_0[55] = x273*(x390*x501 + x391*x503 + x392*x502);
    res_0[56] = x273*(x393*x502 + x394*x503 + x395*x501);
    res_0[57] = x273*(x362*x507 + x379*x508 + x389*x509);
    res_0[58] = x273*(x390*x507 + x391*x509 + x392*x508);
    res_0[59] = x273*(x393*x508 + x394*x509 + x395*x507);
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
    
    double x0 = s*s;
    double x1 = (1.0/4.0)*x0;
    double x2 = x1 - 1.0/4.0;
    double x3 = (1.0/4.0)*r;
    double x4 = -x3;
    double x5 = x0*x3;
    double x6 = x4 + x5;
    double x7 = x2 + x6;
    double *x8 = coord;
    double x9 = x8[27];
    double x10 = r*r;
    double x11 = (1.0/4.0)*x10;
    double x12 = x11 - 1.0/4.0;
    double x13 = (1.0/4.0)*s;
    double x14 = -x13;
    double x15 = x10*x13;
    double x16 = x14 + x15;
    double x17 = x12 + x16;
    double x18 = x8[30];
    double x19 = x3 - x5;
    double x20 = x19 + x2;
    double x21 = x8[33];
    double x22 = x13 - x15;
    double x23 = x12 + x22;
    double x24 = x8[24];
    double x25 = 1.0/4.0 - x1;
    double x26 = x25 + x6;
    double x27 = x8[45];
    double x28 = 1.0/4.0 - x11;
    double x29 = x16 + x28;
    double x30 = x8[36];
    double x31 = x22 + x28;
    double x32 = x8[42];
    double x33 = x19 + x25;
    double x34 = x8[39];
    double x35 = (1.0/2.0)*t;
    double x36 = -x35;
    double x37 = r*x35;
    double x38 = x36 + x37;
    double x39 = s*x35;
    double x40 = s*x37;
    double x41 = -x40;
    double x42 = x39 + x41;
    double x43 = x38 + x42;
    double x44 = x8[48];
    double x45 = -x39;
    double x46 = x40 + x45;
    double x47 = x38 + x46;
    double x48 = x8[57];
    double x49 = -x37;
    double x50 = x36 + x49;
    double x51 = x39 + x40;
    double x52 = x50 + x51;
    double x53 = x8[51];
    double x54 = x41 + x45;
    double x55 = x50 + x54;
    double x56 = x8[54];
    double x57 = (1.0/4.0)*t;
    double x58 = t*x13;
    double x59 = -x58;
    double x60 = t*x3;
    double x61 = s*x60;
    double x62 = x59 + x61;
    double x63 = x57 + x62;
    double x64 = (1.0/8.0)*x10;
    double x65 = s*x64;
    double x66 = -x65;
    double x67 = (1.0/8.0)*r;
    double x68 = s*x67;
    double x69 = x66 + x68;
    double x70 = (1.0/8.0)*x0;
    double x71 = r*x70;
    double x72 = -x71;
    double x73 = x64 - 1.0/8.0;
    double x74 = -x60;
    double x75 = x70 + x74;
    double x76 = x72 + x73 + x75;
    double x77 = x63 + x69 + x76;
    double x78 = x8[12];
    double x79 = -x68;
    double x80 = x65 + x79;
    double x81 = -x61;
    double x82 = x58 + x81;
    double x83 = x57 + x82;
    double x84 = x76 + x80 + x83;
    double x85 = x8[21];
    double x86 = x59 + x81;
    double x87 = x73 + x86;
    double x88 = x57 + x66 + x79;
    double x89 = x60 + x70;
    double x90 = x71 + x89;
    double x91 = x87 + x88 + x90;
    double x92 = x8[15];
    double x93 = x57 + x65 + x68;
    double x94 = x58 + x61;
    double x95 = x73 + x94;
    double x96 = x90 + x93 + x95;
    double x97 = x8[18];
    double x98 = -x70;
    double x99 = 1.0/8.0 - x64;
    double x100 = x98 + x99;
    double x101 = x100 + x71 + x74;
    double x102 = x101 + x63 + x80;
    double x103 = x8[0];
    double x104 = x101 + x69 + x83;
    double x105 = x8[9];
    double x106 = x100 + x60 + x72;
    double x107 = x106 + x86 + x93;
    double x108 = x8[3];
    double x109 = x106 + x88 + x94;
    double x110 = x8[6];
    double x111 = x102*x103 + x104*x105 + x107*x108 + x109*x110 + x17*x18 + x20*x21 + x23*x24 + x26*x27 + x29*x30 + x31*x32 + x33*x34 + x43*x44 + x47*x48 + x52*x53 + x55*x56 + x7*x9 + x77*x78 + x84*x85 + x91*x92 + x96*x97;
    double x112 = t*t;
    double x113 = (1.0/4.0)*x112;
    double x114 = x113 - 1.0/4.0;
    double x115 = x112*x3;
    double x116 = x115 + x4;
    double x117 = x114 + x116;
    double x118 = x8[53];
    double x119 = -x57;
    double x120 = t*x11;
    double x121 = x119 + x120;
    double x122 = x12 + x121;
    double x123 = x8[38];
    double x124 = -x115 + x3;
    double x125 = x114 + x124;
    double x126 = x8[50];
    double x127 = -x120 + x57;
    double x128 = x12 + x127;
    double x129 = x8[26];
    double x130 = 1.0/4.0 - x113;
    double x131 = x116 + x130;
    double x132 = x8[59];
    double x133 = x121 + x28;
    double x134 = x8[32];
    double x135 = x127 + x28;
    double x136 = x8[44];
    double x137 = x124 + x130;
    double x138 = x8[56];
    double x139 = (1.0/2.0)*s;
    double x140 = -x139;
    double x141 = r*x139;
    double x142 = x140 + x141;
    double x143 = x142 + x42;
    double x144 = x8[35];
    double x145 = x142 + x46;
    double x146 = x8[47];
    double x147 = -x141;
    double x148 = x140 + x147;
    double x149 = x148 + x51;
    double x150 = x8[29];
    double x151 = x148 + x54;
    double x152 = x8[41];
    double x153 = t*x64;
    double x154 = x13 - x153;
    double x155 = (1.0/8.0)*x112;
    double x156 = r*x155;
    double x157 = -x156;
    double x158 = t*x67;
    double x159 = x157 + x158;
    double x160 = s*x3;
    double x161 = -x160;
    double x162 = x155 + x161;
    double x163 = x162 + x73;
    double x164 = x154 + x159 + x163 + x62;
    double x165 = x8[11];
    double x166 = x13 + x153;
    double x167 = -x158;
    double x168 = x157 + x167;
    double x169 = x163 + x166 + x168 + x82;
    double x170 = x8[23];
    double x171 = x156 + x167;
    double x172 = x154 + x160;
    double x173 = x155 + x171 + x172 + x87;
    double x174 = x8[8];
    double x175 = x156 + x158;
    double x176 = x160 + x166;
    double x177 = x155 + x175 + x176 + x95;
    double x178 = x8[20];
    double x179 = -x155;
    double x180 = x179 + x99;
    double x181 = x161 + x180;
    double x182 = x166 + x171 + x181 + x62;
    double x183 = x8[2];
    double x184 = x154 + x175 + x181 + x82;
    double x185 = x8[14];
    double x186 = x159 + x176 + x180 + x86;
    double x187 = x8[5];
    double x188 = x168 + x172 + x180 + x94;
    double x189 = x8[17];
    double x190 = x117*x118 + x122*x123 + x125*x126 + x128*x129 + x131*x132 + x133*x134 + x135*x136 + x137*x138 + x143*x144 + x145*x146 + x149*x150 + x151*x152 + x164*x165 + x169*x170 + x173*x174 + x177*x178 + x182*x183 + x184*x185 + x186*x187 + x188*x189;
    double x191 = x111*x190;
    double x192 = x102*x183 + x104*x165 + x107*x187 + x109*x174 + x118*x52 + x123*x29 + x126*x43 + x129*x23 + x132*x47 + x134*x17 + x136*x31 + x138*x55 + x144*x20 + x146*x26 + x150*x7 + x152*x33 + x170*x84 + x178*x96 + x185*x77 + x189*x91;
    double x193 = x103*x182 + x105*x164 + x108*x186 + x110*x173 + x117*x53 + x122*x30 + x125*x44 + x128*x24 + x131*x48 + x133*x18 + x135*x32 + x137*x56 + x143*x21 + x145*x27 + x149*x9 + x151*x34 + x169*x85 + x177*x97 + x184*x78 + x188*x92;
    double x194 = x192*x193;
    double x195 = x191 - x194;
    double x196 = x112*x13;
    double x197 = x14 + x196;
    double x198 = x114 + x197;
    double x199 = x8[58];
    double x200 = t*x1;
    double x201 = x119 + x200;
    double x202 = x2 + x201;
    double x203 = x8[46];
    double x204 = x13 - x196;
    double x205 = x114 + x204;
    double x206 = x8[49];
    double x207 = -x200 + x57;
    double x208 = x2 + x207;
    double x209 = x8[34];
    double x210 = x130 + x197;
    double x211 = x8[52];
    double x212 = x201 + x25;
    double x213 = x8[28];
    double x214 = x207 + x25;
    double x215 = x8[40];
    double x216 = x130 + x204;
    double x217 = x8[55];
    double x218 = -1.0/2.0*r;
    double x219 = x141 + x218;
    double x220 = x219 + x37 + x41;
    double x221 = x8[25];
    double x222 = x219 + x40 + x49;
    double x223 = x8[37];
    double x224 = x147 + x218;
    double x225 = x224 + x37 + x40;
    double x226 = x8[31];
    double x227 = x224 + x41 + x49;
    double x228 = x8[43];
    double x229 = x3 + x61;
    double x230 = t*x70;
    double x231 = -x230;
    double x232 = (1.0/8.0)*s*t;
    double x233 = x231 + x232;
    double x234 = s*x155;
    double x235 = -x234;
    double x236 = x162 + x235 - 1.0/8.0;
    double x237 = x229 + x233 + x236 + x75;
    double x238 = x8[4];
    double x239 = x230 + x89;
    double x240 = -x232;
    double x241 = x3 + x81;
    double x242 = x240 + x241;
    double x243 = x236 + x239 + x242;
    double x244 = x8[16];
    double x245 = x160 + x231;
    double x246 = x155 + x234 - 1.0/8.0;
    double x247 = x242 + x245 + x246 + x75;
    double x248 = x8[7];
    double x249 = x160 + x232;
    double x250 = x229 + x239 + x246 + x249;
    double x251 = x8[19];
    double x252 = x161 + x179 + x234 + x98 + 1.0/8.0;
    double x253 = x230 + x74;
    double x254 = x229 + x240;
    double x255 = x252 + x253 + x254;
    double x256 = x8[1];
    double x257 = x233 + x241 + x252 + x60;
    double x258 = x8[13];
    double x259 = x179 + x235 + x98 + 1.0/8.0;
    double x260 = x241 + x249 + x253 + x259;
    double x261 = x8[10];
    double x262 = x245 + x254 + x259 + x60;
    double x263 = x8[22];
    double x264 = x198*x199 + x202*x203 + x205*x206 + x208*x209 + x210*x211 + x212*x213 + x214*x215 + x216*x217 + x220*x221 + x222*x223 + x225*x226 + x227*x228 + x237*x238 + x243*x244 + x247*x248 + x250*x251 + x255*x256 + x257*x258 + x260*x261 + x262*x263;
    double x265 = x102*x256 + x104*x261 + x107*x238 + x109*x248 + x17*x226 + x199*x47 + x20*x209 + x203*x26 + x206*x43 + x211*x52 + x213*x7 + x215*x33 + x217*x55 + x221*x23 + x223*x29 + x228*x31 + x244*x91 + x251*x96 + x258*x77 + x263*x84;
    double x266 = x118*x210 + x123*x222 + x126*x205 + x129*x220 + x132*x198 + x134*x225 + x136*x227 + x138*x216 + x144*x208 + x146*x202 + x150*x212 + x152*x214 + x165*x260 + x170*x262 + x174*x247 + x178*x250 + x183*x255 + x185*x257 + x187*x237 + x189*x243;
    double x267 = x193*x266;
    double x268 = x117*x211 + x122*x223 + x125*x206 + x128*x221 + x131*x199 + x133*x226 + x135*x228 + x137*x217 + x143*x209 + x145*x203 + x149*x213 + x151*x215 + x164*x261 + x169*x263 + x173*x248 + x177*x251 + x182*x256 + x184*x258 + x186*x238 + x188*x244;
    double x269 = x103*x255 + x105*x260 + x108*x237 + x110*x247 + x18*x225 + x198*x48 + x202*x27 + x205*x44 + x208*x21 + x210*x53 + x212*x9 + x214*x34 + x216*x56 + x220*x24 + x222*x30 + x227*x32 + x243*x92 + x250*x97 + x257*x78 + x262*x85;
    double x270 = x192*x269;
    double x271 = x111*x266;
    double x272 = x190*x269;
    double x273 = x191*x264 - x194*x264 + x265*x267 - x265*x272 + x268*x270 - x268*x271;
    double x274 = 1.0/x273;
    double *x275 = U;
    double x276 = x275[59];
    double x277 = x275[47];
    double x278 = x275[50];
    double x279 = x275[35];
    double x280 = x275[53];
    double x281 = x275[29];
    double x282 = x275[41];
    double x283 = x275[56];
    double x284 = x275[26];
    double x285 = x275[38];
    double x286 = x275[32];
    double x287 = x275[44];
    double x288 = x275[5];
    double x289 = x275[17];
    double x290 = x275[8];
    double x291 = x275[20];
    double x292 = x275[2];
    double x293 = x275[14];
    double x294 = x275[11];
    double x295 = x275[23];
    double x296 = x274*(x198*x276 + x202*x277 + x205*x278 + x208*x279 + x210*x280 + x212*x281 + x214*x282 + x216*x283 + x220*x284 + x222*x285 + x225*x286 + x227*x287 + x237*x288 + x243*x289 + x247*x290 + x250*x291 + x255*x292 + x257*x293 + x260*x294 + x262*x295);
    double x297 = x270 - x271;
    double x298 = x274*(x117*x280 + x122*x285 + x125*x278 + x128*x284 + x131*x276 + x133*x286 + x135*x287 + x137*x283 + x143*x279 + x145*x277 + x149*x281 + x151*x282 + x164*x294 + x169*x295 + x173*x290 + x177*x291 + x182*x292 + x184*x293 + x186*x288 + x188*x289);
    double x299 = x267 - x272;
    double x300 = x274*(x102*x292 + x104*x294 + x107*x288 + x109*x290 + x17*x286 + x20*x279 + x23*x284 + x26*x277 + x276*x47 + x278*x43 + x280*x52 + x281*x7 + x282*x33 + x283*x55 + x285*x29 + x287*x31 + x289*x91 + x291*x96 + x293*x77 + x295*x84);
    double x301 = x195*x296 + x297*x298 + x299*x300;
    double x302 = -x192*x264 + x265*x266;
    double x303 = x275[52];
    double x304 = x275[37];
    double x305 = x275[49];
    double x306 = x275[25];
    double x307 = x275[58];
    double x308 = x275[31];
    double x309 = x275[43];
    double x310 = x275[55];
    double x311 = x275[34];
    double x312 = x275[46];
    double x313 = x275[28];
    double x314 = x275[40];
    double x315 = x275[10];
    double x316 = x275[22];
    double x317 = x275[7];
    double x318 = x275[19];
    double x319 = x275[1];
    double x320 = x275[13];
    double x321 = x275[4];
    double x322 = x275[16];
    double x323 = x274*(x117*x303 + x122*x304 + x125*x305 + x128*x306 + x131*x307 + x133*x308 + x135*x309 + x137*x310 + x143*x311 + x145*x312 + x149*x313 + x151*x314 + x164*x315 + x169*x316 + x173*x317 + x177*x318 + x182*x319 + x184*x320 + x186*x321 + x188*x322);
    double x324 = -x190*x265 + x192*x268;
    double x325 = x274*(x198*x307 + x202*x312 + x205*x305 + x208*x311 + x210*x303 + x212*x313 + x214*x314 + x216*x310 + x220*x306 + x222*x304 + x225*x308 + x227*x309 + x237*x321 + x243*x322 + x247*x317 + x250*x318 + x255*x319 + x257*x320 + x260*x315 + x262*x316);
    double x326 = x190*x264 - x266*x268;
    double x327 = x274*(x102*x319 + x104*x315 + x107*x321 + x109*x317 + x17*x308 + x20*x311 + x23*x306 + x26*x312 + x29*x304 + x303*x52 + x305*x43 + x307*x47 + x309*x31 + x310*x55 + x313*x7 + x314*x33 + x316*x84 + x318*x96 + x320*x77 + x322*x91);
    double x328 = x302*x323 + x324*x325 + x326*x327;
    double x329 = x296*x324 + x298*x302 + x300*x326;
    double x330 = x195*x325 + x297*x323 + x299*x327;
    double x331 = x330 + 1.0;
    double x332 = x301*x328 - x329*x331;
    double x333 = x275[57];
    double x334 = x275[45];
    double x335 = x275[48];
    double x336 = x275[33];
    double x337 = x275[51];
    double x338 = x275[27];
    double x339 = x275[39];
    double x340 = x275[54];
    double x341 = x275[24];
    double x342 = x275[36];
    double x343 = x275[30];
    double x344 = x275[42];
    double x345 = x275[3];
    double x346 = x275[15];
    double x347 = x275[6];
    double x348 = x275[18];
    double x349 = x275[0];
    double x350 = x275[12];
    double x351 = x275[9];
    double x352 = x275[21];
    double x353 = x274*(x198*x333 + x202*x334 + x205*x335 + x208*x336 + x210*x337 + x212*x338 + x214*x339 + x216*x340 + x220*x341 + x222*x342 + x225*x343 + x227*x344 + x237*x345 + x243*x346 + x247*x347 + x250*x348 + x255*x349 + x257*x350 + x260*x351 + x262*x352);
    double x354 = x274*(x117*x337 + x122*x342 + x125*x335 + x128*x341 + x131*x333 + x133*x343 + x135*x344 + x137*x340 + x143*x336 + x145*x334 + x149*x338 + x151*x339 + x164*x351 + x169*x352 + x173*x347 + x177*x348 + x182*x349 + x184*x350 + x186*x345 + x188*x346);
    double x355 = x274*(x102*x349 + x104*x351 + x107*x345 + x109*x347 + x17*x343 + x20*x336 + x23*x341 + x26*x334 + x29*x342 + x31*x344 + x33*x339 + x333*x47 + x335*x43 + x337*x52 + x338*x7 + x340*x55 + x346*x91 + x348*x96 + x350*x77 + x352*x84);
    double x356 = x195*x353 + x297*x354 + x299*x355;
    double x357 = x302*x354 + x324*x353 + x326*x355;
    double x358 = x357 + 1.0;
    double x359 = -x301*x358 + x329*x356;
    double x360 = -x328*x356 + x331*x358;
    double x361 = S11*x332 + S12*x359 + S13*x360;
    double x362 = S12*x332 + S22*x359 + S23*x360;
    double x363 = S13*x332 + S23*x359 + S33*x360;
    double x364 = -x328*x362 - x329*x363 - x357*x361;
    double x365 = x111*x264 - x265*x269;
    double x366 = x182*x274;
    double x367 = -x111*x268 + x193*x265;
    double x368 = x255*x274;
    double x369 = -x193*x264 + x268*x269;
    double x370 = x102*x274;
    double x371 = x365*x366 + x367*x368 + x369*x370;
    double x372 = x323*x365 + x325*x367 + x327*x369;
    double x373 = x296*x367 + x298*x365 + x300*x369;
    double x374 = x373 + 1.0;
    double x375 = -x328*x374 + x329*x372;
    double x376 = x353*x367 + x354*x365 + x355*x369;
    double x377 = -x329*x376 + x358*x374;
    double x378 = x328*x376 - x358*x372;
    double x379 = S11*x375 + S12*x377 + S13*x378;
    double x380 = S12*x375 + S22*x377 + S23*x378;
    double x381 = S13*x375 + S23*x377 + S33*x378;
    double x382 = -x328*x380 - x329*x381 - x357*x379;
    double x383 = x195*x368 + x297*x366 + x299*x370;
    double x384 = x302*x366 + x324*x368 + x326*x370;
    double x385 = -x301*x372 + x331*x374;
    double x386 = x301*x376 - x356*x374;
    double x387 = -x331*x376 + x356*x372;
    double x388 = S11*x385 + S12*x386 + S13*x387;
    double x389 = S12*x385 + S22*x386 + S23*x387;
    double x390 = S13*x385 + S23*x386 + S33*x387;
    double x391 = 1.0*PENER + 1.0*SENER;
    double x392 = -x328*x389 - x329*x390 - x357*x388 + x391;
    double x393 = -x301*x363 - x330*x362 - x356*x361;
    double x394 = -x301*x390 - x330*x389 - x356*x388;
    double x395 = -x301*x381 - x330*x380 - x356*x379 + x391;
    double x396 = -x372*x380 - x373*x381 - x376*x379;
    double x397 = -x372*x389 - x373*x390 - x376*x388;
    double x398 = -x361*x376 - x362*x372 - x363*x373 + x391;
    double x399 = x186*x274;
    double x400 = x237*x274;
    double x401 = x107*x274;
    double x402 = x365*x399 + x367*x400 + x369*x401;
    double x403 = x195*x400 + x297*x399 + x299*x401;
    double x404 = x302*x399 + x324*x400 + x326*x401;
    double x405 = x173*x274;
    double x406 = x247*x274;
    double x407 = x109*x274;
    double x408 = x365*x405 + x367*x406 + x369*x407;
    double x409 = x195*x406 + x297*x405 + x299*x407;
    double x410 = x302*x405 + x324*x406 + x326*x407;
    double x411 = x164*x274;
    double x412 = x260*x274;
    double x413 = x104*x274;
    double x414 = x365*x411 + x367*x412 + x369*x413;
    double x415 = x195*x412 + x297*x411 + x299*x413;
    double x416 = x302*x411 + x324*x412 + x326*x413;
    double x417 = x184*x274;
    double x418 = x257*x274;
    double x419 = x274*x77;
    double x420 = x365*x417 + x367*x418 + x369*x419;
    double x421 = x195*x418 + x297*x417 + x299*x419;
    double x422 = x302*x417 + x324*x418 + x326*x419;
    double x423 = x188*x274;
    double x424 = x243*x274;
    double x425 = x274*x91;
    double x426 = x365*x423 + x367*x424 + x369*x425;
    double x427 = x195*x424 + x297*x423 + x299*x425;
    double x428 = x302*x423 + x324*x424 + x326*x425;
    double x429 = x177*x274;
    double x430 = x250*x274;
    double x431 = x274*x96;
    double x432 = x365*x429 + x367*x430 + x369*x431;
    double x433 = x195*x430 + x297*x429 + x299*x431;
    double x434 = x302*x429 + x324*x430 + x326*x431;
    double x435 = x169*x274;
    double x436 = x262*x274;
    double x437 = x274*x84;
    double x438 = x365*x435 + x367*x436 + x369*x437;
    double x439 = x195*x436 + x297*x435 + x299*x437;
    double x440 = x302*x435 + x324*x436 + x326*x437;
    double x441 = x128*x274;
    double x442 = x220*x274;
    double x443 = x23*x274;
    double x444 = x365*x441 + x367*x442 + x369*x443;
    double x445 = x195*x442 + x297*x441 + x299*x443;
    double x446 = x302*x441 + x324*x442 + x326*x443;
    double x447 = x149*x274;
    double x448 = x212*x274;
    double x449 = x274*x7;
    double x450 = x365*x447 + x367*x448 + x369*x449;
    double x451 = x195*x448 + x297*x447 + x299*x449;
    double x452 = x302*x447 + x324*x448 + x326*x449;
    double x453 = x133*x274;
    double x454 = x225*x274;
    double x455 = x17*x274;
    double x456 = x365*x453 + x367*x454 + x369*x455;
    double x457 = x195*x454 + x297*x453 + x299*x455;
    double x458 = x302*x453 + x324*x454 + x326*x455;
    double x459 = x143*x274;
    double x460 = x208*x274;
    double x461 = x20*x274;
    double x462 = x365*x459 + x367*x460 + x369*x461;
    double x463 = x195*x460 + x297*x459 + x299*x461;
    double x464 = x302*x459 + x324*x460 + x326*x461;
    double x465 = x122*x274;
    double x466 = x222*x274;
    double x467 = x274*x29;
    double x468 = x365*x465 + x367*x466 + x369*x467;
    double x469 = x195*x466 + x297*x465 + x299*x467;
    double x470 = x302*x465 + x324*x466 + x326*x467;
    double x471 = x151*x274;
    double x472 = x214*x274;
    double x473 = x274*x33;
    double x474 = x365*x471 + x367*x472 + x369*x473;
    double x475 = x195*x472 + x297*x471 + x299*x473;
    double x476 = x302*x471 + x324*x472 + x326*x473;
    double x477 = x135*x274;
    double x478 = x227*x274;
    double x479 = x274*x31;
    double x480 = x365*x477 + x367*x478 + x369*x479;
    double x481 = x195*x478 + x297*x477 + x299*x479;
    double x482 = x302*x477 + x324*x478 + x326*x479;
    double x483 = x145*x274;
    double x484 = x202*x274;
    double x485 = x26*x274;
    double x486 = x365*x483 + x367*x484 + x369*x485;
    double x487 = x195*x484 + x297*x483 + x299*x485;
    double x488 = x302*x483 + x324*x484 + x326*x485;
    double x489 = x125*x274;
    double x490 = x205*x274;
    double x491 = x274*x43;
    double x492 = x365*x489 + x367*x490 + x369*x491;
    double x493 = x195*x490 + x297*x489 + x299*x491;
    double x494 = x302*x489 + x324*x490 + x326*x491;
    double x495 = x117*x274;
    double x496 = x210*x274;
    double x497 = x274*x52;
    double x498 = x365*x495 + x367*x496 + x369*x497;
    double x499 = x195*x496 + x297*x495 + x299*x497;
    double x500 = x302*x495 + x324*x496 + x326*x497;
    double x501 = x137*x274;
    double x502 = x216*x274;
    double x503 = x274*x55;
    double x504 = x365*x501 + x367*x502 + x369*x503;
    double x505 = x195*x502 + x297*x501 + x299*x503;
    double x506 = x302*x501 + x324*x502 + x326*x503;
    double x507 = x131*x274;
    double x508 = x198*x274;
    double x509 = x274*x47;
    double x510 = x365*x507 + x367*x508 + x369*x509;
    double x511 = x195*x508 + x297*x507 + x299*x509;
    double x512 = x302*x507 + x324*x508 + x326*x509;
    
    res_0[0] = x273*(x364*x371 + x382*x383 + x384*x392);
    res_0[1] = x273*(x371*x393 + x383*x395 + x384*x394);
    res_0[2] = x273*(x371*x398 + x383*x396 + x384*x397);
    res_0[3] = x273*(x364*x402 + x382*x403 + x392*x404);
    res_0[4] = x273*(x393*x402 + x394*x404 + x395*x403);
    res_0[5] = x273*(x396*x403 + x397*x404 + x398*x402);
    res_0[6] = x273*(x364*x408 + x382*x409 + x392*x410);
    res_0[7] = x273*(x393*x408 + x394*x410 + x395*x409);
    res_0[8] = x273*(x396*x409 + x397*x410 + x398*x408);
    res_0[9] = x273*(x364*x414 + x382*x415 + x392*x416);
    res_0[10] = x273*(x393*x414 + x394*x416 + x395*x415);
    res_0[11] = x273*(x396*x415 + x397*x416 + x398*x414);
    res_0[12] = x273*(x364*x420 + x382*x421 + x392*x422);
    res_0[13] = x273*(x393*x420 + x394*x422 + x395*x421);
    res_0[14] = x273*(x396*x421 + x397*x422 + x398*x420);
    res_0[15] = x273*(x364*x426 + x382*x427 + x392*x428);
    res_0[16] = x273*(x393*x426 + x394*x428 + x395*x427);
    res_0[17] = x273*(x396*x427 + x397*x428 + x398*x426);
    res_0[18] = x273*(x364*x432 + x382*x433 + x392*x434);
    res_0[19] = x273*(x393*x432 + x394*x434 + x395*x433);
    res_0[20] = x273*(x396*x433 + x397*x434 + x398*x432);
    res_0[21] = x273*(x364*x438 + x382*x439 + x392*x440);
    res_0[22] = x273*(x393*x438 + x394*x440 + x395*x439);
    res_0[23] = x273*(x396*x439 + x397*x440 + x398*x438);
    res_0[24] = x273*(x364*x444 + x382*x445 + x392*x446);
    res_0[25] = x273*(x393*x444 + x394*x446 + x395*x445);
    res_0[26] = x273*(x396*x445 + x397*x446 + x398*x444);
    res_0[27] = x273*(x364*x450 + x382*x451 + x392*x452);
    res_0[28] = x273*(x393*x450 + x394*x452 + x395*x451);
    res_0[29] = x273*(x396*x451 + x397*x452 + x398*x450);
    res_0[30] = x273*(x364*x456 + x382*x457 + x392*x458);
    res_0[31] = x273*(x393*x456 + x394*x458 + x395*x457);
    res_0[32] = x273*(x396*x457 + x397*x458 + x398*x456);
    res_0[33] = x273*(x364*x462 + x382*x463 + x392*x464);
    res_0[34] = x273*(x393*x462 + x394*x464 + x395*x463);
    res_0[35] = x273*(x396*x463 + x397*x464 + x398*x462);
    res_0[36] = x273*(x364*x468 + x382*x469 + x392*x470);
    res_0[37] = x273*(x393*x468 + x394*x470 + x395*x469);
    res_0[38] = x273*(x396*x469 + x397*x470 + x398*x468);
    res_0[39] = x273*(x364*x474 + x382*x475 + x392*x476);
    res_0[40] = x273*(x393*x474 + x394*x476 + x395*x475);
    res_0[41] = x273*(x396*x475 + x397*x476 + x398*x474);
    res_0[42] = x273*(x364*x480 + x382*x481 + x392*x482);
    res_0[43] = x273*(x393*x480 + x394*x482 + x395*x481);
    res_0[44] = x273*(x396*x481 + x397*x482 + x398*x480);
    res_0[45] = x273*(x364*x486 + x382*x487 + x392*x488);
    res_0[46] = x273*(x393*x486 + x394*x488 + x395*x487);
    res_0[47] = x273*(x396*x487 + x397*x488 + x398*x486);
    res_0[48] = x273*(x364*x492 + x382*x493 + x392*x494);
    res_0[49] = x273*(x393*x492 + x394*x494 + x395*x493);
    res_0[50] = x273*(x396*x493 + x397*x494 + x398*x492);
    res_0[51] = x273*(x364*x498 + x382*x499 + x392*x500);
    res_0[52] = x273*(x393*x498 + x394*x500 + x395*x499);
    res_0[53] = x273*(x396*x499 + x397*x500 + x398*x498);
    res_0[54] = x273*(x364*x504 + x382*x505 + x392*x506);
    res_0[55] = x273*(x393*x504 + x394*x506 + x395*x505);
    res_0[56] = x273*(x396*x505 + x397*x506 + x398*x504);
    res_0[57] = x273*(x364*x510 + x382*x511 + x392*x512);
    res_0[58] = x273*(x393*x510 + x394*x512 + x395*x511);
    res_0[59] = x273*(x396*x511 + x397*x512 + x398*x510);
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
    
    double *x0 = coord;
    double x1 = -x0[2];
    double x2 = x1 + x0[11];
    double x3 = -x0[0];
    double x4 = x3 + x0[3];
    double x5 = -x0[1];
    double x6 = x5 + x0[7];
    double x7 = x4*x6;
    double x8 = x1 + x0[5];
    double x9 = x3 + x0[6];
    double x10 = x5 + x0[10];
    double x11 = x10*x9;
    double x12 = x1 + x0[8];
    double x13 = x3 + x0[9];
    double x14 = x5 + x0[4];
    double x15 = x13*x14;
    double x16 = x10*x4;
    double x17 = x14*x9;
    double x18 = x13*x6;
    double x19 = x11*x8 + x12*x15 - x12*x16 - x17*x2 - x18*x8 + x2*x7;
    double x20 = 1.0/x19;
    double x21 = x20*(-x17 + x7);
    double x22 = x20*(x11 - x18);
    double x23 = x20*(x15 - x16);
    double x24 = -x21 - x22 - x23;
    double *x25 = U;
    double x26 = -x25[1];
    double x27 = x26 + x25[4];
    double x28 = x20*(-x10*x12 + x2*x6);
    double x29 = x26 + x25[7];
    double x30 = x20*(x10*x8 - x14*x2);
    double x31 = x26 + x25[10];
    double x32 = x20*(x12*x14 - x6*x8);
    double x33 = x27*x28 + x29*x30 + x31*x32;
    double x34 = -x25[2];
    double x35 = x34 + x25[5];
    double x36 = x20*(x12*x13 - x2*x9);
    double x37 = x34 + x25[8];
    double x38 = x20*(-x13*x8 + x2*x4);
    double x39 = x34 + x25[11];
    double x40 = x20*(-x12*x4 + x8*x9);
    double x41 = x35*x36 + x37*x38 + x39*x40;
    double x42 = x28*x35 + x30*x37 + x32*x39;
    double x43 = x27*x36 + x29*x38 + x31*x40;
    double x44 = x43 + 1.0;
    double x45 = x33*x41 - x42*x44;
    double x46 = -x25[0];
    double x47 = x46 + x25[3];
    double x48 = x46 + x25[6];
    double x49 = x46 + x25[9];
    double x50 = x36*x47 + x38*x48 + x40*x49;
    double x51 = x28*x47 + x30*x48 + x32*x49;
    double x52 = x51 + 1.0;
    double x53 = -x41*x52 + x42*x50;
    double x54 = -x33*x50 + x44*x52;
    double x55 = S11*x45 + S12*x53 + S13*x54;
    double x56 = S12*x45 + S22*x53 + S23*x54;
    double x57 = S13*x45 + S23*x53 + S33*x54;
    double x58 = -x33*x56 - x42*x57 - x51*x55;
    double x59 = -x36 - x38 - x40;
    double x60 = x21*x31 + x22*x27 + x23*x29;
    double x61 = x21*x39 + x22*x35 + x23*x37;
    double x62 = x61 + 1.0;
    double x63 = -x33*x62 + x42*x60;
    double x64 = x21*x49 + x22*x47 + x23*x48;
    double x65 = -x42*x64 + x52*x62;
    double x66 = x33*x64 - x52*x60;
    double x67 = S11*x63 + S12*x65 + S13*x66;
    double x68 = S12*x63 + S22*x65 + S23*x66;
    double x69 = S13*x63 + S23*x65 + S33*x66;
    double x70 = -x33*x68 - x42*x69 - x51*x67;
    double x71 = -x28 - x30 - x32;
    double x72 = -x41*x60 + x44*x62;
    double x73 = x41*x64 - x50*x62;
    double x74 = -x44*x64 + x50*x60;
    double x75 = S11*x72 + S12*x73 + S13*x74;
    double x76 = S12*x72 + S22*x73 + S23*x74;
    double x77 = S13*x72 + S23*x73 + S33*x74;
    double *x78 = V;
    double x79 = x78[3];
    double x80 = x78[6];
    double x81 = x78[9];
    double x82 = -r - s - t + 1;
    double x83 = x78[0];
    double x84 = r*x79 + s*x80 + t*x81 + x82*x83;
    double x85 = x78[4];
    double x86 = x78[7];
    double x87 = x78[10];
    double x88 = x78[1];
    double x89 = r*x85 + s*x86 + t*x87 + x82*x88;
    double x90 = x78[5];
    double x91 = x78[8];
    double x92 = x78[11];
    double x93 = x78[2];
    double x94 = r*x90 + s*x91 + t*x92 + x82*x93;
    double x95 = 1.0*PENER + 1.0*SENER - 0.5*rho*(x84*x84 + x89*x89 + x94*x94);
    double x96 = -x33*x76 - x42*x77 - x51*x75 + x95;
    double *x97 = A;
    double x98 = r*x97[3] + s*x97[6] + t*x97[9] + x82*x97[0];
    double x99 = r*x97[4] + s*x97[7] + t*x97[10] + x82*x97[1];
    double x100 = r*x97[5] + s*x97[8] + t*x97[11] + x82*x97[2];
    double x101 = -x83;
    double x102 = x101 + x79;
    double x103 = x101 + x80;
    double x104 = x101 + x81;
    double x105 = -x88;
    double x106 = x105 + x85;
    double x107 = x105 + x86;
    double x108 = x105 + x87;
    double x109 = -x93;
    double x110 = x109 + x90;
    double x111 = x109 + x91;
    double x112 = x109 + x92;
    double x113 = x100*x42 + x33*x99 + x51*x98 + x84*(x102*x28 + x103*x30 + x104*x32) + x89*(x106*x28 + x107*x30 + x108*x32) + x94*(x110*x28 + x111*x30 + x112*x32);
    double x114 = rho*x82;
    double x115 = -x41*x57 - x43*x56 - x50*x55;
    double x116 = -x41*x69 - x43*x68 - x50*x67 + x95;
    double x117 = -x41*x77 - x43*x76 - x50*x75;
    double x118 = x100*x41 + x43*x99 + x50*x98 + x84*(x102*x36 + x103*x38 + x104*x40) + x89*(x106*x36 + x107*x38 + x108*x40) + x94*(x110*x36 + x111*x38 + x112*x40);
    double x119 = -x55*x64 - x56*x60 - x57*x61 + x95;
    double x120 = -x60*x68 - x61*x69 - x64*x67;
    double x121 = -x60*x76 - x61*x77 - x64*x75;
    double x122 = x100*x61 + x60*x99 + x64*x98 + x84*(x102*x22 + x103*x23 + x104*x21) + x89*(x106*x22 + x107*x23 + x108*x21) + x94*(x110*x22 + x111*x23 + x112*x21);
    double x123 = r*rho;
    double x124 = rho*s;
    double x125 = rho*t;
    
    res_0[0] = x19*(-x113*x114 + x24*x58 + x59*x70 + x71*x96);
    res_0[1] = x19*(-x114*x118 + x115*x24 + x116*x59 + x117*x71);
    res_0[2] = x19*(-x114*x122 + x119*x24 + x120*x59 + x121*x71);
    res_0[3] = x19*(-x113*x123 + x22*x58 + x28*x96 + x36*x70);
    res_0[4] = x19*(x115*x22 + x116*x36 + x117*x28 - x118*x123);
    res_0[5] = x19*(x119*x22 + x120*x36 + x121*x28 - x122*x123);
    res_0[6] = x19*(-x113*x124 + x23*x58 + x30*x96 + x38*x70);
    res_0[7] = x19*(x115*x23 + x116*x38 + x117*x30 - x118*x124);
    res_0[8] = x19*(x119*x23 + x120*x38 + x121*x30 - x122*x124);
    res_0[9] = x19*(-x113*x125 + x21*x58 + x32*x96 + x40*x70);
    res_0[10] = x19*(x115*x21 + x116*x40 + x117*x32 - x118*x125);
    res_0[11] = x19*(x119*x21 + x120*x40 + x121*x32 - x122*x125);
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
    
    double *x0 = coord;
    double x1 = -x0[2];
    double x2 = x1 + x0[11];
    double x3 = -x0[0];
    double x4 = x3 + x0[3];
    double x5 = -x0[1];
    double x6 = x5 + x0[7];
    double x7 = x4*x6;
    double x8 = x1 + x0[5];
    double x9 = x3 + x0[6];
    double x10 = x5 + x0[10];
    double x11 = x10*x9;
    double x12 = x1 + x0[8];
    double x13 = x3 + x0[9];
    double x14 = x5 + x0[4];
    double x15 = x13*x14;
    double x16 = x10*x4;
    double x17 = x14*x9;
    double x18 = x13*x6;
    double x19 = x11*x8 + x12*x15 - x12*x16 - x17*x2 - x18*x8 + x2*x7;
    double x20 = 1.0/x19;
    double x21 = x20*(-x17 + x7);
    double x22 = x20*(x11 - x18);
    double x23 = x20*(x15 - x16);
    double x24 = -x21 - x22 - x23;
    double *x25 = U;
    double x26 = -x25[1];
    double x27 = x26 + x25[4];
    double x28 = x20*(-x10*x12 + x2*x6);
    double x29 = x26 + x25[7];
    double x30 = x20*(x10*x8 - x14*x2);
    double x31 = x26 + x25[10];
    double x32 = x20*(x12*x14 - x6*x8);
    double x33 = x27*x28 + x29*x30 + x31*x32;
    double x34 = -x25[2];
    double x35 = x34 + x25[5];
    double x36 = x20*(x12*x13 - x2*x9);
    double x37 = x34 + x25[8];
    double x38 = x20*(-x13*x8 + x2*x4);
    double x39 = x34 + x25[11];
    double x40 = x20*(-x12*x4 + x8*x9);
    double x41 = x35*x36 + x37*x38 + x39*x40;
    double x42 = x28*x35 + x30*x37 + x32*x39;
    double x43 = x27*x36 + x29*x38 + x31*x40 + 1.0;
    double x44 = x33*x41 - x42*x43;
    double x45 = -x25[0];
    double x46 = x45 + x25[3];
    double x47 = x45 + x25[6];
    double x48 = x45 + x25[9];
    double x49 = x36*x46 + x38*x47 + x40*x48;
    double x50 = x28*x46 + x30*x47 + x32*x48 + 1.0;
    double x51 = -x41*x50 + x42*x49;
    double x52 = -x33*x49 + x43*x50;
    double x53 = S11*x44 + S12*x51 + S13*x52;
    double x54 = S12*x44 + S22*x51 + S23*x52;
    double x55 = S13*x44 + S23*x51 + S33*x52;
    double x56 = -x33*x54 - x42*x55 - x50*x53;
    double x57 = -x36 - x38 - x40;
    double x58 = x21*x31 + x22*x27 + x23*x29;
    double x59 = x21*x39 + x22*x35 + x23*x37 + 1.0;
    double x60 = -x33*x59 + x42*x58;
    double x61 = x21*x48 + x22*x46 + x23*x47;
    double x62 = -x42*x61 + x50*x59;
    double x63 = x33*x61 - x50*x58;
    double x64 = S11*x60 + S12*x62 + S13*x63;
    double x65 = S12*x60 + S22*x62 + S23*x63;
    double x66 = S13*x60 + S23*x62 + S33*x63;
    double x67 = -x33*x65 - x42*x66 - x50*x64;
    double x68 = -x28 - x30 - x32;
    double x69 = -x41*x58 + x43*x59;
    double x70 = x41*x61 - x49*x59;
    double x71 = -x43*x61 + x49*x58;
    double x72 = S11*x69 + S12*x70 + S13*x71;
    double x73 = S12*x69 + S22*x70 + S23*x71;
    double x74 = S13*x69 + S23*x70 + S33*x71;
    double x75 = 1.0*PENER + 1.0*SENER;
    double x76 = -x33*x73 - x42*x74 - x50*x72 + x75;
    double x77 = -x41*x55 - x43*x54 - x49*x53;
    double x78 = -x41*x66 - x43*x65 - x49*x64 + x75;
    double x79 = -x41*x74 - x43*x73 - x49*x72;
    double x80 = -x53*x61 - x54*x58 - x55*x59 + x75;
    double x81 = -x58*x65 - x59*x66 - x61*x64;
    double x82 = -x58*x73 - x59*x74 - x61*x72;
    
    res_0[0] = x19*(x24*x56 + x57*x67 + x68*x76);
    res_0[1] = x19*(x24*x77 + x57*x78 + x68*x79);
    res_0[2] = x19*(x24*x80 + x57*x81 + x68*x82);
    res_0[3] = x19*(x22*x56 + x28*x76 + x36*x67);
    res_0[4] = x19*(x22*x77 + x28*x79 + x36*x78);
    res_0[5] = x19*(x22*x80 + x28*x82 + x36*x81);
    res_0[6] = x19*(x23*x56 + x30*x76 + x38*x67);
    res_0[7] = x19*(x23*x77 + x30*x79 + x38*x78);
    res_0[8] = x19*(x23*x80 + x30*x82 + x38*x81);
    res_0[9] = x19*(x21*x56 + x32*x76 + x40*x67);
    res_0[10] = x19*(x21*x77 + x32*x79 + x40*x78);
    res_0[11] = x19*(x21*x80 + x32*x82 + x40*x81);
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
    
    double *x0 = coord;
    double x1 = -x0[2];
    double x2 = x1 + x0[11];
    double x3 = -x0[0];
    double x4 = x3 + x0[3];
    double x5 = -x0[1];
    double x6 = x5 + x0[7];
    double x7 = x4*x6;
    double x8 = x1 + x0[5];
    double x9 = x3 + x0[6];
    double x10 = x5 + x0[10];
    double x11 = x10*x9;
    double x12 = x1 + x0[8];
    double x13 = x3 + x0[9];
    double x14 = x5 + x0[4];
    double x15 = x13*x14;
    double x16 = x10*x4;
    double x17 = x14*x9;
    double x18 = x13*x6;
    double x19 = x11*x8 + x12*x15 - x12*x16 - x17*x2 - x18*x8 + x2*x7;
    double x20 = 1.0/x19;
    double x21 = x20*(-x17 + x7);
    double x22 = x20*(x11 - x18);
    double x23 = x20*(x15 - x16);
    double x24 = -x21 - x22 - x23;
    double *x25 = U;
    double x26 = -x25[1];
    double x27 = x26 + x25[4];
    double x28 = x20*(-x10*x12 + x2*x6);
    double x29 = x26 + x25[7];
    double x30 = x20*(x10*x8 - x14*x2);
    double x31 = x26 + x25[10];
    double x32 = x20*(x12*x14 - x6*x8);
    double x33 = x27*x28 + x29*x30 + x31*x32;
    double x34 = -x25[2];
    double x35 = x34 + x25[5];
    double x36 = x20*(x12*x13 - x2*x9);
    double x37 = x34 + x25[8];
    double x38 = x20*(-x13*x8 + x2*x4);
    double x39 = x34 + x25[11];
    double x40 = x20*(-x12*x4 + x8*x9);
    double x41 = x35*x36 + x37*x38 + x39*x40;
    double x42 = x28*x35 + x30*x37 + x32*x39;
    double x43 = x27*x36 + x29*x38 + x31*x40;
    double x44 = x43 + 1.0;
    double x45 = x33*x41 - x42*x44;
    double x46 = -x25[0];
    double x47 = x46 + x25[3];
    double x48 = x46 + x25[6];
    double x49 = x46 + x25[9];
    double x50 = x36*x47 + x38*x48 + x40*x49;
    double x51 = x28*x47 + x30*x48 + x32*x49;
    double x52 = x51 + 1.0;
    double x53 = -x41*x52 + x42*x50;
    double x54 = -x33*x50 + x44*x52;
    double x55 = S11*x45 + S12*x53 + S13*x54;
    double x56 = S12*x45 + S22*x53 + S23*x54;
    double x57 = S13*x45 + S23*x53 + S33*x54;
    double x58 = -x33*x56 - x42*x57 - x51*x55;
    double x59 = -x36 - x38 - x40;
    double x60 = x21*x31 + x22*x27 + x23*x29;
    double x61 = x21*x39 + x22*x35 + x23*x37;
    double x62 = x61 + 1.0;
    double x63 = -x33*x62 + x42*x60;
    double x64 = x21*x49 + x22*x47 + x23*x48;
    double x65 = -x42*x64 + x52*x62;
    double x66 = x33*x64 - x52*x60;
    double x67 = S11*x63 + S12*x65 + S13*x66;
    double x68 = S12*x63 + S22*x65 + S23*x66;
    double x69 = S13*x63 + S23*x65 + S33*x66;
    double x70 = -x33*x68 - x42*x69 - x51*x67;
    double x71 = -x28 - x30 - x32;
    double x72 = -x41*x60 + x44*x62;
    double x73 = x41*x64 - x50*x62;
    double x74 = -x44*x64 + x50*x60;
    double x75 = S11*x72 + S12*x73 + S13*x74;
    double x76 = S12*x72 + S22*x73 + S23*x74;
    double x77 = S13*x72 + S23*x73 + S33*x74;
    double x78 = 1.0*PENER + 1.0*SENER;
    double x79 = -x33*x76 - x42*x77 - x51*x75 + x78;
    double x80 = -x41*x57 - x43*x56 - x50*x55;
    double x81 = -x41*x69 - x43*x68 - x50*x67 + x78;
    double x82 = -x41*x77 - x43*x76 - x50*x75;
    double x83 = -x55*x64 - x56*x60 - x57*x61 + x78;
    double x84 = -x60*x68 - x61*x69 - x64*x67;
    double x85 = -x60*x76 - x61*x77 - x64*x75;
    
    res_0[0] = x19*(x24*x58 + x59*x70 + x71*x79);
    res_0[1] = x19*(x24*x80 + x59*x81 + x71*x82);
    res_0[2] = x19*(x24*x83 + x59*x84 + x71*x85);
    res_0[3] = x19*(x22*x58 + x28*x79 + x36*x70);
    res_0[4] = x19*(x22*x80 + x28*x82 + x36*x81);
    res_0[5] = x19*(x22*x83 + x28*x85 + x36*x84);
    res_0[6] = x19*(x23*x58 + x30*x79 + x38*x70);
    res_0[7] = x19*(x23*x80 + x30*x82 + x38*x81);
    res_0[8] = x19*(x23*x83 + x30*x85 + x38*x84);
    res_0[9] = x19*(x21*x58 + x32*x79 + x40*x70);
    res_0[10] = x19*(x21*x80 + x32*x82 + x40*x81);
    res_0[11] = x19*(x21*x83 + x32*x85 + x40*x84);
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
    double *x2 = coord;
    double x3 = 4*s;
    double x4 = -x3;
    double x5 = 4*t;
    double x6 = 4 - x5;
    double x7 = -8*r + x4 + x6;
    double x8 = x2[12];
    double x9 = x2[21];
    double x10 = -x5*x9;
    double x11 = x2[15];
    double x12 = x2[24];
    double x13 = x0 + x3 + x5 - 3;
    double x14 = x13*x2[0];
    double x15 = x2[18];
    double x16 = x14 - x15*x3;
    double x17 = x1*x2[3] + x10 + x11*x3 + x12*x5 + x16 + x7*x8;
    double x18 = x3 - 1;
    double x19 = -x0;
    double x20 = -8*s + x19 + x6;
    double x21 = x2[19];
    double x22 = x2[13];
    double x23 = -x0*x22;
    double x24 = x2[16];
    double x25 = x2[28];
    double x26 = x13*x2[1];
    double x27 = x2[22];
    double x28 = x26 - x27*x5;
    double x29 = x0*x24 + x18*x2[7] + x20*x21 + x23 + x25*x5 + x28;
    double x30 = x17*x29;
    double x31 = -x21*x3;
    double x32 = x2[25];
    double x33 = x1*x2[4] + x22*x7 + x24*x3 + x28 + x31 + x32*x5;
    double x34 = -x0*x8;
    double x35 = x2[27];
    double x36 = x0*x11 + x10 + x14 + x15*x20 + x18*x2[6] + x34 + x35*x5;
    double x37 = x33*x36;
    double x38 = x30 - x37;
    double x39 = x5 - 1;
    double x40 = -8*t + x19 + x4 + 4;
    double x41 = x2[23];
    double x42 = x2[20];
    double x43 = -x3*x42;
    double x44 = x2[26];
    double x45 = x2[29];
    double x46 = x13*x2[2];
    double x47 = x2[14];
    double x48 = -x0*x47 + x46;
    double x49 = x0*x44 + x3*x45 + x39*x2[11] + x40*x41 + x43 + x48;
    double x50 = -x41*x5;
    double x51 = x2[17];
    double x52 = x0*x51 + x18*x2[8] + x20*x42 + x45*x5 + x48 + x50;
    double x53 = x0*x12 + x16 + x3*x35 + x34 + x39*x2[9] + x40*x9;
    double x54 = x33*x53;
    double x55 = x1*x2[5] + x3*x51 + x43 + x44*x5 + x46 + x47*x7 + x50;
    double x56 = x0*x32 + x23 + x25*x3 + x26 + x27*x40 + x31 + x39*x2[10];
    double x57 = x36*x56;
    double x58 = x17*x56;
    double x59 = x29*x53;
    double x60 = x30*x49 - x37*x49 + x52*x54 - x52*x58 + x55*x57 - x55*x59;
    double x61 = 1.0/x60;
    double *x62 = U;
    double x63 = x62[22];
    double x64 = x62[19];
    double x65 = -x3*x64;
    double x66 = x62[25];
    double x67 = x62[28];
    double x68 = x13*x62[1];
    double x69 = x62[13];
    double x70 = -x0*x69 + x68;
    double x71 = x61*(x0*x66 + x3*x67 + x39*x62[10] + x40*x63 + x65 + x70);
    double x72 = x54 - x58;
    double x73 = -x5*x63;
    double x74 = x62[16];
    double x75 = x61*(x0*x74 + x18*x62[7] + x20*x64 + x5*x67 + x70 + x73);
    double x76 = x57 - x59;
    double x77 = x61*(x1*x62[4] + x3*x74 + x5*x66 + x65 + x68 + x69*x7 + x73);
    double x78 = x38*x71 + x72*x75 + x76*x77;
    double x79 = -x29*x55 + x33*x52;
    double x80 = x62[23];
    double x81 = x62[20];
    double x82 = -x3*x81;
    double x83 = x62[26];
    double x84 = x62[29];
    double x85 = x13*x62[2];
    double x86 = x62[14];
    double x87 = -x0*x86 + x85;
    double x88 = x61*(x0*x83 + x3*x84 + x39*x62[11] + x40*x80 + x82 + x87);
    double x89 = -x33*x49 + x55*x56;
    double x90 = -x5*x80;
    double x91 = x62[17];
    double x92 = x61*(x0*x91 + x18*x62[8] + x20*x81 + x5*x84 + x87 + x90);
    double x93 = x29*x49 - x52*x56;
    double x94 = x61*(x1*x62[5] + x3*x91 + x5*x83 + x7*x86 + x82 + x85 + x90);
    double x95 = x79*x88 + x89*x92 + x93*x94;
    double x96 = x71*x79 + x75*x89 + x77*x93;
    double x97 = x38*x88 + x72*x92 + x76*x94;
    double x98 = x97 + 1.0;
    double x99 = x78*x95 - x96*x98;
    double x100 = x62[21];
    double x101 = x62[18];
    double x102 = -x101*x3;
    double x103 = x62[24];
    double x104 = x62[27];
    double x105 = x13*x62[0];
    double x106 = x62[12];
    double x107 = -x0*x106 + x105;
    double x108 = x61*(x0*x103 + x100*x40 + x102 + x104*x3 + x107 + x39*x62[9]);
    double x109 = -x100*x5;
    double x110 = x62[15];
    double x111 = x61*(x0*x110 + x101*x20 + x104*x5 + x107 + x109 + x18*x62[6]);
    double x112 = x61*(x1*x62[3] + x102 + x103*x5 + x105 + x106*x7 + x109 + x110*x3);
    double x113 = x108*x79 + x111*x89 + x112*x93;
    double x114 = x113 + 1.0;
    double x115 = x108*x38 + x111*x72 + x112*x76;
    double x116 = x114*x98 - x115*x95;
    double x117 = -x114*x78 + x115*x96;
    double x118 = S11*x99 + S12*x116 + S13*x117;
    double x119 = S12*x99 + S22*x116 + S23*x117;
    double x120 = S13*x99 + S23*x116 + S33*x117;
    double x121 = -x113*x118 - x119*x96 - x120*x95;
    double x122 = x17*x49 - x53*x55;
    double x123 = x13*x61;
    double x124 = -x17*x52 + x36*x55;
    double x125 = -x36*x49 + x52*x53;
    double x126 = x122*x123 + x123*x124 + x123*x125;
    double x127 = x122*x92 + x124*x88 + x125*x94;
    double x128 = x122*x75 + x124*x71 + x125*x77;
    double x129 = x128 + 1.0;
    double x130 = x127*x96 - x129*x95;
    double x131 = x108*x124 + x111*x122 + x112*x125;
    double x132 = -x114*x127 + x131*x95;
    double x133 = x114*x129 - x131*x96;
    double x134 = S11*x130 + S12*x132 + S13*x133;
    double x135 = S12*x130 + S22*x132 + S23*x133;
    double x136 = S13*x130 + S23*x132 + S33*x133;
    double x137 = -x113*x134 - x135*x96 - x136*x95;
    double x138 = x123*x38 + x123*x72 + x123*x76;
    double x139 = x123*x79 + x123*x89 + x123*x93;
    double x140 = -x127*x78 + x129*x98;
    double x141 = x115*x127 - x131*x98;
    double x142 = -x115*x129 + x131*x78;
    double x143 = S11*x140 + S12*x141 + S13*x142;
    double x144 = S12*x140 + S22*x141 + S23*x142;
    double x145 = S13*x140 + S23*x141 + S33*x142;
    double x146 = r*r;
    double x147 = 2*x146;
    double x148 = -r + x147;
    double *x149 = V;
    double x150 = x149[3];
    double x151 = s*s;
    double x152 = 2*x151;
    double x153 = -s + x152;
    double x154 = x149[6];
    double x155 = t*t;
    double x156 = 2*x155;
    double x157 = -t + x156;
    double x158 = x149[9];
    double x159 = s*x0;
    double x160 = -x159;
    double x161 = t*x0;
    double x162 = -x161;
    double x163 = x0 - 4*x146 + x160 + x162;
    double x164 = x149[12];
    double x165 = t*x3;
    double x166 = -x165;
    double x167 = -4*x151 + x160 + x166 + x3;
    double x168 = x149[18];
    double x169 = -4*x155 + x162 + x166 + x5;
    double x170 = x149[21];
    double x171 = -3*r - 3*s - 3*t + x147 + x152 + x156 + x159 + x161 + x165 + 1;
    double x172 = x149[0];
    double x173 = x149[15];
    double x174 = x149[24];
    double x175 = x149[27];
    double x176 = x148*x150 + x153*x154 + x157*x158 + x159*x173 + x161*x174 + x163*x164 + x165*x175 + x167*x168 + x169*x170 + x171*x172;
    double x177 = x149[4];
    double x178 = x149[7];
    double x179 = x149[10];
    double x180 = x149[13];
    double x181 = x149[19];
    double x182 = x149[22];
    double x183 = x149[1];
    double x184 = x149[16];
    double x185 = x149[25];
    double x186 = x149[28];
    double x187 = x148*x177 + x153*x178 + x157*x179 + x159*x184 + x161*x185 + x163*x180 + x165*x186 + x167*x181 + x169*x182 + x171*x183;
    double x188 = x149[5];
    double x189 = x149[8];
    double x190 = x149[11];
    double x191 = x149[14];
    double x192 = x149[20];
    double x193 = x149[23];
    double x194 = x149[2];
    double x195 = x149[17];
    double x196 = x149[26];
    double x197 = x149[29];
    double x198 = x148*x188 + x153*x189 + x157*x190 + x159*x195 + x161*x196 + x163*x191 + x165*x197 + x167*x192 + x169*x193 + x171*x194;
    double x199 = 1.0*PENER + 1.0*SENER - 0.5*rho*(x176*x176 + x187*x187 + x198*x198);
    double x200 = -x113*x143 - x144*x96 - x145*x95 + x199;
    double *x201 = A;
    double x202 = x148*x201[3] + x153*x201[6] + x157*x201[9] + x159*x201[15] + x161*x201[24] + x163*x201[12] + x165*x201[27] + x167*x201[18] + x169*x201[21] + x171*x201[0];
    double x203 = x148*x201[4] + x153*x201[7] + x157*x201[10] + x159*x201[16] + x161*x201[25] + x163*x201[13] + x165*x201[28] + x167*x201[19] + x169*x201[22] + x171*x201[1];
    double x204 = x148*x201[5] + x153*x201[8] + x157*x201[11] + x159*x201[17] + x161*x201[26] + x163*x201[14] + x165*x201[29] + x167*x201[20] + x169*x201[23] + x171*x201[2];
    double x205 = -x168*x3;
    double x206 = x13*x172;
    double x207 = -x0*x164 + x206;
    double x208 = x0*x174 + x158*x39 + x170*x40 + x175*x3 + x205 + x207;
    double x209 = x61*x79;
    double x210 = -x170*x5;
    double x211 = x0*x173 + x154*x18 + x168*x20 + x175*x5 + x207 + x210;
    double x212 = x61*x89;
    double x213 = x1*x150 + x164*x7 + x173*x3 + x174*x5 + x205 + x206 + x210;
    double x214 = x61*x93;
    double x215 = -x181*x3;
    double x216 = x13*x183;
    double x217 = -x0*x180 + x216;
    double x218 = x0*x185 + x179*x39 + x182*x40 + x186*x3 + x215 + x217;
    double x219 = -x182*x5;
    double x220 = x0*x184 + x178*x18 + x181*x20 + x186*x5 + x217 + x219;
    double x221 = x1*x177 + x180*x7 + x184*x3 + x185*x5 + x215 + x216 + x219;
    double x222 = -x192*x3;
    double x223 = x13*x194;
    double x224 = -x0*x191 + x223;
    double x225 = x0*x196 + x190*x39 + x193*x40 + x197*x3 + x222 + x224;
    double x226 = -x193*x5;
    double x227 = x0*x195 + x18*x189 + x192*x20 + x197*x5 + x224 + x226;
    double x228 = x1*x188 + x191*x7 + x195*x3 + x196*x5 + x222 + x223 + x226;
    double x229 = x113*x202 + x176*(x208*x209 + x211*x212 + x213*x214) + x187*(x209*x218 + x212*x220 + x214*x221) + x198*(x209*x225 + x212*x227 + x214*x228) + x203*x96 + x204*x95;
    double x230 = rho*x171;
    double x231 = -x127*x136 - x128*x135 - x131*x134;
    double x232 = -x127*x145 - x128*x144 - x131*x143;
    double x233 = -x118*x131 - x119*x128 - x120*x127 + x199;
    double x234 = x122*x61;
    double x235 = x124*x61;
    double x236 = x125*x61;
    double x237 = x127*x204 + x128*x203 + x131*x202 + x176*(x208*x235 + x211*x234 + x213*x236) + x187*(x218*x235 + x220*x234 + x221*x236) + x198*(x225*x235 + x227*x234 + x228*x236);
    double x238 = -x115*x118 - x119*x78 - x120*x97;
    double x239 = -x115*x143 - x144*x78 - x145*x97;
    double x240 = -x115*x134 - x135*x78 - x136*x97 + x199;
    double x241 = x38*x61;
    double x242 = x61*x72;
    double x243 = x61*x76;
    double x244 = x115*x202 + x176*(x208*x241 + x211*x242 + x213*x243) + x187*(x218*x241 + x220*x242 + x221*x243) + x198*(x225*x241 + x227*x242 + x228*x243) + x203*x78 + x204*x97;
    double x245 = rho*x148;
    double x246 = x1*x243;
    double x247 = x1*x214;
    double x248 = x1*x236;
    double x249 = rho*x153;
    double x250 = x18*x234;
    double x251 = x18*x242;
    double x252 = x18*x212;
    double x253 = rho*x157;
    double x254 = x241*x39;
    double x255 = x209*x39;
    double x256 = x235*x39;
    double x257 = x0*x234;
    double x258 = x0*x235;
    double x259 = x236*x7 - x257 - x258;
    double x260 = x0*x241;
    double x261 = x0*x242;
    double x262 = x243*x7 - x260 - x261;
    double x263 = x0*x209;
    double x264 = x0*x212;
    double x265 = x214*x7 - x263 - x264;
    double x266 = rho*x163;
    double x267 = x236*x3;
    double x268 = x257 + x267;
    double x269 = x243*x3;
    double x270 = x261 + x269;
    double x271 = x214*x3;
    double x272 = x264 + x271;
    double x273 = rho*x159;
    double x274 = x235*x3;
    double x275 = x20*x234 - x267 - x274;
    double x276 = x241*x3;
    double x277 = x20*x242 - x269 - x276;
    double x278 = x209*x3;
    double x279 = x20*x212 - x271 - x278;
    double x280 = rho*x167;
    double x281 = x234*x5;
    double x282 = x236*x5;
    double x283 = x235*x40 - x281 - x282;
    double x284 = x242*x5;
    double x285 = x243*x5;
    double x286 = x241*x40 - x284 - x285;
    double x287 = x212*x5;
    double x288 = x214*x5;
    double x289 = x209*x40 - x287 - x288;
    double x290 = rho*x169;
    double x291 = x260 + x285;
    double x292 = x263 + x288;
    double x293 = x258 + x282;
    double x294 = rho*x161;
    double x295 = x276 + x284;
    double x296 = x278 + x287;
    double x297 = x274 + x281;
    double x298 = rho*x165;
    
    res_0[0] = x60*(x121*x126 + x137*x138 + x139*x200 - x229*x230);
    res_0[1] = x60*(x126*x233 + x138*x231 + x139*x232 - x230*x237);
    res_0[2] = x60*(x126*x238 + x138*x240 + x139*x239 - x230*x244);
    res_0[3] = x60*(x121*x248 + x137*x246 + x200*x247 - x229*x245);
    res_0[4] = x60*(x231*x246 + x232*x247 + x233*x248 - x237*x245);
    res_0[5] = x60*(x238*x248 + x239*x247 + x240*x246 - x244*x245);
    res_0[6] = x60*(x121*x250 + x137*x251 + x200*x252 - x229*x249);
    res_0[7] = x60*(x231*x251 + x232*x252 + x233*x250 - x237*x249);
    res_0[8] = x60*(x238*x250 + x239*x252 + x240*x251 - x244*x249);
    res_0[9] = x60*(x121*x256 + x137*x254 + x200*x255 - x229*x253);
    res_0[10] = x60*(x231*x254 + x232*x255 + x233*x256 - x237*x253);
    res_0[11] = x60*(x238*x256 + x239*x255 + x240*x254 - x244*x253);
    res_0[12] = x60*(x121*x259 + x137*x262 + x200*x265 - x229*x266);
    res_0[13] = x60*(x231*x262 + x232*x265 + x233*x259 - x237*x266);
    res_0[14] = x60*(x238*x259 + x239*x265 + x240*x262 - x244*x266);
    res_0[15] = x60*(x121*x268 + x137*x270 + x200*x272 - x229*x273);
    res_0[16] = x60*(x231*x270 + x232*x272 + x233*x268 - x237*x273);
    res_0[17] = x60*(x238*x268 + x239*x272 + x240*x270 - x244*x273);
    res_0[18] = x60*(x121*x275 + x137*x277 + x200*x279 - x229*x280);
    res_0[19] = x60*(x231*x277 + x232*x279 + x233*x275 - x237*x280);
    res_0[20] = x60*(x238*x275 + x239*x279 + x240*x277 - x244*x280);
    res_0[21] = x60*(x121*x283 + x137*x286 + x200*x289 - x229*x290);
    res_0[22] = x60*(x231*x286 + x232*x289 + x233*x283 - x237*x290);
    res_0[23] = x60*(x238*x283 + x239*x289 + x240*x286 - x244*x290);
    res_0[24] = x60*(x121*x293 + x137*x291 + x200*x292 - x229*x294);
    res_0[25] = x60*(x231*x291 + x232*x292 + x233*x293 - x237*x294);
    res_0[26] = x60*(x238*x293 + x239*x292 + x240*x291 - x244*x294);
    res_0[27] = x60*(x121*x297 + x137*x295 + x200*x296 - x229*x298);
    res_0[28] = x60*(x231*x295 + x232*x296 + x233*x297 - x237*x298);
    res_0[29] = x60*(x238*x297 + x239*x296 + x240*x295 - x244*x298);
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
    double *x2 = coord;
    double x3 = 4*s;
    double x4 = -x3;
    double x5 = 4*t;
    double x6 = 4 - x5;
    double x7 = -8*r + x4 + x6;
    double x8 = x2[12];
    double x9 = x2[21];
    double x10 = -x5*x9;
    double x11 = x2[15];
    double x12 = x2[24];
    double x13 = x0 + x3 + x5 - 3;
    double x14 = x13*x2[0];
    double x15 = x2[18];
    double x16 = x14 - x15*x3;
    double x17 = x1*x2[3] + x10 + x11*x3 + x12*x5 + x16 + x7*x8;
    double x18 = x3 - 1;
    double x19 = -x0;
    double x20 = -8*s + x19 + x6;
    double x21 = x2[19];
    double x22 = x2[13];
    double x23 = -x0*x22;
    double x24 = x2[16];
    double x25 = x2[28];
    double x26 = x13*x2[1];
    double x27 = x2[22];
    double x28 = x26 - x27*x5;
    double x29 = x0*x24 + x18*x2[7] + x20*x21 + x23 + x25*x5 + x28;
    double x30 = x17*x29;
    double x31 = -x21*x3;
    double x32 = x2[25];
    double x33 = x1*x2[4] + x22*x7 + x24*x3 + x28 + x31 + x32*x5;
    double x34 = -x0*x8;
    double x35 = x2[27];
    double x36 = x0*x11 + x10 + x14 + x15*x20 + x18*x2[6] + x34 + x35*x5;
    double x37 = x33*x36;
    double x38 = x30 - x37;
    double x39 = x5 - 1;
    double x40 = -8*t + x19 + x4 + 4;
    double x41 = x2[23];
    double x42 = x2[20];
    double x43 = -x3*x42;
    double x44 = x2[26];
    double x45 = x2[29];
    double x46 = x13*x2[2];
    double x47 = x2[14];
    double x48 = -x0*x47 + x46;
    double x49 = x0*x44 + x3*x45 + x39*x2[11] + x40*x41 + x43 + x48;
    double x50 = -x41*x5;
    double x51 = x2[17];
    double x52 = x0*x51 + x18*x2[8] + x20*x42 + x45*x5 + x48 + x50;
    double x53 = x0*x12 + x16 + x3*x35 + x34 + x39*x2[9] + x40*x9;
    double x54 = x33*x53;
    double x55 = x1*x2[5] + x3*x51 + x43 + x44*x5 + x46 + x47*x7 + x50;
    double x56 = x0*x32 + x23 + x25*x3 + x26 + x27*x40 + x31 + x39*x2[10];
    double x57 = x36*x56;
    double x58 = x17*x56;
    double x59 = x29*x53;
    double x60 = x30*x49 - x37*x49 + x52*x54 - x52*x58 + x55*x57 - x55*x59;
    double x61 = 1.0/x60;
    double *x62 = U;
    double x63 = x62[22];
    double x64 = x62[19];
    double x65 = -x3*x64;
    double x66 = x62[25];
    double x67 = x62[28];
    double x68 = x13*x62[1];
    double x69 = x62[13];
    double x70 = -x0*x69 + x68;
    double x71 = x61*(x0*x66 + x3*x67 + x39*x62[10] + x40*x63 + x65 + x70);
    double x72 = x54 - x58;
    double x73 = -x5*x63;
    double x74 = x62[16];
    double x75 = x61*(x0*x74 + x18*x62[7] + x20*x64 + x5*x67 + x70 + x73);
    double x76 = x57 - x59;
    double x77 = x61*(x1*x62[4] + x3*x74 + x5*x66 + x65 + x68 + x69*x7 + x73);
    double x78 = x38*x71 + x72*x75 + x76*x77;
    double x79 = -x29*x55 + x33*x52;
    double x80 = x62[23];
    double x81 = x62[20];
    double x82 = -x3*x81;
    double x83 = x62[26];
    double x84 = x62[29];
    double x85 = x13*x62[2];
    double x86 = x62[14];
    double x87 = -x0*x86 + x85;
    double x88 = x61*(x0*x83 + x3*x84 + x39*x62[11] + x40*x80 + x82 + x87);
    double x89 = -x33*x49 + x55*x56;
    double x90 = -x5*x80;
    double x91 = x62[17];
    double x92 = x61*(x0*x91 + x18*x62[8] + x20*x81 + x5*x84 + x87 + x90);
    double x93 = x29*x49 - x52*x56;
    double x94 = x61*(x1*x62[5] + x3*x91 + x5*x83 + x7*x86 + x82 + x85 + x90);
    double x95 = x79*x88 + x89*x92 + x93*x94;
    double x96 = x71*x79 + x75*x89 + x77*x93;
    double x97 = x38*x88 + x72*x92 + x76*x94 + 1.0;
    double x98 = x78*x95 - x96*x97;
    double x99 = x62[21];
    double x100 = x62[18];
    double x101 = -x100*x3;
    double x102 = x62[24];
    double x103 = x62[27];
    double x104 = x13*x62[0];
    double x105 = x62[12];
    double x106 = -x0*x105 + x104;
    double x107 = x61*(x0*x102 + x101 + x103*x3 + x106 + x39*x62[9] + x40*x99);
    double x108 = -x5*x99;
    double x109 = x62[15];
    double x110 = x61*(x0*x109 + x100*x20 + x103*x5 + x106 + x108 + x18*x62[6]);
    double x111 = x61*(x1*x62[3] + x101 + x102*x5 + x104 + x105*x7 + x108 + x109*x3);
    double x112 = x107*x79 + x110*x89 + x111*x93 + 1.0;
    double x113 = x107*x38 + x110*x72 + x111*x76;
    double x114 = x112*x97 - x113*x95;
    double x115 = -x112*x78 + x113*x96;
    double x116 = S11*x98 + S12*x114 + S13*x115;
    double x117 = S12*x98 + S22*x114 + S23*x115;
    double x118 = S13*x98 + S23*x114 + S33*x115;
    double x119 = -x112*x116 - x117*x96 - x118*x95;
    double x120 = x17*x49 - x53*x55;
    double x121 = x13*x61;
    double x122 = -x17*x52 + x36*x55;
    double x123 = -x36*x49 + x52*x53;
    double x124 = x120*x121 + x121*x122 + x121*x123;
    double x125 = x120*x92 + x122*x88 + x123*x94;
    double x126 = x120*x75 + x122*x71 + x123*x77 + 1.0;
    double x127 = x125*x96 - x126*x95;
    double x128 = x107*x122 + x110*x120 + x111*x123;
    double x129 = -x112*x125 + x128*x95;
    double x130 = x112*x126 - x128*x96;
    double x131 = S11*x127 + S12*x129 + S13*x130;
    double x132 = S12*x127 + S22*x129 + S23*x130;
    double x133 = S13*x127 + S23*x129 + S33*x130;
    double x134 = -x112*x131 - x132*x96 - x133*x95;
    double x135 = x121*x38 + x121*x72 + x121*x76;
    double x136 = x121*x79 + x121*x89 + x121*x93;
    double x137 = -x125*x78 + x126*x97;
    double x138 = x113*x125 - x128*x97;
    double x139 = -x113*x126 + x128*x78;
    double x140 = S11*x137 + S12*x138 + S13*x139;
    double x141 = S12*x137 + S22*x138 + S23*x139;
    double x142 = S13*x137 + S23*x138 + S33*x139;
    double x143 = 1.0*PENER + 1.0*SENER;
    double x144 = -x112*x140 - x141*x96 - x142*x95 + x143;
    double x145 = -x125*x133 - x126*x132 - x128*x131;
    double x146 = -x125*x142 - x126*x141 - x128*x140;
    double x147 = -x116*x128 - x117*x126 - x118*x125 + x143;
    double x148 = -x113*x116 - x117*x78 - x118*x97;
    double x149 = -x113*x140 - x141*x78 - x142*x97;
    double x150 = -x113*x131 - x132*x78 - x133*x97 + x143;
    double x151 = x1*x61;
    double x152 = x151*x76;
    double x153 = x151*x93;
    double x154 = x123*x151;
    double x155 = x18*x61;
    double x156 = x120*x155;
    double x157 = x155*x72;
    double x158 = x155*x89;
    double x159 = x39*x61;
    double x160 = x159*x38;
    double x161 = x159*x79;
    double x162 = x122*x159;
    double x163 = x61*x7;
    double x164 = x0*x61;
    double x165 = x120*x164;
    double x166 = x122*x164;
    double x167 = x123*x163 - x165 - x166;
    double x168 = x164*x38;
    double x169 = x164*x72;
    double x170 = x163*x76 - x168 - x169;
    double x171 = x164*x79;
    double x172 = x164*x89;
    double x173 = x163*x93 - x171 - x172;
    double x174 = x3*x61;
    double x175 = x123*x174;
    double x176 = x165 + x175;
    double x177 = x174*x76;
    double x178 = x169 + x177;
    double x179 = x174*x93;
    double x180 = x172 + x179;
    double x181 = x20*x61;
    double x182 = x122*x174;
    double x183 = x120*x181 - x175 - x182;
    double x184 = x174*x38;
    double x185 = -x177 + x181*x72 - x184;
    double x186 = x174*x79;
    double x187 = -x179 + x181*x89 - x186;
    double x188 = x40*x61;
    double x189 = x5*x61;
    double x190 = x120*x189;
    double x191 = x123*x189;
    double x192 = x122*x188 - x190 - x191;
    double x193 = x189*x72;
    double x194 = x189*x76;
    double x195 = x188*x38 - x193 - x194;
    double x196 = x189*x89;
    double x197 = x189*x93;
    double x198 = x188*x79 - x196 - x197;
    double x199 = x168 + x194;
    double x200 = x171 + x197;
    double x201 = x166 + x191;
    double x202 = x184 + x193;
    double x203 = x186 + x196;
    double x204 = x182 + x190;
    
    res_0[0] = x60*(x119*x124 + x134*x135 + x136*x144);
    res_0[1] = x60*(x124*x147 + x135*x145 + x136*x146);
    res_0[2] = x60*(x124*x148 + x135*x150 + x136*x149);
    res_0[3] = x60*(x119*x154 + x134*x152 + x144*x153);
    res_0[4] = x60*(x145*x152 + x146*x153 + x147*x154);
    res_0[5] = x60*(x148*x154 + x149*x153 + x150*x152);
    res_0[6] = x60*(x119*x156 + x134*x157 + x144*x158);
    res_0[7] = x60*(x145*x157 + x146*x158 + x147*x156);
    res_0[8] = x60*(x148*x156 + x149*x158 + x150*x157);
    res_0[9] = x60*(x119*x162 + x134*x160 + x144*x161);
    res_0[10] = x60*(x145*x160 + x146*x161 + x147*x162);
    res_0[11] = x60*(x148*x162 + x149*x161 + x150*x160);
    res_0[12] = x60*(x119*x167 + x134*x170 + x144*x173);
    res_0[13] = x60*(x145*x170 + x146*x173 + x147*x167);
    res_0[14] = x60*(x148*x167 + x149*x173 + x150*x170);
    res_0[15] = x60*(x119*x176 + x134*x178 + x144*x180);
    res_0[16] = x60*(x145*x178 + x146*x180 + x147*x176);
    res_0[17] = x60*(x148*x176 + x149*x180 + x150*x178);
    res_0[18] = x60*(x119*x183 + x134*x185 + x144*x187);
    res_0[19] = x60*(x145*x185 + x146*x187 + x147*x183);
    res_0[20] = x60*(x148*x183 + x149*x187 + x150*x185);
    res_0[21] = x60*(x119*x192 + x134*x195 + x144*x198);
    res_0[22] = x60*(x145*x195 + x146*x198 + x147*x192);
    res_0[23] = x60*(x148*x192 + x149*x198 + x150*x195);
    res_0[24] = x60*(x119*x201 + x134*x199 + x144*x200);
    res_0[25] = x60*(x145*x199 + x146*x200 + x147*x201);
    res_0[26] = x60*(x148*x201 + x149*x200 + x150*x199);
    res_0[27] = x60*(x119*x204 + x134*x202 + x144*x203);
    res_0[28] = x60*(x145*x202 + x146*x203 + x147*x204);
    res_0[29] = x60*(x148*x204 + x149*x203 + x150*x202);
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
    double *x2 = coord;
    double x3 = 4*s;
    double x4 = -x3;
    double x5 = 4*t;
    double x6 = 4 - x5;
    double x7 = -8*r + x4 + x6;
    double x8 = x2[12];
    double x9 = x2[21];
    double x10 = -x5*x9;
    double x11 = x2[15];
    double x12 = x2[24];
    double x13 = x0 + x3 + x5 - 3;
    double x14 = x13*x2[0];
    double x15 = x2[18];
    double x16 = x14 - x15*x3;
    double x17 = x1*x2[3] + x10 + x11*x3 + x12*x5 + x16 + x7*x8;
    double x18 = x3 - 1;
    double x19 = -x0;
    double x20 = -8*s + x19 + x6;
    double x21 = x2[19];
    double x22 = x2[13];
    double x23 = -x0*x22;
    double x24 = x2[16];
    double x25 = x2[28];
    double x26 = x13*x2[1];
    double x27 = x2[22];
    double x28 = x26 - x27*x5;
    double x29 = x0*x24 + x18*x2[7] + x20*x21 + x23 + x25*x5 + x28;
    double x30 = x17*x29;
    double x31 = -x21*x3;
    double x32 = x2[25];
    double x33 = x1*x2[4] + x22*x7 + x24*x3 + x28 + x31 + x32*x5;
    double x34 = -x0*x8;
    double x35 = x2[27];
    double x36 = x0*x11 + x10 + x14 + x15*x20 + x18*x2[6] + x34 + x35*x5;
    double x37 = x33*x36;
    double x38 = x30 - x37;
    double x39 = x5 - 1;
    double x40 = -8*t + x19 + x4 + 4;
    double x41 = x2[23];
    double x42 = x2[20];
    double x43 = -x3*x42;
    double x44 = x2[26];
    double x45 = x2[29];
    double x46 = x13*x2[2];
    double x47 = x2[14];
    double x48 = -x0*x47 + x46;
    double x49 = x0*x44 + x3*x45 + x39*x2[11] + x40*x41 + x43 + x48;
    double x50 = -x41*x5;
    double x51 = x2[17];
    double x52 = x0*x51 + x18*x2[8] + x20*x42 + x45*x5 + x48 + x50;
    double x53 = x0*x12 + x16 + x3*x35 + x34 + x39*x2[9] + x40*x9;
    double x54 = x33*x53;
    double x55 = x1*x2[5] + x3*x51 + x43 + x44*x5 + x46 + x47*x7 + x50;
    double x56 = x0*x32 + x23 + x25*x3 + x26 + x27*x40 + x31 + x39*x2[10];
    double x57 = x36*x56;
    double x58 = x17*x56;
    double x59 = x29*x53;
    double x60 = x30*x49 - x37*x49 + x52*x54 - x52*x58 + x55*x57 - x55*x59;
    double x61 = 1.0/x60;
    double *x62 = U;
    double x63 = x62[22];
    double x64 = x62[19];
    double x65 = -x3*x64;
    double x66 = x62[25];
    double x67 = x62[28];
    double x68 = x13*x62[1];
    double x69 = x62[13];
    double x70 = -x0*x69 + x68;
    double x71 = x61*(x0*x66 + x3*x67 + x39*x62[10] + x40*x63 + x65 + x70);
    double x72 = x54 - x58;
    double x73 = -x5*x63;
    double x74 = x62[16];
    double x75 = x61*(x0*x74 + x18*x62[7] + x20*x64 + x5*x67 + x70 + x73);
    double x76 = x57 - x59;
    double x77 = x61*(x1*x62[4] + x3*x74 + x5*x66 + x65 + x68 + x69*x7 + x73);
    double x78 = x38*x71 + x72*x75 + x76*x77;
    double x79 = -x29*x55 + x33*x52;
    double x80 = x62[23];
    double x81 = x62[20];
    double x82 = -x3*x81;
    double x83 = x62[26];
    double x84 = x62[29];
    double x85 = x13*x62[2];
    double x86 = x62[14];
    double x87 = -x0*x86 + x85;
    double x88 = x61*(x0*x83 + x3*x84 + x39*x62[11] + x40*x80 + x82 + x87);
    double x89 = -x33*x49 + x55*x56;
    double x90 = -x5*x80;
    double x91 = x62[17];
    double x92 = x61*(x0*x91 + x18*x62[8] + x20*x81 + x5*x84 + x87 + x90);
    double x93 = x29*x49 - x52*x56;
    double x94 = x61*(x1*x62[5] + x3*x91 + x5*x83 + x7*x86 + x82 + x85 + x90);
    double x95 = x79*x88 + x89*x92 + x93*x94;
    double x96 = x71*x79 + x75*x89 + x77*x93;
    double x97 = x38*x88 + x72*x92 + x76*x94;
    double x98 = x97 + 1.0;
    double x99 = x78*x95 - x96*x98;
    double x100 = x62[21];
    double x101 = x62[18];
    double x102 = -x101*x3;
    double x103 = x62[24];
    double x104 = x62[27];
    double x105 = x13*x62[0];
    double x106 = x62[12];
    double x107 = -x0*x106 + x105;
    double x108 = x61*(x0*x103 + x100*x40 + x102 + x104*x3 + x107 + x39*x62[9]);
    double x109 = -x100*x5;
    double x110 = x62[15];
    double x111 = x61*(x0*x110 + x101*x20 + x104*x5 + x107 + x109 + x18*x62[6]);
    double x112 = x61*(x1*x62[3] + x102 + x103*x5 + x105 + x106*x7 + x109 + x110*x3);
    double x113 = x108*x79 + x111*x89 + x112*x93;
    double x114 = x113 + 1.0;
    double x115 = x108*x38 + x111*x72 + x112*x76;
    double x116 = x114*x98 - x115*x95;
    double x117 = -x114*x78 + x115*x96;
    double x118 = S11*x99 + S12*x116 + S13*x117;
    double x119 = S12*x99 + S22*x116 + S23*x117;
    double x120 = S13*x99 + S23*x116 + S33*x117;
    double x121 = -x113*x118 - x119*x96 - x120*x95;
    double x122 = x17*x49 - x53*x55;
    double x123 = x13*x61;
    double x124 = -x17*x52 + x36*x55;
    double x125 = -x36*x49 + x52*x53;
    double x126 = x122*x123 + x123*x124 + x123*x125;
    double x127 = x122*x92 + x124*x88 + x125*x94;
    double x128 = x122*x75 + x124*x71 + x125*x77;
    double x129 = x128 + 1.0;
    double x130 = x127*x96 - x129*x95;
    double x131 = x108*x124 + x111*x122 + x112*x125;
    double x132 = -x114*x127 + x131*x95;
    double x133 = x114*x129 - x131*x96;
    double x134 = S11*x130 + S12*x132 + S13*x133;
    double x135 = S12*x130 + S22*x132 + S23*x133;
    double x136 = S13*x130 + S23*x132 + S33*x133;
    double x137 = -x113*x134 - x135*x96 - x136*x95;
    double x138 = x123*x38 + x123*x72 + x123*x76;
    double x139 = x123*x79 + x123*x89 + x123*x93;
    double x140 = -x127*x78 + x129*x98;
    double x141 = x115*x127 - x131*x98;
    double x142 = -x115*x129 + x131*x78;
    double x143 = S11*x140 + S12*x141 + S13*x142;
    double x144 = S12*x140 + S22*x141 + S23*x142;
    double x145 = S13*x140 + S23*x141 + S33*x142;
    double x146 = 1.0*PENER + 1.0*SENER;
    double x147 = -x113*x143 - x144*x96 - x145*x95 + x146;
    double x148 = -x127*x136 - x128*x135 - x131*x134;
    double x149 = -x127*x145 - x128*x144 - x131*x143;
    double x150 = -x118*x131 - x119*x128 - x120*x127 + x146;
    double x151 = -x115*x118 - x119*x78 - x120*x97;
    double x152 = -x115*x143 - x144*x78 - x145*x97;
    double x153 = -x115*x134 - x135*x78 - x136*x97 + x146;
    double x154 = x1*x61;
    double x155 = x154*x76;
    double x156 = x154*x93;
    double x157 = x125*x154;
    double x158 = x18*x61;
    double x159 = x122*x158;
    double x160 = x158*x72;
    double x161 = x158*x89;
    double x162 = x39*x61;
    double x163 = x162*x38;
    double x164 = x162*x79;
    double x165 = x124*x162;
    double x166 = x61*x7;
    double x167 = x0*x61;
    double x168 = x122*x167;
    double x169 = x124*x167;
    double x170 = x125*x166 - x168 - x169;
    double x171 = x167*x38;
    double x172 = x167*x72;
    double x173 = x166*x76 - x171 - x172;
    double x174 = x167*x79;
    double x175 = x167*x89;
    double x176 = x166*x93 - x174 - x175;
    double x177 = x3*x61;
    double x178 = x125*x177;
    double x179 = x168 + x178;
    double x180 = x177*x76;
    double x181 = x172 + x180;
    double x182 = x177*x93;
    double x183 = x175 + x182;
    double x184 = x20*x61;
    double x185 = x124*x177;
    double x186 = x122*x184 - x178 - x185;
    double x187 = x177*x38;
    double x188 = -x180 + x184*x72 - x187;
    double x189 = x177*x79;
    double x190 = -x182 + x184*x89 - x189;
    double x191 = x40*x61;
    double x192 = x5*x61;
    double x193 = x122*x192;
    double x194 = x125*x192;
    double x195 = x124*x191 - x193 - x194;
    double x196 = x192*x72;
    double x197 = x192*x76;
    double x198 = x191*x38 - x196 - x197;
    double x199 = x192*x89;
    double x200 = x192*x93;
    double x201 = x191*x79 - x199 - x200;
    double x202 = x171 + x197;
    double x203 = x174 + x200;
    double x204 = x169 + x194;
    double x205 = x187 + x196;
    double x206 = x189 + x199;
    double x207 = x185 + x193;
    
    res_0[0] = x60*(x121*x126 + x137*x138 + x139*x147);
    res_0[1] = x60*(x126*x150 + x138*x148 + x139*x149);
    res_0[2] = x60*(x126*x151 + x138*x153 + x139*x152);
    res_0[3] = x60*(x121*x157 + x137*x155 + x147*x156);
    res_0[4] = x60*(x148*x155 + x149*x156 + x150*x157);
    res_0[5] = x60*(x151*x157 + x152*x156 + x153*x155);
    res_0[6] = x60*(x121*x159 + x137*x160 + x147*x161);
    res_0[7] = x60*(x148*x160 + x149*x161 + x150*x159);
    res_0[8] = x60*(x151*x159 + x152*x161 + x153*x160);
    res_0[9] = x60*(x121*x165 + x137*x163 + x147*x164);
    res_0[10] = x60*(x148*x163 + x149*x164 + x150*x165);
    res_0[11] = x60*(x151*x165 + x152*x164 + x153*x163);
    res_0[12] = x60*(x121*x170 + x137*x173 + x147*x176);
    res_0[13] = x60*(x148*x173 + x149*x176 + x150*x170);
    res_0[14] = x60*(x151*x170 + x152*x176 + x153*x173);
    res_0[15] = x60*(x121*x179 + x137*x181 + x147*x183);
    res_0[16] = x60*(x148*x181 + x149*x183 + x150*x179);
    res_0[17] = x60*(x151*x179 + x152*x183 + x153*x181);
    res_0[18] = x60*(x121*x186 + x137*x188 + x147*x190);
    res_0[19] = x60*(x148*x188 + x149*x190 + x150*x186);
    res_0[20] = x60*(x151*x186 + x152*x190 + x153*x188);
    res_0[21] = x60*(x121*x195 + x137*x198 + x147*x201);
    res_0[22] = x60*(x148*x198 + x149*x201 + x150*x195);
    res_0[23] = x60*(x151*x195 + x152*x201 + x153*x198);
    res_0[24] = x60*(x121*x204 + x137*x202 + x147*x203);
    res_0[25] = x60*(x148*x202 + x149*x203 + x150*x204);
    res_0[26] = x60*(x151*x204 + x152*x203 + x153*x202);
    res_0[27] = x60*(x121*x207 + x137*x205 + x147*x206);
    res_0[28] = x60*(x148*x205 + x149*x206 + x150*x207);
    res_0[29] = x60*(x151*x207 + x152*x206 + x153*x205);
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
    double *x2 = coord;
    double x3 = 4*s;
    double x4 = -x3;
    double x5 = 4*t;
    double x6 = 4 - x5;
    double x7 = -8*r + x4 + x6;
    double x8 = x2[12];
    double x9 = x2[21];
    double x10 = -x5*x9;
    double x11 = x2[15];
    double x12 = x2[24];
    double x13 = x0 + x3 + x5 - 3;
    double x14 = x13*x2[0];
    double x15 = x2[18];
    double x16 = x14 - x15*x3;
    double x17 = x1*x2[3] + x10 + x11*x3 + x12*x5 + x16 + x7*x8;
    double x18 = x3 - 1;
    double x19 = -x0;
    double x20 = -8*s + x19 + x6;
    double x21 = x2[19];
    double x22 = x2[13];
    double x23 = -x0*x22;
    double x24 = x2[16];
    double x25 = x2[28];
    double x26 = x13*x2[1];
    double x27 = x2[22];
    double x28 = x26 - x27*x5;
    double x29 = x0*x24 + x18*x2[7] + x20*x21 + x23 + x25*x5 + x28;
    double x30 = x17*x29;
    double x31 = -x21*x3;
    double x32 = x2[25];
    double x33 = x1*x2[4] + x22*x7 + x24*x3 + x28 + x31 + x32*x5;
    double x34 = -x0*x8;
    double x35 = x2[27];
    double x36 = x0*x11 + x10 + x14 + x15*x20 + x18*x2[6] + x34 + x35*x5;
    double x37 = x33*x36;
    double x38 = x30 - x37;
    double x39 = x5 - 1;
    double x40 = -8*t + x19 + x4 + 4;
    double x41 = x2[23];
    double x42 = x2[20];
    double x43 = -x3*x42;
    double x44 = x2[26];
    double x45 = x2[29];
    double x46 = x13*x2[2];
    double x47 = x2[14];
    double x48 = -x0*x47 + x46;
    double x49 = x0*x44 + x3*x45 + x39*x2[11] + x40*x41 + x43 + x48;
    double x50 = -x41*x5;
    double x51 = x2[17];
    double x52 = x0*x51 + x18*x2[8] + x20*x42 + x45*x5 + x48 + x50;
    double x53 = x0*x12 + x16 + x3*x35 + x34 + x39*x2[9] + x40*x9;
    double x54 = x33*x53;
    double x55 = x1*x2[5] + x3*x51 + x43 + x44*x5 + x46 + x47*x7 + x50;
    double x56 = x0*x32 + x23 + x25*x3 + x26 + x27*x40 + x31 + x39*x2[10];
    double x57 = x36*x56;
    double x58 = x17*x56;
    double x59 = x29*x53;
    double x60 = x30*x49 - x37*x49 + x52*x54 - x52*x58 + x55*x57 - x55*x59;
    double x61 = 1.0/x60;
    double *x62 = U;
    double x63 = x62[22];
    double x64 = x62[19];
    double x65 = -x3*x64;
    double x66 = x62[25];
    double x67 = x62[28];
    double x68 = x13*x62[1];
    double x69 = x62[13];
    double x70 = -x0*x69 + x68;
    double x71 = x61*(x0*x66 + x3*x67 + x39*x62[10] + x40*x63 + x65 + x70);
    double x72 = x54 - x58;
    double x73 = -x5*x63;
    double x74 = x62[16];
    double x75 = x61*(x0*x74 + x18*x62[7] + x20*x64 + x5*x67 + x70 + x73);
    double x76 = x57 - x59;
    double x77 = x61*(x1*x62[4] + x3*x74 + x5*x66 + x65 + x68 + x69*x7 + x73);
    double x78 = x38*x71 + x72*x75 + x76*x77;
    double x79 = -x29*x55 + x33*x52;
    double x80 = x62[23];
    double x81 = x62[20];
    double x82 = -x3*x81;
    double x83 = x62[26];
    double x84 = x62[29];
    double x85 = x13*x62[2];
    double x86 = x62[14];
    double x87 = -x0*x86 + x85;
    double x88 = x61*(x0*x83 + x3*x84 + x39*x62[11] + x40*x80 + x82 + x87);
    double x89 = -x33*x49 + x55*x56;
    double x90 = -x5*x80;
    double x91 = x62[17];
    double x92 = x61*(x0*x91 + x18*x62[8] + x20*x81 + x5*x84 + x87 + x90);
    double x93 = x29*x49 - x52*x56;
    double x94 = x61*(x1*x62[5] + x3*x91 + x5*x83 + x7*x86 + x82 + x85 + x90);
    double x95 = x79*x88 + x89*x92 + x93*x94;
    double x96 = x71*x79 + x75*x89 + x77*x93;
    double x97 = x38*x88 + x72*x92 + x76*x94;
    double x98 = x97 + 1.0;
    double x99 = x78*x95 - x96*x98;
    double x100 = x62[21];
    double x101 = x62[18];
    double x102 = -x101*x3;
    double x103 = x62[24];
    double x104 = x62[27];
    double x105 = x13*x62[0];
    double x106 = x62[12];
    double x107 = -x0*x106 + x105;
    double x108 = x61*(x0*x103 + x100*x40 + x102 + x104*x3 + x107 + x39*x62[9]);
    double x109 = -x100*x5;
    double x110 = x62[15];
    double x111 = x61*(x0*x110 + x101*x20 + x104*x5 + x107 + x109 + x18*x62[6]);
    double x112 = x61*(x1*x62[3] + x102 + x103*x5 + x105 + x106*x7 + x109 + x110*x3);
    double x113 = x108*x79 + x111*x89 + x112*x93;
    double x114 = x113 + 1.0;
    double x115 = x108*x38 + x111*x72 + x112*x76;
    double x116 = x114*x98 - x115*x95;
    double x117 = -x114*x78 + x115*x96;
    double x118 = S11*x99 + S12*x116 + S13*x117;
    double x119 = S12*x99 + S22*x116 + S23*x117;
    double x120 = S13*x99 + S23*x116 + S33*x117;
    double x121 = -x113*x118 - x119*x96 - x120*x95;
    double x122 = x17*x49 - x53*x55;
    double x123 = x13*x61;
    double x124 = -x17*x52 + x36*x55;
    double x125 = -x36*x49 + x52*x53;
    double x126 = x122*x123 + x123*x124 + x123*x125;
    double x127 = x122*x92 + x124*x88 + x125*x94;
    double x128 = x122*x75 + x124*x71 + x125*x77;
    double x129 = x128 + 1.0;
    double x130 = x127*x96 - x129*x95;
    double x131 = x108*x124 + x111*x122 + x112*x125;
    double x132 = -x114*x127 + x131*x95;
    double x133 = x114*x129 - x131*x96;
    double x134 = S11*x130 + S12*x132 + S13*x133;
    double x135 = S12*x130 + S22*x132 + S23*x133;
    double x136 = S13*x130 + S23*x132 + S33*x133;
    double x137 = -x113*x134 - x135*x96 - x136*x95;
    double x138 = x123*x38 + x123*x72 + x123*x76;
    double x139 = x123*x79 + x123*x89 + x123*x93;
    double x140 = -x127*x78 + x129*x98;
    double x141 = x115*x127 - x131*x98;
    double x142 = -x115*x129 + x131*x78;
    double x143 = S11*x140 + S12*x141 + S13*x142;
    double x144 = S12*x140 + S22*x141 + S23*x142;
    double x145 = S13*x140 + S23*x141 + S33*x142;
    double x146 = r*r;
    double x147 = 2*x146;
    double x148 = -r + x147;
    double *x149 = V;
    double x150 = x149[3];
    double x151 = s*s;
    double x152 = 2*x151;
    double x153 = -s + x152;
    double x154 = x149[6];
    double x155 = t*t;
    double x156 = 2*x155;
    double x157 = -t + x156;
    double x158 = x149[9];
    double x159 = s*x0;
    double x160 = -x159;
    double x161 = t*x0;
    double x162 = -x161;
    double x163 = x0 - 4*x146 + x160 + x162;
    double x164 = x149[12];
    double x165 = t*x3;
    double x166 = -x165;
    double x167 = -4*x151 + x160 + x166 + x3;
    double x168 = x149[18];
    double x169 = -4*x155 + x162 + x166 + x5;
    double x170 = x149[21];
    double x171 = -3*r - 3*s - 3*t + x147 + x152 + x156 + x159 + x161 + x165 + 1;
    double x172 = x149[0];
    double x173 = x149[15];
    double x174 = x149[24];
    double x175 = x149[27];
    double x176 = x148*x150 + x153*x154 + x157*x158 + x159*x173 + x161*x174 + x163*x164 + x165*x175 + x167*x168 + x169*x170 + x171*x172;
    double x177 = x149[4];
    double x178 = x149[7];
    double x179 = x149[10];
    double x180 = x149[13];
    double x181 = x149[19];
    double x182 = x149[22];
    double x183 = x149[1];
    double x184 = x149[16];
    double x185 = x149[25];
    double x186 = x149[28];
    double x187 = x148*x177 + x153*x178 + x157*x179 + x159*x184 + x161*x185 + x163*x180 + x165*x186 + x167*x181 + x169*x182 + x171*x183;
    double x188 = x149[5];
    double x189 = x149[8];
    double x190 = x149[11];
    double x191 = x149[14];
    double x192 = x149[20];
    double x193 = x149[23];
    double x194 = x149[2];
    double x195 = x149[17];
    double x196 = x149[26];
    double x197 = x149[29];
    double x198 = x148*x188 + x153*x189 + x157*x190 + x159*x195 + x161*x196 + x163*x191 + x165*x197 + x167*x192 + x169*x193 + x171*x194;
    double x199 = 1.0*PENER + 1.0*SENER - 0.5*rho*(x176*x176 + x187*x187 + x198*x198);
    double x200 = -x113*x143 - x144*x96 - x145*x95 + x199;
    double *x201 = A;
    double x202 = x148*x201[3] + x153*x201[6] + x157*x201[9] + x159*x201[15] + x161*x201[24] + x163*x201[12] + x165*x201[27] + x167*x201[18] + x169*x201[21] + x171*x201[0];
    double x203 = x148*x201[4] + x153*x201[7] + x157*x201[10] + x159*x201[16] + x161*x201[25] + x163*x201[13] + x165*x201[28] + x167*x201[19] + x169*x201[22] + x171*x201[1];
    double x204 = x148*x201[5] + x153*x201[8] + x157*x201[11] + x159*x201[17] + x161*x201[26] + x163*x201[14] + x165*x201[29] + x167*x201[20] + x169*x201[23] + x171*x201[2];
    double x205 = -x168*x3;
    double x206 = x13*x172;
    double x207 = -x0*x164 + x206;
    double x208 = x0*x174 + x158*x39 + x170*x40 + x175*x3 + x205 + x207;
    double x209 = x61*x79;
    double x210 = -x170*x5;
    double x211 = x0*x173 + x154*x18 + x168*x20 + x175*x5 + x207 + x210;
    double x212 = x61*x89;
    double x213 = x1*x150 + x164*x7 + x173*x3 + x174*x5 + x205 + x206 + x210;
    double x214 = x61*x93;
    double x215 = -x181*x3;
    double x216 = x13*x183;
    double x217 = -x0*x180 + x216;
    double x218 = x0*x185 + x179*x39 + x182*x40 + x186*x3 + x215 + x217;
    double x219 = -x182*x5;
    double x220 = x0*x184 + x178*x18 + x181*x20 + x186*x5 + x217 + x219;
    double x221 = x1*x177 + x180*x7 + x184*x3 + x185*x5 + x215 + x216 + x219;
    double x222 = -x192*x3;
    double x223 = x13*x194;
    double x224 = -x0*x191 + x223;
    double x225 = x0*x196 + x190*x39 + x193*x40 + x197*x3 + x222 + x224;
    double x226 = -x193*x5;
    double x227 = x0*x195 + x18*x189 + x192*x20 + x197*x5 + x224 + x226;
    double x228 = x1*x188 + x191*x7 + x195*x3 + x196*x5 + x222 + x223 + x226;
    double x229 = x113*x202 + x176*(x208*x209 + x211*x212 + x213*x214) + x187*(x209*x218 + x212*x220 + x214*x221) + x198*(x209*x225 + x212*x227 + x214*x228) + x203*x96 + x204*x95;
    double x230 = rho*x171;
    double x231 = -x127*x136 - x128*x135 - x131*x134;
    double x232 = -x127*x145 - x128*x144 - x131*x143;
    double x233 = -x118*x131 - x119*x128 - x120*x127 + x199;
    double x234 = x122*x61;
    double x235 = x124*x61;
    double x236 = x125*x61;
    double x237 = x127*x204 + x128*x203 + x131*x202 + x176*(x208*x235 + x211*x234 + x213*x236) + x187*(x218*x235 + x220*x234 + x221*x236) + x198*(x225*x235 + x227*x234 + x228*x236);
    double x238 = -x115*x118 - x119*x78 - x120*x97;
    double x239 = -x115*x143 - x144*x78 - x145*x97;
    double x240 = -x115*x134 - x135*x78 - x136*x97 + x199;
    double x241 = x38*x61;
    double x242 = x61*x72;
    double x243 = x61*x76;
    double x244 = x115*x202 + x176*(x208*x241 + x211*x242 + x213*x243) + x187*(x218*x241 + x220*x242 + x221*x243) + x198*(x225*x241 + x227*x242 + x228*x243) + x203*x78 + x204*x97;
    double x245 = rho*x148;
    double x246 = x1*x243;
    double x247 = x1*x214;
    double x248 = x1*x236;
    double x249 = rho*x153;
    double x250 = x18*x234;
    double x251 = x18*x242;
    double x252 = x18*x212;
    double x253 = rho*x157;
    double x254 = x241*x39;
    double x255 = x209*x39;
    double x256 = x235*x39;
    double x257 = x0*x234;
    double x258 = x0*x235;
    double x259 = x236*x7 - x257 - x258;
    double x260 = x0*x241;
    double x261 = x0*x242;
    double x262 = x243*x7 - x260 - x261;
    double x263 = x0*x209;
    double x264 = x0*x212;
    double x265 = x214*x7 - x263 - x264;
    double x266 = rho*x163;
    double x267 = x236*x3;
    double x268 = x257 + x267;
    double x269 = x243*x3;
    double x270 = x261 + x269;
    double x271 = x214*x3;
    double x272 = x264 + x271;
    double x273 = rho*x159;
    double x274 = x235*x3;
    double x275 = x20*x234 - x267 - x274;
    double x276 = x241*x3;
    double x277 = x20*x242 - x269 - x276;
    double x278 = x209*x3;
    double x279 = x20*x212 - x271 - x278;
    double x280 = rho*x167;
    double x281 = x234*x5;
    double x282 = x236*x5;
    double x283 = x235*x40 - x281 - x282;
    double x284 = x242*x5;
    double x285 = x243*x5;
    double x286 = x241*x40 - x284 - x285;
    double x287 = x212*x5;
    double x288 = x214*x5;
    double x289 = x209*x40 - x287 - x288;
    double x290 = rho*x169;
    double x291 = x260 + x285;
    double x292 = x263 + x288;
    double x293 = x258 + x282;
    double x294 = rho*x161;
    double x295 = x276 + x284;
    double x296 = x278 + x287;
    double x297 = x274 + x281;
    double x298 = rho*x165;
    
    res_0[0] = x60*(x121*x126 + x137*x138 + x139*x200 - x229*x230);
    res_0[1] = x60*(x126*x233 + x138*x231 + x139*x232 - x230*x237);
    res_0[2] = x60*(x126*x238 + x138*x240 + x139*x239 - x230*x244);
    res_0[3] = x60*(x121*x248 + x137*x246 + x200*x247 - x229*x245);
    res_0[4] = x60*(x231*x246 + x232*x247 + x233*x248 - x237*x245);
    res_0[5] = x60*(x238*x248 + x239*x247 + x240*x246 - x244*x245);
    res_0[6] = x60*(x121*x250 + x137*x251 + x200*x252 - x229*x249);
    res_0[7] = x60*(x231*x251 + x232*x252 + x233*x250 - x237*x249);
    res_0[8] = x60*(x238*x250 + x239*x252 + x240*x251 - x244*x249);
    res_0[9] = x60*(x121*x256 + x137*x254 + x200*x255 - x229*x253);
    res_0[10] = x60*(x231*x254 + x232*x255 + x233*x256 - x237*x253);
    res_0[11] = x60*(x238*x256 + x239*x255 + x240*x254 - x244*x253);
    res_0[12] = x60*(x121*x259 + x137*x262 + x200*x265 - x229*x266);
    res_0[13] = x60*(x231*x262 + x232*x265 + x233*x259 - x237*x266);
    res_0[14] = x60*(x238*x259 + x239*x265 + x240*x262 - x244*x266);
    res_0[15] = x60*(x121*x268 + x137*x270 + x200*x272 - x229*x273);
    res_0[16] = x60*(x231*x270 + x232*x272 + x233*x268 - x237*x273);
    res_0[17] = x60*(x238*x268 + x239*x272 + x240*x270 - x244*x273);
    res_0[18] = x60*(x121*x275 + x137*x277 + x200*x279 - x229*x280);
    res_0[19] = x60*(x231*x277 + x232*x279 + x233*x275 - x237*x280);
    res_0[20] = x60*(x238*x275 + x239*x279 + x240*x277 - x244*x280);
    res_0[21] = x60*(x121*x283 + x137*x286 + x200*x289 - x229*x290);
    res_0[22] = x60*(x231*x286 + x232*x289 + x233*x283 - x237*x290);
    res_0[23] = x60*(x238*x283 + x239*x289 + x240*x286 - x244*x290);
    res_0[24] = x60*(x121*x293 + x137*x291 + x200*x292 - x229*x294);
    res_0[25] = x60*(x231*x291 + x232*x292 + x233*x293 - x237*x294);
    res_0[26] = x60*(x238*x293 + x239*x292 + x240*x291 - x244*x294);
    res_0[27] = x60*(x121*x297 + x137*x295 + x200*x296 - x229*x298);
    res_0[28] = x60*(x231*x295 + x232*x296 + x233*x297 - x237*x298);
    res_0[29] = x60*(x238*x297 + x239*x296 + x240*x295 - x244*x298);
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
    double *x2 = coord;
    double x3 = 4*s;
    double x4 = -x3;
    double x5 = 4*t;
    double x6 = 4 - x5;
    double x7 = -8*r + x4 + x6;
    double x8 = x2[12];
    double x9 = x2[21];
    double x10 = -x5*x9;
    double x11 = x2[15];
    double x12 = x2[24];
    double x13 = x0 + x3 + x5 - 3;
    double x14 = x13*x2[0];
    double x15 = x2[18];
    double x16 = x14 - x15*x3;
    double x17 = x1*x2[3] + x10 + x11*x3 + x12*x5 + x16 + x7*x8;
    double x18 = x3 - 1;
    double x19 = -x0;
    double x20 = -8*s + x19 + x6;
    double x21 = x2[19];
    double x22 = x2[13];
    double x23 = -x0*x22;
    double x24 = x2[16];
    double x25 = x2[28];
    double x26 = x13*x2[1];
    double x27 = x2[22];
    double x28 = x26 - x27*x5;
    double x29 = x0*x24 + x18*x2[7] + x20*x21 + x23 + x25*x5 + x28;
    double x30 = x17*x29;
    double x31 = -x21*x3;
    double x32 = x2[25];
    double x33 = x1*x2[4] + x22*x7 + x24*x3 + x28 + x31 + x32*x5;
    double x34 = -x0*x8;
    double x35 = x2[27];
    double x36 = x0*x11 + x10 + x14 + x15*x20 + x18*x2[6] + x34 + x35*x5;
    double x37 = x33*x36;
    double x38 = x30 - x37;
    double x39 = x5 - 1;
    double x40 = -8*t + x19 + x4 + 4;
    double x41 = x2[23];
    double x42 = x2[20];
    double x43 = -x3*x42;
    double x44 = x2[26];
    double x45 = x2[29];
    double x46 = x13*x2[2];
    double x47 = x2[14];
    double x48 = -x0*x47 + x46;
    double x49 = x0*x44 + x3*x45 + x39*x2[11] + x40*x41 + x43 + x48;
    double x50 = -x41*x5;
    double x51 = x2[17];
    double x52 = x0*x51 + x18*x2[8] + x20*x42 + x45*x5 + x48 + x50;
    double x53 = x0*x12 + x16 + x3*x35 + x34 + x39*x2[9] + x40*x9;
    double x54 = x33*x53;
    double x55 = x1*x2[5] + x3*x51 + x43 + x44*x5 + x46 + x47*x7 + x50;
    double x56 = x0*x32 + x23 + x25*x3 + x26 + x27*x40 + x31 + x39*x2[10];
    double x57 = x36*x56;
    double x58 = x17*x56;
    double x59 = x29*x53;
    double x60 = x30*x49 - x37*x49 + x52*x54 - x52*x58 + x55*x57 - x55*x59;
    double x61 = 1.0/x60;
    double *x62 = U;
    double x63 = x62[22];
    double x64 = x62[19];
    double x65 = -x3*x64;
    double x66 = x62[25];
    double x67 = x62[28];
    double x68 = x13*x62[1];
    double x69 = x62[13];
    double x70 = -x0*x69 + x68;
    double x71 = x61*(x0*x66 + x3*x67 + x39*x62[10] + x40*x63 + x65 + x70);
    double x72 = x54 - x58;
    double x73 = -x5*x63;
    double x74 = x62[16];
    double x75 = x61*(x0*x74 + x18*x62[7] + x20*x64 + x5*x67 + x70 + x73);
    double x76 = x57 - x59;
    double x77 = x61*(x1*x62[4] + x3*x74 + x5*x66 + x65 + x68 + x69*x7 + x73);
    double x78 = x38*x71 + x72*x75 + x76*x77;
    double x79 = -x29*x55 + x33*x52;
    double x80 = x62[23];
    double x81 = x62[20];
    double x82 = -x3*x81;
    double x83 = x62[26];
    double x84 = x62[29];
    double x85 = x13*x62[2];
    double x86 = x62[14];
    double x87 = -x0*x86 + x85;
    double x88 = x61*(x0*x83 + x3*x84 + x39*x62[11] + x40*x80 + x82 + x87);
    double x89 = -x33*x49 + x55*x56;
    double x90 = -x5*x80;
    double x91 = x62[17];
    double x92 = x61*(x0*x91 + x18*x62[8] + x20*x81 + x5*x84 + x87 + x90);
    double x93 = x29*x49 - x52*x56;
    double x94 = x61*(x1*x62[5] + x3*x91 + x5*x83 + x7*x86 + x82 + x85 + x90);
    double x95 = x79*x88 + x89*x92 + x93*x94;
    double x96 = x71*x79 + x75*x89 + x77*x93;
    double x97 = x38*x88 + x72*x92 + x76*x94 + 1.0;
    double x98 = x78*x95 - x96*x97;
    double x99 = x62[21];
    double x100 = x62[18];
    double x101 = -x100*x3;
    double x102 = x62[24];
    double x103 = x62[27];
    double x104 = x13*x62[0];
    double x105 = x62[12];
    double x106 = -x0*x105 + x104;
    double x107 = x61*(x0*x102 + x101 + x103*x3 + x106 + x39*x62[9] + x40*x99);
    double x108 = -x5*x99;
    double x109 = x62[15];
    double x110 = x61*(x0*x109 + x100*x20 + x103*x5 + x106 + x108 + x18*x62[6]);
    double x111 = x61*(x1*x62[3] + x101 + x102*x5 + x104 + x105*x7 + x108 + x109*x3);
    double x112 = x107*x79 + x110*x89 + x111*x93 + 1.0;
    double x113 = x107*x38 + x110*x72 + x111*x76;
    double x114 = x112*x97 - x113*x95;
    double x115 = -x112*x78 + x113*x96;
    double x116 = S11*x98 + S12*x114 + S13*x115;
    double x117 = S12*x98 + S22*x114 + S23*x115;
    double x118 = S13*x98 + S23*x114 + S33*x115;
    double x119 = -x112*x116 - x117*x96 - x118*x95;
    double x120 = x17*x49 - x53*x55;
    double x121 = x13*x61;
    double x122 = -x17*x52 + x36*x55;
    double x123 = -x36*x49 + x52*x53;
    double x124 = x120*x121 + x121*x122 + x121*x123;
    double x125 = x120*x92 + x122*x88 + x123*x94;
    double x126 = x120*x75 + x122*x71 + x123*x77 + 1.0;
    double x127 = x125*x96 - x126*x95;
    double x128 = x107*x122 + x110*x120 + x111*x123;
    double x129 = -x112*x125 + x128*x95;
    double x130 = x112*x126 - x128*x96;
    double x131 = S11*x127 + S12*x129 + S13*x130;
    double x132 = S12*x127 + S22*x129 + S23*x130;
    double x133 = S13*x127 + S23*x129 + S33*x130;
    double x134 = -x112*x131 - x132*x96 - x133*x95;
    double x135 = x121*x38 + x121*x72 + x121*x76;
    double x136 = x121*x79 + x121*x89 + x121*x93;
    double x137 = -x125*x78 + x126*x97;
    double x138 = x113*x125 - x128*x97;
    double x139 = -x113*x126 + x128*x78;
    double x140 = S11*x137 + S12*x138 + S13*x139;
    double x141 = S12*x137 + S22*x138 + S23*x139;
    double x142 = S13*x137 + S23*x138 + S33*x139;
    double x143 = 1.0*PENER + 1.0*SENER;
    double x144 = -x112*x140 - x141*x96 - x142*x95 + x143;
    double x145 = -x125*x133 - x126*x132 - x128*x131;
    double x146 = -x125*x142 - x126*x141 - x128*x140;
    double x147 = -x116*x128 - x117*x126 - x118*x125 + x143;
    double x148 = -x113*x116 - x117*x78 - x118*x97;
    double x149 = -x113*x140 - x141*x78 - x142*x97;
    double x150 = -x113*x131 - x132*x78 - x133*x97 + x143;
    double x151 = x1*x61;
    double x152 = x151*x76;
    double x153 = x151*x93;
    double x154 = x123*x151;
    double x155 = x18*x61;
    double x156 = x120*x155;
    double x157 = x155*x72;
    double x158 = x155*x89;
    double x159 = x39*x61;
    double x160 = x159*x38;
    double x161 = x159*x79;
    double x162 = x122*x159;
    double x163 = x61*x7;
    double x164 = x0*x61;
    double x165 = x120*x164;
    double x166 = x122*x164;
    double x167 = x123*x163 - x165 - x166;
    double x168 = x164*x38;
    double x169 = x164*x72;
    double x170 = x163*x76 - x168 - x169;
    double x171 = x164*x79;
    double x172 = x164*x89;
    double x173 = x163*x93 - x171 - x172;
    double x174 = x3*x61;
    double x175 = x123*x174;
    double x176 = x165 + x175;
    double x177 = x174*x76;
    double x178 = x169 + x177;
    double x179 = x174*x93;
    double x180 = x172 + x179;
    double x181 = x20*x61;
    double x182 = x122*x174;
    double x183 = x120*x181 - x175 - x182;
    double x184 = x174*x38;
    double x185 = -x177 + x181*x72 - x184;
    double x186 = x174*x79;
    double x187 = -x179 + x181*x89 - x186;
    double x188 = x40*x61;
    double x189 = x5*x61;
    double x190 = x120*x189;
    double x191 = x123*x189;
    double x192 = x122*x188 - x190 - x191;
    double x193 = x189*x72;
    double x194 = x189*x76;
    double x195 = x188*x38 - x193 - x194;
    double x196 = x189*x89;
    double x197 = x189*x93;
    double x198 = x188*x79 - x196 - x197;
    double x199 = x168 + x194;
    double x200 = x171 + x197;
    double x201 = x166 + x191;
    double x202 = x184 + x193;
    double x203 = x186 + x196;
    double x204 = x182 + x190;
    
    res_0[0] = x60*(x119*x124 + x134*x135 + x136*x144);
    res_0[1] = x60*(x124*x147 + x135*x145 + x136*x146);
    res_0[2] = x60*(x124*x148 + x135*x150 + x136*x149);
    res_0[3] = x60*(x119*x154 + x134*x152 + x144*x153);
    res_0[4] = x60*(x145*x152 + x146*x153 + x147*x154);
    res_0[5] = x60*(x148*x154 + x149*x153 + x150*x152);
    res_0[6] = x60*(x119*x156 + x134*x157 + x144*x158);
    res_0[7] = x60*(x145*x157 + x146*x158 + x147*x156);
    res_0[8] = x60*(x148*x156 + x149*x158 + x150*x157);
    res_0[9] = x60*(x119*x162 + x134*x160 + x144*x161);
    res_0[10] = x60*(x145*x160 + x146*x161 + x147*x162);
    res_0[11] = x60*(x148*x162 + x149*x161 + x150*x160);
    res_0[12] = x60*(x119*x167 + x134*x170 + x144*x173);
    res_0[13] = x60*(x145*x170 + x146*x173 + x147*x167);
    res_0[14] = x60*(x148*x167 + x149*x173 + x150*x170);
    res_0[15] = x60*(x119*x176 + x134*x178 + x144*x180);
    res_0[16] = x60*(x145*x178 + x146*x180 + x147*x176);
    res_0[17] = x60*(x148*x176 + x149*x180 + x150*x178);
    res_0[18] = x60*(x119*x183 + x134*x185 + x144*x187);
    res_0[19] = x60*(x145*x185 + x146*x187 + x147*x183);
    res_0[20] = x60*(x148*x183 + x149*x187 + x150*x185);
    res_0[21] = x60*(x119*x192 + x134*x195 + x144*x198);
    res_0[22] = x60*(x145*x195 + x146*x198 + x147*x192);
    res_0[23] = x60*(x148*x192 + x149*x198 + x150*x195);
    res_0[24] = x60*(x119*x201 + x134*x199 + x144*x200);
    res_0[25] = x60*(x145*x199 + x146*x200 + x147*x201);
    res_0[26] = x60*(x148*x201 + x149*x200 + x150*x199);
    res_0[27] = x60*(x119*x204 + x134*x202 + x144*x203);
    res_0[28] = x60*(x145*x202 + x146*x203 + x147*x204);
    res_0[29] = x60*(x148*x204 + x149*x203 + x150*x202);
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
    double *x2 = coord;
    double x3 = 4*s;
    double x4 = -x3;
    double x5 = 4*t;
    double x6 = 4 - x5;
    double x7 = -8*r + x4 + x6;
    double x8 = x2[12];
    double x9 = x2[21];
    double x10 = -x5*x9;
    double x11 = x2[15];
    double x12 = x2[24];
    double x13 = x0 + x3 + x5 - 3;
    double x14 = x13*x2[0];
    double x15 = x2[18];
    double x16 = x14 - x15*x3;
    double x17 = x1*x2[3] + x10 + x11*x3 + x12*x5 + x16 + x7*x8;
    double x18 = x3 - 1;
    double x19 = -x0;
    double x20 = -8*s + x19 + x6;
    double x21 = x2[19];
    double x22 = x2[13];
    double x23 = -x0*x22;
    double x24 = x2[16];
    double x25 = x2[28];
    double x26 = x13*x2[1];
    double x27 = x2[22];
    double x28 = x26 - x27*x5;
    double x29 = x0*x24 + x18*x2[7] + x20*x21 + x23 + x25*x5 + x28;
    double x30 = x17*x29;
    double x31 = -x21*x3;
    double x32 = x2[25];
    double x33 = x1*x2[4] + x22*x7 + x24*x3 + x28 + x31 + x32*x5;
    double x34 = -x0*x8;
    double x35 = x2[27];
    double x36 = x0*x11 + x10 + x14 + x15*x20 + x18*x2[6] + x34 + x35*x5;
    double x37 = x33*x36;
    double x38 = x30 - x37;
    double x39 = x5 - 1;
    double x40 = -8*t + x19 + x4 + 4;
    double x41 = x2[23];
    double x42 = x2[20];
    double x43 = -x3*x42;
    double x44 = x2[26];
    double x45 = x2[29];
    double x46 = x13*x2[2];
    double x47 = x2[14];
    double x48 = -x0*x47 + x46;
    double x49 = x0*x44 + x3*x45 + x39*x2[11] + x40*x41 + x43 + x48;
    double x50 = -x41*x5;
    double x51 = x2[17];
    double x52 = x0*x51 + x18*x2[8] + x20*x42 + x45*x5 + x48 + x50;
    double x53 = x0*x12 + x16 + x3*x35 + x34 + x39*x2[9] + x40*x9;
    double x54 = x33*x53;
    double x55 = x1*x2[5] + x3*x51 + x43 + x44*x5 + x46 + x47*x7 + x50;
    double x56 = x0*x32 + x23 + x25*x3 + x26 + x27*x40 + x31 + x39*x2[10];
    double x57 = x36*x56;
    double x58 = x17*x56;
    double x59 = x29*x53;
    double x60 = x30*x49 - x37*x49 + x52*x54 - x52*x58 + x55*x57 - x55*x59;
    double x61 = 1.0/x60;
    double *x62 = U;
    double x63 = x62[22];
    double x64 = x62[19];
    double x65 = -x3*x64;
    double x66 = x62[25];
    double x67 = x62[28];
    double x68 = x13*x62[1];
    double x69 = x62[13];
    double x70 = -x0*x69 + x68;
    double x71 = x61*(x0*x66 + x3*x67 + x39*x62[10] + x40*x63 + x65 + x70);
    double x72 = x54 - x58;
    double x73 = -x5*x63;
    double x74 = x62[16];
    double x75 = x61*(x0*x74 + x18*x62[7] + x20*x64 + x5*x67 + x70 + x73);
    double x76 = x57 - x59;
    double x77 = x61*(x1*x62[4] + x3*x74 + x5*x66 + x65 + x68 + x69*x7 + x73);
    double x78 = x38*x71 + x72*x75 + x76*x77;
    double x79 = -x29*x55 + x33*x52;
    double x80 = x62[23];
    double x81 = x62[20];
    double x82 = -x3*x81;
    double x83 = x62[26];
    double x84 = x62[29];
    double x85 = x13*x62[2];
    double x86 = x62[14];
    double x87 = -x0*x86 + x85;
    double x88 = x61*(x0*x83 + x3*x84 + x39*x62[11] + x40*x80 + x82 + x87);
    double x89 = -x33*x49 + x55*x56;
    double x90 = -x5*x80;
    double x91 = x62[17];
    double x92 = x61*(x0*x91 + x18*x62[8] + x20*x81 + x5*x84 + x87 + x90);
    double x93 = x29*x49 - x52*x56;
    double x94 = x61*(x1*x62[5] + x3*x91 + x5*x83 + x7*x86 + x82 + x85 + x90);
    double x95 = x79*x88 + x89*x92 + x93*x94;
    double x96 = x71*x79 + x75*x89 + x77*x93;
    double x97 = x38*x88 + x72*x92 + x76*x94;
    double x98 = x97 + 1.0;
    double x99 = x78*x95 - x96*x98;
    double x100 = x62[21];
    double x101 = x62[18];
    double x102 = -x101*x3;
    double x103 = x62[24];
    double x104 = x62[27];
    double x105 = x13*x62[0];
    double x106 = x62[12];
    double x107 = -x0*x106 + x105;
    double x108 = x61*(x0*x103 + x100*x40 + x102 + x104*x3 + x107 + x39*x62[9]);
    double x109 = -x100*x5;
    double x110 = x62[15];
    double x111 = x61*(x0*x110 + x101*x20 + x104*x5 + x107 + x109 + x18*x62[6]);
    double x112 = x61*(x1*x62[3] + x102 + x103*x5 + x105 + x106*x7 + x109 + x110*x3);
    double x113 = x108*x79 + x111*x89 + x112*x93;
    double x114 = x113 + 1.0;
    double x115 = x108*x38 + x111*x72 + x112*x76;
    double x116 = x114*x98 - x115*x95;
    double x117 = -x114*x78 + x115*x96;
    double x118 = S11*x99 + S12*x116 + S13*x117;
    double x119 = S12*x99 + S22*x116 + S23*x117;
    double x120 = S13*x99 + S23*x116 + S33*x117;
    double x121 = -x113*x118 - x119*x96 - x120*x95;
    double x122 = x17*x49 - x53*x55;
    double x123 = x13*x61;
    double x124 = -x17*x52 + x36*x55;
    double x125 = -x36*x49 + x52*x53;
    double x126 = x122*x123 + x123*x124 + x123*x125;
    double x127 = x122*x92 + x124*x88 + x125*x94;
    double x128 = x122*x75 + x124*x71 + x125*x77;
    double x129 = x128 + 1.0;
    double x130 = x127*x96 - x129*x95;
    double x131 = x108*x124 + x111*x122 + x112*x125;
    double x132 = -x114*x127 + x131*x95;
    double x133 = x114*x129 - x131*x96;
    double x134 = S11*x130 + S12*x132 + S13*x133;
    double x135 = S12*x130 + S22*x132 + S23*x133;
    double x136 = S13*x130 + S23*x132 + S33*x133;
    double x137 = -x113*x134 - x135*x96 - x136*x95;
    double x138 = x123*x38 + x123*x72 + x123*x76;
    double x139 = x123*x79 + x123*x89 + x123*x93;
    double x140 = -x127*x78 + x129*x98;
    double x141 = x115*x127 - x131*x98;
    double x142 = -x115*x129 + x131*x78;
    double x143 = S11*x140 + S12*x141 + S13*x142;
    double x144 = S12*x140 + S22*x141 + S23*x142;
    double x145 = S13*x140 + S23*x141 + S33*x142;
    double x146 = 1.0*PENER + 1.0*SENER;
    double x147 = -x113*x143 - x144*x96 - x145*x95 + x146;
    double x148 = -x127*x136 - x128*x135 - x131*x134;
    double x149 = -x127*x145 - x128*x144 - x131*x143;
    double x150 = -x118*x131 - x119*x128 - x120*x127 + x146;
    double x151 = -x115*x118 - x119*x78 - x120*x97;
    double x152 = -x115*x143 - x144*x78 - x145*x97;
    double x153 = -x115*x134 - x135*x78 - x136*x97 + x146;
    double x154 = x1*x61;
    double x155 = x154*x76;
    double x156 = x154*x93;
    double x157 = x125*x154;
    double x158 = x18*x61;
    double x159 = x122*x158;
    double x160 = x158*x72;
    double x161 = x158*x89;
    double x162 = x39*x61;
    double x163 = x162*x38;
    double x164 = x162*x79;
    double x165 = x124*x162;
    double x166 = x61*x7;
    double x167 = x0*x61;
    double x168 = x122*x167;
    double x169 = x124*x167;
    double x170 = x125*x166 - x168 - x169;
    double x171 = x167*x38;
    double x172 = x167*x72;
    double x173 = x166*x76 - x171 - x172;
    double x174 = x167*x79;
    double x175 = x167*x89;
    double x176 = x166*x93 - x174 - x175;
    double x177 = x3*x61;
    double x178 = x125*x177;
    double x179 = x168 + x178;
    double x180 = x177*x76;
    double x181 = x172 + x180;
    double x182 = x177*x93;
    double x183 = x175 + x182;
    double x184 = x20*x61;
    double x185 = x124*x177;
    double x186 = x122*x184 - x178 - x185;
    double x187 = x177*x38;
    double x188 = -x180 + x184*x72 - x187;
    double x189 = x177*x79;
    double x190 = -x182 + x184*x89 - x189;
    double x191 = x40*x61;
    double x192 = x5*x61;
    double x193 = x122*x192;
    double x194 = x125*x192;
    double x195 = x124*x191 - x193 - x194;
    double x196 = x192*x72;
    double x197 = x192*x76;
    double x198 = x191*x38 - x196 - x197;
    double x199 = x192*x89;
    double x200 = x192*x93;
    double x201 = x191*x79 - x199 - x200;
    double x202 = x171 + x197;
    double x203 = x174 + x200;
    double x204 = x169 + x194;
    double x205 = x187 + x196;
    double x206 = x189 + x199;
    double x207 = x185 + x193;
    
    res_0[0] = x60*(x121*x126 + x137*x138 + x139*x147);
    res_0[1] = x60*(x126*x150 + x138*x148 + x139*x149);
    res_0[2] = x60*(x126*x151 + x138*x153 + x139*x152);
    res_0[3] = x60*(x121*x157 + x137*x155 + x147*x156);
    res_0[4] = x60*(x148*x155 + x149*x156 + x150*x157);
    res_0[5] = x60*(x151*x157 + x152*x156 + x153*x155);
    res_0[6] = x60*(x121*x159 + x137*x160 + x147*x161);
    res_0[7] = x60*(x148*x160 + x149*x161 + x150*x159);
    res_0[8] = x60*(x151*x159 + x152*x161 + x153*x160);
    res_0[9] = x60*(x121*x165 + x137*x163 + x147*x164);
    res_0[10] = x60*(x148*x163 + x149*x164 + x150*x165);
    res_0[11] = x60*(x151*x165 + x152*x164 + x153*x163);
    res_0[12] = x60*(x121*x170 + x137*x173 + x147*x176);
    res_0[13] = x60*(x148*x173 + x149*x176 + x150*x170);
    res_0[14] = x60*(x151*x170 + x152*x176 + x153*x173);
    res_0[15] = x60*(x121*x179 + x137*x181 + x147*x183);
    res_0[16] = x60*(x148*x181 + x149*x183 + x150*x179);
    res_0[17] = x60*(x151*x179 + x152*x183 + x153*x181);
    res_0[18] = x60*(x121*x186 + x137*x188 + x147*x190);
    res_0[19] = x60*(x148*x188 + x149*x190 + x150*x186);
    res_0[20] = x60*(x151*x186 + x152*x190 + x153*x188);
    res_0[21] = x60*(x121*x195 + x137*x198 + x147*x201);
    res_0[22] = x60*(x148*x198 + x149*x201 + x150*x195);
    res_0[23] = x60*(x151*x195 + x152*x201 + x153*x198);
    res_0[24] = x60*(x121*x204 + x137*x202 + x147*x203);
    res_0[25] = x60*(x148*x202 + x149*x203 + x150*x204);
    res_0[26] = x60*(x151*x204 + x152*x203 + x153*x202);
    res_0[27] = x60*(x121*x207 + x137*x205 + x147*x206);
    res_0[28] = x60*(x148*x205 + x149*x206 + x150*x207);
    res_0[29] = x60*(x151*x207 + x152*x206 + x153*x205);
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

/******************************************************************************

Welcome to GDB Online.
GDB online is an online compiler and debugger tool for C, C++, Python, Java, PHP, Ruby, Perl,
C#, OCaml, VB, Swift, Pascal, Fortran, Haskell, Objective-C, Assembly, HTML, CSS, JS, SQLite, Prolog.
Code, Compile, Run and Debug online from anywhere in world.

*******************************************************************************/
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/SparseLU>
#include <Eigen/Sparse>
#include <set>
#include <numeric>
using namespace std;
//const int n = 6; // The size of the matrix

struct Problem {
	
	int n, m;
	std::vector<double> x0, xmin, xmax;
	double x_density_percell;
	double *cell_array;
	double *solver_vec;
	double *diff_solver_vec;
	int dofs_all_num;
	int dofs_per_cell;
	int tol_dimen;
	double norm_value;
	double x_density_coff;
	double norm2=0.0;
	int num_cells=0;
	double obj_c;
	double *dc_vector;
	double *dv_vector;
	Eigen::MatrixXd dense_matrix ;
        Eigen::MatrixXd inv_dense_matrix ;
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver;
        Eigen::MatrixXd  dev_cell_matrix;
	Problem()
		: n(1)
		, m(2)
		, x0(1,0.1e6)
		, xmin(1,0.1e6)
		, xmax(1,1e6)
	{ 
                // for(int i;i<x0.size();i++){

		// }
	
	}
  void Obj(double *x, double *f0x, double *fx) {

	int nelx=16;
	int nely=8;
	int nelz=8;

	std::vector<int> iH, jH;
        int rmin=2;
	//int nelx=solid_3d.nelx;
	//int nely=solid_3d.nely;
	//int nelz=solid_3d.nelz;
	int nele=nelx*nely*nelz;
	//tol_dimen=nele;
	int iH_size = nele * (2 * (std::ceil(rmin) - 1) + 1) * (2 * (std::ceil(rmin) - 1) + 1);
//        int k = 0;
    // Define the matrix A and nu (adjust values as needed)
    std::vector<std::vector<double>> A = {
        {32, 6, -8, 6, -6, 4, 3, -6, -10, 3, -3, -3, -4, -8},
        {-48, 0, 0, -24, 24, 0, 0, 0, 12, -12, 0, 12, 12, 12}
    };
    double nu = 0.3; // Update nu with the desired value





    Eigen::MatrixXd K1(n, n);
    Eigen::MatrixXd K2(n, n);
    Eigen::MatrixXd K3(n, n);
    Eigen::MatrixXd K4(n, n);

    // You need to provide the actual content of these matrices here.
    // For example:
    // K1 << 1, 2, 3, 4, 5, 6,
    //       7, 8, 9, 10, 11, 12,
    //       ... and so on

    // Define the matrices K5 and K6
    Eigen::MatrixXd K5(n, n);
    Eigen::MatrixXd K6(n, n);



    // Compute the k vector
    std::vector<double> k(A[0].size());
    for (int i = 0; i < n; ++i) {
        k[i] = 1.0 / 144.0 * (A[0][i] + nu * A[1][i]);
    }

    // Define and compute K matrices (K1, K2, K3, K4, K5, K6)


    K1<<k[0], k[1], k[1], k[2], k[4], k[4],
        k[1], k[0], k[1], k[3], k[5], k[6],
	k[1], k[1], k[0], k[3], k[6], k[5],
	k[2], k[3], k[3], k[0], k[7], k[7],
	k[4], k[5], k[6], k[7], k[0], k[1],
	k[4], k[6], k[5], k[7], k[1], k[0];



	K2<<k[8],  k[7],  k[11], k[5],  k[3],  k[6],
	    k[7],  k[8],  k[11], k[4],  k[2],  k[5],
	    k[9],  k[9],  k[12], k[6],  k[3],  k[4],
	    k[5],  k[4],  k[10], k[8],  k[1],  k[9],
	    k[3],  k[2],  k[4],  k[1],  k[8],  k[11],
	    k[10], k[3],  k[4],  k[11], k[9],  k[12];

    K3<<k[5],  k[6],  k[3],  k[8],  k[11], k[7],
        k[6],  k[5],  k[3],  k[9],  k[12], k[9],
	k[4],  k[4],  k[2],  k[7],  k[11], k[8],
	k[8],  k[9],  k[1],  k[5],  k[10], k[4],
	k[11], k[12], k[10], k[11], k[6],  k[3],
	k[2],  k[12], k[9],  k[4],  k[5],  k[3];

   K4<<k[13], k[10], k[10], k[12], k[9],  k[9],
       k[10], k[13], k[10], k[11], k[8],  k[7],
       k[10], k[10], k[13], k[11], k[7],  k[8],
       k[12], k[11], k[11], k[13], k[6],  k[6],
       k[9],  k[8],  k[7],  k[6],  k[13], k[10],
       k[9],  k[7],  k[8],  k[6],  k[10], k[13];



K5<<k[0], k[1],  k[7],  k[2], k[4],  k[3],
    k[1], k[0],  k[7],  k[3], k[5],  k[10],
    k[7], k[7],  k[0],  k[4], k[10], k[5],
    k[2], k[3],  k[4],  k[0], k[7],  k[1],
    k[4], k[5],  k[10], k[7], k[0],  k[7],
    k[3], k[10], k[5],  k[1], k[7],  k[0];



    K6<<k[13], k[10], k[6],  k[12], k[9],  k[11],
        k[10], k[13], k[6],  k[11], k[8],  k[2],
	k[6],  k[6],  k[13], k[10], k[2],  k[9],
	k[12], k[11], k[10], k[13], k[7],  k[10],
	k[9],  k[8],  k[2],  k[7],  k[13], k[6],
	k[11], k[2],  k[9],  k[10], k[6],  k[13];

Eigen::MatrixXd KE(K1.rows()+K2.rows()+K3.rows()+K4.rows(),K1.cols()+K2.cols()+K3.cols()+K4.cols()); 
	
	KE<< K1, K2, K3, K4,
          K2.transpose(), K5, K6, K3.transpose(),
          K3.transpose(), K6, K5.transpose(), K2.transpose(),
          K4, K3, K2, K1;

    std::vector<std::vector<std::vector<int>>> il(1, std::vector<std::vector<int>>(1, std::vector<int>(nelz + 1)));
    std::vector<std::vector<std::vector<int>>> jl(1, std::vector<std::vector<int>>(1, std::vector<int>(nelz + 1)));
    std::vector<std::vector<std::vector<int>>> kl(1, std::vector<std::vector<int>>(1, std::vector<int>(nelz + 1)));
    std::vector<std::vector<std::vector<int>>> loadnid(1, std::vector<std::vector<int>>(1, std::vector<int>(nelz + 1)));
    // Fill the 1x1x9 coordinate grids using nested loops
    Eigen::VectorXd loaddof(nelz+1);
    for (int k = 0; k <= nelz; ++k) {
        il[0][0][k] = nelx;
        jl[0][0][k] = 0;
        kl[0][0][k] = k;
        loadnid[0][0][k]=kl[0][0][k]*(nelx+1)*(nely+1)+il[0][0][k]*(nely+1)+(nely+1-jl[0][0][k]);     
	loaddof[k]=3*loadnid[0][0][k]-1;
    }


    std::vector<std::vector<std::vector<int>>> iif(nely + 1, std::vector<std::vector<int>>(1, std::vector<int>(nelz + 1)));
    std::vector<std::vector<std::vector<int>>> jf(nely + 1, std::vector<std::vector<int>>(1, std::vector<int>(nelz + 1)));
    std::vector<std::vector<std::vector<int>>> kf(nely + 1, std::vector<std::vector<int>>(1, std::vector<int>(nelz + 1)));
    std::vector<std::vector<std::vector<int>>> fixednid(nely + 1, std::vector<std::vector<int>>(1, std::vector<int>(nelz + 1))); 


    Eigen::VectorXd fixednid_1((nelz+1)*(nely+1));
    Eigen::VectorXd fixednid_2((nelz+1)*(nely+1));
    Eigen::VectorXd fixednid_3((nelz+1)*(nely+1));
    Eigen::VectorXd fixeddof(3*(nelz+1)*(nely+1));
    int fix_counter=0;
    #pragma omp parallel for sum(+:fix_counter)
    for (int k = 0; k <= nelz; ++k) {
        for (int j = 0; j <= nely; ++j) {
            iif[j][0][k] = 0;
            jf[j][0][k] = j;
            kf[j][0][k] = k;
	    fixednid[j][0][k]=kf[j][0][k]*(nelx+1)*(nely+1)+iif[j][0][k]*(nely+1)+(nely+1-jf[j][0][k]);
            fixednid_1[fix_counter]=3*fixednid[j][0][k];
            fixednid_2[fix_counter]=3*fixednid[j][0][k]-1;
	    fixednid_3[fix_counter]=3*fixednid[j][0][k]-2;
            fix_counter++;

        }
    }


fixeddof<<fixednid_1,fixednid_2,fixednid_3;

Eigen::VectorXi fixeddof_xi(3*(nelz+1)*(nely+1));
fixeddof_xi<<fixednid_1.cast<int>(),fixednid_2.cast<int>(),fixednid_3.cast<int>();
//int nele = nelx*nely*nelz;
int ndof = 3*(nelx+1)*(nely+1)*(nelz+1);


    //Eigen::VectorXi loaddof(243); // Example data with 243 elements

    // Assuming ndof is the total number of degrees of freedom in the finite element model

    // Create the sparse matrix F
    Eigen::SparseMatrix<int> F(ndof, 1);
    F.reserve(loaddof.size());
    for (int i = 0; i < loaddof.size(); ++i) {
        F.insert(loaddof(i), 0) = -1;
    }
    F.makeCompressed();

    // Access and print the result (optional)
  //  std::cout << F << std::endl;
  Eigen::VectorXd U = Eigen::VectorXd::Zero(ndof);

  Eigen::VectorXi freedofs(ndof);
  freedofs = Eigen::VectorXi::LinSpaced(ndof, 1, ndof);
  std::set<int> set1(freedofs.data(), freedofs.data() + freedofs.size());
  std::set<int> set2(fixeddof_xi.data(), fixeddof_xi.data() + fixeddof_xi.size());

  std::vector<int>freedof_result;
  std::set_difference(set1.begin(), set1.end(), set2.begin(), set2.end(), std::back_inserter(freedof_result));
  std::cout<<freedof_result.size()<<std::endl;
  Eigen::VectorXi freedof_xi(freedof_result.size());
  for (int i = 0; i < freedof_result.size(); ++i) {
        freedof_xi(i) = freedof_result[i];
    }
//  std::cout<<freedof_xi.transpose()<<std::endl;

//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
/*std::vector<std::vector<int>> nodegrd(nely + 1, std::vector<int>(nelx + 1));

    int counter = 1;
    #pragma omp for
    for (int i = 0; i <= nelx; i++) {
        for (int j = 0; j <= nely; j++) {
            nodegrd[j][i] = counter;
            counter++;
        }
     }
       

   std::vector<int> nodeids(nelx * nely);
   Eigen::MatrixXi eigen_nodeids(nely * nelx, 1);
    int index = 0;
    #pragma omp for
    for (int i = 0; i < nelx; i++) {
        for (int j = 0; j < nely; j++) {
            nodeids[index] = nodegrd[j][i];
            eigen_nodeids(index)=nodeids[index];
            index++;
        }
    }

    std::vector<int> nodeidz;

    int step = (nely + 1) * (nelx + 1);
     for (int i = 0; i <= (nelz - 1) * step; i += step) {
        nodeidz.push_back(i);

    }//
  
Eigen::MatrixXi eigen_nodeidz(nodeidz.size(),1);

   #pragma omp for 
   for (int i = 0; i < nodeidz.size(); i++) {
       eigen_nodeidz(i)=nodeidz[i];

    }

x_density_percell=x[0];
std::vector<std::vector<std::vector<double>>> xPhys(nely, std::vector<std::vector<double>>(nelx, std::vector<double>(nelz)));
std::cout<<"-------------------mu is:"<<x_density_percell<<"------------------"<<std::endl;
  #pragma omp for
  for (int i = 0; i < nely; ++i) {
        for (int j = 0; j < nelx; ++j) {
            for (int k = 0; k < nelz; ++k) {
              xPhys[i][j][k]=x_density_percell;
            }
        }
    }

*/


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
/*int eigen_nodeidz_cols = eigen_nodeidz.transpose().cols();
    Eigen::MatrixXi replicated_nodeids = eigen_nodeids.replicate(1,eigen_nodeidz_cols);

    // Column replication of nodeidz
    int eigen_nodeids_rows = eigen_nodeids.rows();
    Eigen::MatrixXi replicated_nodeidz = eigen_nodeidz.transpose().replicate(eigen_nodeids_rows,1);
    // Element-wise addition of replicated matrices
   Eigen::MatrixXi result = replicated_nodeids + replicated_nodeidz;

   Eigen::Map<Eigen::VectorXi>result_vec(result.data(), result.size());


   Eigen::VectorXi edofVec(result_vec.rows());
   #pragma omp for
   for(int i=0;i<result_vec.rows();i++){

      edofVec(i) = result_vec(i)*3+1;
      //result_vec(i) = edofVec(i);
   }
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   Eigen::MatrixXi edofMat1(nele,dofs_per_cell);
   edofMat1 = edofVec.replicate(1, dofs_per_cell);
   //std::cout<<edofMat1<<std::endl;
    Eigen::VectorXd ele1(3);
    ele1 << 0, 1, 2;

    Eigen::VectorXd ele2(6);
    ele2 << 3, 4, 5, 0, 1, 2;
    ele2 = ele2.array() + 3 * nely; // Add 3*nely to each element of ele2

    Eigen::VectorXd ele3(3);
    ele3 << -3, -2, -1;

    // Concatenate ele1, ele2, and ele3 into test_ele
    Eigen::VectorXd test_ele(ele1.size() + ele2.size() + ele3.size());
    test_ele.block(0, 0, ele1.size(), 1) = ele1;
    test_ele.block(ele1.size(), 0, ele2.size(), 1) = ele2;
    test_ele.block(ele1.size() + ele2.size(), 0, ele3.size(), 1) = ele3;


    Eigen::VectorXd test_ele1(test_ele.size()*2);
    Eigen::VectorXd temp_ele(test_ele.size());

    temp_ele=3*(nely+1)*(nelx+1)+test_ele.array();

   // cout<<"-------------------"<<endl;
    test_ele1.block(0,0,test_ele.size(),1)=test_ele;
    test_ele1.block(test_ele.size(),0, temp_ele.size(),1)=temp_ele;
    //std::cout << "test_ele: \n" << test_ele1.transpose() << std::endl;
    Eigen::MatrixXi edofMat2(nele,dofs_per_cell);
    Eigen::VectorXi edofVec2(test_ele1.rows());
   
    edofVec2=test_ele1.cast<int>();
    edofMat2=edofVec2.transpose().replicate(nele,1);
   // std::cout<<edofMat2<<std::endl;
   // std::cout<<"-----------------------------------"<<std::endl;

    Eigen::MatrixXi edofMat=edofMat1+edofMat2 ;
    
    
    //>>>>>>>>>>>>>
    //Eigen::VectorXi onesVector(24);
    //onesVector.setOnes();
    //int rows = edofMat.rows();
    //int cols = edofMat.cols();
//    int kronRows = rows * onesVector.rows();
  //  int kronCols = cols * onesVector.cols();
   // Eigen::MatrixXi kronResult(kronRows, kronCols);
*/
/*
     for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            kronResult.block(i * onesVector.rows(), j * onesVector.cols(), onesVector.rows(), onesVector.cols()) = edofMat(i, j) * onesVector;
        }
    }
*/
   // Eigen::VectorXi iK = Eigen::Map<Eigen::VectorXi>(kronResult.data(), kronResult.size());
    /*
    Eigen::MatrixXi kronResult = Eigen::kroneckerProduct(edofMat, onesVector).transpose();
    Eigen::VectorXi iK = Eigen::Map<Eigen::VectorXi>(kronResult.data(), kronResult.size());
    */
  //  std::cout<<iK.size()<<std::endl;


    //std::cout<<edofMat<<std::endl;
/*
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  std::set<std::pair<int, int>> added_indices;
  //Eigen::SparseMatrix<double> H(nele, nele);
  std::vector<Eigen::Triplet<double>> sH;
  for (int k1 = 1; k1 <= nelz; ++k1) {
        for (int i1 = 1; i1 <= nelx; ++i1) {
            for (int j1 = 1; j1 <= nely; ++j1) {
                int e1 = (k1 - 1) * nelx * nely + (i1 - 1) * nely + j1;
                for (int k2 = std::max(k1 - static_cast<int>(std::ceil(rmin) - 1), 1); k2 <= std::min(k1 + static_cast<int>(std::ceil(rmin) - 1), nelz); ++k2) {
                    for (int i2 = std::max(i1 - static_cast<int>(std::ceil(rmin) - 1), 1); i2 <= std::min(i1 + static_cast<int>(std::ceil(rmin) - 1), nelx); ++i2) {
                        for (int j2 = std::max(j1 - static_cast<int>(std::ceil(rmin) - 1), 1); j2 <= std::min(j1 + static_cast<int>(std::ceil(rmin) - 1), nely); ++j2) {
                            int e2 = (k2 - 1) * nelx * nely + (i2 - 1) * nely + j2;

                            // Check if (e1, e2) pair is already added
                            if (added_indices.find(std::make_pair(e1, e2)) == added_indices.end()) {
                                iH.push_back(e1);
                                jH.push_back(e2);
                                //double distance_squared = (i1 - i2) * (i1 - i2) + (j1 - j2) * (j1 - j2) + (k1 - k2) * (k1 - k2);
                                //double distance = std::sqrt(distance_squared);
                                //sH.push_back(std::max(0.0, rmin - distance));
                                sH.push_back(Eigen::Triplet<double>(e1 - 1, e2 - 1, std::max(0.0, rmin - sqrt((i1 - i2) * (i1 - i2) + (j1 - j2) * (j1 - j2) + (k1 - k2) * (k1 - k2)))));

                                added_indices.insert(std::make_pair(e1, e2)); // Add the pair to the set
                            }
                        }
                    }
                }
            }
        }
    }
//<<<<<<<<<<<<<<<<<<<<<

auto maxih_element=std::max_element(iH.begin(),iH.end());
auto maxjh_element=std::max_element(jH.begin(),jH.end());
std::cout<<"+++++++++++++++++iH max:"<<*maxih_element<<std::endl;
std::cout<<"+++++++++++++++++jH max:"<<*maxjh_element<<std::endl;
Eigen::MatrixXd H_dense(*maxih_element,*maxjh_element);
Eigen::SparseMatrix<double> H(*maxih_element,*maxjh_element);
H.setFromTriplets(sH.begin(), sH.end());
Eigen::VectorXd Hs(*maxih_element);
//std::cout<<"-------------"<<*maxih_element<<","<<*maxjh_element<<std::endl;
//std::cout<<H.rows()<<","<<H.cols<<std::endl;
//std::cout<<Hs.rows()<<std::endl;
double sum_temp=0.0;
#pragma omp for
for(int i=0;i<*maxih_element;i++){
  for(int j=0;j<*maxih_element;j++){
     
     sum_temp+=H.coeff(i, j);
     H_dense(i,j)=H.coeff(i, j);
   }

   Hs[i]=sum_temp;
   sum_temp=0.0;

//   std::cout<<"############################("<<i<<",1): "<<Hs[i]<<std::endl;
}

std::cout<<"++++++++++++++++++++"<<Hs.rows()<<std::endl;

//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 // const int n_dofs = dof_handler.n_dofs();active_cells
//    const int dofs_per_cell = fe.dofs_per_cell;dofs_per_cell
std::cout<<"+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl;

int num_cells=edofMat.rows();
int num_dofs=edofMat.cols();
Eigen::MatrixXd U_edofMat(num_cells, num_dofs);
#pragma omp for
for (int i = 0; i < num_cells; ++i) {
        for (int j = 0; j < num_dofs; ++j) {
            int dof_index = edofMat(i, j);  // Get dof index from edofMat
            if(dof_index<dofs_all_num){
	    U_edofMat(i, j) = solid_3d.solution_n[dof_index];
	    } // Extract corresponding U value
        }
    }
//std::cout << "U(edofMat):\n" << U_edofMat << std::endl;


Eigen::MatrixXd U_KE=U_edofMat*KE;


Eigen::MatrixXd KE_U_multi_U=U_KE.array()*U_edofMat.array();

 Eigen::VectorXd sumVector =sumAlongRows(KE_U_multi_U);

// std::cout<<sumVector.rows()<<std::endl;


//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


dv_vector=new double[nely*nelx*nelz]; 
 
 
int penal=3;
double E0=1.0;
double Emin=1e-9;
std::vector<std::vector<std::vector<double>>> ce(nely, std::vector<std::vector<double>>(nelx, std::vector<double>(nelz)));
//std::vector<std::vector<std::vector<double>>> xPhys(nely, std::vector<std::vector<double>>(nelx, std::vector<double>(nelz)));
std::vector<std::vector<std::vector<double>>> c_temp(nely, std::vector<std::vector<double>>(nelx, std::vector<double>(nelz)));
std::vector<std::vector<std::vector<double>>> dc(nely, std::vector<std::vector<double>>(nelx, std::vector<double>(nelz)));
std::vector<std::vector<std::vector<double>>> dv(nely, std::vector<std::vector<double>>(nelx, std::vector<double>(nelz,1)));
    int idx = 0;
    #pragma omp for
    for (int i = 0; i < nely; ++i) {
        for (int j = 0; j < nelx; ++j) {
            for (int k = 0; k < nelz; ++k) {
                ce[i][j][k] =sumVector[idx];
		//xPhys[i][j][k]=1.3;
		dc[i][j][k]=-1*penal*(E0-Emin)*pow(xPhys[i][j][k],penal-1)*ce[i][j][k];
		xPhys[i][j][k]=Emin+pow(xPhys[i][j][k],penal-1)*(E0-Emin);
                c_temp[i][j][k]=ce[i][j][k]*xPhys[i][j][k];

                ++idx;
            }
        }
    }




  std::vector<std::vector<double>> sum_nelx_nelz(nelx, std::vector<double>(nelz, 0.0));



  for (int i = 0; i < nelx; ++i) {
      for (int j = 0; j < nelz; ++j) {
         for (int k = 0; k < nely; ++k) {
             sum_nelx_nelz[i][j] += c_temp[k][i][j];
         }
     }
  }


    std::vector<double> sum_nelz(nelz, 0.0);

    for (int i = 0; i < nelz; ++i) {
        for (int k = 0; k < nelx; ++k) {
           sum_nelz[i]=sum_nelx_nelz[k][i];
        }
    }



  double c=std::accumulate(sum_nelz.begin(),sum_nelz.end(),0.0);
  std::cout<<"-------****************************************************************************--------------------------"<<std::endl;
  std::cout<<c<<std::endl;
  std::cout<<"-------****************************************************************************--------------------------"<<std::endl;
         *f0x=c;

  int dc_num=0;
  Eigen::VectorXd dc_vec(nely*nelx*nelz);
  Eigen::VectorXd dv_vec(nely*nelx*nelz);
 
   for (int k = 0; k < nelz; ++k) {
        for (int j = 0; j < nelx; ++j) {
            for (int i = 0; i < nely; ++i) {
               dc_vec(dc_num)=dc[i][j][k];
	       dv_vec(dc_num)=dv[i][j][k];
	       dc_num++;
            }
        }
    }
   
  dc_vector=new double[nely*nelx*nelz];
  dc_vec=dc_vec.array()/Hs.array();
  dc_vec=H_dense*dc_vec;
  dv_vec=dv_vec.array()/Hs.array();
  dv_vec=H_dense*dv_vec;

  #pragma omp for
  for(int i=0;i<dc_vec.rows();i++){
   dc_vector[i] = dc_vec[i];
   dv_vector[i] = dv_vec[i];
  }
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
*/
 }

    Eigen::VectorXd sumAlongRows(const Eigen::MatrixXd& matrix) {
    Eigen::VectorXd sumVector(matrix.rows());
    #pragma omp for
    for (int i = 0; i < matrix.rows(); ++i) {
        sumVector(i) = matrix.row(i).sum();
     }
      return sumVector;
     
    }

};



int main(){
double f=0.0;
Problem problem;
vector<double>x={0.1e6};
vector<double>g={0,0};
problem.Obj(x.data(),&f,g.data());


}




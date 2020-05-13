#include <fstream>
#include <stdlib.h>
#include <iostream>
#include "gurobi_c++.h"

#define BLOCK_SIZE (16)
#define ROUNDS (15)
using namespace std;


static void XORR (GRBEnv env, GRBModel& model,int num_rounds, GRBVar A[][BLOCK_SIZE], int a, GRBVar B[][BLOCK_SIZE],GRBVar C[][BLOCK_SIZE], GRBVar D[][BLOCK_SIZE])
{
	
		for (int j=0;  j<BLOCK_SIZE; j++)
		{
		model.addConstr(A[a-1][j] + B[num_rounds-1][j] + C[num_rounds-1][j]>=2*D[num_rounds-1][j]);
		model.addConstr(D[num_rounds-1][j]>=A[a-1][j]);
		model.addConstr(D[num_rounds-1][j]>=B[num_rounds-1][j]);
		model.addConstr(D[num_rounds-1][j]>=C[num_rounds-1][j]);
		model.addConstr(A[a-1][j] + B[num_rounds-1][j] + C[num_rounds-1][j]<=2);
		}

}

static void HREPRESENTATION (GRBEnv env, GRBModel& model,int num_rounds, GRBVar T1[][BLOCK_SIZE], GRBVar T3[][BLOCK_SIZE],
							 GRBVar L[][BLOCK_SIZE],int r)
{
	
		for (int j=0;  j<BLOCK_SIZE; j++)
		{
			
			model.addConstr(L[r-1][(j+0)%BLOCK_SIZE]+ L[r-1][(j+5)%BLOCK_SIZE] - T1[num_rounds-1][j]>=0);
			model.addConstr(L[r-1][(j+1)%BLOCK_SIZE]==T3[num_rounds-1][j]);
		
		}

}

static void
		inOputmask (GRBEnv env, GRBModel& model, GRBVar L[][BLOCK_SIZE],GRBVar R[][BLOCK_SIZE],GRBVar L1[][BLOCK_SIZE],
		GRBConstr m1[BLOCK_SIZE],GRBConstr m2[BLOCK_SIZE],GRBConstr c1[BLOCK_SIZE],GRBConstr c2[BLOCK_SIZE])
{			
			for (int i=0;i<BLOCK_SIZE;i++)
		{
		
			m1[i] = model.addConstr(L[0][i]==0);
	
			m2[i] = model.addConstr(R[0][i]==0);

			c1[i] = model.addConstr(L1[ROUNDS-1][i]==0);

		    c2[i] =  model.addConstr(L1[ROUNDS-2][i]==0);
		}

}


// static void
// 		inOputmask11 (GRBEnv env, GRBModel& model, GRBVar L[][BLOCK_SIZE],GRBVar R[][BLOCK_SIZE],GRBVar L1[][BLOCK_SIZE],
// 		GRBConstr m1[BLOCK_SIZE],GRBConstr m2[BLOCK_SIZE],GRBConstr c1[BLOCK_SIZE],GRBConstr c2[BLOCK_SIZE])
// {			
// 			for (int i=0;i<BLOCK_SIZE;i++)
// 		{

// 				 m1[i] = model.addConstr(L[0][i]==0);}

// 			for (int i=0;i<BLOCK_SIZE;i++){
// 			if(i==20 )
// 			m2[i] = model.addConstr(R[0][i]==1);
// 			else
// 				m2[i] = model.addConstr(R[0][i]==0);}

// 			for (int i=0;i<BLOCK_SIZE;i++)
// 				if(i==8 )
// 			c1[i] = model.addConstr(L1[ROUNDS-1][i]==0);

// 			for (int i=0;i<BLOCK_SIZE;i++){
// 				 if(i==20 )
// 		    c2[i] =  model.addConstr(L1[ROUNDS-2][i]==1);
// 				 else
//              c2[i] =  model.addConstr(L1[ROUNDS-2][i]==0);
// 		}

// }

int main(int argc, char *argv[])
{
  GRBVar *vars = 0;
	
		GRBEnv env = GRBEnv();
		GRBModel model = GRBModel(env);

		// Set variables
		GRBVar L[ROUNDS][BLOCK_SIZE] ;
		GRBVar R[ROUNDS][BLOCK_SIZE] ;
		GRBVar T1[ROUNDS][BLOCK_SIZE] ;
		GRBVar T2[ROUNDS][BLOCK_SIZE] ;
		GRBVar T3[ROUNDS][BLOCK_SIZE] ;
		GRBVar R1[ROUNDS][BLOCK_SIZE] ;
		GRBVar DT1[ROUNDS][BLOCK_SIZE] ;
		GRBVar DT3[ROUNDS][BLOCK_SIZE] ;
		GRBConstr m1[BLOCK_SIZE];
		GRBConstr m2[BLOCK_SIZE];
		GRBConstr c1[BLOCK_SIZE];
		GRBConstr c2[BLOCK_SIZE];
		
		for(int i=0;i<ROUNDS; i++)
			for(int j=0;j<BLOCK_SIZE; j++)
			{
				L[i][j]= model.addVar(0.0,1.0,0.0,GRB_BINARY);
		        R[i][j]= model.addVar(0.0,1.0,0.0,GRB_BINARY);
				T1[i][j]= model.addVar(0.0,1.0,0.0,GRB_BINARY);
				T2[i][j]= model.addVar(0.0,1.0,0.0,GRB_BINARY);
				T3[i][j]= model.addVar(0.0,1.0,0.0,GRB_BINARY);
				R1[i][j]= model.addVar(0.0,1.0,0.0,GRB_BINARY);
				DT1[i][j]= model.addVar(0.0,1.0,0.0,GRB_BINARY);
				DT3[i][j]= model.addVar(0.0,1.0,0.0,GRB_BINARY);
			}

		model.update();
		
		
		// ROUND 1
		XORR(env,model,1,R,1,T1,T2,DT1);
		XORR(env,model,1,T2,1,T3,R1,DT3);
		HREPRESENTATION(env,model,1,T1,T3,L,1);

		// ROUND 2
		if(ROUNDS!=1)
		{
		XORR(env,model,2,L,1,T1,T2,DT1);
		XORR(env,model,2,T2,2,T3,R1,DT3);
		HREPRESENTATION(env,model,2,T1,T3,R1,1);
		}
		//OTHER ROUNDS
		for(int i=3;i<=ROUNDS;i++)
		{
		XORR(env,model,i,R1,i-2,T1,T2,DT1);
		XORR(env,model,i,T2,i,T3,R1,DT3);
		HREPRESENTATION(env,model,i,T1,T3,R1,i-1);
		}

		inOputmask(env,model,L,R,R1,m1,m2,c1,c2);

		//objective function
		// GRBQuadExpr OBJfun1,OBJfun2,OBJfun;
		GRBLinExpr OBJfun1,OBJfun2,OBJfun;
		for(int j=0;j<BLOCK_SIZE;j++)
			OBJfun1 += L[0][j];
		for(int i=2;i<=ROUNDS;i++)
			for(int j=0;j<BLOCK_SIZE;j++)
			OBJfun2+=R1[i-2][j];
		OBJfun=OBJfun1+OBJfun2;
		model.addConstr(OBJfun>=1);


		model.setObjective(OBJfun,GRB_MINIMIZE);
		model.write("model.lp");

int numvars = model.get(GRB_IntAttr_NumVars);
vars = model.getVars();


		model.optimize();
		
std::cout << "----------------------\n ";


	for(int i=0;i<BLOCK_SIZE;i++){
				m2[i].set(GRB_DoubleAttr_RHS, 1.0);
				for(int j=0;j<BLOCK_SIZE;j++){
					c1[j].set(GRB_DoubleAttr_RHS, 1.0);
			
				model.update();
		
		model.optimize();
			if (model.get(GRB_IntAttr_Status) == 3){
			std::cout<< "model is infeasible"<<"\n";
					
			cout <<"m2 "<< i <<"\t"<<"c1 "<< j <<"\n";
std::cout << "---------------------\n ";
	//			}
				c1[j].set(GRB_DoubleAttr_RHS, 0.0);
				model.update();}
			
				m2[i].set(GRB_DoubleAttr_RHS, 0.0);
				model.update();
			}


	
}

	return 0;
}

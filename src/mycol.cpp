#include "mycol.h"
 
int **userCol;

void initColors(){
	userCol = new int* [6];
	for(int i=0; i<6; i++) {
		userCol[i] = new int [3];
	}
	
	userCol[0][0] = 255;
	userCol[0][1] = 255;
	userCol[0][2] = 255;
	
	userCol[1][0] = 0;
	userCol[1][1] = 0;
	userCol[1][2] = 0;
	
	userCol[2][0] = 255;
	userCol[2][1] = 0;
	userCol[2][2] = 0;
	
	userCol[3][0] = 0;
	userCol[3][1] = 0;
	userCol[3][2] = 255;
	
	userCol[4][0] = 255;
	userCol[4][1] = 255;
	userCol[4][2] = 0;
	
	userCol[5][0] = 0;
	userCol[5][1] = 255;
	userCol[5][2] = 255;
}

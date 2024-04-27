#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* 
    Dona't que operarem amb matrius
    definirem una classe que ens permeti 
    accedir rapidament a les files i les 
    columnes
*/

typedef struct matrix{
    int row;
    int col;
    double** data;
}matrix;

/*
    Necessitem un operador per inicialitzar
    les matrius donades 2 dimensions
*/

void StartMatrix2(matrix *X, int n, int m){
  X->data= (double**)calloc(n,sizeof(double*));
  int i;
  for( i = 0 ; i<n ; i++){
    X->data[i] = (double*)calloc(m,sizeof(double));
  }
  X->row = n;
  X->col = m;
  return;
}
void ReadMatrixFile( char* name, matrix *C){
  FILE *fp;
  fp = fopen(name,"r");
  if(fp == NULL) return;
  int i,j;

  fscanf(fp,"%d",&i);
  fscanf(fp,"%d",&j);
 
  StartMatrix2(C,i,j);
  
  
  // LEE la dimension, modificar por scanf si el usuario la introduce manualmente

  
  for(i=0;i<C->row;i++){
    for(j=0;j<C->col;j++){
      fscanf(fp,"%lf",&C->data[i][j]);
    }
  }

  fclose (fp);
  return;
}
/*
    Fara falta imprimir matrius
*/
void PrintMatrixFile( matrix A, char* name){
  FILE *fp;
  int i,j;
  fp = fopen(name,"w");
  if(fp== NULL) return;
  for(i=0;i<A.row;i++){
    for(j=0;j<A.col;j++){
      fprintf(fp,"%e\t",A.data[i][j]);
    }
    fprintf(fp,"\n");
  }
  fclose (fp);
  return;
}
//            Imprimir Matrius al temrinal
void PrintMatrix(matrix X){
  int i,j;
  printf("\n");
  for(i=0;i<X.row;i++){
    for(j=0;j<X.col;j++){
      printf("%0.15e\t",X.data[i][j]);
    }
    printf("\n");
  }
  return;
}

/*
    Fara falta multiplicar matrius, 
    A*B = X
    el operador assumeix que les
    dimensions son correctes
*/

void ProductMatrix( matrix A,  matrix B, matrix X){
     int i,j,k;
     matrix C;
     StartMatrix2(&C,A.row,B.col);
     // Calculem la multiplicacio
     for(i=0;i<A.row;i++){
          for(j=0;j<B.col;j++){
               for(k=0;k<A.col;k++){
                    C.data[i][j]+=A.data[i][k]*B.data[k][j];
               }
          }
     }
     // Copiem els valors a X
     for(i=0;i<A.row;i++){
          for(j=0;j<B.col;j++){
               for(k=0;k<A.col;k++){
                    X.data[i][j]=C.data[i][j];
               }
          }
     }  
    // Alliberem C
    for( i =0; i<C.row ; i++){
      free(C.data[i]);
    }
    free(C.data);
     return;
}

double norm_inf(matrix X){
    int i,j;
    double result = 0.0, aux = 0.0;

    for(i=0; i<X.row; i++){
        for(j=0; j<X.col; j++){
            aux += fabs(X.data[i][j]);
        }
        
        if(aux >= result) result = aux;
        
        aux = 0;
    }
    return result;
}


void swap_rows(matrix *A, int row1, int row2) {
    double *temp = A->data[row1];
    A->data[row1] = A->data[row2];
    A->data[row2] = temp;
}

void gauss_jordan(matrix C, matrix b, matrix x) {
    int i, j, k, max_row;
    double pivot, temp;
    
    matrix A;
    StartMatrix2(&A,C.row,C.col);
    for ( i = 0; i < C.row; i++){
       
        for ( j = 0; j < C.col; j++){
           A.data[i][j] = C.data[i][j];
        } 
    }

    
    // Copy the solution to matrix x
    for (i = 0; i < x.row; i++) {
        x.data[i][0] = b.data[i][0];
    }
    // Forward elimination with partial pivoting
    for (i = 0; i < A.row; i++) {
        max_row = i;
        for (j = i + 1; j < A.row; j++) {
            if (abs(A.data[j][i]) > abs(A.data[max_row][i])) {
                max_row = j;
            }
        }
        swap_rows(&A, i, max_row);
        swap_rows(&b, i, max_row);

        pivot = A.data[i][i];
        if(pivot == 0){
            printf("Singular matrix");
            return;
        }
        for (j = i; j < A.col; j++) {
            A.data[i][j] /= pivot;
        }
        b.data[i][0] /= pivot;

        for (k = 0; k < A.row; k++) {
            if (k != i) {
                temp = A.data[k][i];
                for (j = i; j < A.col; j++) {
                    A.data[k][j] -= temp * A.data[i][j];
                }
                b.data[k][0] -= temp * b.data[i][0];
            }
        }
    }


for( i =0; i<A.row ; i++){
    free(A.data[i]);
}

free(A.data);

    return;
}

void F(matrix z, matrix f){
    double x,y;
    x = z.data[0][0];
    y = z.data[1][0];
    double u = 0.5;
    f.data[0][0] = x+u*cos(x+y)-u*u;
    f.data[1][0] = x + u*0.5*y+u*x*x+y*y;
    return;
}
void DF(matrix z, matrix D){
    double x,y;
    x = z.data[0][0];
    y = z.data[1][0];
    double u = 0.5;
    D.data[0][0] = 1-u*sin(x+y);
    D.data[0][1] = -u*sin(x+y);
    D.data[1][0] = 1+2*u*x;
    D.data[1][1] = y*y + 0.5*u;

}

void Newton(matrix x){
    matrix A;
    StartMatrix2(&A,x.row,x.row);
    matrix f;
    StartMatrix2(&f,x.row,1);
    DF(x,A);
    F(x,f);
    int i;
    for ( i = 0; i < f.row; i++)
    {
        f.data[i][0] *= -1;
    }
    
    gauss_jordan(A,f,f);

    for ( i = 0; i <x.row; i++){
        x.data[i][0] +=  f.data[i][0];
    }


    for( i =0; i<A.row ; i++){
    free(A.data[i]);
    }
    for( i =0; i<f.row ; i++){
    free(f.data[i]);
    }

    free(A.data);

return;
}



int main() {

    matrix x;
    matrix y;
    StartMatrix2(&y,2,1);
    StartMatrix2(&x,2,1);


    int i;
    
    for ( i = 0; i < 10000; i++)
    {
          Newton(x);
          F(x,y);
          if (norm_inf(y)<10e-10)
          {
            break;
          }
    }
    PrintMatrix(x);
    printf("Tenim error: %.16g",norm_inf(y));
   

    // Print the solution
    
    // Free memory
    for (int i = 0; i < x.row; i++) {
        free(x.data[i]);
    }
    free(x.data);

    return 0;
}
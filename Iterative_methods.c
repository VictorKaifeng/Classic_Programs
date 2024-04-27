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
      printf("%lf\t",X.data[i][j]);
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

/*
    Treballarem amb la norma infinit
*/

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

/* 
    Computa un iterar del 
    metode Gauss-Seidel de Ax = b
*/
void  Gauss_Seidel(matrix A, matrix b, matrix x){
  // x_i = 1/a_ii (b_i - \sum a_ij x_j)
  int i,j;
  for(i=0; i<x.row; i++ ){
    x.data[i][0] = b.data[i][0] / A.data[i][i];
    for(j=0; (j< x.row) ; j++){
      if(j != i){ 
        x.data[i][0] -= A.data[i][j] / A.data[i][i] * x.data[j][0];
      }
    }
  }
  return;

}

/* 
    Computa un iterar del 
    metode Jacobi de Ax = b
*/
void Jacobi(matrix A, matrix b, matrix x) {
  int i,j;
  matrix  y;
  StartMatrix2(&y,x.row,x.col);
  for( i = 0 ; i< x.row ; i++){
    for( j = 0; j<x.col; j++){
      y.data[i][j] = x.data[i][j];
    }
  }
  for(i=0; i<A.col; i++ ){
    x.data[i][0] = b.data[i][0] / A.data[i][i];
    for(j=0; j< A.col; j++){
      if(j!= i){
        x.data[i][0] -= A.data[i][j] / A.data[i][i] * y.data[j][0];
      }
    }
  }
  //Alliberem y
  for( i =0; i<y.row ; i++){
    free(y.data[i]);
  }
    free(y.data);
  return ;

}

/* 
    Computa un iterar del 
    metode SOR de Ax = b
    amb par. de relaxacio w
*/
void  SOR(matrix A, matrix b, matrix x, double w){
  // x_i = x_i + w* (1/a_ii (b_i - \sum a_ij x_j)) [=: r_i]
  int i,j;
  double r = 0;
  for(i=0; i<x.row; i++ ){
    r = b.data[i][0] / A.data[i][i];
    for(j=0; j< x.row ; j++){
      r -= A.data[i][j] / A.data[i][i] * x.data[j][0];
    }
    x.data[i][0] += w*r;
  }
  return;

}

/*
    Computa el error absolut de Ax-b,
    amb la norma infinit
*/ 
double absolut_err_inf(matrix A, matrix b, matrix x){
  matrix C;
  StartMatrix2(&C,A.row,x.col);
  ProductMatrix(A,x,C);
  int i;
  for(i=0; i< C.row; i++){
    C.data[i][0] -= b.data[i][0];
  }
  double value = norm_inf(C);
  //Alliberem C
  for( i =0; i<C.row ; i++){
    free(C.data[i]);
  }
    free(C.data);
  return value; 
}

/*
  Computa la norma infinit de la diferencia de dos vectors
*/
double max_err_inf( matrix b, matrix x){
  matrix C;
  StartMatrix2(&C,b.row,b.col);
  int i;
  for(i=0; i< C.row; i++){
    C.data[i][0] = b.data[i][0]- x.data[i][0];
  }
  double value = norm_inf(C);
  //Alliberem C
  for( i =0; i<C.row ; i++){
    free(C.data[i]);
  }
    free(C.data);
  return value; 
}

/*
    Tindrem una funció f
*/
double function_f(double x){
    double pi = 3.141531415926535897932384;
    return 4*(pi*pi+1)*x*sin(2*pi*x)-4*pi*cos(2*pi*x);
}
/*
    Tindrem una funció u
*/
double function_u(double x){
    double pi = 3.141531415926535897932384;
    return x*sin(2*pi*x);
}

// MAIN
int  main(void){
    // Cal demanar al usuari la dimensio del problema
    int n=30;
    printf("\n Siusplau escriu la dimensio del problema, com a minim 2, amb un nombre enter  positiu i prem ENTER \n ");
    //scanf("%d", &n);
    // Ara cal fer la matriu A del problema.
    matrix A;
    StartMatrix2(&A,n-1,n-1);
    int i;
    int j;
    for( i=0; i<A.row; i++){
        for( j=0; j<A.col;j++){
            if( i == j ) A.data[i][j] = 2+4.0/(n*n);
            if( i == j+1) A.data[i][j] = -1;
            if( i+1 == j) A.data[i][j] = -1;
        }
    }
    // Ara necessitem computar les b_i=f_i
    matrix b;
    StartMatrix2(&b,n-1,1);

    for(i=0; i<b.row; i++){
        b.data[i][0]= function_f(1.0*(i+1)/n)/(n*n);
    }

    // Necessitem epsilon
    double epsilon = pow(10,-6);
    printf("\n Siusplau escriu l'error maxim acceptat i prem ENTER\n");
    //scanf("%lf",&epsilon);
    // Imposarem un nombre maxim d'iterats n_max
    int n_max = pow(10,5);
    printf("\n Siusplau escriu el maxim nombre d'iterats (enter positiu) i prem ENTER \n ");
    //scanf("%d", &n_max);
    // Necessitem un vector auxiliar per computar els radis espectrals
    matrix y;
    StartMatrix2(&y,n-1,1);
    double delta;
    double delta_aux;
    double theta;

    // Metode Jacobi
    matrix x_Jacobi;
    StartMatrix2(&x_Jacobi,n-1,1);

    delta = 0;
    delta_aux = 0;

    for( i=0; i<n_max; i++){
      // y = x_Jacobi
        for(j=0; j<y.row; j++){
          y.data[j][0] = x_Jacobi.data[j][0];
        }
      // Iterem Jacobi
        Jacobi(A,b,x_Jacobi);
      // Copiem delta al pas k
        delta_aux = delta;
      // Computem delta al pas k+1
        delta = max_err_inf(x_Jacobi,y);
      // Mirem si es compleix la condicio d'aturada
        if(absolut_err_inf(A,b,x_Jacobi)<epsilon){
          if(delta_aux !=0){
            theta = delta/delta_aux;
          }else{
            printf("\n Error delta==0 en J\n");
          }
          printf("\n Jacobi ha tardat %d  amb radi espectral %e",i+1,theta);
          break;
        }
    }
    // Guardem els resultats
    printf("\nSha creat el fitxer amb la solucio de Jacobi: Entrega1_Jacobi_Resposta.dat\n");
    PrintMatrixFile(x_Jacobi, "Entrega1_Jacobi_Resposta.dat");
    for( i =0; i<x_Jacobi.row ; i++){
      free(x_Jacobi.data[i]);
    }
      free(x_Jacobi.data);
    // Metode Gauss-Seidel
    matrix x_GS;
    StartMatrix2(&x_GS,n-1,1);
    // Resetejem delta
    theta = 0;
    delta = 0;
    delta_aux = 0;
    for( i=0; i<n_max; i++){
        // y = x_GS
        for(j=0; j<y.row; j++){
          y.data[j][0] = x_GS.data[j][0];
        }
        // Iterem GS
        Gauss_Seidel(A,b,x_GS);
        // Copiem delta al pas k
        delta_aux = delta;
        // Computem delta al pas k+1
        delta = max_err_inf(x_GS,y);
        // Mirem si es compleix la condicio d'aturada
        if(absolut_err_inf(A,b,x_GS)<epsilon){
          if(delta_aux !=0){
            theta = delta/delta_aux;
          }else{
            printf("\nError delta==0 en GS\n");
          }
          printf("\n Gauss-Seidel ha tardat %d amb radi espectral %e",i+1,theta);
          break;
        }
    }
    PrintMatrixFile(x_GS, "Entrega1_Gauss_Seidel_Resposta.dat");
    printf("\nSha creat el fitxer amb la solucio de Gauss-Seidel: Entrega1_Gauss_Seidel_Resposta.dat\n");

    for(i=0; i<x_GS.row ; i++){
      free(x_GS.data[i]);
    }
    free(x_GS.data);

    // Metode SOR
    
    FILE *fp;
    fp = fopen("Entrega1_w-i.dat","w");
    // El problema ens demana les u(x_i)
    matrix u;
    StartMatrix2(&u,n-1,1);
    for(i=0; i<u.row; i++){
        u.data[i][0]= function_u(1.0*(i+1)/n);
    }
    matrix x_SOR;
    StartMatrix2(&x_SOR,n-1,1);
    fprintf(fp,"w\t iterat \t e.max\n");
    double w;
    double w_min;
    int i_min = pow(10,3);
    for(j=0;j<1500;j++){
      // MODIFICAR w aqui i la j al for
        w=0.5 + 1.0*(j)/1000.0;
        printf("\n Actualment estem al pas %d de SOR amb w = %lf",j,w);
        for(i=0; i<n_max; i++){
            SOR(A,b,x_SOR,w);
            if(absolut_err_inf(A,b,x_SOR)<epsilon){
              if(i<= i_min){
                w_min = w;
                i_min = i;
              }
              fprintf(fp,"  %e \t %d\t %e\n",w,i+1, max_err_inf(x_SOR,u));
             break;
              }
        }
        for(i=0; i<x_SOR.row ; i++){
            x_SOR.data[i][0] = 0;
        }
    }
  printf("\nSha creat el fitxer amb els diversos valors iterats-w-error: Entrega1_w-i.dat\n");
  PrintMatrixFile(u,"Entrega1_u.dat");
  printf("\nEl valor de i_min = %d i el de  w_min = %e \n",i_min,w_min);
  for( i =0; i<A.row ; i++){
    free(A.data[i]);
    free(b.data[i]);
    free(x_SOR.data[i]);
  }
  free(A.data);
  free(b.data);
  free(x_SOR.data);

  return 0;
}